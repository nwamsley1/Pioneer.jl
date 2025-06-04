"""
Memory-efficient protein group scoring using top N precursor scores as features with cross-validation

This module implements out-of-memory machine learning-based protein group scoring that:
1. Samples protein groups for training (similar to PSM sampling)
2. Trains XGBoost random forest models on the sample
3. Applies models file-by-file to avoid loading all data into memory
4. Maintains cross-validation consistency with precursor scoring
"""

using DataFrames
using XGBoost
using StatsBase
using Statistics
using Random
using Arrow

"""
    count_protein_groups(protein_groups_folder::String)::Int

Count total protein groups across all Arrow files in the folder.
"""
function count_protein_groups(protein_groups_folder::String)::Int
    file_paths = [path for path in readdir(protein_groups_folder, join=true) if endswith(path, ".arrow")]
    total_count = 0
    
    for file_path in file_paths
        if isfile(file_path)
            total_count += length(Arrow.Table(file_path)[1])
        end
    end
    
    return total_count
end

"""
    sample_protein_groups_for_ml(
        protein_groups_folder::String,
        passing_psms_paths::Vector{String},
        precursors::LibraryPrecursors,
        max_proteins::Int;
        n_top_precursors::Int = 5,
        random_seed::Int = 1776
    )

Sample protein groups from multiple files for ML model training.
Uses proportional sampling similar to PSM sampling.

# Arguments
- `protein_groups_folder`: Folder containing protein group Arrow files
- `passing_psms_paths`: Paths to PSM files for extracting precursor scores
- `precursors`: Library precursor information for CV fold assignment
- `max_proteins`: Maximum number of protein groups to sample
- `n_top_precursors`: Number of top precursor scores to use as features
- `random_seed`: Random seed for reproducible sampling

# Returns
- DataFrame with sampled protein group features ready for ML training
"""
function sample_protein_groups_for_ml(
    protein_groups_folder::String,
    passing_psms_paths::Vector{String},
    precursors::LibraryPrecursors,
    max_proteins::Int;
    n_top_precursors::Int = 5,
    random_seed::Int = 1776
)
    @info "Sampling protein groups for ML training (max: $max_proteins)"
    
    # Get all protein group files
    pg_file_paths = [path for path in readdir(protein_groups_folder, join=true) if endswith(path, ".arrow")]
    
    if isempty(pg_file_paths)
        @error "No protein group files found in $protein_groups_folder"
        return DataFrame()
    end
    
    # Count total protein groups
    total_proteins = count_protein_groups(protein_groups_folder)
    @info "Total protein groups: $total_proteins"
    
    if total_proteins == 0
        @error "No protein groups found"
        return DataFrame()
    end
    
    # Load PSM data for extracting precursor scores
    @info "Loading PSM data for precursor scores..."
    psms_df = DataFrame()
    for psm_path in passing_psms_paths
        if isfile(psm_path)
            psm_data = DataFrame(Arrow.Table(psm_path))
            # Add CV fold information
            psm_data[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psm_data.precursor_idx]
            # Add sequence information
            psm_data[!, :sequence] = [getSequence(precursors)[pid] for pid in psm_data.precursor_idx]
            append!(psms_df, psms_data)
        end
    end
    
    if nrow(psms_df) == 0
        @error "No PSM data found"
        return DataFrame()
    end
    
    # Group PSMs by precursor to get best score per precursor
    gpsms = groupby(psms_df, [:sequence, :charge, :proteins])
    precursor_scores = combine(gpsms, 
        :prob => maximum => :best_prob,
        :cv_fold => first => :cv_fold
    )
    
    @info "Processed $(nrow(precursor_scores)) unique precursors"
    
    # Sample protein groups proportionally from each file
    Random.seed!(random_seed)
    sampled_features = DataFrame()
    
    for pg_file_path in pg_file_paths
        if !isfile(pg_file_path)
            continue
        end
        
        # Load protein groups from this file
        pg_table = Arrow.Table(pg_file_path)
        num_proteins = length(pg_table[1])
        
        if num_proteins == 0
            continue
        end
        
        # Calculate sample size proportional to file size
        sample_size = min(
            ceil(Int, (num_proteins / total_proteins) * max_proteins),
            num_proteins
        )
        
        if sample_size == 0
            continue
        end
        
        # Generate random indices for sampling
        sampled_indices = sort!(sample(1:num_proteins, sample_size, replace=false))
        
        # Sample and convert to DataFrame
        pg_df = DataFrame(pg_table)[sampled_indices, :]
        
        @info "Sampled $sample_size protein groups from $(basename(pg_file_path))"
        
        # Prepare features for sampled protein groups
        features_df = prepare_protein_features_from_pg_data(
            pg_df, precursor_scores, precursors, n_top_precursors
        )
        
        append!(sampled_features, features_df)
    end
    
    @info "Total sampled protein groups: $(nrow(sampled_features))"
    
    # Check target/decoy distribution
    if nrow(sampled_features) > 0
        target_dist = countmap(sampled_features.is_target)
        @info "Sampled target distribution: $target_dist"
        
        # Ensure we have both targets and decoys
        if length(target_dist) < 2 || get(target_dist, false, 0) == 0
            @warn "Insufficient decoys in sample for ML training"
            return DataFrame()
        end
    end
    
    return sampled_features
end

"""
    prepare_protein_features_from_pg_data(
        pg_df::DataFrame,
        precursor_scores::DataFrame,
        precursors::LibraryPrecursors,
        n_top_precursors::Int
    )

Prepare ML features from protein group data and precursor scores.
"""
function prepare_protein_features_from_pg_data(
    pg_df::DataFrame,
    precursor_scores::DataFrame,
    precursors::LibraryPrecursors,
    n_top_precursors::Int
)
    if nrow(pg_df) == 0
        return DataFrame()
    end
    
    # Extract unique protein groups and their associated peptides
    # This requires knowledge of how protein groups are stored in the files
    # For now, assume we have protein_name, target, entrap_id columns
    
    features_df = DataFrame()
    
    for row in eachrow(pg_df)
        # Extract protein group information
        protein_name = row.protein_name
        is_target = row.target
        entrap_id = row.entrap_id
        current_score = row.global_pg_score  # or pg_score
        
        # For this implementation, we'll need to infer peptides from PSM data
        # In a real implementation, we'd need the peptides stored with protein groups
        # or reconstruct them from the PSM data
        
        # Find peptides for this protein from precursor scores
        protein_peptides = precursor_scores[
            precursor_scores.proteins .== protein_name, :
        ]
        
        if nrow(protein_peptides) == 0
            continue
        end
        
        # Extract scores and sort
        peptide_scores = sort(protein_peptides.best_prob, rev=true)
        
        # Get top N scores (pad with zeros if needed)
        top_scores = zeros(Float32, n_top_precursors)
        for i in 1:min(length(peptide_scores), n_top_precursors)
            top_scores[i] = peptide_scores[i]
        end
        
        # Calculate protein-level features
        n_peptides = length(unique(protein_peptides.sequence))
        n_precursors = nrow(protein_peptides)
        mean_score = length(peptide_scores) > 0 ? mean(peptide_scores) : 0.0f0
        median_score = length(peptide_scores) > 0 ? median(peptide_scores) : 0.0f0
        std_score = length(peptide_scores) > 1 ? std(peptide_scores) : 0.0f0
        
        # Determine CV fold (use most common among peptides)
        cv_folds = protein_peptides.cv_fold
        if !isempty(cv_folds)
            fold_counts = countmap(cv_folds)
            max_count = maximum(values(fold_counts))
            most_common_folds = [f for (f, c) in fold_counts if c == max_count]
            cv_fold = length(most_common_folds) == 1 ? most_common_folds[1] : rand(most_common_folds)
        else
            cv_fold = rand(UInt8[0, 1])
        end
        
        # Create feature row
        feature_row = Dict{Symbol, Any}(
            :protein_name => protein_name,
            :is_target => is_target,
            :entrap_id => entrap_id,
            :n_peptides => n_peptides,
            :n_precursors => n_precursors,
            :mean_precursor_score => mean_score,
            :median_precursor_score => median_score,
            :std_precursor_score => std_score,
            :current_score => current_score,
            :cv_fold => cv_fold
        )
        
        # Add top score features
        for i in 1:n_top_precursors
            feature_row[Symbol("top_score_", i)] = top_scores[i]
        end
        
        push!(features_df, feature_row)
    end
    
    return features_df
end

"""
    train_protein_group_models_on_sample(
        sampled_features::DataFrame;
        n_rounds::Int = 1,
        n_top_precursors::Int = 5,
        num_parallel_tree::Int = 100,
        subsample::Float64 = 0.5,
        colsample_bynode::Float64 = 0.3
    )

Train XGBoost random forest models on sampled protein group features.
"""
function train_protein_group_models_on_sample(
    sampled_features::DataFrame;
    n_rounds::Int = 1,
    n_top_precursors::Int = 5,
    num_parallel_tree::Int = 100,
    subsample::Float64 = 0.5,
    colsample_bynode::Float64 = 0.3
)
    @info "Training protein group models on $(nrow(sampled_features)) sampled protein groups"
    
    if nrow(sampled_features) == 0
        @warn "No protein features available for training"
        return Dict{UInt8, Booster}()
    end
    
    # Check target distribution
    target_dist = countmap(sampled_features.is_target)
    @info "Training target distribution: $target_dist"
    
    if length(target_dist) < 2 || get(target_dist, false, 0) == 0
        @warn "Cannot train binary classifier: insufficient decoys in sample"
        return Dict{UInt8, Booster}()
    end
    
    unique_folds = unique(sampled_features.cv_fold)
    @info "Training on CV folds: $unique_folds"
    
    # Define feature columns
    base_features = [
        :n_peptides, :n_precursors,
        :mean_precursor_score, :median_precursor_score, :std_precursor_score
    ]
    for i in 1:n_top_precursors
        push!(base_features, Symbol("top_score_", i))
    end
    
    # Random forest parameters
    rf_params = Dict{String, Any}(
        "objective" => "binary:logistic",
        "eval_metric" => "logloss",
        "num_parallel_tree" => num_parallel_tree,
        "subsample" => subsample,
        "colsample_bynode" => colsample_bynode,
        "colsample_bytree" => 1.0,
        "learning_rate" => 1.0,
        "max_depth" => 6,
        "min_child_weight" => 1,
        "gamma" => 0.0,
        "reg_alpha" => 0.0,
        "reg_lambda" => 1.0
    )
    
    # Train models for each fold
    models = Dict{UInt8, Booster}()
    
    for test_fold in unique_folds
        train_mask = sampled_features.cv_fold .!= test_fold
        test_mask = sampled_features.cv_fold .== test_fold
        
        train_data = sampled_features[train_mask, :]
        test_data = sampled_features[test_mask, :]
        
        if nrow(train_data) < 10 || nrow(test_data) < 1
            @warn "Skipping fold $test_fold due to insufficient data (train: $(nrow(train_data)), test: $(nrow(test_data)))"
            continue
        end
        
        @info "Training fold $test_fold with $(nrow(train_data)) samples"
        
        # Prepare training data
        X_train = Matrix(train_data[:, base_features])
        y_train = Float32.(train_data.is_target)
        
        X_test = Matrix(test_data[:, base_features])
        y_test = Float32.(test_data.is_target)
        
        # Check for issues
        if any(isnan, X_train) || any(ismissing, X_train)
            @error "Training data contains NaN or missing values for fold $test_fold"
            continue
        end
        
        try
            # Create DMatrix and train model
            dtrain = DMatrix(X_train, label=y_train)
            dtest = DMatrix(X_test, label=y_test)
            
            booster = xgboost(
                dtrain,
                num_round=n_rounds,
                param=rf_params,
                verbose=0
            )
            
            models[test_fold] = booster
            @info "Successfully trained model for fold $test_fold"
            
        catch e
            @error "Failed to train model for fold $test_fold: $e"
            continue
        end
    end
    
    @info "Trained $(length(models)) protein group models"
    return models
end

"""
    apply_ml_scoring_to_protein_files!(
        protein_groups_paths::Vector{String},
        passing_psms_paths::Vector{String},
        models::Dict{UInt8, Booster},
        precursors::LibraryPrecursors;
        n_top_precursors::Int = 5
    )

Apply trained ML models to protein group files one by one to avoid memory issues.
"""
function apply_ml_scoring_to_protein_files!(
    protein_groups_paths::Vector{String},
    passing_psms_paths::Vector{String},
    models::Dict{UInt8, Booster},
    precursors::LibraryPrecursors;
    n_top_precursors::Int = 5
)
    @info "Applying ML scoring to $(length(protein_groups_paths)) protein group files"
    
    if isempty(models)
        @info "No ML models available, skipping ML scoring"
        return
    end
    
    # Load PSM data once for all files (this should be manageable in memory)
    @info "Loading PSM data for precursor scores..."
    psms_df = DataFrame()
    for psm_path in passing_psms_paths
        if isfile(psm_path)
            psm_data = DataFrame(Arrow.Table(psm_path))
            psm_data[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psm_data.precursor_idx]
            psm_data[!, :sequence] = [getSequence(precursors)[pid] for pid in psm_data.precursor_idx]
            append!(psms_df, psm_data)
        end
    end
    
    # Group PSMs by precursor
    gpsms = groupby(psms_df, [:sequence, :charge, :proteins])
    precursor_scores = combine(gpsms, 
        :prob => maximum => :best_prob,
        :cv_fold => first => :cv_fold
    )
    
    # Define feature columns
    feature_names = [
        :n_peptides, :n_precursors,
        :mean_precursor_score, :median_precursor_score, :std_precursor_score
    ]
    for i in 1:n_top_precursors
        push!(feature_names, Symbol("top_score_", i))
    end
    
    # Process each protein group file
    for (file_idx, pg_path) in enumerate(protein_groups_paths)
        if !isfile(pg_path)
            @warn "Protein group file not found: $pg_path"
            continue
        end
        
        try
            # Load protein groups for this file
            pg_df = DataFrame(Arrow.Table(pg_path))
            
            if nrow(pg_df) == 0
                @info "No protein groups in file $(basename(pg_path))"
                continue
            end
            
            @info "Processing $(nrow(pg_df)) protein groups in file $(basename(pg_path))"
            
            # Prepare features for this file
            features_df = prepare_protein_features_from_pg_data(
                pg_df, precursor_scores, precursors, n_top_precursors
            )
            
            if nrow(features_df) == 0
                @info "No features prepared for file $(basename(pg_path))"
                continue
            end
            
            # Apply ML scoring
            ml_scores = Vector{Float32}(undef, nrow(features_df))
            
            for (i, row) in enumerate(eachrow(features_df))
                fold = row.cv_fold
                
                if haskey(models, fold)
                    # Prepare features for this protein group
                    X = Matrix(features_df[i:i, feature_names])
                    
                    if !any(isnan, X) && !any(ismissing, X)
                        dtest = DMatrix(X)
                        pred = XGBoost.predict(models[fold], dtest)
                        ml_scores[i] = pred[1]
                    else
                        ml_scores[i] = row.current_score
                    end
                else
                    # Fallback to original score
                    ml_scores[i] = row.current_score
                end
            end
            
            # Update protein group scores in DataFrame
            if hasproperty(pg_df, :global_pg_score)
                pg_df[!, :global_pg_score] = ml_scores
            end
            if hasproperty(pg_df, :pg_score)
                # Update run-specific scores too
                pg_df[!, :pg_score] = ml_scores
            end
            
            # Write back the updated protein groups
            Arrow.write(pg_path, pg_df)
            @info "Updated ML scores for file $(basename(pg_path))"
            
        catch e
            @error "Failed to process protein group file $pg_path: $e"
            continue
        end
    end
    
    @info "Completed ML scoring application"
end

"""
    apply_ml_protein_scoring_oom!(
        protein_groups_folder::String,
        passing_psms_paths::Vector{String},
        passing_pg_paths::Vector{String},
        precursors::LibraryPrecursors;
        max_proteins_for_training::Int = 50000,
        n_top_precursors::Int = 5,
        num_parallel_tree::Int = 100,
        n_rounds::Int = 1,
        subsample::Float64 = 0.5,
        colsample_bynode::Float64 = 0.3
    )

Main function for out-of-memory protein group ML scoring.
"""
function apply_ml_protein_scoring_oom!(
    protein_groups_folder::String,
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    precursors::LibraryPrecursors;
    max_proteins_for_training::Int = 50000,
    n_top_precursors::Int = 5,
    num_parallel_tree::Int = 100,
    n_rounds::Int = 1,
    subsample::Float64 = 0.5,
    colsample_bynode::Float64 = 0.3
)
    @info "Starting out-of-memory protein group ML scoring"
    
    # Step 1: Sample protein groups for training
    sampled_features = sample_protein_groups_for_ml(
        protein_groups_folder,
        passing_psms_paths,
        precursors,
        max_proteins_for_training,
        n_top_precursors=n_top_precursors
    )
    
    if nrow(sampled_features) == 0
        @warn "No protein groups sampled for training, skipping ML scoring"
        return
    end
    
    # Step 2: Train models on sample
    models = train_protein_group_models_on_sample(
        sampled_features,
        n_rounds=n_rounds,
        n_top_precursors=n_top_precursors,
        num_parallel_tree=num_parallel_tree,
        subsample=subsample,
        colsample_bynode=colsample_bynode
    )
    
    if isempty(models)
        @warn "No models trained, skipping ML scoring application"
        return
    end
    
    # Step 3: Apply models to all files
    apply_ml_scoring_to_protein_files!(
        passing_pg_paths,
        passing_psms_paths,
        models,
        precursors,
        n_top_precursors=n_top_precursors
    )
    
    @info "Completed out-of-memory protein group ML scoring"
end