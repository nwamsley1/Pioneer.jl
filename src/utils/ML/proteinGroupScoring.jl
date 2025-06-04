"""
Protein group scoring using top N precursor scores as features with cross-validation

This module implements machine learning-based protein group scoring that uses
the scores of the top N precursors as features, trained with the same 
cross-validation scheme as the precursor scoring models using XGBoost random forests.
"""

using DataFrames
using XGBoost
using StatsBase
using Statistics
using Random
using Arrow
using Dictionaries

"""
    prepare_protein_group_features(
        psms::DataFrame, 
        protein_groups::Dict, 
        lib_precursors::LibraryPrecursors,
        n_top_precursors::Int=5
    )

Prepare features for protein group scoring by extracting top N precursor scores
and additional protein-level features.

# Arguments
- `psms`: DataFrame containing PSM data with scores and cv_fold assignments
- `protein_groups`: Dictionary mapping protein groups to their peptides
- `lib_precursors`: LibraryPrecursors object for accessing CV fold information
- `n_top_precursors`: Number of top scoring precursors to use as features

# Returns
- DataFrame with protein group features including CV fold assignments
"""
function prepare_protein_group_features(
    psms::DataFrame,
    protein_groups::Dictionary{
        @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, 
        @NamedTuple{pg_score::Float32, peptides::Set{String}}
    },
    lib_precursors::LibraryPrecursors,
    n_top_precursors::Int=5
)
    @info "Preparing protein group features with top $n_top_precursors precursor scores"
    # Group PSMs by precursor to get best score per precursor
    gpsms = groupby(psms, [:sequence, :charge, :proteins])
    precursor_scores = combine(gpsms, 
        :prob => maximum => :best_prob,
        :cv_fold => first => :cv_fold  # Keep CV fold assignment
    )
    
    # Pre-allocate vectors for efficient DataFrame construction
    n_proteins = length(protein_groups)
    
    # Initialize vectors for all columns
    protein_group_keys = Vector{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}}(undef, n_proteins)
    is_target_vec = Vector{Bool}(undef, n_proteins)
    n_peptides_vec = Vector{Int64}(undef, n_proteins)
    n_precursors_vec = Vector{Int64}(undef, n_proteins)
    mean_precursor_score_vec = Vector{Float32}(undef, n_proteins)
    median_precursor_score_vec = Vector{Float32}(undef, n_proteins)
    std_precursor_score_vec = Vector{Float32}(undef, n_proteins)
    current_score_vec = Vector{Float32}(undef, n_proteins)
    cv_fold_vec = Vector{UInt8}(undef, n_proteins)
    
    # Initialize top score vectors
    top_score_vecs = [Vector{Float32}(undef, n_proteins) for _ in 1:n_top_precursors]
    
    # Process each protein group
    for (idx, (protein_key, pg_data)) in enumerate(pairs(protein_groups))
        # Extract protein group data
        pg_score, peptides = pg_data
        protein_name, is_target, entrap_id = protein_key
        
        # Get all precursor scores and CV folds for this protein group
        peptide_scores = Float32[]
        peptide_cv_folds = UInt8[]
        
        for peptide in peptides
            matching_rows = precursor_scores[
                precursor_scores.sequence .== peptide, :
            ]
            if nrow(matching_rows) > 0
                append!(peptide_scores, matching_rows.best_prob)
                append!(peptide_cv_folds, matching_rows.cv_fold)
            end
        end
        
        # Sort scores in descending order
        if !isempty(peptide_scores)
            sort_idx = sortperm(peptide_scores, rev=true)
            peptide_scores = peptide_scores[sort_idx]
        end
        
        # Extract top N scores (pad with zeros if fewer than N)
        for i in 1:n_top_precursors
            if i <= length(peptide_scores)
                top_score_vecs[i][idx] = peptide_scores[i]
            else
                top_score_vecs[i][idx] = 0.0f0
            end
        end
        
        # Calculate additional protein-level features
        n_peptides = length(unique(peptides))
        n_precursors = length(peptide_scores)
        mean_score = length(peptide_scores) > 0 ? mean(peptide_scores) : 0.0f0
        median_score = length(peptide_scores) > 0 ? median(peptide_scores) : 0.0f0
        std_score = length(peptide_scores) > 1 ? std(peptide_scores) : 0.0f0
        
        # Determine CV fold for protein group
        # Use the most common fold among its peptides (consistent with LibraryIon.jl)
        if !isempty(peptide_cv_folds)
            fold_counts = countmap(peptide_cv_folds)
            max_count = maximum(values(fold_counts))
            most_common_folds = [f for (f, c) in fold_counts if c == max_count]
            cv_fold = length(most_common_folds) == 1 ? most_common_folds[1] : rand(most_common_folds)
        else
            # Fallback: randomly assign fold
            cv_fold = rand(UInt8[0, 1])
        end
        
        # Store values in vectors
        protein_group_keys[idx] = protein_key
        is_target_vec[idx] = is_target
        n_peptides_vec[idx] = n_peptides
        n_precursors_vec[idx] = n_precursors
        mean_precursor_score_vec[idx] = mean_score
        median_precursor_score_vec[idx] = median_score
        std_precursor_score_vec[idx] = std_score
        current_score_vec[idx] = pg_score
        cv_fold_vec[idx] = cv_fold
    end
    
    # Construct DataFrame efficiently from pre-allocated vectors
    protein_features = DataFrame(
        protein_group = protein_group_keys,
        is_target = is_target_vec,
        n_peptides = n_peptides_vec,
        n_precursors = n_precursors_vec,
        mean_precursor_score = mean_precursor_score_vec,
        median_precursor_score = median_precursor_score_vec,
        std_precursor_score = std_precursor_score_vec,
        current_score = current_score_vec,
        cv_fold = cv_fold_vec
    )
    
    # Add top score columns
    for i in 1:n_top_precursors
        protein_features[!, Symbol("top_score_", i)] = top_score_vecs[i]
    end
    
    return protein_features
end

"""
    train_protein_group_models(
        protein_features::DataFrame;
        n_rounds::Int = 1,  # For random forest, typically use 1 round
        n_top_precursors::Int = 5,
        num_parallel_tree::Int = 100,
        subsample::Float64 = 0.5,
        colsample_bynode::Float64 = 0.3
    )

Train XGBoost random forest models for protein group scoring using cross-validation.

# Arguments
- `protein_features`: DataFrame with protein group features and cv_fold assignments
- `n_rounds`: Number of boosting rounds (typically 1 for random forest)
- `n_top_precursors`: Number of top precursor scores used as features
- `num_parallel_tree`: Number of trees in the random forest
- `subsample`: Fraction of data to sample per tree
- `colsample_bynode`: Fraction of features to sample per node

# Returns
- Dictionary mapping CV fold indices to trained models
"""
function train_protein_group_models(
    protein_features::DataFrame;
    n_rounds::Int = 1,
    n_top_precursors::Int = 5,
    num_parallel_tree::Int = 100,
    subsample::Float64 = 0.5,
    colsample_bynode::Float64 = 0.3
)
    @info "Training protein group models with $(nrow(protein_features)) protein groups"
    @info "Protein features columns: $(names(protein_features))"
    
    target_dist = countmap(protein_features.is_target)
    @info "Target distribution: $target_dist"
    
    # Check if we have both targets and decoys
    if length(target_dist) < 2 || get(target_dist, false, 0) == 0
        @warn "Cannot train binary classifier: only targets found, no decoys. Skipping ML scoring."
        return Dict{UInt8, Booster}()
    end
    
    unique_folds = unique(protein_features.cv_fold)
    @info "CV folds found: $unique_folds"
    
    # Initialize model storage
    models = Dict{UInt8, Booster}()
    
    # Feature columns for training
    base_features = [
        :n_peptides, :n_precursors, 
        :mean_precursor_score, :median_precursor_score, :std_precursor_score
    ]
    
    # Add top score features
    for i in 1:n_top_precursors
        push!(base_features, Symbol("top_score_", i))
    end
    
    # Random forest parameters (manually set based on XGBoost.randomforest defaults)
    rf_params = Dict{String, Any}(
        "objective" => "binary:logistic",
        "eval_metric" => "logloss",
        "num_parallel_tree" => num_parallel_tree,
        "subsample" => subsample,
        "colsample_bynode" => colsample_bynode,
        # Random forest specific parameters
        "colsample_bytree" => 1.0,  # Use all features per tree
        "learning_rate" => 1.0,     # No shrinkage for random forest
        "max_depth" => 6,           # Default depth
        "min_child_weight" => 1,    # Default
        "gamma" => 0.0,             # No regularization
        "reg_alpha" => 0.0,         # No L1 regularization
        "reg_lambda" => 1.0         # Default L2 regularization
    )
    
    # Train models for each fold
    for test_fold in unique_folds
        # Split data
        train_mask = protein_features.cv_fold .!= test_fold
        test_mask = protein_features.cv_fold .== test_fold
        
        train_data = protein_features[train_mask, :]
        test_data = protein_features[test_mask, :]
        
        # Skip if not enough data
        if nrow(train_data) < 10 || nrow(test_data) < 1
            @warn "Skipping fold $test_fold due to insufficient data (train: $(nrow(train_data)), test: $(nrow(test_data)))"
            continue
        end
        
        @info "Fold $test_fold: Training with $(nrow(train_data)) samples, testing with $(nrow(test_data)) samples"
        
        # Prepare training data
        X_train = Matrix(train_data[:, base_features])
        y_train = Float32.(train_data.is_target)
        
        X_test = Matrix(test_data[:, base_features])
        y_test = Float32.(test_data.is_target)
        
        @info "Data shapes - X_train: $(size(X_train)), y_train: $(length(y_train)), X_test: $(size(X_test)), y_test: $(length(y_test))"
        @info "Features being used: $base_features"
        
        # Check for data consistency
        if size(X_train, 1) != length(y_train)
            @error "Mismatch: X_train has $(size(X_train, 1)) rows but y_train has $(length(y_train)) elements"
            continue
        end
        
        if size(X_test, 1) != length(y_test)
            @error "Mismatch: X_test has $(size(X_test, 1)) rows but y_test has $(length(y_test)) elements"
            continue
        end
        
        # Write training data to desktop for inspection
        train_debug_df = copy(train_data)
        train_debug_df[!, :fold] .= test_fold
        train_debug_df[!, :split] .= "train"
        test_debug_df = copy(test_data)
        test_debug_df[!, :fold] .= test_fold
        test_debug_df[!, :split] .= "test"
        debug_df = vcat(train_debug_df, test_debug_df)
        
        try
            @info size(train_debug_df)
            Arrow.write("/Users/nathanwamsley/Desktop/protein_training_data_fold_$(test_fold).arrow", train_debug_df)
            @info "Training data written to desktop for fold $test_fold"
        catch e
            @warn "Could not write debug data to desktop: $e"
        end
        
        # Check for missing values
        if any(ismissing, X_train) || any(isnan, X_train)
            @error "X_train contains missing or NaN values"
            @info "X_train summary: $(describe(DataFrame(X_train, base_features)))"
            continue
        end
        
        if any(ismissing, y_train) || any(isnan, y_train)
            @error "y_train contains missing or NaN values"
            continue
        end
        
        @info "Before DMatrix creation: X_train size = $(size(X_train)), y_train length = $(length(y_train))"
        
        # Initialize dtrain and dtest outside try-catch blocks
        dtrain = nothing
        dtest = nothing
        
        # Try creating DMatrix without transposition first (XGBoost.jl might handle this automatically)
        try
            @info "Attempting DMatrix creation without transposition..."
            dtrain = DMatrix(X_train, label=y_train)
            dtest = DMatrix(X_test, label=y_test)
            @info "DMatrix created successfully without transposition"
        catch e1
            @info "Failed without transposition: $e1"
            try
                @info "Attempting DMatrix creation with transposition..."
                dtrain = DMatrix(X_train', label=y_train)
                dtest = DMatrix(X_test', label=y_test)
                @info "DMatrix created successfully with transposition"
            catch e2
                @error "Both DMatrix creation methods failed!"
                @error "Without transposition: $e1"
                @error "With transposition: $e2"
                continue
            end
        end
        
        # Check if DMatrix creation was successful
        if dtrain === nothing || dtest === nothing
            @error "Failed to create DMatrix objects"
            continue
        end
        
        # Train random forest model
        # Remove watchlist for now to avoid the tuple format issue
        booster = xgboost(
            dtrain, 
            num_round=n_rounds,  # Typically 1 for random forest
            param=rf_params,
            verbose=0
        )
        
        models[test_fold] = booster
    end
    
    return models
end

"""
    score_protein_groups(
        protein_features::DataFrame,
        models::Dict{UInt8, Booster},
        feature_names::Vector{Symbol}
    )

Score protein groups using trained random forest models.

# Returns
- Vector of protein group scores
"""
function score_protein_groups(
    protein_features::DataFrame,
    models::Dict{UInt8, Booster},
    feature_names::Vector{Symbol}
)
    scores = zeros(Float32, nrow(protein_features))
    
    # If no models were trained, return current scores
    if isempty(models)
        @info "No ML models available, using original scores"
        return protein_features.current_score
    end
    
    for i in 1:nrow(protein_features)
        fold = protein_features[i, :cv_fold]
        
        if haskey(models, fold)
            # Prepare features
            X = Matrix(protein_features[i:i, feature_names])
            dtest = DMatrix(X')
            
            # Predict
            pred = XGBoost.predict(models[fold], dtest)
            scores[i] = pred[1]
        else
            # Fallback to current score if no model available
            scores[i] = protein_features[i, :current_score]
        end
    end
    
    return scores
end

"""
    integrate_protein_scoring!(
        psms::DataFrame,
        protein_groups::Dict,
        lib_precursors::LibraryPrecursors;
        n_top_precursors::Int = 5,
        num_parallel_tree::Int = 100,
        n_rounds::Int = 1,
        subsample::Float64 = 0.5,
        colsample_bynode::Float64 = 0.3
    )

Main function to integrate protein group scoring into the workflow using random forests.

This function should be called after initial PSM scoring but before final
protein group selection and FDR calculation.

# Arguments
- `psms`: DataFrame with scored PSMs including cv_fold column
- `protein_groups`: Dictionary of protein groups
- `lib_precursors`: LibraryPrecursors object containing CV fold information
- `n_top_precursors`: Number of top precursor scores to use as features
- `num_parallel_tree`: Number of trees in the random forest
- `n_rounds`: Boosting rounds (typically 1 for random forest)
- `subsample`: Fraction of data to sample per tree
- `colsample_bynode`: Fraction of features to sample per node

# Returns
- `protein_features`: DataFrame with protein features and scores
- `models`: Dictionary of trained models per CV fold
"""
function integrate_protein_scoring!(
    psms::DataFrame,
    protein_groups::Dictionary{
        @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, 
        @NamedTuple{pg_score::Float32, peptides::Set{String}}
    },
    lib_precursors::LibraryPrecursors;
    n_top_precursors::Int = 5,
    num_parallel_tree::Int = 100,
    n_rounds::Int = 1,
    subsample::Float64 = 0.5,
    colsample_bynode::Float64 = 0.3
)
    # Prepare protein features
    protein_features = prepare_protein_group_features(
        psms, protein_groups, lib_precursors, n_top_precursors
    )
    
    # Train random forest models
    models = train_protein_group_models(
        protein_features,
        n_rounds=n_rounds,
        n_top_precursors=n_top_precursors,
        num_parallel_tree=num_parallel_tree,
        subsample=subsample,
        colsample_bynode=colsample_bynode
    )
    
    # Score protein groups
    feature_names = [
        :n_peptides, :n_precursors, 
        :mean_precursor_score, :median_precursor_score, :std_precursor_score
    ]
    for i in 1:n_top_precursors
        push!(feature_names, Symbol("top_score_", i))
    end
    
    scores = score_protein_groups(protein_features, models, feature_names)
    
    # Add ML scores to protein features
    protein_features[!, :ml_score] = scores
    
    # Update protein groups with new scores
    # Keep the original structure - only update pg_score
    for (i, row) in enumerate(eachrow(protein_features))
        protein_key = row.protein_group
        pg_score, peptides = protein_groups[protein_key]
        
        # Update with ML score (or original if no ML model was used)
        protein_groups[protein_key] = (
            pg_score = scores[i],
            peptides = peptides
        )
    end
    
    return protein_features, models
end

"""
    export_protein_features(protein_features::DataFrame, output_path::String)

Export protein features for analysis or debugging.
"""
function export_protein_features(protein_features::DataFrame, output_path::String)
    Arrow.write(output_path, protein_features)
end

# ============================================================================
# Memory-Efficient Out-of-Memory (OOM) Protein Group Scoring Functions
# ============================================================================

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
This is a simplified version that uses the existing in-memory approach
but with better memory management.
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
    
    # Step 1: Sample protein groups for training from across files (similar to PSM sampling)
    sampled_protein_features = sample_protein_groups_for_training(
        passing_psms_paths,
        passing_pg_paths,
        precursors,
        max_proteins_for_training,
        n_top_precursors
    )
    
    if nrow(sampled_protein_features) == 0
        @warn "No protein groups sampled for training, keeping original scores"
        return
    end
    
    # Step 2: Train ML models on the sample
    models = train_protein_ml_models(
        sampled_protein_features,
        n_top_precursors,
        num_parallel_tree,
        n_rounds,
        subsample,
        colsample_bynode
    )
    
    if isempty(models)
        @warn "No models trained, keeping original scores"
        return
    end
    
    # Step 3: Apply models file-by-file (true OOM approach)
    apply_protein_ml_models_oom!(
        passing_psms_paths,
        passing_pg_paths,
        models,
        precursors,
        n_top_precursors
    )
    
    @info "Completed out-of-memory protein group ML scoring"
end

"""
Sample protein groups for ML training using proportional sampling across files.
This ensures we get a representative sample without loading all data into memory.
"""
function sample_protein_groups_for_training(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    precursors::LibraryPrecursors,
    max_proteins::Int,
    n_top_precursors::Int
)
    @info "Sampling protein groups for ML training (max: $max_proteins)"
    
    # Count total protein groups across all files
    total_proteins = 0
    @info passing_pg_paths "Found $(length(passing_pg_paths)) protein group files"
    for pg_path in passing_pg_paths
        if isfile(pg_path)
            try
                pg_table = Arrow.Table(pg_path)
                total_proteins += length(pg_table.protein_name)
            catch
                continue
            end
        end
    end
    
    if total_proteins == 0
        @warn "No protein groups found for sampling"
        return DataFrame()
    end
    
    @info "Found $total_proteins total protein groups across $(length(passing_pg_paths)) files"
    
    # Sample proportionally from each file
    sampled_features = DataFrame()
    Random.seed!(1776)
    
    for (file_idx, pg_path) in enumerate(passing_pg_paths)
        psm_path = passing_psms_paths[file_idx]
        
        if !isfile(pg_path) || !isfile(psm_path)
            continue
        end
        
        try
            # Load protein groups for this file
            pg_table = Arrow.Table(pg_path)
            num_proteins = length(pg_table.protein_name)
            
            if num_proteins == 0
                continue
            end
            
            # Calculate proportional sample size
            sample_size = min(
                ceil(Int, (num_proteins / total_proteins) * max_proteins),
                num_proteins
            )
            
            if sample_size == 0
                continue
            end
            
            # Sample indices
            sampled_indices = sort!(sample(1:num_proteins, sample_size, replace=false))
            
            # Load PSMs for this file to extract features
            psm_table = Arrow.Table(psm_path)
            psm_df = DataFrame(psm_table)
            
            # Add necessary columns for feature extraction
            psm_df[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psm_df.precursor_idx]
            psm_df[!, :sequence] = [getSequence(precursors)[pid] for pid in psm_df.precursor_idx]
            psm_df[!, :charge] = [getCharge(precursors)[pid] for pid in psm_df.precursor_idx]
            psm_df[!, :proteins] = [getAccessionNumbers(precursors)[pid] for pid in psm_df.precursor_idx]
            
            # Extract features for sampled protein groups
            file_features = extract_protein_features_from_file(
                pg_table, sampled_indices, psm_df, precursors, n_top_precursors
            )
            
            append!(sampled_features, file_features)
            
            @info "Sampled $sample_size protein groups from $(basename(pg_path))"
            
        catch e
            @warn "Failed to sample from file $pg_path: $e"
            continue
        end
    end
    
    @info "Total sampled protein groups: $(nrow(sampled_features))"
    return sampled_features
end

"""
Extract protein group features from a single file's PSMs and protein groups.
"""
function extract_protein_features_from_file(
    pg_table::Arrow.Table,
    sampled_indices::Vector{Int},
    psm_df::DataFrame,
    precursors::LibraryPrecursors,
    n_top_precursors::Int
)
    # Initialize with proper schema instead of empty DataFrame
    features_df = DataFrame(
        protein_name = Vector{String}(undef, length(sampled_indices)),
        is_target = Vector{Bool}(undef, length(sampled_indices)),
        entrap_id = Vector{UInt8}(undef, length(sampled_indices)),
        top_score_1 = Vector{Float32}(undef, length(sampled_indices)),
        top_score_2 = Vector{Float32}(undef, length(sampled_indices)),
        top_score_3 = Vector{Float32}(undef, length(sampled_indices)),
        top_score_4 = Vector{Float32}(undef, length(sampled_indices)),
        top_score_5 = Vector{Float32}(undef, length(sampled_indices)),
        n_peptides = Vector{Int}(undef, length(sampled_indices)),
        n_precursors = Vector{Int}(undef, length(sampled_indices)),
        mean_precursor_score = Vector{Float32}(undef, length(sampled_indices)),
        median_precursor_score = Vector{Float32}(undef, length(sampled_indices)),
        std_precursor_score = Vector{Float32}(undef, length(sampled_indices)),
        current_score = Vector{Float32}(undef, length(sampled_indices)),
        cv_fold = Vector{UInt8}(undef, length(sampled_indices))
    )
    for (i, idx) in enumerate(sampled_indices)
        protein_name = pg_table.protein_name[idx]
        features_df[i,:is_target] = pg_table.target[idx]
        features_df[i,:entrap_id]= pg_table.entrap_id[idx]
        features_df[i,:current_score] = pg_table.pg_score[idx]
        
        # Find PSMs for this protein
        protein_psms = filter(row -> row.proteins == protein_name, psm_df)
        
        if nrow(protein_psms) == 0
            continue
        end
        
        # Extract precursor scores and sort
        precursor_scores = sort(protein_psms.prob, rev=true)
        
        # Get top N scores (pad with zeros if needed)
        top_scores = zeros(Float32, n_top_precursors)
        for i in 1:min(length(precursor_scores), n_top_precursors)
            top_scores[i] = precursor_scores[i]
        end
        
        # Calculate protein-level features
        features_df[i,:n_peptides] = length(unique(protein_psms.sequence))
        features_df[i,:n_precursors] = nrow(protein_psms)
        features_df[i,:mean_precursor_score] = mean(precursor_scores)
        features_df[i,:median_precursor_score] = median(precursor_scores)
        features_df[i,:std_precursor_score] = length(precursor_scores) > 1 ? std(precursor_scores) : 0.0f0
        
        # Get CV fold (most common among peptides)
        cv_folds = protein_psms.cv_fold
        features_df[i,:cv_fold] = isempty(cv_folds) ? UInt8(0) : mode(cv_folds)
        
        # Add top score features
        for i in 1:n_top_precursors
            features_df[i,Symbol("top_score_", i)] = top_scores[i]
        end
        
        #push!(features_df, feature_row)
    end
    
    return features_df
end

"""
Train ML models on sampled protein features using cross-validation.
"""
function train_protein_ml_models(
    sampled_features::DataFrame,
    n_top_precursors::Int,
    num_parallel_tree::Int,
    n_rounds::Int,
    subsample::Float64,
    colsample_bynode::Float64
)
    @info "Training protein ML models on $(nrow(sampled_features)) samples"
    
    if nrow(sampled_features) < 100
        @warn "Insufficient data for ML training"
        return Dict{UInt8, Any}()
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
    
    # XGBoost random forest parameters
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
    models = Dict{UInt8, Any}()
    
    for test_fold in unique_folds
        train_mask = sampled_features.cv_fold .!= test_fold
        train_data = sampled_features[train_mask, :]
        
        if nrow(train_data) < 50
            @warn "Skipping fold $test_fold due to insufficient training data"
            continue
        end
        
        # Prepare training data
        X_train = Matrix(train_data[:, base_features])
        y_train = Float32.(train_data.is_target)
        
        try
            dtrain = DMatrix(X_train, label=y_train)
            booster = xgboost(dtrain, num_round=n_rounds, param=rf_params, verbose=0)
            models[test_fold] = booster
            @info "Successfully trained model for fold $test_fold"
        catch e
            @warn "Failed to train model for fold $test_fold: $e"
        end
    end
    
    @info "Trained $(length(models)) protein ML models"
    return models
end

"""
Apply trained ML models to protein group files one by one (true OOM approach).
"""
function apply_protein_ml_models_oom!(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    models::Dict{UInt8, Any},
    precursors::LibraryPrecursors,
    n_top_precursors::Int
)
    @info "Applying ML models to $(length(passing_pg_paths)) protein group files"
    
    if isempty(models)
        @info "No models available, keeping original scores"
        return
    end
    
    # Define feature columns
    feature_names = [
        :n_peptides, :n_precursors,
        :mean_precursor_score, :median_precursor_score, :std_precursor_score
    ]
    for i in 1:n_top_precursors
        push!(feature_names, Symbol("top_score_", i))
    end
    
    # Process each file individually
    for (file_idx, pg_path) in enumerate(passing_pg_paths)
        psm_path = passing_psms_paths[file_idx]
        
        if !isfile(pg_path) || !isfile(psm_path)
            @warn "Missing files for file $file_idx, skipping"
            continue
        end
        
        try
            # Load protein groups for this file
            pg_df = DataFrame(Arrow.Table(pg_path))
            
            if nrow(pg_df) == 0
                @info "No protein groups in $(basename(pg_path)), skipping"
                continue
            end
            
            # Load PSMs for this file
            psm_table = Arrow.Table(psm_path)
            psm_df = DataFrame(psm_table)
            
            # Add necessary columns
            psm_df[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psm_df.precursor_idx]
            psm_df[!, :sequence] = [getSequence(precursors)[pid] for pid in psm_df.precursor_idx]
            psm_df[!, :proteins] = [getAccessionNumbers(precursors)[pid] for pid in psm_df.precursor_idx]
            
            # Extract features and apply models
            ml_scores = Vector{Float32}(undef, nrow(pg_df))
            
            for (i, row) in enumerate(eachrow(pg_df))
                protein_name = row.protein_name
                
                # Find PSMs for this protein
                protein_psms = filter(psm_row -> psm_row.proteins == protein_name, psm_df)
                
                if nrow(protein_psms) == 0
                    ml_scores[i] = row.pg_score  # Keep original score
                    continue
                end
                
                # Extract features
                precursor_scores = sort(protein_psms.prob, rev=true)
                top_scores = zeros(Float32, n_top_precursors)
                for j in 1:min(length(precursor_scores), n_top_precursors)
                    top_scores[j] = precursor_scores[j]
                end
                
                n_peptides = length(unique(protein_psms.sequence))
                n_precursors = nrow(protein_psms)
                mean_score = mean(precursor_scores)
                median_score = median(precursor_scores)
                std_score = length(precursor_scores) > 1 ? std(precursor_scores) : 0.0f0
                
                cv_fold = isempty(protein_psms.cv_fold) ? UInt8(0) : mode(protein_psms.cv_fold)
                
                # Apply model if available for this fold
                if haskey(models, cv_fold)
                    # Prepare features
                    feature_values = [
                        n_peptides, n_precursors, mean_score, median_score, std_score,
                        top_scores...
                    ]
                    
                    X = reshape(Float32.(feature_values), 1, :)
                    
                    try
                        dtest = DMatrix(X)
                        pred = XGBoost.predict(models[cv_fold], dtest)
                        ml_scores[i] = pred[1]
                    catch e
                        @warn "Model prediction failed for protein $protein_name: $e"
                        ml_scores[i] = row.pg_score
                    end
                else
                    ml_scores[i] = row.pg_score  # Keep original score
                end
            end
            
            # Update protein group scores
            pg_df[!, :pg_score] = ml_scores
            if hasproperty(pg_df, :global_pg_score)
                pg_df[!, :global_pg_score] = ml_scores
            end
            
            # Write back updated protein groups
            Arrow.write(pg_path, pg_df)
            @info "Updated ML scores for $(nrow(pg_df)) protein groups in $(basename(pg_path))"
            
            # Free memory
            pg_df = nothing
            psm_df = nothing
            psm_table = nothing
            GC.gc()
            
        catch e
            @warn "Failed to process file $(basename(pg_path)): $e"
            continue
        end
    end
    
    @info "Completed OOM ML scoring application"
end