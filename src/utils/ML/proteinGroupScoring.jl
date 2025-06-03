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
    protein_groups::Dict,
    lib_precursors::LibraryPrecursors,
    n_top_precursors::Int=5
)
    # Group PSMs by precursor to get best score per precursor
    gpsms = groupby(psms, [:sequence, :charge, :proteins])
    precursor_scores = combine(gpsms, 
        :prob => maximum => :best_prob,
        :cv_fold => first => :cv_fold  # Keep CV fold assignment
    )
    
    # Initialize protein feature DataFrame
    protein_features = DataFrame()
    
    for (protein_key, pg_data) in protein_groups
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
        sort_idx = sortperm(peptide_scores, rev=true)
        peptide_scores = peptide_scores[sort_idx]
        
        # Extract top N scores (pad with zeros if fewer than N)
        top_scores = zeros(Float32, n_top_precursors)
        for i in 1:min(length(peptide_scores), n_top_precursors)
            top_scores[i] = peptide_scores[i]
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
        
        # Create feature row
        row = Dict{Symbol, Any}(
            :protein_group => protein_key,
            :is_target => is_target,
            :n_peptides => n_peptides,
            :n_precursors => n_precursors,
            :mean_precursor_score => mean_score,
            :median_precursor_score => median_score,
            :std_precursor_score => std_score,
            :current_score => pg_score,  # Current aggregated score
            :cv_fold => cv_fold
        )
        
        # Add top N scores as individual features
        for i in 1:n_top_precursors
            row[Symbol("top_score_", i)] = top_scores[i]
        end
        
        push!(protein_features, row)
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
    unique_folds = unique(protein_features.cv_fold)
    
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
    
    # Random forest parameters
    rf_params = merge(
        Dict(
            "objective" => "binary:logistic",
            "eval_metric" => "logloss",
            "num_parallel_tree" => num_parallel_tree,
            "subsample" => subsample,
            "colsample_bynode" => colsample_bynode
        ),
        XGBoost.randomforest()...  # Get default random forest parameters
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
            @warn "Skipping fold $test_fold due to insufficient data"
            continue
        end
        
        # Prepare training data
        X_train = Matrix(train_data[:, base_features])
        y_train = Float32.(train_data.is_target)
        
        X_test = Matrix(test_data[:, base_features])
        y_test = Float32.(test_data.is_target)
        
        # Create DMatrix
        dtrain = DMatrix(X_train', label=y_train)
        dtest = DMatrix(X_test', label=y_test)
        
        # Train random forest model
        booster = xgboost(
            dtrain, 
            num_round=n_rounds,  # Typically 1 for random forest
            param=rf_params,
            watchlist=[(dtrain, "train"), (dtest, "test")],
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
    protein_groups::Dict,
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
    for (i, row) in enumerate(eachrow(protein_features))
        protein_key = row.protein_group
        pg_score, peptides = protein_groups[protein_key]
        
        # Store both original and ML scores
        protein_groups[protein_key] = (
            pg_score = scores[i],  # Use ML score as primary
            peptides = peptides,
            ml_score = scores[i],
            original_score = pg_score
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