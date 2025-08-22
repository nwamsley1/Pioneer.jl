# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#==========================================================
Model Comparison Framework for ScoringSearch
==========================================================#

using DataFrames, Random, StatsBase, CSV, FileIO, JLD2
using EvoTrees
using ..Pioneer: @user_warn, @user_info, get_qvalues!, sort_of_percolator_in_memory!, probit_regression_scoring_cv!

"""
Configuration for a single model in the comparison framework.

# Fields
- `name`: Model identifier (e.g., "SimpleXGBoost", "ProbitRegression")
- `model_type`: Algorithm type (:xgboost or :probit)
- `features`: Vector of feature symbols to use
- `hyperparams`: Dictionary of hyperparameters for the model
"""
struct ModelConfig
    name::String
    model_type::Symbol  # :xgboost or :probit
    features::Vector{Symbol}
    hyperparams::Dict{Symbol, Any}
end

"""
Performance metrics for a trained model on validation data.

# Fields
- `model_name`: Name of the model
- `n_targets_passing_qval`: Number of targets passing q-value threshold (primary metric)
- `validation_auc`: Area under ROC curve on validation set
- `validation_accuracy`: Classification accuracy on validation set
- `validation_sensitivity`: True positive rate
- `validation_specificity`: True negative rate
- `training_time`: Time taken to train the model (seconds)
- `n_features`: Number of features used by the model
"""
struct ModelPerformance
    model_name::String
    n_targets_passing_qval::Int64
    validation_auc::Float64
    validation_accuracy::Float64
    validation_sensitivity::Float64
    validation_specificity::Float64
    training_time::Float64
    n_features::Int
end

"""
Result container for a trained model including metadata.

# Fields
- `model`: Trained model object (EvoTrees models or probit coefficients)
- `model_config`: Configuration used to train this model
- `training_time`: Time taken for training
- `n_features`: Number of features actually used
- `trained_features`: List of features used (may differ from config if some unavailable)
"""
struct ModelResult
    model::Any
    model_config::ModelConfig
    training_time::Float64
    n_features::Int
    trained_features::Vector{Symbol}
end

# Feature set definitions as per the implementation plan
const REDUCED_FEATURE_SET = [
    # Core peptide properties
    :missed_cleavage, :Mox, :prec_mz, :sequence_length, :charge,
    # RT features
    :irt_pred, :irt_error, :irt_diff, :rt_diff, :ms1_irt_diff,
    # Spectral features
    :max_y_ions, :y_ions_sum, :longest_y, :y_count, :b_count, :isotope_count,
    :total_ions, :best_rank, :best_rank_iso, :topn, :topn_iso, :gof,
    # Quality metrics
    :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_residual, :max_unmatched_residual, :max_gof,
    :fitted_spectral_contrast, :spectral_contrast, :max_matched_ratio,
    :err_norm, :poisson, :weight, :log2_intensity_explained, :tic, :num_scans,
    :smoothness, :percent_theoretical_ignored, :scribe, :max_scribe,
    # MS1 features
    :weight_ms1, :gof_ms1, :max_matched_residual_ms1, :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1, :error_ms1, :m0_error_ms1, :n_iso_ms1,
    :big_iso_ms1, :ms1_features_missing
    # MBR features added automatically if match_between_runs=true
]

const MINIMAL_FEATURE_SET = [
    :fitted_spectral_contrast,
    :max_matched_residual,
    :max_unmatched_residual,
    :err_norm,
    :log2_intensity_explained
]

"""
    create_model_configurations() -> Vector{ModelConfig}

Creates the three model configurations for comparison as defined in the implementation plan.

# Returns
- Vector of ModelConfig objects for SimpleXGBoost, ProbitRegression, and SuperSimplified models
"""
function create_model_configurations()
    return [
        # Model 1: Simple XGBoost (Current Default)
        ModelConfig(
            "SimpleXGBoost",
            :xgboost,
            REDUCED_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.8,
                :colsample_bynode => 0.8,
                :min_child_weight => 20,
                :gamma => 0.1,
                :subsample => 0.8,
                :max_depth => 4,
                :eta => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        ),
        
        # Model 2: Probit Regression
        ModelConfig(
            "ProbitRegression",
            :probit,
            REDUCED_FEATURE_SET,
            Dict(
                :n_folds => 3,
                :max_iter => 30
            )
        ),
        
        # Model 3: Super Simplified Model
        ModelConfig(
            "SuperSimplified",
            :xgboost,
            MINIMAL_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.8,
                :colsample_bynode => 0.8,
                :min_child_weight => 20,
                :gamma => 0.1,
                :subsample => 0.8,
                :max_depth => 4,
                :eta => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        )
    ]
end

#==========================================================
Data Splitting Framework
==========================================================#

"""
    create_train_validation_split(psms::DataFrame, validation_ratio::Float64 = 0.2; seed::Int = 1776) -> (Vector{Int}, Vector{Int})

Creates simple 80/20 split with basic stratification by target/decoy status.

# Arguments
- `psms`: DataFrame containing PSMs with target column
- `validation_ratio`: Fraction of data to use for validation (default: 0.2)
- `seed`: Random seed for reproducible splits (default: 1776)

# Returns
- `train_indices`: Indices for training set
- `validation_indices`: Indices for validation set
"""
function create_train_validation_split(psms::DataFrame, 
                                     validation_ratio::Float64 = 0.2;
                                     seed::Int = 1776)
    Random.seed!(seed)
    
    # Simple stratification by target/decoy only
    target_indices = findall(psms.target)
    decoy_indices = findall(.!psms.target)
    
    # Split targets
    n_target_val = max(1, round(Int, length(target_indices) * validation_ratio))
    target_val = sample(target_indices, n_target_val, replace=false)
    target_train = setdiff(target_indices, target_val)
    
    # Split decoys
    n_decoy_val = max(1, round(Int, length(decoy_indices) * validation_ratio))
    decoy_val = sample(decoy_indices, n_decoy_val, replace=false)
    decoy_train = setdiff(decoy_indices, decoy_val)
    
    train_indices = vcat(target_train, decoy_train)
    validation_indices = vcat(target_val, decoy_val)
    
    # Basic validation warnings
    if length(target_train) < 100
        @user_warn "Few target PSMs in training set: $(length(target_train))"
    end
    if length(decoy_train) < 100
        @user_warn "Few decoy PSMs in training set: $(length(decoy_train))"
    end
    if length(target_val) < 100
        @user_warn "Few target PSMs in validation set: $(length(target_val))"
    end
    if length(decoy_val) < 100
        @user_warn "Few decoy PSMs in validation set: $(length(decoy_val))"
    end
    
    return train_indices, validation_indices
end

#==========================================================
Model Training Framework
==========================================================#

"""
    train_model(model_config::ModelConfig, psms_train::DataFrame, match_between_runs::Bool) -> ModelResult

Trains a single model configuration and returns trained model + metadata.

# Arguments
- `model_config`: Configuration for the model to train
- `psms_train`: Training DataFrame with PSMs
- `match_between_runs`: Whether match between runs is enabled

# Returns
- `ModelResult`: Trained model with metadata
"""
function train_model(model_config::ModelConfig, 
                    psms_train::DataFrame,
                    match_between_runs::Bool) 
    
    start_time = time()
    
    if model_config.model_type == :xgboost
        result = train_xgboost_model(model_config, psms_train, match_between_runs)
    elseif model_config.model_type == :probit
        result = train_probit_model(model_config, psms_train, match_between_runs)
    else
        error("Unknown model type: $(model_config.model_type)")
    end
    
    training_time = time() - start_time
    
    return ModelResult(
        result.model,
        model_config,
        training_time,
        length(result.trained_features),
        result.trained_features
    )
end

"""
    train_xgboost_model(config::ModelConfig, psms_train::DataFrame, match_between_runs::Bool) -> NamedTuple

Wrapper around existing sort_of_percolator_in_memory! for model comparison.

# Arguments
- `config`: Model configuration with hyperparameters
- `psms_train`: Training PSMs DataFrame
- `match_between_runs`: Whether MBR is enabled

# Returns
- NamedTuple with model and trained_features
"""
function train_xgboost_model(config::ModelConfig, 
                           psms_train::DataFrame, 
                           match_between_runs::Bool)
    
    # Filter features to those available and specified
    available_features = filter(f -> f in propertynames(psms_train), config.features)
    if match_between_runs
        mbr_features = [f for f in propertynames(psms_train) if startswith(String(f), "MBR_")]
        append!(available_features, mbr_features)
    end
    
    # Create a copy to avoid modifying original during training
    psms_copy = copy(psms_train)
    
    # Add required columns that may be missing
    if !("accession_numbers" in propertynames(psms_copy))
        psms_copy[!, :accession_numbers] = ["dummy_accession" for _ in 1:nrow(psms_copy)]
    end
    if !("q_value" in propertynames(psms_copy))
        psms_copy[!, :q_value] = zeros(Float32, nrow(psms_copy))
    end
    if !("decoy" in propertynames(psms_copy))
        psms_copy[!, :decoy] = .!psms_copy[!, :target]
    end
    
    # Create temporary file paths (needed by existing function)
    temp_file_paths = ["temp_training_file_$(i).arrow" for i in 1:3]  # Dummy paths for training
    
    # Call existing training function with config hyperparameters
    models = sort_of_percolator_in_memory!(
        psms_copy,
        temp_file_paths,
        available_features,
        match_between_runs;
        max_q_value_xgboost_rescore = Float32(0.01),
        max_q_value_xgboost_mbr_rescore = Float32(0.20),
        min_PEP_neg_threshold_xgboost_rescore = Float32(0.90),
        colsample_bytree = Float64(config.hyperparams[:colsample_bytree]),
        colsample_bynode = Float64(config.hyperparams[:colsample_bynode]),
        min_child_weight = Int(config.hyperparams[:min_child_weight]),
        gamma = Float64(config.hyperparams[:gamma]),
        subsample = Float64(config.hyperparams[:subsample]),
        max_depth = Int(config.hyperparams[:max_depth]),
        eta = Float64(config.hyperparams[:eta]),
        iter_scheme = Vector{Int}(config.hyperparams[:iter_scheme]),
        print_importance = false
    )
    
    return (model = models, trained_features = available_features)
end

"""
    train_probit_model(config::ModelConfig, psms_train::DataFrame, match_between_runs::Bool) -> NamedTuple

Wrapper around existing probit_regression_scoring_cv! for model comparison.

# Arguments
- `config`: Model configuration with hyperparameters
- `psms_train`: Training PSMs DataFrame
- `match_between_runs`: Whether MBR is enabled

# Returns
- NamedTuple with model (nothing for probit) and trained_features
"""
function train_probit_model(config::ModelConfig, 
                          psms_train::DataFrame,
                          match_between_runs::Bool)
    
    # Create copy to avoid modifying original
    psms_copy = copy(psms_train)
    
    # Filter features
    available_features = filter(f -> f in propertynames(psms_copy), config.features)
    if match_between_runs
        mbr_features = [f for f in propertynames(psms_copy) if startswith(String(f), "MBR_")]
        append!(available_features, mbr_features)
    end
    
    # Create temporary file paths
    temp_file_paths = ["temp_training_file_$(i).arrow" for i in 1:3]
    
    # Call existing probit function
    probit_regression_scoring_cv!(
        psms_copy,
        temp_file_paths, 
        available_features,
        match_between_runs;
        n_folds = Int(config.hyperparams[:n_folds])
    )
    
    # Probit doesn't return models - extract coefficients if needed
    return (model = nothing, trained_features = available_features)
end

#==========================================================
Model Evaluation Framework
==========================================================#

"""
    evaluate_model_performance(model_result::ModelResult, psms_val::DataFrame, match_between_runs::Bool, qvalue_threshold::Float64 = 0.01) -> ModelPerformance

Evaluates trained model on validation set and computes performance metrics.
Primary metric: number of targets passing q-value threshold.

# Arguments
- `model_result`: Trained model result
- `psms_val`: Validation DataFrame
- `match_between_runs`: Whether MBR is enabled
- `qvalue_threshold`: Q-value threshold for counting targets (default: 0.01)

# Returns
- `ModelPerformance`: Performance metrics for the model
"""
function evaluate_model_performance(model_result::ModelResult,
                                  psms_val::DataFrame,
                                  match_between_runs::Bool,
                                  qvalue_threshold::Float64 = 0.01)
    
    # Generate predictions on validation set
    if model_result.model_config.model_type == :xgboost
        val_predictions = predict_xgboost_validation(model_result, psms_val)
    elseif model_result.model_config.model_type == :probit
        val_predictions = predict_probit_validation(model_result, psms_val, match_between_runs)
    end
    
    # Compute q-values for validation predictions
    val_qvalues = zeros(Float32, length(val_predictions))
    get_qvalues!(val_predictions, psms_val.target, val_qvalues)
    
    # Primary selection metric: number of targets passing q-value threshold
    n_targets_passing = sum((psms_val.target) .& (val_qvalues .<= qvalue_threshold))
    
    # Secondary metrics for logging
    auc = compute_auc(val_predictions, psms_val.target)
    accuracy = compute_accuracy(val_predictions, psms_val.target)
    sensitivity, specificity = compute_sensitivity_specificity(val_predictions, psms_val.target)
    
    return ModelPerformance(
        model_result.model_config.name,
        n_targets_passing,  # Primary selection metric
        auc,
        accuracy, 
        sensitivity,
        specificity,
        model_result.training_time,
        model_result.n_features
    )
end

"""
    predict_xgboost_validation(model_result::ModelResult, psms_val::DataFrame) -> Vector{Float32}

Generates XGBoost predictions on validation set using CV models.

# Arguments
- `model_result`: Trained XGBoost model result
- `psms_val`: Validation DataFrame

# Returns
- Vector of prediction scores
"""
function predict_xgboost_validation(model_result::ModelResult, 
                                  psms_val::DataFrame)
    
    # Extract CV models from result
    cv_models = model_result.model
    features = model_result.trained_features
    
    # Predict using CV fold assignment
    predictions = zeros(Float32, nrow(psms_val))
    for (fold_idx, model) in cv_models
        fold_mask = psms_val.cv_fold .== fold_idx
        if any(fold_mask)
            # Use final iteration model for prediction
            X_val = Matrix{Float32}(psms_val[fold_mask, features])
            predictions[fold_mask] = predict(model[end], X_val)
        end
    end
    
    return predictions
end

"""
    predict_probit_validation(model_result::ModelResult, psms_val::DataFrame, match_between_runs::Bool) -> Vector{Float32}

For probit, extract predictions from probability column of validation PSMs.

# Arguments
- `model_result`: Trained probit model result
- `psms_val`: Validation DataFrame
- `match_between_runs`: Whether MBR is enabled

# Returns
- Vector of prediction scores
"""
function predict_probit_validation(model_result::ModelResult,
                                 psms_val::DataFrame, 
                                 match_between_runs::Bool)
    
    # For probit, we need to re-apply the trained model
    # Since probit doesn't store coefficients, use the probabilities from the trained PSMs
    if "prob" in propertynames(psms_val)
        return Vector{Float32}(psms_val.prob)
    else
        error("Probit validation prediction requires 'prob' column in validation data")
    end
end

#==========================================================
Performance Metrics
==========================================================#

"""
    compute_auc(predictions::Vector{Float32}, targets::Vector{Bool}) -> Float64

Computes Area Under the ROC Curve.

# Arguments
- `predictions`: Model prediction scores
- `targets`: True labels (true = target, false = decoy)

# Returns
- AUC value between 0.0 and 1.0
"""
function compute_auc(predictions::Vector{Float32}, targets::Vector{Bool})
    n_pos = sum(targets)
    n_neg = length(targets) - n_pos
    
    if n_pos == 0 || n_neg == 0
        return 0.5  # No discrimination possible
    end
    
    # Sort by prediction score (descending)
    sorted_indices = sortperm(predictions, rev=true)
    sorted_targets = targets[sorted_indices]
    
    # Calculate AUC using trapezoidal rule
    tp = 0
    fp = 0
    auc = 0.0
    
    for i in 1:length(sorted_targets)
        if sorted_targets[i]  # True positive
            tp += 1
        else  # False positive
            fp += 1
            # Add area of rectangle
            auc += tp
        end
    end
    
    # Normalize by total area
    return auc / (n_pos * n_neg)
end

"""
    compute_accuracy(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0) -> Float64

Computes classification accuracy.

# Arguments
- `predictions`: Model prediction scores
- `targets`: True labels
- `threshold`: Classification threshold (default: 0.5)

# Returns
- Accuracy between 0.0 and 1.0
"""
function compute_accuracy(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0)
    predicted_labels = predictions .> threshold
    return mean(predicted_labels .== targets)
end

"""
    compute_sensitivity_specificity(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0) -> Tuple{Float64, Float64}

Computes sensitivity (true positive rate) and specificity (true negative rate).

# Arguments
- `predictions`: Model prediction scores
- `targets`: True labels
- `threshold`: Classification threshold (default: 0.5)

# Returns
- Tuple of (sensitivity, specificity)
"""
function compute_sensitivity_specificity(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0)
    predicted_labels = predictions .> threshold
    
    tp = sum(predicted_labels .& targets)
    fn = sum(.!predicted_labels .& targets) 
    tn = sum(.!predicted_labels .& .!targets)
    fp = sum(predicted_labels .& .!targets)
    
    sensitivity = tp > 0 ? tp / (tp + fn) : 0.0  # True positive rate
    specificity = tn > 0 ? tn / (tn + fp) : 0.0  # True negative rate
    
    return sensitivity, specificity
end