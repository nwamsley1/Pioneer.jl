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
using ..Pioneer: @user_warn, @user_info

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