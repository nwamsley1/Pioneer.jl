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

# We need to access the original function for fallback
# This will be available from the parent module after include
# score_precursor_isotope_traces_in_memory!

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

#==========================================================
Model Selection and Integration
==========================================================#

"""
    select_best_model(performances::Vector{ModelPerformance}) -> String

Selects best model based on number of targets passing q-value threshold.
Tie-breaking: highest AUC, then fastest training time.

# Arguments
- `performances`: Vector of performance results for all models

# Returns
- Name of the best performing model
"""
function select_best_model(performances::Vector{ModelPerformance})
    # Primary metric: number of targets passing q-value threshold
    max_targets = maximum([p.n_targets_passing_qval for p in performances])
    best_candidates = filter(p -> p.n_targets_passing_qval == max_targets, performances)
    
    if length(best_candidates) == 1
        best_model = best_candidates[1]
    else
        # Tie-breaking: highest AUC
        max_auc = maximum([p.validation_auc for p in best_candidates])
        auc_candidates = filter(p -> p.validation_auc == max_auc, best_candidates)
        
        if length(auc_candidates) == 1
            best_model = auc_candidates[1]
        else
            # Final tie-breaking: fastest training time
            best_model = auc_candidates[argmin([p.training_time for p in auc_candidates])]
        end
    end
    
    @user_info "Selected $(best_model.model_name) with $(best_model.n_targets_passing_qval) targets passing q-value threshold"
    
    return best_model.model_name
end

"""
    train_selected_model_full_dataset(best_model_name::String, model_configs::Vector{ModelConfig}, psms_full::DataFrame, file_paths::Vector{String}, match_between_runs::Bool, max_q_value_xgboost_rescore::Float32, max_q_value_xgboost_mbr_rescore::Float32, min_PEP_neg_threshold_xgboost_rescore::Float32) -> Any

Trains the selected best model on the full dataset (100% of PSMs) using existing procedures.

# Arguments
- `best_model_name`: Name of the selected best model
- `model_configs`: Vector of all model configurations
- `psms_full`: Full PSM DataFrame
- `file_paths`: File paths for training
- `match_between_runs`: Whether MBR is enabled
- `max_q_value_xgboost_rescore`: Max q-value for XGBoost rescoring
- `max_q_value_xgboost_mbr_rescore`: Max q-value for MBR rescoring
- `min_PEP_neg_threshold_xgboost_rescore`: Min PEP threshold for negative relabeling

# Returns
- Trained model(s) or nothing for probit
"""
function train_selected_model_full_dataset(best_model_name::String,
                                         model_configs::Vector{ModelConfig},
                                         psms_full::DataFrame,
                                         file_paths::Vector{String},
                                         match_between_runs::Bool,
                                         max_q_value_xgboost_rescore::Float32,
                                         max_q_value_xgboost_mbr_rescore::Float32,
                                         min_PEP_neg_threshold_xgboost_rescore::Float32)
    
    selected_config_idx = findfirst(c -> c.name == best_model_name, model_configs)
    config = model_configs[selected_config_idx]
    
    @user_info "Training selected model $(config.name) on full dataset ($(nrow(psms_full)) PSMs)"
    
    # Filter features to those available and specified
    available_features = filter(f -> f in propertynames(psms_full), config.features)
    if match_between_runs
        mbr_features = [f for f in propertynames(psms_full) if startswith(String(f), "MBR_")]
        append!(available_features, mbr_features)
    end
    
    # Use existing training procedures exactly as implemented
    if config.model_type == :xgboost
        # Call existing sort_of_percolator_in_memory! with selected config
        return sort_of_percolator_in_memory!(
            psms_full,
            file_paths,
            available_features,
            match_between_runs;
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore,
            colsample_bytree = Float64(config.hyperparams[:colsample_bytree]),
            colsample_bynode = Float64(config.hyperparams[:colsample_bynode]),
            min_child_weight = Int(config.hyperparams[:min_child_weight]),
            gamma = Float64(config.hyperparams[:gamma]),
            subsample = Float64(config.hyperparams[:subsample]),
            max_depth = Int(config.hyperparams[:max_depth]),
            eta = Float64(config.hyperparams[:eta]),
            iter_scheme = Vector{Int}(config.hyperparams[:iter_scheme]),
            print_importance = true
        )
    elseif config.model_type == :probit
        # Call existing probit_regression_scoring_cv! procedure
        probit_regression_scoring_cv!(
            psms_full,
            file_paths,
            available_features,
            match_between_runs;
            n_folds = Int(config.hyperparams[:n_folds])
        )
        return nothing  # Probit doesn't return models
    end
end

"""
    log_model_comparison_results(performances::Vector{ModelPerformance})

Logs detailed comparison results for all models.

# Arguments
- `performances`: Vector of performance results for all models
"""
function log_model_comparison_results(performances::Vector{ModelPerformance})
    @user_info "Model Comparison Results:"
    @user_info "========================"
    
    # Sort by primary metric (targets passing q-value threshold)
    sorted_perfs = sort(performances, by=p->p.n_targets_passing_qval, rev=true)
    
    for (i, perf) in enumerate(sorted_perfs)
        @user_info "$(i). $(perf.model_name):"
        @user_info "   Targets Passing Q≤0.01: $(perf.n_targets_passing_qval)"
        @user_info "   AUC: $(round(perf.validation_auc, digits=4))"
        @user_info "   Accuracy: $(round(perf.validation_accuracy, digits=4))"  
        @user_info "   Sensitivity: $(round(perf.validation_sensitivity, digits=4))"
        @user_info "   Specificity: $(round(perf.validation_specificity, digits=4))"
        @user_info "   Training Time: $(round(perf.training_time, digits=2))s"
        @user_info "   Features: $(perf.n_features)"
        @user_info ""
    end
end

"""
    write_model_comparison_report(performances::Vector{ModelPerformance}, output_dir::String)

Writes detailed CSV report of model comparison results.

# Arguments
- `performances`: Vector of performance results
- `output_dir`: Directory to write the report
"""
function write_model_comparison_report(performances::Vector{ModelPerformance}, 
                                     output_dir::String)
    
    report_path = joinpath(output_dir, "model_comparison_report.csv")
    
    df = DataFrame(
        model_name = [p.model_name for p in performances],
        n_targets_passing_qval = [p.n_targets_passing_qval for p in performances],
        validation_auc = [p.validation_auc for p in performances],
        validation_accuracy = [p.validation_accuracy for p in performances],
        validation_sensitivity = [p.validation_sensitivity for p in performances], 
        validation_specificity = [p.validation_specificity for p in performances],
        training_time_seconds = [p.training_time for p in performances],
        n_features = [p.n_features for p in performances]
    )
    
    CSV.write(report_path, df)
    @user_info "Model comparison report written to: $report_path"
end

#==========================================================
Main Entry Point
==========================================================#

"""
    score_precursor_isotope_traces_in_memory_with_comparison!(best_psms::DataFrame, file_paths::Vector{String}, precursors, match_between_runs::Bool, max_q_value_xgboost_rescore::Float32, max_q_value_xgboost_mbr_rescore::Float32, min_PEP_neg_threshold_xgboost_rescore::Float32, enable_model_comparison::Bool = false, validation_split_ratio::Float64 = 0.2, qvalue_threshold::Float64 = 0.01, output_dir::String = ".")

Enhanced version of score_precursor_isotope_traces_in_memory! with model comparison.
Only applies to in-memory approach (<100k PSMs).

# Arguments
- `best_psms`: PSM DataFrame for training
- `file_paths`: File paths for training
- `precursors`: Library precursors (for compatibility)
- `match_between_runs`: Whether MBR is enabled
- `max_q_value_xgboost_rescore`: Max q-value for XGBoost rescoring
- `max_q_value_xgboost_mbr_rescore`: Max q-value for MBR rescoring
- `min_PEP_neg_threshold_xgboost_rescore`: Min PEP threshold
- `enable_model_comparison`: Whether to enable model comparison
- `validation_split_ratio`: Fraction for validation split
- `qvalue_threshold`: Q-value threshold for target counting
- `output_dir`: Output directory for reports

# Returns
- Trained models from best performing approach
"""
function score_precursor_isotope_traces_in_memory_with_comparison!(
    best_psms::DataFrame,
    file_paths::Vector{String}, 
    precursors,  # Keep for compatibility, not used in comparison logic
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    enable_model_comparison::Bool = false,
    validation_split_ratio::Float64 = 0.2,
    qvalue_threshold::Float64 = 0.01,
    output_dir::String = "."
)
    
    n_psms = size(best_psms, 1)
    
    # Check if model comparison should be enabled (in-memory only)
    if !enable_model_comparison || n_psms < 1000 || n_psms >= 100000
        @user_info "Model comparison disabled - falling back to existing logic"
        # Call the original score_precursor_isotope_traces_in_memory! function
        return score_precursor_isotope_traces_in_memory!(
            best_psms,
            file_paths,
            precursors,
            match_between_runs,
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore
        )
    end
    
    @user_info "Starting model comparison with $(n_psms) PSMs (in-memory approach)"
    
    # Phase 1: Create train/validation split
    train_indices, val_indices = create_train_validation_split(best_psms, validation_split_ratio)
    psms_train = best_psms[train_indices, :]
    psms_val = best_psms[val_indices, :]
    
    @user_info "Split: $(nrow(psms_train)) training, $(nrow(psms_val)) validation PSMs"
    
    # Phase 2: Define model configurations  
    model_configs = create_model_configurations()
    
    # Phase 3: Train all models on training set
    model_results = ModelResult[]
    for config in model_configs
        @user_info "Training model: $(config.name)"
        try
            result = train_model(config, psms_train, match_between_runs)
            push!(model_results, result)
        catch e
            @user_warn "Failed to train model $(config.name): $e"
        end
    end
    
    if isempty(model_results)
        @user_warn "No models trained successfully - falling back to existing logic"
        return nothing
    end
    
    # Phase 4: Evaluate all models on validation set
    performances = ModelPerformance[]
    for result in model_results
        @user_info "Evaluating model: $(result.model_config.name)"
        try
            perf = evaluate_model_performance(result, psms_val, match_between_runs, qvalue_threshold)
            push!(performances, perf)
        catch e
            @user_warn "Failed to evaluate model $(result.model_config.name): $e"
        end
    end
    
    if isempty(performances)
        @user_warn "No models evaluated successfully - falling back to existing logic"
        return nothing
    end
    
    # Phase 5: Select best model and log results
    log_model_comparison_results(performances)
    write_model_comparison_report(performances, output_dir)
    best_model_name = select_best_model(performances)
    
    # Phase 6: Train selected model on full dataset using existing procedures
    @user_info "Training $(best_model_name) on full dataset"
    models = train_selected_model_full_dataset(
        best_model_name, model_configs, best_psms, file_paths, 
        match_between_runs, max_q_value_xgboost_rescore,
        max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
    )
    
    return models
end