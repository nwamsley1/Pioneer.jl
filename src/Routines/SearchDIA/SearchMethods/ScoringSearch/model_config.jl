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
Model Configuration for ScoringSearch
This file contains only the minimal functionality needed for model selection.
The full model_comparison.jl file is deprecated and should not be used.
==========================================================#

"""
Configuration for a single model in the scoring framework.

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

# Feature set definitions

# Full feature set used for advanced XGBoost model (matches out-of-memory case)
# Note: :target is excluded as it's the label, not a feature
const ADVANCED_FEATURE_SET = [
    :missed_cleavage,
    :Mox,
    :prec_mz,
    :sequence_length,
    :charge,
    :irt_pred,
    :irt_error,
    :irt_diff,
    :max_y_ions,
    :y_ions_sum,
    :longest_y,
    :y_count,
    :b_count,
    :isotope_count,
    :total_ions,
    :best_rank,
    :best_rank_iso,
    :topn,
    :topn_iso,
    :gof,
    :max_fitted_manhattan_distance,
    :max_fitted_spectral_contrast,
    :max_matched_residual,
    :max_unmatched_residual,
    :max_gof,
    :fitted_spectral_contrast,
    :spectral_contrast,
    :max_matched_ratio,
    :err_norm,
    :poisson,
    :weight,
    :log2_intensity_explained,
    :tic,
    :num_scans,
    :smoothness,
    :rt_diff,
    :ms1_irt_diff,
    :weight_ms1,
    :gof_ms1,
    :max_matched_residual_ms1,
    :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1,
    :error_ms1,
    :m0_error_ms1,
    :n_iso_ms1,
    :big_iso_ms1,
    :ms1_features_missing,
    :percent_theoretical_ignored,
    :scribe,
    :max_scribe
    # MBR features added automatically if match_between_runs=true
]

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

#=
const PROBIT_FEATURE_SET = [
    # Core peptide properties
    :missed_cleavage, :Mox, :prec_mz, :sequence_length, :charge,
    # RT features
    :irt_error,
    # Spectral features
    :y_count, :b_count, :isotope_count,
    :total_ions, :best_rank, :topn, :gof,
    # Quality metrics
    :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_residual, :max_unmatched_residual, :max_gof,
    :err_norm, :poisson, :weight, :log2_intensity_explained, :tic, :num_scans,
    :smoothness, :percent_theoretical_ignored, :scribe, :max_scribe,
    # MS1 features
    :weight_ms1, :gof_ms1, :error_ms1, :ms1_features_missing
    # MBR features added automatically if match_between_runs=true
]
=#


const MINIMAL_FEATURE_SET = [
    :fitted_spectral_contrast,
    :max_matched_residual,
    :max_unmatched_residual,
    :err_norm,
    :log2_intensity_explained
]

"""
    create_model_configurations() -> Vector{ModelConfig}

Creates the model configurations for comparison.

# Returns
- Vector of ModelConfig objects for SimpleXGBoost, AdvancedXGBoost, ProbitRegression,
- ProbitRegressionSimple, and SuperSimplified models
"""
function create_model_configurations()
    return [
        # Model 1: Simple XGBoost (Default for small datasets)
        ModelConfig(
            "SimpleXGBoost",
            :xgboost,
            REDUCED_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.8,
                :min_child_weight => 20,
                :gamma => 0.1,
                :subsample => 0.8,
                :max_depth => 4,
                :eta => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        ),
        
        # Model 2: Advanced XGBoost (Same as used for >100K PSMs)
        ModelConfig(
            "AdvancedXGBoost",
            :xgboost,
            ADVANCED_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.5,
                :min_child_weight => 5,
                :gamma => 1.0,
                :subsample => 0.25,
                :max_depth => 10,
                :eta => 0.05,
                :iter_scheme => [100, 200, 200]
            )
        ),
        
        # Model 3: Probit Regression
        ModelConfig(
            "ProbitRegression",
            :probit,
            vcat(REDUCED_FEATURE_SET, [:intercept]),
            Dict(
                :max_iter => 30
            )
        ),
        
        # Model 4: Probit Regression (Simplified feature set)
        ModelConfig(
            "ProbitRegressionSimple",
            :probit,
            vcat([
                # Keep core, broadly available features
                :spectral_contrast,
                :entropy_score,
                :scribe,
                :irt_error,
                :err_norm,
                :y_count,
                :tic
            ], [:intercept]),
            Dict(
                :max_iter => 30
            )
        ),

        # Model 5: Super Simplified Model
        ModelConfig(
            "SuperSimplified",
            :xgboost,
            MINIMAL_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.8,
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
