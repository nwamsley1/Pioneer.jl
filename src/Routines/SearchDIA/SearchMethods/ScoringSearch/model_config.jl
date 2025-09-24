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
- `name`: Model identifier (e.g., "SimpleLightGBM", "ProbitRegression")
- `model_type`: Algorithm type (:lightgbm or :probit)
- `features`: Vector of feature symbols to use
- `hyperparams`: Dictionary of hyperparameters for the model
"""
struct ModelConfig
    name::String
    model_type::Symbol  # :lightgbm or :probit
    features::Vector{Symbol}
    hyperparams::Dict{Symbol, Any}
end

# Feature set definitions

# Full feature set used for advanced LightGBM model (matches out-of-memory case)
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
    :ms1_ms2_rt_diff,  # MS1-MS2 RT difference in iRT space
    #:ms1_irt_diff,
    #:weight_ms1,
    
    :gof_ms1,
    :max_matched_residual_ms1,
    :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1,
    :error_ms1,
    :m0_error_ms1,
    :n_iso_ms1,
    :big_iso_ms1,
    :rt_max_intensity_ms1,
    :rt_diff_max_intensity_ms1,
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
    :irt_pred, :irt_error, :irt_diff,
    :ms1_ms2_rt_diff,  # MS1-MS2 RT difference in iRT space
    #:ms1_irt_diff,
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
    :weight_ms1, 
    :gof_ms1, :max_matched_residual_ms1, :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1, :error_ms1, :m0_error_ms1, :n_iso_ms1,
    :big_iso_ms1, :rt_max_intensity_ms1, 
    :rt_diff_max_intensity_ms1, 
    :ms1_features_missing
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
- Vector of ModelConfig objects for SimpleLightGBM, AdvancedLightGBM, ProbitRegression,
- ProbitRegressionSimple, and SuperSimplified models
"""
function create_model_configurations()
    return [
        # Model 1: Simple LightGBM (Default for small datasets)
        ModelConfig(
            "SimpleLightGBM",
            :lightgbm,
            REDUCED_FEATURE_SET,
            Dict(
                :feature_fraction => 0.8,
                :min_data_in_leaf => 20,
                :min_gain_to_split => 0.1,
                :bagging_fraction => 0.8,
                :max_depth => 4,
                :num_leaves => 15,
                :learning_rate => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        ),
        
        # Model 2: Advanced LightGBM (Same as used for >100K PSMs)
        ModelConfig(
            "AdvancedLightGBM",
            :lightgbm,
            ADVANCED_FEATURE_SET,
            Dict(
                :feature_fraction => 0.5,
                :min_data_in_leaf => 200,
                :min_gain_to_split => 0.0,
                :bagging_fraction => 0.25,
                :max_depth => -1,
                :num_leaves => 63,
                :learning_rate => 0.05,
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

        # Model 5: Super Simplified LightGBM Model
        ModelConfig(
            "SuperSimplified",
            :lightgbm,
            MINIMAL_FEATURE_SET,
            Dict(
                :feature_fraction => 0.8,
                :min_data_in_leaf => 10,
                :min_gain_to_split => 0.1,
                :bagging_fraction => 0.8,
                :max_depth => 4,
                :num_leaves => 15,
                :learning_rate => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        )
    ]
end

"""
    apply_ms1_filtering!(features::Vector{Symbol}, ms1_scoring::Bool)

Remove MS1 features from feature list if MS1 scoring is disabled.
This is applied to model configurations to prevent using MS1 features
when they would have zero variance (ms1_scoring=false).

# Arguments
- `features`: Vector of feature symbols to filter
- `ms1_scoring`: Whether MS1 scoring is enabled

# Returns
- Modified feature vector with MS1 features removed if ms1_scoring=false
"""
function apply_ms1_filtering!(features::Vector{Symbol}, ms1_scoring::Bool)
    if !ms1_scoring
        # MS1 features to exclude when ms1_scoring=false
        ms1_features = Set([
            :ms1_irt_diff, :ms1_ms2_rt_diff, :weight_ms1, :gof_ms1, :max_matched_residual_ms1,
            :max_unmatched_residual_ms1, :fitted_spectral_contrast_ms1, :error_ms1,
            :m0_error_ms1, :n_iso_ms1, :big_iso_ms1, :rt_max_intensity_ms1,
            :rt_diff_max_intensity_ms1, :ms1_features_missing
        ])

        original_count = length(features)
        filter!(feature -> !(feature in ms1_features), features)
        removed_count = original_count - length(features)

        if removed_count > 0
            @user_info "Filtered $removed_count MS1 features from model (ms1_scoring=false)"
        end
    end

    return features
end

"""
    create_filtered_model_configurations(ms1_scoring::Bool = true) -> Vector{ModelConfig}

Create model configurations with optional MS1 feature filtering.
When ms1_scoring=false, removes MS1 features from all feature sets to prevent
zero-variance columns and singular matrices in ML training.

# Arguments
- `ms1_scoring`: Whether MS1 scoring is enabled

# Returns
- Vector of ModelConfig objects with appropriately filtered features
"""
function create_filtered_model_configurations(ms1_scoring::Bool = true)
    configs = create_model_configurations()

    # Apply MS1 filtering to all model configurations
    for config in configs
        apply_ms1_filtering!(config.features, ms1_scoring)
    end

    return configs
end
