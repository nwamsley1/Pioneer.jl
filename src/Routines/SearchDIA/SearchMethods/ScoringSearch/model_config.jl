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
    :prec_mz_qbin,     # Quantile-binned version of :prec_mz
    :sequence_length,
    :charge,
    :irt_pred_qbin,    # Quantile-binned version of :irt_pred
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
    :weight_qbin,      # Quantile-binned version of :weight
    :log2_intensity_explained,
    :tic_qbin,         # Quantile-binned version of :tic
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
    :missed_cleavage, :Mox, :prec_mz_qbin, :sequence_length, :charge,
    # RT features
    :irt_pred_qbin, :irt_error, :irt_diff,
    :ms1_ms2_rt_diff,  # MS1-MS2 RT difference in iRT space
    #:ms1_irt_diff,
    # Spectral features
    :max_y_ions, :y_ions_sum, :longest_y, :y_count, :b_count, :isotope_count,
    :total_ions, :best_rank, :best_rank_iso, :topn, :topn_iso, :gof,
    # Quality metrics
    :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_residual, :max_unmatched_residual, :max_gof,
    :fitted_spectral_contrast, :spectral_contrast, :max_matched_ratio,
    :err_norm, :poisson, :weight_qbin, :log2_intensity_explained, :tic_qbin, :num_scans,
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
    create_model_configurations(ms1_scoring::Bool = true) -> Vector{ModelConfig}

Creates the model configurations for comparison.

# Arguments
- `ms1_scoring`: Whether MS1 scoring is enabled (default: true)

# Returns
- Vector of ModelConfig objects for SimpleLightGBM, AdvancedLightGBM, ProbitRegression,
- ProbitRegressionSimple, and SuperSimplified models
"""
function create_model_configurations(ms1_scoring::Bool = true)
    # Apply MS1 filtering to feature sets if needed
    reduced_features = copy(REDUCED_FEATURE_SET)
    apply_ms1_filtering!(reduced_features, ms1_scoring)

    advanced_features = copy(ADVANCED_FEATURE_SET)
    apply_ms1_filtering!(advanced_features, ms1_scoring)

    minimal_features = copy(MINIMAL_FEATURE_SET)
    apply_ms1_filtering!(minimal_features, ms1_scoring)

    return [
        # Model 1: Simple LightGBM (Default for small datasets)
        ModelConfig(
            "SimpleLightGBM",
            :lightgbm,
            reduced_features,
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
            advanced_features,
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
            vcat(reduced_features, [:intercept]),
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
            minimal_features,
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

"""
    add_quantile_binned_features!(df::DataFrame, features::Vector{Symbol}, n_bins::Int=100)

Creates quantile-binned versions of specified features.
Adds new columns with suffix `_qbin` containing bin indices (1 to n_bins).

This transformation discretizes continuous features into quantile bins, which can:
- Reduce overfitting by limiting model complexity
- Make models more robust to outliers
- Improve generalization on small datasets
- Reduce information leakage from high-precision numeric features

# Arguments
- `df`: DataFrame containing features to bin
- `features`: Vector of feature column names to quantize
- `n_bins`: Number of quantile bins (default: 100). Must be ≤ 255 for UInt8 storage.

# Behavior
- Missing values are preserved as missing in binned columns
- Bins are numbered 1 to n_bins (not 0-indexed)
- Uses UInt8 storage for n_bins ≤ 255, UInt16 otherwise
- Bin edges are computed on non-missing data only

# Example
```julia
add_quantile_binned_features!(psms, [:prec_mz, :irt_pred, :weight, :tic], 100)
# Creates: :prec_mz_qbin, :irt_pred_qbin, :weight_qbin, :tic_qbin
```

# Notes
- For n_bins=100, each bin contains ~1% of the data
- Quantile binning adapts to the actual data distribution
- More robust than fixed-width binning for skewed distributions
"""
function add_quantile_binned_features!(df::DataFrame, features::Vector{Symbol}, n_bins::Int=100)
    if n_bins > 65535
        error("n_bins must be ≤ 65535 (UInt16 limit)")
    end

    for feature in features
        if !hasproperty(df, feature)
            @user_warn "Feature $feature not found in DataFrame, skipping quantile binning"
            continue
        end

        col_data = df[!, feature]

        # Handle missing values
        non_missing_mask = .!ismissing.(col_data)
        if !any(non_missing_mask)
            @user_warn "Feature $feature has all missing values, skipping quantile binning"
            continue
        end

        # Calculate quantiles on non-missing data
        non_missing_data = col_data[non_missing_mask]
        quantiles = range(0, 1, length=n_bins+1)

        bin_edges = try
            StatsBase.quantile(non_missing_data, quantiles)
        catch e
            @user_warn "Failed to compute quantiles for $feature: $e"
            continue
        end

        # Create binned column (UInt8 for n_bins ≤ 255, UInt16 otherwise)
        bin_type = n_bins <= 255 ? UInt8 : UInt16
        binned_col = Vector{Union{Missing, bin_type}}(missing, length(col_data))

        # Assign bin indices (1 to n_bins)
        for i in eachindex(col_data)
            if non_missing_mask[i]
                val = col_data[i]
                # searchsortedfirst gives index in bin_edges, subtract 1 for bin index
                bin_idx = searchsortedfirst(bin_edges, val) - 1
                # Clamp to valid range [1, n_bins]
                bin_idx = max(bin_idx, 1)
                bin_idx = min(bin_idx, n_bins)
                binned_col[i] = bin_type(bin_idx)
            end
        end

        # Add binned column with _qbin suffix
        binned_feature_name = Symbol(string(feature) * "_qbin")
        df[!, binned_feature_name] = binned_col

        @user_info "Created quantile-binned feature $binned_feature_name with $n_bins bins"
    end

    return nothing
end
