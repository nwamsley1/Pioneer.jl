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
Feature set for the experiment-wide LightGBM classifier in
PrecursorScoringSearch. Hyperparameters and the low-data probit fallback
live in MainSearch/scoring.jl (`SHARED_LGBM_HP`,
`train_psm_classifier_with_fallback`) — both stages share that single
classifier definition. This file is just the feature list + an MS1 filter.
==========================================================#

# Full feature set used for the experiment-wide LightGBM classifier.
# Note: `:target` is excluded as it's the label, not a feature. The `_qbin`
# variants are computed by `add_quantile_binned_features!` before training and
# are kept commented out — the un-binned versions are currently used.
const ADVANCED_FEATURE_SET = [
    :missed_cleavage,
    :Mox,
    #:prec_mz_qbin,     # Quantile-binned version of :prec_mz
    :prec_mz,
    :sequence_length,
    :charge,
    #:irt_pred_qbin,    # Quantile-binned version of :irt_pred
    :irt_pred,
    :irt_error,
    :irt_diff,
    :max_y_ions,
    :y_ions_sum,
    :longest_y,
    :y_count,
    :b_count,
    :total_ions,
    :total_ions_iso,
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
    #:weight_qbin,      # Quantile-binned version of :weight
    :weight,
    :log2_intensity_explained,
    #:tic_qbin,         # Quantile-binned version of :tic
    :tic,
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
    :max_scribe,
    :max_weight,
    :fitted_hellinger,
    :n_scans,
    # MBR Phase 2 features (auto-skipped on first pass when columns are absent;
    # used in the second pass after compute_mbr_features! populates them).
    :MBR_max_pair_prob,
    :MBR_log2_weight_ratio,
    :MBR_log2_explained_ratio,
    :MBR_best_irt_diff,
    :MBR_is_missing,
]

"""
    apply_ms1_filtering!(features::Vector{Symbol}, ms1_scoring::Bool)

Remove MS1 features from `features` in-place when `ms1_scoring=false`.
This prevents using MS1 features when they would have zero variance.

# Returns
The mutated feature vector (for chaining).
"""
function apply_ms1_filtering!(features::Vector{Symbol}, ms1_scoring::Bool)
    if !ms1_scoring
        ms1_features = Set([
            :ms1_irt_diff, :ms1_ms2_rt_diff, :weight_ms1, :gof_ms1,
            :max_matched_residual_ms1, :max_unmatched_residual_ms1,
            :fitted_spectral_contrast_ms1, :error_ms1, :m0_error_ms1,
            :n_iso_ms1, :big_iso_ms1, :rt_max_intensity_ms1,
            :rt_diff_max_intensity_ms1, :ms1_features_missing,
        ])
        filter!(f -> !(f in ms1_features), features)
    end
    return features
end
