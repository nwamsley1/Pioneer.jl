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
    #:irt_pred_qbin,    # Quantile-binned version of :irt_pred
    :irt_pred,
    :irt_error,
    :irt_diff,
    :longest_y,
    :y_count,
    # `isotope_count` was redundant with `total_ions_iso` (the new name in
    # MainSearchScoredPSM) — removed 2026-05-11.
    :total_ions,
    # best_rank / best_rank_iso / topn / topn_iso re-added to MainSearchScoredPSM
    # via Score! on 2026-05-11 so they now propagate to second_pass_psms.
    :best_rank,
    :best_rank_iso,
    :topn,
    :topn_iso,
    :gof,
    :max_fitted_manhattan_distance,   # populated by select_best_per_precursor!
    :max_matched_residual,
    :max_unmatched_residual,
    :max_gof,                         # populated by select_best_per_precursor!
    # Removed 2026-05-11 — never populated by the current pipeline:
    # :max_y_ions, :y_ions_sum, :max_fitted_spectral_contrast,
    # :fitted_spectral_contrast, :spectral_contrast, :max_matched_ratio,
    # :tic, :scribe, :max_scribe. Add back when the corresponding feature-
    # computation code is wired up.
    :err_norm,
    :poisson,
    #:weight_qbin,      # Quantile-binned version of :weight
    :weight,
    :log2_intensity_explained,
    :smoothness,
    :num_scans,
    :irt_fwhm,
    :rt_fwhm,
    :n_above_hm,

    # ============================================================================
    # PRESCORE-shared features. Every feature in PRESCORE_FEATURES (the per-file
    # MainSearch LightGBM input) is also listed here so the experiment-wide
    # LightGBM trains on the same signal. Previous *_ms1 family
    # (gof_ms1/error_ms1/m0_error_ms1/etc.) was dead — no code populates those
    # columns and they were silently dropped by hasproperty filtering. Replaced
    # 2026-05-11 with the live set; cross-checked against the second_pass_psms
    # Arrow schema.
    # ============================================================================
    :fitted_manhattan_distance,
    :total_ions_iso,
    :spectrum_peak_count,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :best_gof_3scan, :best_manhattan_3scan, :best_max_residual_3scan,
    :irt_dist_best_gof_3scan, :irt_dist_best_manhattan_3scan, :irt_dist_best_max_residual_3scan,
    :worst_max_residual_3scan, :worst_manhattan_3scan,
    :best_gof_5scan, :best_manhattan_5scan, :best_max_residual_5scan,
    :irt_dist_best_gof_5scan, :irt_dist_best_manhattan_5scan, :irt_dist_best_max_residual_5scan,
    :worst_max_residual_11scan, :worst_manhattan_11scan,
    :irt_dist_to_weight_apex,

    # MS1 point-lookup + chromatogram features.
    :ms1_m0_intensity, :ms1_m1_intensity, :ms1_m0_mass_err_ppm,
    :ms1_corr_weight_m0, :ms1_corr_m0_m1,
    :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,

    # Per-rank MS2-fragment intensities + per-precursor chromatogram features.
    :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int,
    :frag_corr_top1_top2, :frag_corr_top1_top3, :frag_corr_top1_weight,
    :frag_corr_mean_pairwise, :frag_corr_min_pairwise, :frag_corr_top3_weight,
    :frag_apex_dispersion_irt, :n_correlated_fragments,

    :max_weight,                      # populated by select_best_per_precursor!
    :fitted_hellinger,
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
        # Live MS1 features (populated by add_ms1_features! and the per-precursor
        # MS1 chromatogram pass). Updated 2026-05-11 to match the live set in
        # ADVANCED_FEATURE_SET (the previous *_ms1 names were dead).
        ms1_features = Set([
            :ms1_m0_intensity, :ms1_m1_intensity, :ms1_m0_mass_err_ppm,
            :ms1_corr_weight_m0, :ms1_corr_m0_m1,
            :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,
        ])
        filter!(f -> !(f in ms1_features), features)
    end
    return features
end
