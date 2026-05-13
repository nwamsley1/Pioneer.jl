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
    # Batch A/B/C/D peak-shape features — DISABLED. A/B test 2026-05-12 showed
    # −667 IDs when stacked on top of Batch E (heavy overlap).
    # :relative_position, :dist_to_relative_center,
    # :weight_apex_relative_pos, :weight_chrom_skewness, :weight_chrom_kurtosis,
    # :apex_to_edge_weight_ratio, :n_above_hm_right_of_apex,
    # :weight_chrom_gaussian_r2, :weight_chrom_gaussian_sigma, :weight_chrom_gaussian_apex_irt,
    # :shape_m2, :shape_m1, :shape_0, :shape_p1, :shape_p2,
    # :gof_minus_max_gof_precursor,

    # MS1 point-lookup + chromatogram features.
    :ms1_m0_intensity, :ms1_m1_intensity, :ms1_m0_mass_err_ppm,
    :ms1_corr_weight_m0, :ms1_corr_m0_m1, :ms1_corr_weight_m1,
    :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,
    # Isotope envelope features (2026-05-11)
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred, :ms1_envelope_dev_log2,
    # Sequence composition
    :his_count, :pro_count, :lys_count, :arg_count,
    # Cross-precursor competition
    :weight_ratio_to_2nd_best, :weight_ratio_to_3rd_best, :n_competitors_50pct,

    # Per-rank MS2-fragment intensities + per-precursor chromatogram features.
    :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int,
    :frag_corr_top1_top2, :frag_corr_top1_top3, :frag_corr_top1_weight,
    :frag_corr_mean_pairwise, :frag_corr_min_pairwise, :frag_corr_top3_weight,
    :frag_corr_top5_weight,
    :frag_apex_dispersion_irt, :n_correlated_fragments,
    :n_correlated_fragments_50, :n_correlated_fragments_90,
    :frag_corr_best_weight, :frag_corr_best_m0,
    # M+1 fragment intensities + per-fragment M0/M+1 correlations (2026-05-12)
    :frag1_int_m1, :frag2_int_m1, :frag3_int_m1, :frag4_int_m1, :frag5_int_m1, :frag6_int_m1,
    :frag_corr_m0_m1_top1, :frag_corr_m0_m1_top2, :frag_corr_m0_m1_top3,
    :frag_corr_m0_m1_mean,
    :n_correlated_m0_m1_fragments, :n_correlated_m0_m1_fragments_50,
    # Batch E isolated test — E7 + E1/E2 (cosine confirmed non-redundant via A/B drop test)
    :top3_ms2_mass_error_mean,
    :pred_obs_max_spectral_contrast, :pred_obs_max_scribe,
    :pred_obs_area_spectral_contrast, :pred_obs_area_scribe,
    # E10 dropped (A/B test 2026-05-12: −38 PGs, redundant with frag_corr_best_*)
    # :frag_corr_to_median_mean,
    # E14 (delta median-apex vs scan center) + E6 M0 (log b/y ratio)
    :delta_frame_peak_center, :log_by_ratio_m0,
    # E3 (matched_ratio) DROPPED — A/B test 2026-05-12: −104 PGs vs E7+E1+E2+E14+E6 M0
    # E6 M01 DROPPED — redundant with log_by_ratio_m0
    # :matched_ratio,
    # :log_by_ratio_m01,

    :max_weight,                      # populated by select_best_per_precursor!
    :fitted_hellinger,
    :n_scans,
    # MBR Batch F: only the _true variants appear here (Pass-1 LGBM auto-skips
    # them since they're computed after Pass 1). The _false variants live only
    # in the FTR controller's doubled training frame.
    :MBR_max_pair_prob_true,
    :MBR_log2_weight_ratio_true,
    :MBR_log2_explained_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_is_missing_true,
    # 2026-05-12 evening: consensus aggregates + b/y diff.
    :MBR_top_n_median_score_true,
    :MBR_top_n_irt_diff_true,
    :MBR_log_by_diff_true,
    # 2026-05-12 night: fragment-pattern similarity (RV-coefficient analog).
    :MBR_frag_pattern_cosine_true,
    :MBR_frag_pattern_scribe_true,
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

"""
    apply_feature_blacklist!(features::Vector{Symbol})

Remove features named in env var `PIONEER_FEATURE_BLACKLIST` (comma-separated
list of symbols, with or without leading `:`). Used by the overnight ablation
campaign to test feature minimization without recompiling.
"""
function apply_feature_blacklist!(features::Vector{Symbol})
    bl_str = get(ENV, "PIONEER_FEATURE_BLACKLIST", "")
    isempty(bl_str) && return features
    blset = Set{Symbol}()
    for raw in split(bl_str, ',')
        s = strip(raw)
        isempty(s) && continue
        s = startswith(s, ":") ? s[2:end] : s
        push!(blset, Symbol(s))
    end
    n_before = length(features)
    filter!(f -> !(f in blset), features)
    n_after = length(features)
    if n_before != n_after
        @user_info "  PIONEER_FEATURE_BLACKLIST removed $(n_before - n_after) features; $n_after remain"
    end
    return features
end
