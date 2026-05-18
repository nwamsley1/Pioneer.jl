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

# ADVANCED_FEATURE_SET = PRESCORE_FEATURES (the smart-composite reduced set)
# + 5 MBR features that are computed AFTER MainSearch's Pass 1 prepass
# and are therefore only available to the experiment-wide LightGBM here.
# By construction, the two LGBMs are trained on the same features apart
# from those 5 MBR-specific signals that depend on cross-run pooling.
# Listed as a literal here because Julia precompile struggled with
# `vcat(collect(PRESCORE_FEATURES), …)` evaluated at const time.
const ADVANCED_FEATURE_SET = [
    # ── PRESCORE_FEATURES (50 features, kept in sync manually) ──
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,
    # frag3_int / frag4_int dropped 2026-05-18 (A/B test).
    :frag1_int, :frag2_int,
    # frag_corr_mean_pairwise (Spearman) dropped 2026-05-13.
    :frag_apex_dispersion_irt,
    :n_correlated_fragments,
    :frag_corr_best_m0,
    :top3_ms2_mass_error_mean,
    :delta_frame_peak_center,
    :log_by_ratio_m0,
    :n_scans,
    :prec_mz,
    :irt_fwhm,
    :smoothness,
    :log2_intensity_explained, :longest_y,
    # Tier-2 drop-all-5 (2026-05-13): rt_fwhm, num_scans, irt_pred,
    # best_rank_iso, total_ions_iso removed.
    # Tier-4 drop-top-10 (2026-05-14): max_matched_residual, lys_count,
    # frag_corr_min_pairwise, ms1_corr_weight_m0, irt_dist_to_weight_apex,
    # best_rank, ms1_envelope_dev_log2, n_above_hm, topn, topn_iso removed —
    # 8-file Olsen: −476 IDs (−0.11%), +203 PGs (+0.39%), pf −864 (−0.24%).
    # ── MBR features (only available post-MainSearch) ──
    :MBR_max_pair_prob_true,
    :MBR_log2_weight_ratio_true,
    :MBR_log2_explained_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_is_missing_true,
    # rtv1 (2026-05-13)
    :MBR_best_rt_diff_true,
]

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

