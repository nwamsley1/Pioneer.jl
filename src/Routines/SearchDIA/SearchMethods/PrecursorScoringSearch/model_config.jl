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

# ADVANCED_FEATURE_SET drives the ScoringSearch Pass-1 LGBM
# (in score_psms.jl::_score_precursor_isotope_traces_{mbr,no_mbr}).
#
# Important: the MBR features (MBR_max_pair_prob_true, etc.) are NOT in
# this list. They are computed AFTER Pass-1 LGBM trains (see steps 5-6
# in _score_precursor_isotope_traces_mbr), so they wouldn't exist on the
# DataFrame at Pass-1 training time anyway. MBR features are consumed
# only by the FTR controller (mbr_ftr.jl::FTR_FEATURES_F_TRUE), which is
# a separate LGBM trained later in the pipeline.
#
# Pre-2026-05-19 this list redundantly included the 6 MBR features —
# they were always silently filtered out by the hasproperty guard in
# pass1_oom.jl. Removed to make the design intent explicit.
#
# ADVANCED_FEATURE_SET mostly tracks PRESCORE_FEATURES, but it can use
# best-per-precursor features that MainSearch cannot use for scan selection.
# Lists are kept separate due to load-order constraint in importScripts.jl:
# PrecursorScoringSearch loads before MainSearch.
const SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD = 0.03f0
const SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD = 0.01f0
const SCORING_SEMISUPERVISED_MIN_TARGET_GAIN = 0.01f0
const SCORING_SEMISUPERVISED_MAX_ITERATIONS = 8

const FLANKING_WINDOW_FEATURES = [
    :flanking_ms1_m0_candidate_fraction,
    :flanking_frag_candidate_fraction,
    :flanking_ms1_frag_sum_corr,
    :flanking_frag_corr_mean,
    :flanking_frag_corr_strength,
    :flanking_frag_corr_effective_n,
    :flanking_frag_corr_best_m0,
    :flanking_signal_support,
    :n_contiguous_scans,
    :frag_apex_gt2x_flank_bitvec_rank,
]

const SMOOTHED_FRAGMENT_INTENSITY_COLUMNS = (
    :frag1_smoothed_intensity,
    :frag2_smoothed_intensity,
    :frag3_smoothed_intensity,
    :frag4_smoothed_intensity,
    :frag5_smoothed_intensity,
    :frag6_smoothed_intensity,
    :frag7_smoothed_intensity,
    :frag8_smoothed_intensity,
)

const FITTED_FRAGMENT_INTENSITY_COLUMNS = (
    :fitted_frag1_int, :fitted_frag2_int, :fitted_frag3_int, :fitted_frag4_int,
    :fitted_frag5_int, :fitted_frag6_int, :fitted_frag7_int, :fitted_frag8_int,
)

const SHADOW_FRAGMENT_INTENSITY_COLUMNS = (
    :shadow_frag1_int, :shadow_frag2_int, :shadow_frag3_int, :shadow_frag4_int,
    :shadow_frag5_int, :shadow_frag6_int, :shadow_frag7_int, :shadow_frag8_int,
)

const ADVANCED_FEATURE_SET = [
    # Kept close to PRESCORE_FEATURES (MainSearch/features.jl), with
    # cross-run-only substitutions where the selected best PSM is known.
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,
    :ms1_isotope_dotp_m0_m1_m2,
    :ms1_m0_m1_m2_window_fraction, :ms1_ms2_explained_delta,
    :ms1_m0_m1_m2_window_fraction_pc, :ms1_ms2_explained_delta_pc,
    SMOOTHED_FRAGMENT_INTENSITY_COLUMNS...,
    :smoothed_2d_shadow_hellinger,
    # frag_corr_mean_pairwise (Spearman) dropped 2026-05-13.
    :frag_apex_dispersion_irt,
    :n_correlated_fragments,
    :n_correlated_fragments_bitvec_rank,
    :frag_corr_strength,
    :frag_corr_effective_n,
    :frag_corr_best_m0,
    # Sciex ZT within-metascan (shape) fragment features — present only when the
    # meta-scan collapse ran (PIONEER_ZT_METASCAN_K>0); filtered out by hasproperty
    # on all other datasets.
    :frag_corr_strength_shape,
    :frag_corr_effective_n_shape,
    :frag_corr_best_shape,
    :frag_apex_dispersion_shape,
    :n_correlated_fragments_shape,
    :n_correlated_fragments_bitvec_rank_shape,
    # Sciex ZT across-cycle (elution) fragment/MS1 features + real cycle count.
    :ms1_corr_weight_m0_elution,
    :ms1_corr_m0_m1_elution,
    :ms1_apex_offset_irt_elution,
    :ms1_weight_apex_to_m0_apex_irt_elution,
    :frag_corr_strength_elution,
    :frag_corr_effective_n_elution,
    :frag_corr_best_elution,
    :frag_apex_dispersion_elution,
    :n_correlated_fragments_elution,
    :n_correlated_fragments_bitvec_rank_elution,
    :n_scans_metascan,
    :n_frags_detected_union,
    :n_frags_detected_intersection,
    :n_frags_detected_union_bitvec_rank,
    :n_frags_detected_intersection_bitvec_rank,
    :precursor_fraction_transmitted,
    :n_scans_other_windows,
    :other_window_weight_corr,
    :other_window_apex_delta_irt,
    FLANKING_WINDOW_FEATURES...,
    :ms1_m0_peak_frag_intensity_fraction,
    :ms1_m0_peak_n_precursors,
    :scan_prec_mz_n_precursors,
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
    # 2026-05-21: removed 11 features added during the 2026-05-20 experiment
    # (5 max_*, 3 Tier-2 re-adds, 3 min_*). 8-file Olsen showed only ~+1% ID
    # gain (commit 8f2a1583) — not worth the compute + disk cost. Reverted
    # along with their computation in select_best_per_precursor!.
]
