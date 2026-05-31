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

# Raw full-window co-elution features computed only for the experiment-wide
# PrecursorScoringSearch LightGBM. MainSearch's per-run PRESCORE_FEATURES does
# not consume them.
const WIDE_WINDOW_FEATURES = [
    :wide_ms1_m0_candidate_fraction,
    :wide_frag_candidate_fraction,
    :wide_ms1_frag_sum_corr,
    :wide_frag_corr_mean,
    :wide_n_correlated_fragments,
    :wide_n_correlated_fragments_bitvec_rank,
    :wide_frag_corr_best_m0,
    :wide_signal_support,
    :wide_core_n_scans,
]

# Empirical observed-fragment shape features. These are cross-run only:
# they compare each per-run best PSM's deconvolved observed ("shadow") top-8
# fragment intensities against a top-5 leave-one-out empirical consensus for
# the same precursor, weighted by per-run confidence (1 - PEP).
# empirical_frag_ref_pep is computed as a diagnostic/reference-confidence column
# but is intentionally excluded from ADVANCED_FEATURE_SET to avoid inflating
# false transfer rate through per-run PEP feedback.
const EMPIRICAL_FRAGMENT_FEATURES = [
    :empirical_frag_best_hellinger,
    :empirical_frag_ref_pep,
]

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
# ADVANCED_FEATURE_SET intentionally extends PRESCORE_FEATURES with features
# that only exist, or should only be learned, after per-run filtering (for
# example wide-window and scan-local competition features). Lists are kept
# separate due to load-order constraint in importScripts.jl —
# PrecursorScoringSearch loads before MainSearch.
const ADVANCED_FEATURE_SET = [
    # Core prescore-compatible features.
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :charge, :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger, :spectral_contrast,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,
    :frag5_int, :frag6_int, :frag7_int, :frag8_int,
    # frag_corr_mean_pairwise (Spearman) dropped 2026-05-13.
    :frag_apex_dispersion_irt,
    :n_correlated_fragments,
    :n_correlated_fragments_bitvec_rank,
    :frag_corr_best_m0,
    :ms1_m0_peak_frag_intensity_fraction,
    :ms1_m0_peak_n_precursors,
    :scan_prec_mz_n_precursors,
    :frag_competition_num_unique_fragments,
    :frag_competition_mean_candidates,
    :precursor_fraction_transmitted,
    :n_scans_other_traces,
    :n_frags_detected_union,
    :n_frags_detected_intersection,
    :n_frags_detected_union_bitvec_rank,
    :n_frags_detected_intersection_bitvec_rank,
    :frag_observed_sum_hellinger,
    :empirical_frag_best_hellinger,
    :trace_other_weight_corr,
    :trace_other_frag_sum_corr,
    :trace_other_apex_delta_irt,
    :wide_ms1_m0_candidate_fraction,
    :wide_frag_candidate_fraction,
    :wide_ms1_frag_sum_corr,
    :wide_frag_corr_mean,
    :wide_n_correlated_fragments,
    :wide_n_correlated_fragments_bitvec_rank,
    :wide_frag_corr_best_m0,
    :wide_signal_support,
    :wide_core_n_scans,
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
