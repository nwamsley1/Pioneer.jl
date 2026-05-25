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
# Conceptually ADVANCED_FEATURE_SET == PRESCORE_FEATURES at training
# time (modulo a few per-precursor aggregate features like :smoothness
# and :irt_fwhm that exist on best-per-precursor rows but not per-scan
# rows). Lists kept separate due to load-order constraint in
# importScripts.jl — PrecursorScoringSearch loads before MainSearch.
const ADVANCED_FEATURE_SET = [
    # Kept in sync manually with PRESCORE_FEATURES (MainSearch/features.jl).
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :charge, :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,
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
    # 2026-05-21: removed 11 features added during the 2026-05-20 experiment
    # (5 max_*, 3 Tier-2 re-adds, 3 min_*). 8-file Olsen showed only ~+1% ID
    # gain (commit 8f2a1583) — not worth the compute + disk cost. Reverted
    # along with their computation in select_best_per_precursor!.
]

