# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Types for MBRRecoverySearch.
#
# The recovery stage targets B-bucket (Pioneer "PSM never scored") cells for
# donor pids that were identified at q ≤ MBR_RECOVERY_DONOR_QVAL in ≥1 file.
# Each cell is re-extracted via targeted PMM deconvolution in a per-file
# calibrated RT window. The emitted RecoveredSeedRow carries the seven MBR_*
# features inline so PrecursorScoringSearch can vertical-concat it into the
# unified PSM table without further per-row feature compute.

"""
Per-recovery-row payload. One row per (donor pid × receiver file) cell where:
  - pid is in the global donor pool (qval ≤ MBR_RECOVERY_DONOR_QVAL in ≥1 file)
  - pid is NOT in receiver file's main_search_psms
  - PMM weight is finite at the donor-projected scan window

The seven MBR_*_true features mirror mbr_streaming.jl outputs and are
computed inline from (receiver result, best donor entry) at extraction time
— no donor-dict search is needed downstream.
"""
struct RecoveredSeedRow
    # Identity + provenance
    precursor_idx::UInt32
    target::Bool
    scan_idx::UInt32                          # Max-weight scan inside the window
    cv_fold::UInt8                            # From library; preserves FTR fold-split

    # Receiver-side raw features (mirror of _MBRDonorEntry's payload schema)
    weight::Float32
    log2_intensity_explained::Float32
    irt_obs::Float32                          # Receiver's iRT at the chosen scan
    irt_pred::Float32                         # Library iRT for this pid
    rt_obs::Float32                           # Receiver's RT (min) at the chosen scan
    log_by_ratio::Float32
    smoothed_frag_sqrt::NTuple{8, Float32}
    prec_mz::Float32

    # Precomputed MBR_*_true features (the _false counterfactual is filled
    # downstream by the existing partner-pool machinery).
    MBR_max_pair_prob_true::Float32           # = seed donor's trace_prob
    MBR_log2_weight_ratio_true::Float32
    MBR_log2_explained_ratio_true::Float32
    MBR_best_irt_diff_true::Float32
    MBR_best_rt_diff_true::Float32
    MBR_log_by_diff_true::Float32
    MBR_smoothed_frag_hellinger_true::Float32
    MBR_is_missing_true::Bool                 # Always false for a seeded row

    # Markers
    mbr_recovered_seed::Bool                  # Always true
    trace_prob::Float32                       # Always 0 — signal to FTR that this row
                                              # has no intrinsic LightGBM identification.
end

"""
Parameters for MBRRecoverySearch.

Donor pool gating, window safety floor, and quad-window/isotope settings that
match MainSearch so the design-matrix candidate sets are comparable.
"""
struct MBRRecoverySearchParameters <: SearchParameters
    donor_qval_threshold::Float32             # 0.01 (mirrors MBR_DONOR_Q_THRESHOLD)
    max_donors_per_pid::Int                   # 2 (mirrors MBR_MAX_DONOR_FILES_PER_PRECURSOR)
    rt_tol_floor_min::Float32                 # 0.05 — applied if calibrated tol returns 0
    isotope_err_bounds::Tuple{UInt8, UInt8}   # (1, 0) — same as MainSearch
    min_fraction_transmitted::Float32         # quad transmission cutoff
end

"""
Results container for MBRRecoverySearch.

Per-file recovery is written to disk as a sidecar (see sidecar.jl); the
in-memory results track totals for the post-run summary log plus the
once-built, file-independent state the per-file `process_file!` calls reuse.

The donor pool, donor-pid set, and fragment-annotation keys are built on the
first `process_file!` call (`initialized` guards the build) and persist across
files because `reset_results!` is a no-op for this method.
"""
mutable struct MBRRecoverySearchResults <: SearchResults
    n_donor_pids::Int                         # Distinct pids in the donor pool
    n_cells_attempted::Int                    # Sum of "cells we tried to recover"
    n_cells_emitted::Int                      # Sum of "cells with positive-weight seed"
    elapsed_sec::Float64
    t_start::Float64                          # time() at init, for elapsed reporting

    # Built once on first process_file!; reused for every file.
    initialized::Bool
    donor_pool::Dict{UInt32, Vector{_MBRDonorEntry}}
    donor_pid_set::Set{UInt32}
    fragment_keys::_MBRFragmentAnnotationKeys
end

# Constants shared with the rest of the recovery code.
const MBR_RECOVERY_DONOR_QVAL = Float32(0.01)
const MBR_RECOVERY_MAX_DONORS = 2
const MBR_RECOVERY_RT_TOL_FLOOR_MIN = Float32(0.05)
const RECOVERY_SEED_SIDECAR_SUFFIX = ".recovery_seed.arrow"
