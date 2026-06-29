# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

#==========================================================
MBRRecoverySearch — donor-driven targeted PSM recovery.

For every (donor pid × file) cell where the donor was identified at
q ≤ MBR_RECOVERY_DONOR_QVAL in ≥1 file but is missing from this file's
main_search_psms, run a targeted PMM deconvolution in the file's
calibrated RT window around the donor's projected RT. Emit a sidecar
of recovered PSM-like rows with the seven MBR_*_true features
precomputed inline.

PrecursorScoringSearch later vertical-concatenates the sidecar into the
unified PSM table, fills missing columns with sentinels, and feeds the
augmented set through FTR (retrained to use the recovery rows + the
explicit :mbr_recovered_seed boolean feature).

Pipeline position: AFTER PrecursorScoringSearch (which already computes per-file
qval, trace_prob, and writes the unified PSM table) and BEFORE
IntegrateChromatogramSearch. Source-of-truth files for donor selection are the
117-column main_search_psms/{file}.arrow that PrecursorScoringSearch writes.

The seven MBR_*_true features are computed inline at extraction time; the
_false counterfactual + the receiver-side FTR scoring is deferred to V1.5 (V0
emits recovery rows that pass through to downstream stages with
`mbr_recovered_seed = true` as the only quality gate).
==========================================================#

struct MBRRecoverySearch <: SearchMethod end

#==========================================================
Interface
==========================================================#

get_search_method_name(::MBRRecoverySearch) = "MBR Recovery"

function get_parameters(::MBRRecoverySearch, params::Any)
    return MBRRecoverySearchParameters(
        MBR_RECOVERY_DONOR_QVAL,
        MBR_RECOVERY_MAX_DONORS,
        MBR_RECOVERY_RT_TOL_FLOOR_MIN,
        (UInt8(1), UInt8(0)),
        Float32(params.global.isotope_settings.min_fraction_transmitted),
    )
end

function init_search_results(::MBRRecoverySearch, ::MBRRecoverySearchParameters, ::SearchContext)
    return MBRRecoverySearchResults(0, 0, 0, 0.0)
end

reset_results!(::MBRRecoverySearchResults) = nothing

#==========================================================
Core processing
==========================================================#

"""
    perform_recovery_search!(method, results, params, search_context)

Top-level entry. Builds the donor pool once from the (post-prescore) per-file
PSM tables, then loops over files for the per-cell extraction. Each file
gets a `recovery_seed_sidecar.arrow` written to `temp_data/main_search_psms/`.

V0: stubbed. Real implementation comes in donor_pool.jl + extraction.jl +
sidecar.jl.
"""
function perform_recovery_search!(
    method::MBRRecoverySearch,
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
)
    t0 = time()
    @user_info "MBRRecoverySearch: stub — donor pool + extraction not yet wired"
    # TODO: see donor_pool.jl + extraction.jl + sidecar.jl
    results.elapsed_sec = time() - t0
    return nothing
end

#==========================================================
SearchMethod hooks
==========================================================#

# MBRRecoverySearch doesn't iterate per-file via process_file!. It overrides
# performSearch! to coordinate the (donor pool build) → (per-file extraction)
# flow itself.

function process_file!(
    ::MBRRecoverySearch, ::MBRRecoverySearchResults, ::MBRRecoverySearchParameters,
    ::SearchContext, ::Int64, ::Any
)
    return nothing
end

function process_search_results!(
    ::MBRRecoverySearchResults, ::MBRRecoverySearchParameters,
    ::SearchContext, ::Int64, ::Any
)
    return nothing
end

function summarize_results!(
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
)
    perform_recovery_search!(MBRRecoverySearch(), results, params, search_context)
    @user_info "MBRRecoverySearch done — pids=$(results.n_donor_pids) " *
               "cells_attempted=$(results.n_cells_attempted) " *
               "cells_emitted=$(results.n_cells_emitted) " *
               "elapsed=$(round(results.elapsed_sec, digits=2))s"
    return nothing
end
