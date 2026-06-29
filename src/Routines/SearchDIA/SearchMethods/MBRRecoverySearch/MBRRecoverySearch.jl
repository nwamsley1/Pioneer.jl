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

Top-level entry. Builds the donor pool once from the post-PrecursorScoringSearch
per-file PSM tables, then per file:
  - identifies pids missing from this receiver file's main_search_psms
  - (V1.0) runs extract_max_weight_in_window! per cell to compute the recovered
    weight + receiver-side payload
  - (V1.0) computes the 7 MBR_*_true features inline via compute_seed_mbr_features
  - writes recovery_seed_sidecar.arrow next to the main file

V0 (this commit): builds the donor pool + identifies missing cells per file
and writes an EMPTY sidecar per file so the on-disk shape lines up. The
actual extraction call is gated behind the per-cell scratch-buffer wiring,
which depends on SearchData lifecycle access from PrecursorScoringSearch's
output. The extraction.jl function is in place and unit-testable in isolation
once that wiring lands.
"""
function perform_recovery_search!(
    method::MBRRecoverySearch,
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
)
    t0 = time()

    spec_lib = getSpecLib(search_context)
    psm_paths = list_post_scoring_psm_files(search_context)
    if isempty(psm_paths)
        @user_warn "MBRRecoverySearch: no post-scoring PSM files found — skipping"
        results.elapsed_sec = time() - t0
        return nothing
    end

    donor_pool = build_recovery_donor_pool(
        psm_paths, spec_lib;
        qval_threshold = params.donor_qval_threshold,
        max_donors = params.max_donors_per_pid,
    )
    results.n_donor_pids = length(donor_pool)
    donor_pid_set = Set{UInt32}(keys(donor_pool))

    for (file_idx, psm_path) in enumerate(psm_paths)
        main_pids = main_search_pid_set(psm_path)
        missing_pids = setdiff(donor_pid_set, main_pids)
        results.n_cells_attempted += length(missing_pids)

        # V0: write an empty sidecar to lock the on-disk format.
        # V1.0 will replace this with the per-cell extraction loop:
        #
        #   seeds = RecoveredSeedRow[]
        #   for donor_pid in missing_pids
        #       donor = first(donor_pool[donor_pid])  # best by trace_prob
        #       center_rt = inverse(rt_irt_model)(donor.irt_obs)
        #       receiver = extract_max_weight_in_window!(
        #           search_data, spectra, donor_pid, main_pids, rt_index,
        #           center_rt, rt_tol, rt_irt_model, mass_err, quad, nce,
        #           precursors, frag_lookup, iso_splines,
        #           params.isotope_err_bounds, params.min_fraction_transmitted,
        #           irt_tol, rt_binned_tol, max_iter_outer, max_diff, deconv_solver,
        #       )
        #       receiver === nothing && continue
        #       features = compute_seed_mbr_features(
        #           receiver, donor, getIrt(precursors)[donor_pid], fragment_keys,
        #       )
        #       push!(seeds, _build_recovered_seed_row(donor_pid, donor, receiver, features))
        #   end
        #   results.n_cells_emitted += length(seeds)
        seeds = RecoveredSeedRow[]
        sidecar_path = recovery_seed_sidecar_path(psm_path)
        write_recovery_seed_sidecar(sidecar_path, seeds)
        @user_info "  $(basename(psm_path)): $(length(missing_pids)) missing-pid cells, " *
                   "$(length(seeds)) seeds emitted (V0 stub)"
    end

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
