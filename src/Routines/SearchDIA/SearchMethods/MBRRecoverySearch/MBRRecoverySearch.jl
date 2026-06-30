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
    # min_fraction_transmitted is hardcoded to 0.25 to match MainSearch and
    # IntegrateChromatogramSearch (the former global.isotope_settings knob was
    # removed — every shipping config used this value), so the per-scan
    # candidate sets are comparable across stages.
    return MBRRecoverySearchParameters(
        MBR_RECOVERY_DONOR_QVAL,
        MBR_RECOVERY_MAX_DONORS,
        MBR_RECOVERY_RT_TOL_FLOOR_MIN,
        (UInt8(1), UInt8(0)),
        Float32(0.25),
    )
end

function init_search_results(::MBRRecoverySearch, ::MBRRecoverySearchParameters, ::SearchContext)
    return MBRRecoverySearchResults(
        0, 0, 0, 0.0, time(),
        false,
        Dict{UInt32, Vector{_MBRDonorEntry}}(),
        Set{UInt32}(),
        _MBRFragmentAnnotationKeys(UInt16[]),
    )
end

# No-op: the cached donor pool / pid set / fragment keys must survive between
# per-file process_file! calls, so we deliberately keep results across files.
reset_results!(::MBRRecoverySearchResults) = nothing

#==========================================================
Core processing
==========================================================#

"""
    ensure_recovery_initialized!(results, params, search_context)

Builds the file-independent recovery state once and caches it on `results`:
the donor pool (pid → top donor entries by trace_prob, from every
post-PrecursorScoringSearch PSM table), the donor-pid set, and the library
fragment-annotation keys used by the inline Hellinger feature. Guarded by
`results.initialized`, so it runs on the first `process_file!` call only.
"""
function ensure_recovery_initialized!(
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
)
    results.initialized && return nothing
    results.initialized = true

    spec_lib = getSpecLib(search_context)
    psm_paths = list_post_scoring_psm_files(search_context)
    if isempty(psm_paths)
        @user_warn "MBRRecoverySearch: no post-scoring PSM files found — skipping"
        return nothing
    end

    pool = build_recovery_donor_pool(
        psm_paths, spec_lib;
        qval_threshold = params.donor_qval_threshold,
        max_donors = params.max_donors_per_pid,
    )
    results.donor_pool = pool
    results.donor_pid_set = Set{UInt32}(keys(pool))
    results.n_donor_pids = length(pool)
    results.fragment_keys = build_mbr_fragment_annotation_keys(getFragmentLookupTable(spec_lib))
    return nothing
end

"""
    _build_recovered_seed_row(pid, donor, receiver, feats, precursors, prec_mz_arr, prec_irt_arr) -> RecoveredSeedRow

Assembles one recovery row from the per-cell PMM extraction (`receiver`), the
seeding donor entry, and the inline MBR_*_true features (`feats`).
"""
@inline function _build_recovered_seed_row(
    pid::UInt32,
    donor::_MBRDonorEntry,
    receiver::ReceiverExtraction,
    feats,
    precursors,
    prec_mz_arr,
    prec_irt_arr,
)
    return RecoveredSeedRow(
        pid,
        !donor.is_decoy,                       # target (donor shares this pid)
        receiver.scan_idx,
        getCvFold(precursors, pid),
        receiver.weight,
        receiver.log2_intensity_explained,
        receiver.irt_obs,
        prec_irt_arr[pid],                     # irt_pred (library)
        receiver.rt_obs,
        receiver.log_by_ratio,
        receiver.smoothed_frag_sqrt,
        prec_mz_arr[pid],
        feats.MBR_max_pair_prob_true,
        feats.MBR_log2_weight_ratio_true,
        feats.MBR_log2_explained_ratio_true,
        feats.MBR_best_irt_diff_true,
        feats.MBR_best_rt_diff_true,
        feats.MBR_log_by_diff_true,
        feats.MBR_smoothed_frag_hellinger_true,
        feats.MBR_is_missing_true,
        true,                                  # mbr_recovered_seed
        0.0f0,                                 # trace_prob (no intrinsic ID)
    )
end

"""
    extract_file_recovery_seeds(results, params, search_context, ms_file_idx, spectra) -> Vector{RecoveredSeedRow}

Per-file extraction. Partitions this file's missing donor pids across
`nthreads()` tasks, each using its own `getSearchData(...)[thread_id]` scratch,
and runs `extract_max_weight_in_window!` + inline MBR feature compute per cell.
"""
function extract_file_recovery_seeds(
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
    rt_index::retentionTimeIndex,
    main_pids::Set{UInt32},
    missing_pids::Vector{UInt32},
)
    n = length(missing_pids)
    n == 0 && return RecoveredSeedRow[]

    # Per-file calibrated models + library handles (read-only, shared).
    rt_irt_model     = getRtIrtModel(search_context, ms_file_idx)
    nce_model        = getNceModel(search_context, ms_file_idx)
    mass_error_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model       = getQuadTransmissionModel(search_context, ms_file_idx)
    irt_tol          = Float32(getIrtErrors(search_context)[ms_file_idx])
    rt_binned_tol    = haskey(getRtTolerances(search_context), ms_file_idx) ?
                       getRtTolerance(search_context, ms_file_idx) : nothing

    spec_lib     = getSpecLib(search_context)
    precursors   = getPrecursors(spec_lib)
    frag_lookup  = getFragmentLookupTable(spec_lib)
    prec_mz_arr  = getMz(precursors)
    prec_irt_arr = getIrt(precursors)

    donor_pool    = results.donor_pool
    fragment_keys = results.fragment_keys
    deconv_solver = PoissonMMSolver()

    # Contiguous partition of the missing-pid vector across ≤ nthreads tasks.
    nt = max(min(Threads.nthreads(), n), 1)
    bounds = [(k * n) ÷ nt for k in 0:nt]   # task t covers bounds[t]+1 : bounds[t+1]

    tasks = map(1:nt) do thread_id
        Threads.@spawn begin
            search_data = getSearchData(search_context)[thread_id]
            iso_splines = getIsoSplines(search_data)
            local_seeds = RecoveredSeedRow[]
            @inbounds for j in (bounds[thread_id] + 1):bounds[thread_id + 1]
                pid = missing_pids[j]
                entries = donor_pool[pid]
                isempty(entries) && continue
                donor = entries[1]   # best donor by trace_prob

                receiver = extract_max_weight_in_window!(
                    search_data, spectra, pid, main_pids, rt_index,
                    donor.irt_obs, rt_irt_model,
                    mass_error_model, quad_model, nce_model,
                    precursors, frag_lookup, iso_splines,
                    params.isotope_err_bounds, params.min_fraction_transmitted,
                    irt_tol, rt_binned_tol,
                    DECONV_MAX_ITER, DECONV_CONVERGENCE_TOL, deconv_solver,
                )
                receiver === nothing && continue

                feats = compute_seed_mbr_features(receiver, donor, prec_irt_arr[pid], fragment_keys)
                push!(local_seeds, _build_recovered_seed_row(
                    pid, donor, receiver, feats, precursors, prec_mz_arr, prec_irt_arr,
                ))
            end
            local_seeds
        end
    end

    return reduce(vcat, fetch.(tasks); init = RecoveredSeedRow[])
end

#==========================================================
SearchMethod hooks
==========================================================#

"""
Per-file recovery. The framework hands us the receiver file's `spectra`; we
build the donor pool once (lazily), find the donor pids missing from this
file's main_search_psms, extract each in a targeted PMM window, and write the
recovery seed sidecar next to the main file.
"""
function process_file!(
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
)
    ensure_recovery_initialized!(results, params, search_context)
    isempty(results.donor_pid_set) && return nothing

    # Map this file index to its post-scoring PSM table.
    file_name = getParsedFileName(search_context, ms_file_idx)
    psm_path = joinpath(getDataOutDir(search_context), "temp_data",
                        "main_search_psms", file_name * ".arrow")
    if !isfile(psm_path)
        @debug_l2 "MBRRecoverySearch: no main_search_psms for $file_name — skipping"
        return nothing
    end

    main_pids = main_search_pid_set(psm_path)
    missing_pids = collect(setdiff(results.donor_pid_set, main_pids))
    results.n_cells_attempted += length(missing_pids)

    sidecar_path = recovery_seed_sidecar_path(psm_path)

    # Need the per-file RT index (empirical RT space) to enumerate competitors,
    # plus a calibrated per-file iRT tolerance to size the recovery window. A
    # file lacking either was not fully calibrated upstream (e.g. too few PSMs);
    # skip it with an empty sidecar.
    rt_index_path = getRtIndex(getMSData(search_context), ms_file_idx)
    if isempty(rt_index_path) || !isfile(rt_index_path) ||
       !haskey(getIrtErrors(search_context), ms_file_idx)
        @debug_l2 "MBRRecoverySearch: no rt_index / iRT tol for $file_name — emitting empty sidecar"
        write_recovery_seed_sidecar(sidecar_path, RecoveredSeedRow[])
        return nothing
    end
    rt_index = buildRtIndex(DataFrame(Arrow.Table(rt_index_path)), bin_rt_size = 0.1)

    seeds = extract_file_recovery_seeds(
        results, params, search_context, ms_file_idx,
        spectra, rt_index, main_pids, missing_pids,
    )
    results.n_cells_emitted += length(seeds)

    write_recovery_seed_sidecar(sidecar_path, seeds)
    @user_info "  $file_name: $(length(missing_pids)) missing-pid cells, " *
               "$(length(seeds)) seeds emitted"
    return nothing
end

function process_search_results!(
    ::MBRRecoverySearchResults, ::MBRRecoverySearchParameters,
    ::SearchContext, ::Int64, ::MassSpecData,
)
    return nothing
end

function summarize_results!(
    results::MBRRecoverySearchResults,
    params::MBRRecoverySearchParameters,
    search_context::SearchContext,
)
    results.elapsed_sec = time() - results.t_start
    @user_info "MBRRecoverySearch done — pids=$(results.n_donor_pids) " *
               "cells_attempted=$(results.n_cells_attempted) " *
               "cells_emitted=$(results.n_cells_emitted) " *
               "elapsed=$(round(results.elapsed_sec, digits=2))s"
    return nothing
end
