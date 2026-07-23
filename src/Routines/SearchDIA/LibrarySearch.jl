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

# Toggle to profile steps with PProf during MainSearch only.
# Saves a flamegraph .pb.gz per file to the output directory.
const PROFILE_FRAG_INDEX = false   # fragment index search (searchFragmentIndexPartitionMajorHinted)
const PROFILE_DECONV = false        # deconvolution (process_scans!)

"""
    library_search(spectra, search_context, params, ms_file_idx) -> DataFrame

Single entry point for all fragment-index-based searches (ParameterTuning,
QuadTuning, MainSearch, etc.).

Three dispatch helpers (defined in process_scans.jl) handle the behavioral differences:

| Helper              | ParameterTuning          | Default / MainSearch         |
|---------------------|--------------------------|------------------------------|
| `get_fragment_index`| presearch (small) index  | full partitioned index       |
| `get_irt_tolerance` | `params.irt_tol`         | calibrated per-file or `Inf` |
| `get_nce_models`    | 1 calibrated model       | 1 calibrated model           |

Pipeline within this function:
1. Extract per-file models (mass error, RT, quad transmission) from SearchContext
2. Fragment index search — find candidate precursors for each scan (runs once)
3. For each NCE model, spawn threads calling `process_scans!` to score candidates
4. Concatenate results across threads and NCE models into a single DataFrame
"""
function library_search(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64;
    scan_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    fragment_index = nothing,
    max_peaks::Int = 0,
) where {P<:FragmentIndexSearchParameters}

    # --- 1. Extract per-file models and library data ---
    spec_lib = getSpecLib(search_context)
    search_data = getSearchData(search_context)

    irt_tol = get_irt_tolerance(search_context, params, ms_file_idx)
    # Use presearch index (no iRT bins) when iRT tolerance is infinite (RT model not fitted).
    # This avoids iterating all fine iRT bins for no benefit.
    partitioned_index = if fragment_index !== nothing
        fragment_index
    elseif irt_tol >= typemax(Float32)
        getPresearchPartitionedIndex(spec_lib)
    else
        get_fragment_index(spec_lib, params)
    end
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    # Sciex ZT scanning DIA (PIONEER_ZT_METASCAN_K>0): use two square quad models.
    # The fragment index gets a tight 1-bin box (overhang 0) for center-bin matching;
    # deconvolution gets a wide box spanning the ±k-bin meta-scan (overhang = k·S) so
    # every bin passes full transmission for any in-meta-scan precursor and the deconv
    # weights reflect the actual (triangular) ion current instead of being masked by a
    # narrow fitted quad. qtm_deconv is finalized after S is measured (below).
    _zt_k = something(tryparse(Int, get(ENV, "PIONEER_ZT_METASCAN_K", "")), 0)
    _zt = _zt_k > 0
    # Only the MAIN search expands to meta-scans + uses the wide deconv box.
    # Parameter/quad/NCE tuning must stay on the clean center-bin (1 m/z) path
    # so calibrations aren't fit on the metascan-expanded, wide-box deconvolution.
    _zt_main = _zt && (params isa MainSearchParameters)
    # ZT fragment-index candidacy box. Default overhang 0 = tight 1-bin (isolation-width) box:
    # a precursor is a candidate only where its M0 falls in a single bin. PIONEER_ZT_FRAG_OVERHANG
    # (m/z) widens it so a precursor whose transmission is split across bins (offset near a bin
    # boundary) — or whose M0-bin signal alone is below the index-score threshold — gets candidacy
    # chances in the neighboring bins too. Center assignment (filter_to_center_bin!) stays tight.
    # Default 0.5 m/z overhang (≈±1 bin) — a 3-file Condition-A sweep showed the M0-only tight
    # box (0.0) was dropping split/marginal-M0 precursors at candidacy: 0.0→77,549, 0.5→78,432
    # (+883 prec, +1.14%, 1% FDR), 1.0→78,362 (plateau). Center assignment stays tight.
    _zt_frag_overhang = something(tryparse(Float32, get(ENV, "PIONEER_ZT_FRAG_OVERHANG", "")), 0.5f0)
    qtm_frag = _zt ? SquareQuadModel(_zt_frag_overhang) : qtm
    # Experiment: ZT ParameterTuning is PSM-starved because each ~1 Da reconstructed bin carries
    # fragments from the FULL ~5 Da swept-quad isolation, but the tight candidacy box only matches
    # precursors whose M0 sits in that bin. PIONEER_ZT_TUNING_OVERHANG (m/z) widens the tuning
    # candidacy box to the physical quad width (overhang≈2 → ±2.5 Da = 5 Da window) so tuning finds
    # enough PSMs for a well-constrained mass/RT/NCE calibration. ParameterTuning only; 0 = off.
    _zt_tuning = _zt && (params isa ParameterTuningSearchParameters)
    _zt_tuning_overhang = something(tryparse(Float32, get(ENV, "PIONEER_ZT_TUNING_OVERHANG", "")), 0.0f0)
    if _zt_tuning && _zt_tuning_overhang > 0f0
        qtm_frag = SquareQuadModel(_zt_tuning_overhang)
    end
    # Sciex ZT wide-emit candidacy (MAIN search): PIONEER_ZT_CANDIDACY_TOL (Da half-width) evaluates
    # each precursor across ~2·tol bins so it gets a bitvec emission chance in EVERY metascan bin,
    # not just its center. A precursor that clears in ANY bin survives (map_any_hit_to_center! below,
    # then the usual ±k fill) — recovering diluted precursors without lowering the bitvec threshold.
    # 0 = off (current filter_to_center path). isoWidth/2 ≈ 0.51, so overhang = tol − 0.51 → half-width ≈ tol.
    _zt_cand_tol = something(tryparse(Float32, get(ENV, "PIONEER_ZT_CANDIDACY_TOL", "")), 0.0f0)
    _zt_wide_emit = _zt_main && _zt_cand_tol > 0.51f0
    if _zt_wide_emit
        qtm_frag = SquareQuadModel(_zt_cand_tol - 0.51f0)
    end
    # Sciex ZT candidacy merge (no-reset counter accumulation): PIONEER_ZT_MERGE_K (bins).
    # For each center MS2 scan, the fragment-index scores its ±K same-cycle neighbors'
    # peaks into the SAME bitvec counter before emitting — bit-identical to merging the
    # ±K bins' peaks (idempotent OR) but with no buffer/sort and no precursor-window
    # widening. Recovers swept-quad-diluted precursors at the default bitvec threshold.
    # MAIN search only; keeps the tight center-bin path (filter_to_center_bin! + ±k
    # expand + baseline deconv). Independent of PIONEER_ZT_CANDIDACY_TOL (wide-emit).
    _zt_merge_k = _zt_main ? something(tryparse(Int, get(ENV, "PIONEER_ZT_MERGE_K", "")), 0) : 0
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt = getRtIrtModel(search_context, ms_file_idx)
    precursors = getPrecursors(spec_lib)
    ion_list = getFragmentLookupTable(spec_lib)

    # NCE models to iterate: [(calibrated_model, nothing)]
    nce_entries = get_nce_models(search_context, params, ms_file_idx)
    isempty(nce_entries) && return DataFrame()

    # --- 2. Fragment index search (runs once, shared across all NCE models) ---
    if scan_indices === nothing
        # Default: partition all scans across threads
        thread_tasks = partition_scans(spectra, Threads.nthreads())
        all_scan_idxs = Int[]
        for tt in thread_tasks
            append!(all_scan_idxs, last(tt))
        end
        filter!(si -> si > 0 && si <= length(spectra) &&
            getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)
    else
        # Use provided scan indices, filter to valid MS2 scans
        all_scan_idxs = filter(si -> si > 0 && si <= length(spectra) &&
            getMsOrder(spectra, si) ∈ getSpecOrder(params), Int.(scan_indices))
        # Partition provided scans evenly across threads
        n_threads = Threads.nthreads()
        thread_tasks = [(i, Int[]) for i in 1:n_threads]
        for (idx, si) in enumerate(all_scan_idxs)
            push!(thread_tasks[mod1(idx, n_threads)][2], si)
        end
        filter!(tt -> !isempty(last(tt)), thread_tasks)
    end

    isempty(all_scan_idxs) && return DataFrame()

    # Build score filter: use learned LUT if available, otherwise fall back to count_ones
    bitvec_filter = getBitVecFilter(search_context, ms_file_idx)
    score_filter = if bitvec_filter !== nothing
        LUTFilter(bitvec_filter)
    else
        CountFilter(getMinIndexSearchScore(params))
    end

    # scan_to_prec_idx[scan] = range into precursors_passed (or missing if no candidates)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    t_frag_start = time()
    if PROFILE_FRAG_INDEX && params isa MainSearchParameters
        Profile.clear()
        Profile.init(n=50_000_000, delay=0.0005)
        Profile.@profile precursors_passed, scores_passed = searchFragmentIndexPartitionMajorHinted(
            scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
            Threads.nthreads(), params, qtm_frag, mem, rt_to_irt, irt_tol,
            getMz(precursors);
            score_filter = score_filter, max_peaks = max_peaks,
            scratch = getFragIndexScratch(search_context),
            iso_bounds_override = _zt_main ? (UInt8(0), UInt8(0)) : nothing,
            zt_merge_k = _zt_merge_k)
        out_dir = getDataOutDir(search_context)
        prof_path = joinpath(out_dir, "frag_index_profile_$(ms_file_idx).pb.gz")
        pprof(out=prof_path, web=false)
        @user_info "Fragment index profile saved to $prof_path\n"
    else
        precursors_passed, scores_passed = searchFragmentIndexPartitionMajorHinted(
            scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
            Threads.nthreads(), params, qtm_frag, mem, rt_to_irt, irt_tol,
            getMz(precursors);
            score_filter = score_filter, max_peaks = max_peaks,
            scratch = getFragIndexScratch(search_context),
            iso_bounds_override = _zt_main ? (UInt8(0), UInt8(0)) : nothing,
            zt_merge_k = _zt_merge_k)
    end
    t_frag = time() - t_frag_start

    # --- DEBUG: dump fragment index bitmask scores to Arrow and bail ---
    # Only dump during MainSearch, not tuning stages.
    # Applies RT and precursor m/z filtering (same as selectTransitions!) for realistic counts.
    # --- 2a''. Sciex ZT: filter the emitted candidates down to the center bin.
    # The partitioned index can emit precursors slightly outside a bin's ±half-bin
    # window (partition granularity); keep only |prec_mz - centerMz| ≤ isoWidth/2 so
    # each (precursor, scan) is a clean single-center-bin assignment before the ±k
    # meta-scan expansion. ---
    if _zt_main
        n_before_cb = length(precursors_passed)
        if _zt_wide_emit
            # Wide-emit: any precursor that cleared the bitvec in ANY bin → mapped to its center
            # bin (expand_to_metascans! then fills its ±k, same as the tight path). search_halfbins
            # spans the ±tol emission range so the true nearest bin is found.
            _S_hb = _zt_bin_step(spectra, all_scan_idxs)
            _search_hb = ceil(Int, _zt_cand_tol / max(_S_hb, 1f-3)) + 1
            precursors_passed = map_any_hit_to_center!(
                scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, getMz(precursors), _search_hb)
            @debug_l1 "ZT wide-emit (tol=$(_zt_cand_tol) Da, hb=$_search_hb): $n_before_cb emissions → $(length(precursors_passed)) center candidates"
        else
            precursors_passed = filter_to_center_bin!(
                scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, getMz(precursors))
            @debug_l1 "ZT center-bin filter: $n_before_cb → $(length(precursors_passed)) candidates"
        end
    end

    # --- 2a. Pre-filter: require candidates to appear in ≥ N scans ---
    prefilter_n = getPrefilterMinScanCount(params)
    if prefilter_n > 1
        n_before = length(precursors_passed)
        precursors_passed = filter_low_scan_candidates!(
            scan_to_prec_idx, precursors_passed, prefilter_n)
        @debug_l1 "Pre-filter: $n_before → $(length(precursors_passed)) candidates " *
              "($(round(100*(1 - length(precursors_passed)/max(1,n_before)), digits=1))% removed)"
    end

    # --- 2a'. Sciex ZT scanning DIA: expand each candidate to its ±k adjacent
    # Q1 bins (same cycle). A precursor's ions spread across a meta-scan of
    # 2k+1 bins; the fragment index only matched its center bin, so we add it
    # as a candidate to the neighboring bins before deconvolution estimates a
    # per-bin weight. Off (k=0) by default. ---
    # Tuning stages (ZT but not main search) deconvolve the clean center bins with
    # a 1 m/z square box; only the main search expands + widens the box.
    qtm_deconv = _zt ? SquareQuadModel(0.0f0) : qtm
    if _zt_tuning && _zt_tuning_overhang > 0f0
        qtm_deconv = SquareQuadModel(_zt_tuning_overhang)   # keep tuning deconv consistent with the wide candidacy box
    end
    if _zt_main
        n_before_ms = length(precursors_passed)
        precursors_passed = expand_to_metascans!(
            scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, _zt_k)
        # Deconv transmission: prefer the QuadTuning-fitted empirical LUT (the real
        # swept-quad shape); fall back to the wide square box (overhang = k·S spanning
        # the whole (2k+1)-bin meta-scan) if QuadTuning didn't produce one.
        S_est = _zt_bin_step(spectra, all_scan_idxs)
        qtm_deconv = qtm isa EmpiricalQuadModel ? qtm : SquareQuadModel(Float32(_zt_k) * S_est)
        @debug_l1 "ZT meta-scan expansion (k=$_zt_k): $n_before_ms → " *
                  "$(length(precursors_passed)) candidates " *
                  "($(round(length(precursors_passed)/max(1,n_before_ms), digits=2))×); " *
                  "S=$(round(S_est,digits=3)); deconv=$(qtm_deconv isa EmpiricalQuadModel ? "EmpiricalQuadModel(LUT)" : "SquareQuadModel($(round(S_est*(2*_zt_k+1),digits=1)) m/z)")"
    end

    # --- 2b. Build precursor index ---
    prec_index = PerScanPrecursorIndex(scan_to_prec_idx, precursors_passed)

    # --- 3zt. Sciex ZT main search (5b): batched deconvolve-once + per-batch collapse.
    # Deconvolve CONTIGUOUS scan batches (each scan once), run the batch-safe pre-collapse
    # passes (scan_competition, ms1_lookup), then collapse each batch borrowing ±k boundary
    # rows from neighbors — vcat ~3M meta-PSMs instead of ~35M raw rows. The caller
    # (process_file!/process_search_results!) skips its own scan_competition/ms1/permute
    # and the global collapse for this path. Byte-identical to the global collapse. ---
    # ZT disk-spill (chosen approach): round-robin deconv UNCHANGED, spill raw to on-disk
    # scan-range bins, collapse one bin at a time. Byte-identical, low peak. PIONEER_ZT_SPILL=1.
    if _zt_main && get(ENV, "PIONEER_ZT_SPILL", "0") == "1"
        t_zt_start = time()
        n_bins = something(tryparse(Int, get(ENV, "PIONEER_ZT_SPILL_BINS", "")), 16)
        result = zt_disk_spill_deconv_collapse(
            thread_tasks, nce_entries, spectra, prec_index, ms_file_idx, search_context,
            params, precursors, ion_list, qtm_deconv, mem, rt_to_irt, irt_tol, _zt_k,
            all_scan_idxs, n_bins)
        @debug_l1 "  library_search (ZT disk-spill, bins=$n_bins): frag_index=$(round(t_frag, digits=2))s  " *
                   "deconv+collapse=$(round(time()-t_zt_start, digits=2))s  meta_psms=$(nrow(result))"
        return result
    end

    if _zt_main && get(ENV, "PIONEER_ZT_BATCHED", "0") == "1"
        t_zt_start = time()
        zt_tasks = partition_scans_contiguous_zt(all_scan_idxs, Threads.nthreads())
        result = zt_batched_deconv_collapse(
            zt_tasks, nce_entries, spectra, prec_index, ms_file_idx, search_context,
            params, precursors, ion_list, qtm_deconv, mem, rt_to_irt, irt_tol, _zt_k)
        @debug_l1 "  library_search (ZT batched): frag_index=$(round(t_frag, digits=2))s  " *
                   "deconv+collapse=$(round(time()-t_zt_start, digits=2))s  " *
                   "meta_psms=$(nrow(result))  candidates=$(length(get_precursors(prec_index)))"
        return result
    end

    # --- 3. Threaded scan processing, once per NCE model ---
    # All search methods use the fused per-precursor scan loop. The classic
    # `process_scans!` path was deleted along with the multi-pass pipeline
    # (selectTransitions! + matchPeaks! + buildDesignMatrix! + sortSparse!).
    # When nce_tag is not nothing (NCE tuning), tag each result with the NCE value.
    t_deconv_start = time()
    all_results = map(nce_entries) do (nce_model, nce_tag)
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn process_scans_fused!(
                last(thread_task), spectra, prec_index,
                ms_file_idx,
                search_data[first(thread_task)], params, precursors, ion_list,
                nce_model, qtm_deconv, mem, rt_to_irt, irt_tol)
        end
        # Unwrap TaskFailedException so the real error surfaces instead of
        # being buried inside a Task wrapper.
        fetched = map(tasks) do t
            try
                fetch(t)
            catch e
                while e isa TaskFailedException
                    e = e.task.exception
                end
                rethrow(e)
            end
        end
        result = vcat(fetched...)
        if nce_tag !== nothing && !isempty(result)
            result[!, :nce] .= nce_tag
        end
        return result
    end
    t_deconv = time() - t_deconv_start

    t_post_start = time()
    # Single NCE model is the common case; `vcat(all_results...)` would copy the
    # whole (large) PSM table for nothing. all_results[1] is already the
    # materialized per-NCE table (from the inner vcat(fetched...)), so return it
    # directly. Only concatenate when there are multiple NCE models.
    result = length(all_results) == 1 ? all_results[1] : vcat(all_results...)
    t_vcat = time() - t_post_start

    if params isa MainSearchParameters
        @debug_l1 "  library_search breakdown: frag_index=$(round(t_frag, digits=2))s  " *
                   "deconv=$(round(t_deconv, digits=2))s  " *
                   "vcat=$(round(t_vcat, digits=2))s  " *
                   "candidates=$(length(get_precursors(prec_index)))"
        # General min-memory win: the per-thread scored scratch (~6 GB SoA) was just copied
        # into `result` by the vcat and is now stale. The copycols=false view DataFrames that
        # referenced it (the per-NCE `fetched` locals) are out of scope, so reassigning the
        # field to an empty StructArray drops the last reference — GC reclaims it under memory
        # pressure BEFORE process_search_results!'s prepare pass (the peak instant), on every
        # DIA run incl. non-ZT. Regrown next file by growScoredPSMs!. Main-search only: the
        # tuning stages call library_search many times, where regrow churn would not pay off.
        # Env gate (default on) for A/B benchmarking: PIONEER_FREE_SCORED_SCRATCH=0 disables.
        if get(ENV, "PIONEER_FREE_SCORED_SCRATCH", "1") != "0"
            for sd in search_data
                sd.main_search_scored_psms = StructArray(MainSearchScoredPSM{Float32,Float16}[])
            end
        end
    end

    return result
end

"""
    filter_low_scan_candidates!(scan_to_prec_idx, precursors_passed, min_scan_count)

Remove fragment index candidates whose precursor appears in fewer than `min_scan_count` scans.
Modifies `scan_to_prec_idx` in-place and returns a new `precursors_passed` vector.
"""
function filter_low_scan_candidates!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    min_scan_count::Int
)
    # Step 1: Count scans per precursor
    scan_counts = Dict{UInt32, Int}()
    for scan_idx in eachindex(scan_to_prec_idx)
        range = scan_to_prec_idx[scan_idx]
        ismissing(range) && continue
        for i in range
            pid = precursors_passed[i]
            scan_counts[pid] = get(scan_counts, pid, 0) + 1
        end
    end

    # Step 2: Rebuild filtered vectors
    new_passed = UInt32[]
    sizehint!(new_passed, length(precursors_passed))
    for scan_idx in eachindex(scan_to_prec_idx)
        range = scan_to_prec_idx[scan_idx]
        ismissing(range) && continue
        start = length(new_passed) + 1
        for i in range
            pid = precursors_passed[i]
            if get(scan_counts, pid, 0) >= min_scan_count
                push!(new_passed, pid)
            end
        end
        scan_to_prec_idx[scan_idx] = length(new_passed) >= start ?
            (start:length(new_passed)) : missing
    end
    return new_passed
end

"""
    filter_by_bitvec!(scan_to_prec_idx, precursors_passed, scores_passed, filter_table)

Remove fragment index candidates whose bitmask pattern does not pass the BitVec filter.
`filter_table` is a `Vector{Bool}` of length 256, indexed by `Int(score_mask) + 1`.
Modifies `scan_to_prec_idx` in-place and returns a new `precursors_passed` vector.
"""
function filter_by_bitvec!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    scores_passed::Vector{UInt8},
    filter_table::Vector{Bool}
)
    new_passed = UInt32[]
    sizehint!(new_passed, length(precursors_passed))
    for scan_idx in eachindex(scan_to_prec_idx)
        range = scan_to_prec_idx[scan_idx]
        ismissing(range) && continue
        start = length(new_passed) + 1
        for i in range
            if filter_table[Int(scores_passed[i]) + 1]
                push!(new_passed, precursors_passed[i])
            end
        end
        scan_to_prec_idx[scan_idx] = length(new_passed) >= start ?
            (start:length(new_passed)) : missing
    end
    return new_passed
end

"""
    filter_to_center_bin!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, prec_mzs)

Sciex ZT: keep only candidates whose precursor m/z lies within the bin's own
±(isolationWidth/2) window, i.e. its single center bin. Removes the frag index's
partition-granularity over-emission so the ±k meta-scan expansion is well-defined.
Rebuilds `precursors_passed` / reindexes `scan_to_prec_idx`; returns the new vector.
"""
# Sciex ZT wide-emit candidacy. The fragment index (run with a wide ±tol precursor box) emits a
# precursor in every metascan bin whose bitmask pattern cleared the threshold. Here we map each such
# emission to the precursor's CENTER bin (the same-cycle MS2 bin whose centerMz is nearest its m/z),
# deduped per (cycle, precursor). Net effect: a precursor survives if it cleared in ANY of its bins,
# vs filter_to_center_bin! which requires the center bin itself to clear. expand_to_metascans! then
# fills each survivor's ±k as usual, so the meta-scan collapse + shape features are unchanged.
# `search_halfbins` bounds the nearest-bin scan (covers the ±tol emission span).
function map_any_hit_to_center!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    prec_mzs::AbstractVector{Float32},
    search_halfbins::Int,
)
    nspec = length(spectra)
    center_sets = Dict{Int, Set{UInt32}}()          # center scan_idx -> precursors filled there
    seen = Set{Tuple{UInt32, UInt32}}()             # (cycle, precursor) dedup
    @inbounds for si in all_scan_idxs
        isassigned(scan_to_prec_idx, si) || continue
        rng = scan_to_prec_idx[si]
        ismissing(rng) && continue
        cyc = getCycleIdx(spectra, si)
        for r in rng
            p = precursors_passed[r]
            key = (UInt32(cyc), p)
            (key in seen) && continue
            push!(seen, key)
            pmz = Float32(prec_mzs[p])
            # nearest same-cycle MS2 bin to the precursor m/z, within ±search_halfbins of si
            best = 0
            bestd = Inf32
            for d in -search_halfbins:search_halfbins
                sj = si + d
                (sj < 1 || sj > nspec) && continue
                getMsOrder(spectra, sj) == 2 || continue
                getCycleIdx(spectra, sj) == cyc || continue
                cm = getCenterMz(spectra, sj)
                ismissing(cm) && continue
                dd = abs(Float32(cm) - pmz)
                if dd < bestd
                    bestd = dd
                    best = sj
                end
            end
            best == 0 && continue
            push!(get!(center_sets, best, Set{UInt32}()), p)
        end
    end
    new_precursors = UInt32[]
    sizehint!(new_precursors, sum(length, values(center_sets); init = 0))
    @inbounds for si in all_scan_idxs
        s = get(center_sets, si, nothing)
        if s === nothing || isempty(s)
            scan_to_prec_idx[si] = missing
        else
            start = length(new_precursors) + 1
            for p in s
                push!(new_precursors, p)
            end
            scan_to_prec_idx[si] = start:length(new_precursors)
        end
    end
    return new_precursors
end

function filter_to_center_bin!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    prec_mzs::AbstractVector{Float32},
)
    new_precursors = UInt32[]
    sizehint!(new_precursors, length(precursors_passed))
    @inbounds for si in all_scan_idxs
        rng = scan_to_prec_idx[si]
        ismissing(rng) && continue
        cv = getCenterMz(spectra, si)
        wv = getIsolationWidthMz(spectra, si)
        start = length(new_precursors) + 1
        if ismissing(cv) || ismissing(wv)
            for r in rng
                push!(new_precursors, precursors_passed[r])
            end
        else
            c = Float32(cv); hw = Float32(wv) / 2
            for r in rng
                p = precursors_passed[r]
                abs(prec_mzs[p] - c) <= hw && push!(new_precursors, p)
            end
        end
        scan_to_prec_idx[si] = length(new_precursors) >= start ?
            (start:length(new_precursors)) : missing
    end
    return new_precursors
end

# Median MS2 isolation width (≈ Q1 bin step S) — sizes the ZT deconv box.
function _zt_bin_step(spectra::MassSpecData, all_scan_idxs::Vector{Int})
    ws = Float32[]
    @inbounds for si in all_scan_idxs
        w = getIsolationWidthMz(spectra, si)
        ismissing(w) && continue
        push!(ws, Float32(w))
        length(ws) >= 5000 && break
    end
    isempty(ws) && return 1.0f0
    sort!(ws)
    return ws[cld(length(ws), 2)]
end

"""
    expand_to_metascans!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, k)

Sciex ZT scanning-DIA meta-scan expansion. The scanning quadrupole spreads a
precursor's ions across ~`2k+1` adjacent Q1 bins (consecutive MS2 scans within a
cycle), peaking at the center bin. The fragment index only matched a precursor to
its center bin (±half-bin tolerance); this replaces each scan's candidate set with
the UNION of the candidates of all scans within ±`k` of it in the SAME cycle, so
downstream deconvolution estimates a per-bin weight for the precursor across its
whole meta-scan (which the triangular-reduction step then exploits).

Bins are consecutive scan indices within a cycle, so neighbors are `si±j` guarded
by `getMsOrder == 2` and equal `getCycleIdx` (never crosses a cycle/MS1 boundary).
Rebuilds `precursors_passed` and reindexes `scan_to_prec_idx` in place (mirrors
`merge_precursors_by_window!`). Returns the new `precursors_passed`.
"""
function expand_to_metascans!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    k::Int,
)
    (k <= 0 || isempty(all_scan_idxs)) && return precursors_passed
    n = length(all_scan_idxs)
    nspec = length(spectra)

    # Pass 1: per searched scan, union candidates from neighbors si±j (same cycle,
    # MS2). Reads the OLD scan_to_prec_idx / precursors_passed; no mutation yet.
    new_sets = Vector{Set{UInt32}}(undef, n)
    @inbounds for i in 1:n
        si = all_scan_idxs[i]
        ci = getCycleIdx(spectra, si)
        s = Set{UInt32}()
        for j in -k:k
            sj = si + j
            (sj < 1 || sj > nspec) && continue
            isassigned(scan_to_prec_idx, sj) || continue
            (getMsOrder(spectra, sj) == 2) || continue
            (getCycleIdx(spectra, sj) == ci) || continue
            rng = scan_to_prec_idx[sj]
            ismissing(rng) && continue
            for r in rng
                push!(s, precursors_passed[r])
            end
        end
        new_sets[i] = s
    end

    # Pass 2: build the new flat candidate vector and reindex scan_to_prec_idx.
    new_precursors = UInt32[]
    sizehint!(new_precursors, sum(length, new_sets; init = 0))
    @inbounds for i in 1:n
        si = all_scan_idxs[i]
        s = new_sets[i]
        if isempty(s)
            scan_to_prec_idx[si] = missing
        else
            start = length(new_precursors) + 1
            append!(new_precursors, s)
            scan_to_prec_idx[si] = start:length(new_precursors)
        end
    end
    return new_precursors
end

"""
    merge_precursors_by_window!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs)

Group consecutive scans with the same isolation window (center_mz ± 0.1 Da,
isolation_width ± 0.1 Da) and replace each scan's precursor range with the
union of all precursor IDs in the window group.

This allows `process_scans!` to build the transition list once per window
and reuse it for all scans in the group. Returns a new `precursors_passed`
vector with the merged ranges.
"""
function merge_precursors_by_window!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int};
    max_group_size::Int = 3
)
    n_scans = length(all_scan_idxs)
    n_scans == 0 && return precursors_passed

    # Pass 1: identify window groups (consecutive scans with same isolation window)
    # group_id[i] = which group scan i belongs to
    # Groups are capped at max_group_size to limit precursor list inflation.
    group_id = Vector{Int}(undef, n_scans)
    current_group = 1
    group_size = 1
    group_id[1] = 1
    prev_center = getCenterMz(spectra, all_scan_idxs[1])
    prev_width = getIsolationWidthMz(spectra, all_scan_idxs[1])

    for i in 2:n_scans
        scan_idx = all_scan_idxs[i]
        center = getCenterMz(spectra, scan_idx)
        width = getIsolationWidthMz(spectra, scan_idx)
        if abs(center - prev_center) > 0.1f0 || abs(width - prev_width) > 0.1f0 ||
           group_size >= max_group_size
            current_group += 1
            group_size = 0
            prev_center = center
            prev_width = width
        end
        group_size += 1
        group_id[i] = current_group
    end
    n_groups = current_group

    # Pass 2: for each group, collect the union of precursor IDs (deduplicated)
    # Use a Set per group to deduplicate
    group_prec_sets = [Set{UInt32}() for _ in 1:n_groups]
    for i in 1:n_scans
        scan_idx = all_scan_idxs[i]
        range = scan_to_prec_idx[scan_idx]
        ismissing(range) && continue
        gid = group_id[i]
        for j in range
            push!(group_prec_sets[gid], precursors_passed[j])
        end
    end

    # Pass 3: build new flat precursors vector and update scan_to_prec_idx
    new_precursors = UInt32[]
    sizehint!(new_precursors, length(precursors_passed))
    group_ranges = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_groups)

    for gid in 1:n_groups
        precs = group_prec_sets[gid]
        if isempty(precs)
            group_ranges[gid] = missing
        else
            start = length(new_precursors) + 1
            append!(new_precursors, precs)
            group_ranges[gid] = start:length(new_precursors)
        end
    end

    # Capture original per-scan counts BEFORE updating scan_to_prec_idx
    orig_counts = Dict{Int, Int}()  # scan_idx → original precursor count
    total_original = 0
    total_merged = 0
    n_active_scans = 0
    for i in 1:n_scans
        scan_idx = all_scan_idxs[i]
        orig_range = scan_to_prec_idx[scan_idx]
        merged_range = group_ranges[group_id[i]]
        if !ismissing(orig_range)
            orig_counts[scan_idx] = length(orig_range)
            total_original += length(orig_range)
            n_active_scans += 1
        else
            orig_counts[scan_idx] = 0
        end
        if !ismissing(merged_range)
            total_merged += length(merged_range)
        end
    end

    # Now update scan_to_prec_idx
    for i in 1:n_scans
        scan_idx = all_scan_idxs[i]
        scan_to_prec_idx[scan_idx] = group_ranges[group_id[i]]
    end

    if n_active_scans > 0
        avg_orig = round(total_original / n_active_scans, digits=1)
        avg_merged = round(total_merged / n_active_scans, digits=1)
        inflation = round(100 * (total_merged - total_original) / max(total_original, 1), digits=1)

        # Detailed: sample 10 groups from data-dense region (~15 min, ~600 m/z)
        # Find groups near RT=15 min and mz=600
        target_rt = 15.0f0
        target_mz = 600.0f0
        best_gid = 1
        best_dist = Inf
        for gid in 1:n_groups
            members = [i for i in 1:n_scans if group_id[i] == gid]
            isempty(members) && continue
            s = all_scan_idxs[first(members)]
            rt = Float64(getRetentionTime(spectra, s))
            mz = Float64(getCenterMz(spectra, s))
            dist = abs(rt - target_rt) + abs(mz - target_mz) * 0.01
            if dist < best_dist
                best_dist = dist
                best_gid = gid
            end
        end
        sample_start = max(1, best_gid - 5)
        sample_end = min(n_groups, sample_start + 9)
        group_details = String[]
        for gid in sample_start:sample_end
            members = [i for i in 1:n_scans if group_id[i] == gid]
            if !isempty(members)
                scan_idxs = [all_scan_idxs[m] for m in members]
                rts = [Float64(getRetentionTime(spectra, s)) for s in scan_idxs]
                mzs = [Float64(getCenterMz(spectra, s)) for s in scan_idxs]
                per_scan_counts = [get(orig_counts, s, 0) for s in scan_idxs]
                union_count = ismissing(group_ranges[gid]) ? 0 : length(group_ranges[gid])
                push!(group_details, "    g$gid: $(length(members)) scans, " *
                    "RT=$(round.(rts, digits=2)), mz=$(round.(mzs, digits=1)), " *
                    "per_scan=$(per_scan_counts), union=$union_count")
            end
        end

        @debug_l1 "  window merge: $n_groups groups from $n_scans scans (max_group=$max_group_size), " *
                   "avg precs/scan: $avg_orig → $avg_merged (+$(inflation)%)\n"
        for d in group_details
            @debug_l1 "$d\n"
        end
    end

    return new_precursors
end

"""
    write_fragment_index_matches(scan_to_prec_idx, precursors_passed, output_path)

Flatten (scan_idx, precursor_idx) pairs from fragment index search and write to Arrow.
"""
function write_fragment_index_matches(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    output_path::String
)
    scan_idxs = UInt32[]
    prec_idxs = UInt32[]
    for (scan_idx, range) in enumerate(scan_to_prec_idx)
        ismissing(range) && continue
        for i in range
            push!(scan_idxs, UInt32(scan_idx))
            push!(prec_idxs, precursors_passed[i])
        end
    end
    Arrow.write(output_path, (scan_idx=scan_idxs, precursor_idx=prec_idxs))
end

function getRTWindow(irt::U, irt_tol::T) where {T,U<:AbstractFloat}
    return Float32(irt - irt_tol), Float32(irt + irt_tol)
end
