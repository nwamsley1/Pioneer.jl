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

# An empty PSM frame carrying the columns the search would have produced.
#
# `process_scans_fused!` builds its result from the per-thread scored-PSM
# struct-of-arrays, so a zero-length view of that same store gives the exact
# schema with no rows — and cannot drift from it the way a hand-written column
# list would.
#
# The early exits in `library_search` returned a bare `DataFrame()` instead, so a
# file that matched nothing produced a frame with *no columns at all*, and the
# first downstream reader of `:precursor_idx` threw `ArgumentError: column name
# :precursor_idx not found` rather than seeing zero PSMs. Every such reader
# already handles zero rows; none of them could handle zero columns.
_empty_scored_psms(search_data, params) =
    DataFrame(@view(get_scored_psms(first(search_data), params)[1:0]))

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
    zt_chunk_candidates::Int = 0,
    zt_reduce = nothing,
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
    # Scanning-quad (ZT) two-width quad model. The recorded isolation width is only the Q1
    # step; the physical window is several Da swept across m/z. So the fragment index keeps a
    # NARROW ~1 m/z box (candidacy stays tight and specific) while deconvolution uses a box
    # wide enough to span the whole meta-scan. `qtm_deconv` is finalized once k is known below.
    zt_geom = getZTGeometry(search_context, ms_file_idx)
    zt_on   = zt_geom !== nothing
    zt_k    = zt_on ? Int(zt_geom.metascan_k) : 0
    # Candidacy expansion and the wide deconv box are MAIN-search only: the tuning searches
    # must not be calibrated on the metascan-expanded, wide-box deconvolution.
    zt_main = zt_on && (params isa MainSearchParameters)
    # Scanning-quad (ZT) quad tuning needs the SAME meta-scan view as the main search: the
    # expansion, so a precursor is seen in all 2k+1 bins, and the WIDE square deconvolution box,
    # so its fitted weight tracks true transmission instead of being divided by an assumed
    # model. Without both, every (precursor, cycle) group holds 1-2 bins and the triangle fit
    # has nothing to regress against. Only the fit differs from the main search, not the view.
    zt_qtune = zt_on && (params isa QuadTuningSearchParameters)
    zt_meta  = zt_main || zt_qtune
    qtm_frag = zt_on ? SquareQuadModel(zt_candidacy_overhang(zt_geom)) : qtm

    # DIAGNOSTIC (PIONEER_ZT_QUAD_PROBE=<Da>): widen the quad-tuning candidacy AND deconv box
    # so the isotope-pair measurement can observe offsets far enough off bin-center to resolve
    # the swept transmission. QuadTuning measures yt = log T(x0) - log T(x1) across the isotope
    # spacing, so mapping a several-Da-wide profile needs x0 to span several Da — which it
    # cannot if candidacy stops at the recorded ~1 Da step. Inert unless the env var is set.
    zt_probe = (zt_on && params isa QuadTuningSearchParameters) ?
        something(tryparse(Float32, get(ENV, "PIONEER_ZT_QUAD_PROBE", "")), 0f0) : 0f0
    if zt_probe > 0f0
        qtm_frag = SquareQuadModel(zt_probe)
    end

    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt = getRtIrtModel(search_context, ms_file_idx)
    precursors = getPrecursors(spec_lib)
    ion_list = getFragmentLookupTable(spec_lib)

    # NCE models to iterate: [(calibrated_model, nothing)]
    nce_entries = get_nce_models(search_context, params, ms_file_idx)
    isempty(nce_entries) && return _empty_scored_psms(search_data, params)

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
        thread_tasks = _round_robin_tasks(all_scan_idxs, Threads.nthreads())
    end

    isempty(all_scan_idxs) && return _empty_scored_psms(search_data, params)

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
    # DEVELOPER PROFILING (disabled by default; deliberately not behind a runtime flag).
    #
    # PProf is NOT a dependency of Pioneer. It used to be, purely to support the block below --
    # which is guarded by a `const false`, so it could never execute in a release build, yet the
    # package was still compiled into every shipped binary. It also pulled FlameGraphs, which caps
    # FixedPointNumbers below 0.9 and makes local `incremental=false` builds fragile.
    #
    # To profile the fragment index search: uncomment this block, add `using Profile, PProf` to
    # src/Pioneer.jl, and run from the `dev` environment (which has PProf). Do not commit either.
    #
    #     Profile.clear()
    #     Profile.init(n = 50_000_000, delay = 0.0005)
    #     Profile.@profile precursors_passed, scores_passed = searchFragmentIndexPartitionMajorHinted(
    #         scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
    #         Threads.nthreads(), params, qtm, mem, rt_to_irt, irt_tol,
    #         getMz(precursors);
    #         score_filter = score_filter, max_peaks = max_peaks,
    #         scratch = getFragIndexScratch(search_context))
    #     prof_path = joinpath(getDataOutDir(search_context),
    #                          "frag_index_profile_$(ms_file_idx).pb.gz")
    #     pprof(out = prof_path, web = false)
    #     @user_info "Fragment index profile saved to $prof_path\n"
    #
    precursors_passed, scores_passed = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        Threads.nthreads(), params, qtm_frag, mem, rt_to_irt, irt_tol,
        getMz(precursors);
        score_filter = score_filter, max_peaks = max_peaks,
        scratch = getFragIndexScratch(search_context))
    t_frag = time() - t_frag_start

    # --- DEBUG: dump fragment index bitmask scores to Arrow and bail ---
    # Only dump during MainSearch, not tuning stages.
    # Applies RT and precursor m/z filtering (same as selectTransitions!) for realistic counts.
    # --- 2a. Pre-filter: require candidates to appear in ≥ N scans ---
    prefilter_n = getPrefilterMinScanCount(params)
    if prefilter_n > 1
        n_before = length(precursors_passed)
        precursors_passed = filter_low_scan_candidates!(
            scan_to_prec_idx, precursors_passed, prefilter_n)
        @debug_l1 "Pre-filter: $n_before → $(length(precursors_passed)) candidates " *
              "($(round(100*(1 - length(precursors_passed)/max(1,n_before)), digits=1))% removed)"
    end

    # --- 2b. ZT: anchor candidacy on the precursor's own bin, then span the meta-scan ---
    # Non-ZT files and the tuning searches keep the tuned quad model untouched.
    # Deconvolution uses the file's installed transmission model. On ZT that is the flat box
    # spanning the meta-scan, set once in ensure_zt_geometry! — so there is no per-call-site
    # patching here, and every other consumer sees the same model.
    # Tuning searches must NOT calibrate on the wide meta-scan deconvolution box: it admits
    # many off-center precursors whose interference degrades the mass-error and NCE fits. Only
    # the MAIN search deconvolves across the meta-scan; everything else stays on the bin.
    qtm_deconv = (zt_on && !zt_meta) ? SquareQuadModel(0.0f0) : qtm
    zt_probe > 0f0 && (qtm_deconv = SquareQuadModel(zt_probe))
    if zt_meta
        n_emitted = length(precursors_passed)
        # Wide-emit re-anchors every emission to the precursor's own bin (survives if it cleared
        # in ANY bin); the narrow path requires the center bin itself to clear.
        # Re-anchoring is useful even at the DEFAULT box: with isotope_err_bounds (1,0) the
        # effective search window is ~1.52 Da against a ~1.02 Da bin step, so a precursor in the
        # upper half of its bin is also a candidate in the bin above. filter_to_center_bin!
        # discards that emission; if it cleared the bitvec there but not in its own bin, the
        # precursor is lost. Re-anchoring keeps it at no extra emission cost.
        _reanchor = zt_candidacy_tol() > zt_geom.nominal_width / 2 ||
                    get(ENV, "PIONEER_ZT_REANCHOR", "0") != "0"
        precursors_passed = if _reanchor
            map_any_hit_to_center!(scan_to_prec_idx, precursors_passed, spectra,
                                   all_scan_idxs, getMz(precursors), zt_geom)
        else
            filter_to_center_bin!(scan_to_prec_idx, precursors_passed, spectra,
                                  all_scan_idxs, getMz(precursors))
        end
        n_center = length(precursors_passed)
        if zt_k > 0
            precursors_passed = expand_to_metascans!(
                scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, zt_k)
        end
        # Report the boxes by INTERROGATING the models actually in use, never by recomputing
        # what they were meant to be — a log that restates intent cannot catch a model that
        # something else overwrote.
        _c = Float32(500)
        _fq = getQuadTransmissionFunction(qtm_frag,   _c, zt_geom.nominal_width)
        _dq = getQuadTransmissionFunction(qtm_deconv, _c, zt_geom.nominal_width)
        _fw = (getPrecMaxBound(_fq) - getPrecMinBound(_fq)) / 2
        _dw = (getPrecMaxBound(_dq) - getPrecMinBound(_dq)) / 2
        @debug_l1 "ZT candidacy (k=$zt_k): $n_emitted emitted -> $n_center center -> " *
                   "$(length(precursors_passed)) expanded; candidacy box +/-$(round(_fw, digits=2)) Da, " *
                   "deconv box +/-$(round(_dw, digits=2)) Da (expansion span +/-" *
                   "$(round(Float32(zt_k) * zt_geom.bin_step, digits=2)) Da)"
        _dw < Float32(zt_k) * zt_geom.bin_step &&
            @user_warn "ZT: deconv box (+/-$(round(_dw,digits=2)) Da) is NARROWER than the " *
                       "expansion span (+/-$(round(Float32(zt_k)*zt_geom.bin_step,digits=2)) Da) — " *
                       "expanded candidates in outer bins will get zero transmission." 
    end

    # DIAGNOSTIC (PIONEER_ZT_OVERLAP=1, main search): candidate-set overlap between consecutive
    # bins of a cycle after expansion — the template-reuse opportunity. Candidates are sorted
    # within each scan, so overlap is a linear merge.
    if zt_meta && zt_k > 0 && get(ENV, "PIONEER_ZT_OVERLAP", "0") != "0" && params isa MainSearchParameters
        _zt_report_adjacent_overlap(scan_to_prec_idx, precursors_passed, spectra)
    end

    if zt_meta && zt_k > 0 && zt_chunk_candidates > 0 && zt_reduce !== nothing
        # EXPERIMENT (PIONEER_ZT_EVEN_BINS=1): after expansion, deconvolve only every second MS2
        # scan of each cycle (even position within the cycle). Halves the solves; the collapse
        # still sees ~half the bins of every meta-scan. Candidate sets are untouched.
        # PIONEER_ZT_EVEN_BINS=1 drops every odd-position scan. PIONEER_ZT_EVEN_BINS=core:<c>
        # keeps every scan whose cycle position is within c bins of ANY precursor's own bin —
        # approximated here by keeping the scan if it is within c of a scan that a candidate
        # anchors on. Simpler proxy used: keep odd scans too when (pos mod 2k+1) is within c of
        # the meta-scan centre. Since meta-scans are not aligned to cycle positions, the
        # practical variant is: drop odd scans only where the bin is >= c from the precursor
        # m/z for ALL its candidates — too costly to evaluate here. So "core:<c>" instead thins
        # the OUTER bins by candidate: for each dropped odd scan, candidates whose precursor m/z
        # lies within c bins of that scan's centre are moved back in (kept). See _zt_thin_bins!.
        # Default on scanning-quad files: thin the outer bins (core ±2). Measured on 5 Da A_REP1:
        # −20% main search for −0.7% precursors (29,635 -> 29,420); the outer bins of a
        # meta-scan carry <40% transmission and every second one is enough for the collapse.
        # PIONEER_ZT_EVEN_BINS=0 restores all bins; =1 or core:<c> select other variants.
        _mode = get(ENV, "PIONEER_ZT_EVEN_BINS", "core:2")
        begin
            if _mode == "0"
                # all bins
            elseif _mode == "1"
                _dropped = 0
                for r in zt_cycle_scan_ranges(spectra), (pos, si) in enumerate(r)
                    if isodd(pos) && !ismissing(scan_to_prec_idx[si])
                        scan_to_prec_idx[si] = missing; _dropped += 1
                    end
                end
                @user_info "ZT even-bins experiment: dropped $_dropped odd-position MS2 scans from deconvolution"
            elseif startswith(_mode, "core:")
                _c = parse(Int, _mode[6:end])
                precursors_passed = _zt_thin_outer_bins!(scan_to_prec_idx, precursors_passed, spectra,
                                                         getMz(precursors), zt_geom, _c)
            end
        end
    end

    prec_index = PerScanPrecursorIndex(scan_to_prec_idx, precursors_passed)

    # --- 3. Threaded scan processing, once per NCE model ---
    # All search methods use the fused per-precursor scan loop. The classic
    # `process_scans!` path was deleted along with the multi-pass pipeline
    # (selectTransitions! + matchPeaks! + buildDesignMatrix! + sortSparse!).
    # When nce_tag is not nothing (NCE tuning), tag each result with the NCE value.
    t_deconv_start = time()
    if params isa MainSearchParameters
        DECONV_SOLVES[] = 0; DECONV_ITERS[] = 0
        PMM_STATS_ON[] = haskey(ENV, "PIONEER_PMM_STATS")
        PMM_STAT_VISITS[] = 0; PMM_STAT_ZERO_VISITS[] = 0; PMM_STAT_NNZ[] = 0; PMM_STAT_COLS[] = 0; PMM_STAT_ZERO_END[] = 0
    end
    # DIAGNOSTIC (PIONEER_PROFILE_DECONV=1, main search only): sample the threaded deconv with
    # Julia's Profile and write a flat self-time profile to the output dir, so the run_fused!
    # match/design-matrix build can be separated from the solver. Adds sampling overhead when on.
    _prof_deconv = (params isa MainSearchParameters) && get(ENV, "PIONEER_PROFILE_DECONV", "0") != "0"
    if _prof_deconv
        Profile.clear()
        Profile.init(n = 200_000_000, delay = 0.001)
    end
    _deconv_body = (tt) -> map(nce_entries) do (nce_model, nce_tag)
        intensity_model = prepare_fragment_intensity_model(ion_list, nce_model)
        tasks = map(tt) do thread_task
            Threads.@spawn process_scans_fused!(
                last(thread_task), spectra, prec_index,
                ms_file_idx,
                search_data[first(thread_task)], params, precursors, ion_list,
                intensity_model, qtm_deconv, mem, rt_to_irt, irt_tol)
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
    # --- 3b. Scanning-quad chunked deconvolution (MainSearch only, zt_chunk_candidates > 0) ---
    # The fragment index + expansion above ran ONCE over the whole file, so every scan's
    # candidate count is exact. Cycles are grouped into chunks of ~zt_chunk_candidates
    # candidates (cycle-aligned: expansion and the meta-scan collapse never cross a cycle),
    # each chunk is deconvolved with segment-major thread tasks and handed to `zt_reduce`,
    # which returns the reduced (collapsed) table. Only the reduced tables accumulate, so the
    # raw per-bin rows of one chunk are the peak, never the whole file's.
    if zt_meta && zt_k > 0 && zt_chunk_candidates > 0 && zt_reduce !== nothing
        chunks = zt_cycle_chunks_by_candidates(zt_cycle_scan_ranges(spectra), scan_to_prec_idx,
                                               zt_chunk_candidates)
        reduced = DataFrame(); n_raw_total = 0; max_raw = 0
        for (ci, chunk) in enumerate(chunks)
            t_c = time()
            n_cand = zt_candidate_count(chunk, scan_to_prec_idx)
            tt = zt_thread_tasks(chunk, zt_k, Threads.nthreads())
            raw_all = _prof_deconv ? (Profile.@profile _deconv_body(tt)) : _deconv_body(tt)
            raw = length(raw_all) == 1 ? raw_all[1] : vcat(raw_all...)
            t_d = time() - t_c
            n_raw = nrow(raw); n_raw_total += n_raw; max_raw = max(max_raw, n_raw)
            part = zt_reduce(raw, ci)
            raw = nothing
            reduced = isempty(reduced) ? part : (append!(reduced, part); reduced)
            @user_info "ZT chunk $ci/$(length(chunks)): $(length(chunk)) cycles, " *
                       "$(sum(length, chunk)) scans, $n_cand candidates -> $n_raw raw rows " *
                       "(deconv $(round(t_d; digits=1))s) -> $(nrow(part)) reduced " *
                       "(reduce $(round(time() - t_c - t_d; digits=1))s); cumulative $(nrow(reduced))"
        end
        t_deconv = time() - t_deconv_start
        if _prof_deconv
            _pp = joinpath(getDataOutDir(search_context), "deconv_profile_$(ms_file_idx).txt")
            open(_pp, "w") do io
                Profile.print(IOContext(io, :displaysize => (100000, 320));
                              format = :flat, sortedby = :count, mincount = 20)
            end
            @user_info "Deconv flat profile (all chunks) saved to $_pp"
        end
        @user_info "ZT chunked main search: $(length(chunks)) chunks, largest $max_raw raw rows, " *
                   "$n_raw_total raw -> $(nrow(reduced)) reduced; frag_index=$(round(t_frag, digits=1))s " *
                   "deconv+reduce=$(round(t_deconv, digits=1))s candidates=$(length(precursors_passed)); " *
                   "solver: $(DECONV_SOLVES[]) solves, mean $(round(DECONV_ITERS[] / max(DECONV_SOLVES[], 1); digits=2)) iters"
        if PMM_STATS_ON[]
            @user_info "PMM stats: $(PMM_STAT_COLS[]) columns over $(DECONV_SOLVES[]) solves " *
                       "(mean $(round(PMM_STAT_COLS[] / max(DECONV_SOLVES[],1); digits=0))/solve), " *
                       "$(round(100 * PMM_STAT_ZERO_END[] / max(PMM_STAT_COLS[],1); digits=1))% at zero at exit; " *
                       "$(PMM_STAT_VISITS[]) column-visits, $(round(100 * PMM_STAT_ZERO_VISITS[] / max(PMM_STAT_VISITS[],1); digits=1))% zero->zero; " *
                       "mean nnz/column $(round(PMM_STAT_NNZ[] / max(PMM_STAT_VISITS[],1); digits=1))"
        end
        return reduced
    end

    all_results = _prof_deconv ? (Profile.@profile _deconv_body(thread_tasks)) : _deconv_body(thread_tasks)
    if _prof_deconv
        _pp = joinpath(getDataOutDir(search_context), "deconv_profile_$(ms_file_idx).txt")
        open(_pp, "w") do io
            Profile.print(IOContext(io, :displaysize => (100000, 320));
                          format = :flat, sortedby = :count, mincount = 20)
        end
        @user_info "Deconv flat profile saved to $_pp"
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
        @user_info "library_search breakdown: frag_index=$(round(t_frag, digits=2))s  " *
                   "deconv=$(round(t_deconv, digits=2))s  " *
                   "vcat=$(round(t_vcat, digits=2))s  " *
                   "candidates=$(length(get_precursors(prec_index))); " *
                   "solver: $(DECONV_SOLVES[]) solves, mean $(round(DECONV_ITERS[] / max(DECONV_SOLVES[], 1); digits=2)) iters"
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
    filter_to_center_bin!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, prec_mzs)

Scanning-quad (ZT) center-bin candidacy. The fragment index runs with a widened box so a
precursor can be emitted from neighbouring Q1 bins; this keeps only the emissions whose scan is
the precursor's own bin, i.e. `|prec_mz - centerMz| <= isolationWidth/2`. `expand_to_metascans!`
then refills the ±k neighbours, so the meta-scan is anchored on the precursor's true bin rather
than on wherever it happened to be emitted.

Rebuilds `precursors_passed` and reindexes `scan_to_prec_idx` in place (mirrors
`filter_low_scan_candidates!`). Returns the new `precursors_passed`.
"""
function filter_to_center_bin!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    prec_mzs::AbstractVector{Float32},
)
    new_passed = UInt32[]
    sizehint!(new_passed, length(precursors_passed))
    @inbounds for si in all_scan_idxs
        rng = scan_to_prec_idx[si]
        ismissing(rng) && continue
        cv = getCenterMz(spectra, si)
        wv = getIsolationWidthMz(spectra, si)
        start = length(new_passed) + 1
        if ismissing(cv) || ismissing(wv)
            # No window metadata: keep everything rather than silently dropping candidates.
            for r in rng
                push!(new_passed, precursors_passed[r])
            end
        else
            c = Float32(cv); hw = Float32(wv) / 2
            for r in rng
                p = precursors_passed[r]
                abs(prec_mzs[p] - c) <= hw && push!(new_passed, p)
            end
        end
        scan_to_prec_idx[si] = length(new_passed) >= start ?
            (start:length(new_passed)) : missing
    end
    return new_passed
end

"""
    map_any_hit_to_center!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs,
                           prec_mzs, geom) -> Vector{UInt32}

Scanning-quad (ZT) wide-emit candidacy. The fragment index runs with a WIDENED box, so a
precursor can be emitted from any bin of its meta-scan; this maps every emission back to the
precursor's OWN bin, deduped per (cycle, precursor). Net effect: a precursor survives if it
cleared the bitvec in ANY of its bins, where `filter_to_center_bin!` requires the center bin
itself to clear. Different bins expose different fragment subsets, so those are real second
looks rather than noise admission.

`expand_to_metascans!` then fills each survivor's +/-k as usual, so collapse and the shape
features are unchanged.

Implementation notes — the reference version was the single fattest serial step in the search
(~181 s over ~60 M emissions), for three separable reasons, all addressed here:

  * It found each precursor's bin by scanning `+/-search_halfbins` neighbours with three
    accessor calls apiece. The lattice is uniform to Float32 granularity, so the bin is
    `si + round((prec_mz - centerMz[si]) / S)` in O(1).
  * It deduped with a `Set{Tuple{UInt32,UInt32}}` and accumulated into a
    `Dict{Int,Set{UInt32}}`. Center scan is bounded by `length(spectra)`, so this is a COUNTING
    sort: count per center scan, prefix-sum, scatter by index. No hashing, no `push!` in the
    hot loop, one exactly-sized allocation, and dedup becomes a sort of each scan's ~160-element
    slice rather than one global sort of tens of millions.
  * It was serial. Both passes and the per-scan dedup are threaded.

`cs` is deliberately recomputed in the scatter pass rather than stored: the arithmetic is a
multiply, a round and a clamp, against hundreds of MB of stores and loads.
"""
function map_any_hit_to_center!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    prec_mzs::AbstractVector{Float32},
    geom::ZTGeometry,
)
    (isempty(all_scan_idxs) || isempty(precursors_passed)) && return precursors_passed
    nspec = length(spectra)
    inv_S = 1.0f0 / geom.bin_step

    # Per-cycle MS2 scan bounds, so a re-anchored center never crosses a ramp boundary.
    cyc = Vector{UInt32}(getCycleIdxs(spectra))
    ncyc = 0
    @inbounds for si in 1:nspec
        getMsOrder(spectra, si) == 2 || continue
        c = Int(cyc[si]); c > ncyc && (ncyc = c)
    end
    ncyc == 0 && return precursors_passed
    cyc_lo = fill(typemax(Int), ncyc)
    cyc_hi = zeros(Int, ncyc)
    @inbounds for si in 1:nspec
        getMsOrder(spectra, si) == 2 || continue
        c = Int(cyc[si])
        si < cyc_lo[c] && (cyc_lo[c] = si)
        si > cyc_hi[c] && (cyc_hi[c] = si)
    end

    cmzs = Float32.(coalesce.(getCenterMzs(spectra), NaN32))

    n = length(all_scan_idxs)
    nchunks = min(Threads.nthreads(), n)
    bounds = [(n * (t - 1)) ÷ nchunks + 1 for t in 1:(nchunks + 1)]
    bounds[end] = n + 1

    # --- Pass A: count emissions per center scan (parallel; nothing allocates in the loop) ---
    counts = [zeros(Int64, nspec) for _ in 1:nchunks]
    Threads.@threads for t in 1:nchunks
        ct = counts[t]
        @inbounds for i in bounds[t]:(bounds[t + 1] - 1)
            si = all_scan_idxs[i]
            rng = scan_to_prec_idx[si]
            ismissing(rng) && continue
            cm = cmzs[si]; isnan(cm) && continue
            c = Int(cyc[si]); (c < 1 || c > ncyc) && continue
            lo, hi = cyc_lo[c], cyc_hi[c]
            for r in rng
                cs = si + round(Int, (prec_mzs[precursors_passed[r]] - cm) * inv_S)
                ct[clamp(cs, lo, hi)] += 1
            end
        end
    end

    # --- Prefix sums: each chunk gets its own write cursor per scan, so scatters never race ---
    scan_start = Vector{Int}(undef, nspec)
    scan_len   = zeros(Int32, nspec)
    total = 0
    @inbounds for cs in 1:nspec
        scan_start[cs] = total
        s = 0
        for t in 1:nchunks
            ct = counts[t][cs]
            counts[t][cs] = total + s      # this chunk's start slot for this scan
            s += ct
        end
        scan_len[cs] = Int32(s)
        total += s
    end
    total == 0 && return UInt32[]

    # --- Pass B: scatter by index (parallel; no push!, one exactly-sized allocation) ---
    scratch = Vector{UInt32}(undef, total)
    Threads.@threads for t in 1:nchunks
        pos = counts[t]
        @inbounds for i in bounds[t]:(bounds[t + 1] - 1)
            si = all_scan_idxs[i]
            rng = scan_to_prec_idx[si]
            ismissing(rng) && continue
            cm = cmzs[si]; isnan(cm) && continue
            c = Int(cyc[si]); (c < 1 || c > ncyc) && continue
            lo, hi = cyc_lo[c], cyc_hi[c]
            for r in rng
                p = precursors_passed[r]
                cs = clamp(si + round(Int, (prec_mzs[p] - cm) * inv_S), lo, hi)
                k = pos[cs] + 1
                scratch[k] = p
                pos[cs] = k
            end
        end
    end

    # --- Dedup within each scan's slice (parallel; ~160 elements each, not one global sort) ---
    new_len = zeros(Int32, nspec)
    Threads.@threads for cs in 1:nspec
        L = Int(scan_len[cs]); L == 0 && continue
        st = scan_start[cs]
        v = view(scratch, (st + 1):(st + L))
        sort!(v)
        m = 1
        @inbounds for i in 2:L
            if v[i] != v[m]
                m += 1; v[m] = v[i]
            end
        end
        new_len[cs] = Int32(m)
    end

    # --- Compact and reindex ---
    outn = 0
    @inbounds for cs in 1:nspec
        outn += new_len[cs]
    end
    new_passed = Vector{UInt32}(undef, outn)
    @inbounds for si in 1:nspec
        scan_to_prec_idx[si] = missing
    end
    off = 0
    @inbounds for cs in 1:nspec
        L = Int(new_len[cs]); L == 0 && continue
        copyto!(new_passed, off + 1, scratch, scan_start[cs] + 1, L)
        scan_to_prec_idx[cs] = (off + 1):(off + L)
        off += L
    end
    return new_passed
end

"""
    _sort_dedup!(v) -> Int

Sort `v` and compact duplicates to the front in place. Returns the count of unique elements;
`v[1:n]` holds them. Replaces a per-scan `Set{UInt32}` in `expand_to_metascans!`.
"""
@inline function _sort_dedup!(v::Vector{UInt32})
    isempty(v) && return 0
    sort!(v)
    m = 1
    @inbounds for i in 2:length(v)
        if v[i] != v[m]
            m += 1
            v[m] = v[i]
        end
    end
    return m
end

"""
    expand_to_metascans!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs, k)

Scanning-quad (ZT) meta-scan expansion. The swept quadrupole spreads a precursor's ions across
~`2k+1` adjacent Q1 bins (consecutive MS2 scans within a cycle), so replace each searched scan's
candidate set with the UNION of the candidates of every scan within ±`k` of it in the SAME cycle.
Deconvolution then estimates a per-bin weight across the whole meta-scan, which the collapse step
exploits.

Neighbours are `si±j` guarded by `getMsOrder == 2` and equal `getCycleIdx`, so expansion never
crosses a cycle or MS1 boundary.

Rebuilds `precursors_passed` and reindexes `scan_to_prec_idx` in place (mirrors
`filter_low_scan_candidates!`). Returns the new `precursors_passed`. Candidates are sorted within
each scan.

Threaded over contiguous chunks of `all_scan_idxs`. **The parallel phase is strictly read-only on
`scan_to_prec_idx`; every write happens in the serial phase after it.** A scan's neighbours must
be read as they were BEFORE any expansion, so the two phases cannot be fused.
"""
function expand_to_metascans!(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    spectra::MassSpecData,
    all_scan_idxs::Vector{Int},
    k::Int,
)
    (k <= 0 || isempty(all_scan_idxs)) && return precursors_passed
    n     = length(all_scan_idxs)
    nspec = length(spectra)

    # Contiguous chunks, assigned in order, so concatenating buffers in chunk order
    # reproduces scan order.
    nchunks = min(Threads.nthreads(), n)
    bounds  = [(n * (t - 1)) ÷ nchunks + 1 for t in 1:(nchunks + 1)]
    bounds[end] = n + 1

    counts   = Vector{Int32}(undef, n)
    out_bufs = Vector{Vector{UInt32}}(undef, nchunks)

    # --- Phase 1: gather + dedup (parallel, READ-ONLY on scan_to_prec_idx) ---
    Threads.@threads for t in 1:nchunks
        lo, hi  = bounds[t], bounds[t + 1] - 1
        buf     = UInt32[]
        scratch = UInt32[]
        @inbounds for i in lo:hi
            si = all_scan_idxs[i]
            ci = getCycleIdx(spectra, si)
            empty!(scratch)
            for j in -k:k
                sj = si + j
                (sj < 1 || sj > nspec) && continue
                getMsOrder(spectra, sj) == 2 || continue
                getCycleIdx(spectra, sj) == ci || continue
                rng = scan_to_prec_idx[sj]
                ismissing(rng) && continue
                for r in rng
                    push!(scratch, precursors_passed[r])
                end
            end
            m = _sort_dedup!(scratch)
            counts[i] = Int32(m)
            m > 0 && append!(buf, @view scratch[1:m])
        end
        out_bufs[t] = buf
    end

    # --- Phase 2: concatenate and reindex (serial) ---
    total = 0
    @inbounds for i in 1:n
        total += counts[i]
    end
    new_precursors = Vector{UInt32}(undef, total)
    pos = 1
    @inbounds for t in 1:nchunks
        buf = out_bufs[t]
        isempty(buf) && continue
        copyto!(new_precursors, pos, buf, 1, length(buf))
        pos += length(buf)
    end

    off = 0
    @inbounds for i in 1:n
        si = all_scan_idxs[i]
        c  = Int(counts[i])
        if c == 0
            scan_to_prec_idx[si] = missing
        else
            scan_to_prec_idx[si] = (off + 1):(off + c)
            off += c
        end
    end
    return new_precursors
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


"""Strided round-robin scan partition, exactly sized (no push!)."""
function _round_robin_tasks(all_scan_idxs::Vector{Int}, n_threads::Int)
    n = length(all_scan_idxs)
    T = min(n_threads, max(n, 1))
    tasks = [(t, Vector{Int}(undef, length(t:T:n))) for t in 1:T]
    @inbounds for t in 1:T
        v = last(tasks[t]); p = 1
        for i in t:T:n
            v[p] = all_scan_idxs[i]; p += 1
        end
    end
    filter!(tt -> !isempty(last(tt)), tasks)
    return tasks
end


"""
    _zt_report_adjacent_overlap(scan_to_prec_idx, precursors_passed, spectra)

For every pair of consecutive MS2 scans in the same cycle, count |A ∩ B| against |A|, |B|.
Reports the mean fraction of a bin's candidates already present in the previous bin (what a
per-thread template cache could reuse) and the mean Jaccard. Diagnostic only.
"""
function _zt_report_adjacent_overlap(scan_to_prec_idx, precursors_passed::Vector{UInt32},
                                     spectra::MassSpecData)
    cyc = getCycleIdxs(spectra); n = length(spectra)
    tot_b = 0; tot_inter = 0; tot_union = 0; npairs = 0
    reuse_fracs = Float64[]
    @inbounds for si in 2:n
        (getMsOrder(spectra, si) == 2 && getMsOrder(spectra, si - 1) == 2 && cyc[si] == cyc[si - 1]) || continue
        ra = scan_to_prec_idx[si - 1]; rb = scan_to_prec_idx[si]
        (ismissing(ra) || ismissing(rb)) && continue
        ia = first(ra); ib = first(rb); inter = 0
        while ia <= last(ra) && ib <= last(rb)
            a = precursors_passed[ia]; b = precursors_passed[ib]
            if a == b; inter += 1; ia += 1; ib += 1
            elseif a < b; ia += 1
            else; ib += 1
            end
        end
        nb = length(rb); na = length(ra)
        tot_b += nb; tot_inter += inter; tot_union += na + nb - inter; npairs += 1
        push!(reuse_fracs, inter / nb)
    end
    npairs == 0 && return
    sort!(reuse_fracs)
    q(p) = reuse_fracs[clamp(round(Int, p * length(reuse_fracs)), 1, length(reuse_fracs))]
    @user_info "ZT adjacent-bin candidate overlap: $npairs pairs; reusable fraction of a bin's " *
               "candidates (present in previous bin) = $(round(100 * tot_inter / tot_b; digits=1))% " *
               "(median $(round(100 * q(0.5); digits=1))%, p10 $(round(100 * q(0.1); digits=1))%); " *
               "Jaccard = $(round(100 * tot_inter / tot_union; digits=1))%"
end


"""
    _zt_thin_outer_bins!(scan_to_prec_idx, precursors_passed, spectra, prec_mzs, geom, c)

EXPERIMENT: on odd-position scans of each cycle, keep only the candidates whose precursor m/z
lies within `c` bins of the scan centre (the core of their meta-scan); drop the rest. Even scans
are untouched. So every meta-scan keeps all bins within ±c of its centre and every second bin
outside. Rebuilds `precursors_passed` and reindexes in place (mirrors filter_low_scan_candidates!).
"""
function _zt_thin_outer_bins!(scan_to_prec_idx, precursors_passed::Vector{UInt32}, spectra::MassSpecData,
                              prec_mzs::AbstractVector{Float32}, geom::ZTGeometry, c::Int)
    cmzs = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
    lim = Float32(c) * geom.bin_step + geom.bin_step / 2
    odd = falses(length(spectra))
    for r in zt_cycle_scan_ranges(spectra), (pos, si) in enumerate(r); isodd(pos) && (odd[si] = true); end
    new_passed = UInt32[]; sizehint!(new_passed, length(precursors_passed))
    n_before = length(precursors_passed)
    @inbounds for si in eachindex(scan_to_prec_idx)
        rng = scan_to_prec_idx[si]; ismissing(rng) && continue
        start = length(new_passed) + 1
        if odd[si]
            cm = cmzs[si]
            for i in rng
                pid = precursors_passed[i]
                abs(prec_mzs[pid] - cm) <= lim && push!(new_passed, pid)
            end
        else
            for i in rng; push!(new_passed, precursors_passed[i]); end
        end
        scan_to_prec_idx[si] = length(new_passed) >= start ? (start:length(new_passed)) : missing
    end
    @user_info "ZT thin-outer-bins experiment (core ±$c): candidates $n_before -> $(length(new_passed))"
    return new_passed
end
