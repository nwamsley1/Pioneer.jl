# Scanning-quad (ZT) main search: deconvolve in cycle-aligned chunks and collapse each chunk to
# meta-PSMs before the next, so the per-bin deconvolved table is never resident for the whole file.

"""
Target candidates (precursor x scan pairs after expansion) per main-search chunk on a
scanning-quad file. Exact, not predicted: the fragment index + expansion run once over the whole
file before chunking. Rows are 0.13-0.26 of candidates on the ZT files seen, so 60M candidates
is ~10-15M raw rows per chunk. nano15 (634M candidates, 164M rows in one pass) swapped a 48 GB
machine; 60M -> 30M (2026-09-17) cut the chunk plateau 26.9 -> 22.5 GB on EV1109 at no time cost.
"""
const ZT_CHUNK_CANDIDATES = 30_000_000

"""
    zt_cycle_scan_ranges(spectra) -> Vector{UnitRange{Int}}

Contiguous MS2 scan ranges, one per acquisition cycle, in scan order.
"""
function zt_cycle_scan_ranges(spectra::MassSpecData)
    cycles = getCycleIdxs(spectra)
    n = length(spectra)
    out = UnitRange{Int}[]
    i = 1
    while i <= n
        if getMsOrder(spectra, i) != 2
            i += 1; continue
        end
        c = cycles[i]; j = i
        while j + 1 <= n && getMsOrder(spectra, j + 1) == 2 && cycles[j + 1] == c
            j += 1
        end
        push!(out, i:j)
        i = j + 1
    end
    return out
end

"""
    zt_candidate_count(cycle_ranges, scan_to_prec_idx) -> Int

Exact number of (precursor, scan) candidates in the given cycles, from the per-scan ranges the
fragment index + expansion produced.
"""
function zt_candidate_count(cycle_ranges::Vector{UnitRange{Int}},
                            scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}})
    n = 0
    @inbounds for r in cycle_ranges, si in r
        rng = scan_to_prec_idx[si]
        ismissing(rng) || (n += length(rng))
    end
    return n
end

"""
    zt_cycle_chunks_by_candidates(ranges, scan_to_prec_idx, target) -> Vector{Vector{UnitRange{Int}}}

Group consecutive cycles into chunks of >= `target` candidates. Every chunk boundary is a cycle
boundary, so no meta-scan is split: candidate expansion and the collapse are both confined to
one cycle. Dense elution regions get many small chunks, empty regions one large chunk.
"""
function zt_cycle_chunks_by_candidates(ranges::Vector{UnitRange{Int}},
                                       scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                                       target::Int)
    chunks = Vector{Vector{UnitRange{Int}}}()
    cur = UnitRange{Int}[]; acc = 0
    for r in ranges
        push!(cur, r); acc += zt_candidate_count([r], scan_to_prec_idx)
        if acc >= target
            push!(chunks, cur); cur = UnitRange{Int}[]; acc = 0
        end
    end
    isempty(cur) || push!(chunks, cur)
    return chunks
end

"""
    zt_thread_tasks(cycle_ranges, k, n_threads) -> Vector{Tuple{Int, Vector{Int}}}

Deal one chunk's cycles to threads so that at any moment every thread is working the SAME m/z
segment of the ramp, on different cycles. Segment width is one meta-scan (2k+1 bins). Thread t
owns cycles t, t+T, t+2T, ... and walks them segment-major: segment 1 of each of its cycles, then
segment 2, and so on. All threads therefore touch the same precursors' library entries at the
same time (cache locality), and each thread processes a full meta-scan width consecutively (what
a warm-started solver and per-precursor template reuse will need). Per-thread vectors are exactly
sized up front; nothing grows.
"""
function zt_thread_tasks(cycle_ranges::Vector{UnitRange{Int}}, k::Int, n_threads::Int)
    C = length(cycle_ranges)
    W = 2k + 1
    T = min(n_threads, C)
    counts = zeros(Int, T)
    @inbounds for c in 1:C
        counts[mod1(c, T)] += length(cycle_ranges[c])
    end
    tasks = [(t, Vector{Int}(undef, counts[t])) for t in 1:T]
    pos = ones(Int, T)
    L = maximum(length, cycle_ranges)
    S = cld(L, W)
    @inbounds for s in 1:S, c in 1:C
        r = cycle_ranges[c]
        lo = first(r) + (s - 1) * W
        hi = min(lo + W - 1, last(r))
        lo > hi && continue
        t = mod1(c, T); v = last(tasks[t]); p = pos[t]
        for si in lo:hi
            v[p] = si; p += 1
        end
        pos[t] = p
    end
    return tasks
end

"""
    zt_chunked_deconvolution(geom, params, deconv, spectra, search_context, ms_file_idx,
                             scan_to_prec_idx) -> Union{Nothing, DataFrame}

`library_search` hook. On a ZT file in the main search, group the cycles into chunks of
~`ZT_CHUNK_CANDIDATES` candidates, deconvolve each (`deconv(thread_tasks)`, segment-major tasks
from `zt_thread_tasks`) and reduce it with `_zt_reduce_chunk` before the next. Returns the
concatenated meta-PSM table, or `nothing` when the file and stage are searched as usual.

The fragment index and expansion ran once over the whole file, so each scan's candidate count is
exact. Chunks are cycle-aligned because expansion and the collapse never cross a cycle.
"""
zt_chunked_deconvolution(::Nothing, params, deconv, spectra, search_context, ms_file_idx,
                         scan_to_prec_idx) = nothing

function zt_chunked_deconvolution(g::ZTGeometry, params, deconv, spectra::MassSpecData,
                                  search_context::SearchContext, ms_file_idx::Int64,
                                  scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}})
    (params isa MainSearchParameters && g.metascan_k > 0) || return nothing
    k = Int(g.metascan_k)
    precursors = getPrecursors(getSpecLib(search_context))
    bitvec_rank_table = getBitVecExcessRanks(search_context, ms_file_idx)
    lookups = ZTCollapseLookups(spectra, precursors)
    t0 = time()
    chunks = zt_cycle_chunks_by_candidates(zt_cycle_scan_ranges(spectra), scan_to_prec_idx,
                                           ZT_CHUNK_CANDIDATES)
    # More than one chunk: each chunk's meta-PSMs go to disk, sorted, and MainSearch merges them into
    # precursor-complete parts, so no more than one chunk is held in memory (ZT/partitioned_scoring.jl).
    on_disk = length(chunks) > 1
    psm_dir = zt_psm_dir(search_context, ms_file_idx)
    if on_disk
        rm(psm_dir; recursive = true, force = true); mkpath(psm_dir)
    end
    chunk_paths = String[]; n_meta = 0
    reduced = DataFrame(); n_raw_total = 0; max_raw = 0
    t_deconv = 0.0; t_reduce = 0.0
    for (ci, chunk) in enumerate(chunks)
        t_c = time()
        raw_all = deconv(zt_thread_tasks(chunk, k, Threads.nthreads()))
        raw = length(raw_all) == 1 ? raw_all[1] : vcat(raw_all...)
        t_d = time() - t_c
        n_raw = nrow(raw); n_raw_total += n_raw; max_raw = max(max_raw, n_raw)
        part = _zt_reduce_chunk(raw, spectra, search_context, ms_file_idx, precursors, g,
                                bitvec_rank_table, lookups)
        raw = nothing; raw_all = nothing
        n_meta += nrow(part)
        if on_disk
            if nrow(part) > 0
                _assert_precursor_scan_sorted(part)
                path = joinpath(psm_dir, "chunk_$(lpad(ci, 4, '0')).arrow")
                Arrow.write(path, part); push!(chunk_paths, path)
            end
        else
            reduced = isempty(reduced) ? part : (append!(reduced, part); reduced)
        end
        t_deconv += t_d; t_reduce += time() - t_c - t_d
        @debug_l1 "ZT chunk $ci/$(length(chunks)): $(length(chunk)) cycles, " *
                  "$(zt_candidate_count(chunk, scan_to_prec_idx)) candidates -> $n_raw raw rows " *
                  "(deconv $(round(t_d; digits=1))s) -> $(nrow(part)) meta-PSMs"
    end
    # The whole-file candidate index and its expansion scratch are dead once the chunks are done;
    # collect them before the meta-PSM table is permuted and featurised, which is where the
    # resident set peaks on ZT files (measured: EV1109 29.8 GB at that point).
    GC.gc()
    @user_info "ZT main search: $(length(chunks)) chunks, $n_raw_total raw rows " *
               "(largest chunk $max_raw) -> $n_meta meta-PSMs; " *
               "deconv $(round(t_deconv; digits=1))s, reduce $(round(t_reduce; digits=1))s"
    on_disk || return reduced
    # The chunk files are merged into precursor-complete parts while MainSearch scores them, after
    # this call has returned and the whole-file candidate index is freed.
    setZTPsmPartitions!(search_context, ms_file_idx, chunk_paths)
    return DataFrame()          # the meta-PSMs are on disk
end

"""
Per-chunk reduce: the two per-scan feature passes MainSearch normally runs after `library_search`
(they need the raw per-bin rows, contiguous by scan), then the meta-scan collapse.
"""
function _zt_reduce_chunk(raw::DataFrame, spectra::MassSpecData, search_context::SearchContext,
                          ms_file_idx::Int64, precursors, g::ZTGeometry, bitvec_rank_table,
                          lookups::ZTCollapseLookups)
    nrow(raw) == 0 && return raw
    @alloc_bucket "scan_competition_features" add_scan_competition_features!(raw)
    @alloc_bucket "ms1_lookup_features" add_ms1_lookup_features!(raw, spectra, search_context, ms_file_idx)
    return @alloc_bucket "metascan_collapse" collapse_to_metascans(
        raw, spectra, precursors, g; bitvec_rank_table = bitvec_rank_table, lookups = lookups)
end

"""
    zt_collapsed_in_search(geom) -> Bool

Whether `library_search` already ran the per-scan feature passes and the meta-scan collapse for
this file (the chunked ZT main search), so MainSearch must not run them again.
"""
zt_collapsed_in_search(::Nothing) = false
zt_collapsed_in_search(g::ZTGeometry) = g.metascan_k > 0
