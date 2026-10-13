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
Candidates per main-search chunk on a file that is not scanning-quad. A file with more candidates
than this is deconvolved one chunk of scans at a time and its PSMs go to disk (see
`chunked_main_deconvolution`), so the per-thread scored-PSM buffers and the gathered table hold one
chunk, not the file. Same budget as `ZT_CHUNK_CANDIDATES`.
"""
const MAIN_CHUNK_CANDIDATES = 30_000_000

"""
Peak disk use of a chunked MainSearch file per fragment-index candidate, with margin: the round
files (~195 B per PSM) and the featured parts (~280 B per PSM) are on disk together while the parts
are built, and files keep 0.13-0.24 PSMs per candidate (65.5M PSMs from 279M candidates measured
33 GB at peak).
"""
const MAIN_SPILL_BYTES_PER_CANDIDATE = 150

"""Disk space the chunked path of a file with `n_candidates` candidates needs, with 2 GB to spare."""
main_spill_bytes(n_candidates::Integer) = n_candidates * MAIN_SPILL_BYTES_PER_CANDIDATE + 2_000_000_000

"""Bytes available to this process on the volume holding `path` (`typemax(Int)` if unknown)."""
function disk_free_bytes(path::AbstractString)
    try
        return Int(diskstat(isdir(path) ? path : dirname(path)).available)
    catch
        return typemax(Int)
    end
end

"""
    main_thread_chunks(thread_tasks, scan_to_prec_idx, n_chunks) -> Vector{Vector{Tuple{Int, Vector{Int}}}}

Split a MainSearch file's deconvolution tasks into `n_chunks` rounds. Each task (thread `t` and its
scans, as `partition_scans` dealt them) is cut into `n_chunks` consecutive segments holding about
equal numbers of candidates; round `c` is segment `c` of every task.

Every thread therefore deconvolves exactly the scans it would in one pass, in the same order. That
matters: a thread warm-starts each fit from the weights its previous scans found for the same
precursors (`initialize_weights!` / `update_precursor_weights!`), so a scan's result depends on the
scans its thread fitted before it.
"""
function main_thread_chunks(thread_tasks::AbstractVector,
                            scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}}, n_chunks::Int)
    rounds = [Tuple{Int, Vector{Int}}[] for _ in 1:n_chunks]
    for (t, scans) in thread_tasks
        cand = [(r = scan_to_prec_idx[si]; ismissing(r) ? 0 : length(r)) for si in scans]
        total = sum(cand; init = 0)
        c = 1; acc = 0; start = 1
        for (i, n) in enumerate(cand)
            acc += n
            # close segment c once it holds its share of this task's candidates
            if c < n_chunks && acc >= c * total / n_chunks
                push!(rounds[c], (t, scans[start:i])); start = i + 1; c += 1
            end
        end
        push!(rounds[c], (t, scans[start:end]))
    end
    return [filter(tk -> !isempty(last(tk)), r) for r in rounds]
end

"""
    chunked_main_deconvolution(params, deconv, spectra, search_context, ms_file_idx,
                               scan_to_prec_idx, thread_tasks, n_candidates) -> Union{Nothing, DataFrame}

`library_search` hook for MainSearch on files that are not scanning-quad. With more than
`MAIN_CHUNK_CANDIDATES` candidates, the file is deconvolved in rounds (`main_thread_chunks`); each
round gets the two per-scan feature passes MainSearch otherwise runs after `library_search` (scan
competition and the MS1 lookup; both read only a row's own scan), is sorted by
(precursor_idx, scan_idx) and is written to disk. MainSearch then merges the files into
precursor-complete parts and scores them as it does a multi-chunk ZT file
(`zt_mainsearch_best_partitioned!`), so the file's PSM table is never in memory whole. Returns
`nothing` (search in memory as usual) or a table with no rows (the PSMs are on disk, listed by
`getZTPsmPartitions`; or there were none).

Each thread scores the same scans in the same order as in one pass, and the merged files hold the
rows in (precursor_idx, scan_idx) order, the order the in-memory table is sorted into; so the PSM
rows and their order are the same.
"""
function chunked_main_deconvolution(params, deconv, spectra::MassSpecData,
                                    search_context::SearchContext, ms_file_idx::Int64,
                                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                                    thread_tasks::AbstractVector, n_candidates::Int)
    params isa MainSearchParameters || return nothing
    getZTGeometry(search_context, ms_file_idx) === nothing || return nothing
    n_candidates > MAIN_CHUNK_CANDIDATES || return nothing
    psm_dir = zt_psm_dir(search_context, ms_file_idx)
    # Not enough free disk for the round files and the featured parts (both are on disk while the
    # parts are built): search in memory instead, which gives the same results.
    need = main_spill_bytes(n_candidates)
    free = disk_free_bytes(getDataOutDir(search_context))
    if free < need
        @user_warn "Main search: $(round(free / 1e9, digits = 1)) GB free on the results disk, about " *
                   "$(round(need / 1e9, digits = 1)) GB needed to keep this file's PSMs on disk; " *
                   "searching it in memory (same results, more memory)"
        return nothing
    end
    chunks = main_thread_chunks(thread_tasks, scan_to_prec_idx, cld(n_candidates, MAIN_CHUNK_CANDIDATES))
    rm(psm_dir; recursive = true, force = true); mkpath(psm_dir)
    workspace = Int32SortPermWorkspace()
    chunk_paths = String[]; n_rows = 0; max_rows = 0; empty_schema = nothing
    t_deconv = 0.0; t_features = 0.0; t_write = 0.0
    for (ci, tasks) in enumerate(chunks)
        isempty(tasks) && continue
        t_c = time()
        parts = deconv(tasks)
        raw = length(parts) == 1 ? parts[1] : vcat(parts...)
        parts = nothing
        t_d = time() - t_c
        empty_schema === nothing && (empty_schema = raw[1:0, :])
        if nrow(raw) > 0
            t_features += @elapsed begin
                @alloc_bucket "scan_competition_features" add_scan_competition_features!(raw)
                @alloc_bucket "ms1_lookup_features" add_ms1_lookup_features!(raw, spectra, search_context, ms_file_idx)
                permute_psms_by_precursor_idx!(raw, workspace)
            end
            _assert_precursor_scan_sorted(raw)
            path = joinpath(psm_dir, "chunk_$(lpad(ci, 4, '0')).arrow")
            t_write += @elapsed Arrow.write(path, raw)
            push!(chunk_paths, path)
            # Diagnostic: keep a copy of each round file (e.g. to benchmark the merge on it)
            keep = get(ENV, "PIONEER_DIAG_KEEP_CHUNKS", "")
            isempty(keep) || (mkpath(keep); cp(path, joinpath(keep, basename(path)); force = true))
        end
        n_rows += nrow(raw); max_rows = max(max_rows, nrow(raw))
        t_deconv += t_d
        @debug_l1 "Main chunk $ci/$(length(chunks)): $(sum(t -> length(last(t)), tasks)) scans -> $(nrow(raw)) PSMs " *
                  "(deconv $(round(t_d; digits=1))s)"
        raw = nothing
    end
    @user_info "Main search: $(length(chunks)) chunks, $n_rows PSMs (largest chunk $max_rows) on disk; " *
               "deconv $(round(t_deconv; digits=1))s, per-scan features $(round(t_features; digits=1))s, " *
               "write $(round(t_write; digits=1))s"
    if isempty(chunk_paths)
        rm(psm_dir; recursive = true, force = true)
        return empty_schema === nothing ? nothing : empty_schema
    end
    setZTPsmPartitions!(search_context, ms_file_idx, chunk_paths)
    return DataFrame()                 # the PSMs are on disk
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
