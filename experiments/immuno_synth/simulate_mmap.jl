# EXPERIMENT: replay a search's detailed-fragment accesses against candidate precursor orderings, to judge how a
# memory-mapped detailed_fragments file would perform. Inputs: a synthetic library's precursors.arrow (mz, irt),
# the harness's exact-candidate capture (scan_idx, precursor_idx) and the raw file (scan RT / center m/z, ordered as
# partitionScansToThreads orders them: 1-min RT bins, m/z within; threads move through that order together).
# Fixed-size records (`rec_bytes` per precursor) laid out in each ordering, 4 KB pages.
#
# usage: julia --project=PIONEER_DIR simulate_mmap.jl PRECURSORS.arrow CAPTURE.arrow RAW.arrow OUT.csv
using Pioneer, Arrow, DataFrames, CSV, Random, Printf, Statistics
const P = Pioneer
const PAGE = 4096

logmsg(s) = (println(round(time() - T0, digits = 1), "s  ", s); flush(stdout))
const T0 = time()

"Position of each precursor when sorted by iRT, cut into blocks of `w` iRT (as chronologer_parse does), m/z within."
function block_order(irt::Vector{Float32}, mz::Vector{Float32}, w::Float32)
    perm = sortperm(irt; alg = Base.Sort.DEFAULT_STABLE)
    if w > 0
        start = 1; start_irt = irt[perm[1]]
        for i in 2:length(perm)
            if irt[perm[i]] - start_irt > w
                sort!(view(perm, start:i-1), by = p -> mz[p]); start = i; start_irt = irt[perm[i]]
            end
        end
        sort!(view(perm, start:length(perm)), by = p -> mz[p])
    end
    return invperm32(perm)
end
function invperm32(perm::Vector{Int})
    pos = Vector{UInt32}(undef, length(perm))
    @inbounds for (r, p) in enumerate(perm); pos[p] = UInt32(r - 1); end
    return pos
end

"Exact LRU over page ids: number of misses and how many misses landed within +32 pages of the previous miss."
function lru_misses(pages::Vector{Int32}, n_pages::Int, cap::Int)
    prev = zeros(Int32, n_pages); nxt = zeros(Int32, n_pages); incache = falses(n_pages)
    head = Int32(0); tail = Int32(0); size = 0; misses = 0; near = 0; last_miss = Int32(-1000)
    @inbounds for pg in pages
        if incache[pg]
            pg == head && continue
            # unlink, then push front
            p, n = prev[pg], nxt[pg]
            p != 0 && (nxt[p] = n); n != 0 && (prev[n] = p); pg == tail && (tail = p)
            prev[pg] = 0; nxt[pg] = head; head != 0 && (prev[head] = pg); head = pg
        else
            misses += 1
            (0 <= pg - last_miss <= 32) && (near += 1); last_miss = pg
            if size == cap  # evict tail
                t = tail; incache[t] = false; tail = prev[t]; tail != 0 && (nxt[tail] = 0); size -= 1
                t == head && (head = Int32(0))
            end
            incache[pg] = true; prev[pg] = 0; nxt[pg] = head; head != 0 && (prev[head] = pg); head = pg
            tail == 0 && (tail = pg); size += 1
        end
    end
    return misses, near
end

"Bytes read when the unique pages are fetched in order, merging runs whose gap is at most `gap` pages."
function gather_bytes(uniq::Vector{Int32}, gap::Int)
    isempty(uniq) && return 0
    total = 0; run_start = uniq[1]; run_end = uniq[1]
    for pg in @view uniq[2:end]
        if pg - run_end - 1 <= gap
            run_end = pg
        else
            total += run_end - run_start + 1; run_start = run_end = pg
        end
    end
    return (total + run_end - run_start + 1) * PAGE
end

function main(prec_path, cap_path, raw_path, out_path)
    pt = Arrow.Table(prec_path)
    mz = Vector{Float32}(pt.mz); irt = Vector{Float32}(pt.irt); N = length(mz)
    logmsg("precursors: $N")

    # scan processing order (as partitionScansToThreads)
    spectra = P.BasicMassSpecData(raw_path)
    rts = P.getRetentionTimes(spectra); cmz = P.getCenterMzs(spectra); mso = P.getMsOrders(spectra)
    order = Int64[i for i in eachindex(mso) if mso[i] == 2]
    P._sort_scans_by_mz_in_rt_bins!(order, rts, cmz)
    scan_rank = zeros(Int32, length(rts)); for (r, s) in enumerate(order); scan_rank[s] = r; end
    rt_bin = zeros(Int32, length(rts)); for s in order; rt_bin[s] = Int32(floor(rts[s])); end

    cap = Arrow.Table(cap_path)
    cs = Vector{Int32}(cap.scan_idx); cp = Vector{UInt32}(cap.precursor_idx)
    acc = sortperm(scan_rank[cs]; alg = Base.Sort.DEFAULT_STABLE)   # accesses in processing order
    cs = cs[acc]; cp = cp[acc]
    n_unique_prec = length(unique(cp))
    logmsg("accesses: $(length(cp)), unique precursors: $n_unique_prec, scans: $(length(order))")

    orderings = [("iRT blocks 3.0 (current)", () -> block_order(irt, mz, 3.0f0)),
                 ("iRT blocks 1.0", () -> block_order(irt, mz, 1.0f0)),
                 ("iRT blocks 0.5", () -> block_order(irt, mz, 0.5f0)),
                 ("iRT only", () -> block_order(irt, mz, 0.0f0)),
                 ("m/z only", () -> invperm32(sortperm(mz))),
                 ("random", () -> invperm32(randperm(MersenneTwister(1), N)))]
    rows = DataFrame()
    for (name, mk) in orderings
        pos = mk(); GC.gc()
        for rec in (280, 200)
            n_pages = cld(N * rec, PAGE)
            # page accesses in processing order (a record can straddle two pages)
            pages = Int32[]; sizehint!(pages, 2 * length(cp))
            bins = Int32[]; sizehint!(bins, 2 * length(cp))
            for (s, p) in zip(cs, cp)
                b0 = Int(pos[p]) * rec
                for pg in (b0 ÷ PAGE):((b0 + rec - 1) ÷ PAGE)
                    push!(pages, Int32(pg + 1)); push!(bins, rt_bin[s])
                end
            end
            uniq = sort!(unique(pages))
            # active set per 1-min RT bin of scans
            active = Int[]; i = 1
            while i <= length(pages)
                j = i; while j < length(pages) && bins[j + 1] == bins[i]; j += 1; end
                push!(active, length(unique(view(pages, i:j)))); i = j + 1
            end
            row = Dict{Symbol, Any}(:ordering => name, :rec_bytes => rec, :file_gb => N * rec / 1e9,
                :needed_gb => n_unique_prec * rec / 1e9, :unique_pages_gb => length(uniq) * PAGE / 1e9,
                :page_fraction => length(uniq) / n_pages,
                :active_max_gb => maximum(active) * PAGE / 1e9, :active_median_gb => median(active) * PAGE / 1e9,
                :gather_gb_gap0 => gather_bytes(uniq, 0) / 1e9, :gather_gb_gap64k => gather_bytes(uniq, 16) / 1e9,
                :gather_gb_gap1m => gather_bytes(uniq, 256) / 1e9)
            for cache_gb in (2, 8, 16)
                m, near = lru_misses(pages, n_pages, round(Int, cache_gb * 1e9 / PAGE))
                row[Symbol("lru$(cache_gb)g_read_gb")] = m * PAGE / 1e9
                row[Symbol("lru$(cache_gb)g_seq_frac")] = near / max(m, 1)
            end
            push!(rows, row; cols = :union)
            logmsg(@sprintf("%-26s rec=%d: file %.1f GB, unique pages %.1f GB (%.0f%%), active max %.2f GB, LRU 8 GB reads %.1f GB",
                name, rec, row[:file_gb], row[:unique_pages_gb], 100 * row[:page_fraction], row[:active_max_gb], row[:lru8g_read_gb]))
        end
        pos = nothing; GC.gc()
    end
    CSV.write(out_path, rows)
    show(stdout, rows; allcols = true); println()
end

main(ARGS...)
