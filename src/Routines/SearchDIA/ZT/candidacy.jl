# Scanning-quad (ZT) candidacy: re-anchor fragment-index emissions on the precursor's own Q1
# bin and expand each candidate across its meta-scan.

"""
    map_any_hit_to_center!(scan_to_prec_idx, precursors_passed, spectra, all_scan_idxs,
                           prec_mzs, geom) -> Vector{UInt32}

Scanning-quad (ZT) wide-emit candidacy. The fragment index runs with a WIDENED box, so a
precursor can be emitted from any bin of its meta-scan; this maps every emission back to the
precursor's OWN bin, deduped per (cycle, precursor). Net effect: a precursor survives if it
cleared the bitvec in ANY of its bins, rather than requiring the center bin itself to clear. Different bins expose different fragment subsets, so those are real second
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

Three parallel phases over scans, each writing only its own scan's slot; `scan_to_prec_idx` is
read as it was BEFORE expansion and rewritten only in the final serial phase. Exact per-scan
counts are computed first so the output is allocated ONCE at its final size: the previous
per-thread append!-grown buffers plus a concatenation copy allocated 14.8 GB and held ~7 GB live
to produce a 2.6 GB result on EV1109 (647M candidates).
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

    # --- Phase 1 (parallel, read-only): exact per-scan counts. After re-anchoring each precursor
    # occupies exactly one bin per cycle, so the neighbour union has no duplicates and its size
    # is the plain sum of neighbour range lengths. (Dedup below is kept as a guard; if it ever
    # removes anything the scan's range is simply shorter — gaps in the array are never read.)
    counts = Vector{Int64}(undef, n)
    Threads.@threads for i in 1:n
        si = all_scan_idxs[i]
        ci = getCycleIdx(spectra, si)
        c = 0
        @inbounds for j in -k:k
            sj = si + j
            (sj < 1 || sj > nspec) && continue
            getMsOrder(spectra, sj) == 2 || continue
            getCycleIdx(spectra, sj) == ci || continue
            rng = scan_to_prec_idx[sj]
            ismissing(rng) || (c += length(rng))
        end
        counts[i] = c
    end

    # --- Phase 2 (serial): exact prefix offsets, one allocation of the final size ---
    offsets = Vector{Int64}(undef, n + 1)
    offsets[1] = 0
    @inbounds for i in 1:n
        offsets[i + 1] = offsets[i] + counts[i]
    end
    total = offsets[n + 1]
    new_precursors = Vector{UInt32}(undef, total)

    # --- Phase 3 (parallel, disjoint slots): gather each scan's union in place, sort, dedup ---
    # Reads scan_to_prec_idx as it was BEFORE expansion (not yet rewritten), writes only its
    # own slot of new_precursors. Sorting is in place (no scratch).
    new_counts = Vector{Int64}(undef, n)
    Threads.@threads for i in 1:n
        si = all_scan_idxs[i]
        ci = getCycleIdx(spectra, si)
        lo = offsets[i] + 1
        p = lo
        @inbounds for j in -k:k
            sj = si + j
            (sj < 1 || sj > nspec) && continue
            getMsOrder(spectra, sj) == 2 || continue
            getCycleIdx(spectra, sj) == ci || continue
            rng = scan_to_prec_idx[sj]
            ismissing(rng) && continue
            for r in rng
                new_precursors[p] = precursors_passed[r]; p += 1
            end
        end
        m = p - lo
        if m > 1
            v = @view new_precursors[lo:(p - 1)]
            sort!(v; alg = QuickSort)
            w = 1
            @inbounds for t in 2:m
                if v[t] != v[w]
                    w += 1
                    v[w] = v[t]
                end
            end
            m = w
        end
        new_counts[i] = m
    end

    # --- Phase 4 (serial): reindex ---
    @inbounds for i in 1:n
        si = all_scan_idxs[i]
        c  = new_counts[i]
        scan_to_prec_idx[si] = c == 0 ? missing : (offsets[i] + 1):(offsets[i] + c)
    end
    return new_precursors
end

# --- library_search hooks. Each has a `::Nothing` method (not a ZT file) that returns develop's
# --- behaviour unchanged.

"""
Stages that search a ZT file across the meta-scan: the main search, and quad tuning, whose
triangle fit needs each precursor in all 2k+1 bins with a weight that tracks transmission. The
other tuning stages stay on the bin, because the meta-scan's off-center precursors degrade the
mass-error and NCE fits.
"""
_zt_metascan_stage(params) = params isa MainSearchParameters || params isa QuadTuningSearchParameters

"""
    zt_candidacy_quad_model(geom, qtm) -> QuadTransmissionModel

The model the fragment index uses to admit candidates. On a ZT file the recorded isolation width
is only the Q1 step, so candidacy uses a narrow square box around the precursor's own bin.
"""
zt_candidacy_quad_model(::Nothing, qtm::QuadTransmissionModel) = qtm
zt_candidacy_quad_model(g::ZTGeometry, ::QuadTransmissionModel) = SquareQuadModel(zt_candidacy_overhang(g))

"""
    zt_deconv_quad_model(geom, params, qtm) -> QuadTransmissionModel

The model deconvolution uses. On a ZT file the meta-scan stages use the installed flat box that
spans the meta-scan (`ensure_zt_geometry!`); the other stages deconvolve on the bin alone.
"""
zt_deconv_quad_model(::Nothing, params, qtm::QuadTransmissionModel) = qtm
zt_deconv_quad_model(::ZTGeometry, params, qtm::QuadTransmissionModel) =
    _zt_metascan_stage(params) ? qtm : SquareQuadModel(0.0f0)

"""
    zt_expand_candidates!(geom, params, scan_to_prec_idx, precursors_passed, spectra,
                          all_scan_idxs, prec_mzs) -> Vector{UInt32}

On a ZT file in a meta-scan stage: re-anchor every emission on the precursor's own bin
(`map_any_hit_to_center!`) and expand it across ±k bins (`expand_to_metascans!`). Returns
`precursors_passed` unchanged otherwise.
"""
zt_expand_candidates!(::Nothing, params, scan_to_prec_idx, precursors_passed::Vector{UInt32},
                      spectra, all_scan_idxs, prec_mzs) = precursors_passed

function zt_expand_candidates!(g::ZTGeometry, params, scan_to_prec_idx,
                               precursors_passed::Vector{UInt32}, spectra::MassSpecData,
                               all_scan_idxs::Vector{Int}, prec_mzs::AbstractVector{Float32})
    _zt_metascan_stage(params) || return precursors_passed
    k = Int(g.metascan_k)
    n_emitted = length(precursors_passed)
    precursors_passed = map_any_hit_to_center!(scan_to_prec_idx, precursors_passed, spectra,
                                               all_scan_idxs, prec_mzs, g)
    n_center = length(precursors_passed)
    k > 0 || return precursors_passed
    precursors_passed = expand_to_metascans!(scan_to_prec_idx, precursors_passed, spectra,
                                             all_scan_idxs, k)
    @debug_l1 "ZT candidacy (k=$k): $n_emitted emitted -> $n_center center -> " *
              "$(length(precursors_passed)) expanded"
    return precursors_passed
end
