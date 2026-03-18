# ── SIMD primitives for SoA linear scan ──────────────────────────────────────

const F32x8 = NTuple{8, Core.VecElement{Float32}}

@inline function _vbroadcast8(x::Float32)::F32x8
    ntuple(_ -> Core.VecElement(x), Val(8))
end

@inline function _vload8(arr::Vector{Float32}, i::Int)::F32x8
    unsafe_load(Ptr{F32x8}(pointer(arr, i)))
end

@inline function _vcmpge_mask(a::F32x8, b::F32x8)::UInt8
    Core.Intrinsics.llvmcall("""
        %cmp = fcmp oge <8 x float> %0, %1
        %mask = bitcast <8 x i1> %cmp to i8
        ret i8 %mask
    """, UInt8, Tuple{F32x8, F32x8}, a, b)
end

"""
    _find_first_ge(highs, start, stop, threshold) -> UInt32

SIMD-accelerated scan: find first index i in start:stop where highs[i] >= threshold.
Returns stop + 1 if no match found.
"""
@inline function _find_first_ge(highs::Vector{Float32}, start::UInt32, stop::UInt32, threshold::Float32)::UInt32
    thr_vec = _vbroadcast8(threshold)
    i = Int(start)
    stop_i = Int(stop)
    # SIMD: 8 elements per iteration
    while i + 7 <= stop_i
        mask = _vcmpge_mask(_vload8(highs, i), thr_vec)
        mask != 0x00 && return UInt32(i + trailing_zeros(mask))
        i += 8
    end
    # Scalar tail
    while i <= stop_i
        @inbounds highs[i] >= threshold && return UInt32(i)
        i += 1
    end
    return stop + one(UInt32)
end

"""
    _findFirstFragBin_soa(highs, lb, ub, frag_min) -> UInt32

Branchless binary search on the SoA `highs` array.
Equivalent to Pioneer.findFirstFragmentBin but operates on a contiguous Float32 array.
"""
@inline function _findFirstFragBin_soa(highs::Vector{Float32}, lb::UInt32, ub::UInt32, frag_min::Float32)
    @inbounds @fastmath begin
        len = ub - lb + one(UInt32)
        mid = len >>> 0x01
        base = lb
        while len > 1
            base += (highs[base + mid - one(UInt32)] < frag_min) * mid
            len -= mid
            mid = len >>> 0x01
        end
    end
    return base
end

"""
    _findFirstFragBin_hybrid(highs, lb, ub, frag_min, simd_cutoff) -> UInt32

Hybrid binary+SIMD search. Runs branchless binary search iterations until the
remaining range is ≤ simd_cutoff, then finishes with a SIMD linear scan.
Gets the best of both: O(log n) narrowing for large ranges, cache-friendly
SIMD scan for the final stretch.
"""
@inline function _findFirstFragBin_hybrid(highs::Vector{Float32}, lb::UInt32, ub::UInt32,
                                           frag_min::Float32, simd_cutoff::UInt32)
    @inbounds @fastmath begin
        len = ub - lb + one(UInt32)
        mid = len >>> 0x01
        base = lb
        # Binary search until range is small enough for SIMD
        while len > simd_cutoff
            base += (highs[base + mid - one(UInt32)] < frag_min) * mid
            len -= mid
            mid = len >>> 0x01
        end
    end
    # Finish with SIMD linear scan over the remaining ≤simd_cutoff elements
    return _find_first_ge(highs, base, base + len - one(UInt32), frag_min)
end

# ── Scoring functions ────────────────────────────────────────────────────────

"""
    searchFragmentBinUnconditional!(counter, fragments, frag_id_range)

Score all fragments in the range unconditionally (no binary search on prec_mz).
This is the key performance change — partitioning guarantees all fragments in the
partition have prec_mz within ~partition_width of the query window.

Works with any fragment type that supports `getPrecID` and `getScore`.
"""
@inline function searchFragmentBinUnconditional!(
        prec_id_to_score::Pioneer.Counter{I, UInt8},
        fragments::AbstractVector,
        frag_id_range::UnitRange{UInt32}) where {I}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        Pioneer.inc!(prec_id_to_score, Pioneer.getPrecID(frag), Pioneer.getScore(frag))
    end
    return nothing
end

@inline function searchFragmentBinUnconditional!(
        prec_id_to_score::LocalCounter{I, UInt8},
        fragments::AbstractVector,
        frag_id_range::UnitRange{UInt32}) where {I}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        inc!(prec_id_to_score, Pioneer.getPrecID(frag), Pioneer.getScore(frag))
    end
    return nothing
end

# ── Bitmask scoring functions ─────────────────────────────────────────────

@inline function searchFragmentBinBitmask!(
        counter::BitmaskCounter{I},
        fragments::AbstractVector,
        frag_id_range::UnitRange{UInt32}) where {I}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        bitmask_inc!(counter, Pioneer.getPrecID(frag), Pioneer.getScore(frag))
    end
    return nothing
end

function queryFragmentBitmask!(
        counter::BitmaskCounter{I},
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{T}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32) where {I, T<:AbstractFloat}

    lower_bound_guess, upper_bound_guess = Pioneer.exponentialFragmentBinSearch(
        frag_bins, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
        frag_mz_min, frag_mz_max, UInt32(4))

    frag_bin_idx = Pioneer.findFirstFragmentBin(
        frag_bins, lower_bound_guess, upper_bound_guess, frag_mz_min)

    @inbounds @fastmath begin
        if iszero(frag_bin_idx)
            return lower_bound_guess, upper_bound_guess
        end
        while frag_bin_idx <= frag_bin_max_idx
            frag_bin = frag_bins[frag_bin_idx]
            if Pioneer.getLow(frag_bin) > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if Pioneer.getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                frag_id_range = Pioneer.getSubBinRange(frag_bin)
                searchFragmentBinBitmask!(counter, fragments, frag_id_range)
                frag_bin_idx += 1
            end
        end
    end
    return lower_bound_guess, upper_bound_guess
end

@inline function _score_partition_bitmask!(
        counter::BitmaskCounter{UInt16},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    p_rt_bins   = getRTBins(partition)
    p_frag_bins = getFragBins(partition)
    p_fragments = getFragments(partition)

    isempty(p_frag_bins) && return nothing
    n_rt = length(p_rt_bins)

    local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
    local_rt_bin_idx > n_rt && return nothing

    @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            lower_bound_guess = min_frag_bin
            upper_bound_guess = min_frag_bin

            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                lower_bound_guess, upper_bound_guess = queryFragmentBitmask!(
                    counter,
                    max_frag_bin,
                    lower_bound_guess,
                    upper_bound_guess,
                    p_frag_bins,
                    p_fragments,
                    frag_min,
                    frag_max
                )
            end
        end

        local_rt_bin_idx += 1
        if local_rt_bin_idx > n_rt
            break
        end
    end

    return nothing
end

# ── Standard scoring functions ───────────────────────────────────────────

"""
    queryFragmentPartitioned!(counter, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
                              frag_bins, fragments, frag_mz_min, frag_mz_max)

Same frag-bin lookup as `queryFragment!` (exponential search + binary search to find
matching frag bins), but calls `searchFragmentBinUnconditional!` instead of
`searchFragmentBin!`. No prec_mz_min/prec_mz_max needed.
"""
function queryFragmentPartitioned!(
        prec_id_to_score::Pioneer.Counter{I, UInt8},
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{T}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32) where {I, T<:AbstractFloat}

    lower_bound_guess, upper_bound_guess = Pioneer.exponentialFragmentBinSearch(
        frag_bins,
        frag_bin_max_idx,
        lower_bound_guess,
        upper_bound_guess,
        frag_mz_min,
        frag_mz_max,
        UInt32(4)
    )

    frag_bin_idx = Pioneer.findFirstFragmentBin(
        frag_bins,
        lower_bound_guess,
        upper_bound_guess,
        frag_mz_min
    )

    @inbounds @fastmath begin
        if iszero(frag_bin_idx)
            return lower_bound_guess, upper_bound_guess
        end

        while frag_bin_idx <= frag_bin_max_idx
            frag_bin = frag_bins[frag_bin_idx]
            if Pioneer.getLow(frag_bin) > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if Pioneer.getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                frag_id_range = Pioneer.getSubBinRange(frag_bin)
                searchFragmentBinUnconditional!(prec_id_to_score, fragments, frag_id_range)
                frag_bin_idx += 1
            end
        end
    end

    return lower_bound_guess, upper_bound_guess
end

# LocalCounter version of queryFragmentPartitioned!
function queryFragmentPartitioned!(
        prec_id_to_score::LocalCounter{I, UInt8},
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{T}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32) where {I, T<:AbstractFloat}

    lower_bound_guess, upper_bound_guess = Pioneer.exponentialFragmentBinSearch(
        frag_bins, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
        frag_mz_min, frag_mz_max, UInt32(4))

    frag_bin_idx = Pioneer.findFirstFragmentBin(
        frag_bins, lower_bound_guess, upper_bound_guess, frag_mz_min)

    @inbounds @fastmath begin
        if iszero(frag_bin_idx)
            return lower_bound_guess, upper_bound_guess
        end
        while frag_bin_idx <= frag_bin_max_idx
            frag_bin = frag_bins[frag_bin_idx]
            if Pioneer.getLow(frag_bin) > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if Pioneer.getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                frag_id_range = Pioneer.getSubBinRange(frag_bin)
                searchFragmentBinUnconditional!(prec_id_to_score, fragments, frag_id_range)
                frag_bin_idx += 1
            end
        end
    end
    return lower_bound_guess, upper_bound_guess
end

"""
    searchScanPartitioned!(counter, pfi, irt_low, irt_high, masses, mass_err_model, quad_func, isotope_err_bounds)

For each scan, determine which partitions overlap the quad isolation window, then
search each relevant partition using the unconditional scoring loop.

Each partition has its own RT bins, so we find the correct starting RT bin per
partition via binary search on `irt_low`.

Works with both PartitionedFragmentIndex and CompactPartitionedFragmentIndex.
"""
function searchScanPartitioned!(
        prec_id_to_score::Pioneer.Counter{UInt32, UInt8},
        pfi::Union{PartitionedFragmentIndex{T}, CompactPartitionedFragmentIndex{T}},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        quad_transmission_func::Pioneer.QuadTransmissionFunction,
        isotope_err_bounds::Tuple{UInt8, UInt8}
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    prec_min = U(Pioneer.getPrecMinBound(quad_transmission_func) - Pioneer.NEUTRON * first(isotope_err_bounds) / 2)
    prec_max = U(Pioneer.getPrecMaxBound(quad_transmission_func) + Pioneer.NEUTRON * last(isotope_err_bounds) / 2)

    # Determine which partitions overlap the precursor window
    first_k = get_partition_idx(pfi, T(prec_min))
    last_k  = get_partition_idx(pfi, T(prec_max))

    for k in first_k:last_k
        partition = getPartition(pfi, k)
        p_rt_bins   = getRTBins(partition)
        p_frag_bins = getFragBins(partition)
        p_fragments = getFragments(partition)

        isempty(p_frag_bins) && continue
        n_rt = length(p_rt_bins)

        # Find the first RT bin whose upper bound >= irt_low (binary search)
        local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
        local_rt_bin_idx > n_rt && continue

        @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
            sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
            min_frag_bin = first(sub_bin_range)
            max_frag_bin = last(sub_bin_range)

            # Skip empty RT bins (min > max means empty range)
            if min_frag_bin <= max_frag_bin
                lower_bound_guess = min_frag_bin
                upper_bound_guess = min_frag_bin

                for mass in masses
                    corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                    frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                    lower_bound_guess, upper_bound_guess = queryFragmentPartitioned!(
                        prec_id_to_score,
                        max_frag_bin,
                        lower_bound_guess,
                        upper_bound_guess,
                        p_frag_bins,
                        p_fragments,
                        frag_min,
                        frag_max
                    )
                end
            end

            local_rt_bin_idx += 1
            if local_rt_bin_idx > n_rt
                break
            end
        end
    end

    return nothing
end

"""
Binary search for the first RT bin whose upper bound (getHigh) >= irt_low.
Returns index in 1:length(rt_bins), or length(rt_bins)+1 if none qualifies.
"""
@inline function _find_rt_bin_start(rt_bins::Vector{Pioneer.FragIndexBin{T}}, irt_low::Float32) where {T}
    lo, hi = 1, length(rt_bins)
    result = hi + 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if Pioneer.getHigh(rt_bins[mid]) >= irt_low
            result = mid
            hi = mid - 1
        else
            lo = mid + 1
        end
    end
    return result
end

"""
    _score_partition!(local_counter, partition, irt_low, irt_high, masses, mass_err_model)

Score one partition's fragments into the local counter. Does not reset the counter.
"""
@inline function _score_partition!(
        local_counter::LocalCounter{UInt16, UInt8},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    p_rt_bins   = getRTBins(partition)
    p_frag_bins = getFragBins(partition)
    p_fragments = getFragments(partition)

    isempty(p_frag_bins) && return nothing
    n_rt = length(p_rt_bins)

    local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
    local_rt_bin_idx > n_rt && return nothing

    @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            lower_bound_guess = min_frag_bin
            upper_bound_guess = min_frag_bin

            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                lower_bound_guess, upper_bound_guess = queryFragmentPartitioned!(
                    local_counter,
                    max_frag_bin,
                    lower_bound_guess,
                    upper_bound_guess,
                    p_frag_bins,
                    p_fragments,
                    frag_min,
                    frag_max
                )
            end
        end

        local_rt_bin_idx += 1
        if local_rt_bin_idx > n_rt
            break
        end
    end

    return nothing
end

"""
    searchFragmentIndexPartitioned(scan_to_prec_idx, pfi, spectra, thread_task,
                                   search_data, params, qtm, mem, rt_to_irt_spline,
                                   irt_tol, precursor_mzs)

Per-thread entry point for partitioned fragment index search. Same scan loop as the
current `searchFragmentIndex`, but calls `searchScanPartitioned!` and applies a
post-filter to remove false-positive precursors outside the actual quad window.
"""
function searchFragmentIndexPartitioned(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::Union{PartitionedFragmentIndex{Float32}, CompactPartitionedFragmentIndex{Float32}},
        spectra::Pioneer.MassSpecData,
        thread_task::Vector{Int64},
        search_data::S,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32}
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel,
                 S<:Pioneer.SearchDataStructures, P<:Pioneer.FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)

    for scan_idx in thread_task
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, scan_idx) ∉ Pioneer.getSpecOrder(params) && continue

        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)

        counter = Pioneer.getPrecursorScores(search_data)

        # Partitioned search — each partition finds its own RT bin via binary search
        searchScanPartitioned!(
            counter,
            pfi,
            irt_lo,
            irt_hi,
            Pioneer.getMzArray(spectra, scan_idx),
            mem,
            Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx)),
            Pioneer.getIsotopeErrBounds(params)
        )

        # Post-filter: zero out scores for precursors outside the actual quad window
        quad_func = Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx))
        iso_bounds = Pioneer.getIsotopeErrBounds(params)
        actual_prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        actual_prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)

        @inbounds for idx in 1:(Pioneer.getSize(counter) - 1)
            pid = Pioneer.getID(counter, idx)
            pid == 0 && continue
            pmz = precursor_mzs[pid]
            if pmz < actual_prec_min || pmz > actual_prec_max
                counter.counts[pid] = zero(UInt8)
            end
        end

        # Filter and collect results (same as searchFragmentIndex)
        _match_count, _prec_count = Pioneer.filterPrecursorMatches!(counter, Pioneer.getMinIndexSearchScore(params))

        if Pioneer.getID(counter, 1) > 0
            start_idx = prec_id + 1
            n = 1
            while n <= counter.matches
                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring,
                            Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring)))
                end
                precursors_passed_scoring[prec_id] = Pioneer.getID(counter, n)
                n += 1
            end
            scan_to_prec_idx[scan_idx] = start_idx:prec_id
        else
            scan_to_prec_idx[scan_idx] = missing
        end

        Pioneer.reset!(counter)
    end

    return precursors_passed_scoring[1:prec_id]
end

"""
Per-thread entry point for LocalPartitionedFragmentIndex search.

Scores each partition into a small LocalCounter{UInt16, UInt8}, applies the
score threshold and post-filter on local IDs, then converts only passing
precursors to global UInt32 IDs. No global counter needed — the local counter
fits in L1/L2 cache.
"""
function searchFragmentIndexPartitioned(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        thread_task::Vector{Int64},
        search_data::S,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32}
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel,
                 S<:Pioneer.SearchDataStructures, P<:Pioneer.FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)

    # Small counter for partition-local scoring (reused across scans+partitions)
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    local_counter = LocalCounter(UInt16, UInt8, max_local + 1)
    min_score = Pioneer.getMinIndexSearchScore(params)

    for scan_idx in thread_task
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, scan_idx) ∉ Pioneer.getSpecOrder(params) && continue

        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        masses = Pioneer.getMzArray(spectra, scan_idx)
        quad_func = Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx))
        iso_bounds = Pioneer.getIsotopeErrBounds(params)
        actual_prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        actual_prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)

        first_k, last_k = get_partition_range(pfi, Float32(actual_prec_min), Float32(actual_prec_max))
        scan_start = prec_id + 1

        for k in first_k:last_k
            partition = getPartition(pfi, k)

            # Score this partition's fragments into the local counter
            _score_partition!(local_counter, partition, irt_lo, irt_hi, masses, mem)

            # Filter + collect: only convert passing local IDs to global
            l2g = partition.local_to_global
            @inbounds for i in 1:(local_counter.size - 1)
                lid = local_counter.ids[i]
                score = local_counter.counts[lid]
                score < min_score && continue

                # Post-filter: check prec_mz against quad window
                global_pid = l2g[lid]
                pmz = precursor_mzs[global_pid]
                (pmz < actual_prec_min || pmz > actual_prec_max) && continue

                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring,
                            Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring)))
                end
                precursors_passed_scoring[prec_id] = global_pid
            end

            reset!(local_counter)
        end

        if prec_id >= scan_start
            scan_to_prec_idx[scan_idx] = scan_start:prec_id
        else
            scan_to_prec_idx[scan_idx] = missing
        end
    end

    return precursors_passed_scoring[1:prec_id]
end

"""
    searchFragmentIndexPartitionMajor(scan_to_prec_idx, pfi, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

Partition-major search: outer loop over partitions, inner loop fans MS2 scans
across threads. All threads working at any given time read the same partition's
fragment/bin arrays, maximizing shared cache utilization.

Returns a flat vector of global precursor IDs (concatenated across all scans).
"""
function searchFragmentIndexPartitionMajor(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        all_scan_idxs::Vector{Int},
        n_threads::Int,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32}
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel,
                 P<:Pioneer.FragmentIndexSearchParameters}

    min_score = Pioneer.getMinIndexSearchScore(params)
    iso_bounds = Pioneer.getIsotopeErrBounds(params)
    n_scans = length(all_scan_idxs)

    # ── Pre-compute per-scan properties ──────────────────────────────────────
    scan_irt_lo  = Vector{Float32}(undef, n_scans)
    scan_irt_hi  = Vector{Float32}(undef, n_scans)
    scan_prec_min = Vector{Float32}(undef, n_scans)
    scan_prec_max = Vector{Float32}(undef, n_scans)

    for (si, scan_idx) in enumerate(all_scan_idxs)
        irt_lo, irt_hi = Pioneer.getRTWindow(
            rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        quad_func = Pioneer.getQuadTransmissionFunction(
            qtm, Pioneer.getCenterMz(spectra, scan_idx),
            Pioneer.getIsolationWidthMz(spectra, scan_idx))
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        scan_irt_lo[si]   = irt_lo
        scan_irt_hi[si]   = irt_hi
        scan_prec_min[si] = prec_min
        scan_prec_max[si] = prec_max
    end

    # ── Per-thread result buffers and local counters ─────────────────────────
    thread_results = [Vector{Tuple{Int, UInt32}}(undef, 0) for _ in 1:n_threads]
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    thread_counters = [LocalCounter(UInt16, UInt8, max_local + 1) for _ in 1:n_threads]

    # ── Pre-compute partition → scan mapping ────────────────────────────────
    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end

    # ── Build flat work list: (partition_idx, scan_si) pairs ────────────────
    # Group by partition so all threads working on the same partition share cache
    work_items = Tuple{Int, Int}[]  # (partition_idx, scan_si)
    for k in 1:pfi.n_partitions
        isempty(getFragBins(getPartition(pfi, k))) && continue
        for si in partition_to_scans[k]
            push!(work_items, (k, si))
        end
    end
    n_work = length(work_items)

    # ── Partition-major parallel execution ───────────────────────────────────
    # Spawn n_threads persistent workers that iterate all partitions. Within
    # each partition, threads take interleaved slices of scans. All threads
    # process the same partition before moving to the next, sharing cache.
    # Barrier synchronization between partitions via a shared atomic counter.
    partition_barrier = Threads.Atomic{Int}(0)

    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            lc = thread_counters[tid]
            results = thread_results[tid]

            for k in 1:pfi.n_partitions
                relevant = partition_to_scans[k]
                partition = getPartition(pfi, k)
                n_relevant = length(relevant)
                (n_relevant == 0 || isempty(getFragBins(partition))) && continue

                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition!(lc, partition,
                        scan_irt_lo[si], scan_irt_hi[si],
                        Pioneer.getMzArray(spectra, scan_idx), mem)

                    prec_lo = scan_prec_min[si]
                    prec_hi = scan_prec_max[si]
                    @inbounds for i in 1:(lc.size - 1)
                        lid = lc.ids[i]
                        score = lc.counts[lid]
                        score < min_score && continue
                        global_pid = l2g[lid]
                        pmz = precursor_mzs[global_pid]
                        (pmz < prec_lo || pmz > prec_hi) && continue
                        push!(results, (si, global_pid))
                    end

                    reset!(lc)
                    scan_i += n_threads
                end
            end
        end
    end
    fetch.(tasks)

    # ── Collect results into scan_to_prec_idx ────────────────────────────────
    scan_results = [UInt32[] for _ in 1:n_scans]
    for tid in 1:n_threads
        for (si, gpid) in thread_results[tid]
            push!(scan_results[si], gpid)
        end
        empty!(thread_results[tid])
    end

    precursors_passed = UInt32[]
    for si in 1:n_scans
        scan_idx = all_scan_idxs[si]
        pids = scan_results[si]
        if isempty(pids)
            scan_to_prec_idx[scan_idx] = missing
        else
            start = length(precursors_passed) + 1
            append!(precursors_passed, pids)
            scan_to_prec_idx[scan_idx] = start:length(precursors_passed)
        end
    end

    return precursors_passed
end

"""
Partition-major search using BitmaskCounter. Fragments store bit positions
(1, 2, 4) and scoring uses OR. Precursors pass when all bits are set.
"""
function searchFragmentIndexPartitionMajorBitmask(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        all_scan_idxs::Vector{Int},
        n_threads::Int,
        pass_mask::UInt8,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32}
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel}

    iso_bounds = (UInt8(1), UInt8(0))  # match FirstPassSearch defaults
    n_scans = length(all_scan_idxs)

    scan_irt_lo  = Vector{Float32}(undef, n_scans)
    scan_irt_hi  = Vector{Float32}(undef, n_scans)
    scan_prec_min = Vector{Float32}(undef, n_scans)
    scan_prec_max = Vector{Float32}(undef, n_scans)

    for (si, scan_idx) in enumerate(all_scan_idxs)
        irt_lo, irt_hi = Pioneer.getRTWindow(
            rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        quad_func = Pioneer.getQuadTransmissionFunction(
            qtm, Pioneer.getCenterMz(spectra, scan_idx),
            Pioneer.getIsolationWidthMz(spectra, scan_idx))
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        scan_irt_lo[si]   = irt_lo
        scan_irt_hi[si]   = irt_hi
        scan_prec_min[si] = prec_min
        scan_prec_max[si] = prec_max
    end

    thread_results = [Vector{Tuple{Int, UInt32}}(undef, 0) for _ in 1:n_threads]
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    thread_counters = [BitmaskCounter(UInt16, max_local + 1) for _ in 1:n_threads]

    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end

    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            bc = thread_counters[tid]
            results = thread_results[tid]

            for k in 1:pfi.n_partitions
                relevant = partition_to_scans[k]
                partition = getPartition(pfi, k)
                n_relevant = length(relevant)
                (n_relevant == 0 || isempty(getFragBins(partition))) && continue

                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition_bitmask!(bc, partition,
                        scan_irt_lo[si], scan_irt_hi[si],
                        Pioneer.getMzArray(spectra, scan_idx), mem)

                    prec_lo = scan_prec_min[si]
                    prec_hi = scan_prec_max[si]
                    @inbounds for i in 1:(bc.size - 1)
                        lid = bc.ids[i]
                        bc.counts[lid] == pass_mask || continue
                        global_pid = l2g[lid]
                        pmz = precursor_mzs[global_pid]
                        (pmz < prec_lo || pmz > prec_hi) && continue
                        push!(results, (si, global_pid))
                    end

                    reset!(bc)
                    scan_i += n_threads
                end
            end
        end
    end
    fetch.(tasks)

    scan_results = [UInt32[] for _ in 1:n_scans]
    for tid in 1:n_threads
        for (si, gpid) in thread_results[tid]
            push!(scan_results[si], gpid)
        end
        empty!(thread_results[tid])
    end

    precursors_passed = UInt32[]
    for si in 1:n_scans
        scan_idx = all_scan_idxs[si]
        pids = scan_results[si]
        if isempty(pids)
            scan_to_prec_idx[scan_idx] = missing
        else
            start = length(precursors_passed) + 1
            append!(precursors_passed, pids)
            scan_to_prec_idx[scan_idx] = start:length(precursors_passed)
        end
    end

    return precursors_passed
end

# ── Hint-based search functions ──────────────────────────────────────────

const HINT_LINEAR_THRESHOLD = UInt32(32)

"""
    queryFragmentHinted!(counter, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
                          frag_bins, fragments, frag_mz_min, frag_mz_max,
                          hints, prev_mz, linear_threshold)

5-Da direct hint search with advancing lb, SoA layout + SIMD linear scan.

Hint semantics: hints[j] = k where frag_bins.lows[j+k] >= frag_bins.lows[j] + 5 Da.
est_step = hint * (delta_mz / 5.0).

Algorithm:
1. Hint-based lb advancement (provably safe: hint measures exactly 5 Da in lows,
   so advancing by hint when delta >= 5 Da cannot overshoot).
   For delta < 5 Da, linear interpolation with factor=1.0 (validated zero-overshoot).
2. Hint-based UB guess: new_lb + est_step * 1.5.
3. Exponential doubling only if UB guess insufficient (covers ~10% of cases).
4. SIMD _find_first_ge on highs array for small ranges, binary search for large ranges.
5. Return (first_match, ub) so lb advances for next peak.
"""
@inline function queryFragmentHinted!(
        counter,
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::SoAFragBins{T},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32,
        hints::Vector{UInt16},
        prev_mz::Float32,
        linear_threshold::UInt32) where {T<:AbstractFloat}

    fb_lows = frag_bins.lows
    fb_highs = frag_bins.highs
    fb_first = frag_bins.first_bins
    fb_last = frag_bins.last_bins

    # ── Hint-based lb advancement ────────────────────────────────────
    new_lb = lower_bound_guess
    est_step_f = 0.0f0

    if prev_mz > 0.0f0 && lower_bound_guess <= frag_bin_max_idx
        delta_mz = frag_mz_min - @inbounds fb_lows[lower_bound_guess]
        if delta_mz > 0.0f0
            hint = @inbounds hints[lower_bound_guess]
            est_step_f = Float32(hint) * (delta_mz / 5.0f0)

            if delta_mz > 5.0f0
                new_lb = min(lower_bound_guess + UInt32(hint), frag_bin_max_idx)
            else
                advance = max(unsafe_trunc(UInt32, est_step_f), one(UInt32))
                new_lb = min(lower_bound_guess + advance, frag_bin_max_idx)
            end
        end
    end

    # ── Hint-based UB guess ──────────────────────────────────────────
    ub_guess = new_lb + max(one(UInt32), unsafe_trunc(UInt32, est_step_f * 1.5f0) + one(UInt32))
    ub_guess = min(ub_guess, frag_bin_max_idx)

    # ── Check if UB guess is sufficient, exponential doubling if not ─
    ub_final = ub_guess
    if @inbounds fb_highs[ub_final] < frag_mz_max
        step = one(UInt32)
        while @inbounds fb_highs[ub_final] < frag_mz_max
            ub_final += step
            step = step << one(UInt8)
            if ub_final > frag_bin_max_idx
                ub_final = frag_bin_max_idx
                break
            end
        end
    end

    # ── Find first matching bin (SIMD or binary→SIMD) ──────────────
    range_size = ub_final - new_lb + one(UInt32)
    if range_size <= linear_threshold
        # Small range: direct SIMD scan
        frag_bin_idx = _find_first_ge(fb_highs, new_lb, ub_final, frag_mz_min)
    else
        # Large range: binary search narrows to ≤threshold, then SIMD finishes
        frag_bin_idx = _findFirstFragBin_hybrid(fb_highs, new_lb, ub_final,
                                                 frag_mz_min, linear_threshold)
    end

    first_match = frag_bin_idx

    # ── Score matching bins ──────────────────────────────────────────
    @inbounds @fastmath begin
        while frag_bin_idx <= frag_bin_max_idx
            if fb_lows[frag_bin_idx] > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if fb_highs[frag_bin_idx] < frag_mz_min
                        break
                    end
                end
                frag_id_range = fb_first[frag_bin_idx]:fb_last[frag_bin_idx]
                searchFragmentBinUnconditional!(counter, fragments, frag_id_range)
                frag_bin_idx += one(UInt32)
            end
        end
    end

    # Return (first_match, ub) — lb advances for next peak
    return first_match, ub_final
end

"""
    _score_partition_hinted!(local_counter, partition, irt_low, irt_high, masses, mass_err_model;
                              linear_threshold=HINT_LINEAR_THRESHOLD)

Score one partition's fragments using hint-informed exponential search.
Uses queryFragmentHinted! which computes a smarter initial step_size from
per-bin density hints, then falls back to the proven exponential + binary search.
"""
@inline function _score_partition_hinted!(
        local_counter::LocalCounter{UInt16, UInt8},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel;
        linear_threshold::UInt32 = HINT_LINEAR_THRESHOLD,
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    p_rt_bins   = getRTBins(partition)
    p_frag_bins = getFragBins(partition)
    p_fragments = getFragments(partition)
    p_hints     = getSkipHints(partition)

    isempty(p_frag_bins) && return nothing
    n_rt = length(p_rt_bins)

    local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
    local_rt_bin_idx > n_rt && return nothing

    @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            lower_bound_guess = min_frag_bin
            upper_bound_guess = min_frag_bin
            prev_mz = Float32(0)

            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                lower_bound_guess, upper_bound_guess = queryFragmentHinted!(
                    local_counter, max_frag_bin,
                    lower_bound_guess, upper_bound_guess,
                    p_frag_bins, p_fragments,
                    frag_min, frag_max,
                    p_hints, prev_mz, linear_threshold)
                prev_mz = frag_min
            end
        end

        local_rt_bin_idx += 1
        if local_rt_bin_idx > n_rt
            break
        end
    end

    return nothing
end

"""
Partition-major search using hint-informed exponential search.
Same as searchFragmentIndexPartitionMajor but calls _score_partition_hinted!.
"""
function searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        all_scan_idxs::Vector{Int},
        n_threads::Int,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32};
        linear_threshold::UInt32 = HINT_LINEAR_THRESHOLD,
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel,
                 P<:Pioneer.FragmentIndexSearchParameters}

    min_score = Pioneer.getMinIndexSearchScore(params)
    iso_bounds = Pioneer.getIsotopeErrBounds(params)
    n_scans = length(all_scan_idxs)

    # ── Pre-compute per-scan properties ──────────────────────────────────────
    scan_irt_lo  = Vector{Float32}(undef, n_scans)
    scan_irt_hi  = Vector{Float32}(undef, n_scans)
    scan_prec_min = Vector{Float32}(undef, n_scans)
    scan_prec_max = Vector{Float32}(undef, n_scans)

    for (si, scan_idx) in enumerate(all_scan_idxs)
        irt_lo, irt_hi = Pioneer.getRTWindow(
            rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        quad_func = Pioneer.getQuadTransmissionFunction(
            qtm, Pioneer.getCenterMz(spectra, scan_idx),
            Pioneer.getIsolationWidthMz(spectra, scan_idx))
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        scan_irt_lo[si]   = irt_lo
        scan_irt_hi[si]   = irt_hi
        scan_prec_min[si] = prec_min
        scan_prec_max[si] = prec_max
    end

    # ── Per-thread result buffers and local counters ─────────────────────────
    thread_results = [Vector{Tuple{Int, UInt32}}(undef, 0) for _ in 1:n_threads]
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    thread_counters = [LocalCounter(UInt16, UInt8, max_local + 1) for _ in 1:n_threads]

    # ── Pre-compute partition → scan mapping ────────────────────────────────
    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end

    # ── Partition-major parallel execution (hinted) ──────────────────────────
    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            lc = thread_counters[tid]
            results = thread_results[tid]

            for k in 1:pfi.n_partitions
                relevant = partition_to_scans[k]
                partition = getPartition(pfi, k)
                n_relevant = length(relevant)
                (n_relevant == 0 || isempty(getFragBins(partition))) && continue

                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition_hinted!(lc, partition,
                        scan_irt_lo[si], scan_irt_hi[si],
                        Pioneer.getMzArray(spectra, scan_idx), mem;
                        linear_threshold=linear_threshold)

                    prec_lo = scan_prec_min[si]
                    prec_hi = scan_prec_max[si]
                    @inbounds for i in 1:(lc.size - 1)
                        lid = lc.ids[i]
                        score = lc.counts[lid]
                        score < min_score && continue
                        global_pid = l2g[lid]
                        pmz = precursor_mzs[global_pid]
                        (pmz < prec_lo || pmz > prec_hi) && continue
                        push!(results, (si, global_pid))
                    end

                    reset!(lc)
                    scan_i += n_threads
                end
            end
        end
    end
    fetch.(tasks)

    # ── Collect results into scan_to_prec_idx ────────────────────────────────
    scan_results = [UInt32[] for _ in 1:n_scans]
    for tid in 1:n_threads
        for (si, gpid) in thread_results[tid]
            push!(scan_results[si], gpid)
        end
        empty!(thread_results[tid])
    end

    precursors_passed = UInt32[]
    for si in 1:n_scans
        scan_idx = all_scan_idxs[si]
        pids = scan_results[si]
        if isempty(pids)
            scan_to_prec_idx[scan_idx] = missing
        else
            start = length(precursors_passed) + 1
            append!(precursors_passed, pids)
            scan_to_prec_idx[scan_idx] = start:length(precursors_passed)
        end
    end

    return precursors_passed
end
