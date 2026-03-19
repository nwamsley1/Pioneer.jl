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
    _findFirstFragBin_hybrid(highs, lb, ub, frag_min, simd_cutoff) -> UInt32

Hybrid binary+SIMD search. Runs branchless binary search iterations until the
remaining range is ≤ simd_cutoff, then finishes with a SIMD linear scan.
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

# ── Scoring function ─────────────────────────────────────────────────────────

"""
    searchFragmentBinUnconditional!(counter, fragments, frag_id_range)

Score all fragments in the range unconditionally (no binary search on prec_mz).
Partitioning guarantees all fragments in the partition have prec_mz within
~partition_width of the query window.
"""
@inline function searchFragmentBinUnconditional!(
        prec_id_to_score::LocalCounter{I, UInt8},
        fragments::AbstractVector,
        frag_id_range::UnitRange{UInt32}) where {I}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        inc!(prec_id_to_score, getPrecID(frag), getScore(frag))
    end
    return nothing
end

# ── Hinted search functions ─────────────────────────────────────────────────

"""
    queryFragmentHinted!(counter, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
                          frag_bins, fragments, frag_mz_min, frag_mz_max,
                          hints, prev_mz, linear_threshold)

5-Da direct hint search with advancing lb, SoA layout + SIMD linear scan.

Algorithm:
1. Hint-based lb advancement (provably safe within 5 Da).
2. Hint-based UB guess: new_lb + est_step * 1.5.
3. Exponential doubling only if UB guess insufficient (~10% of cases).
4. Hybrid binary→SIMD search for first matching bin.
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

    # ── Lower bound: previous first_match (always safe, peaks sorted by m/z) ──
    new_lb = lower_bound_guess

    # ── Hint-based UB guess (safe: exponential doubling corrects if too low) ──
    est_step = one(UInt32)
    if prev_mz > 0.0f0 && lower_bound_guess <= frag_bin_max_idx
        delta_mz = frag_mz_min - @inbounds fb_lows[lower_bound_guess]
        if delta_mz > 0.0f0
            hint = @inbounds hints[lower_bound_guess]
            est_step = max(one(UInt32), unsafe_trunc(UInt32, Float32(hint) * (delta_mz / 5.0f0)))
        end
    end
    ub_guess = min(new_lb + est_step + est_step, frag_bin_max_idx)

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
        frag_bin_idx = _find_first_ge(fb_highs, new_lb, ub_final, frag_mz_min)
    else
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

    return first_match, ub_final
end

"""
    _score_partition_hinted!(local_counter, partition, irt_low, irt_high, masses, mass_err_model)

Score one partition's fragments using hint-informed search with SIMD acceleration.
"""
@inline function _score_partition_hinted!(
        local_counter::LocalCounter{UInt16, UInt8},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::MassErrorModel;
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

    @inbounds @fastmath while getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            lower_bound_guess = min_frag_bin
            upper_bound_guess = min_frag_bin
            prev_mz = Float32(0)

            for mass in masses
                corrected_mz = getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = getMzBoundsReverse(mass_err_model, corrected_mz)

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
Binary search for the first RT bin whose upper bound (getHigh) >= irt_low.
"""
@inline function _find_rt_bin_start(rt_bins::Vector{FragIndexBin{T}}, irt_low::Float32) where {T}
    lo, hi = 1, length(rt_bins)
    result = hi + 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if getHigh(rt_bins[mid]) >= irt_low
            result = mid
            hi = mid - 1
        else
            lo = mid + 1
        end
    end
    return result
end

# ── Partition-major hinted search (orchestrator) ─────────────────────────────

"""
    searchFragmentIndexPartitionMajorHinted(scan_to_prec_idx, pfi, spectra,
        all_scan_idxs, n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol,
        precursor_mzs)

Partition-major search: outer loop over partitions, inner loop fans MS2 scans
across threads. All threads working at any given time read the same partition's
fragment/bin arrays, maximizing shared cache utilization.

Returns a flat vector of global precursor IDs (concatenated across all scans).
"""
function searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::MassSpecData,
        all_scan_idxs::Vector{Int},
        n_threads::Int,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32};
        linear_threshold::UInt32 = HINT_LINEAR_THRESHOLD,
        ) where {M<:MassErrorModel, Q<:QuadTransmissionModel,
                 P<:FragmentIndexSearchParameters}

    min_score = getMinIndexSearchScore(params)
    iso_bounds = getIsotopeErrBounds(params)
    n_scans = length(all_scan_idxs)

    # ── Pre-compute per-scan properties ──────────────────────────────────────
    scan_irt_lo  = Vector{Float32}(undef, n_scans)
    scan_irt_hi  = Vector{Float32}(undef, n_scans)
    scan_prec_min = Vector{Float32}(undef, n_scans)
    scan_prec_max = Vector{Float32}(undef, n_scans)

    for (si, scan_idx) in enumerate(all_scan_idxs)
        irt_lo, irt_hi = getRTWindow(
            rt_to_irt_spline(getRetentionTime(spectra, scan_idx)), irt_tol)
        quad_func = getQuadTransmissionFunction(
            qtm, getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx))
        prec_min = Float32(getPrecMinBound(quad_func) - NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(getPrecMaxBound(quad_func) + NEUTRON * last(iso_bounds) / 2)
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
                        getMzArray(spectra, scan_idx), mem;
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
