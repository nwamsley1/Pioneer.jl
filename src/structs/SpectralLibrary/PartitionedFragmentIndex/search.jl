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
        prec_id_to_score::Counter{I, UInt8},
        fragments::AbstractVector,
        frag_id_range::UnitRange{UInt32}) where {I}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        or!(prec_id_to_score, getPrecID(frag), getScore(frag))
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
    _score_partition_hinted!(local_counter, partition, irt_low, irt_high,
                             masses, intensities, mass_err_model)

Score one partition's fragments using hint-informed search with SIMD acceleration.
When `mass_err_model` is an IntensityMassErrorModel, uses per-peak intensity-aware
bias correction and tolerance windows via the 3-arg getCorrectedMz / getMzBoundsReverse.
"""
@inline function _score_partition_hinted!(
        local_counter::Counter{UInt16, UInt8},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        intensities::AbstractArray{Union{Missing, V}},
        mass_err_model::AbstractMassErrorModel;
        linear_threshold::UInt32 = HINT_LINEAR_THRESHOLD,
        intensity_threshold::Float32 = 0.0f0,
        ) where {T<:AbstractFloat, U<:AbstractFloat, V<:AbstractFloat}

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
            n_peaks = length(masses)

            for peak_i in 1:n_peaks
                @inbounds int_f32 = intensities[peak_i]::Float32
                int_f32 < intensity_threshold && continue
                @inbounds mass_f32 = masses[peak_i]::Float32
                corrected_mz, frag_min, frag_max = getCorrectedMzAndBounds(mass_err_model, mass_f32, int_f32)

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

# ── Emit strategies for fragment index scoring ─────────────────────────────
#
# Two modes for processing the counter after scoring a scan's fragments:
#
# EmitToBuffer: MainSearch path. Writes (scan_idx, precursor_id) pairs to flat
#   buffers. Only applies the bitmask LUT filter — precise m/z and iRT filtering
#   is deferred to selectTransitions! in the deconvolution step.
#
# EmitToAccumulator: BitVecCalibration path. Tallies 256-bin target/decoy pattern
#   counts directly. Applies prec_mz AND iRT filtering because this is the final
#   consumer (no selectTransitions! downstream).

abstract type FragIndexEmitStrategy end

"""
Emit to flat buffers for downstream deconvolution. Only bitmask filter applied.
Parameterized on filter type F for compile-time dispatch.
"""
struct EmitToBuffer{F<:AbstractBitVecFilter} <: FragIndexEmitStrategy
    sf::F
end

"""
Emit to pattern accumulator for BitVecCalibration. Applies prec_mz + iRT filtering.
"""
struct EmitToAccumulator{A<:PatternAccumulator, V<:AbstractVector{Float32}} <: FragIndexEmitStrategy
    acc::A
    precursor_irts::V               # concrete element/array type: avoids per-candidate boxing in emit_candidates!
    irt_tol::Float32
end

"""
    emit_candidates!(strategy, lc, l2g, si, scan_irt, prec_lo, prec_hi,
                      precursor_mzs, tid, si_buf, pid_buf, wp) -> wp

Process the scored counter for one scan. Dispatches on strategy type.
Returns updated write position (for EmitToBuffer) or 0 (for EmitToAccumulator).
"""
@inline function emit_candidates!(s::EmitToBuffer{F}, lc::Counter{UInt16, UInt8},
        l2g::Vector{UInt32}, si::Int, scan_irt::Float32,
        prec_lo::Float32, prec_hi::Float32,
        precursor_mzs::AbstractVector{Float32},
        tid::Int, si_buf::Vector{Int32}, pid_buf::Vector{UInt32}, wp::Int) where {F<:AbstractBitVecFilter}
    sf = s.sf
    @inbounds for i in 1:(lc.size - 1)
        lid = lc.ids[i]
        score = lc.counts[lid]
        passes_filter(sf, score) || continue
        global_pid = l2g[lid]
        wp += 1
        if wp > length(pid_buf)
            new_len = length(pid_buf) * 2
            resize!(si_buf, new_len)
            resize!(pid_buf, new_len)
        end
        si_buf[wp] = Int32(si)
        pid_buf[wp] = global_pid
    end
    return wp
end

@inline function emit_candidates!(s::EmitToAccumulator, lc::Counter{UInt16, UInt8},
        l2g::Vector{UInt32}, si::Int, scan_irt::Float32,
        prec_lo::Float32, prec_hi::Float32,
        precursor_mzs::AbstractVector{Float32},
        tid::Int, si_buf::Vector{Int32}, pid_buf::Vector{UInt32}, wp::Int)
    acc = s.acc
    acc_min = acc.min_score
    irt_tol = s.irt_tol
    prec_irts = s.precursor_irts
    @inbounds for i in 1:(lc.size - 1)
        lid = lc.ids[i]
        score = lc.counts[lid]
        count_ones(score) < acc_min && continue
        global_pid = l2g[lid]
        # Prec m/z filter (partition edges)
        pmz = precursor_mzs[global_pid]
        (pmz < prec_lo || pmz > prec_hi) && continue
        # iRT filter (precise per-precursor check)
        abs(prec_irts[global_pid] - scan_irt) > irt_tol && continue
        accumulate!(acc, tid, global_pid, score)
    end
    return wp  # unchanged for accumulator mode
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
        score_filter::AbstractBitVecFilter = CountFilter(getMinIndexSearchScore(params)),
        pattern_accumulator::Union{Nothing, PatternAccumulator} = nothing,
        precursor_irts::AbstractVector{Float32} = Float32[],
        max_peaks::Int = 0,
        scratch::Union{Nothing, FragIndexScratch} = nothing,
        iso_bounds_override::Union{Nothing, Tuple{UInt8, UInt8}} = nothing,
        ) where {M<:AbstractMassErrorModel, Q<:QuadTransmissionModel,
                 P<:FragmentIndexSearchParameters}

    # `iso_bounds_override` lets a caller narrow the precursor-selection window
    # independently of the params (used by Sciex ZT main search: the center-bin
    # filter deletes isotope-error catches anyway, so the (1,0) down-extension is
    # wasted frag-index work — override to (0,0)).
    iso_bounds = iso_bounds_override === nothing ? getIsotopeErrBounds(params) : iso_bounds_override
    n_scans = length(all_scan_idxs)

    # ── 1. Pre-compute per-scan properties ─────────────────────────────────
    scan_irt_lo, scan_irt_hi, scan_prec_min, scan_prec_max, scan_irts =
        _precompute_scan_properties(spectra, all_scan_idxs, rt_to_irt_spline, irt_tol, qtm, iso_bounds, n_scans)

    # ── 2. Map partitions to scans ─────────────────────────────────────────
    partition_to_scans = _build_partition_scan_mapping(pfi, scan_prec_min, scan_prec_max, n_scans)

    # ── 3. Per-thread buffers (reused via FragIndexScratch when supplied) ──
    # Per-thread initial guess: n_scans*200 estimates TOTAL emitted candidates,
    # but each thread only handles a fraction, so divide by n_threads. The emit
    # buffers grow in place (doubling) if a thread exceeds this, and FragIndexScratch
    # is grow-only across files, so a low first-file guess self-corrects to the
    # true high-water-mark. Avoids ~3-20x per-thread over-provisioning (worst on SCP).
    est_per_thread = max(div(n_scans * 200, n_threads), 100_000)
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    int_buf_size = max_peaks > 0 ?
        maximum(si -> length(getMzArray(spectra, all_scan_idxs[si])),
                1:n_scans; init=0) : 0

    if scratch === nothing
        thread_si_bufs  = [Vector{Int32}(undef, est_per_thread) for _ in 1:n_threads]
        thread_pid_bufs = [Vector{UInt32}(undef, est_per_thread) for _ in 1:n_threads]
        thread_counters = [Counter(UInt16, UInt8, max_local + 1) for _ in 1:n_threads]
        thread_int_bufs = max_peaks > 0 ?
            [Vector{Float32}(undef, int_buf_size) for _ in 1:n_threads] :
            [Float32[] for _ in 1:n_threads]
    else
        prepare!(scratch;
            n_threads = n_threads,
            est_per_thread = est_per_thread,
            counter_size = max_local + 1,
            int_buf_size = int_buf_size)
        thread_si_bufs  = scratch.si
        thread_pid_bufs = scratch.pid
        thread_counters = scratch.counters
        thread_int_bufs = scratch.int_bufs
    end
    thread_counts = zeros(Int, n_threads)

    # ── 4. Build emit strategy (compile-time dispatch) ─────────────────────
    emit_strategy = if pattern_accumulator !== nothing
        EmitToAccumulator(pattern_accumulator, precursor_irts, Float32(irt_tol))
    else
        EmitToBuffer(score_filter)
    end

    # ── 6. Partition-major parallel execution ──────────────────────────────
    thread_times = zeros(Float64, n_threads)
    t_parallel_start = time()
    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            t_thread = time()
            thread_counts[tid] = _run_thread(tid, emit_strategy,
                thread_counters[tid],
                thread_si_bufs[tid], thread_pid_bufs[tid],
                pfi, partition_to_scans, all_scan_idxs, spectra,
                scan_irt_lo, scan_irt_hi, scan_prec_min, scan_prec_max,
                scan_irts, precursor_mzs,
                mem, linear_threshold, n_threads,
                thread_int_bufs[tid], max_peaks)
            thread_times[tid] = time() - t_thread
        end
    end
    fetch.(tasks)
    t_parallel = time() - t_parallel_start

    # ── 6. Collect results via counting sort ───────────────────────────────
    precursors_passed, total = _collect_frag_index_results!(
        scan_to_prec_idx, thread_si_bufs, thread_pid_bufs, thread_counts,
        all_scan_idxs, n_scans, n_threads)

    if params isa MainSearchParameters
        @debug_l1 "  frag_index: parallel=$(round(t_parallel, digits=2))s, results=$total"
    end

    scores_passed = UInt8[]
    return precursors_passed, scores_passed
end

# ── Helper functions ─────────────────────────────────────────────────────────

"""Pre-compute per-scan iRT windows, prec m/z bounds, and center iRTs."""
function _precompute_scan_properties(spectra, all_scan_idxs, rt_to_irt_spline,
        irt_tol, qtm, iso_bounds, n_scans)
    scan_irt_lo  = Vector{Float32}(undef, n_scans)
    scan_irt_hi  = Vector{Float32}(undef, n_scans)
    scan_prec_min = Vector{Float32}(undef, n_scans)
    scan_prec_max = Vector{Float32}(undef, n_scans)
    scan_irts    = Vector{Float32}(undef, n_scans)

    for (si, scan_idx) in enumerate(all_scan_idxs)
        scan_irt = Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx)))
        irt_lo, irt_hi = getRTWindow(scan_irt, irt_tol)
        quad_func = getQuadTransmissionFunction(
            qtm, getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx))
        scan_irt_lo[si] = irt_lo
        scan_irt_hi[si] = irt_hi
        scan_irts[si] = scan_irt
        scan_prec_min[si] = Float32(getPrecMinBound(quad_func) - C13_C12_MASS_DIFF * first(iso_bounds) / 2)
        scan_prec_max[si] = Float32(getPrecMaxBound(quad_func) + C13_C12_MASS_DIFF * last(iso_bounds) / 2)
    end
    return scan_irt_lo, scan_irt_hi, scan_prec_min, scan_prec_max, scan_irts
end

"""Build partition → scan mapping for partition-major traversal."""
function _build_partition_scan_mapping(pfi, scan_prec_min, scan_prec_max, n_scans)
    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end
    return partition_to_scans
end

"""Typed inner function: Julia specializes on E (emit strategy) + M (mass error model)."""
function _run_thread(tid::Int, emit::E,
        lc::Counter{UInt16, UInt8},
        si_buf::Vector{Int32}, pid_buf::Vector{UInt32},
        pfi_ref::LocalPartitionedFragmentIndex{Float32},
        partition_to_scans::Vector{Vector{Int}},
        all_scan_idxs::Vector{Int},
        spectra::MassSpecData,
        scan_irt_lo::Vector{Float32}, scan_irt_hi::Vector{Float32},
        scan_prec_min::Vector{Float32}, scan_prec_max::Vector{Float32},
        scan_irts::Vector{Float32},
        precursor_mzs::AbstractVector{Float32},
        mem::M,
        linear_threshold::UInt32, n_threads::Int,
        int_buf::Vector{Float32}, max_peaks::Int
        ) where {E<:FragIndexEmitStrategy, M<:AbstractMassErrorModel}

    wp = 0

    for k in 1:pfi_ref.n_partitions
        relevant = partition_to_scans[k]
        partition = getPartition(pfi_ref, k)
        n_relevant = length(relevant)
        (n_relevant == 0 || isempty(getFragBins(partition))) && continue

        l2g = partition.local_to_global
        scan_i = tid
        while scan_i <= n_relevant
            si = relevant[scan_i]
            scan_idx = all_scan_idxs[si]

            # Compute intensity threshold for top-N peak filtering (wide scout only)
            intensity_threshold = 0.0f0
            if max_peaks > 0
                scan_intensities = getIntensityArray(spectra, scan_idx)
                n_scan_peaks = length(scan_intensities)
                if n_scan_peaks > max_peaks
                    @inbounds for j in 1:n_scan_peaks
                        int_buf[j] = scan_intensities[j]::Float32
                    end
                    k = n_scan_peaks - max_peaks + 1
                    partialsort!(view(int_buf, 1:n_scan_peaks), k)
                    intensity_threshold = int_buf[k]
                end
            end

            _score_partition_hinted!(lc, partition,
                scan_irt_lo[si], scan_irt_hi[si],
                getMzArray(spectra, scan_idx),
                getIntensityArray(spectra, scan_idx),
                mem;
                linear_threshold=linear_threshold,
                intensity_threshold=intensity_threshold)

            wp = emit_candidates!(emit, lc, l2g, si, scan_irts[si],
                scan_prec_min[si], scan_prec_max[si],
                precursor_mzs, tid, si_buf, pid_buf, wp)

            reset!(lc)
            scan_i += n_threads
        end
    end
    return wp
end

"""Counting-sort merge of per-thread flat buffers into scan-indexed output."""
function _collect_frag_index_results!(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        thread_si_bufs::Vector{Vector{Int32}},
        thread_pid_bufs::Vector{Vector{UInt32}},
        thread_counts::Vector{Int},
        all_scan_idxs::Vector{Int},
        n_scans::Int, n_threads::Int)

    # 1. Count results per scan
    scan_counts = zeros(Int, n_scans)
    for tid in 1:n_threads
        si_buf = thread_si_bufs[tid]
        @inbounds for j in 1:thread_counts[tid]
            scan_counts[si_buf[j]] += 1
        end
    end

    # 2. Prefix sum for offsets
    total = sum(scan_counts)
    offsets = Vector{Int}(undef, n_scans)
    offsets[1] = 1
    @inbounds for si in 2:n_scans
        offsets[si] = offsets[si-1] + scan_counts[si-1]
    end
    write_pos = copy(offsets)

    # 3. Scatter into pre-allocated output
    precursors_passed = Vector{UInt32}(undef, total)
    for tid in 1:n_threads
        si_buf = thread_si_bufs[tid]
        pid_buf = thread_pid_bufs[tid]
        @inbounds for j in 1:thread_counts[tid]
            si = si_buf[j]
            precursors_passed[write_pos[si]] = pid_buf[j]
            write_pos[si] += 1
        end
    end

    # 4. Build scan_to_prec_idx
    @inbounds for si in 1:n_scans
        if scan_counts[si] == 0
            scan_to_prec_idx[all_scan_idxs[si]] = missing
        else
            scan_to_prec_idx[all_scan_idxs[si]] = offsets[si]:(offsets[si] + scan_counts[si] - 1)
        end
    end

    return precursors_passed, total
end
