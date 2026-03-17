#
# Approach B: Delta-based skip with adaptive linear/jump scan
#
# Since peaks are sorted by m/z, use the delta between consecutive peaks
# to estimate how many frag bins to skip. Multiply delta_mz by a precomputed
# inv_avg_bin_width (one multiply, no division).
#
# - Small skip (≤ threshold): linear scan forward (cache-friendly, SIMD-able)
# - Large skip: jump ahead by estimated amount, then short linear refine
#
# No exponential search, no binary search.
#

const LINEAR_SCAN_THRESHOLD = 32  # switch to jump if estimated skip > this

"""
    DeltaSearchState

Mutable state carried across consecutive peaks within an RT bin.
Tracks the current frag bin position so the next peak can estimate
how far to jump.
"""
mutable struct DeltaSearchState
    current_idx::UInt32      # current frag bin position
    current_mz::Float32      # m/z of current position (for delta calc)
end

"""
Query frag bins using delta-based skip from previous peak position.
No exponential search, no binary search — just estimate + linear scan.
"""
@inline function queryFragmentDelta!(
        prec_id_to_score::LocalCounter{UInt16, UInt8},
        frag_bin_max_idx::UInt32,
        frag_bin_min_idx::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{Float32}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32,
        state::DeltaSearchState,
        inv_avg_bin_width::Float32)

    # Estimate how many bins to skip based on m/z delta from last position
    delta_mz = frag_mz_min - state.current_mz
    est_skip = floor(Int32, delta_mz * inv_avg_bin_width)

    # Compute starting position
    if est_skip <= Int32(0)
        # Going backward or same position — start from current
        start_idx = state.current_idx
    elseif est_skip <= LINEAR_SCAN_THRESHOLD
        # Small jump — just start from current and linear scan will cover it
        start_idx = state.current_idx
    else
        # Large jump — skip ahead by estimated amount (with safety margin)
        jump_target = state.current_idx + UInt32(max(0, est_skip - 2))
        start_idx = min(jump_target, frag_bin_max_idx)
    end

    # Clamp to RT bin range
    start_idx = clamp(start_idx, frag_bin_min_idx, frag_bin_max_idx)

    # Linear scan backward if we overshot
    @inbounds while start_idx > frag_bin_min_idx && Pioneer.getLow(frag_bins[start_idx]) > frag_mz_max
        start_idx -= UInt32(1)
    end

    # Linear scan forward to find first matching bin
    @inbounds while start_idx <= frag_bin_max_idx && Pioneer.getHigh(frag_bins[start_idx]) < frag_mz_min
        start_idx += UInt32(1)
    end

    # Score matching bins
    frag_bin_idx = start_idx
    @inbounds @fastmath while frag_bin_idx <= frag_bin_max_idx
        frag_bin = frag_bins[frag_bin_idx]
        if Pioneer.getLow(frag_bin) > frag_mz_max
            break
        end
        frag_id_range = Pioneer.getSubBinRange(frag_bin)
        searchFragmentBinUnconditional!(prec_id_to_score, fragments, frag_id_range)
        frag_bin_idx += UInt32(1)
    end

    # Update state for next peak
    state.current_idx = frag_bin_idx > frag_bin_max_idx ? frag_bin_max_idx : frag_bin_idx
    state.current_mz = frag_mz_min

    return nothing
end

"""
Score a partition using delta-based skip.
"""
@inline function _score_partition_delta!(
        local_counter::LocalCounter{UInt16, UInt8},
        partition,
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        inv_avg_bin_width::Float32,
        ) where {U<:AbstractFloat}

    p_rt_bins   = getRTBins(partition)
    p_frag_bins = getFragBins(partition)
    p_fragments = getFragments(partition)

    isempty(p_frag_bins) && return nothing
    n_rt = length(p_rt_bins)

    local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
    local_rt_bin_idx > n_rt && return nothing

    state = DeltaSearchState(UInt32(1), 0.0f0)

    @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            # Reset state for each RT bin
            state.current_idx = min_frag_bin
            state.current_mz = Pioneer.getLow(p_frag_bins[min_frag_bin])

            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                queryFragmentDelta!(
                    local_counter,
                    max_frag_bin,
                    min_frag_bin,
                    p_frag_bins,
                    p_fragments,
                    frag_min,
                    frag_max,
                    state,
                    inv_avg_bin_width
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
Compute inverse average bin width for a partition (for delta estimation).
"""
function compute_inv_avg_bin_width(partition)
    fb = getFragBins(partition)
    isempty(fb) && return 1.0f0
    total_range = Pioneer.getHigh(fb[end]) - Pioneer.getLow(fb[1])
    total_range <= 0 && return 1.0f0
    avg_width = total_range / length(fb)
    return 1.0f0 / avg_width
end

"""
Partition-major search using delta-based skip.
"""
function searchPartitionMajorDelta(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::LocalPartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        all_scan_idxs::Vector{Int},
        n_threads::Int,
        min_score::UInt8,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32};
        iso_bounds_val::Tuple{UInt8, UInt8} = (UInt8(1), UInt8(0))
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel}
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
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds_val) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds_val) / 2)
        scan_irt_lo[si]   = irt_lo
        scan_irt_hi[si]   = irt_hi
        scan_prec_min[si] = prec_min
        scan_prec_max[si] = prec_max
    end

    # Precompute inv_avg_bin_width per partition
    inv_bin_widths = Float32[compute_inv_avg_bin_width(getPartition(pfi, k)) for k in 1:pfi.n_partitions]

    thread_results = [Vector{Tuple{Int, UInt32}}(undef, 0) for _ in 1:n_threads]
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    thread_counters = [LocalCounter(UInt16, UInt8, max_local + 1) for _ in 1:n_threads]

    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end

    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            lc = thread_counters[tid]
            results = thread_results[tid]

            for k in 1:pfi.n_partitions
                relevant = partition_to_scans[k]
                partition = getPartition(pfi, k)
                n_relevant = length(relevant)
                (n_relevant == 0 || isempty(getFragBins(partition))) && continue

                inv_bw = inv_bin_widths[k]
                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition_delta!(lc, partition,
                        scan_irt_lo[si], scan_irt_hi[si],
                        Pioneer.getMzArray(spectra, scan_idx), mem, inv_bw)

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
