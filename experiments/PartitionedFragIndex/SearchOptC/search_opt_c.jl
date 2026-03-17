#
# O(1) fragment lookup using fixed-width m/z bins.
#
# For each peak, compute the bin index with one multiply, load the fragment
# range from the CSR pointer array, and score all fragments. No search.
#

"""
Score all fragments in a fixed m/z bin range into the local counter.
"""
@inline function _score_fixed_bins!(
        lc::LocalCounter{UInt16, UInt8},
        rt_data::FixedBinRTData,
        bin_lo::Int32,
        bin_hi::Int32,
        n_mz_bins::Int32)

    bin_lo = max(bin_lo, Int32(1))
    bin_hi = min(bin_hi, n_mz_bins)
    bin_lo > bin_hi && return nothing

    frag_ptrs = rt_data.frag_ptrs
    fragments = rt_data.fragments

    @inbounds for bin in bin_lo:bin_hi
        frag_start = UInt32(frag_ptrs[bin])
        frag_end = UInt32(frag_ptrs[bin + 1]) - UInt32(1)
        for fi in frag_start:frag_end
            frag = fragments[fi]
            inc!(lc, frag.local_id, frag.score)
        end
    end

    return nothing
end

"""
Score one partition's fragments using fixed-width bin O(1) lookup.
"""
@inline function _score_partition_fixed!(
        lc::LocalCounter{UInt16, UInt8},
        partition::FixedBinPartition,
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        ) where {U<:AbstractFloat}

    rt_data_list = partition.rt_data
    isempty(rt_data_list) && return nothing

    mz_min = partition.mz_min
    n_mz_bins = partition.n_mz_bins
    n_mz_bins <= 0 && return nothing

    # Find first RT bin with irt_hi >= irt_low
    n_rt = length(rt_data_list)
    rt_idx = 1
    @inbounds while rt_idx <= n_rt && rt_data_list[rt_idx].irt_hi < irt_low
        rt_idx += 1
    end

    @inbounds while rt_idx <= n_rt && rt_data_list[rt_idx].irt_lo < irt_high
        rd = rt_data_list[rt_idx]

        for mass in masses
            corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
            frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

            # O(1) bin lookup — one multiply each, no search
            bin_lo = floor(Int32, (frag_min - mz_min) * INV_FIXED_BIN_WIDTH) + Int32(1)
            bin_hi = floor(Int32, (frag_max - mz_min) * INV_FIXED_BIN_WIDTH) + Int32(1)

            _score_fixed_bins!(lc, rd, bin_lo, bin_hi, n_mz_bins)
        end

        rt_idx += 1
    end

    return nothing
end

"""
Partition-major search using fixed-width bin O(1) lookup.
"""
function searchPartitionMajorFixed(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::FixedBinPartitionedIndex,
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

    # Pre-compute per-scan properties
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

    # Per-thread resources
    thread_results = [Vector{Tuple{Int, UInt32}}(undef, 0) for _ in 1:n_threads]
    max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pfi); init=0)
    thread_counters = [LocalCounter(UInt16, UInt8, max_local + 1) for _ in 1:n_threads]

    # Pre-compute partition → scan mapping
    partition_to_scans = [Int[] for _ in 1:pfi.n_partitions]
    for si in 1:n_scans
        first_k, last_k = get_partition_range(pfi, scan_prec_min[si], scan_prec_max[si])
        for k in first_k:last_k
            push!(partition_to_scans[k], si)
        end
    end

    # Partition-major parallel execution
    tasks = map(1:n_threads) do tid
        Threads.@spawn begin
            lc = thread_counters[tid]
            results = thread_results[tid]

            for k in 1:pfi.n_partitions
                relevant = partition_to_scans[k]
                partition = getPartition(pfi, k)
                n_relevant = length(relevant)
                (n_relevant == 0 || isempty(partition.rt_data)) && continue

                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition_fixed!(lc, partition,
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

    # Collect results
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
