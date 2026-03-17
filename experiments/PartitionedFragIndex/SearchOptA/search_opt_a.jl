#
# Approach A: Per-partition m/z lookup table
#
# Pre-build a coarse lookup table per partition mapping m/z → approximate
# frag bin index (as fraction × n_bins). Given a peak m/z, one reciprocal
# multiply gives the bucket, one array load gives the approximate bin index.
# Then a short linear scan finds the exact bin.
#
# No division at runtime — uses precomputed inv_bucket_width.
#

const LOOKUP_RESOLUTION = 0.05f0  # ~37K buckets for 150-2000 Da range
const LOOKUP_MZ_MIN = 100.0f0
const LOOKUP_MZ_MAX = 2100.0f0
const N_LOOKUP_BUCKETS = ceil(Int, (LOOKUP_MZ_MAX - LOOKUP_MZ_MIN) / LOOKUP_RESOLUTION)
const INV_LOOKUP_RESOLUTION = 1.0f0 / LOOKUP_RESOLUTION

"""
    MzLookupTable

Maps m/z → approximate frag bin index within an RT bin's range.
One table per partition (shared across RT bins since frag m/z distribution
is similar). Values are stored as Float32 fractions (0.0 = first bin,
1.0 = last bin) which are scaled to each RT bin's actual range.

Actually, store absolute frag bin indices for the median RT bin and use
as a starting guess for any RT bin. The error is bounded by how much
the bin layout varies across RT bins (usually small).
"""
struct MzLookupTable
    # For each bucket: the approximate frag bin INDEX (within this partition's
    # frag_bins array) where fragments near this m/z start. UInt32 index.
    bin_idx::Vector{UInt32}
end

"""
Build a lookup table for a partition by scanning its frag bins.
For each m/z bucket, store the index of the first frag bin whose
lb falls into or after that bucket.
"""
function build_mz_lookup(partition_frag_bins::Vector{Pioneer.FragIndexBin{Float32}})
    table = zeros(UInt32, N_LOOKUP_BUCKETS + 1)
    n_fb = length(partition_frag_bins)
    n_fb == 0 && return MzLookupTable(table)

    # Walk frag bins and fill buckets
    fb_idx = UInt32(1)
    for bucket in 1:(N_LOOKUP_BUCKETS + 1)
        bucket_mz = LOOKUP_MZ_MIN + (bucket - 1) * LOOKUP_RESOLUTION
        # Advance fb_idx to the first frag bin whose ub >= bucket_mz
        while fb_idx < n_fb && Pioneer.getHigh(partition_frag_bins[fb_idx]) < bucket_mz
            fb_idx += 1
        end
        table[bucket] = fb_idx
    end

    return MzLookupTable(table)
end

"""
Look up the approximate frag bin index for a given m/z.
Uses reciprocal multiply (no division).
"""
@inline function lookup_frag_bin(table::MzLookupTable, mz::Float32)
    bucket = floor(Int, (mz - LOOKUP_MZ_MIN) * INV_LOOKUP_RESOLUTION) + 1
    bucket = clamp(bucket, 1, N_LOOKUP_BUCKETS + 1)
    @inbounds return table.bin_idx[bucket]
end

"""
Query frag bins using lookup table for the starting position,
then linear scan forward/backward to find exact matching bins.
Replaces exponentialFragmentBinSearch + findFirstFragmentBin.
"""
function queryFragmentLookup!(
        prec_id_to_score::LocalCounter{UInt16, UInt8},
        frag_bin_max_idx::UInt32,
        frag_bin_min_idx::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{Float32}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32,
        lookup::MzLookupTable)

    # Get approximate starting position from lookup
    approx_idx = lookup_frag_bin(lookup, frag_mz_min)

    # Clamp to this RT bin's range
    approx_idx = clamp(approx_idx, frag_bin_min_idx, frag_bin_max_idx)

    # Linear scan backward to find the true first bin
    # (lookup might overshoot slightly)
    @inbounds while approx_idx > frag_bin_min_idx && Pioneer.getHigh(frag_bins[approx_idx - UInt32(1)]) >= frag_mz_min
        approx_idx -= UInt32(1)
    end

    # Now scan forward, scoring matching bins
    frag_bin_idx = approx_idx
    @inbounds @fastmath while frag_bin_idx <= frag_bin_max_idx
        frag_bin = frag_bins[frag_bin_idx]
        if Pioneer.getLow(frag_bin) > frag_mz_max
            break
        end
        if Pioneer.getHigh(frag_bin) >= frag_mz_min
            frag_id_range = Pioneer.getSubBinRange(frag_bin)
            searchFragmentBinUnconditional!(prec_id_to_score, fragments, frag_id_range)
        end
        frag_bin_idx += UInt32(1)
    end

    return nothing
end

"""
Score a partition using the lookup table approach.
"""
@inline function _score_partition_lookup!(
        local_counter::LocalCounter{UInt16, UInt8},
        partition,
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        lookup::MzLookupTable,
        ) where {U<:AbstractFloat}

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
            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                queryFragmentLookup!(
                    local_counter,
                    max_frag_bin,
                    min_frag_bin,
                    p_frag_bins,
                    p_fragments,
                    frag_min,
                    frag_max,
                    lookup
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
Partition-major search using lookup tables.
"""
function searchPartitionMajorLookup(
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

    # Build lookup tables per partition
    lookups = [build_mz_lookup(getFragBins(getPartition(pfi, k))) for k in 1:pfi.n_partitions]

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

                lookup = lookups[k]
                l2g = partition.local_to_global
                scan_i = tid
                while scan_i <= n_relevant
                    si = relevant[scan_i]
                    scan_idx = all_scan_idxs[si]

                    _score_partition_lookup!(lc, partition,
                        scan_irt_lo[si], scan_irt_hi[si],
                        Pioneer.getMzArray(spectra, scan_idx), mem, lookup)

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
