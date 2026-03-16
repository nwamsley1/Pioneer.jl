"""
    build_partitioned_index(nfi::NativeFragmentIndex{T}; partition_width=5.0f0)

Reorganize a NativeFragmentIndex into a PartitionedFragmentIndex by splitting fragments
according to their precursor m/z into partitions of width `partition_width` Da.

Every partition gets the full set of RT bins (same iRT boundaries as the original),
so partition-local searches can reuse the same rt_bin_idx advancement logic.
"""
function build_partitioned_index(nfi::Pioneer.NativeFragmentIndex{T};
                                  partition_width::T = T(5.0)) where {T<:AbstractFloat}
    rt_bins_orig  = Pioneer.getRTBins(nfi)
    frag_bins_orig = Pioneer.getFragBins(nfi)
    fragments_orig = Pioneer.getFragments(nfi)
    n_rt_bins = length(rt_bins_orig)

    # Step 1: Find prec_mz range across all fragments
    min_prec_mz = typemax(T)
    max_prec_mz = typemin(T)
    for frag in fragments_orig
        pmz = Pioneer.getPrecMZ(frag)
        min_prec_mz = min(min_prec_mz, pmz)
        max_prec_mz = max(max_prec_mz, pmz)
    end
    n_partitions = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / partition_width))
    println("  Prec m/z range: [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))]")
    println("  Partition width: $partition_width Da → $n_partitions partitions")

    # Step 2: Collect fragments per partition, preserving RT-bin and frag-bin structure.
    #
    # For each partition k, we build:
    #   partition_fragments[k]  — Vector{IndexFragment{T}}
    #   partition_frag_bins[k]  — Vector{FragIndexBin{T}}  (within each RT bin)
    #   partition_rt_bins[k]    — Vector{FragIndexBin{T}}
    #
    # We walk the original hierarchy: RT bin → frag bin → fragment.

    # Pre-allocate per-partition fragment accumulators
    partition_fragments = [Pioneer.IndexFragment{T}[] for _ in 1:n_partitions]

    # We'll build frag_bins and rt_bins in a second pass after collecting fragments
    # per (partition, rt_bin, frag_bin).
    #
    # Strategy: iterate RT bins → frag bins → fragments.
    # For each fragment, assign to partition. Record (rt_bin_idx, frag_bin_idx, partition_k)
    # groupings so we can reconstruct per-partition frag_bin boundaries.

    # Per-partition, per-RT-bin: list of (frag_bin_lb, frag_bin_ub, fragment_indices_in_partition)
    # We'll use a simpler approach: for each partition, accumulate fragments grouped by
    # (rt_bin_idx, frag_bin_idx) and then build bins from those groups.

    # Data structure: partition → rt_bin → Vector of (frag_bin_lb, frag_bin_ub, frag_start, frag_count)
    # where frag_start/count index into partition_fragments[k].

    FragBinInfo = @NamedTuple{lb::T, ub::T, start::UInt32, count::UInt32}

    # partition → rt_bin → Vector{FragBinInfo}
    part_rt_fbin = [[FragBinInfo[] for _ in 1:n_rt_bins] for _ in 1:n_partitions]

    # Walk the hierarchy
    for rt_idx in 1:n_rt_bins
        rt_bin = rt_bins_orig[rt_idx]
        frag_bin_range = Pioneer.getSubBinRange(rt_bin)
        for fb_idx in frag_bin_range
            fbin = frag_bins_orig[fb_idx]
            frag_range = Pioneer.getSubBinRange(fbin)
            fb_lb = Pioneer.getLow(fbin)
            fb_ub = Pioneer.getHigh(fbin)

            # Temporarily collect per-partition fragments for this frag_bin
            # Use a dictionary keyed by partition index
            per_part_start = Dict{Int, UInt32}()
            per_part_count = Dict{Int, UInt32}()

            for frag_idx in frag_range
                frag = fragments_orig[frag_idx]
                pmz = Pioneer.getPrecMZ(frag)
                k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1, 1, n_partitions)

                if !haskey(per_part_start, k)
                    per_part_start[k] = UInt32(length(partition_fragments[k]) + 1)
                    per_part_count[k] = UInt32(0)
                end
                push!(partition_fragments[k], frag)
                per_part_count[k] += UInt32(1)
            end

            # Record frag bin info for each partition that got fragments
            for (k, start) in per_part_start
                push!(part_rt_fbin[k][rt_idx],
                      (lb=fb_lb, ub=fb_ub, start=start, count=per_part_count[k]))
            end
        end
    end

    # Step 3: Build NativeFragmentIndex for each partition
    partitions = Vector{Pioneer.NativeFragmentIndex{T}}(undef, n_partitions)

    for k in 1:n_partitions
        frags = partition_fragments[k]

        # Build frag_bins and rt_bins for this partition
        p_frag_bins = Pioneer.FragIndexBin{T}[]
        p_rt_bins   = Pioneer.FragIndexBin{T}[]

        for rt_idx in 1:n_rt_bins
            orig_rt = rt_bins_orig[rt_idx]
            rt_lb = Pioneer.getLow(orig_rt)
            rt_ub = Pioneer.getHigh(orig_rt)

            fb_infos = part_rt_fbin[k][rt_idx]

            if isempty(fb_infos)
                # Empty RT bin: point to an empty frag_bin range
                fb_start = UInt32(length(p_frag_bins) + 1)
                fb_end   = UInt32(length(p_frag_bins))  # empty range
            else
                fb_start = UInt32(length(p_frag_bins) + 1)
                for info in fb_infos
                    push!(p_frag_bins, Pioneer.FragIndexBin{T}(
                        info.lb, info.ub, info.start, info.start + info.count - UInt32(1)
                    ))
                end
                fb_end = UInt32(length(p_frag_bins))
            end

            push!(p_rt_bins, Pioneer.FragIndexBin{T}(rt_lb, rt_ub, fb_start, fb_end))
        end

        partitions[k] = Pioneer.NativeFragmentIndex{T}(p_frag_bins, p_rt_bins, frags)
    end

    return PartitionedFragmentIndex{T}(partitions, min_prec_mz, partition_width, n_partitions)
end
