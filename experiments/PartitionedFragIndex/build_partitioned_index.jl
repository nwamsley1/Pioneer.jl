"""
    build_partitioned_index(nfi::NativeFragmentIndex{T}; partition_width=5.0f0)

Reorganize a NativeFragmentIndex into a PartitionedFragmentIndex by splitting fragments
according to their precursor m/z into partitions of width `partition_width` Da.

Every partition preserves the **full** frag bin and RT bin structure from the original
index. Frag bins that have no fragments for a given partition get empty ranges
(first > last), so the exponential + binary search navigates the identical "landscape"
as the baseline. RT bins are copied verbatim since they index into frag bins at the
same positions.
"""
function build_partitioned_index(nfi::Pioneer.NativeFragmentIndex{T};
                                  partition_width::T = T(5.0)) where {T<:AbstractFloat}
    rt_bins_orig   = Pioneer.getRTBins(nfi)
    frag_bins_orig = Pioneer.getFragBins(nfi)
    fragments_orig = Pioneer.getFragments(nfi)
    n_frag_bins    = length(frag_bins_orig)

    # Step 1: Find prec_mz range across all fragments → compute n_partitions
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

    # Step 2: Distribute fragments preserving full frag bin structure.
    # Walk each frag bin in the original and distribute its fragments across partitions.
    # Each partition gets the same number of frag bins (n_frag_bins) in the same order,
    # but with partition-local first_bin/last_bin ranges.
    partition_fragments = [Pioneer.IndexFragment{T}[] for _ in 1:n_partitions]
    partition_fb_first  = [Vector{UInt32}(undef, n_frag_bins) for _ in 1:n_partitions]
    partition_fb_last   = [Vector{UInt32}(undef, n_frag_bins) for _ in 1:n_partitions]
    starts = Vector{UInt32}(undef, n_partitions)  # reusable per-frag-bin buffer

    for j in 1:n_frag_bins
        fbin = frag_bins_orig[j]
        frag_range = Pioneer.getSubBinRange(fbin)

        # Snapshot current end of each partition's fragments vector
        for k in 1:n_partitions
            starts[k] = UInt32(length(partition_fragments[k]) + 1)
        end

        # Distribute this frag bin's fragments to their respective partitions
        for fi in frag_range
            frag = fragments_orig[fi]
            pmz = Pioneer.getPrecMZ(frag)
            k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1,
                       1, n_partitions)
            push!(partition_fragments[k], frag)
        end

        # Record first/last for each partition's copy of frag bin j
        # When no fragments were added, new_end < starts[k] → empty UnitRange
        for k in 1:n_partitions
            new_end = UInt32(length(partition_fragments[k]))
            partition_fb_first[k][j] = starts[k]
            partition_fb_last[k][j]  = new_end
        end
    end

    # Step 3: Assemble NativeFragmentIndex per partition.
    # Each partition gets the full frag bin array (same lb/ub, partition-local ranges)
    # and a verbatim copy of the original RT bins (indices into frag bins are unchanged).
    partitions = Vector{Pioneer.NativeFragmentIndex{T}}(undef, n_partitions)
    for k in 1:n_partitions
        p_frag_bins = Vector{Pioneer.FragIndexBin{T}}(undef, n_frag_bins)
        for j in 1:n_frag_bins
            orig = frag_bins_orig[j]
            p_frag_bins[j] = Pioneer.FragIndexBin{T}(
                Pioneer.getLow(orig), Pioneer.getHigh(orig),
                partition_fb_first[k][j], partition_fb_last[k][j])
        end
        partitions[k] = Pioneer.NativeFragmentIndex{T}(
            p_frag_bins, copy(rt_bins_orig), partition_fragments[k])
    end

    return PartitionedFragmentIndex{T}(partitions, min_prec_mz, partition_width, n_partitions)
end
