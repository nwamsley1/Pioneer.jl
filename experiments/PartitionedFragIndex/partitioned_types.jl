"""
    PartitionedFragmentIndex{T<:AbstractFloat}

A precursor m/z-partitioned fragment index. Each partition contains only fragments
from precursors within a narrow m/z range (partition_width Da). Within each partition,
ALL fragments unconditionally match the precursor m/z criterion (up to partition_width-wide
false positives), eliminating the binary search over prec_mz entirely.

Each partition reuses `NativeFragmentIndex{T}` and shares the same RT bin boundaries
as the original index.
"""
struct PartitionedFragmentIndex{T<:AbstractFloat}
    partitions::Vector{Pioneer.NativeFragmentIndex{T}}
    prec_mz_min::T
    partition_width::T
    n_partitions::Int
end

getPartitions(pfi::PartitionedFragmentIndex) = pfi.partitions
getPartition(pfi::PartitionedFragmentIndex, k::Int) = pfi.partitions[k]
getNPartitions(pfi::PartitionedFragmentIndex) = pfi.n_partitions

"""
    get_partition_idx(pfi, prec_mz) -> Int

Compute which partition a precursor m/z falls into.
"""
@inline function get_partition_idx(pfi::PartitionedFragmentIndex{T}, prec_mz::T) where {T}
    return clamp(floor(Int, (prec_mz - pfi.prec_mz_min) / pfi.partition_width) + 1, 1, pfi.n_partitions)
end
