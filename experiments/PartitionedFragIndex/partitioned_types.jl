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

# ── Compact fragment type for partitioned index ──────────────────────────────

"""
    CompactFragment

Minimal fragment for the partitioned index. Since the partition guarantees
precursor m/z proximity, we only need the precursor ID (for scoring) and
the score (rank-based weight). No prec_mz or charge needed.

Layout: 5 bytes payload, padded to 8 bytes by Julia.
vs IndexFragment{Float32}: 10 bytes payload, padded to 12 bytes.
→ 33% memory reduction on the fragment array.
"""
struct CompactFragment
    prec_id::UInt32
    score::UInt8
end

# Make CompactFragment quack like IndexFragment for the scoring loop
Pioneer.getPrecID(f::CompactFragment) = f.prec_id
Pioneer.getScore(f::CompactFragment) = f.score

"""
    CompactPartition{T}

A single partition's index using CompactFragment instead of IndexFragment.
Same frag_bins and rt_bins as NativeFragmentIndex, but smaller fragments.
"""
struct CompactPartition{T<:AbstractFloat}
    fragment_bins::Vector{Pioneer.FragIndexBin{T}}
    rt_bins::Vector{Pioneer.FragIndexBin{T}}
    fragments::Vector{CompactFragment}
end

# Accessors for CompactPartition
getFragBins(cp::CompactPartition) = cp.fragment_bins
getRTBins(cp::CompactPartition) = cp.rt_bins
getFragments(cp::CompactPartition) = cp.fragments

# Forward to Pioneer accessors for NativeFragmentIndex (so searchScanPartitioned! works
# with both PartitionedFragmentIndex and CompactPartitionedFragmentIndex)
getFragBins(nfi::Pioneer.NativeFragmentIndex) = Pioneer.getFragBins(nfi)
getRTBins(nfi::Pioneer.NativeFragmentIndex) = Pioneer.getRTBins(nfi)
getFragments(nfi::Pioneer.NativeFragmentIndex) = Pioneer.getFragments(nfi)

"""
    CompactPartitionedFragmentIndex{T}

Like PartitionedFragmentIndex but uses CompactFragment for ~33% less memory.
"""
struct CompactPartitionedFragmentIndex{T<:AbstractFloat}
    partitions::Vector{CompactPartition{T}}
    prec_mz_min::T
    partition_width::T
    n_partitions::Int
end

getPartitions(pfi::CompactPartitionedFragmentIndex) = pfi.partitions
getPartition(pfi::CompactPartitionedFragmentIndex, k::Int) = pfi.partitions[k]
getNPartitions(pfi::CompactPartitionedFragmentIndex) = pfi.n_partitions

@inline function get_partition_idx(pfi::CompactPartitionedFragmentIndex{T}, prec_mz::T) where {T}
    return clamp(floor(Int, (prec_mz - pfi.prec_mz_min) / pfi.partition_width) + 1, 1, pfi.n_partitions)
end
