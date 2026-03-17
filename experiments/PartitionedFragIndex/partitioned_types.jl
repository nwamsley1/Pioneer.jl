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

# ── Local-ID fragment type for partitioned index ─────────────────────────────

const MAX_LOCAL_PRECS = 65535  # UInt16 max

"""
    LocalFragment

Ultra-compact fragment using partition-local UInt16 precursor IDs.
4 bytes total (UInt16 local_id + UInt8 score + 1 byte padding).
vs CompactFragment: 8 bytes. vs IndexFragment{Float32}: 12 bytes.

Requires a per-partition `local_to_global::Vector{UInt32}` to map back
to global precursor IDs after scoring.
"""
struct LocalFragment
    local_id::UInt16
    score::UInt8
end

Pioneer.getPrecID(f::LocalFragment) = f.local_id
Pioneer.getScore(f::LocalFragment) = f.score

"""
    LocalPartition{T}

A single partition's index using LocalFragment with partition-local UInt16 IDs.
Includes `local_to_global` mapping to recover global UInt32 precursor IDs.
"""
struct LocalPartition{T<:AbstractFloat}
    fragment_bins::Vector{Pioneer.FragIndexBin{T}}
    rt_bins::Vector{Pioneer.FragIndexBin{T}}
    fragments::Vector{LocalFragment}
    local_to_global::Vector{UInt32}  # local_id → global prec_id
    n_local_precs::UInt16
end

getFragBins(lp::LocalPartition) = lp.fragment_bins
getRTBins(lp::LocalPartition) = lp.rt_bins
getFragments(lp::LocalPartition) = lp.fragments

"""
    LocalPartitionedFragmentIndex{T}

Partitioned index using LocalFragment (UInt16 local IDs) for maximum cache
efficiency. Each partition has ≤ 65535 unique precursors.

The search loop uses a small Counter{UInt16, UInt8} (~65K slots) that fits
in L1/L2 cache, then translates results back to global UInt32 IDs.
"""
struct LocalPartitionedFragmentIndex{T<:AbstractFloat}
    partitions::Vector{LocalPartition{T}}
    partition_bounds::Vector{Tuple{T, T}}  # (prec_mz_min, prec_mz_max) per partition
    n_partitions::Int
end

getPartitions(pfi::LocalPartitionedFragmentIndex) = pfi.partitions
getPartition(pfi::LocalPartitionedFragmentIndex, k::Int) = pfi.partitions[k]
getNPartitions(pfi::LocalPartitionedFragmentIndex) = pfi.n_partitions

"""
Find the range of partitions whose prec_mz bounds overlap [query_min, query_max].
Partition bounds are sorted by prec_mz_min, so we use binary search.
"""
# ── LocalCounter: type-correct Counter for non-Float32 use ───────────────────

"""
    LocalCounter{I, C}

A type-correct reimplementation of Pioneer's Counter that works with any
integer key/count types. Pioneer's Counter has hardcoded `zero(T)` and
`zero(Float32)` in reset!/countFragMatches which breaks for non-UInt32/UInt8.

Same branchless `inc!` algorithm as Pioneer's Counter.
"""
mutable struct LocalCounter{I, C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
end

function LocalCounter(::Type{I}, ::Type{C}, n::Int) where {I, C<:Unsigned}
    LocalCounter{I, C}(zeros(I, n), zeros(C, n), 1)
end

@inline function inc!(c::LocalCounter{I, C}, id::I, score::C) where {I, C}
    @inbounds @fastmath begin
        no_previous_encounter = c.counts[id] === zero(C)
        c.ids[c.size] = id
        c.size += no_previous_encounter
        c.counts[id] += score
    end
    return nothing
end

@inline function reset!(c::LocalCounter{I, C}) where {I, C}
    @inbounds for i in 1:(c.size - 1)
        c.counts[c.ids[i]] = zero(C)
        c.ids[i] = zero(I)
    end
    c.size = 1
    return nothing
end

@inline function get_partition_range(pfi::LocalPartitionedFragmentIndex{T}, query_min::T, query_max::T) where {T}
    bounds = pfi.partition_bounds
    n = pfi.n_partitions

    # Find first partition whose prec_mz_max >= query_min
    lo, hi = 1, n
    first_k = n + 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if bounds[mid][2] >= query_min
            first_k = mid
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    # Find last partition whose prec_mz_min <= query_max
    lo, hi = first_k, n
    last_k = first_k - 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if bounds[mid][1] <= query_max
            last_k = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    return first_k, last_k
end
