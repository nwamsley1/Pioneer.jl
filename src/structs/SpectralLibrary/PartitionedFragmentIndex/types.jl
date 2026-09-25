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

# ── Fragment index bin types (shared by partitioned index and Arrow I/O) ──────

abstract type FragmentIndexBin{T<:AbstractFloat} end

getLow(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.lb
getHigh(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.ub
getSubBinRange(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.first_bin:fb.last_bin

"""
    FragIndexBin{T}

A bin in the fragment index with m/z bounds and a range of sub-bins (fragment IDs).
Used as RT bins and fragment m/z bins in the partitioned index.
"""
struct FragIndexBin{T<:AbstractFloat} <: FragmentIndexBin{T}
    lb::T
    ub::T
    first_bin::UInt32
    last_bin::UInt32
end
ArrowTypes.arrowname(::Type{FragIndexBin{Float32}}) = :FragIndexBin
ArrowTypes.JuliaType(::Val{:FragIndexBin}) = FragIndexBin

"""
    IndexFragment{T}

A fragment in the non-partitioned index with global precursor ID.
Used during index building and by the Arrow serialization layer.
"""
struct IndexFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32
    prec_mz::T
    score::UInt8
    charge::UInt8
end
ArrowTypes.arrowname(::Type{IndexFragment{Float32}}) = :IndexFragment
ArrowTypes.JuliaType(::Val{:IndexFragment}) = IndexFragment

getScore(ind_frag::IndexFragment{T}) where {T<:AbstractFloat} = ind_frag.score

# ── SoA layout for fragment bins (cache-friendly for field-specific scans) ────

"""
    SoAFragBins{T}

Struct-of-arrays layout for fragment bin data. Each field-specific scan
(e.g., scanning `highs` to find first bin ≥ threshold) touches only one
contiguous array instead of striding through 16-byte AoS records.

Enables SIMD "find first ≥ threshold" over the `highs` array.
"""
struct SoAFragBins{T<:AbstractFloat}
    lows::Vector{T}
    highs::Vector{T}
    first_bins::Vector{UInt32}
    last_bins::Vector{UInt32}
end

@inline Base.length(s::SoAFragBins) = length(s.lows)
@inline Base.isempty(s::SoAFragBins) = isempty(s.lows)

# ── Local-ID fragment type for partitioned index ─────────────────────────────

const MAX_LOCAL_PRECS = 65535  # UInt16 max

# Two index variants that differ only in the width of the partition-local precursor ID:
#   UInt16: LocalFragment / LocalPartition / LocalPartitionedFragmentIndex (the original types, unchanged so that
#           every existing serialized library still loads; partitions split at 65,535 precursors)
#   UInt32: LocalFragment32 / LocalPartition32 / LocalPartitionedFragmentIndex32 (no practical split limit)
# The search is written against the abstract supertypes and the Counter's ID type, so each variant compiles to its
# own specialized code; the UInt16 path is the same machine code as before the UInt32 variant existed.
abstract type AbstractLocalFragment end
abstract type AbstractLocalPartition{T<:AbstractFloat} end
abstract type AbstractLocalPartitionedFragmentIndex{T<:AbstractFloat} end

"""
    LocalFragment

Ultra-compact fragment using partition-local UInt16 precursor IDs.
4 bytes total (UInt16 local_id + UInt8 score + 1 byte padding).
vs IndexFragment{Float32}: 12 bytes.

Requires a per-partition `local_to_global::Vector{UInt32}` to map back
to global precursor IDs after scoring.
"""
struct LocalFragment <: AbstractLocalFragment
    local_id::UInt16
    score::UInt8
end

getPrecID(f::LocalFragment) = f.local_id
getScore(f::LocalFragment) = f.score

"""
    LocalFragment32

`LocalFragment` with a UInt32 partition-local precursor ID (8 bytes with padding), for partitions of more than
65,535 precursors.
"""
struct LocalFragment32 <: AbstractLocalFragment
    local_id::UInt32
    score::UInt8
end

getPrecID(f::LocalFragment32) = f.local_id
getScore(f::LocalFragment32) = f.score

"Partition-local precursor ID type of a fragment type."
local_id_type(::Type{LocalFragment}) = UInt16
local_id_type(::Type{LocalFragment32}) = UInt32
"Largest number of precursors a partition may hold with local IDs of type `I` (Counter slot 0 is unused)."
max_local_precs(::Type{UInt16}) = MAX_LOCAL_PRECS
max_local_precs(::Type{UInt32}) = Int(typemax(UInt32)) - 1
"Fragment type for local IDs of type `I`."
local_fragment_type(::Type{UInt16}) = LocalFragment
local_fragment_type(::Type{UInt32}) = LocalFragment32

"""
    LocalPartition{T}

A single partition's index using LocalFragment with partition-local UInt16 IDs.
Includes `local_to_global` mapping to recover global UInt32 precursor IDs.
"""
struct LocalPartition{T<:AbstractFloat} <: AbstractLocalPartition{T}
    fragment_bins::SoAFragBins{T}
    rt_bins::Vector{FragIndexBin{T}}
    fragments::Vector{LocalFragment}
    local_to_global::Vector{UInt32}  # local_id → global prec_id
    n_local_precs::UInt16
    skip_hints::Vector{UInt16}  # per-frag-bin: bins to +5 Da in getLow
end

"""
    LocalPartition32{T}

`LocalPartition` with `LocalFragment32` fragments and a UInt32 precursor count.
"""
struct LocalPartition32{T<:AbstractFloat} <: AbstractLocalPartition{T}
    fragment_bins::SoAFragBins{T}
    rt_bins::Vector{FragIndexBin{T}}
    fragments::Vector{LocalFragment32}
    local_to_global::Vector{UInt32}  # local_id → global prec_id
    n_local_precs::UInt32
    skip_hints::Vector{UInt16}  # per-frag-bin: bins to +5 Da in getLow
end

getFragBins(lp::AbstractLocalPartition) = lp.fragment_bins
getRTBins(lp::AbstractLocalPartition) = lp.rt_bins
getFragments(lp::AbstractLocalPartition) = lp.fragments
getSkipHints(lp::AbstractLocalPartition) = lp.skip_hints

"Partition type for local IDs of type `I`."
local_partition_type(::Type{UInt16}) = LocalPartition
local_partition_type(::Type{UInt32}) = LocalPartition32

"""
    LocalPartitionedFragmentIndex{T}

Partitioned index using LocalFragment (UInt16 local IDs) for maximum cache
efficiency. Each partition has ≤ 65535 unique precursors.

The search loop uses a small LocalCounter{UInt16, UInt8} (~65K slots) that fits
in L1/L2 cache, then translates results back to global UInt32 IDs.
"""
struct LocalPartitionedFragmentIndex{T<:AbstractFloat} <: AbstractLocalPartitionedFragmentIndex{T}
    partitions::Vector{LocalPartition{T}}
    partition_bounds::Vector{Tuple{T, T}}  # (prec_mz_min, prec_mz_max) per partition
    n_partitions::Int
end

"""
    LocalPartitionedFragmentIndex32{T}

`LocalPartitionedFragmentIndex` of `LocalPartition32`s (UInt32 local IDs): partitions are not split at 65,535
precursors, so they keep the nominal `partition_width`. The search counter grows to the largest partition.
"""
struct LocalPartitionedFragmentIndex32{T<:AbstractFloat} <: AbstractLocalPartitionedFragmentIndex{T}
    partitions::Vector{LocalPartition32{T}}
    partition_bounds::Vector{Tuple{T, T}}  # (prec_mz_min, prec_mz_max) per partition
    n_partitions::Int
end

getPartitions(pfi::AbstractLocalPartitionedFragmentIndex) = pfi.partitions
getPartition(pfi::AbstractLocalPartitionedFragmentIndex, k::Int) = pfi.partitions[k]
getNPartitions(pfi::AbstractLocalPartitionedFragmentIndex) = pfi.n_partitions

"Partition-local precursor ID type of an index (UInt16 or UInt32); a compile-time constant for a concrete index type."
local_id_type(::LocalPartitionedFragmentIndex) = UInt16
local_id_type(::LocalPartitionedFragmentIndex32) = UInt32
"Index type for local IDs of type `I`."
local_index_type(::Type{UInt16}) = LocalPartitionedFragmentIndex
local_index_type(::Type{UInt32}) = LocalPartitionedFragmentIndex32

"""
Find the range of partitions whose prec_mz bounds overlap [query_min, query_max].
Partition bounds are sorted by prec_mz_min, so we use binary search.
"""
@inline function get_partition_range(pfi::AbstractLocalPartitionedFragmentIndex{T}, query_min::T, query_max::T) where {T}
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

# ── SIMD threshold constant ─────────────────────────────────────────────────

const HINT_LINEAR_THRESHOLD = UInt32(128)
