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

"""
    LocalFragment

Ultra-compact fragment using partition-local UInt16 precursor IDs.
4 bytes total (UInt16 local_id + UInt8 score + 1 byte padding).
vs IndexFragment{Float32}: 12 bytes.

Requires a per-partition `local_to_global::Vector{UInt32}` to map back
to global precursor IDs after scoring.
"""
struct LocalFragment
    local_id::UInt16
    score::UInt8
end

getPrecID(f::LocalFragment) = f.local_id
getScore(f::LocalFragment) = f.score

"""
    LocalPartition{T}

A single partition's index using LocalFragment with partition-local UInt16 IDs.
Includes `local_to_global` mapping to recover global UInt32 precursor IDs.
"""
struct LocalPartition{T<:AbstractFloat}
    fragment_bins::SoAFragBins{T}
    rt_bins::Vector{FragIndexBin{T}}
    fragments::Vector{LocalFragment}
    local_to_global::Vector{UInt32}  # local_id → global prec_id
    n_local_precs::UInt16
    skip_hints::Vector{UInt16}  # per-frag-bin: bins to +5 Da in getLow
end

getFragBins(lp::LocalPartition) = lp.fragment_bins
getRTBins(lp::LocalPartition) = lp.rt_bins
getFragments(lp::LocalPartition) = lp.fragments
getSkipHints(lp::LocalPartition) = lp.skip_hints

"""
    LocalPartitionedFragmentIndex{T}

Partitioned index using LocalFragment (UInt16 local IDs) for maximum cache
efficiency. Each partition has ≤ 65535 unique precursors.

The search loop uses a small LocalCounter{UInt16, UInt8} (~65K slots) that fits
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

# ── LocalCounter: type-correct Counter for UInt16/UInt8 use ──────────────────

"""
    LocalCounter{I, C}

A type-correct reimplementation of Pioneer's Counter that works with any
integer key/count types. Pioneer's Counter has hardcoded `zero(T)` in
reset!/countFragMatches which doesn't generalize cleanly.

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

# ── SIMD threshold constant ─────────────────────────────────────────────────

const HINT_LINEAR_THRESHOLD = UInt32(128)
