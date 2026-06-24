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

"""
    Counter{I, C}

Array-backed sparse counter for accumulating scores by integer key.
Used in fragment index search to count matching fragments per precursor.

The branchless `inc!` algorithm avoids branch misprediction in the hot loop:
instead of `if counts[id]==0 then track id`, it unconditionally writes the id
and increments size by the boolean `counts[id]===zero(C)`.

Fields:
- `ids`: insertion-order list of keys that have been incremented
- `counts`: array indexed by key, stores accumulated scores
- `size`: next insertion slot (1-indexed, so `size-1` = number of unique keys)
"""
mutable struct Counter{I,C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
    function Counter(I::DataType, C::DataType, size::Int)
        new{I, C}(zeros(I, size), zeros(C, size), 1)
    end
end

# ── Bitmask filter types (used by Counter and fragment index search) ─────
#
# Two strategies for deciding which bitmask patterns pass:
#   CountFilter: count_ones(score) >= min_score (training / calibration)
#   LUTFilter:   learned lookup table (main search after calibration)

abstract type AbstractBitVecFilter end

struct CountFilter <: AbstractBitVecFilter
    min_score::UInt8
end

struct LUTFilter <: AbstractBitVecFilter
    bits::NTuple{4, UInt64}  # 256-bit bitset (32 bytes, fits one cache line)
    function LUTFilter(v::Vector{Bool})
        length(v) == 256 || error("LUTFilter requires exactly 256 entries")
        b = ntuple(4) do w
            word = UInt64(0)
            for bit in 0:63
                idx = (w - 1) * 64 + bit + 1
                v[idx] && (word |= UInt64(1) << bit)
            end
            word
        end
        new(b)
    end
end

@inline passes_filter(f::CountFilter, score::UInt8) = count_ones(score) >= f.min_score
@inline function passes_filter(f::LUTFilter, score::UInt8)
    word = (score >> 6) + 1        # which UInt64 (1-4)
    bit = score & UInt8(0x3f)      # bit position within word
    @inbounds return (f.bits[word] >> bit) & UInt64(1) != UInt64(0)
end

"""
    PatternAccumulator

Accumulates bitmask pattern counts (target/decoy × 256 bins) directly in the
fragment index scoring loop. One pair of vectors per thread to avoid contention.
Created once, reused across files — no per-scan allocation.
"""
struct PatternAccumulator{D<:AbstractVector{Bool}}
    thread_tc::Vector{Vector{Int}}  # thread_tc[tid][score+1]
    thread_dc::Vector{Vector{Int}}
    is_decoy::D                     # indexed by global precursor_idx (concrete: avoids boxing in accumulate!)
    min_score::UInt8

    function PatternAccumulator(n_threads::Int, is_decoy::AbstractVector{Bool}; min_score::UInt8 = UInt8(2))
        new{typeof(is_decoy)}(
            [zeros(Int, 256) for _ in 1:n_threads],
            [zeros(Int, 256) for _ in 1:n_threads],
            is_decoy, min_score)
    end
end

@inline function accumulate!(acc::PatternAccumulator, tid::Int, global_pid::UInt32, score::UInt8)
    idx = Int(score) + 1
    @inbounds if acc.is_decoy[global_pid]
        acc.thread_dc[tid][idx] += 1
    else
        acc.thread_tc[tid][idx] += 1
    end
    return nothing
end

function merge_accumulator(acc::PatternAccumulator)
    tc = zeros(Int, 256); dc = zeros(Int, 256)
    for tid in eachindex(acc.thread_tc)
        tc .+= acc.thread_tc[tid]; dc .+= acc.thread_dc[tid]
    end
    return tc, dc
end

"""
    FragIndexScratch

Reusable per-thread scratch buffers for `searchFragmentIndexPartitionMajorHinted`.
A single instance lives on `SearchContext` and is reused across all calls (across
files, across BitVec batches), so the per-call ~50–100 MB transient allocation
stops happening once the buffers have grown to fit. Fields are sized lazily by
`prepare!` — the constructor allocates empty placeholders.

Fields (all length `n_threads`):
- `si`, `pid`        : per-thread output buffers (Int32 / UInt32, sized to est_per_thread)
- `counters`         : per-thread `Counter{UInt16,UInt8}` (sized to `max_local+1`)
- `int_bufs`         : per-thread intensity buffers for top-N peak filtering
"""
mutable struct FragIndexScratch
    si::Vector{Vector{Int32}}
    pid::Vector{Vector{UInt32}}
    counters::Vector{Counter{UInt16,UInt8}}
    int_bufs::Vector{Vector{Float32}}
    function FragIndexScratch(n_threads::Int)
        new(
            [Int32[]    for _ in 1:n_threads],
            [UInt32[]   for _ in 1:n_threads],
            [Counter(UInt16, UInt8, 1) for _ in 1:n_threads],
            [Float32[]  for _ in 1:n_threads],
        )
    end
end

"""
    prepare!(s::FragIndexScratch; n_threads, est_per_thread, counter_size, int_buf_size)

Grow each per-thread buffer to at least the requested size. Counters are
replaced (not resized) when `counter_size` exceeds their current capacity —
the freshly-allocated Counter is already zeroed.
"""
function prepare!(s::FragIndexScratch;
        n_threads::Int,
        est_per_thread::Int,
        counter_size::Int,
        int_buf_size::Int)
    @assert length(s.si) == n_threads
    @inbounds for tid in 1:n_threads
        length(s.si[tid])  < est_per_thread && resize!(s.si[tid],  est_per_thread)
        length(s.pid[tid]) < est_per_thread && resize!(s.pid[tid], est_per_thread)
        if length(s.counters[tid].counts) < counter_size
            s.counters[tid] = Counter(UInt16, UInt8, counter_size)
        end
        length(s.int_bufs[tid]) < int_buf_size && resize!(s.int_bufs[tid], int_buf_size)
    end
    return s
end

"""
    inc!(c::Counter, id, score)

Accumulate `score` for `id`. Branchless: tracks new ids without branching.
"""
function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned}
    @inbounds @fastmath begin
        no_previous_encounter = c.counts[id]===zero(C)
        c.ids[c.size] = id
        c.size += no_previous_encounter
        c.counts[id] += pred_intensity;
    end
    return nothing
end

"""
    or!(c::Counter, id, mask)

Bitwise-OR `mask` into `id`. Branchless: tracks new ids without branching.
Used for bitmask fragment index scoring where each bit represents a fragment rank.
"""
function or!(c::Counter{I,C}, id::I, mask::C) where {I,C<:Unsigned}
    @inbounds @fastmath begin
        no_previous_encounter = c.counts[id] === zero(C)
        c.ids[c.size] = id
        c.size += no_previous_encounter
        c.counts[id] |= mask
    end
    return nothing
end

"""
    filter_counter!(c::Counter, f::AbstractBitVecFilter) -> Int

Compact passing IDs to the front of `c.ids[1:n_pass]` in-place.
Returns the number of passing IDs. After this call, `c.ids[1:n_pass]`
contains only IDs whose bitmask passes the filter, and `c.counts` still
holds the original scores for lookup.

Does NOT reset the counter — call `reset!` after reading the results.
"""
function filter_counter!(c::Counter{I,C}, f::AbstractBitVecFilter) where {I,C<:Unsigned}
    n_pass = 0
    @inbounds for i in 1:(c.size - 1)
        id = c.ids[i]
        if passes_filter(f, c.counts[id])
            n_pass += 1
            c.ids[n_pass] = id
        end
    end
    return n_pass
end

"""
    reset!(c::Counter)

Zero out all tracked ids and counts, reset size to 1.
"""
function reset!(c::Counter{I,C}) where {I,C<:Unsigned}
    @inbounds for i in 1:(c.size - 1)
        c.counts[c.ids[i]] = zero(C)
        c.ids[i] = zero(I)
    end
    c.size = 1
    return nothing
end
