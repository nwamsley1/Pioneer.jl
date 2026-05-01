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

#==========================================================
AbstractPrecursorMap{V} — interface for per-thread keyed-by-precursor-id
storage used during DIA search.

Two on-disk callers carry different lifecycle expectations:

  * `id_to_col` (UInt16) — RESET between every scan. Live keys at any
    moment are the precursors actually matched in the current scan
    (~5k for OlsenAstral). Dense storage wastes most of its
    n_precursors-sized backing array.

  * `precursor_weights` (Float32) — PERSISTENT across scans within a
    file. Live keys grow monotonically over the run and saturate at
    "ever-matched" count. For typical DIA libraries this is well below
    n_precursors, so sparse storage saves memory.

Both can be implemented either with a dense array indexed by precursor_id
(`DensePrecMap`, mirrors classic `ArrayDict`/`Vector` behavior) or with a
hash-based sparse table (`SparsePrecMap`, future Phase 2).

# Required methods on a concrete subtype
  Base.getindex(m, k::Integer)        :: V
  Base.setindex!(m, v::V, k::Integer)
  reset!(m)                            # zero out all live entries
  active_keys(m)                       # iterator yielding (k, v::V)
  n_active(m) :: Int                   # count of live entries
  Base.sizehint!(m, n)                 # backing-specific
==========================================================#

abstract type AbstractPrecursorMap{V} end

"""
    DensePrecMap{V}(n_precursors::Int)

Dense backing: a length-`n_precursors` `Vector{V}` indexed directly by
precursor id, plus an auxiliary `active::Vector{UInt32}` listing the
keys touched since the last `reset!`. Mirrors the original `ArrayDict`
semantics — every read/write is O(1), reset cost is O(n_active).
"""
mutable struct DensePrecMap{V} <: AbstractPrecursorMap{V}
    vals::Vector{V}
    active::Vector{UInt32}        # touched keys since last reset
    n_active::Int
end

function DensePrecMap{V}(n_precursors::Integer) where V
    DensePrecMap{V}(zeros(V, n_precursors), Vector{UInt32}(undef, 0), 0)
end

@inline function Base.getindex(m::DensePrecMap{V}, k::Integer) where V
    @inbounds m.vals[k]
end

@inline function Base.setindex!(m::DensePrecMap{V}, v::V, k::Integer) where V
    @inbounds begin
        # Track the key the first time it transitions away from zero. We
        # only care about reset! efficiency; tracking every set would make
        # `active` grow forever within one scan.
        if iszero(m.vals[k])
            m.n_active += 1
            if m.n_active > length(m.active)
                push!(m.active, UInt32(k))
            else
                m.active[m.n_active] = UInt32(k)
            end
        end
        m.vals[k] = v
    end
    return v
end

function reset!(m::DensePrecMap{V}) where V
    @inbounds for i in 1:m.n_active
        m.vals[m.active[i]] = zero(V)
    end
    m.n_active = 0
    return
end

n_active(m::DensePrecMap) = m.n_active

# Iterator over (key, value) for currently-live entries.
struct DensePrecMapKVIter{V}
    m::DensePrecMap{V}
end

active_keys(m::DensePrecMap) = DensePrecMapKVIter(m)

Base.length(it::DensePrecMapKVIter) = it.m.n_active
Base.eltype(::Type{DensePrecMapKVIter{V}}) where V = Tuple{UInt32, V}

function Base.iterate(it::DensePrecMapKVIter{V}, state::Int = 1) where V
    state > it.m.n_active && return nothing
    @inbounds k = it.m.active[state]
    @inbounds v = it.m.vals[k]
    return ((k, v), state + 1)
end

# Dense map is fixed-size — sizehint! is a no-op (backing is already
# sized to n_precursors at construction).
Base.sizehint!(::DensePrecMap, ::Integer) = nothing

"""
    resize_dense!(m::DensePrecMap{V}, new_n::Integer)

Used by QuadTuning's artificial-id encoding (which needs `3 *
n_precursors + 1` slots). No-op on sparse backings — they grow on demand.
"""
function resize_dense!(m::DensePrecMap{V}, new_n::Integer) where V
    if length(m.vals) < new_n
        old = length(m.vals)
        resize!(m.vals, new_n)
        @inbounds for i in (old + 1):new_n
            m.vals[i] = zero(V)
        end
    end
    return
end

# Catch-all for non-DensePrecMap impls — they don't need explicit resize.
resize_dense!(::AbstractPrecursorMap, ::Integer) = nothing

#==========================================================
SparsePrecMap{V} — Dict-backed alternative.

Memory grows with the live key count rather than n_precursors. Best for:
  * id_to_col (reset every scan → ~5k entries × ~28 B each = ~140 KB
                       vs DensePrecMap's ~38 MB per thread on Olsen)
  * precursor_weights (cumulative matched precursors over a run, typically
                       well below n_precursors for DIA libraries)

Uses Julia stdlib `Dict{UInt32, V}` internally; consumers access only
through the `AbstractPrecursorMap` interface so the dense vs sparse
choice is a one-line swap in `initSimpleSearchContext`.

`Base.getindex(m, k)` returns `zero(V)` when `k` is absent, mirroring
the dense vector's "all unset → zero" semantics. This is essential for
downstream code that expects e.g. `weights[col] = precursor_weights[id]`
to write zero into `weights` for unseen precursors.
==========================================================#

struct SparsePrecMap{V} <: AbstractPrecursorMap{V}
    d::Dict{UInt32, V}
end

function SparsePrecMap{V}(; sizehint::Integer = 0) where V
    d = Dict{UInt32, V}()
    sizehint > 0 && sizehint!(d, sizehint)
    SparsePrecMap{V}(d)
end

@inline function Base.getindex(m::SparsePrecMap{V}, k::Integer) where V
    get(m.d, UInt32(k), zero(V))
end

@inline function Base.setindex!(m::SparsePrecMap{V}, v::V, k::Integer) where V
    m.d[UInt32(k)] = v
    return v
end

reset!(m::SparsePrecMap) = (empty!(m.d); nothing)

n_active(m::SparsePrecMap) = length(m.d)

# Dict iteration already yields (k::UInt32, v::V) pairs.
active_keys(m::SparsePrecMap) = m.d

Base.sizehint!(m::SparsePrecMap, n::Integer) = (sizehint!(m.d, n); m)
