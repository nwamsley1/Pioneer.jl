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
    AbstractSparseDesignMatrix{Ti, T}

Supertype for column-sparse design matrices used by the deconvolution
solvers and distance-metric code. Currently only one concrete subtype
ships: `SparseArrayFused` (CSC order written directly by `run_fused!`).

Solver and metric helpers (`initResiduals!`, `initMu!`, `initObserved!`,
`updateMu!`, `solveOLS!`, `solvePoissonMM_fast!`, `getDistanceMetrics`)
dispatch on this supertype so additional layouts can be added without
re-typing every helper.
"""
abstract type AbstractSparseDesignMatrix{Ti<:Integer, T<:AbstractFloat} end

"""
    matched_at(H, i) -> Bool
    isotope_at(H, i) -> UInt8
    rank_at(H, i) -> UInt8

Per-entry accessors that abstract over storage layout. `SparseArrayFused`
packs `matched + isotope` into a single `UInt8`, so the methods on it
unpack the byte. Methods here are generic forwarders for any future
layout that stores them as separate arrays. Layouts without rank storage
return 0, which is outside the top-8 fragment-rank summaries.
"""
@inline matched_at(H::AbstractSparseDesignMatrix, i::Integer) = H.matched[i]
@inline isotope_at(H::AbstractSparseDesignMatrix, i::Integer) = H.isotope[i]
@inline rank_at(::AbstractSparseDesignMatrix, ::Integer) = UInt8(0)

"""
    initResiduals!(r, sa, w)

Initialize residual vector r = Hw - x, where H is the design matrix (sa),
w is the weight vector, and x is the observed intensities. Used at the
start of deconvolution iterations.
"""
function initResiduals!(r::Vector{T}, sa::AbstractSparseDesignMatrix{Ti,T}, w::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    if length(r) < sa.m
        append!(r, zeros(T, sa.m - length(r)))
    end

    fill!(@view(r[1:sa.m]), zero(T))

    for n in range(1, sa.n_vals)
        if iszero(r[sa.rowval[n]])
            r[sa.rowval[n]] = -sa.x[n]
        end
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        for n in start:stop
            r[sa.rowval[n]] += w[col]*sa.nzval[n]
        end
    end
end
