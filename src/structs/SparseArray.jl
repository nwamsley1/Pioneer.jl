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

mutable struct SparseArray{Ti<:Integer,T<:AbstractFloat}
    n_vals::Int64
    m::Int64
    n::Int64
    rowval::Vector{Ti}
    colval::Vector{UInt16}
    nzval::Vector{T}
    matched::Vector{Bool}
    isotope::Vector{UInt8}
    x::Vector{T}
    colptr::Vector{Ti}
    cursors::Vector{Ti}   # pre-allocated, reused per scan for countingSortByCol!
end
SparseArray(N::I) where {I<:Integer} = SparseArray(
                    0,
                    0,
                    0,
                    zeros(I, N), #rowval
                    zeros(UInt16, N), #colval
                    zeros(Float32, N), #nzval
                    ones(Bool, N), #matched,
                    zeros(UInt8, N), #isotope
                    zeros(Float32, N), #x
                    zeros(I, N), #colptr
                    zeros(I, 200), #cursors
)

@inline function selectpivot!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds begin
        mi = Base.midpoint(lo, hi)
        # sort v[mi] <= v[lo] <= v[hi] such that the pivot is immediately in place
        if lt(o, v.colval[lo], v.colval[mi])
            v.colval[mi], v.colval[lo] = v.colval[lo], v.colval[mi]
            v.rowval[mi], v.rowval[lo] = v.rowval[lo], v.rowval[mi]
            v.nzval[mi], v.nzval[lo] = v.nzval[lo], v.nzval[mi]
            v.x[mi], v.x[lo] = v.x[lo], v.x[mi]
            v.matched[mi], v.matched[lo] = v.matched[lo], v.matched[mi]
            v.isotope[mi], v.isotope[lo] = v.isotope[lo], v.isotope[mi]
        end

        if lt(o, v.colval[hi], v.colval[lo])
            if lt(o, v.colval[hi], v.colval[mi])
                #v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]
                v.colval[hi], v.colval[lo], v.colval[mi] = v.colval[lo], v.colval[mi], v.colval[hi]
                v.rowval[hi], v.rowval[lo], v.rowval[mi] = v.rowval[lo], v.rowval[mi], v.rowval[hi]
                v.nzval[hi], v.nzval[lo], v.nzval[mi] = v.nzval[lo], v.nzval[mi], v.nzval[hi]
                v.x[hi], v.x[lo], v.x[mi] = v.x[lo], v.x[mi], v.x[hi]
                v.matched[hi], v.matched[lo], v.matched[mi] = v.matched[lo], v.matched[mi], v.matched[hi]
                v.isotope[hi], v.isotope[lo], v.isotope[mi] = v.isotope[lo], v.isotope[mi], v.isotope[hi]

            else
                #v[hi], v[lo] = v[lo], v[hi]
                v.colval[hi], v.colval[lo] = v.colval[lo], v.colval[hi]
                v.rowval[hi], v.rowval[lo] = v.rowval[lo], v.rowval[hi]
                v.nzval[hi], v.nzval[lo] = v.nzval[lo], v.nzval[hi] 
                v.x[hi], v.x[lo] = v.x[lo], v.x[hi] 
                v.matched[hi], v.matched[lo] = v.matched[lo], v.matched[hi] 
                v.isotope[hi], v.isotope[lo] = v.isotope[lo], v.isotope[hi] 
            end
        end

        # return the pivot
        return v.colval[lo], v.rowval[lo], v.nzval[lo], v.x[lo], v.matched[lo], v.isotope[lo]
    end
end

function partition!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    pivot = selectpivot!(v, lo, hi, o)
    # pivot == v[lo], v[hi] > pivot
    i, j = lo, hi
    @inbounds while true
        i += 1; j -= 1
        while lt(o, v.colval[i], pivot[1]); i += 1; end;
        while lt(o, pivot[1], v.colval[j]); j -= 1; end;
        i >= j && break
        v.colval[i], v.colval[j] = v.colval[j], v.colval[i]
        v.rowval[i], v.rowval[j] = v.rowval[j], v.rowval[i]
        v.nzval[i], v.nzval[j] = v.nzval[j], v.nzval[i]
        v.x[i], v.x[j] = v.x[j], v.x[i]
        v.matched[i], v.matched[j] = v.matched[j], v.matched[i]
        v.isotope[i], v.isotope[j] = v.isotope[j], v.isotope[i]
    end
    v.colval[j], v.colval[lo] = pivot[1], v.colval[j]
    v.rowval[j], v.rowval[lo] = pivot[2], v.rowval[j]
    v.nzval[j], v.nzval[lo] = pivot[3], v.nzval[j]
    v.x[j], v.x[lo] = pivot[4], v.x[j]
    v.matched[j], v.matched[lo] = pivot[5], v.matched[j]
    v.isotope[j], v.isotope[lo] = pivot[6], v.isotope[j]
    # v[j] == pivot
    # v[k] >= pivot for k > j
    # v[i] <= pivot for i < j
    return j
end

function specialsort!(v::SparseArray{Ti, T}, lo::Integer, hi::Integer, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds while lo < hi
        hi-lo <= 20 && return smallsort!(v, lo, hi, o)
        j = partition!(v, lo, hi, o)
        if j-lo < hi-j
            # recurse on the smaller chunk
            # this is necessary to preserve O(log(n))
            # stack space in the worst case (rather than O(n))
            lo < (j-1) && specialsort!(v, lo, j-1, o)
            lo = j+1
        else
            j+1 < hi && specialsort!(v, j+1, hi, o)
            hi = j-1
        end
    end
    return v
end

function smallsort!(v::SparseArray{Ti, T}, lo::Int64, hi::Int64, o::Ordering) where {Ti<:Integer, T<:AbstractFloat}
    #getkw lo hi scratch
    lo_plus_1 = (lo + 1)::Int64
    @inbounds for i = lo_plus_1:hi
        j = i
        col_x = v.colval[i]
        row_x = v.rowval[i]
        nzval_x = v.nzval[i]
        x_x = v.x[i]
        matched_x = v.matched[i]
        isotope_x = v.isotope[i]
        while j > lo
            #y = v[j-1]
            col_y = v.colval[j - 1]
            row_y = v.rowval[j - 1]
            nzval_y = v.nzval[j - 1]
            x_y = v.x[j - 1]
            matched_y = v.matched[j - 1]
            isotope_y = v.isotope[j - 1]
            if !(lt(o, col_x, col_y)::Bool)
                break
            end
            v.colval[j] = col_y
            v.rowval[j] = row_y
            v.nzval[j] = nzval_y
            v.x[j] = x_y
            v.matched[j] = matched_y
            v.isotope[j] = isotope_y
            j -= 1
        end
        v.colval[j] = col_x
        v.rowval[j] = row_x
        v.nzval[j] = nzval_x
        v.x[j] = x_x
        v.matched[j] = matched_x
        v.isotope[j] = isotope_x
    end
    #scratch
end



function reset!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    @turbo for i in range(1, sa.n_vals)
    #for i in range(1, sa.n_vals)
        sa.colval[i] = zero(UInt16)
        sa.rowval[i] = zero(Ti)
        sa.x[i] = zero(T)
        sa.nzval[i] = zero(T)
        sa.matched[i] = true
        sa.colptr[i] = zero(Ti)
    end
    @turbo for i in range(1, sa.n_vals)
        sa.isotope[i] = zero(eltype(sa.isotope))
    end
    sa.n_vals = 0
    sa.m = 0
    sa.n = 0
    return 
end



"""
    countingSortByCol!(sa::SparseArray)

O(n) counting sort by column index, replacing the O(n log n) quicksort.
Uses slack space after n_vals in existing arrays as temp storage — zero allocation.
Builds colptr as a natural byproduct of the sort.
"""
function countingSortByCol!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    n = sa.n_vals
    n == 0 && return

    # Find max column and max row
    max_col = zero(UInt16)
    max_row = zero(Ti)
    @inbounds for i in 1:n
        if sa.colval[i] > max_col; max_col = sa.colval[i]; end
        if sa.rowval[i] > max_row; max_row = sa.rowval[i]; end
    end

    # Ensure cursors buffer is large enough
    if Int(max_col) > length(sa.cursors)
        resize!(sa.cursors, Int(max_col) + 64)
    end

    # Ensure arrays have capacity >= 2n for temp space
    capacity = length(sa.colval)
    if capacity < 2 * n
        grow = 2 * n - capacity
        append!(sa.colval, zeros(eltype(sa.colval), grow))
        append!(sa.rowval, zeros(eltype(sa.rowval), grow))
        append!(sa.nzval, zeros(eltype(sa.nzval), grow))
        append!(sa.x, zeros(eltype(sa.x), grow))
        append!(sa.matched, zeros(eltype(sa.matched), grow))
        append!(sa.isotope, zeros(eltype(sa.isotope), grow))
    end

    # Count per column → build colptr via prefix sum
    mc = Int(max_col)
    @inbounds for c in 1:mc+1; sa.colptr[c] = zero(Ti); end
    @inbounds for i in 1:n; sa.colptr[Int(sa.colval[i])] += one(Ti); end
    total = one(Ti)
    @inbounds for c in 1:mc
        count = sa.colptr[c]
        sa.colptr[c] = total
        total += count
    end
    sa.colptr[mc + 1] = total

    # Copy to temp space (slack after n_vals — no allocation)
    @inbounds for i in 1:n
        sa.colval[n + i] = sa.colval[i]
        sa.rowval[n + i] = sa.rowval[i]
        sa.nzval[n + i]  = sa.nzval[i]
        sa.x[n + i]      = sa.x[i]
        sa.matched[n + i] = sa.matched[i]
        sa.isotope[n + i] = sa.isotope[i]
    end

    # Scatter back in column order
    @inbounds for c in 1:mc; sa.cursors[c] = sa.colptr[c]; end
    @inbounds for i in 1:n
        c = Int(sa.colval[n + i])
        dest = Int(sa.cursors[c])
        sa.colval[dest] = sa.colval[n + i]
        sa.rowval[dest] = sa.rowval[n + i]
        sa.nzval[dest]  = sa.nzval[n + i]
        sa.x[dest]      = sa.x[n + i]
        sa.matched[dest] = sa.matched[n + i]
        sa.isotope[dest] = sa.isotope[n + i]
        sa.cursors[c] += one(Ti)
    end

    sa.n = Int(max_col)
    sa.m = Int(max_row)
end

function sortSparse!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    countingSortByCol!(sa)
end

function initResiduals!( r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    #resize if necessary
    if length(r) < sa.m
        append!(r, zeros(T, sa.m - length(r)))
    end

    @turbo for i in range(1, sa.m)
        r[i] = zero(T)
    end

    for n in range(1, sa.n_vals)
        if iszero(r[sa.rowval[n]])
            r[sa.rowval[n]] = -sa.x[n]
        end
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        #@turbo for n in start:stop
        for n in start:stop
            r[sa.rowval[n]] += w[col]*sa.nzval[n]
        end
    end
end

function multiply!(sa::SparseArray{Ti,T}, w::Vector{T}, r::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    #resize if necessary
    if length(r) < length(sa.n_vals)
        append!(r, zeros(T, length(sa.n_vals) - length(r)))
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
            r[sa.rowval[n]] += w[col]*sa.nzval[n]
        end
    end
end