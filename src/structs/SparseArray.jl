mutable struct SparseArray{Ti<:Integer,T<:AbstractFloat}
    n_vals::Int64
    m::Int64
    n::Int64
    rowval::Vector{Ti}
    colval::Vector{Ti}
    nzval::Vector{T}
    mask::Vector{Bool}
    matched::Vector{Bool}
    x::Vector{T}
    colptr::Vector{Ti}
end

function setMask!(sa::SparseArray{Ti, T}, i::Int64, v::Bool) where {Ti<:Integer,T<:AbstractFloat}
    sa.mask[i] = v
end

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
        end

        if lt(o, v.colval[hi], v.colval[lo])
            if lt(o, v.colval[hi], v.colval[mi])
                #v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]

                v.colval[hi], v.colval[lo], v.colval[mi] = v.colval[lo], v.colval[mi], v.colval[hi]
                v.rowval[hi], v.rowval[lo], v.rowval[mi] = v.rowval[lo], v.rowval[mi], v.rowval[hi]
                v.nzval[hi], v.nzval[lo], v.nzval[mi] = v.nzval[lo], v.nzval[mi], v.nzval[hi]
                v.x[hi], v.x[lo], v.x[mi] = v.x[lo], v.x[mi], v.x[hi]
                v.matched[hi], v.matched[lo], v.matched[mi] = v.matched[lo], v.matched[mi], v.matched[hi]

            else
                #v[hi], v[lo] = v[lo], v[hi]

                v.colval[hi], v.colval[lo] = v.colval[lo], v.colval[hi]
                v.rowval[hi], v.rowval[lo] = v.rowval[lo], v.rowval[hi]
                v.nzval[hi], v.nzval[lo] = v.nzval[lo], v.nzval[hi] 
                v.x[hi], v.x[lo] = v.x[hi], v.x[lo] 
                v.matched[hi], v.matched[lo] = v.matched[hi], v.matched[lo] 
            end
        end

        # return the pivot
        return v.colval[lo], v.rowval[lo], v.nzval[lo], v.x[lo], v.matched[lo]
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
    end
    v.colval[j], v.colval[lo] = pivot[1], v.colval[j]
    v.rowval[j], v.rowval[lo] = pivot[2], v.rowval[j]
    v.nzval[j], v.nzval[lo] = pivot[3], v.nzval[j]
    v.x[j], v.x[lo] = pivot[4], v.x[j]
    v.matched[j], v.matched[lo] = pivot[5], v.matched[j]

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
        while j > lo
            #y = v[j-1]
            col_y = v.colval[j - 1]
            row_y = v.rowval[j - 1]
            nzval_y = v.nzval[j - 1]
            x_y = v.x[j - 1]
            matched_y = v.matched[j - 1]
            if !(lt(o, col_x, col_y)::Bool)
                break
            end
            v.colval[j] = col_y
            v.rowval[j] = row_y
            v.nzval[j] = nzval_y
            v.x[j] = x_y
            v.matched[j] = matched_y
            j -= 1
        end
        v.colval[j] = col_x
        v.rowval[j] = row_x
        v.nzval[j] = nzval_x
        v.x[j] = x_x
        v.matched[j] = matched_x
    end
    #scratch
end

function getRowVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.rowval[i]
end

function getColVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.colval[i]
end

function getNzVal(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.nzval[i]
end

function getRowColNZ(sa::SparseArray{Ti, T}, i::Int64) where {Ti<:Integer,T<:AbstractFloat}
    return sa.row_col_nzval_x[i]
end

SparseArray(N::Int) = SparseArray(
                    0,
                    0,
                    0,
                    Vector{Int64}(undef, N), #rowval 
                    Vector{Int64}(undef, N), #colval 
                    Vector{Float32}(undef, N), #nzval
                    ones(Bool, N), #mask
                    ones(Bool, N), #matched,
                    Vector{Float32}(undef, N), #x
                    Vector{Int64}(undef, N), #colptr
                    
)

function reset!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    @turbo for i in range(1, sa.n_vals)
        sa.colval[i] = 0
        sa.rowval[i] = 0
        sa.x[i] = zero(T)
        sa.nzval[i] = zero(T)
        sa.mask[i] = true
        sa.matched[i] = true
        sa.colptr[i] = 0
    end
    sa.n_vals = 0
    sa.m = 0
    sa.n = 0
    return 
end

function sortSparse!(sa::SparseArray{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    #Get sorted indices by column
    #partialsort!(sa.row_col_nzval[1:sa.n_vals], 
    #                    1:sa.n_vals, 
    #                    by = x -> x[2])
    specialsort!(sa, 1, sa.n_vals, Base.Order.Forward);
    max_col = 1
    max_row  = 1
    sa.colptr[1] = 1
    for i in range(1, sa.n_vals - 1)
        #If row greater than max row
        if sa.rowval[i + 1] > max_row
            max_row = sa.rowval[i + 1][1] 
        end 

        if sa.colval[i + 1] == sa.colval[i]
            continue
        else
            max_col += 1
            sa.colptr[max_col] = i + 1
        end
    end
    sa.colptr[max_col + 1] = sa.n_vals
    #Number of columns
    sa.n = max_col
    sa.m = max_row
end

function initResiduals!( r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer,T<:AbstractFloat}
    #resize if necessary
    if length(r) < sa.n_vals
        append!(r, zeros(T, sa.n_vals - length(r)))
    end


    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
            r[sa.rowval[n]] = -sa.x[n]
        end
    end

    for col in range(1, sa.n)
        start = sa.colptr[col]
        stop = sa.colptr[col+1] - 1
        @turbo for n in start:stop
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