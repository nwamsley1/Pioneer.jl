# Minimal self-contained SparseArray for testing.
# Mirrors the fields from Pioneer.jl's src/structs/SparseArray.jl
# but strips sorting/reset helpers that depend on @turbo and Base.Order.

mutable struct SparseArray{Ti<:Integer,T<:AbstractFloat}
    n_vals::Int64
    m::Int64        # number of rows
    n::Int64        # number of columns
    rowval::Vector{Ti}
    colval::Vector{UInt16}
    nzval::Vector{T}
    matched::Vector{Bool}
    isotope::Vector{UInt8}
    x::Vector{T}           # observed values (one per nonzero entry)
    colptr::Vector{Ti}
end

"""
    build_SparseArray(A, y)

Build a SparseArray from dense matrix `A` (m×n, Float32) and observed vector `y` (length m, Float32).
Only nonzero entries of `A` are stored.  `sa.x` stores the observed value `y[row]` for each nonzero.
"""
function build_SparseArray(A::Matrix{Float32}, y::Vector{Float32})
    m, n = size(A)
    @assert length(y) == m

    # Collect nonzeros column-major
    rows   = Int64[]
    cols   = UInt16[]
    vals   = Float32[]
    obs    = Float32[]
    colptr = Int64[1]

    for j in 1:n
        for i in 1:m
            if A[i,j] != 0f0
                push!(rows, i)
                push!(cols, UInt16(j))
                push!(vals, A[i,j])
                push!(obs,  y[i])
            end
        end
        push!(colptr, length(rows) + 1)
    end

    n_vals = length(rows)
    sa = SparseArray{Int64, Float32}(
        n_vals, m, n,
        rows, cols, vals,
        ones(Bool, n_vals),
        zeros(UInt8, n_vals),
        obs,
        colptr
    )
    return sa
end
