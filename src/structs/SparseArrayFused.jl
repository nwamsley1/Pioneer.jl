# Fused-path design matrix.
#
# Written in CSC order directly by `run_fused!` — no post-hoc sortSparse!
# needed. Drops `colval` (redundant once CSC is built). Keeps the same
# `rowval/nzval/x/matched/isotope/colptr` fields as `SparseArray` so all
# downstream solver + metric code works unchanged via the shared
# `AbstractSparseDesignMatrix` supertype.

"""
    SparseArrayFused{Ti, T} <: AbstractSparseDesignMatrix{Ti, T}

Column-sparse design matrix produced by the fused per-precursor scan path.

# Fields
- `n_vals::Int64`: number of non-zero entries currently stored.
- `m::Int64`: number of rows (unique spectrum peaks + unique miss rows).
- `n::Int64`: number of columns (precursors that produced ≥1 match).
- `rowval::Vector{Ti}`: row index per entry. Real peak indices for matched
  entries; synthetic rows > n_peaks for missed fragments (one per miss,
  never collide with peak indices or with each other).
- `nzval::Vector{T}`: predicted fragment intensity.
- `x::Vector{T}`: observed intensity (zero for misses).
- `meta::Vector{UInt8}`: packed match-flag + isotope index. Bit 0 = matched
  (1 = matched, 0 = miss). Bits 1-7 = isotope index (0=mono, 1=M+1, …, max
  127). Replaces classic SparseArray's separate `matched::Vector{Bool}` and
  `isotope::Vector{UInt8}` arrays — one byte per entry rather than two.
  Read via the `matched_at(H, i)` and `isotope_at(H, i)` accessors defined
  on `AbstractSparseDesignMatrix`.
- `rank::Vector{UInt8}`: fragment rank for top-rank shadow/fitted summaries.
- `colptr::Vector{Ti}`: CSC column pointer. `colptr[c]:colptr[c+1]-1` gives
  the entry range for column `c`. Written directly by `finalize_column!`;
  no sortSparse! step.

# Lifecycle
1. `reset!(sa)` — clear counters + zero touched entries.
2. `run_fused!` — fills all fields in CSC order while matching/scoring.
3. Solver (`solve_deconvolution!`) + metrics (`getDistanceMetrics`) read.
4. `reset!` again for the next scan.
"""
mutable struct SparseArrayFused{Ti<:Integer, T<:AbstractFloat} <: AbstractSparseDesignMatrix{Ti, T}
    n_vals::Int64
    m::Int64
    n::Int64
    rowval::Vector{Ti}
    nzval::Vector{T}
    x::Vector{T}
    meta::Vector{UInt8}
    rank::Vector{UInt8}
    colptr::Vector{Ti}
end

# Packed-meta layout: bit 0 = matched, bits 1-7 = isotope index (0..127).
@inline pack_meta(matched::Bool, iso::UInt8) =
    UInt8(matched) | (iso << UInt8(1))
@inline unpack_matched(m::UInt8) = (m & UInt8(0x01)) != UInt8(0)
@inline unpack_isotope(m::UInt8) = m >> UInt8(1)

# Override the abstract accessors for the packed layout.
@inline matched_at(H::SparseArrayFused, i::Integer) = unpack_matched(H.meta[i])
@inline isotope_at(H::SparseArrayFused, i::Integer) = unpack_isotope(H.meta[i])
@inline rank_at(H::SparseArrayFused, i::Integer) = H.rank[i]

SparseArrayFused(N::I) where {I<:Integer} = SparseArrayFused(
    0, 0, 0,
    zeros(I, N),        # rowval
    zeros(Float32, N),  # nzval
    zeros(Float32, N),  # x
    fill(UInt8(1), N),  # meta — default matched=true, iso=0 (mirrors SparseArray)
    zeros(UInt8, N),    # rank
    zeros(I, N),        # colptr
)

function reset!(sa::SparseArrayFused{Ti,T}) where {Ti<:Integer,T<:AbstractFloat}
    @inbounds for i in 1:sa.n_vals
        sa.rowval[i] = zero(Ti)
        sa.x[i]      = zero(T)
        sa.nzval[i]  = zero(T)
        sa.meta[i]   = UInt8(1)  # matched=true, iso=0
        sa.rank[i]   = UInt8(0)
    end
    # colptr is indexed by column count, not entry count; zero up to (n + 1).
    @inbounds for i in 1:(sa.n + 1)
        sa.colptr[i] = zero(Ti)
    end
    sa.n_vals = 0
    sa.m = 0
    sa.n = 0
    return
end

"""
    grow_if_needed!(sa, needed_entries)

Expand the entry-indexed vectors (`rowval/nzval/x/meta`) so they can hold
at least `needed_entries` entries. No-op if already large enough.
"""
function grow_if_needed!(sa::SparseArrayFused{Ti,T}, needed_entries::Integer) where {Ti,T}
    cap = length(sa.rowval)
    if needed_entries > cap
        new_cap = max(needed_entries, cap * 2)
        extra = new_cap - cap
        append!(sa.rowval, zeros(Ti, extra))
        append!(sa.nzval,  zeros(T,  extra))
        append!(sa.x,      zeros(T,  extra))
        append!(sa.meta,   fill(UInt8(1), extra))
        append!(sa.rank,   zeros(UInt8, extra))
    end
    return
end

"""
    grow_colptr_if_needed!(sa, needed_cols)

Expand `colptr` so it can index at least `needed_cols` columns (i.e. has
length ≥ `needed_cols + 1`).
"""
function grow_colptr_if_needed!(sa::SparseArrayFused{Ti,T}, needed_cols::Integer) where {Ti,T}
    if length(sa.colptr) < needed_cols + 1
        append!(sa.colptr, zeros(Ti, needed_cols + 1 - length(sa.colptr)))
    end
    return
end
