"""
    fast_df_sort!(df::DataFrame, cols::NTuple{N, Symbol};
                  rev::NTuple{N, Bool}=ntuple(_ -> false, Val(N))) where N

Sort a DataFrame in-place using zip+sortperm, which avoids DataFrames' DFPerm
comparison wrapper overhead. ~3.5x faster than `sort!` on large tables (25M rows:
1.36s vs 4.85s, 1.1GB vs 6GB alloc).

Dispatches on `NTuple{N}` so the compiler knows the tuple length at compile time,
producing type-stable `zip` keys without manual branching per column count.
For reversed columns, negates numeric values or flips Bools so ascending sort
produces descending order. Falls back to `sort!` only for unsupported rev types.
"""
function fast_df_sort!(df::DataFrame, cols::NTuple{N, Symbol};
                       rev::NTuple{N, Bool}=ntuple(_ -> false, Val(N))) where N
    nrow(df) == 0 && return df

    # Fall back if a reversed column has an unsupported type (e.g. String)
    if !_can_fast_rev(df, cols, rev)
        sort!(df, collect(cols), rev=collect(rev))
        return df
    end

    perm = if N == 1
        # Single column: sortperm supports rev directly, no zip needed
        sortperm(df[!, cols[1]], rev=rev[1])
    else
        # ntuple with Val(N) builds a concrete-typed tuple at compile time,
        # so zip(...) produces Vector{Tuple{T1,T2,...}} with known types
        key_cols = ntuple(i -> _maybe_rev_col(df[!, cols[i]], rev[i]), Val(N))
        sortperm(collect(zip(key_cols...)))
    end

    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
    return df
end

# Convenience: accept Vector{Symbol} from callers that build column lists dynamically.
# Julia's ntuple(f, n::Int) specializes for small n (≤10), returning a concrete
# NTuple{n} that dispatches to the typed method above.
function fast_df_sort!(df::DataFrame, cols::Vector{Symbol};
                       rev::Vector{Bool}=fill(false, length(cols)))
    col_t = ntuple(i -> cols[i], length(cols))
    rev_t = ntuple(i -> rev[i], length(cols))
    fast_df_sort!(df, col_t; rev=rev_t)
end

# --- Helpers for reverse-sort negation ---

_supports_rev(::Type{<:Signed}) = true
_supports_rev(::Type{<:AbstractFloat}) = true
_supports_rev(::Type{<:Unsigned}) = true
_supports_rev(::Type{Bool}) = true
_supports_rev(::Type) = false

function _can_fast_rev(df::DataFrame, cols::NTuple{N, Symbol}, rev::NTuple{N, Bool}) where N
    for i in 1:N
        rev[i] || continue
        _supports_rev(eltype(df[!, cols[i]])) || return false
    end
    return true
end

# Return column as-is if not reversed, or negate for reverse sorting
_maybe_rev_col(col, do_rev::Bool) = do_rev ? _rev_val.(col) : col

# Negate values so ascending sort produces descending order
_rev_val(x::Union{Signed, AbstractFloat}) = -x
_rev_val(x::Unsigned) = ~x  # bitwise complement reverses unsigned ordering
_rev_val(x::Bool) = !x
