"""
    fast_df_sort!(df::DataFrame, cols::Vector{Symbol};
                  rev::Vector{Bool}=fill(false, length(cols)))

Sort a DataFrame in-place using zip+sortperm, which avoids DataFrames' DFPerm
comparison wrapper overhead. ~3.5x faster than `sort!` on large tables (25M rows:
1.36s vs 4.85s, 1.1GB vs 6GB alloc).

For 1-4 sort keys, builds concrete tuple keys via `zip` and uses `sortperm`.
For reversed columns, negates numeric values or flips Bools so ascending sort
produces descending order. Falls back to `sort!` for >4 keys or unsupported types.
"""
function fast_df_sort!(df::DataFrame, cols::Vector{Symbol};
                       rev::Vector{Bool}=fill(false, length(cols)))
    ncols = length(cols)
    nrow(df) == 0 && return df

    # Fallback for >4 keys
    if ncols > 4
        sort!(df, cols, rev=rev)
        return df
    end

    # Check that reversed columns have types we can negate
    if !_can_fast_rev(df, cols, rev)
        sort!(df, cols, rev=rev)
        return df
    end

    perm = if ncols == 1
        # Single column: sortperm supports rev directly
        sortperm(df[!, cols[1]], rev=rev[1])
    elseif ncols == 2
        c1 = _maybe_rev_col(df[!, cols[1]], rev[1])
        c2 = _maybe_rev_col(df[!, cols[2]], rev[2])
        sortperm(collect(zip(c1, c2)))
    elseif ncols == 3
        c1 = _maybe_rev_col(df[!, cols[1]], rev[1])
        c2 = _maybe_rev_col(df[!, cols[2]], rev[2])
        c3 = _maybe_rev_col(df[!, cols[3]], rev[3])
        sortperm(collect(zip(c1, c2, c3)))
    else  # ncols == 4
        c1 = _maybe_rev_col(df[!, cols[1]], rev[1])
        c2 = _maybe_rev_col(df[!, cols[2]], rev[2])
        c3 = _maybe_rev_col(df[!, cols[3]], rev[3])
        c4 = _maybe_rev_col(df[!, cols[4]], rev[4])
        sortperm(collect(zip(c1, c2, c3, c4)))
    end

    # Apply permutation to all columns
    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
    return df
end

# --- Helpers for reverse-sort negation ---

# Types that support negation for reverse sorting
_supports_rev(::Type{<:Signed}) = true
_supports_rev(::Type{<:AbstractFloat}) = true
_supports_rev(::Type{<:Unsigned}) = true
_supports_rev(::Type{Bool}) = true
_supports_rev(::Type) = false

function _can_fast_rev(df::DataFrame, cols::Vector{Symbol}, rev::Vector{Bool})
    for (col, r) in zip(cols, rev)
        r || continue  # non-reversed columns always work
        _supports_rev(eltype(df[!, col])) || return false
    end
    return true
end

# Return column as-is if not reversed, or negate for reverse sorting
_maybe_rev_col(col, do_rev::Bool) = do_rev ? _rev_val.(col) : col

# Negate values so ascending sort produces descending order
_rev_val(x::Union{Signed, AbstractFloat}) = -x
_rev_val(x::Unsigned) = ~x  # bitwise complement reverses unsigned ordering
_rev_val(x::Bool) = !x
