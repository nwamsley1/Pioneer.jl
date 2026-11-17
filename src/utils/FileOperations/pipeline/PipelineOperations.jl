# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
Pipeline operation builders.

Provides convenient builders for common file operations
that can be composed into transformation pipelines.
"""

using DataFrames
using Interpolations

#==========================================================
Pipeline Operation Builders
==========================================================#

"""
    add_column(name::Symbol, compute_fn::Function)

Add a new column by applying compute_fn to each row.
"""
function add_column(name::Symbol, compute_fn::Function)
    desc = "add_column_$name"
    op = function(df)
        # Pre-allocate result vector for performance
        n_rows = nrow(df)
        result = Vector{Any}(undef, n_rows)
        
        # Use indexed iteration instead of eachrow for better type stability
        for i in 1:n_rows
            # Create a NamedTuple row for the compute function
            row_data = NamedTuple{Tuple(propertynames(df))}(Tuple(df[i, col] for col in propertynames(df)))
            result[i] = compute_fn(row_data)
        end
        
        df[!, name] = result
        return df
    end
    return desc => op
end


"""
    filter_rows(predicate::Function; desc::String="filter")

Filter rows based on a predicate function.
"""
function filter_rows(predicate::Function; desc::String="filter")
    op = function(df)
        filter!(predicate, df)
        return df
    end
    return desc => op
end

# `column => predicate` form: the predicate sees the typed column's values, not a DataFrameRow.
function filter_rows(column_predicate::Pair{Symbol, <:Function}; desc::String="filter")
    op = function(df)
        filter!(column_predicate, df)
        return df
    end
    return desc => op
end

"""
    sort_by(cols::Vector{Symbol}; rev::Vector{Bool}=fill(false, length(cols)))

Sort DataFrame by specified columns. Includes post-action to update FileReference sort state.
"""
function sort_by(cols::Vector{Symbol}; rev::Vector{Bool}=fill(false, length(cols)))
    desc = "sort_by_$(join(cols, '_'))"
    op = function(df)
        fast_df_sort!(df, cols, rev=rev)
        return df
    end
    # Create post-action to mark the file as sorted
    post_action = ref -> mark_sorted!(ref, cols...)
    return PipelineOperation(desc => op, post_action)
end

"""
    add_interpolated_column(new_col::Symbol, source_col::Symbol, 
                           interpolator)

Add a new column by applying a callable interpolation or score-calibration mapping to an existing column.
Commonly used for adding q-value columns based on probability scores.

Example:
```julia
pipeline = TransformPipeline() |>
    add_interpolated_column(:global_qval, :global_prob, global_qval_interp)
```
"""
# Function barrier: `df[!, col]` is type-unstable (returns an abstract column type), so indexing it
# inside a per-row loop boxes every element (~100 B/row measured). Passing the column to a method
# specializes on its concrete runtime type, making the loop allocation-free. `1:length` (not
# `eachindex`) so the result Vector and an Arrow-backed source stay index-aligned (see A1).
@inline function _fill_interpolated!(result::Vector{Float32}, source::AbstractVector, interp)
    @inbounds for i in 1:length(source)
        result[i] = Float32(interp(source[i]))
    end
end

function add_interpolated_column(new_col::Symbol, source_col::Symbol,
                               interpolator)
    desc = "add_interpolated_column($new_col from $source_col)"
    op = function(df)
        source_data = df[!, source_col]
        result = Vector{Float32}(undef, length(source_data))
        _fill_interpolated!(result, source_data, interpolator)
        df[!, new_col] = result
        return df
    end
    return desc => op
end

"""
    add_dict_column_composite_key(new_col::Symbol, key_cols::Vector{Symbol}, lookup_dict::Dict{Tuple{K,Bool,UInt8},V}) where {K,V}

Add a new column by looking up values in a dictionary using multiple key columns to form a composite key.
Specifically designed for protein groups which are identified by (protein_name, target, entrap_id).

# Arguments
- `new_col`: Name of the new column to create
- `key_cols`: Vector of column names to use for creating composite keys (must be [:protein_name, :target, :entrap_id])
- `lookup_dict`: Dictionary mapping tuple keys (protein name or pg_id, Bool, UInt8) to values

# Example
```julia
pipeline = TransformPipeline() |>
    add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], pg_qval_dict)
```
"""
function add_dict_column_composite_key(new_col::Symbol, key_cols::Vector{Symbol},
                                      lookup_dict::Dict{Tuple{K,Bool,UInt8},V}) where {K,V}
    desc = "add_dict_column_composite_key($new_col from $(join(key_cols, ", ")))"
    op = function(df)
        # Verify key columns exist
        for col in key_cols
            if !hasproperty(df, col)
                error("Column $col not found in DataFrame for composite key lookup")
            end
        end

        df[!, new_col] = _lookup_composite_keys(lookup_dict, V, df.protein_name, df.target, df.entrap_id)
        return df
    end
    return desc => op
end

"""
    filter_by_threshold(col::Symbol, threshold::Real; comparison::Symbol = :<=)

Filter rows where column meets threshold condition.

Supported comparisons: :<=, :<, :>=, :>, :==, :!=

Example:
```julia
pipeline = TransformPipeline() |>
    filter_by_threshold(:qval, 0.01)  # Keep rows where qval <= 0.01
```
"""
function filter_by_threshold(col::Symbol, threshold::Real; comparison::Symbol = :<=)
    desc = "filter_by_threshold($col $comparison $threshold)"
    
    # Create comparison function
    comp_fn = if comparison == :<=
        (x) -> x <= threshold
    elseif comparison == :<
        (x) -> x < threshold
    elseif comparison == :>=
        (x) -> x >= threshold
    elseif comparison == :>
        (x) -> x > threshold
    elseif comparison == :(==)
        (x) -> x == threshold
    elseif comparison == :(!=)
        (x) -> x != threshold
    else
        error("Unsupported comparison: $comparison")
    end
    
    op = function(df)
        # Extract column with type assertion for performance
        col_data = df[!, col]::AbstractVector{Float32}
        
        # Create boolean mask for filtering
        keep_mask = Vector{Bool}(undef, length(col_data))
        for i in eachindex(col_data)
            keep_mask[i] = comp_fn(col_data[i])
        end
        
        # Filter in place (no copy of the kept rows)
        keepat!(df, keep_mask)
        return df
    end
    
    return desc => op
end

# Type-stable helper functions for AND logic
function update_mask_and!(mask::BitVector, col::AbstractVector{T}, threshold::T, op::Function) where T
    @inbounds for i in eachindex(mask, col)
        mask[i] = mask[i] && op(col[i], threshold)
    end
end

function update_mask_and!(mask::BitVector, col::AbstractVector{Union{Missing,T}}, threshold::T, op::Function) where T
    @inbounds for i in eachindex(mask, col)
        if mask[i]
            val = col[i]
            mask[i] = !ismissing(val) && op(val, threshold)
        end
    end
end

# String specialization for AND
function update_mask_and!(mask::BitVector, col::AbstractVector{<:AbstractString}, threshold::AbstractString, op::Function)
    if op === (==)
        threshold_lower = lowercase(string(threshold))
        @inbounds for i in eachindex(mask, col)
            mask[i] = mask[i] && lowercase(string(col[i])) == threshold_lower
        end
    elseif op === (!=)
        threshold_lower = lowercase(string(threshold))
        @inbounds for i in eachindex(mask, col)
            mask[i] = mask[i] && lowercase(string(col[i])) != threshold_lower
        end
    else
        error("Comparison operator $op not supported for String columns")
    end
end

# Type-stable helper functions for OR logic
function update_mask_or!(mask::BitVector, col::AbstractVector{T}, threshold::T, op::Function) where T
    @inbounds for i in eachindex(mask, col)
        mask[i] = mask[i] || op(col[i], threshold)
    end
end

function update_mask_or!(mask::BitVector, col::AbstractVector{Union{Missing,T}}, threshold::T, op::Function) where T
    @inbounds for i in eachindex(mask, col)
        if !mask[i]
            val = col[i]
            mask[i] = !ismissing(val) && op(val, threshold)
        end
    end
end

# String specialization for OR
function update_mask_or!(mask::BitVector, col::AbstractVector{<:AbstractString}, threshold::AbstractString, op::Function)
    if op === (==)
        threshold_lower = lowercase(string(threshold))
        @inbounds for i in eachindex(mask, col)
            mask[i] = mask[i] || lowercase(string(col[i])) == threshold_lower
        end
    elseif op === (!=)
        threshold_lower = lowercase(string(threshold))
        @inbounds for i in eachindex(mask, col)
            mask[i] = mask[i] || lowercase(string(col[i])) != threshold_lower
        end
    else
        error("Comparison operator $op not supported for String columns")
    end
end

"""
    filter_by_multiple_thresholds(conditions::Vector{<:Tuple{Symbol, <:Any}};
                                 comparison::Union{Symbol,Function} = <=,
                                 logic::Symbol = :and)

Apply multiple threshold filters with specified logic (AND/OR).
Uses DataFrames' efficient `filter!` function to avoid intermediate allocations.

# Arguments
- `conditions`: Vector of (column_name, threshold) tuples
- `comparison`: Comparison operator (default: <=). Can be a Symbol or Function
- `logic`: :and (all conditions must be true) or :or (at least one must be true)

Example:
```julia
# AND logic (both conditions must be true)
pipeline = TransformPipeline() |>
    filter_by_multiple_thresholds([
        (:global_qval, 0.01),
        (:qval, 0.01)
    ])

# OR logic (at least one condition must be true)
pipeline = TransformPipeline() |>
    filter_by_multiple_thresholds([
        (:score1, 0.95),
        (:score2, 0.90)
    ], logic = :or)

# String comparisons
pipeline |> filter_by_multiple_thresholds([
    (:category, "protein"),
    (:status, "valid")
], comparison = ==)
```
"""
function filter_by_multiple_thresholds(conditions::Vector{<:Tuple{Symbol, <:Any}};
                                     comparison::Union{Symbol,Function} = <=,
                                     logic::Symbol = :and)
    desc = "filter_by_multiple_thresholds($(length(conditions)) conditions, $logic)"
    
    # Convert symbol to function if needed
    comp_fn = if comparison isa Symbol
        # Use getfield to get the function from Base
        if comparison in (:<=, :<, :>=, :>, :(==), :(!=))
            getfield(Base, comparison)
        else
            error("Unsupported comparison symbol: $comparison")
        end
    else
        comparison  # Already a function
    end
    
    op = function(df)
        n_rows = nrow(df)
        
        # Pre-allocate mask based on logic
        keep_mask = if logic == :and
            trues(n_rows)
        elseif logic == :or
            falses(n_rows)
        else
            error("Unsupported logic: $logic. Use :and or :or")
        end
        
        # Process each condition column-wise for type stability
        for (col, threshold) in conditions
            col_data = df[!, col]
            
            # Type-stable column processing via multiple dispatch
            if logic == :and
                update_mask_and!(keep_mask, col_data, threshold, comp_fn)
            else  # :or
                update_mask_or!(keep_mask, col_data, threshold, comp_fn)
            end
        end
        
        # Filter in place (no copy of the kept rows)
        keepat!(df, keep_mask)

        return df
    end
    
    return desc => op
end

# Look up each row's (protein_name, target, entrap_id) key (a function barrier: the columns are untyped).
function _lookup_composite_keys(lookup_dict, ::Type{V}, names::AbstractVector, targets::AbstractVector,
                                entrap_ids::AbstractVector) where {V}
    result = Vector{Union{V, Missing}}(undef, length(names))
    for i in eachindex(names, targets, entrap_ids)
        result[i] = get(lookup_dict, (names[i], targets[i], entrap_ids[i]), missing)
    end
    return result
end
