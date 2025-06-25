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
    rename_column(old::Symbol, new::Symbol)

Rename a column from old to new name.
"""
function rename_column(old::Symbol, new::Symbol)
    desc = "rename_$(old)_to_$(new)"
    op = function(df)
        rename!(df, old => new)
        return df
    end
    return desc => op
end

"""
    select_columns(cols::Vector{Symbol})

Keep only specified columns, ignoring any that don't exist.
"""
function select_columns(cols::Vector{Symbol})
    desc = "select_columns"
    op = function(df)
        available = intersect(cols, Symbol.(names(df)))
        select!(df, available)
        return df
    end
    return desc => op
end

"""
    remove_columns(cols::Symbol...)

Remove specified columns from the DataFrame.
"""
function remove_columns(cols::Symbol...)
    desc = "remove_columns_$(join(cols, '_'))"
    op = function(df)
        select!(df, Not(collect(cols)))
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

"""
    sort_by(cols::Vector{Symbol}; rev::Vector{Bool}=fill(false, length(cols)))

Sort DataFrame by specified columns. Includes post-action to update FileReference sort state.
"""
function sort_by(cols::Vector{Symbol}; rev::Vector{Bool}=fill(false, length(cols)))
    desc = "sort_by_$(join(cols, '_'))"
    op = function(df)
        sort!(df, cols, rev=rev, alg=QuickSort)
        return df
    end
    # Create post-action to mark the file as sorted
    post_action = ref -> mark_sorted!(ref, cols...)
    return PipelineOperation(desc => op, post_action)
end

"""
    add_interpolated_column(new_col::Symbol, source_col::Symbol, 
                           interpolator::Interpolations.Extrapolation)

Add a new column by applying an interpolation function to an existing column.
Commonly used for adding q-value columns based on probability scores.

Example:
```julia
pipeline = TransformPipeline() |>
    add_interpolated_column(:global_qval, :global_prob, global_qval_interp)
```
"""
function add_interpolated_column(new_col::Symbol, source_col::Symbol, 
                               interpolator::Interpolations.Extrapolation)
    desc = "add_interpolated_column($new_col from $source_col)"
    op = function(df)
        # Extract source column with type assertion for performance
        source_data = df[!, source_col]::AbstractVector{Float32}
        
        # Apply interpolator with pre-allocated result vector
        result = Vector{Float32}(undef, length(source_data))
        for i in eachindex(source_data)
            result[i] = Float32(interpolator(source_data[i]))
        end
        
        df[!, new_col] = result
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
        
        # Filter using the mask (more efficient than row-by-row)
        df_filtered = df[keep_mask, :]
        empty!(df)
        append!(df, df_filtered)
        return df
    end
    
    return desc => op
end

"""
    filter_by_multiple_thresholds(conditions::Vector{<:Tuple{Symbol, <:Real}};
                                 comparison::Symbol = :<=,
                                 logic::Symbol = :and)

Apply multiple threshold filters with specified logic (AND/OR).

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
```
"""
function filter_by_multiple_thresholds(conditions::Vector{<:Tuple{Symbol, <:Real}};
                                     comparison::Symbol = :<=,
                                     logic::Symbol = :and)
    desc = "filter_by_multiple_thresholds($(length(conditions)) conditions, $logic)"
    
    # Create comparison function
    comp_fn = if comparison == :<=
        (x, t) -> x <= t
    elseif comparison == :<
        (x, t) -> x < t
    elseif comparison == :>=
        (x, t) -> x >= t
    elseif comparison == :>
        (x, t) -> x > t
    elseif comparison == :(==)
        (x, t) -> x == t
    elseif comparison == :(!=)
        (x, t) -> x != t
    else
        error("Unsupported comparison: $comparison")
    end
    
    op = function(df)
        n_rows = nrow(df)
        keep_mask = Vector{Bool}(undef, n_rows)
        
        if logic == :and
            # AND logic - all conditions must be true
            fill!(keep_mask, true)  # Start with all true
            for (col, threshold) in conditions
                col_data = df[!, col]::AbstractVector{Float32}
                for i in 1:n_rows
                    if keep_mask[i]  # Only check if still true
                        keep_mask[i] = comp_fn(col_data[i], threshold)
                    end
                end
            end
        elseif logic == :or
            # OR logic - at least one condition must be true
            fill!(keep_mask, false)  # Start with all false
            for (col, threshold) in conditions
                col_data = df[!, col]::AbstractVector{Float32}
                for i in 1:n_rows
                    if !keep_mask[i]  # Only check if still false
                        keep_mask[i] = comp_fn(col_data[i], threshold)
                    end
                end
            end
        else
            error("Unsupported logic: $logic. Use :and or :or")
        end
        
        # Filter using the mask
        df_filtered = df[keep_mask, :]
        empty!(df)
        append!(df, df_filtered)
        return df
    end
    
    return desc => op
end