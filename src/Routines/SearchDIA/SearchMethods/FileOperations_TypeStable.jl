"""
Type-stable, high-performance file merging operations for Arrow files.

This module provides optimized merging functions that:
1. Use parametric types for type stability
2. Implement recursive merging to handle large numbers of files
3. Minimize memory usage and file handle requirements
4. Provide specialized operations for different column types

Key improvements over basic FileOperations:
- Type-stable heap operations
- Specialized fillColumn! methods for each type
- Recursive chunking for scalability
- Memory-efficient batch processing
"""

using Arrow, DataFrames, Tables
using DataStructures: BinaryMinHeap, BinaryMaxHeap

# Import FileReferences if not already loaded
if !@isdefined(FileReference)
    include(joinpath(@__DIR__, "FileReferences.jl"))
end

#==========================================================
Type-Stable Column Filling Operations
==========================================================#

"""
Specialized fillColumn! for numeric types (Float32, Float64, Int32, etc.)
"""
function fillColumn!(
    batch_col::Vector{T},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
) where {T<:Real}
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::T
    end
end

"""
Specialized fillColumn! for nullable numeric types
"""
function fillColumn!(
    batch_col::Vector{Union{Missing, T}},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
) where {T<:Real}
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::Union{Missing, T}
    end
end

"""
Specialized fillColumn! for String types
"""
function fillColumn!(
    batch_col::Vector{String},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
)
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::String
    end
end

"""
Specialized fillColumn! for nullable String types
"""
function fillColumn!(
    batch_col::Vector{Union{Missing, String}},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
)
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::Union{Missing, String}
    end
end

"""
Specialized fillColumn! for Boolean types
"""
function fillColumn!(
    batch_col::Vector{Bool},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
)
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::Bool
    end
end

"""
Specialized fillColumn! for UInt8 types (common for entrapment_group_id)
"""
function fillColumn!(
    batch_col::Vector{UInt8},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
)
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::UInt8
    end
end

#==========================================================
Type-Stable DataFrame Creation
==========================================================#

"""
Create a type-stable empty DataFrame with pre-allocated vectors.
Types are determined from the schema of the first table.
"""
function create_typed_dataframe(reference_table::Arrow.Table, batch_size::Int)
    df = DataFrame()
    col_names = Tables.columnnames(reference_table)
    
    for col_name in col_names
        col_type = eltype(Tables.getcolumn(reference_table, col_name))
        
        # Create appropriately typed vector
        if col_type <: Union && Missing <: col_type
            # Handle Union{Missing, T} types
            non_missing_types = Base.uniontypes(col_type)
            actual_type = first(t for t in non_missing_types if t !== Missing)
            df[!, col_name] = Vector{Union{Missing, actual_type}}(undef, batch_size)
        else
            df[!, col_name] = Vector{col_type}(undef, batch_size)
        end
    end
    
    return df
end

#==========================================================
Type-Stable Heap Operations
==========================================================#

"""
Add entry to heap with compile-time known types for sort keys.
"""
function add_to_typed_heap!(
    heap::H,
    table::Arrow.Table,
    table_idx::Int,
    row_idx::Int,
    sort_key1::Symbol,
    sort_key2::Symbol
) where {T1, T2, H<:Union{BinaryMinHeap{Tuple{T1, T2, Int64}}, BinaryMaxHeap{Tuple{T1, T2, Int64}}}}
    val1 = Tables.getcolumn(table, sort_key1)[row_idx]::T1
    val2 = Tables.getcolumn(table, sort_key2)[row_idx]::T2
    push!(heap, (val1, val2, table_idx))
end


#==========================================================
High-Performance Type-Stable Merge
==========================================================#

"""
    stream_sorted_merge_typed(refs::Vector{<:FileReference}, 
                             output_path::String,
                             sort_key1::Symbol, sort_key2::Symbol;
                             reverse::Vector{Bool}=[false, false],
                             batch_size::Int=1_000_000)

Type-stable merge function that determines types at compile time and uses
specialized operations for maximum performance.

The function automatically detects the types of the sort keys from the first
file and dispatches to type-stable code paths.
"""
function stream_sorted_merge_typed(
    refs::Vector{<:FileReference}, 
    output_path::String,
    sort_key1::Symbol, 
    sort_key2::Symbol;
    reverse::Vector{Bool}=[false, false],
    batch_size::Int=1_000_000
)
    # Validate inputs
    isempty(refs) && error("No files to merge")
    length(reverse) != 2 && error("reverse must have exactly 2 elements for 2-key sort")
    
    # Determine types from first file
    first_table = Arrow.Table(file_path(first(refs)))
    T1 = eltype(Tables.getcolumn(first_table, sort_key1))
    T2 = eltype(Tables.getcolumn(first_table, sort_key2))
    
    # Dispatch to type-stable implementation
    return _stream_sorted_merge_typed_impl(
        refs, output_path, sort_key1, sort_key2, 
        T1, T2, reverse, batch_size
    )
end

"""
Internal type-stable implementation with concrete types.
"""
function _stream_sorted_merge_typed_impl(
    refs::Vector{<:FileReference}, 
    output_path::String,
    sort_key1::Symbol, 
    sort_key2::Symbol,
    ::Type{T1},
    ::Type{T2},
    reverse::Vector{Bool},
    batch_size::Int
) where {T1, T2}
    
    # Validate all files exist and have compatible schemas
    for ref in refs
        validate_exists(ref)
    end
    
    # Load all tables (this is the main bottleneck for many files)
    tables = [Arrow.Table(file_path(ref)) for ref in refs]
    table_sizes = [length(Tables.getcolumn(table, 1)) for table in tables]
    table_indices = ones(Int64, length(tables))
    
    # Create type-stable batch DataFrame and heap
    batch_df = create_typed_dataframe(first(tables), batch_size)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, batch_size)
    
    # Choose heap type based on reverse flags
    # For mixed reverse, we need to handle this more carefully
    use_max_heap = all(reverse)  # Only use max heap if both keys are reversed
    
    if use_max_heap
        heap = BinaryMaxHeap{Tuple{T1, T2, Int64}}()
    else
        heap = BinaryMinHeap{Tuple{T1, T2, Int64}}()
    end
    
    # Initialize heap with first row from each table
    for (i, table) in enumerate(tables)
        if table_sizes[i] > 0
            add_to_typed_heap!(heap, table, i, 1, sort_key1, sort_key2)
        end
    end
    
    # Perform merge
    row_idx = 1
    n_writes = 0
    
    while !isempty(heap)
        # Pop minimum element
        _, _, table_idx = pop!(heap)
        current_row_idx = table_indices[table_idx]
        
        # Store row location
        sorted_tuples[row_idx] = (table_idx, current_row_idx)
        
        # Advance to next row in this table
        table_indices[table_idx] += 1
        next_row_idx = table_indices[table_idx]
        
        # Add next row to heap if available
        if next_row_idx <= table_sizes[table_idx]
            add_to_typed_heap!(heap, tables[table_idx], table_idx, next_row_idx, sort_key1, sort_key2)
        end
        
        row_idx += 1
        
        # Write batch when full
        if row_idx > batch_size
            _fill_batch_columns_typed!(batch_df, tables, sorted_tuples, batch_size)
            _write_batch_typed(output_path, batch_df, n_writes, batch_size)
            n_writes += 1
            row_idx = 1
        end
    end
    
    # Write final partial batch
    if row_idx > 1
        final_rows = row_idx - 1
        _fill_batch_columns_typed!(batch_df, tables, sorted_tuples, final_rows)
        _write_batch_typed(output_path, batch_df, n_writes, final_rows)
    end
    
    # Create output reference
    output_ref = create_reference(output_path, typeof(first(refs)))
    
    # If we have mixed reverse flags, we need to sort the final result
    if !all(reverse) && any(reverse)
        # For mixed reverse, we need to sort the final file
        # This is a more complex case that requires post-processing
        sort_file_by_keys!(output_ref, sort_key1, sort_key2, reverse=reverse)
    elseif !use_max_heap
        # For normal ascending sort, mark as sorted
        mark_sorted!(output_ref, sort_key1, sort_key2)
    else
        # For descending sort with max heap, mark as reverse sorted
        mark_sorted!(output_ref, sort_key1, sort_key2)
    end
    
    return output_ref
end

"""
Fill batch DataFrame columns using type-stable fillColumn! methods.
"""
function _fill_batch_columns_typed!(
    batch_df::DataFrame, 
    tables::Vector{Arrow.Table}, 
    sorted_tuples::Vector{Tuple{Int64, Int64}}, 
    n_rows::Int
)
    for col_name in names(batch_df)
        col_symbol = Symbol(col_name)
        fillColumn!(batch_df[!, col_symbol], col_symbol, sorted_tuples, tables, n_rows)
    end
end

"""
Write batch to file with proper error handling.
"""
function _write_batch_typed(
    output_path::String, 
    batch_df::DataFrame, 
    n_writes::Int, 
    n_rows::Int
)
    # Only write the rows we actually filled
    data_to_write = n_rows < nrow(batch_df) ? batch_df[1:n_rows, :] : batch_df
    
    if n_writes == 0
        # First write
        if isfile(output_path)
            rm(output_path)
        end
        open(output_path, "w") do io
            Arrow.write(io, data_to_write; file=false)
        end
    else
        # Append to existing file
        Arrow.append(output_path, data_to_write)
    end
end

#==========================================================
Recursive Merge Strategy
==========================================================#

"""
    recursive_merge(refs::Vector{<:FileReference}, 
                   output_path::String,
                   sort_key1::Symbol, sort_key2::Symbol;
                   reverse::Vector{Bool}=[false, false],
                   group_size::Int=10,
                   batch_size::Int=1_000_000)

Recursively merge large numbers of files by processing them in small groups.
This prevents file handle exhaustion and reduces memory usage.

Algorithm:
1. If files ≤ group_size: direct merge
2. Otherwise: merge in chunks of group_size → intermediate files
3. Recursively merge intermediate files
4. Clean up intermediate files

Example:
- 300 files → 30 intermediate files (groups of 10)
- 30 files → 3 intermediate files (groups of 10)  
- 3 files → 1 final file
"""
function recursive_merge(
    refs::Vector{<:FileReference}, 
    output_path::String,
    sort_key1::Symbol, 
    sort_key2::Symbol;
    reverse::Vector{Bool}=[false, false],
    group_size::Int=10,
    batch_size::Int=1_000_000,
    cleanup_intermediates::Bool=true
)
    # Base case: small enough for direct merge
    if length(refs) <= group_size
        @info "Direct merge of $(length(refs)) files"
        return stream_sorted_merge_typed(
            refs, output_path, sort_key1, sort_key2; 
            reverse=reverse, batch_size=batch_size
        )
    end
    
    # Recursive case: merge in chunks
    @info "Recursive merge: $(length(refs)) files → groups of $group_size"
    
    # Create temporary directory for intermediate files
    temp_dir = mktempdir(prefix="merge_intermediate_")
    intermediate_refs = FileReference[]
    
    try
        # Process chunks
        chunk_num = 0
        for chunk in Iterators.partition(refs, group_size)
            chunk_num += 1
            intermediate_path = joinpath(temp_dir, "intermediate_$(chunk_num).arrow")
            
            @info "Processing chunk $chunk_num: $(length(chunk)) files"
            
            intermediate_ref = stream_sorted_merge_typed(
                collect(chunk), intermediate_path, sort_key1, sort_key2;
                reverse=reverse, batch_size=batch_size
            )
            
            push!(intermediate_refs, intermediate_ref)
        end
        
        @info "Created $(length(intermediate_refs)) intermediate files"
        
        # Recursively merge intermediate files
        result = recursive_merge(
            intermediate_refs, output_path, sort_key1, sort_key2;
            reverse=reverse, group_size=group_size, batch_size=batch_size,
            cleanup_intermediates=false  # We'll clean up manually
        )
        
        return result
        
    finally
        # Clean up intermediate files
        if cleanup_intermediates && isdir(temp_dir)
            @info "Cleaning up intermediate files in $temp_dir"
            rm(temp_dir, recursive=true)
        end
    end
end

#==========================================================
Convenience Functions
==========================================================#

"""
    merge_protein_groups_high_performance(refs::Vector{ProteinGroupFileReference},
                                         output_path::String;
                                         use_recursive::Bool=true,
                                         group_size::Int=10)

High-performance merge specifically optimized for protein group files.
Uses the standard protein group sort keys: :pg_score, :target with reverse=[true, true].
"""
function merge_protein_groups_high_performance(
    refs::Vector{ProteinGroupFileReference},
    output_path::String;
    use_recursive::Bool=true,
    group_size::Int=10,
    batch_size::Int=1_000_000
)
    sort_keys = (:pg_score, :target)
    reverse = [true, true]  # Descending score, descending target (targets first)
    
    if use_recursive && length(refs) > group_size
        return recursive_merge(
            refs, output_path, sort_keys[1], sort_keys[2];
            reverse=reverse, group_size=group_size, batch_size=batch_size
        )
    else
        return stream_sorted_merge_typed(
            refs, output_path, sort_keys[1], sort_keys[2];
            reverse=reverse, batch_size=batch_size
        )
    end
end

# Export main functions
export stream_sorted_merge_typed, recursive_merge, merge_protein_groups_high_performance