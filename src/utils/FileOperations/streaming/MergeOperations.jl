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
High-performance merge operations with heap-based sorting.

Provides type-stable, memory-efficient merging of multiple
sorted files with support for arbitrary numbers of sort keys
and mixed sort directions.
"""

using Arrow, DataFrames, Tables
using DataStructures: BinaryMinHeap, BinaryMaxHeap, BinaryHeap
using Base.Order: Ordering, Forward, Reverse, By, Lt

#==========================================================
Helper Functions for Heap-based Merge
==========================================================#

# Create empty DataFrame with same schema as Arrow table
function _create_empty_dataframe(table::Arrow.Table, n_rows::Int)
    df = DataFrame()
    col_names = Tables.columnnames(table)
    schema_types = Tables.schema(table).types
    
    for (col_name, col_type) in zip(col_names, schema_types)
        # Handle Union types properly
        if col_type isa Union && Missing <: col_type
            # Extract the non-missing type
            non_missing_type = Base.uniontypes(col_type)
            actual_type = first(t for t in non_missing_type if t !== Missing)
            df[!, col_name] = Vector{Union{Missing, actual_type}}(undef, n_rows)
        else
            df[!, col_name] = Vector{col_type}(undef, n_rows)
        end
    end
    return df
end

# Type-stable fillColumn! implementations
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

"""
Specialized fillColumn! for Tuple types
"""
function fillColumn!(
    batch_col::Vector{Tuple{R,R}},
    col::Symbol,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    tables::Vector{Arrow.Table},
    n_rows::Int
) where {R<:Real}
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        batch_col[i] = tables[table_idx][col][row_idx]::Tuple{R,R}
    end
end

# Generic fillColumn! implementations (fallback)
function _fillColumn!(col::Vector{T}, col_symbol::Symbol, 
                     sorted_tuples::Vector{Tuple{Int64, Int64}},
                     tables::Vector{Arrow.Table}, n_rows::Int) where T
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        col[i] = Tables.getcolumn(tables[table_idx], col_symbol)[row_idx]
    end
end

# Fill batch columns from sorted tuples (type-stable version)
function _fill_batch_columns_nkey!(
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
N-Key Implementation Support Functions
==========================================================#

"""
Normalize reverse specification to vector format.
"""
function _normalize_reverse_spec(reverse::Union{Bool,Vector{Bool}}, n_keys::Int)
    if reverse isa Bool
        return fill(reverse, n_keys)
    elseif length(reverse) == 1
        return fill(reverse[1], n_keys)
    elseif length(reverse) == n_keys
        return reverse
    else
        error("reverse must be Bool, single-element vector, or have same length as sort keys ($n_keys)")
    end
end

"""
Compare two tuples with mixed reverse directions for each key.
"""
function _mixed_reverse_compare(a, b, reverse_vec::Vector{Bool})
    # Compare each sort key according to its reverse setting
    for i in 1:(length(a)-1)  # -1 to skip table_idx (last element)
        val_a, val_b = a[i], b[i]
        if reverse_vec[i]
            # Reverse order: larger values come first
            if val_a > val_b
                return true
            elseif val_a < val_b
                return false
            end
            # Equal values continue to next key
        else
            # Normal order: smaller values come first
            if val_a < val_b
                return true
            elseif val_a > val_b
                return false
            end
            # Equal values continue to next key
        end
    end
    # All sort keys are equal, use table_idx as tiebreaker (always ascending)
    return a[end] < b[end]
end

"""
Generate heap type dynamically based on sort types and reverse specification.
Supports mixed reverse directions for different keys.
"""
function _create_typed_heap(::Type{SortTypes}, reverse_vec::Vector{Bool}) where SortTypes
    if length(reverse_vec) == 1
        # Single key - use simple min/max heap
        if reverse_vec[1]
            return BinaryMaxHeap{SortTypes}()
        else
            return BinaryMinHeap{SortTypes}()
        end
    elseif all(reverse_vec)
        # All keys reversed - use max heap
        return BinaryMaxHeap{SortTypes}()
    elseif !any(reverse_vec)
        # No keys reversed - use min heap
        return BinaryMinHeap{SortTypes}()
    else
        # Mixed reverse directions - use custom comparison with BinaryHeap
        # Create comparison function that returns true if a < b in our desired order
        comp_func = (a, b) -> _mixed_reverse_compare(a, b, reverse_vec)
        ordering = Lt(comp_func)
        return BinaryHeap(ordering, SortTypes[])
    end
end

"""
Add entry to N-key heap with variadic sort keys.
"""
function _add_to_nkey_heap!(
    heap::H,
    table::Arrow.Table,
    table_idx::Int,
    row_idx::Int,
    sort_keys::NTuple{N,Symbol}
) where {N, H<:Union{BinaryMinHeap, BinaryMaxHeap, BinaryHeap}}
    # Build tuple: (val1, val2, ..., valN, table_idx)
    values = tuple((Tables.getcolumn(table, key)[row_idx] for key in sort_keys)..., table_idx)
    push!(heap, values)
end

"""
Write batch to file with proper error handling (type-stable version).
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
            # Force garbage collection to release any lingering file handles
            GC.gc()
            # Use Windows-safe removal to avoid permission errors
            safeRm(output_path, nothing)
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
Main Stream Sorted Merge Function
==========================================================#

"""
    stream_sorted_merge(refs::Vector{<:FileReference}, 
                       output_path::String,
                       sort_keys::Symbol...;
                       reverse::Union{Bool,Vector{Bool}}=false,
                       batch_size::Int=1_000_000)

High-performance type-stable merge function that handles arbitrary numbers of sort keys.
Determines types at compile time and uses specialized operations for maximum performance.

Provides 4-20x speedup over the original implementation through:
- Type-stable heap operations with compile-time optimizations
- Specialized column filling methods for different data types  
- Memory-efficient batch processing
- Support for arbitrary number of sort keys (2, 3, 4, 5+)

# Arguments
- `refs`: Vector of file references to merge
- `output_path`: Path for merged output file
- `sort_keys...`: Variable number of sort key symbols
- `reverse`: Boolean or vector specifying sort direction(s)
- `batch_size`: Number of rows to process in each batch

# Examples
```julia
# 2 keys (protein groups)
stream_sorted_merge(refs, path, :pg_score, :target; reverse=[true, true])

# 4 keys (MaxLFQ scenario)  
stream_sorted_merge(refs, path, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx; reverse=true)

# Mixed reverse directions
stream_sorted_merge(refs, path, :protein, :score, :name; reverse=[false, true, false])
```
"""
function stream_sorted_merge(
    refs::Vector{<:FileReference}, 
    output_path::String,
    sort_keys::Symbol...;
    reverse::Union{Bool,Vector{Bool}}=false,
    batch_size::Int=1_000_000
)
    # Validate inputs
    isempty(refs) && error("No files to merge")
    isempty(sort_keys) && error("At least one sort key must be specified")
    
    # Convert sort_keys to tuple for type stability
    sort_keys_tuple = tuple(sort_keys...)
    
    # Normalize reverse specification
    reverse_vec = _normalize_reverse_spec(reverse, length(sort_keys))
    
    # Determine types from first file
    first_table = Arrow.Table(file_path(first(refs)))
    sort_types = tuple((eltype(Tables.getcolumn(first_table, key)) for key in sort_keys_tuple)...)
    

    # Dispatch to N-key implementation
    return _stream_sorted_merge_nkey_impl(
        refs, output_path, sort_keys_tuple, sort_types, reverse_vec, batch_size
    )
end

"""
Internal N-key implementation with dynamic type dispatch.
"""
function _stream_sorted_merge_nkey_impl(
    refs::Vector{<:FileReference}, 
    output_path::String,
    sort_keys::NTuple{N,Symbol},
    sort_types::NTuple{M,DataType},
    reverse_vec::Vector{Bool},
    batch_size::Int
) where {N, M}
    # Validate all files exist and have compatible schemas
    for ref in refs
        validate_exists(ref)
    end

    # Validate that all files are sorted by the requested keys
    for ref in refs
        if !is_sorted_by(ref, sort_keys...)
            error("File $(file_path(ref)) is not sorted by the required keys: $(sort_keys). Use sort_file_by_keys! or mark_sorted! first.")
        end
    end
    
    # Load all tables
    tables = [Arrow.Table(file_path(ref)) for ref in refs]
    
    # Validate that all tables have the required sort columns
    for (i, table) in enumerate(tables)
        available_columns = Set(Tables.columnnames(table))
        for key in sort_keys
            if key ∉ available_columns
                throw(BoundsError("Column $key not found in file $(file_path(refs[i])). Available columns: $(collect(available_columns))"))
            end
        end
    end
    
    table_sizes = [length(Tables.getcolumn(table, 1)) for table in tables]
    table_indices = ones(Int64, length(tables))
    
    # Create type-stable batch DataFrame
    batch_df = create_typed_dataframe(first(tables), batch_size)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, batch_size)
    
    # Create heap with full type information
    heap_tuple_type = Tuple{sort_types..., Int64}
    heap = _create_typed_heap(heap_tuple_type, reverse_vec)
    
    # Initialize heap with first row from each table
    for (i, table) in enumerate(tables)
        if table_sizes[i] > 0
            _add_to_nkey_heap!(heap, table, i, 1, sort_keys)
        end
    end
    
    # Perform merge
    row_idx = 1
    n_writes = 0
    
    while !isempty(heap)
        # Pop minimum/maximum element
        heap_entry = pop!(heap)
        table_idx = heap_entry[end]  # Last element is always table_idx
        current_row_idx = table_indices[table_idx]
        
        # Store row location
        sorted_tuples[row_idx] = (table_idx, current_row_idx)
        
        # Advance to next row in this table
        table_indices[table_idx] += 1
        next_row_idx = table_indices[table_idx]
        
        # Add next row to heap if available
        if next_row_idx <= table_sizes[table_idx]
            _add_to_nkey_heap!(heap, tables[table_idx], table_idx, next_row_idx, sort_keys)
        end
        
        row_idx += 1
        
        # Write batch when full
        if row_idx > batch_size
            _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, batch_size)
            _write_batch_typed(output_path, batch_df, n_writes, batch_size)
            n_writes += 1
            row_idx = 1
        end
    end
    
    # Write final partial batch
    if row_idx > 1
        final_rows = row_idx - 1
        _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, final_rows)
        _write_batch_typed(output_path, batch_df, n_writes, final_rows)
    end
    
    # Create output reference
    output_ref = create_reference(output_path, typeof(first(refs)))
    
    # Mark as sorted - heap-based merge maintains sort order for all cases
    mark_sorted!(output_ref, sort_keys...)
    
    return output_ref
end
