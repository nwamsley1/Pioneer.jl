"""
Column manipulation operations for file processing.

Provides utilities for adding, updating, and modifying columns
in files with memory-efficient streaming support.
"""

using Arrow, DataFrames, Tables

# Dependencies from core modules:
# validate_exists, has_column, schema, row_count, file_path, create_reference, mark_sorted!, sorted_by
# writeArrow (from utils), sort_file_by_keys! (defined later)

#==========================================================
Column Operations
==========================================================#

"""
    add_column_to_file!(ref::FileReference, col_name::Symbol, 
                       compute_fn::Function; batch_size=100_000)

Add a new column to a file by computing it from existing columns.
Updates the file in place and updates the reference's schema.
"""
function add_column_to_file!(ref::FileReference, 
                           col_name::Symbol,
                           compute_fn::Function;
                           batch_size::Int=100_000)
    validate_exists(ref)
    
    # Check if column already exists
    if has_column(schema(ref), col_name)
        @warn "Column $col_name already exists, will be overwritten"
    end
    
    # Create temporary output file
    temp_path = file_path(ref) * ".tmp"
    
    # For small files or if partitioner has issues, process entire file at once
    tbl = Arrow.Table(file_path(ref))
    
    # Check if file is small enough to process at once
    if row_count(ref) <= batch_size
        # Process entire file at once
        df_batch = DataFrame(Tables.columntable(tbl))
        df_batch[!, col_name] = compute_fn(df_batch)
        writeArrow(temp_path, df_batch)
    else
        # Stream through file for larger files
        partitions = Tables.partitioner(tbl, batch_size)
        
        first_write = true
        for partition in partitions
            df_batch = DataFrame(Tables.columntable(partition))
            df_batch[!, col_name] = compute_fn(df_batch)
            
            if first_write
                writeArrow(temp_path, df_batch)
                first_write = false
            else
                Arrow.append(temp_path, df_batch)
            end
        end
    end
    
    # Replace original file
    mv(temp_path, file_path(ref), force=true)
    
    # Update reference metadata
    # We need to re-read the file to update schema
    new_ref = create_reference(file_path(ref), typeof(ref))
    ref.schema = schema(new_ref)
    ref.row_count = row_count(new_ref)
    ref.sorted_by = ()  # Adding column may invalidate sort
    
    return ref
end

"""
    update_column_in_file!(ref::FileReference, col_name::Symbol,
                          update_fn::Function; batch_size=100_000)

Update an existing column in a file.
"""
function update_column_in_file!(ref::FileReference,
                               col_name::Symbol,
                               update_fn::Function;
                               batch_size::Int=100_000)
    validate_exists(ref)
    
    if !has_column(schema(ref), col_name)
        error("Column $col_name does not exist")
    end
    
    # Create temporary output file
    temp_path = file_path(ref) * ".tmp"
    
    # For small files or if partitioner has issues, process entire file at once
    tbl = Arrow.Table(file_path(ref))
    
    # Check if file is small enough to process at once
    if row_count(ref) <= batch_size
        # Process entire file at once
        df_batch = DataFrame(Tables.columntable(tbl))
        df_batch[!, col_name] = update_fn(df_batch[!, col_name])
        writeArrow(temp_path, df_batch)
    else
        # Stream through file for larger files
        partitions = Tables.partitioner(tbl, batch_size)
        
        first_write = true
        for partition in partitions
            df_batch = DataFrame(Tables.columntable(partition))
            df_batch[!, col_name] = update_fn(df_batch[!, col_name])
            
            if first_write
                writeArrow(temp_path, df_batch)
                first_write = false
            else
                Arrow.append(temp_path, df_batch)
            end
        end
    end
    
    # Replace original file
    mv(temp_path, file_path(ref), force=true)
    
    # Sort order may be affected if we updated a sort key
    if col_name in sorted_by(ref)
        mark_sorted!(ref)  # Empty tuple - not sorted
    end
    
    return ref
end

"""
    add_column_and_sort!(ref::FileReference, col_name::Symbol, 
                        compute_fn::Function, sort_keys::Symbol...; 
                        reverse::Bool=false) -> FileReference
    
Add a column and immediately sort by specified keys.
"""
function add_column_and_sort!(ref::FileReference, col_name::Symbol,
                             compute_fn::Function, sort_keys::Symbol...;
                             reverse::Bool=false)
    # Add column
    add_column_to_file!(ref, col_name, compute_fn)
    
    # Sort by keys
    sort_file_by_keys!(ref, sort_keys...; reverse=reverse)
    
    return ref
end

# Export column functions
export add_column_to_file!, update_column_in_file!, add_column_and_sort!