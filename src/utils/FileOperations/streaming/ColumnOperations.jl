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
Column manipulation operations for file processing.

Provides utilities for adding, updating, and modifying columns
in files with memory-efficient streaming support.
"""

using Arrow, DataFrames, Tables

# Dependencies from core modules:
# validate_exists, has_column, schema, row_count, file_path, create_reference, mark_sorted!, sorted_by
# writeArrow (from utils), sort_file_by_keys! (defined later)

#==========================================================
Windows-Safe File Operations
==========================================================#

"""
    safe_replace_file(temp_path::String, target_path::String, file_handle)

Safely replace a file with Windows-specific handling for file locks and permissions.
Uses the same retry and fallback logic as writeArrow and safeRm.

# Arguments
- `temp_path`: Path to temporary file that should replace the target
- `target_path`: Path to file that should be replaced
- `file_handle`: File handle or reference that should be cleared before replacement (e.g., Arrow.Table, IOStream, or nothing)

# Implementation
- Clears the file_handle by setting it to nothing before attempting replacement
- On Windows: Multiple retry attempts, fallback to cmd del, final fallback to rename
- On Unix: Standard mv() call

This function handles common Windows file locking issues that occur with Arrow files.
"""
function safe_replace_file(temp_path::String, target_path::String, file_handle)
    # Clear the file handle before attempting replacement
    file_handle = nothing
    
    if Sys.iswindows()
        # Try to delete the existing file with retries (if it exists)
        if isfile(target_path)
            max_retries = 3  # Reduced retries since we're not relying on GC
            for i in 1:max_retries
                try
                    # Try to delete using Julia's rm with force flag
                    rm(target_path, force=true)
                    break
                catch
                    if i == max_retries
                        # If all retries failed, try Windows-specific deletion
                        try
                            run(`cmd /c del /f /q "$target_path"`)
                        catch
                            # If that also fails, rename the old file instead of deleting
                            backup_path = target_path * ".backup_" * string(time_ns())
                            try
                                mv(target_path, backup_path, force=true)
                                @warn "Could not delete $target_path, renamed to $backup_path"
                            catch
                                error("Unable to remove or rename existing file: $target_path")
                            end
                        end
                    else
                        # Brief wait before retrying (reduced from exponential backoff)
                        sleep(0.05)
                    end
                end
            end
        end
        
        # Move the temporary file to the final location
        try
            mv(temp_path, target_path, force=true)
        catch
            # If move fails, try copy and delete
            try
                cp(temp_path, target_path, force=true)
                rm(temp_path, force=true)
            catch
                error("Unable to replace file: $target_path")
            end
        end
    else
        # For Linux/MacOS, use the simple approach
        mv(temp_path, target_path, force=true)
    end
    return nothing
end

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
    
    # Replace original file using Windows-safe method
    safe_replace_file(temp_path, file_path(ref), tbl)
    
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
    
    # Replace original file using Windows-safe method
    safe_replace_file(temp_path, file_path(ref), tbl)
    
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

