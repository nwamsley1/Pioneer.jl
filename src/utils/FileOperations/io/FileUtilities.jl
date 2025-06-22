"""
General file utilities and operations.

Provides utilities for file path management, temporary files,
and general file operations that support the FileOperations
system.
"""

using DataFrames

#==========================================================
File Path and Management Utilities
==========================================================#

"""
    safe_join_files(left_ref::FileReference, right_ref::FileReference,
                   output_path::String, join_key::Symbol; 
                   join_type=:inner, batch_size=100_000)

Join two files with memory-efficient streaming.
Typically used for joining PSM and protein files.
"""
function safe_join_files(left_ref::FileReference, 
                        right_ref::FileReference,
                        output_path::String,
                        join_key::Symbol;
                        join_type::Symbol=:inner,
                        batch_size::Int=100_000)
    # Validate inputs using validation utilities
    validate_join_compatibility(left_ref, right_ref, join_key)
    
    # For efficient join, right side should be sorted by join key
    if !is_sorted_by(right_ref, join_key)
        @warn "Right side not sorted by join key, join may be slow"
    end
    
    # Load right side (protein groups) - typically smaller
    right_df = DataFrame(Arrow.Table(file_path(right_ref)))
    
    # Stream through left side (PSMs)
    left_tbl = Arrow.Table(file_path(left_ref))
    partitions = Tables.partitioner(left_tbl, batch_size)
    
    first_write = true
    for partition in partitions
        left_batch = DataFrame(partition)
        
        # Perform join
        joined_batch = if join_type == :inner
            innerjoin(left_batch, right_df, on=join_key)
        elseif join_type == :left
            leftjoin(left_batch, right_df, on=join_key)
        else
            error("Unsupported join type: $join_type")
        end
        
        if nrow(joined_batch) > 0
            if first_write
                writeArrow(output_path, joined_batch)
                first_write = false
            else
                Arrow.append(output_path, joined_batch)
            end
        end
    end
    
    # Create reference for output of same type as left input
    output_ref = create_reference(output_path, typeof(left_ref))
    
    return output_ref
end

"""
    ensure_directory_exists(file_path::String)

Ensure the directory for a file path exists, creating it if necessary.
"""
function ensure_directory_exists(file_path::String)
    dir = dirname(file_path)
    !isdir(dir) && mkpath(dir)
    return dir
end

"""
    create_temp_path(base_path::String, suffix::String="_tmp") -> String

Create a temporary file path based on the base path.
"""
function create_temp_path(base_path::String, suffix::String="_tmp")
    return base_path * suffix
end

"""
    cleanup_temp_files(refs::Vector{<:FileReference}; temp_suffix="_tmp")

Clean up temporary files associated with file references.
"""
function cleanup_temp_files(refs::Vector{<:FileReference}; temp_suffix="_tmp")
    for ref in refs
        temp_path = file_path(ref) * temp_suffix
        if isfile(temp_path)
            rm(temp_path, force=true)
        end
    end
end

"""
    safe_move_file(src::String, dst::String; backup::Bool=true)

Safely move a file from src to dst, optionally creating a backup.
"""
function safe_move_file(src::String, dst::String; backup::Bool=true)
    if backup && isfile(dst)
        backup_path = dst * ".backup"
        mv(dst, backup_path, force=true)
    end
    
    mv(src, dst, force=true)
    return dst
end

"""
    get_file_size_mb(file_path::String) -> Float64

Get file size in megabytes.
"""
function get_file_size_mb(file_path::String)
    if !isfile(file_path)
        return 0.0
    end
    
    return filesize(file_path) / (1024 * 1024)
end

"""
    estimate_processing_time(ref::FileReference, operations::Vector{String}) -> Float64

Estimate processing time for a file reference based on size and operations.
Returns estimated time in seconds.
"""
function estimate_processing_time(ref::FileReference, operations::Vector{String})
    # This is a rough estimate - in practice would be calibrated based on benchmarks
    base_time_per_mb = 0.1  # seconds per MB for basic operations
    
    file_size_mb = get_file_size_mb(file_path(ref))
    n_operations = length(operations)
    
    # Apply multipliers for different operation types
    complexity_multiplier = 1.0
    for op in operations
        if contains(op, "sort")
            complexity_multiplier *= 2.0  # Sorting is expensive
        elseif contains(op, "merge")
            complexity_multiplier *= 1.5  # Merging has overhead
        elseif contains(op, "filter")
            complexity_multiplier *= 0.5  # Filtering reduces data
        end
    end
    
    return file_size_mb * base_time_per_mb * n_operations * complexity_multiplier
end

# Export file utilities
export safe_join_files, ensure_directory_exists, create_temp_path, cleanup_temp_files,
       safe_move_file, get_file_size_mb, estimate_processing_time