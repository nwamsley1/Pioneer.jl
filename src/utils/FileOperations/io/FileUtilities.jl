"""
General file utilities and operations.

Provides essential utilities for file path management and basic file operations
that support the FileOperations system.
"""

using DataFrames

#==========================================================
Essential File Utilities
==========================================================#

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
    get_file_size_mb(file_path::String) -> Float64

Get file size in megabytes.
"""
function get_file_size_mb(file_path::String)
    if !isfile(file_path)
        return 0.0
    end
    
    return filesize(file_path) / (1024 * 1024)
end
