# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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

