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

# Export essential utilities
export ensure_directory_exists, get_file_size_mb