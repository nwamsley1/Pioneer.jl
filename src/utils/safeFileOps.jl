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

# Safe file operations for cross-platform compatibility

"""
    safeRm(fpath::String, file_handle; force::Bool=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `file_handle`: File handle or reference that should be cleared before deletion (e.g., Arrow.Table, IOStream, or nothing)
- `force`: Force removal even if file is read-only

# Implementation
- Clears the file_handle by setting it to nothing before attempting deletion
- On Windows: Multiple retry attempts, fallback to cmd del, final fallback to rename
- On Unix: Standard rm() call

This function handles common Windows file locking issues that occur with Arrow files
and other binary formats that may have lingering file handles.
"""
function safeRm(fpath::String, file_handle; force::Bool=false)
    # Clear the file handle before attempting deletion
    file_handle = nothing
    
    if Sys.iswindows()
        # Return early if file doesn't exist
        if !isfile(fpath)
            return nothing
        end
        
        max_retries = 3  # Reduced since we're explicitly clearing handles
        for i in 1:max_retries
            try
                # Try Julia's rm with force flag
                rm(fpath, force=force)
                return nothing
            catch
                if i == max_retries
                    # If all retries failed, try Windows-specific deletion
                    try
                        # Convert to a Windows-style path for the cmd call
                        win_path = replace(abspath(fpath), "/" => "\\")
                        cmd_del = Cmd(["cmd", "/c", "del", "/f", "/q", win_path])
                        run(cmd_del)
                        return nothing
                    catch
                        # If that also fails, rename the old file instead of deleting
                        backup_path = fpath * ".backup_" * string(time_ns())
                        try
                            mv(fpath, backup_path, force=true)
                            @warn "Could not delete $fpath, renamed to $backup_path"
                            return nothing
                        catch
                            error("Unable to remove or rename file: $fpath")
                        end
                    end
                else
                    # Wait with exponential backoff before retrying
                    sleep(0.1 * i)
                end
            end
        end
    else
        # For Linux/MacOS, use the simple approach
        if isfile(fpath)
            rm(fpath, force=force)
        end
    end
    return nothing
end