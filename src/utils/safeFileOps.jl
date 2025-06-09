# Safe file operations for cross-platform compatibility

"""
    safeRm(fpath::String; force::Bool=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `force`: Force removal even if file is read-only

# Implementation
- On Windows: Multiple retry attempts with GC, fallback to cmd del, final fallback to rename
- On Unix: Standard rm() call

This function handles common Windows file locking issues that occur with Arrow files
and other binary formats that may have lingering file handles.
"""
function safeRm(fpath::String; force::Bool=false)
    if Sys.iswindows()
        # Return early if file doesn't exist
        if !isfile(fpath)
            return nothing
        end
        
        max_retries = 5
        for i in 1:max_retries
            try
                # Force garbage collection to release any file handles
                GC.gc()
                
                # Try Julia's rm with force flag
                rm(fpath, force=force)
                return nothing
            catch e
                if i == max_retries
                    # If all retries failed, try Windows-specific deletion
                    try
                        run(`cmd /c del /f /q "$fpath"`)
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