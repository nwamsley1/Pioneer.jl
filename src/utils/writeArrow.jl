#=
function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname()
        Arrow.write(tpath, df)
        if isfile(fpath)
            run(`cmd /c del /f "$fpath"`)
        end
        mv(tpath, fpath, force = true)
    else     #If Linux/MacOS easy
        Arrow.write(fpath, df)
    end
    return nothing
end
=#
function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        # Create a unique temporary file
        tpath = tempname() * ".arrow"
        
        # Write to the temporary file
        Arrow.write(tpath, df)
        
        # Try to delete the existing file with retries
        if isfile(fpath)
            max_retries = 5
            for i in 1:max_retries
                try
                    # Force garbage collection to release any file handles
                    GC.gc()
                    
                    # Try to delete using Julia's rm with force flag
                    rm(fpath, force=true)
                    break
                catch e
                    if i == max_retries
                        # If all retries failed, try Windows-specific deletion
                        try
                            run(`cmd /c del /f /q "$fpath"`)
                        catch
                            # If that also fails, rename the old file instead of deleting
                            backup_path = fpath * ".backup_" * string(time_ns())
                            try
                                mv(fpath, backup_path, force=true)
                            catch
                                error("Unable to remove or rename existing file: $fpath")
                            end
                        end
                    else
                        # Wait a bit before retrying
                        sleep(0.1 * i)
                    end
                end
            end
        end
        
        # Move the temporary file to the final location
        try
            mv(tpath, fpath, force=true)
        catch e
            # If move fails, try copy and delete
            try
                cp(tpath, fpath, force=true)
                rm(tpath, force=true)
            catch
                error("Unable to write to file: $fpath")
            end
        end
    else
        # For Linux/MacOS, use the simple approach
        Arrow.write(fpath, df)
    end
    return nothing
end