# This file initializes the Julia depot path for Windows installations
# It must run before any package operations

if Sys.iswindows() && get(ENV, "JULIA_DEPOT_PATH", "") == ""
    # Check if we're running from Program Files (installed location)
    exe_path = unsafe_string(Base.JLOptions().julia_bin)
    if occursin("Program Files", exe_path)
        # Set depot to user's local app data
        local_appdata = get(ENV, "LOCALAPPDATA", "")
        if !isempty(local_appdata)
            depot_path = joinpath(local_appdata, "Pioneer", "julia")
            
            # Set environment variable
            ENV["JULIA_DEPOT_PATH"] = depot_path
            
            # Create directory if it doesn't exist
            mkpath(depot_path)
            
            # Reset Julia's depot path
            empty!(DEPOT_PATH)
            push!(DEPOT_PATH, depot_path)
            
            # Also set scratch directory
            ENV["JULIA_SCRATCH"] = joinpath(depot_path, "scratchspaces")
            mkpath(ENV["JULIA_SCRATCH"])
        end
    end
end