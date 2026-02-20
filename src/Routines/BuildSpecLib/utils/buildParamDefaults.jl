# No module wrapper - these functions become part of the Pioneer module
# JSON is already imported in Pioneer.jl, no need to import again

"""
    get_build_default_parameters()

Load default parameters for library building from JSON file.

Returns:
- Dict containing default parameters for library building
"""
function get_build_default_parameters()
    # Always load full defaults
    json_path = asset_path("example_config", "defaultBuildLibParams.json")

    if !isfile(json_path)
        error("Default build parameters file not found: $json_path")
    end

    defaults = JSON.parsefile(json_path, dicttype=Dict{String,Any})

    # Remove user-specific sections that must be provided by user
    user_specific_keys = ["fasta_paths", "fasta_names", "library_path", "calibration_raw_file"]
    for key in user_specific_keys
        delete!(defaults, key)
    end

    return defaults
end

"""
    merge_with_build_defaults(user_params::Dict, defaults::Dict)

Recursively merge user parameters over default parameters.
User values always override defaults at any nesting level.
"""
function merge_with_build_defaults(user_params::Dict, defaults::Dict)
    result = deepcopy(defaults)
    
    function recursive_merge!(target::Dict, source::Dict)
        for (key, value) in source
            if haskey(target, key) && isa(target[key], Dict) && isa(value, Dict)
                # Both are dicts, merge recursively
                recursive_merge!(target[key], value)
            else
                # Override with user value
                target[key] = value
            end
        end
    end
    
    recursive_merge!(result, user_params)
    return result
end