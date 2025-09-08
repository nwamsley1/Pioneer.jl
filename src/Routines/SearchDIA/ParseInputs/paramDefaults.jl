# No module wrapper - these functions become part of the Pioneer module
# JSON is already imported in Pioneer.jl, no need to import again

"""
    get_default_parameters(simplified::Bool = false)

Load default parameters from the appropriate JSON file.

Parameters:
- simplified: If true, loads simplified defaults, otherwise loads full defaults

Returns:
- Dict containing default parameters
"""
function get_default_parameters(simplified::Bool = false)
    # Determine which JSON file to load
    json_filename = simplified ? "defaultSearchParamsSimplified.json" : "defaultSearchParams.json"
    
    # Use asset_path to find the JSON file
    json_path = asset_path("example_config", json_filename)
    
    # Load and parse JSON
    if !isfile(json_path)
        error("Default parameters file not found: $json_path")
    end
    
    defaults = JSON.parsefile(json_path, dicttype=Dict{String,Any})
    
    # Remove user-specific sections that should not have defaults
    delete!(defaults, "paths")
    
    return defaults
end

"""
    merge_with_defaults(user_params::Dict, defaults::Dict)

Recursively merges user parameters over default parameters.
User values override defaults at any nesting level.
Missing sections or parameters are filled from defaults.
"""
function merge_with_defaults(user_params::Dict, defaults::Dict)
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