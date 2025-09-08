# Fix Module Warnings in Parameter Loading System

## Executive Summary

**Current State**: Functions are wrapped in modules (ParamDefaults, BuildParamDefaults) causing compilation errors  
**Solution**: Remove module wrappers, make functions part of Pioneer module directly  
**Result**: Clean compilation, no warnings, tests will pass, same functionality

## Problem Summary

The current implementation creates nested modules (ParamDefaults and BuildParamDefaults) that are causing precompilation warnings:
- `invalid using path: "SearchDIA" does not name a module`
- `WARNING: detected unclosed module: Pioneer.BuildParamDefaults`
- `WARNING: detected unclosed module: Pioneer.ParamDefaults`

## Root Cause

1. `SearchDIA` is not a module - it's just a file containing functions in the Pioneer module
2. The nested modules are trying to import `asset_path` using an invalid module path
3. The module errors prevent proper module closure during compilation

## Solution: Convert Modules to Plain Functions

Remove the unnecessary module wrappers and make the functions part of the Pioneer module directly.

## Detailed Implementation Changes

### 1. Update `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl`

**REMOVE the module wrapper and change from:**
```julia
module ParamDefaults

using JSON
using ...SearchDIA: asset_path

export get_default_parameters, merge_with_defaults

# ... functions ...

end # module
```

**TO plain functions:**
```julia
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
    
    # Use asset_path to find the JSON file (asset_path is already in Pioneer module)
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
```

### 2. Update `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl`

**REMOVE the module wrapper and change from:**
```julia
module BuildParamDefaults

using JSON
using ...SearchDIA: asset_path

export get_build_default_parameters, merge_with_build_defaults

# ... functions ...

end # module
```

**TO plain functions:**
```julia
# No module wrapper - these functions become part of the Pioneer module
# JSON is already imported in Pioneer.jl, no need to import again

"""
    get_build_default_parameters(simplified::Bool = false)

Load default parameters for library building from JSON file.

Parameters:
- simplified: If true, loads simplified defaults, otherwise loads full defaults

Returns:
- Dict containing default parameters for library building
"""
function get_build_default_parameters(simplified::Bool = false)
    # Determine which JSON file to load
    json_filename = simplified ? "defaultBuildLibParamsSimplified.json" : "defaultBuildLibParams.json"
    
    # Use asset_path to find the JSON file (asset_path is already in Pioneer module)
    json_path = asset_path("example_config", json_filename)
    
    # Load and parse JSON
    if !isfile(json_path)
        error("Default build parameters file not found: $json_path")
    end
    
    defaults = JSON.parsefile(json_path, dicttype=Dict{String,Any})
    
    # Remove user-specific sections that must be provided by user
    user_specific_keys = ["fasta_paths", "fasta_names", "out_dir", "lib_name", "new_lib_name", "out_name"]
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
```

### 3. Update `src/Routines/SearchDIA/ParseInputs/parseParams.jl`

**REMOVE the module import on line 19:**
```julia
# ParamDefaults module is loaded in importScripts.jl
using .ParamDefaults
```

**CHANGE TO a comment:**
```julia
# Parameter default functions are loaded from paramDefaults.jl
```

**No other changes needed** - the functions `get_default_parameters` and `merge_with_defaults` will be directly available in the Pioneer module.

### 4. Update `src/Routines/BuildSpecLib/utils/check_params.jl`

**REMOVE the module import on line 18:**
```julia
using ..BuildParamDefaults
```

**No replacement needed** - just delete this line. The functions `get_build_default_parameters` and `merge_with_build_defaults` will be directly available in the Pioneer module.

## Key Dependencies Already Available

Since these functions will be part of the Pioneer module, they have access to:
- `JSON` - Already imported in Pioneer.jl (`using JSON, JLD2`)
- `asset_path` - Already defined in SearchDIA.jl which is loaded into Pioneer module

No additional imports are needed!

## Benefits of This Approach

1. **Eliminates module loading errors** - No more invalid module paths
2. **Simplifies code structure** - Fewer nested modules to manage
3. **Direct access to dependencies** - JSON and asset_path are already in the Pioneer module
4. **Maintains all functionality** - Same functions, just without module wrappers
5. **Cleaner precompilation** - No unclosed module warnings
6. **Reduces redundant imports** - No need to re-import JSON

## Impact on Tests

The tests are currently failing because `get_build_default_parameters` is being called but not found. This is expected since we haven't implemented the changes yet. After implementing the changes:

1. **The tests will work automatically** - The functions will be available in the Pioneer module
2. **No test updates needed** - Tests call `Pioneer.BuildSpecLib(params_file)` which internally uses `check_params_bsp`
3. **Function availability** - Once the module wrappers are removed, the functions become part of Pioneer and are accessible

## Testing After Changes

After making these changes:
1. **Restart Julia** to clear any cached module state
   ```bash
   # Exit Julia completely and restart
   ```

2. **Test that Pioneer.jl loads without warnings:**
   ```julia
   using Pkg
   Pkg.activate(".")
   using Pioneer
   ```

3. **Verify parameter loading functions are available:**
   ```julia
   # Test SearchDIA defaults (should work)
   defaults = Pioneer.get_default_parameters()
   println("SearchDIA defaults loaded: ", length(defaults), " sections")
   
   # Test BuildSpecLib defaults (should work)
   build_defaults = Pioneer.get_build_default_parameters()
   println("BuildSpecLib defaults loaded: ", length(build_defaults), " sections")
   ```

4. **Run the BuildSpecLib tests to confirm they pass:**
   ```julia
   using Pkg
   Pkg.test("Pioneer", test_args=["Routines/BuildSpecLib/test_build_spec_lib.jl"])
   ```

5. **Build docs again to ensure no warnings:**
   ```bash
   cd docs
   julia --project=. make.jl
   ```

## Alternative Approach (Not Recommended)

If we wanted to keep the modules, we would need to:
1. Fix the import path by removing the `SearchDIA` reference
2. Ensure asset_path is accessible via proper module hierarchy
3. Add proper module end statements

However, this adds unnecessary complexity for functions that don't need to be in separate modules.