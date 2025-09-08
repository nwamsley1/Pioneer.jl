# Parameter Unification Plan for Pioneer.jl

## Executive Summary

This plan unifies default parameter handling across SearchDIA and BuildSpecLib modules by loading defaults from JSON files (the single source of truth) and consistently merging user parameters to override defaults.

## Core Principle

**JSON files are the single source of truth for all default parameters**
- User parameters always override defaults
- Consistent merge behavior across all modules
- Parameter validation happens after merging

## Implementation Plan

### Phase 1: Extract asset_path Function for Early Loading

#### 1.1 Keep asset_path in SearchDIA.jl but Load it Early

**File**: `src/Routines/SearchDIA.jl`

The `asset_path` function remains in SearchDIA.jl (lines 38-59) as it currently exists. This maintains compatibility and avoids creating a new module.

#### 1.2 Update importScripts.jl to Load SearchDIA.jl Early

**File**: `src/importScripts.jl`

Modify the loading order to load SearchDIA.jl very early (around line 100, before paramDefaults.jl):

```julia
# Line ~100: Load SearchDIA early to make asset_path available
safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA.jl"))

# Line ~101: Now load parameter defaults (can use asset_path)
safe_include!(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "paramDefaults.jl"))
safe_include!(joinpath(package_root, "src", "Routines","BuildSpecLib", "utils", "buildParamDefaults.jl"))
safe_include!(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))

# Remove the duplicate loading of SearchDIA.jl from line ~291
# safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA.jl"))  # DELETE THIS LINE
```

### Phase 2: Update SearchDIA Parameter Loading

#### 2.1 Modify paramDefaults.jl

**File**: `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl`

```julia
module ParamDefaults

using JSON
using ...SearchDIA: asset_path  # Import asset_path from SearchDIA module

export get_default_parameters, merge_with_defaults

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

Recursively merge user parameters over default parameters.
User values always override defaults at any nesting level.
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

end # module
```

#### 2.2 Update parseParams.jl

**File**: `src/Routines/SearchDIA/ParseInputs/parseParams.jl`

Update the function to determine simplified mode and pass it to get_default_parameters:

```julia
function parsePioneerToml(json_path::String; apply_defaults::Bool = true)
    # Read user JSON
    user_params = JSON.parsefile(json_path)
    
    # Apply defaults if requested
    if apply_defaults
        # Determine if we should use simplified defaults based on parameter presence
        is_simplified = !haskey(user_params, "parameter_tuning") && 
                       !haskey(user_params, "first_search") &&
                       !haskey(user_params, "quant_search")
        
        # Get appropriate defaults
        defaults = get_default_parameters(is_simplified)
        
        # Merge user params over defaults
        params = merge_with_defaults(user_params, defaults)
    else
        params = user_params
    end
    
    # Continue with existing code to convert to NamedTuple...
    # [rest of the function remains unchanged]
```

#### 2.3 Update paramsChecks.jl

**File**: `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`

Update to load defaults, merge, then validate:

```julia
function check_params(json_string::String)
    # Parse user parameters
    user_params = JSON.parse(json_string)
    
    # Determine if we should use simplified defaults
    is_simplified = !haskey(user_params, "parameter_tuning") && 
                   !haskey(user_params, "first_search") &&
                   !haskey(user_params, "quant_search")
    
    # Get appropriate defaults
    defaults = get_default_parameters(is_simplified)
    
    # Merge user params over defaults
    params = merge_with_defaults(user_params, defaults)
    
    # Now perform validation on the merged parameters
    # [Continue with all existing validation code but use 'params' instead of parsing json_string]
    
    # Helper function remains the same
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            error("Missing required parameter: $key")
        elseif !(dict[key] isa expected_type)
            error("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))")
        end
    end
    
    # [All the existing validation checks continue here, operating on 'params']
    # ...existing validation code...
    
    return params
end
```

### Phase 3: Update BuildSpecLib Parameter Loading

#### 3.1 Create buildParamDefaults.jl

**New File**: `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl`

```julia
module BuildParamDefaults

using JSON
using ...SearchDIA: asset_path  # Import asset_path from SearchDIA module

export get_build_default_parameters, merge_with_build_defaults

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
    
    # Use asset_path to find the JSON file
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

end # module
```

#### 3.2 Update check_params.jl

**File**: `src/Routines/BuildSpecLib/utils/check_params.jl`

The file already uses `asset_path` (line 176) which will be available from SearchDIA being loaded early.

Update the function to use the new parameter loading:

```julia
# Add at the top of the file after other imports
using ..BuildParamDefaults

function check_params_bsp(json_string::String)
    # Parse user parameters
    user_params = JSON.parse(json_string)
    
    # Determine if simplified based on presence of detailed parameters
    # Simplified version typically omits library_params details and nce_params
    is_simplified = !haskey(user_params, "nce_params") && 
                   (!haskey(user_params, "library_params") || 
                    length(user_params["library_params"]) < 5)
    
    # Get appropriate defaults
    defaults = get_build_default_parameters(is_simplified)
    
    # Merge user params over defaults
    params = merge_with_build_defaults(user_params, defaults)
    
    # IMPORTANT: REMOVE the inline defaults section (lines 72-96 in current file)
    # These defaults now come from defaultBuildLibParams.json:
    # library_defaults = Dict(
    #     "rt_bin_tol" => 1.0,
    #     "frag_bin_tol_ppm" => 10.0,
    #     "rank_to_score" => [8, 4, 4, 2, 2, 1, 1],
    #     "y_start_index" => 4,
    #     "b_start_index" => 3,
    #     "y_start" => 3,
    #     "b_start" => 2,
    #     "include_p_index" => false,
    #     "include_p" => false,
    #     "auto_detect_frag_bounds" => true,
    #     "max_frag_rank" => 255,
    #     "length_to_frag_count_multiple" => 2,
    #     "min_frag_intensity" => 0.0,
    #     "include_isotope" => false,
    #     "include_internal" => false,
    #     "include_immonium" => false,
    #     "include_neutral_diff" => true
    # )
    # DELETE the for loop that merges these hardcoded defaults
    
    # The nce_params default section (lines 42-48) should also be removed since
    # it will come from the JSON defaults
    # Remove:
    # if !haskey(params, "nce_params")
    #     params["nce_params"] = Dict(...)
    # end
    
    # Continue with validation using the merged params
    # Helper function remains the same
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end
    
    # Check required top-level parameters
    # [Continue with all existing validation checks, operating on 'params']
    required_params = ["fasta_digest_params", "library_params", "variable_mods", "fixed_mods", 
                      "isotope_mod_groups", "max_koina_requests", "max_koina_batch", 
                      "match_lib_build_batch", "fasta_paths", "fasta_names", "out_dir", 
                      "lib_name", "new_lib_name", "out_name", "predict_fragments", 
                      "include_contaminants"]
    
    for param in required_params
        check_param(params, param, param in ["predict_fragments", "include_contaminants"] ? Bool : Any)
    end
    
    # [Rest of validation continues unchanged, but operating on 'params']
    # Check fasta_digest_params...
    # Check nce_params...
    # Check library_params...
    # etc.
    
    return params
end
```

#### 3.3 Already Handled in Phase 1

The loading of buildParamDefaults.jl was already covered in Phase 1, step 1.2.

### Phase 4: No Changes Needed for asset_path

Since `asset_path` remains in SearchDIA.jl and SearchDIA is loaded early, all files that use `asset_path` will have access to it:

- `src/Routines/BuildSpecLib/utils/check_params.jl` - Already uses `asset_path`
- `src/Routines/GenerateParams.jl` - Already uses `asset_path`
- `src/Routines/BuildSpecLib.jl` - No direct use of `asset_path`

### Phase 5: Ensure JSON Files Have Complete Defaults

#### 5.1 Add Missing Defaults to defaultBuildLibParams.json

**File**: `assets/example_config/defaultBuildLibParams.json`

Ensure it contains all the defaults currently hardcoded in check_params.jl (lines 72-96):

```json
{
    "library_params": {
        "rt_bin_tol": 1.0,
        "frag_bin_tol_ppm": 10.0,
        "rank_to_score": [8, 4, 4, 2, 2, 1, 1],
        "y_start_index": 4,
        "b_start_index": 3,
        "y_start": 3,
        "b_start": 2,
        "include_p_index": false,
        "include_p": false,
        "auto_detect_frag_bounds": true,
        "max_frag_rank": 255,
        "length_to_frag_count_multiple": 2,
        "min_frag_intensity": 0.0,
        "include_isotope": false,
        "include_internal": false,
        "include_immonium": false,
        "include_neutral_diff": true,
        "frag_mz_min": 100.0,
        "frag_mz_max": 2000.0,
        "prec_mz_min": 400.0,
        "prec_mz_max": 1200.0,
        "max_frag_charge": 2,
        "instrument_type": "ORBIT",
        "prediction_model": "AlphaPept"
    },
    "nce_params": {
        "nce": 25.0,
        "default_charge": 2,
        "dynamic_nce": true
    },
    "fasta_digest_params": {
        "min_length": 7,
        "max_length": 30,
        "min_charge": 2,
        "max_charge": 4,
        "cleavage_regex": "[KR]",
        "missed_cleavages": 2,
        "max_var_mods": 2,
        "add_decoys": true,
        "entrapment_r": 0.01
    },
    "variable_mods": {
        "pattern": [],
        "mass": [],
        "name": []
    },
    "fixed_mods": {
        "pattern": ["C"],
        "mass": [57.021464],
        "name": ["Carbamidomethyl"]
    },
    "isotope_mod_groups": [],
    "max_koina_requests": 10,
    "max_koina_batch": 1000,
    "match_lib_build_batch": true,
    "predict_fragments": true,
    "include_contaminants": false
}
```

#### 5.2 Create defaultBuildLibParamsSimplified.json if it doesn't exist

**File**: `assets/example_config/defaultBuildLibParamsSimplified.json`

Create a simplified version with only essential parameters:

```json
{
    "library_params": {
        "instrument_type": "ORBIT",
        "prediction_model": "AlphaPept"
    },
    "fasta_digest_params": {
        "min_length": 7,
        "max_length": 30,
        "min_charge": 2,
        "max_charge": 4,
        "cleavage_regex": "[KR]",
        "missed_cleavages": 2,
        "max_var_mods": 2,
        "add_decoys": true,
        "entrapment_r": 0.01
    },
    "variable_mods": {
        "pattern": [],
        "mass": [],
        "name": []
    },
    "fixed_mods": {
        "pattern": ["C"],
        "mass": [57.021464],
        "name": ["Carbamidomethyl"]
    },
    "isotope_mod_groups": [],
    "max_koina_requests": 10,
    "max_koina_batch": 1000,
    "match_lib_build_batch": true,
    "predict_fragments": true,
    "include_contaminants": false
}
```

## Verification Steps

After implementation, verify:

1. **SearchDIA works with minimal config**: Create a config with only paths and run SearchDIA
2. **BuildSpecLib works with minimal config**: Create a config with only required fields and run BuildSpecLib
3. **User overrides work**: Verify that user-provided values override JSON defaults
4. **Both simplified and full defaults load**: Test with both parameter generation modes
5. **Error messages are clear**: Test with missing required parameters

## Key Implementation Notes

1. **Module dependencies**: SearchDIA.jl must be loaded early to provide `asset_path` function
2. **Loading order**: SearchDIA → paramDefaults → buildParamDefaults → parseParams
3. **JSON structure**: Ensure JSON files have all parameters currently hardcoded
4. **User-specific fields**: Never include paths or user-specific fields in defaults
5. **Backwards compatibility**: Existing user configs should continue to work
6. **Error handling**: Clear error messages when JSON files are missing
7. **No new modules**: Solution uses existing SearchDIA module for `asset_path`

## Success Criteria

- [x] Single source of truth (JSON files) for all defaults
- [x] No code duplication for path resolution (uses asset_path)
- [x] Consistent merge behavior across modules
- [x] User parameters always override defaults
- [x] Validation happens after merging
- [x] Clear separation of concerns (defaults, merging, validation)