# BuildSpecLib Parameter Simplification Plan

## Executive Summary

This document provides a detailed implementation plan for simplifying BuildSpecLib parameters based on user feedback. The plan addresses specific concerns about path redundancy, parameter complexity, and the need for SearchDIA-style parameter validation and persistence.

---

## Table of Contents

1. [Output Path Simplification](#1-output-path-simplification)
2. [Parameter Restructuring and Reordering](#2-parameter-restructuring-and-reordering)
3. [Parameters to Remove from Simplified Version](#3-parameters-to-remove-from-simplified-version)
4. [Parameter Clarifications](#4-parameter-clarifications)
5. [SearchDIA-Style Parameter Handling](#5-searchdia-style-parameter-handling)
6. [Implementation Checklist](#6-implementation-checklist)

---

## 1. Output Path Simplification

### Current Problem

Currently, BuildSpecLib requires 4 separate parameters for outputs:

```json
{
  "out_dir": "/path/to/output/directory",
  "lib_name": "/path/to/output/library_name",
  "new_lib_name": "/path/to/output/library_name",
  "out_name": "library_name.tsv"
}
```

**Analysis of current usage:**
- `out_dir`: Base output directory
- `lib_name`: Used in code as: `joinpath(out_dir, lib_name * ".poin")`
- `new_lib_name`: **NEVER USED** in BuildSpecLib.jl (only appears in check_params validation)
- `out_name`: **NEVER USED** in BuildSpecLib.jl

### Proposed Solution

**Single parameter:** `library_path`

```json
{
  "library_path": "/path/to/output/my_library"
}
```

**The code intelligently handles the `.poin` extension:**

**Case 1: User omits `.poin` extension** (most common)
```json
"library_path": "/path/to/output/my_library"
```
- Library directory: `/path/to/output/my_library.poin/`
- Library name: `my_library`

**Case 2: User includes `.poin` extension**
```json
"library_path": "/path/to/output/my_library.poin"
```
- Library directory: `/path/to/output/my_library.poin/`
- Library name: `my_library` (`.poin` stripped from name)

**Result:** No matter which format the user provides, the behavior is consistent and correct.

### Code Changes

#### Current Code (BuildSpecLib.jl:72-76)

```julia
# Create output directories
lib_out_dir = params["out_dir"]
mkpath(lib_out_dir)

# Library directory (.poin extension)
lib_dir = joinpath(lib_out_dir, params["lib_name"] * ".poin")
mkpath(lib_dir)
```

#### New Code

```julia
# Parse library path and ensure .poin extension
library_path = params["library_path"]

# Add .poin extension if not already present
if !endswith(library_path, ".poin")
    lib_dir = library_path * ".poin"
    lib_name = basename(library_path)
else
    lib_dir = library_path
    lib_name = basename(library_path[1:end-5])  # Remove .poin for name extraction
end

mkpath(lib_dir)
```

#### Parameter Validation (check_params.jl)

**Remove from required_params list:**
```julia
# OLD:
required_params = [
    "fasta_digest_params", "library_params", "variable_mods", "fixed_mods",
    "isotope_mod_groups", "max_koina_requests", "max_koina_batch",
    "match_lib_build_batch", "fasta_paths", "fasta_names",
    "out_dir", "lib_name", "new_lib_name", "out_name",  # ← Remove these
    "predict_fragments", "include_contaminants"
]

# NEW:
required_params = [
    "fasta_digest_params", "library_params", "variable_mods", "fixed_mods",
    "isotope_mod_groups", "max_koina_requests", "max_koina_batch",
    "match_lib_build_batch", "fasta_paths", "fasta_names",
    "library_path",  # ← Single parameter
    "predict_fragments", "include_contaminants"
]
```

**Add validation:**
```julia
# Validate library_path
if !haskey(params, "library_path")
    throw(InvalidParametersError("Missing required parameter: library_path", params))
end

# Expand home directory
params["library_path"] = expanduser(params["library_path"])

# Extract derived values for internal use
# Handle .poin extension - add if not present, don't duplicate if already there
if !endswith(params["library_path"], ".poin")
    params["_lib_dir"] = params["library_path"] * ".poin"
    params["_lib_name"] = basename(params["library_path"])
else
    params["_lib_dir"] = params["library_path"]
    params["_lib_name"] = basename(params["library_path"][1:end-5])  # Strip .poin
end
```

---

## 2. Parameter Restructuring and Reordering

### 2.1 Move calibration_raw_file to Top Level

**Current structure:**
```json
{
    "library_params": {
        "calibration_raw_file": "/path/to/file.arrow",
        "auto_detect_frag_bounds": true,
        "frag_mz_min": 150.0,
        ...
    }
}
```

**New structure:**
```json
{
    "library_path": "/path/to/output/my_library",
    "calibration_raw_file": "/path/to/calibration/file.arrow",
    "fasta_paths": [...],
    ...
    "library_params": {
        "auto_detect_frag_bounds": true,
        "frag_mz_min": 150.0,
        ...
    }
}
```

**Rationale:**
- `calibration_raw_file` is one of the most important user-provided inputs
- Used for automatic m/z bounds detection (critical feature)
- Deserves top-level prominence alongside library path and FASTA files
- Makes it clear that this is a required input, not a nested library setting

### 2.2 Parameter Ordering in Simplified Version

**New logical order:**
1. **Library output path** - where to save the library
2. **Calibration file** - for automatic bounds detection
3. **FASTA inputs** - source proteins
4. **FASTA regex** - header parsing
5. **Library m/z bounds** - search ranges
6. **Digestion settings** - how to digest proteins
7. **Modifications** - PTMs and labels

This order follows the user's mental model: "Where am I saving? What file do I need? What proteins? How to process them?"

---

## 3. Parameters to Remove from Simplified Version

### 3.1 Instrument and Prediction Model

**Remove:**
```json
"instrument_type": "NONE",
"prediction_model": "altimeter"
```

**Rationale:** Most users use Altimeter (default predictor). Advanced users who need instrument-specific models or alternative predictors (Prosit, MS2PIP) can use the full parameter file.

**Implementation:**

```julia
# In buildParamDefaults.jl, simplified defaults:
function get_simplified_library_params()
    return Dict(
        # ... other params ...
        "instrument_type" => "NONE",      # Hidden default
        "prediction_model" => "altimeter"  # Hidden default
    )
end
```

**JSON before (simplified):**
```json
"library_params": {
    "auto_detect_frag_bounds": true,
    "calibration_raw_file": "/path/to/file.arrow",
    "frag_mz_min": 150.0,
    "frag_mz_max": 2020.0,
    "prec_mz_min": 390.0,
    "prec_mz_max": 1010.0,
    "instrument_type": "NONE",        # ← Remove
    "prediction_model": "altimeter"   # ← Remove
}
```

**JSON after (simplified):**
```json
"library_params": {
    "auto_detect_frag_bounds": true,
    "calibration_raw_file": "/path/to/file.arrow",
    "frag_mz_min": 150.0,
    "frag_mz_max": 2020.0,
    "prec_mz_min": 390.0,
    "prec_mz_max": 1010.0
}
```

### 2.2 Decoy and Entrapment Parameters

**Remove from fasta_digest_params:**
```json
"add_decoys": true,
"decoy_method": "shuffle",
"entrapment_method": "shuffle",
"entrapment_r": 0
```

**Rationale:**
- `add_decoys`: Almost always `true` for FDR control - rare to build library without decoys
- `decoy_method`: "shuffle" is standard (reverse is alternative for specific cases)
- `entrapment_method`: "shuffle" is standard
- `entrapment_r`: 0 means no entrapment groups (most common case)
- Users doing advanced entrapment analysis will use the full parameter file

**Implementation:**
```julia
# In buildParamDefaults.jl:
function get_simplified_fasta_digest_params()
    return Dict(
        "min_length" => 7,
        "max_length" => 40,
        "min_charge" => 2,
        "max_charge" => 3,
        "cleavage_regex" => "[KR][^_|$]",
        "missed_cleavages" => 1,
        "max_var_mods" => 1,
        # Hidden defaults:
        "add_decoys" => true,
        "decoy_method" => "shuffle",
        "entrapment_method" => "shuffle",
        "entrapment_r" => 0
    )
end
```

**JSON before (simplified):**
```json
"fasta_digest_params": {
    "min_length": 7,
    "max_length": 40,
    "min_charge": 2,
    "max_charge": 3,
    "cleavage_regex": "[KR][^_|$]",
    "missed_cleavages": 1,
    "max_var_mods": 1,
    "add_decoys": true,           # ← Remove
    "entrapment_r": 0,            # ← Remove
    "decoy_method": "shuffle",    # ← Remove
    "entrapment_method": "shuffle" # ← Remove
}
```

**JSON after (simplified):**
```json
"fasta_digest_params": {
    "min_length": 7,
    "max_length": 40,
    "min_charge": 2,
    "max_charge": 3,
    "cleavage_regex": "[KR][^_|$]",
    "missed_cleavages": 1,
    "max_var_mods": 1
}
```

### 2.3 NCE Parameters

**Remove entire section:**
```json
"nce_params": {
    "nce": 26.0,
    "default_charge": 2,
    "dynamic_nce": true
}
```

**Rationale:** Collision energy is instrument-specific and usually extracted from raw files. Manual NCE specification is rare.

**Implementation:**
```julia
# Hidden in defaults:
function get_simplified_nce_params()
    return Dict(
        "nce" => 26.0,
        "default_charge" => 2,
        "dynamic_nce" => true
    )
end
```

### 2.4 Processing Parameters

**Remove all processing parameters from simplified version:**

```json
"isotope_mod_groups": [],
"max_koina_requests": 24,
"max_koina_batch": 1000,
"match_lib_build_batch": 100000,
"predict_fragments": true,
"include_contaminants": true
```

**Rationale:**

- **`isotope_mod_groups`**: Only needed for plexDIA/TMT/SILAC experiments (advanced use case)
- **`predict_fragments`**: Almost always `true` - users building libraries want fragment prediction
- **`include_contaminants`**: Almost always `true` - rare to exclude common contaminants
- **`max_koina_*` and `match_lib_build_batch`**: Performance tuning parameters that work well with defaults

**Users who need these settings can use the full parameter file.**

**Implementation:**

```julia
# In buildParamDefaults.jl:
function get_simplified_processing_params()
    return Dict(
        # Hidden defaults:
        "isotope_mod_groups" => [],
        "max_koina_requests" => 24,
        "max_koina_batch" => 1000,
        "match_lib_build_batch" => 100000,
        "predict_fragments" => true,
        "include_contaminants" => true
    )
end
```

**Parameter Explanations (for full version documentation):**

- **`isotope_mod_groups`**: Used for plexDIA, TMT, SILAC, or any isotopic labeling experiment. Empty `[]` for label-free.
- **`predict_fragments`**:
  - `true` (default): Predict fragment intensities using Koina API
  - `false`: Skip prediction (use when library already has fragments)
- **`include_contaminants`**: Auto-adds common contaminant proteins (cRAP database)
- **`max_koina_requests`**: Parallel API requests (increase for faster builds, max ~50)
- **`max_koina_batch`**: Peptides per request (1000 is usually optimal)
- **`match_lib_build_batch`**: Memory management (increase if RAM available)

---

## 4. Parameter Clarifications

### 4.1 Rename "Protein modifications" → "Modifications"

**Current:**
```json
"comment_modifications": "=== Protein modifications ==="
```

**New:**
```json
"comment_modifications": "=== Modifications (variable and fixed) ==="
```

**Rationale:** Clarifies that this includes:
- PTMs (post-translational modifications)
- Peptide N/C-terminal modifications (e.g., TMT labels)
- Fixed modifications (e.g., carbamidomethyl on C)

### 4.2 Add UniProt Regex Comment

**Current:**
```json
"comment_paths": "=== File paths - Edit these first ===",
"fasta_paths": [...],
"fasta_names": [...],
"fasta_header_regex_accessions": [...]
```

**New:**
```json
"comment_paths": "=== File paths - Edit these first ===",
"fasta_paths": [...],
"fasta_names": [...],

"comment_fasta_regex": "=== FASTA header parsing (defaults work for UniProt format) ===",
"fasta_header_regex_accessions": ["^\\w+\\|(\\w+(?:-\\d+)?)\\|"],
"fasta_header_regex_genes": [" GN=(\\S+)"],
"fasta_header_regex_proteins": ["^\\w+\\|\\w+?:-\\d+?\\|[^ ]+ (.*?) [^ ]+="],
"fasta_header_regex_organisms": [" OS=([^ ]+.*?) [^ ]+="
]
```

---

## 5. SearchDIA-Style Parameter Handling

### Current Problems

1. Only user-provided JSON is saved (line 83: `write(params_out_path, params_string)`)
2. No merged user+defaults output
3. No structured parameter object
4. Limited bounds checking

### Proposed Solution: Mirror SearchDIA Architecture

#### 5.1 Create BuildSpecLibParameters Struct

**New file:** `src/Routines/BuildSpecLib/utils/parameter_types.jl`

```julia
"""
Structured parameter container for BuildSpecLib, similar to PioneerParameters for SearchDIA.
"""
struct BuildSpecLibParameters
    paths::NamedTuple
    fasta_digest::NamedTuple
    modifications::NamedTuple
    library::NamedTuple
    processing::NamedTuple
end

"""
Convert BuildSpecLibParameters back to Dict for JSON serialization.
Ensures complete merged configuration (user + defaults) is written to output.
"""
function params_to_dict(params::BuildSpecLibParameters)
    function namedtuple_to_dict(nt::NamedTuple)
        result = Dict{String, Any}()
        for (k, v) in pairs(nt)
            key_str = string(k)
            result[key_str] = if v isa NamedTuple
                namedtuple_to_dict(v)
            elseif v isa Vector{<:NamedTuple}
                [namedtuple_to_dict(x) for x in v]
            else
                v
            end
        end
        return result
    end

    return Dict{String, Any}(
        "paths" => namedtuple_to_dict(params.paths),
        "fasta_digest_params" => namedtuple_to_dict(params.fasta_digest),
        "modifications" => namedtuple_to_dict(params.modifications),
        "library_params" => namedtuple_to_dict(params.library),
        "processing" => namedtuple_to_dict(params.processing)
    )
end
```

#### 5.2 Enhanced Parameter Validation

**New file:** `src/Routines/BuildSpecLib/utils/parameter_bounds.jl`

```julia
"""
Validate parameter bounds and types, similar to SearchDIA parameter checks.
"""
function validate_parameter_bounds!(params::Dict)
    # Validate fasta_digest_params
    digest = params["fasta_digest_params"]

    # Length bounds
    @assert digest["min_length"] >= 5 "min_length must be >= 5"
    @assert digest["max_length"] <= 50 "max_length must be <= 50 (too long peptides rarely observed)"
    @assert digest["min_length"] < digest["max_length"] "min_length must be < max_length"

    # Charge bounds
    @assert digest["min_charge"] >= 1 "min_charge must be >= 1"
    @assert digest["max_charge"] <= 6 "max_charge must be <= 6 (higher charges rare in bottom-up proteomics)"
    @assert digest["min_charge"] <= digest["max_charge"] "min_charge must be <= max_charge"

    # Missed cleavages
    @assert digest["missed_cleavages"] >= 0 "missed_cleavages must be >= 0"
    @assert digest["missed_cleavages"] <= 3 "missed_cleavages must be <= 3 (higher values create huge libraries)"

    # Variable mods
    @assert digest["max_var_mods"] >= 0 "max_var_mods must be >= 0"
    @assert digest["max_var_mods"] <= 5 "max_var_mods must be <= 5 (higher values create huge libraries)"

    # Entrapment
    @assert 0 <= digest["entrapment_r"] <= 1 "entrapment_r must be in [0, 1]"

    # Validate library_params
    lib = params["library_params"]

    # Tolerances
    @assert lib["rt_bin_tol"] > 0 "rt_bin_tol must be > 0"
    @assert lib["frag_bin_tol_ppm"] > 0 "frag_bin_tol_ppm must be > 0"

    # m/z bounds
    @assert 50 <= lib["frag_mz_min"] <= 2000 "frag_mz_min out of reasonable range [50, 2000]"
    @assert 100 <= lib["frag_mz_max"] <= 4000 "frag_mz_max out of reasonable range [100, 4000]"
    @assert lib["frag_mz_min"] < lib["frag_mz_max"] "frag_mz_min must be < frag_mz_max"

    @assert 200 <= lib["prec_mz_min"] <= 2000 "prec_mz_min out of reasonable range [200, 2000]"
    @assert 400 <= lib["prec_mz_max"] <= 4000 "prec_mz_max out of reasonable range [400, 4000]"
    @assert lib["prec_mz_min"] < lib["prec_mz_max"] "prec_mz_min must be < prec_mz_max"

    # Fragment selection
    @assert lib["max_frag_charge"] >= 1 "max_frag_charge must be >= 1"
    @assert lib["max_frag_charge"] <= 4 "max_frag_charge must be <= 4 (higher charges rare for fragments)"

    @assert 0 <= lib["min_frag_intensity"] < 1 "min_frag_intensity must be in [0, 1)"

    # Validate processing params
    @assert params["max_koina_requests"] >= 1 "max_koina_requests must be >= 1"
    @assert params["max_koina_requests"] <= 100 "max_koina_requests must be <= 100 (API rate limiting)"

    @assert params["max_koina_batch"] >= 100 "max_koina_batch must be >= 100"
    @assert params["max_koina_batch"] <= 10000 "max_koina_batch must be <= 10000 (API limits)"

    @assert params["match_lib_build_batch"] >= 10000 "match_lib_build_batch too small, will be very slow"

    return nothing
end
```

#### 5.3 Updated check_params Function

**Modified:** `src/Routines/BuildSpecLib/utils/check_params.jl`

```julia
function check_params_bsp(json_string::String)
    # Parse user parameters
    user_params = JSON.parse(json_string)

    # Determine if simplified based on presence of detailed parameters
    is_simplified = !haskey(user_params, "nce_params") &&
                   (!haskey(user_params, "library_params") ||
                    length(user_params["library_params"]) < 10)

    # Get appropriate defaults
    defaults = get_build_default_parameters(is_simplified)

    # Merge user params over defaults (deep merge)
    params = merge_with_build_defaults(user_params, defaults)

    # Validate parameter bounds
    validate_parameter_bounds!(params)

    # ... (existing validation code) ...

    # Return complete merged parameters
    return params
end
```

#### 5.4 Write Complete Parameters to config.json

**Modified:** `BuildSpecLib.jl`

```julia
# OLD (line 83):
write(params_out_path, params_string)

# NEW:
# Write the COMPLETE merged parameters (user + defaults) to config.json
params_dict = params  # Already a Dict from check_params_bsp
params_json = JSON.json(params_dict, 2)  # Pretty-print with 2-space indent
write(params_out_path, params_json)

dual_println("Complete parameter configuration saved to: ", params_out_path)
dual_println("  (includes all defaults merged with user settings)")
```

---

## 6. Implementation Checklist

### Phase 1: Output Path Simplification
- [ ] Update `BuildSpecLib.jl` to use single `library_path` parameter
- [ ] Modify `check_params.jl` to remove `out_dir`, `lib_name`, `new_lib_name`, `out_name`
- [ ] Add validation for `library_path`
- [ ] Update example JSON files
- [ ] Update documentation

### Phase 2: Simplified Parameters
- [ ] Move parameters to defaults in `buildParamDefaults.jl`:
  - [ ] `instrument_type`, `prediction_model`
  - [ ] `decoy_method`, `entrapment_method`, `entrapment_r`
  - [ ] `nce_params` section
- [ ] Update `defaultBuildLibParamsSimplified.json` to remove hidden params
- [ ] Keep `defaultBuildLibParams.json` with all options for advanced users
- [ ] Add better comments explaining `isotope_mod_groups` and `predict_fragments`

### Phase 3: Parameter Naming and Comments
- [ ] Rename "Protein modifications" → "Modifications"
- [ ] Add UniProt regex comment
- [ ] Improve processing parameter comments

### Phase 4: SearchDIA-Style Parameter Handling
- [ ] Create `parameter_types.jl` with `BuildSpecLibParameters` struct
- [ ] Create `parameter_bounds.jl` with validation functions
- [ ] Update `check_params_bsp` to return merged params
- [ ] Modify `BuildSpecLib.jl` to write complete merged params to config.json
- [ ] Add parameter bounds checking

### Phase 5: Testing
- [ ] Test simplified config with minimal parameters
- [ ] Test full config with all parameters
- [ ] Verify parameter bounds checking catches invalid values
- [ ] Verify config.json contains complete merged parameters
- [ ] Update all test config files in `data/test_build_spec_lib/`

### Phase 6: Documentation
- [ ] Update user guide with new parameter structure
- [ ] Add migration guide for old → new parameter files
- [ ] Document default values for hidden parameters
- [ ] Add examples for common use cases

---

## Example: Before and After

### Before (Current Simplified Config)

```json
{
    "fasta_paths": ["/path/to/file.fasta"],
    "fasta_names": ["HUMAN"],
    "fasta_header_regex_accessions": ["^\\w+\\|(\\w+(?:-\\d+)?)\\|"],
    "fasta_header_regex_genes": [" GN=(\\S+)"],
    "fasta_header_regex_proteins": ["^\\w+\\|\\w+?:-\\d+?\\|[^ ]+ (.*?) [^ ]+="],
    "fasta_header_regex_organisms": [" OS=([^ ]+.*?) [^ ]+="],
    "out_dir": "/path/to/output/directory",
    "lib_name": "library_name",
    "new_lib_name": "library_name",
    "out_name": "library_name.tsv",

    "library_params": {
        "auto_detect_frag_bounds": true,
        "calibration_raw_file": "/path/to/file.arrow",
        "frag_mz_min": 150.0,
        "frag_mz_max": 2020.0,
        "prec_mz_min": 390.0,
        "prec_mz_max": 1010.0,
        "instrument_type": "NONE",
        "prediction_model": "altimeter"
    },

    "fasta_digest_params": {
        "min_length": 7,
        "max_length": 40,
        "min_charge": 2,
        "max_charge": 3,
        "cleavage_regex": "[KR][^_|$]",
        "missed_cleavages": 1,
        "max_var_mods": 1,
        "add_decoys": true,
        "entrapment_r": 0,
        "decoy_method": "shuffle",
        "entrapment_method": "shuffle"
    },

    "variable_mods": {
        "pattern": ["M"],
        "mass": [15.99491],
        "name": ["Unimod:35"]
    },
    "fixed_mods": {
        "pattern": ["C"],
        "mass": [57.021464],
        "name": ["Unimod:4"]
    },

    "nce_params": {
        "nce": 26.0,
        "default_charge": 2,
        "dynamic_nce": true
    },

    "isotope_mod_groups": [],
    "max_koina_requests": 24,
    "max_koina_batch": 1000,
    "match_lib_build_batch": 100000,
    "include_contaminants": true,
    "predict_fragments": true
}
```

**Total: ~40 user-visible parameters**

### After (New Simplified Config)

```json
{
    "comment_output": "=== Output library path ===",
    "library_path": "/path/to/output/my_library",

    "comment_calibration": "=== Calibration file for automatic m/z bounds detection ===\n    Convert a representative RAW file to Arrow format, then provide path here.\n    Used to detect fragment and precursor m/z ranges for your instrument.",
    "calibration_raw_file": "/path/to/calibration/file.arrow",

    "comment_fasta": "=== FASTA files ===",
    "fasta_paths": ["/path/to/file.fasta"],
    "fasta_names": ["HUMAN"],

    "comment_fasta_regex": "=== FASTA header parsing (defaults work for UniProt format) ===",
    "fasta_header_regex_accessions": ["^\\w+\\|(\\w+(?:-\\d+)?)\\|"],
    "fasta_header_regex_genes": [" GN=(\\S+)"],
    "fasta_header_regex_proteins": ["^\\w+\\|\\w+?:-\\d+?\\|[^ ]+ (.*?) [^ ]+="],
    "fasta_header_regex_organisms": [" OS=([^ ]+.*?) [^ ]+="],

    "comment_library": "=== Library m/z bounds ===\n    Note: If auto_detect_frag_bounds=true (default), these bounds are automatically\n    detected from calibration_raw_file. Set to false to use manual bounds below.",
    "library_params": {
        "auto_detect_frag_bounds": true,
        "frag_mz_min": 150.0,
        "frag_mz_max": 2020.0,
        "prec_mz_min": 390.0,
        "prec_mz_max": 1010.0
    },

    "comment_digestion": "=== Protein digestion settings ===",
    "fasta_digest_params": {
        "min_length": 7,
        "max_length": 40,
        "min_charge": 2,
        "max_charge": 3,
        "cleavage_regex": "[KR][^_|$]",
        "missed_cleavages": 1,
        "max_var_mods": 1
    },

    "comment_modifications": "=== Modifications (variable and fixed, including PTMs and peptide-terminal labels) ===",
    "variable_mods": {
        "pattern": ["M"],
        "mass": [15.99491],
        "name": ["Unimod:35"]
    },
    "fixed_mods": {
        "pattern": ["C"],
        "mass": [57.021464],
        "name": ["Unimod:4"]
    }
}
```

**Total: ~20 user-visible parameters** (50% reduction from original!)

**Hidden defaults automatically applied:**
- **Library params:** `instrument_type`: "NONE", `prediction_model`: "altimeter"
- **Digestion params:** `add_decoys`: true, `decoy_method`: "shuffle", `entrapment_method`: "shuffle", `entrapment_r`: 0
- **NCE params:** {nce: 26.0, default_charge: 2, dynamic_nce: true}
- **Processing params:** `isotope_mod_groups`: [], `max_koina_requests`: 24, `max_koina_batch`: 1000, `match_lib_build_batch`: 100000, `predict_fragments`: true, `include_contaminants`: true

### After (Output config.json - Complete with Defaults)

When BuildSpecLib runs, it writes the **complete merged configuration** to `config.json`:

```json
{
    "library_path": "/path/to/output/my_library",
    "_lib_dir": "/path/to/output/my_library.poin",
    "_lib_name": "my_library",

    "calibration_raw_file": "/path/to/calibration/file.arrow",

    "fasta_paths": ["/path/to/file.fasta"],
    "fasta_names": ["HUMAN"],
    "fasta_header_regex_accessions": ["^\\w+\\|(\\w+(?:-\\d+)?)\\|"],
    "fasta_header_regex_genes": [" GN=(\\S+)"],
    "fasta_header_regex_proteins": ["^\\w+\\|\\w+?:-\\d+?\\|[^ ]+ (.*?) [^ ]+="],
    "fasta_header_regex_organisms": [" OS=([^ ]+.*?) [^ ]+="],
    "fasta_digest_params": {
        "min_length": 7,
        "max_length": 40,
        "min_charge": 2,
        "max_charge": 3,
        "cleavage_regex": "[KR][^_|$]",
        "missed_cleavages": 1,
        "max_var_mods": 1,
        "add_decoys": true,
        "decoy_method": "shuffle",
        "entrapment_method": "shuffle",
        "entrapment_r": 0
    },
    "modifications": {
        "variable_mods": {
            "pattern": ["M"],
            "mass": [15.99491],
            "name": ["Unimod:35"]
        },
        "fixed_mods": {
            "pattern": ["C"],
            "mass": [57.021464],
            "name": ["Unimod:4"]
        }
    },
    "library_params": {
        "auto_detect_frag_bounds": true,
        "frag_mz_min": 150.0,
        "frag_mz_max": 2020.0,
        "prec_mz_min": 390.0,
        "prec_mz_max": 1010.0,
        "instrument_type": "NONE",
        "prediction_model": "altimeter",
        "rt_bin_tol": 1.0,
        "frag_bin_tol_ppm": 10.0,
        "rank_to_score": [8, 4, 4, 2, 2, 1, 1],
        "y_start_index": 4,
        "b_start_index": 3,
        "y_start": 3,
        "b_start": 2,
        "include_p_index": false,
        "include_p": false,
        "max_frag_charge": 3,
        "max_frag_rank": 255,
        "length_to_frag_count_multiple": 2,
        "min_frag_intensity": 0.0,
        "include_isotope": false,
        "include_internal": false,
        "include_immonium": false,
        "include_neutral_diff": true
    },
    "nce_params": {
        "nce": 26.0,
        "default_charge": 2,
        "dynamic_nce": true
    },
    "processing": {
        "isotope_mod_groups": [],
        "max_koina_requests": 24,
        "max_koina_batch": 1000,
        "match_lib_build_batch": 100000,
        "include_contaminants": true,
        "predict_fragments": true
    },
    "_metadata": {
        "created_at": "2025-01-15T10:30:00",
        "pioneer_version": "0.1.13",
        "parameter_version": "simplified"
    }
}
```

**Key feature:** Users can see **exactly** what parameters were used, including all defaults.

---

## Implementation Priority

**High Priority (Phase 1-2):**
1. Output path simplification (biggest user pain point)
2. Remove redundant parameters from simplified version
3. Write complete merged parameters to config.json

**Medium Priority (Phase 3-4):**
4. Parameter naming improvements
5. SearchDIA-style parameter validation

**Low Priority (Phase 5-6):**
6. Comprehensive testing
7. Documentation updates

---

## Migration Strategy

### For Existing Users

Provide a migration script: `scripts/migrate_build_params.jl`

```julia
"""
Convert old BuildSpecLib parameter files to new format.
"""
function migrate_build_params(old_json_path::String, new_json_path::String)
    old_params = JSON.parsefile(old_json_path)

    new_params = Dict{String, Any}()

    # Merge old output paths into library_path
    if haskey(old_params, "out_dir") && haskey(old_params, "lib_name")
        new_params["library_path"] = joinpath(old_params["out_dir"], old_params["lib_name"])
    end

    # Copy other sections
    for key in ["fasta_paths", "fasta_names", "fasta_digest_params",
                "variable_mods", "fixed_mods", "library_params",
                "isotope_mod_groups", "max_koina_requests", "max_koina_batch",
                "match_lib_build_batch", "include_contaminants", "predict_fragments"]
        if haskey(old_params, key)
            new_params[key] = old_params[key]
        end
    end

    # Save new format
    open(new_json_path, "w") do io
        JSON.print(io, new_params, 2)
    end

    println("Migrated parameters saved to: $new_json_path")
end
```

---

## Benefits Summary

1. **Dramatically simpler for users:** 50% fewer visible parameters in simplified version (40 → 20)
2. **Less error-prone:** Single output path instead of 4 redundant parameters
3. **Better transparency:** Complete merged parameters written to config.json
4. **Robust validation:** SearchDIA-style bounds checking prevents invalid inputs
5. **Backward compatible:** Keep full parameter file for advanced users
6. **Clearer documentation:** Better comments and parameter grouping
7. **Smart defaults:** Processing, NCE, and decoy parameters have sensible defaults for 95% of use cases

---

*Document prepared for BuildSpecLib parameter simplification project*
*Branch: `simplify-buildspeclib-params`*
*Date: January 2025*
