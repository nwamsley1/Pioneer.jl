# Fix Plan: Missing Parameters in Simplified Build Library Template

## Problem Report

**Issue:** When Emily tried predicting a library using the auto-generated simplified parameters, the system complained about missing "nce_params" and she had to manually add them to the JSON file.

**Additional Request:** Add fragment and precursor min/max m/z to the simplified version, as it will be common to turn off `auto_detect_frag_bounds` when users don't have a converted raw file ready.

## Root Cause Analysis

### 1. The Validation System
**Location:** `src/Routines/BuildSpecLib/utils/check_params.jl`

The `check_params_bsp()` function performs comprehensive validation of all library building parameters. Key validation points:

- **Lines 97-101**: Always validates `nce_params` section (nce, default_charge, dynamic_nce)
- **Lines 117-120**: Always validates fragment/precursor m/z bounds (frag_mz_min, frag_mz_max, prec_mz_min, prec_mz_max)

These validations run **regardless of whether the parameters are simplified or full**.

### 2. The Defaults System
**Location:** `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl`

The defaults system works by:
1. Loading the appropriate template JSON (`defaultBuildLibParamsSimplified.json` or `defaultBuildLibParams.json`)
2. Removing user-specific keys (paths, names)
3. Merging user parameters over the defaults

**Critical Issue**: If a parameter isn't in the template JSON, it won't be in the defaults, and validation will fail.

### 3. The Simplified Template Gap
**File:** `assets/example_config/defaultBuildLibParamsSimplified.json`

The simplified template was missing:

| Parameter | Location | Purpose | Required For |
|-----------|----------|---------|--------------|
| `nce_params` | Top level | Collision energy settings for fragmentation prediction | Always validated |
| `frag_mz_min` | library_params | Minimum fragment m/z bound | Fallback when auto_detect_frag_bounds fails or is disabled |
| `frag_mz_max` | library_params | Maximum fragment m/z bound | Fallback when auto_detect_frag_bounds fails or is disabled |
| `prec_mz_min` | library_params | Minimum precursor m/z bound | Fallback when auto_detect_frag_bounds fails or is disabled |
| `prec_mz_max` | library_params | Maximum precursor m/z bound | Fallback when auto_detect_frag_bounds fails or is disabled |

### 4. Why M/Z Bounds Are Always Required

**Location:** `src/Routines/BuildSpecLib/fragments/get_frag_bounds.jl`

The `get_fragment_bounds()` function signature (lines 58-62):
```julia
function get_fragment_bounds(
    auto_detect_frag_bounds::Bool,
    frag_bounds_detection_raw_file_path::String,
    default_frag_bounds::Tuple{Float32, Float32},
    default_precursor_bounds::Tuple{Float32, Float32})
```

Even when `auto_detect_frag_bounds = true`:
- The bounds serve as **fallback values** if detection fails
- They define the valid m/z range for the instrument
- They're required inputs to the function

When `auto_detect_frag_bounds = false`:
- These values are used directly without any auto-detection
- Essential for users without converted raw files

## Solution Implemented

### Changes Made to `defaultBuildLibParamsSimplified.json`

#### 1. Added NCE Parameters Section
```json
"comment_collision_energy": "=== Collision energy settings ===",
"nce_params": {
    "nce": 26.0,
    "default_charge": 2,
    "dynamic_nce": true
},
```

**Location**: Added after `fixed_mods` section, before `comment_processing`

**Values Chosen**: Match the full template defaults
- `nce: 26.0` - Standard normalized collision energy for typical DDA experiments
- `default_charge: 2` - Most common precursor charge state
- `dynamic_nce: true` - Enable charge-dependent NCE adjustment

#### 2. Added M/Z Bounds to Library Params
```json
"library_params": {
    "auto_detect_frag_bounds": true,
    "calibration_raw_file": "/path/to/calibration/file.arrow",
    "frag_mz_min": 150.0,
    "frag_mz_max": 2020.0,
    "prec_mz_min": 390.0,
    "prec_mz_max": 1010.0,
    "instrument_type": "NONE",
    "prediction_model": "altimeter"
},
```

**Values Chosen**: Match the full template defaults
- `frag_mz_min: 150.0` - Typical low-end cutoff for Orbitrap/QTOF instruments
- `frag_mz_max: 2020.0` - Covers most peptide fragments
- `prec_mz_min: 390.0` - Typical minimum precursor m/z for DDA
- `prec_mz_max: 1010.0` - Covers typical tryptic peptide range (2+ to 4+ charge states)

#### 3. Enhanced Documentation Comment
```json
"comment_library": "=== Library settings ===\n    Note: auto_detect_frag_bounds automatically detects m/z bounds from calibration_raw_file.\n    Set to false if you don't have a converted raw file - will use frag/prec m/z bounds below",
```

**Purpose**: Helps users understand:
- What `auto_detect_frag_bounds` does
- When to disable it
- What the m/z bounds are used for

## Benefits of This Solution

### Immediate Fixes
- Resolves "nce_params missing" error - Users can now use auto-generated simplified params without manual edits
- Prevents validation failures - All required parameters are now present in the template

### User Experience Improvements
- Makes disabling auto_detect_frag_bounds easy - Just change one boolean, bounds are already there
- Provides sensible defaults - Values work for most DDA experiments on modern instruments
- Better documentation - Comments explain when and why to modify settings

### Consistency
- Matches full template structure - Same parameter names and organization
- Uses same default values - Consistent experience across simplified/full templates
- No code changes needed - Pure configuration fix, no risk of introducing bugs

## Use Cases Now Supported

### Case 1: User WITH converted raw file (Original workflow)
```json
"auto_detect_frag_bounds": true,
"calibration_raw_file": "/path/to/real/file.arrow",
```
- System automatically detects optimal m/z bounds from data
- Manual bounds serve as fallback if detection fails

### Case 2: User WITHOUT converted raw file (Emily's case)
```json
"auto_detect_frag_bounds": false,
"calibration_raw_file": "/path/to/calibration/file.arrow",  # Can be dummy/placeholder
```
- System uses the explicit frag/prec m/z bounds directly
- No raw file conversion required
- User can adjust bounds for their specific instrument if needed

## Testing Recommendations

### 1. Basic Validation Test
```julia
# Generate simplified params
params = GetBuildLibParams(
    fasta_paths = ["/path/to/test.fasta"],
    out_dir = "/tmp/test",
    lib_name = "test_lib",
    simplified = true
)

# Should not throw "nce_params missing" error
```

### 2. With Auto-Detection Enabled
- Use actual converted raw file
- Verify bounds are detected and fallback values aren't needed

### 3. With Auto-Detection Disabled
- Set auto_detect_frag_bounds = false
- Verify explicit bounds are used correctly
- Test library building completes without errors

### 4. Edge Cases
- Invalid calibration file path (should fall back to manual bounds gracefully)
- Extreme m/z values (verify validation catches unreasonable bounds)
- Different instrument types (verify bounds are appropriate)

## Files Modified

| File | Location | Changes | Status |
|------|----------|---------|--------|
| `defaultBuildLibParamsSimplified.json` | `assets/example_config/` | Added `nce_params` section and m/z bounds to `library_params`, enhanced comments | Completed |

## Related Code Locations

For future reference, the key code locations involved in this system:

| Component | File | Lines | Purpose |
|-----------|------|-------|---------|
| Parameter validation | `src/Routines/BuildSpecLib/utils/check_params.jl` | 97-101, 117-120 | Validates nce_params and m/z bounds |
| Defaults loading | `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl` | 15-36 | Loads default templates |
| Defaults merging | `src/Routines/BuildSpecLib/utils/buildParamDefaults.jl` | 44-61 | Merges user params with defaults |
| Fragment bounds detection | `src/Routines/BuildSpecLib/fragments/get_frag_bounds.jl` | 58-95 | Auto-detects or uses manual bounds |
| Template usage | `src/Routines/GenerateParams.jl` | 287-463 | Generates parameter files from templates |

## Future Considerations

### Potential Improvements
1. **Make calibration_raw_file truly optional** - When `auto_detect_frag_bounds = false`, the calibration file shouldn't be needed at all
2. **Add instrument-specific presets** - Different default bounds for Orbitrap, QTOF, TOF, etc.
3. **Validation warnings** - Warn users if m/z bounds seem unusual for their instrument type
4. **Better error messages** - When auto-detection fails, explain why and suggest checking the bounds

### Documentation Updates Needed
- User guide explaining when to use simplified vs. full templates
- Explanation of `auto_detect_frag_bounds` feature and its requirements
- Instrument-specific recommendations for m/z bounds
- Troubleshooting guide for parameter validation errors

## Conclusion

This fix resolves the immediate issue Emily encountered and makes Pioneer more user-friendly for researchers who want to start predicting libraries without converting raw files first. The solution maintains backward compatibility while extending functionality in a natural way.

The key insight is that even when auto-detection is enabled, fallback bounds must be specified because:
1. Auto-detection can fail (missing file, corrupted data, etc.)
2. Validation always checks for these parameters
3. They define the instrument's valid operating range

By including sensible defaults in the simplified template, we enable both use cases (with and without raw file conversion) while maintaining a single, clean template structure.
