# Mass Error Model Usage in QuadTuningSearch - Comprehensive Analysis

## Executive Summary

**CRITICAL FINDING**: QuadTuningSearch is correctly using the fitted MassErrorModel from ParameterTuningSearch, NOT the `frag_tol_ppm` parameter. The `frag_tol_ppm` field in QuadTuningSearchParameters appears to be vestigial and unused in the actual search process.

## Analysis Overview

This analysis investigates how mass tolerances are handled in QuadTuningSearch, specifically:
1. Where `frag_tol_ppm` parameter is defined and used
2. How the `collect_psms` function determines mass tolerances 
3. The relationship between parameter values and fitted MassErrorModel
4. Comparison with v0.2.1 behavior

## Key Findings

### 1. frag_tol_ppm Parameter Usage

**Occurrences in QuadTuningSearch codebase:**
- `QuadTuningSearch.jl:35` - Documentation example (30.0 ppm)
- `QuadTuningSearch.jl:97` - Field definition in `QuadTuningSearchParameters`
- `QuadTuningSearch.jl:261` - Logging output only

**CRITICAL**: The `frag_tol_ppm` parameter is **NEVER USED** in the actual search logic. It only appears in:
1. Documentation examples
2. Parameter struct definition
3. Info logging output

### 2. Mass Error Model Flow in collect_psms

The actual mass tolerance determination follows this path:

```
collect_psms() 
  └── library_search(spectra, search_context, params, ms_file_idx)
      └── LibrarySearch(..., getMassErrorModel(search_context, ms_file_idx), ...)
          └── searchFragmentIndex(..., mem, ...)  # mem = MassErrorModel
              └── searchScan!(..., mass_err_model, ...)
                  └── getCorrectedMz(mass_err_model, mass)
                  └── getMzBoundsReverse(mass_err_model, corrected_mz)
```

**Key Function: `searchScan!` in `CommonSearchUtils/queryFragmentIndex.jl:288`**

```julia
function searchScan!(..., mass_err_model::MassErrorModel, ...)
    for mass in masses
        # ACTUAL tolerance calculation using fitted model
        corrected_mz = getCorrectedMz(mass_err_model, mass)
        frag_min, frag_max = getMzBoundsReverse(mass_err_model, corrected_mz)
        # Search fragments within these bounds
    end
end
```

### 3. MassErrorModel vs frag_tol_ppm

**How it ACTUALLY works:**
1. ParameterTuningSearch (runs first) fits a MassErrorModel based on data
2. Model stored in SearchContext via `setMassErrorModel!(search_context, ms_file_idx, fitted_model)`
3. QuadTuningSearch retrieves fitted model via `getMassErrorModel(search_context, ms_file_idx)`
4. `searchScan!` uses the fitted model's `getCorrectedMz()` and `getMzBoundsReverse()` methods
5. **`frag_tol_ppm` parameter is completely ignored**

**What MassErrorModel contains (from ParameterTuningSearch):**
- Mass offset correction (ppm)
- Left tolerance bound (ppm) 
- Right tolerance bound (ppm)
- Asymmetric tolerances allowed
- Data-driven, not hardcoded

**What frag_tol_ppm contains:**
- Single symmetric tolerance value (20 ppm from `init_mass_tol_ppm`)
- Unused in actual search

## Comparison with v0.2.1

### v0.2.1 Behavior:
```julia
# v0.2.1 constructor
Float32(frag_params.tol_ppm),  # Used directly
```
- `frag_tol_ppm` was actually passed to search functions
- Likely used 30 ppm hardcoded value from JSON
- No fitted MassErrorModel integration

### Current Version:
```julia
# Current constructor  
init_tol,  # 20 ppm from init_mass_tol_ppm, but UNUSED in search
```
- `frag_tol_ppm` stored but never used
- Fitted MassErrorModel takes precedence
- Data-driven tolerances (could be 3-50+ ppm depending on instrument)

## Evidence Supporting This Analysis

### 1. Library Search Function Signature
```julia
function library_search(...) where {P<:SearchParameters}
    return vcat(LibrarySearch(
        ...,
        getMassErrorModel(search_context, ms_file_idx),  # ← Fitted model used
        search_parameters,  # ← Contains unused frag_tol_ppm
        ...
    )...)
end
```

### 2. searchScan! Function
```julia
function searchScan!(..., mass_err_model::MassErrorModel, ...)
    # Uses mass_err_model methods directly
    corrected_mz = getCorrectedMz(mass_err_model, mass)
    frag_min, frag_max = getMzBoundsReverse(mass_err_model, corrected_mz)
```

### 3. No Direct Parameter Usage
Search in all QuadTuningSearch files shows `frag_tol_ppm` is never passed to:
- `getCorrectedMz()`
- `getMzBounds()`
- `searchScan!()`
- Any tolerance calculation function

## Current Behavior Analysis

### What the logs show:
```
frag_tol_ppm: 20
```

### What actually happens:
1. ParameterTuningSearch fits model (e.g., offset=1.2 ppm, left_tol=8.3 ppm, right_tol=9.1 ppm)
2. QuadTuningSearch logs "20" but uses the fitted (8.3, 9.1) tolerances
3. All fragment matching uses the fitted asymmetric tolerances
4. The "20" value is never used

## Implications

### 1. Current Implementation is CORRECT
- QuadTuningSearch properly uses data-driven tolerances from ParameterTuningSearch
- Mass error correction and asymmetric bounds are applied
- Better than v0.2.1's hardcoded approach

### 2. Parameter Logging is MISLEADING  
- Logs show "frag_tol_ppm: 20" but system uses fitted model
- Could confuse users about actual tolerances being applied
- Fitted model tolerances could be very different (3-50+ ppm)

### 3. Code Cleanup Opportunity
- `frag_tol_ppm` field could be removed from QuadTuningSearchParameters
- Or replaced with fitted model values for accurate logging
- Parameter extraction code could be simplified

## Recommendations

### 1. Immediate: Fix Misleading Logs
Replace parameter logging with actual fitted model values:

```julia
# Current (misleading)
@info "frag_tol_ppm: $(params.frag_tol_ppm)"

# Proposed (accurate)
fitted_model = getMassErrorModel(search_context, ms_file_idx)  
@info "Using fitted MassErrorModel from ParameterTuningSearch:" *
      "\n  - Mass offset: $(getMassOffset(fitted_model)) ppm" *
      "\n  - Left tolerance: $(getLeftTolerance(fitted_model)) ppm" *  
      "\n  - Right tolerance: $(getRightTolerance(fitted_model)) ppm"
```

### 2. Future: Remove Vestigial Parameter
- Remove `frag_tol_ppm` field from QuadTuningSearchParameters
- Remove `init_mass_tol_ppm` extraction from JSON
- Simplify constructor to not handle mass tolerance

### 3. Documentation Update
- Update QuadTuningSearch documentation to clarify it uses fitted models
- Remove misleading `frag_tol_ppm` examples
- Emphasize dependency on ParameterTuningSearch

## Conclusion

**The QuadTuningSearch is working correctly** - it uses the sophisticated, data-driven MassErrorModel fitted by ParameterTuningSearch rather than crude hardcoded tolerances. The `frag_tol_ppm` parameter is a vestigial field that serves no functional purpose but creates confusion through misleading log output.

This represents a significant improvement over v0.2.1, where QuadTuningSearch used primitive hardcoded tolerances. The current implementation benefits from:

- **Data-driven tolerances** based on actual mass error characteristics  
- **Asymmetric bounds** (left ≠ right tolerance)
- **Mass offset correction** for systematic bias
- **Intensity-dependent corrections** if implemented in the fitted model

The only issue is the misleading logging that suggests the hardcoded 20 ppm value is being used when it's actually ignored in favor of the much more sophisticated fitted model.