# Scan Count Refactoring Plan for ParameterTuningSearch

## Overview
This document outlines the plan to improve scan count parameters for the ParameterTuningSearch convergence loop, making them more intuitive and configurable via JSON.

## Current Issues

1. **Unclear naming**: `expanded_scan_count` doesn't clearly indicate it's the maximum
2. **Fixed increment strategy**: Currently adds fixed chunks instead of doubling
3. **Hardcoded defaults**: Initial scan count defaults to 2500 if not specified
4. **Non-intuitive progression**: Current implementation doesn't follow a clear growth pattern
5. **Hardcoded convergence threshold**: `check_convergence` uses hardcoded 1000 instead of `min_samples` from JSON

## Proposed Solution

### 1. Parameter Renaming
- **Old**: `expanded_scan_count`
- **New**: `max_parameter_tuning_scans`
- **Rationale**: Clearly indicates this is the maximum number of scans for parameter tuning

### 2. Scan Growth Strategy
Implement exponential growth (doubling) strategy:
- **Iteration 0**: Start with `initial_scan_count` scans
- **Iteration 1**: Double to 2 × initial_scan_count
- **Iteration 2**: Double to 4 × initial_scan_count
- **Iteration 3**: Double to 8 × initial_scan_count
- **Continue** until reaching `max_parameter_tuning_scans`

### 3. Use Existing min_samples Parameter
The `min_samples` parameter already exists in JSON and should be used for convergence checking instead of hardcoded 1000.

### 4. Configuration via JSON
All parameters fully configurable:
```json
"search_settings": {
    "min_samples": 3500,                 // Minimum PSMs for convergence check (already exists)
    "initial_scan_count": 500,          // Starting point
    "max_parameter_tuning_scans": 8000  // Maximum limit
}
```

## Implementation Details

### Phase 1: Fix min_samples Usage in check_convergence

Currently, `check_convergence` has a hardcoded threshold of 1000 PSMs:
```julia
if size(psms, 1) < 1000
    @info "size(psms, 1) < 1000, skipping convergence check"
    return false 
end
```

This should use the `min_samples` parameter that already exists in the JSON and ParameterTuningSearchParameters.

#### Files to Modify:

1. **utils.jl** (Line ~972 in check_convergence):
```julia
# Old
function check_convergence(
    psms,
    new_mass_err_model::MassErrorModel,
    old_mass_err_model::MassErrorModel
)
    if size(psms, 1) < 1000
        @info "size(psms, 1) < 1000, skipping convergence check"
        return false 
    end
    # ... rest of function
end

# New - add params parameter
function check_convergence(
    psms,
    new_mass_err_model::MassErrorModel,
    old_mass_err_model::MassErrorModel,
    params::ParameterTuningSearchParameters
)
    min_psms_required = getMinPsms(params)  # Use configurable threshold
    if size(psms, 1) < min_psms_required
        @info "size(psms, 1) < $min_psms_required, skipping convergence check"
        return false 
    end
    # ... rest of function
end
```

2. **ParameterTuningSearch.jl** - Update all calls to check_convergence:
   - Line ~331: Add params argument
   - Line ~438: Add params argument  
   - Line ~464: Add params argument

```julia
# Old
if check_convergence(psms, mass_err_model, results.mass_err_model[])

# New
if check_convergence(psms, mass_err_model, results.mass_err_model[], params)
```

### Phase 2: Rename Parameter

#### Files to Modify:

1. **types.jl** (Line ~145):
```julia
# Old
expanded_scan_count::Int64

# New
max_parameter_tuning_scans::Int64
```

2. **types.jl** (Line ~177):
```julia
# Old
expanded_scan_count = hasproperty(search_params, :expanded_scan_count) ? 
    Int64(search_params.expanded_scan_count) : Int64(10000)

# New
max_parameter_tuning_scans = hasproperty(search_params, :max_parameter_tuning_scans) ? 
    Int64(search_params.max_parameter_tuning_scans) : Int64(8000)
```

3. **types.jl** (Line ~224):
```julia
# Old
getExpandedScanCount(params::ParameterTuningSearchParameters) = params.expanded_scan_count

# New
getMaxParameterTuningScans(params::ParameterTuningSearchParameters) = params.max_parameter_tuning_scans
```

4. **paramsChecks.jl** (Line ~94):
```julia
# Old
if haskey(search_settings, "expanded_scan_count")
    check_param(search_settings, "expanded_scan_count", Integer)
end

# New
if haskey(search_settings, "max_parameter_tuning_scans")
    check_param(search_settings, "max_parameter_tuning_scans", Integer)
end
```

### Phase 2: Implement Doubling Strategy

#### ParameterTuningSearch.jl (Lines 296-314):

**Current Implementation**:
```julia
if n == 0
    additional_scans = params.initial_scan_count
else
    additional_scans = min(params.expanded_scan_count - length(filtered_spectra), 2500)
end
```

**New Implementation**:
```julia
# Calculate scans for this iteration using doubling strategy
if n == 0
    # First iteration: use initial scan count
    additional_scans = getInitialScanCount(params)
else
    # Subsequent iterations: double the current count up to max
    current_scans = length(filtered_spectra)
    target_scans = min(current_scans * 2, getMaxParameterTuningScans(params))
    additional_scans = target_scans - current_scans
    
    # Stop if we've reached the maximum
    if additional_scans <= 0
        @info "Reached maximum scan count of $(getMaxParameterTuningScans(params))"
        break  # Exit the convergence loop
    end
end

@info "Iteration $(n+1): Adding $additional_scans scans (total will be $(length(filtered_spectra) + additional_scans))"
```

### Phase 4: Update JSON Files

#### defaultSearchParams.json & ecoli_test_params.json:
```json
"search_settings": {
    "sample_rate": 0.02,
    "min_samples": 3500,                 // Used for convergence check (not hardcoded 1000)
    "initial_scan_count": 500,           // Reduced from 2500 for faster start
    "max_parameter_tuning_scans": 8000,  // Renamed from expanded_scan_count
    "max_frags_for_mass_err_estimation": 5,
    // ... other settings
}
```

Note: `min_samples` already exists and will now be properly used in `check_convergence` instead of the hardcoded 1000.

## Example Scan Progression

With `initial_scan_count: 500` and `max_parameter_tuning_scans: 8000`:

| Iteration | Scans Added | Total Scans | Notes |
|-----------|-------------|-------------|-------|
| 0 | 500 | 500 | Initial attempt |
| 1 | 500 | 1000 | Doubled |
| 2 | 1000 | 2000 | Doubled |
| 3 | 2000 | 4000 | Doubled |
| 4 | 4000 | 8000 | Hits maximum |

## Benefits

1. **Intuitive progression**: Clear exponential growth pattern
2. **Faster initial attempts**: Can start with fewer scans
3. **Configurable limits**: Both initial and max are tunable
4. **Better naming**: Parameters clearly indicate their purpose
5. **Efficient search**: Balances speed vs thoroughness
6. **Configurable convergence**: min_samples threshold now configurable via JSON instead of hardcoded

## Migration Notes

### For Users:
- Update JSON files to use `max_parameter_tuning_scans` instead of `expanded_scan_count`
- Consider reducing `initial_scan_count` from 2500 to 500-1000 for faster initial attempts
- The maximum can be increased if needed for difficult samples

### Backward Compatibility:
- The code will still check for `expanded_scan_count` for compatibility
- If both are present, `max_parameter_tuning_scans` takes precedence
- Default values ensure existing configs continue to work

## Testing Checklist

- [ ] Verify parameter extraction from JSON
- [ ] Confirm doubling progression: 500 → 1000 → 2000 → 4000 → 8000
- [ ] Test that iteration stops at maximum
- [ ] Ensure convergence still succeeds with new strategy
- [ ] Check backward compatibility with old parameter name
- [ ] Verify logging shows correct scan counts
- [ ] Confirm min_samples is used for convergence check (not hardcoded 1000)
- [ ] Test with different min_samples values (e.g., 2000, 3500, 5000)

## Implementation Order

1. Create this planning document ✓
2. Fix min_samples usage in check_convergence (Phase 1)
3. Rename parameter throughout codebase (Phase 2)
4. Implement doubling strategy (Phase 3)
5. Update JSON files (Phase 4)
6. Test with example data
7. Update documentation

## Notes

- The doubling strategy ensures we quickly explore the scan space
- Starting with fewer scans (500) can significantly speed up initial attempts
- The maximum (8000) provides a safety limit while allowing thorough searches
- This change maintains backward compatibility while improving usability