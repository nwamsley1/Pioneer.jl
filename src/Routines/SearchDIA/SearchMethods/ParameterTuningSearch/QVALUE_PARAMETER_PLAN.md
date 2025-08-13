# Plan: Add Separate Q-Value Parameter for ParameterTuningSearch

## Problem Statement
Currently, ParameterTuningSearch uses the global q-value threshold from `global_params.scoring.q_value_threshold` (line 201 in types.jl). This forces the parameter tuning phase to use the same FDR threshold as the main search phases, which may not be optimal since:
- Parameter tuning needs to be more permissive to collect enough PSMs for robust fitting
- The global threshold (often 0.01) may be too stringent for initial parameter estimation

## Proposed Solution
Add a dedicated `max_q_value` parameter under `parameter_tuning->search_settings` in the JSON configuration, allowing independent control of the FDR threshold for parameter tuning.

## Implementation Plan

### 1. Update JSON Configuration Files

**File**: `data/example_config/defaultSearchParams.json`
**Location**: Under `parameter_tuning->search_settings` (after line 48)
```json
"search_settings": {
    "sample_rate": 0.02,
    "min_samples": 3500,
    "min_quad_tuning_psms": 5000,
    "min_quad_tuning_fragments": 3,
    "max_presearch_iters": 10,
    "frag_err_quantile": 0.01,
    "max_q_value": 0.01,  // NEW: Independent q-value for parameter tuning
    "topn_peaks": 200,
    "max_tolerance_ppm": 50.0,
    "initial_scan_count": 500,
    "max_parameter_tuning_scans": 8000,
    "max_frags_for_mass_err_estimation": 5
}
```

**Also update**: `data/ecoli_test/ecoli_test_params.json` with the same addition

### 2. Update Parameter Extraction in types.jl

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`

**Change at line 201**:
```julia
# OLD (line 201):
Float32(global_params.scoring.q_value_threshold),

# NEW:
# Extract max_q_value with fallback to global threshold
max_q_value = hasproperty(search_params, :max_q_value) ? 
    Float32(search_params.max_q_value) : 
    Float32(global_params.scoring.q_value_threshold)

# Then use in struct construction:
max_q_value,  # Use extracted value instead of global
```

**Add extraction code** (around line 191, after max_frags_for_mass_err_estimation):
```julia
# Extract max q-value for parameter tuning
max_q_value = hasproperty(search_params, :max_q_value) ? 
    Float32(search_params.max_q_value) : 
    Float32(global_params.scoring.q_value_threshold)
```

### 3. Update Parameter Validation

**File**: `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`

**Add validation** (after line 101, with other search_settings checks):
```julia
if haskey(search_settings, "max_q_value")
    check_param(search_settings, "max_q_value", Real)
    # Validate it's between 0 and 1
    if search_settings["max_q_value"] <= 0 || search_settings["max_q_value"] > 1
        error("max_q_value must be between 0 and 1")
    end
end
```

### 4. Add Documentation

**Update CLAUDE.md** in ParameterTuningSearch folder to document:
- Purpose of separate q-value threshold
- Recommended values (0.01 to 0.05 for parameter tuning)
- Relationship to global threshold

## Benefits

1. **Better Parameter Estimation**: Can use more permissive threshold (e.g., 0.05) during parameter tuning while maintaining strict threshold (0.01) for final results
2. **Faster Convergence**: More PSMs available for fitting mass error and RT models
3. **Flexibility**: Different experiments may need different tuning thresholds
4. **Backward Compatibility**: Falls back to global threshold if not specified

## Testing Strategy

1. **Default Behavior**: Verify that omitting the parameter uses global threshold
2. **Custom Value**: Test with max_q_value = 0.05 and verify more PSMs collected
3. **Validation**: Ensure invalid values (negative, >1) are rejected
4. **Integration**: Run full pipeline to verify parameter tuning improves with adjusted threshold

## Risk Assessment

- **Low Risk**: Fully backward compatible with fallback
- **No Breaking Changes**: Existing configurations continue to work
- **Clear Benefits**: More robust parameter estimation

## Recommended Default Values

- **Conservative**: 0.01 (same as global, current behavior)
- **Recommended**: 0.02 (slightly more permissive for parameter tuning)
- **Permissive**: 0.05 (when struggling to collect PSMs)

## Implementation Order

1. Update types.jl to extract and use the parameter
2. Add parameter validation in paramsChecks.jl
3. Update JSON configuration files
4. Test with ecoli_test data
5. Update documentation

## Success Metrics

- Parameter tuning converges in fewer iterations
- More consistent mass error model fits
- Works with both specified and default values
- No regression in existing functionality