# Plan: Add Intensity Filter Quantile Parameter for Mass Error Model Fitting

## Current Issue
In `fit_mass_err_model` function (utils.jl lines 572-573), fragments are filtered by intensity using a hardcoded quantile of 0.25:
```julia
low10 = quantile([fragment.intensity for fragment in fragments], 0.25) #Get median intensity
fragments = [fragment for fragment in fragments if fragment.intensity>low10]
```

This filtering helps remove low-quality/noise matches but the threshold is fixed and not configurable.

## Proposed Solution
Add a configurable `intensity_filter_quantile` parameter under `parameter_tuning->fragment_settings` in JSON configuration.

## Implementation Plan

### 1. Update JSON Configuration Files

**Files to update:**
- `data/example_config/defaultSearchParams.json`
- `data/ecoli_test/ecoli_test_params.json`

**Add to** `parameter_tuning->fragment_settings`:
```json
"fragment_settings": {
    "min_count": 7,
    "max_rank": 25,
    "tol_ppm": 20.0,
    "min_score": 22,
    "min_spectral_contrast": 0.9,
    "relative_improvement_threshold": 1.25,
    "min_log2_ratio": 1.5,
    "min_top_n": [3, 3],
    "n_isotopes": 1,
    "intensity_filter_quantile": 0.25  // NEW: Filter fragments below this intensity quantile
}
```

### 2. Update Parameter Extraction in types.jl

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`

**Add field to struct** (around line 147):
```julia
struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # ... existing fields ...
    intensity_filter_quantile::Float32  # NEW: Quantile for intensity filtering
    # ... rest of fields ...
end
```

**Add extraction logic** (around line 170, where other fragment settings are extracted):
```julia
# Extract intensity filter quantile for mass error model fitting
intensity_filter_quantile = hasproperty(frag_params, :intensity_filter_quantile) ? 
    Float32(frag_params.intensity_filter_quantile) : 
    Float32(0.25)  # Default to 0.25 (25th percentile)
```

**Add to constructor** (in the appropriate position):
```julia
intensity_filter_quantile,  # Add to struct construction
```

**Add accessor function** (around line 233):
```julia
getIntensityFilterQuantile(params::ParameterTuningSearchParameters) = params.intensity_filter_quantile
```

### 3. Update Parameter Validation

**File**: `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`

**Add validation** (in the fragment_settings validation section, around line 77):
```julia
# In the parameter_tuning fragment_settings validation
tuning_frag = tuning_params["fragment_settings"]
# ... existing checks ...
if haskey(tuning_frag, "intensity_filter_quantile")
    check_param(tuning_frag, "intensity_filter_quantile", Real)
    # Validate it's between 0 and 1
    if tuning_frag["intensity_filter_quantile"] < 0 || tuning_frag["intensity_filter_quantile"] >= 1
        error("parameter_tuning.fragment_settings.intensity_filter_quantile must be between 0 and 1 (got $(tuning_frag["intensity_filter_quantile"]))")
    end
end
```

### 4. Update fit_mass_err_model Function

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

**Current code** (lines 553-592):
```julia
function fit_mass_err_model(
    params::P,
    fragments::Vector{FragmentMatch{Float32}}
) where {P<:FragmentIndexSearchParameters}
    # ... calc_ppm_error function ...
    
    # Filter by intensity - HARDCODED
    low10 = quantile([fragment.intensity for fragment in fragments], 0.25)
    fragments = [fragment for fragment in fragments if fragment.intensity>low10]
    
    # ... rest of function ...
```

**Updated code**:
```julia
function fit_mass_err_model(
    params::P,
    fragments::Vector{FragmentMatch{Float32}}
) where {P<:FragmentIndexSearchParameters}
    """
    Fits mass error model from fragment matches with intensity filtering.
    
    # Process
    1. Filters fragments below intensity quantile threshold
    2. Calculates median mass offset
    3. Determines tolerance bounds using error quantiles
    """
    
    # Filter fragments by intensity to remove low-quality matches
    if length(fragments) > 0
        intensity_threshold = getIntensityFilterQuantile(params)
        if intensity_threshold > 0.0
            intensities = [fragment.intensity for fragment in fragments]
            min_intensity = quantile(intensities, intensity_threshold)
            fragments = filter(f -> f.intensity > min_intensity, fragments)
            @info "Filtered to $(length(fragments)) fragments above $(round(intensity_threshold*100, digits=1))th percentile intensity"
        end
    end
    
    # Check if we have enough fragments after filtering
    if length(fragments) == 0
        @warn "No fragments remaining after intensity filtering"
        return MassErrorModel(0.0f0, (30.0f0, 30.0f0)), Float32[]
    end
    
    # Calculate PPM errors
    ppm_errs = [(f.match_mz - f.theoretical_mz)/(f.theoretical_mz/1e6) for f in fragments]
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Calculate error bounds using configured quantile
    frag_err_quantile = getFragErrQuantile(params)
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)
    
    # Diagnostic logging
    @info "Mass error model: offset=$(round(mass_err, digits=2)) ppm, " *
          "tolerance=(-$(round(abs(l_bound), digits=1)), +$(round(abs(r_bound), digits=1))) ppm, " *
          "MAD=$(round(mad(ppm_errs, normalize=false)*4, digits=1)) ppm"
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end
```

### 5. Clean Up Helper Function

Remove the inline `calc_ppm_error` function definition and use a more direct approach or move it to a proper location if needed elsewhere.

## Benefits

1. **Flexibility**: Different datasets may benefit from different filtering thresholds
2. **Robustness**: Can adjust threshold based on data quality
3. **Transparency**: Makes the filtering explicit and configurable
4. **Backward Compatible**: Defaults to current behavior (0.25)

## Recommended Values

- **Conservative** (less filtering): 0.1 (10th percentile)
- **Default**: 0.25 (25th percentile - current behavior)
- **Aggressive** (more filtering): 0.5 (median)
- **No filtering**: 0.0

## Testing Strategy

1. **Default behavior**: Verify 0.25 quantile when not specified
2. **Custom values**: Test with 0.1, 0.5 to see impact on tolerance estimates
3. **Edge cases**: Test with 0.0 (no filtering) and 0.9 (extreme filtering)
4. **Validation**: Ensure invalid values (<0 or >=1) are rejected

## Risk Assessment

- **Low Risk**: Fully backward compatible
- **Tunable**: Can be adjusted per experiment
- **Clear Impact**: Directly affects mass tolerance estimation

## Implementation Order

1. Update types.jl to add field and extraction
2. Add parameter validation
3. Update JSON configuration files
4. Refactor fit_mass_err_model function
5. Test with various quantile values

## Success Metrics

- Cleaner, more readable code
- Configurable intensity filtering
- Better mass tolerance estimates for noisy data
- No regression with default value (0.25)