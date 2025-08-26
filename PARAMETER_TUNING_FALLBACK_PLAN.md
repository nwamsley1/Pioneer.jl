# ParameterTuningSearch Fallback Enhancement Plan

## Executive Summary

This plan details modifications to the ParameterTuningSearch module to improve handling of convergence failures. The current implementation already tracks the best iteration but doesn't fully utilize it. We will enhance the fallback strategy to always use the best iteration's parameters when convergence fails, including testing with 1.5x mass tolerance expansion.

## Current Implementation Analysis

### Existing Best Iteration Tracking

The system already tracks best attempts in `IterationState`:
```julia
struct IterationState
    # ... other fields ...
    best_psm_count::Int64
    best_mass_error_model::Union{Nothing, MassErrorModel}
    best_rt_model::Union{Nothing, SplineRtConversionModel}
    best_ppm_errs::Union{Nothing, Vector{Float32}}
    best_phase::Int64
    best_score::UInt8
    best_iteration::Int64
    best_scan_count::Int64
end
```

### Current Fallback Behavior

When convergence fails (lines 830-889 in ParameterTuningSearch.jl):

1. **Best Attempt Fallback** (lines 835-869):
   - If `best_mass_error_model` exists, uses it
   - Applies best RT model or defaults to identity
   - Logs: "BEST_ATTEMPT_FALLBACK: Used parameters from best attempt"

2. **Conservative Fallback** (lines 871-889):
   - Only triggers when NO valid attempts exist
   - Calls `get_fallback_parameters` to borrow from other files
   - Falls back to ±50 ppm tolerance if no files available

### Issue with Current Implementation

The system doesn't:
1. Test the best iteration with 1.5x tolerance expansion
2. Generate proper diagnostic plots for best iteration parameters
3. Provide detailed warnings about which iteration was used

## Proposed Enhancements

### 1. Enhanced Best Iteration Fallback

#### A. Tolerance Expansion Testing for Best Iteration

**Location**: After line 869 in ParameterTuningSearch.jl

**Implementation**:
```julia
# After applying best attempt models (line 869)
if iteration_state.best_mass_error_model !== nothing
    # Test 1.5x tolerance expansion on best iteration
    expanded_model = MassErrorModel(
        getMassOffset(iteration_state.best_mass_error_model),
        (getLeftTol(iteration_state.best_mass_error_model) * 1.5f0,
         getRightTol(iteration_state.best_mass_error_model) * 1.5f0)
    )
    
    # Apply expanded model and collect PSMs
    setMassErrorModel!(search_context, ms_file_idx, expanded_model)
    expanded_psms, expanded_ppm_errs = collect_psms_with_model(
        filtered_spectra, search_context, params, ms_file_idx
    )
    
    # Check if expansion yields significantly more PSMs
    expansion_threshold = iteration_state.best_psm_count * 1.1
    if size(expanded_psms, 1) > expansion_threshold
        # Use expanded results
        iteration_state.best_mass_error_model = fit_mass_err_model(expanded_ppm_errs, params)
        iteration_state.best_psm_count = size(expanded_psms, 1)
        iteration_state.best_ppm_errs = expanded_ppm_errs
        
        # Update warning message
        push!(warnings, "EXPANDED_BEST_ATTEMPT: Tolerance expansion yielded " *
                       "$(size(expanded_psms, 1)) PSMs ($(round((size(expanded_psms, 1) - iteration_state.best_psm_count) / iteration_state.best_psm_count * 100, digits=1))% increase)")
    end
end
```

#### B. Enhanced Warning Messages

**Current warning**:
```
"BEST_ATTEMPT_FALLBACK: Used parameters from best attempt (X PSMs, Phase Y, Score Z)"
```

**Proposed warning**:
```
"⚠️ Parameter tuning did not converge after X attempts.
Using best iteration: Phase Y (bias=Z ppm), Score threshold W, Iteration V
Best iteration yielded U PSMs with mass offset T ppm and tolerance ±S ppm
[Optional: Tolerance expansion increased PSMs by R%]"
```

### 2. Improved Diagnostic Plotting

#### A. Generate Plots from Best Iteration Data

**Location**: In `process_search_results!` function

**Current behavior**: Generates fallback plots with no data
**Proposed behavior**: Use stored best iteration data

```julia
function generate_best_iteration_plots(
    results::ParameterTuningSearchResults,
    iteration_state::IterationState,
    fname::String
)
    # RT plot with best iteration data
    if length(results.rt) > 0
        rt_plot = generate_rt_plot(results, "$fname (Best Iteration)")
    else
        rt_plot = generate_fallback_rt_plot_with_info(
            results, fname, iteration_state
        )
    end
    
    # Mass error plot with best iteration data
    if length(results.ppm_errs) > 0
        mass_plot = generate_mass_error_plot(
            results, "$fname (Best Iteration)"
        )
    else
        mass_plot = generate_fallback_mass_plot_with_info(
            results, fname, iteration_state
        )
    end
    
    return rt_plot, mass_plot
end
```

#### B. Enhanced Fallback Plot Information

Add iteration details to fallback plots:

```julia
function generate_fallback_rt_plot_with_info(
    results, fname, iteration_state
)
    p = # ... existing plot setup ...
    
    # Add annotation with best iteration details
    annotation_text = "Best Iteration Results:\n" *
                     "Phase $(iteration_state.best_phase), " *
                     "Score $(iteration_state.best_score)\n" *
                     "$(iteration_state.best_psm_count) PSMs found\n" *
                     "Iteration $(iteration_state.best_iteration) of " *
                     "$(iteration_state.total_iterations) total"
    
    Plots.annotate!(60, 20, text(annotation_text, :center, 10))
    return p
end
```

### 3. Helper Function for PSM Collection

Add a helper function to collect PSMs with a specific model:

```julia
function collect_psms_with_model(
    filtered_spectra::FilteredMassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64
)
    # Perform library search with current model
    psms = library_search(
        filtered_spectra,
        search_context,
        params,
        ms_file_idx
    )
    
    # Add columns and score
    add_tuning_search_columns!(psms, ...)
    filter_and_score_psms!(psms, params)
    
    # Get matched fragments for error calculation
    if size(psms, 1) > 0
        matched_frags = get_matched_fragments(...)
        ppm_errs = calculate_ppm_errors(matched_frags)
    else
        ppm_errs = Float32[]
    end
    
    return psms, ppm_errs
end
```

### 4. Update IterationState During Search

Ensure best iteration is properly tracked:

```julia
# In process_file! after each iteration attempt
if psm_count > iteration_state.best_psm_count
    iteration_state.best_psm_count = psm_count
    iteration_state.best_mass_error_model = mass_err_model
    iteration_state.best_rt_model = rt_model
    iteration_state.best_ppm_errs = copy(ppm_errs)
    iteration_state.best_phase = current_phase
    iteration_state.best_score = current_score
    iteration_state.best_iteration = current_iteration
    iteration_state.best_scan_count = current_scan_count
end
```

### 5. Diagnostic Summary Enhancement

Update the diagnostics to show more detail about best iterations:

```julia
struct ParameterTuningStatus
    # ... existing fields ...
    used_best_iteration::Bool  # New field
    best_iteration_details::Union{Nothing, NamedTuple}  # New field
end
```

## Implementation Steps

### Phase 1: Core Fallback Enhancement
1. Modify fallback logic to always test tolerance expansion on best iteration
2. Update warning messages with detailed information
3. Ensure best iteration data is properly stored

### Phase 2: Plotting Improvements
1. Create plot generation functions that use best iteration data
2. Add informative annotations to fallback plots
3. Ensure plots are properly stored in results

### Phase 3: Helper Functions
1. Implement `collect_psms_with_model` helper
2. Add utility functions for tolerance expansion testing
3. Create diagnostic summary functions

### Phase 4: Testing and Validation
1. Test with files that consistently fail to converge
2. Verify tolerance expansion works correctly
3. Ensure plots are generated properly
4. Validate warning messages are informative

## Benefits

1. **Better Parameter Recovery**: Files that fail to converge will use the most successful iteration
2. **Tolerance Optimization**: Automatic testing of expanded tolerance ensures maximum PSM recovery
3. **Improved Diagnostics**: Clear indication of what parameters were used and why
4. **User Transparency**: Detailed warnings help users understand the fallback behavior
5. **Quality Plots**: Diagnostic plots will show actual data when available

## Backward Compatibility

All changes are backward compatible:
- Existing behavior is preserved when convergence succeeds
- Conservative fallback still available when no iterations produce PSMs
- No changes to public API or parameter files

## Testing Strategy

1. **Unit Tests**:
   - Test best iteration tracking
   - Test tolerance expansion logic
   - Test plot generation with partial data

2. **Integration Tests**:
   - Use files known to have convergence issues
   - Verify fallback parameters are applied correctly
   - Check that downstream methods work with fallback parameters

3. **Edge Cases**:
   - No PSMs in any iteration
   - Only one iteration with PSMs
   - Tolerance expansion produces fewer PSMs
   - Missing RT data but valid mass error

## Summary

This plan enhances the existing best iteration tracking to provide a robust fallback strategy when parameter tuning fails to converge. By testing tolerance expansion and providing detailed diagnostics, users will get better results and clearer understanding of the tuning process, even when convergence criteria aren't met.