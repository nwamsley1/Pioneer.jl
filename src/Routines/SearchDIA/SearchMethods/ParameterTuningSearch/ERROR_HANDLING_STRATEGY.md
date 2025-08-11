# Error Handling Strategy for ParameterTuningSearch Convergence Loop

## Problem Statement

The current convergence loop in ParameterTuningSearch can fail when:
1. No PSMs are collected in Strategy 1 (scan expansion)
2. Functions like `fit_mass_err_model` or `fit_irt_model` receive empty data
3. `check_convergence` is called with empty PSMs

When any of these occur, the loop should continue to the next strategy rather than stopping with an error.

## Current Failure Points

### 1. Empty PSM Collection
- **Location**: Line 319-323
- **Issue**: When `collect_psms` returns empty DataFrame, subsequent operations fail:
  - `get_matched_fragments` (line 322) - may return empty fragments
  - `fit_mass_err_model` (line 323) - fails with empty fragments
  - `check_convergence` (line 324) - returns false with <1000 PSMs
  - `fit_irt_model` (line 338) - fails with empty PSMs

### 2. Function Dependencies
These functions expect non-empty input:
- `fit_mass_err_model(params, fragments)` - needs fragments with ppm errors
- `fit_irt_model(params, psms)` - needs PSMs with RT/iRT columns
- `get_matched_fragments(...)` - needs valid PSMs

### 3. Check Convergence Logic
- Returns `false` if PSM count < 1000 (line 972-974)
- Returns `false` if tolerance improvements < 10%
- Never throws errors, so this is safe

## Proposed Error Handling Strategy

### Strategy Overview
1. **Graceful Degradation**: Each strategy should handle empty/insufficient data gracefully
2. **Continue on Failure**: If a strategy fails, proceed to the next one
3. **Maintain Strategy Logic**: Keep the existing convergence strategy unchanged
4. **Final Fallback**: If all strategies fail, use conservative defaults

### Implementation Plan

#### Phase 1: Safe PSM Collection Check
After collecting PSMs (lines 319-321), check if collection was successful before attempting model fitting:

```julia
psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
final_psm_count = size(psms, 1)
@info "Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans"

# Check if we have enough PSMs to proceed with model fitting
if final_psm_count > 0
    # Try to fit models
    fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
    if length(fragments) > 0
        mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
        # Check convergence only if we have valid models
        if check_convergence(psms, mass_err_model, results.mass_err_model[])
            # Convergence success path
            ...
        end
    else
        @info "No fragments matched, proceeding to Strategy 2"
    end
else
    @info "No PSMs collected in Strategy 1, proceeding to Strategy 2"
end
```

#### Phase 2: Add Safety Checks Before Model Fitting

Add simple checks to prevent errors when data is insufficient:

```julia
# Before calling fit_mass_err_model
if length(fragments) == 0
    @info "No fragments available for mass error fitting"
    # Skip to next strategy
else
    mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
end

# Before calling fit_irt_model
if size(psms, 1) == 0
    @info "No PSMs available for RT model fitting"
    # Skip RT model fitting
else
    rt_model_tuple = fit_irt_model(params, psms)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_tuple)
end
```

#### Phase 3: Wrap Critical Sections in Try-Catch

Add minimal error handling to prevent crashes while maintaining the strategy:

```julia
# Strategy 1: Scan Expansion
@info "Strategy 1: Collecting PSMs with current parameters"
psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
final_psm_count = size(psms, 1)

if final_psm_count > 0
    try
        fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
        if length(fragments) > 0
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            if check_convergence(psms, mass_err_model, results.mass_err_model[])
                # Handle convergence
                converged = true
                break
            end
        end
    catch e
        @warn "Error in Strategy 1 model fitting: $e"
        # Continue to Strategy 2
    end
else
    @info "No PSMs in Strategy 1, continuing to Strategy 2"
end

# Continue with Strategy 2 and 3 similarly...
```

### Error Recovery Mechanisms

#### 1. Minimum Data Requirements
Define minimum thresholds for proceeding:
- Minimum PSMs for convergence check: 1000 (existing)
- Minimum PSMs for any model fitting: 1 (prevent empty DataFrame errors)
- Minimum fragments for mass error: 1 (prevent empty vector errors)

#### 2. Fallback Behavior
When fitting fails due to insufficient data:
- Continue to next strategy in the same iteration
- Let the existing convergence loop logic handle progression
- After all iterations, existing fallback mechanism kicks in

#### 3. Iteration Control
- Maintain existing 5-iteration limit
- Each iteration adds more scans as designed
- Existing fallback parameters used if all iterations fail

### Benefits of This Approach

1. **Minimal Changes**: Preserves existing convergence strategy
2. **Robustness**: Prevents crashes from empty data
3. **Diagnostic Clarity**: Clear logging when strategies skip due to no data
4. **Simplicity**: No complex recovery logic needed

### Testing Strategy

1. **Test with edge cases**:
   - Files that yield no initial PSMs
   - Files with very few MS2 scans
   - Files where all PSMs are decoys

2. **Verify behavior**:
   - Loop continues when Strategy 1 yields no PSMs
   - Strategies 2 and 3 execute even if Strategy 1 fails
   - Final fallback parameters are applied correctly

### Implementation Priority

**Immediate fixes needed**:
1. Add check for `final_psm_count > 0` before model fitting in each strategy
2. Add check for `length(fragments) > 0` before `fit_mass_err_model`
3. Skip convergence check when insufficient data
4. Continue to next strategy instead of erroring

## Code Implementation Example

Here's how the main convergence loop would look with simplified error handling:

```julia
while n < N
    # Reset parameters if needed (existing code)
    if n > 0
        results.mass_err_model[] = create_capped_mass_model(
            initial_mass_offset, initial_tolerance, initial_tolerance, max_tolerance
        )
        setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
    end
    
    # Strategy 1: Scan Expansion
    # Add scans
    additional_scans = n == 0 ? params.initial_scan_count : 
                      min(params.expanded_scan_count - length(filtered_spectra), 2500)
    append!(filtered_spectra; max_additional_scans = additional_scans)
    
    # Collect PSMs
    @info "Strategy 1: Collecting PSMs with current parameters"
    psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    final_psm_count = size(psms, 1)
    @info "Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans"
    
    # Check if we can proceed with model fitting
    if final_psm_count > 0
        fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
        
        if length(fragments) > 0
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            
            if check_convergence(psms, mass_err_model, results.mass_err_model[])
                # Convergence success
                results.mass_err_model[] = mass_err_model
                # Store ppm_errs and RT model...
                converged = true
                break
            end
        else
            @info "No fragments matched in Strategy 1, continuing to Strategy 2"
        end
    else
        @info "No PSMs collected in Strategy 1, continuing to Strategy 2"
    end
    
    # Strategy 2: Tolerance Expansion (always runs if Strategy 1 didn't converge)
    @info "Strategy 2: Expanding mass tolerance"
    # Expand tolerance (existing code)
    current_left_tol = getLeftTol(results.mass_err_model[])
    current_right_tol = getRightTol(results.mass_err_model[])
    results.mass_err_model[] = create_capped_mass_model(
        getMassOffset(results.mass_err_model[]),
        current_left_tol * 1.5f0,
        current_right_tol * 1.5f0,
        max_tolerance
    )
    setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
    
    psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    final_psm_count = size(psms, 1)
    @info "Collected $final_psm_count PSMs with expanded tolerance"
    
    # Only try to fit models if we have PSMs
    if final_psm_count > 0
        fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
        if length(fragments) > 0
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            # Don't check convergence here, continue to Strategy 3
        end
    end
    
    # Strategy 3: Bias Adjustment (always runs)
    @info "Strategy 3: Adjusting mass bias"
    # Adjust bias if we have a new model from Strategy 2
    if final_psm_count > 0 && length(fragments) > 0
        new_bias = getMassOffset(mass_err_model)
        results.mass_err_model[] = create_capped_mass_model(
            new_bias,
            getLeftTol(results.mass_err_model[]),
            getRightTol(results.mass_err_model[]),
            max_tolerance
        )
        setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
        
        psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
        final_psm_count = size(psms, 1)
        
        if final_psm_count > 0
            fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            if length(fragments) > 0
                mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
                if check_convergence(psms, mass_err_model, results.mass_err_model[])
                    # Convergence success
                    results.mass_err_model[] = mass_err_model
                    # Store ppm_errs and RT model...
                    converged = true
                    break
                end
            end
        end
    end
    
    n += 1
    @info "Iteration $(n) complete. Moving to next iteration..."
end
```

## Summary

This simplified error handling strategy ensures the ParameterTuningSearch convergence loop:
1. **Continues when no PSMs are collected** - Moves to next strategy instead of crashing
2. **Skips model fitting with empty data** - Prevents errors from empty vectors/DataFrames
3. **Maintains existing convergence logic** - No changes to the overall strategy
4. **Provides clear diagnostics** - Logs when strategies skip due to insufficient data

The implementation requires minimal changes while significantly improving robustness. The existing fallback mechanism handles cases where all iterations fail to converge.