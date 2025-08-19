# Best Attempt Fallback Implementation Plan

## Overview

Currently, when ParameterTuningSearch fails to converge (doesn't reach `min_samples` PSMs), it either borrows parameters from another file or uses conservative defaults. This plan implements a better fallback strategy: use the mass error model and RT alignment from the iteration that produced the most PSMs, even if it didn't meet the convergence criteria.

## Problem Statement

When processing challenging data:
1. None of the iterations may achieve the minimum PSM count for convergence
2. Current fallback to conservative defaults (±50 ppm) is often too broad
3. Borrowing from other files may not be appropriate if files have different characteristics
4. We do generate useful models during iterations, but currently discard them if convergence fails

## Proposed Solution

Track the "best attempt" throughout all iterations and use it as the fallback when convergence fails.

### Definition of "Best Attempt"
The best attempt is the iteration that produced:
1. The highest number of PSMs (primary criterion)
2. With a reasonable mass error model (secondary validation)

## Implementation Design

### 1. Data Structure for Tracking Best Attempt

Add to `IterationState` struct:
```julia
mutable struct IterationState
    # Existing fields...
    
    # New fields for best attempt tracking
    best_psm_count::Int64
    best_mass_error_model::Union{Nothing, MassErrorModel}
    best_rt_model::Union{Nothing, Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32}}
    best_ppm_errs::Union{Nothing, Vector{Float32}}
    best_phase::Int64
    best_score::UInt8
    best_iteration::Int64
    best_scan_count::Int64
end
```

Initialize with:
```julia
best_psm_count = 0
best_mass_error_model = nothing
best_rt_model = nothing
best_ppm_errs = nothing
best_phase = 0
best_score = 0
best_iteration = 0
best_scan_count = 0
```

### 2. Update Best Attempt During Processing

After each PSM collection and model fitting (in `run_single_phase`):

```julia
function update_best_attempt!(
    iteration_state::IterationState,
    psm_count::Int64,
    mass_err_model::Union{Nothing, MassErrorModel},
    rt_model_data::Union{Nothing, Tuple},
    ppm_errs::Union{Nothing, Vector{Float32}},
    phase::Int64,
    score::UInt8,
    iteration::Int64,
    scan_count::Int64
)
    # Only update if we have more PSMs than the previous best
    if psm_count > iteration_state.best_psm_count && mass_err_model !== nothing
        iteration_state.best_psm_count = psm_count
        iteration_state.best_mass_error_model = mass_err_model
        iteration_state.best_rt_model = rt_model_data
        iteration_state.best_ppm_errs = ppm_errs
        iteration_state.best_phase = phase
        iteration_state.best_score = score
        iteration_state.best_iteration = iteration
        iteration_state.best_scan_count = scan_count
        
        @info "    📊 New best attempt: $(psm_count) PSMs " *
              "(Phase $phase, Score $score, Iteration $iteration)"
    end
end
```

### 3. Integration Points

#### A. After Initial Attempt (Before Iteration Loop)
```julia
# After fitting models from initial attempt
if mass_err_model !== nothing && psm_count > 0
    # Fit RT model even if not converged
    rt_model_data = fit_irt_model(params, psms_initial)
    
    # Update best attempt
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
        phase, min_score, 0, current_scan_count
    )
end
```

#### B. After Each Iteration
```julia
# After Step 6: FIT FINAL MODELS
if mass_err_model !== nothing && psm_count > 0
    # Fit RT model
    rt_model_data = fit_irt_model(params, psms_adjusted)
    
    # Update best attempt
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
        phase, min_score, iter, current_scan_count
    )
end
```

### 4. Modified Fallback Logic

In `process_file!`, replace the current fallback with:

```julia
if !converged
    @warn "Failed to converge for file $ms_file_idx after " *
          "$(iteration_state.scan_attempt) attempts"
    
    # Check if we have a best attempt to use
    if iteration_state.best_mass_error_model !== nothing
        @info "Using best attempt as fallback: " *
              "$(iteration_state.best_psm_count) PSMs from " *
              "Phase $(iteration_state.best_phase), " *
              "Score $(iteration_state.best_score), " *
              "Iteration $(iteration_state.best_iteration)"
        
        # Apply best attempt models
        setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
        
        if iteration_state.best_rt_model !== nothing
            set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                iteration_state.best_rt_model)
        end
        
        # Store best ppm errors for plotting
        if iteration_state.best_ppm_errs !== nothing
            resize!(results.ppm_errs, 0)
            append!(results.ppm_errs, iteration_state.best_ppm_errs)
        end
        
        # Update warnings
        push!(warnings, "BEST_ATTEMPT_FALLBACK: Used parameters from best attempt " *
                       "($(iteration_state.best_psm_count) PSMs, " *
                       "Phase $(iteration_state.best_phase), " *
                       "Score $(iteration_state.best_score))")
        
        # Set final PSM count
        final_psm_count = iteration_state.best_psm_count
        
    else
        # No valid attempts at all - use conservative defaults or borrow
        @warn "No valid attempts found, falling back to conservative defaults"
        
        fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
            search_context, params, ms_file_idx, warnings
        )
        
        setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)
        setRtConversionModel!(search_context, ms_file_idx, fallback_rt_model)
        
        push!(warnings, "CONSERVATIVE_FALLBACK: No valid attempts, used defaults")
    end
end
```

### 5. Update store_final_results!

Modify to handle the best attempt case:

```julia
function store_final_results!(results, search_context, params, ms_file_idx, 
                              converged, n_attempts, final_psm_count, warnings,
                              iteration_state)  # Add iteration_state parameter
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Determine convergence type
    convergence_type = if converged
        "CONVERGED"
    elseif iteration_state.best_psm_count > 0
        "BEST_ATTEMPT ($(iteration_state.best_psm_count) PSMs)"
    else
        "FAILED"
    end
    
    # Log final status
    @info "Final parameter tuning status for file $ms_file_idx ($parsed_fname):"
    @info "  - Status: $convergence_type"
    @info "  - RT points: $(length(results.rt))"
    @info "  - iRT points: $(length(results.irt))" 
    @info "  - PPM errors: $(length(results.ppm_errs))"
    @info "  - Final PSM count: $final_psm_count"
    
    if !converged && iteration_state.best_psm_count > 0
        @info "  - Best attempt: Phase $(iteration_state.best_phase), " *
              "Score $(iteration_state.best_score), " *
              "Iteration $(iteration_state.best_iteration)"
    end
    
    # Rest of function remains the same...
end
```

## Benefits

1. **Always uses real data**: Instead of conservative defaults, uses actual calibration from the data
2. **Graceful degradation**: Performance degrades gradually with data quality
3. **Better for difficult samples**: Low-abundance samples still get reasonable parameters
4. **Transparent**: Clear logging shows what parameters were used and why

## Risk Assessment

### Low Risk
- Only affects files that fail to converge
- Always logs clearly when using best attempt
- Falls back to conservative defaults if no valid attempts

### Medium Risk  
- Best attempt might have poor quality models
- Mitigation: Validate that mass offset is reasonable before accepting

## Testing Strategy

1. **Unit Test**: Mock iteration with varying PSM counts, verify best is selected
2. **Integration Test**: Use low-quality data that won't converge, verify best attempt is used
3. **Regression Test**: Ensure files that currently converge still behave the same

## Implementation Steps

1. **Phase 1**: Add IterationState fields and update_best_attempt! function
2. **Phase 2**: Integrate tracking into run_single_phase
3. **Phase 3**: Update fallback logic in process_file!
4. **Phase 4**: Update logging and warnings
5. **Phase 5**: Test with challenging datasets

## Alternative Approaches Considered

1. **Store all attempts**: Rejected - too memory intensive
2. **Score-based selection**: Rejected - PSM count is more reliable indicator
3. **Weighted average of attempts**: Rejected - single best model is simpler and more interpretable

## Success Criteria

- Files that fail convergence use best attempt parameters
- Clear logging indicates when best attempt is used
- Parameters from best attempt produce better results than conservative defaults
- No regression in files that currently converge successfully

## Example Log Output

```
[ Info] Phase 1, Score 22, Iteration 1: 45 PSMs collected
[ Info]     📊 New best attempt: 45 PSMs (Phase 1, Score 22, Iteration 1)
[ Info] Phase 1, Score 17, Iteration 2: 78 PSMs collected  
[ Info]     📊 New best attempt: 78 PSMs (Phase 1, Score 17, Iteration 2)
[ Info] Phase 2, Score 22, Iteration 1: 65 PSMs collected
[ Info] Phase completed without convergence (tried all score thresholds)
[ Warn] Failed to converge for file 1 after 2 attempts
[ Info] Using best attempt as fallback: 78 PSMs from Phase 1, Score 17, Iteration 2
[ Info] Final parameter tuning status for file 1 (sample.raw):
[ Info]   - Status: BEST_ATTEMPT (78 PSMs)
[ Info]   - Best attempt: Phase 1, Score 17, Iteration 2
```