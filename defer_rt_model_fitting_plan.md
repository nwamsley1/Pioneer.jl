# Plan: Defer RT Model Fitting to End of Parameter Tuning

## Executive Summary

Refactor ParameterTuningSearch to **store PSMs during iteration** and **fit RT model only once at the end**, rather than fitting RT models multiple times during the search process.

## Current Behavior (Inefficient)

### RT Model Fitting Locations
Currently, `fit_irt_model()` is called **3 times per file** in multiple scenarios:

1. **Line 550** (`ParameterTuningSearch.jl`): After initial PSM collection attempt
   ```julia
   rt_psms = filter_top_psms_per_precursor(psms_initial, 3)
   rt_model_data = fit_irt_model(params, rt_psms)
   update_best_attempt!(iteration_state, ..., rt_model_data, ...)
   ```

2. **Line 633** (`ParameterTuningSearch.jl`): After each iteration's adjusted PSM collection
   ```julia
   rt_psms = filter_top_psms_per_precursor(psms_adjusted, 3)
   rt_model_data = fit_irt_model(params, rt_psms)
   update_best_attempt!(iteration_state, ..., rt_model_data, ...)
   ```

3. **Line 359** (`ParameterTuningSearch.jl`): At final convergence
   ```julia
   rt_psms = filter_top_psms_per_precursor(final_psms, 3)
   rt_model_data = fit_irt_model(params, rt_psms)
   set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
   ```

### Performance Impact
- **Problem**: For a file with 3 phases × 3 iterations = 9 search attempts, RT model could be fit **9-10 times**
- **Computational cost**: Each fit involves:
  - Filtering top PSMs per precursor
  - DataFrame creation for robust regression
  - Iteratively reweighted least squares (MEstimator)
  - Spline fitting for comparison
  - AIC calculation
- **Why wasteful**: Only the FINAL fitted model is actually used

## Proposed Behavior (Efficient)

### New Strategy
1. **During iteration**: Store the filtered PSMs DataFrame (not the fitted model)
2. **At the end**: Fit RT model exactly once using the best PSM set

### Changes Required

#### 1. Modify `IterationState` Structure
**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`

**Current** (line 155):
```julia
best_rt_model::Union{Nothing, Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}}
```

**New**:
```julia
best_psms::Union{Nothing, DataFrame}  # Store PSMs instead of fitted model
```

**Justification**:
- DataFrames are lightweight references (column vectors)
- PSMs from best attempt typically 20-100 rows
- No need to carry fitted model through iterations

#### 2. Update `update_best_attempt!` Function
**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`

**Current signature** (line 114):
```julia
function update_best_attempt!(
    iteration_state::IterationState,
    psm_count::Int64,
    mass_err_model::Union{Nothing, MassErrorModel},
    rt_model_data::Union{Nothing, Tuple},  # ← Remove this
    ppm_errs::Union{Nothing, Vector{Float64}},
    phase::Int64,
    score::UInt8,
    iteration::Int64,
    scan_count::Int64
)
```

**New signature**:
```julia
function update_best_attempt!(
    iteration_state::IterationState,
    psm_count::Int64,
    mass_err_model::Union{Nothing, MassErrorModel},
    best_psms::Union{Nothing, DataFrame},  # ← Add PSMs
    ppm_errs::Union{Nothing, Vector{Float64}},
    phase::Int64,
    score::UInt8,
    iteration::Int64,
    scan_count::Int64
)
```

**Current body** (line 126-135):
```julia
if psm_count > iteration_state.best_psm_count && mass_err_model !== nothing
    iteration_state.best_psm_count = psm_count
    iteration_state.best_mass_error_model = mass_err_model
    iteration_state.best_rt_model = rt_model_data  # ← Change this
    iteration_state.best_ppm_errs = ppm_errs
    iteration_state.best_phase = phase
    iteration_state.best_score = score
    iteration_state.best_iteration = iteration
    iteration_state.best_scan_count = scan_count
end
```

**New body**:
```julia
if psm_count > iteration_state.best_psm_count && mass_err_model !== nothing
    iteration_state.best_psm_count = psm_count
    iteration_state.best_mass_error_model = mass_err_model
    iteration_state.best_psms = best_psms  # ← Store PSMs not model
    iteration_state.best_ppm_errs = ppm_errs
    iteration_state.best_phase = phase
    iteration_state.best_score = score
    iteration_state.best_iteration = iteration
    iteration_state.best_scan_count = scan_count
end
```

#### 3. Remove RT Fitting from Iteration Loops
**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`

**Location 1** (lines 546-556):
```julia
# REMOVE THIS BLOCK:
if psm_count > 0 && !isempty(psms_initial)
    rt_psms = filter_top_psms_per_precursor(psms_initial, 3)
    rt_model_data = fit_irt_model(params, rt_psms)
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
        phase, min_score, 0, iteration_state.current_scan_count
    )
end
```

**REPLACE WITH**:
```julia
if psm_count > 0 && !isempty(psms_initial)
    # Filter top PSMs but DON'T fit model yet
    rt_psms = filter_top_psms_per_precursor(psms_initial, 3)
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_psms, ppm_errs,
        phase, min_score, 0, iteration_state.current_scan_count
    )
end
```

**Location 2** (lines 629-640):
```julia
# REMOVE THIS BLOCK:
if mass_err_model !== nothing && psm_count > 0 && !isempty(psms_adjusted)
    rt_psms = filter_top_psms_per_precursor(psms_adjusted, 3)
    rt_model_data = fit_irt_model(params, rt_psms)
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
        phase, min_score, iter, iteration_state.current_scan_count
    )
end
```

**REPLACE WITH**:
```julia
if mass_err_model !== nothing && psm_count > 0 && !isempty(psms_adjusted)
    # Filter top PSMs but DON'T fit model yet
    rt_psms = filter_top_psms_per_precursor(psms_adjusted, 3)
    update_best_attempt!(
        iteration_state, psm_count, mass_err_model, rt_psms, ppm_errs,
        phase, min_score, iter, iteration_state.current_scan_count
    )
end
```

#### 4. Fit RT Model Once at End
**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`

**Location**: `process_convergence` function (lines 350-365)

**Current** (line 356-363):
```julia
if !isempty(final_psms)
    rt_psms = filter_top_psms_per_precursor(final_psms, 3)
    rt_model_data = fit_irt_model(params, rt_psms)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
else
    @user_warn "No PSMs available for RT model at convergence - this should not happen"
end
```

**New**:
```julia
# Determine which PSMs to use for RT model fitting
psms_for_rt_model = if !isempty(final_psms)
    # Converged successfully - use final PSMs
    filter_top_psms_per_precursor(final_psms, 3)
elseif iteration_state.best_psms !== nothing
    # Did not converge but have best attempt - use those
    @user_warn "Using best attempt PSMs ($(iteration_state.best_psm_count) PSMs) for RT model fitting"
    iteration_state.best_psms
else
    # No PSMs available at all - this is a failure case
    @user_warn "No PSMs available for RT model fitting"
    nothing
end

# Fit RT model exactly once using best available PSMs
if psms_for_rt_model !== nothing
    rt_model_data = fit_irt_model(params, psms_for_rt_model)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
end
```

#### 5. Handle Fallback Case When No Convergence
**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`

**Location**: After `run_all_phases_with_scan_count` returns false (around line 800-850)

Currently, when convergence fails, the code falls back to best attempt. We need to ensure RT model is fitted from best PSMs in this case.

**Add after fallback logic**:
```julia
# If we didn't converge, process_convergence wasn't called, so fit RT model here
if !converged && iteration_state.best_psms !== nothing
    @user_info "Fitting RT model from best attempt ($(iteration_state.best_psm_count) PSMs)"
    rt_model_data = fit_irt_model(params, iteration_state.best_psms)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
end
```

## Implementation Checklist

### Phase 1: Data Structure Changes
- [ ] Modify `IterationState` in `types.jl`: Replace `best_rt_model` with `best_psms::Union{Nothing, DataFrame}`
- [ ] Update `IterationState` constructor to initialize `best_psms = nothing`

### Phase 2: Function Signature Updates
- [ ] Update `update_best_attempt!` signature: Replace `rt_model_data` parameter with `best_psms`
- [ ] Update function body: Store PSMs instead of fitted model

### Phase 3: Remove Redundant RT Fitting
- [ ] Remove RT model fitting from line 550 (initial attempt tracking)
- [ ] Remove RT model fitting from line 633 (iteration tracking)
- [ ] Update both `update_best_attempt!` calls to pass PSMs instead

### Phase 4: Centralized RT Fitting
- [ ] Update `process_convergence` to handle both convergence and best-attempt cases
- [ ] Fit RT model exactly once using best available PSMs
- [ ] Add proper logging to indicate which PSM set was used

### Phase 5: Fallback Handling
- [ ] Ensure RT model is fitted when falling back to best attempt
- [ ] Add diagnostic logging for fallback cases

### Phase 6: Testing
- [ ] Test convergence case (RT model from final PSMs)
- [ ] Test best-attempt fallback (RT model from best PSMs)
- [ ] Test complete failure case (no PSMs available)
- [ ] Verify performance improvement (fewer `fit_irt_model` calls)

## Expected Benefits

### Performance
- **Reduced computation**: 1 RT model fit instead of 9-10 per file
- **Estimated speedup**: 5-10% faster parameter tuning overall
- **Scalability**: More significant for datasets with many files (100+)

### Code Quality
- **Separation of concerns**: PSM collection vs model fitting
- **Clearer logic**: Explicit single point where RT model is created
- **Easier debugging**: One place to add logging/diagnostics

### Correctness
- **Best data used**: Always use the PSMs from best attempt
- **Consistent behavior**: Same RT fitting logic for convergence and fallback

## Potential Risks & Mitigation

### Risk 1: Memory Usage
**Concern**: Storing DataFrame instead of fitted model might use more memory

**Mitigation**:
- PSM DataFrames are typically small (20-100 rows)
- Only one DataFrame stored per file being processed
- DataFrames are column-oriented and efficient
- Negligible compared to spectral data already in memory

### Risk 2: Breaking Existing Behavior
**Concern**: Code expects `best_rt_model` to exist

**Mitigation**:
- Carefully grep for all uses of `iteration_state.best_rt_model`
- Ensure fallback paths still work correctly
- Add tests for convergence and fallback cases

### Risk 3: Edge Cases
**Concern**: What if PSMs become invalid by the time we fit?

**Mitigation**:
- PSMs are filtered before storage
- No mutations occur during iteration
- Fallback to identity model if fitting fails (already implemented)

## Testing Strategy

### Unit Tests
```julia
# Test 1: Verify IterationState stores PSMs
state = IterationState()
psms = DataFrame(rt=[1.0, 2.0], irt_predicted=[10.0, 20.0])
update_best_attempt!(state, 2, mass_model, psms, nothing, 1, UInt8(22), 0, 100)
@test state.best_psms === psms
@test state.best_psm_count == 2

# Test 2: Verify RT model only fitted once
# Mock fit_irt_model and count calls
```

### Integration Tests
```julia
# Test convergence path
result = run_parameter_tuning(convergent_dataset)
@test result.converged
@test result.rt_model !== nothing

# Test best-attempt fallback
result = run_parameter_tuning(difficult_dataset)
@test !result.converged
@test result.rt_model !== nothing  # Should still have model from best attempt
```

### Performance Tests
```julia
# Measure RT fitting calls
@time SearchDIA("test_params.json")
# Verify fit_irt_model called once per file
```

## Migration Notes

This is a **non-breaking change** from the perspective of downstream code:
- `ParameterTuningSearchResults` structure unchanged
- RT models still stored in `SearchContext` the same way
- External API (`getRtToIrtModel()`, etc.) unchanged

Only internal iteration logic is refactored.

## Appendix: Call Graph

### Before (Multiple Fits)
```
run_single_phase_attempt
├─ collect PSMs (initial)
├─ filter_top_psms_per_precursor
├─ fit_irt_model ← FIT #1
└─ update_best_attempt! (stores fitted model)

run_single_phase (iteration loop)
├─ collect PSMs (adjusted)
├─ filter_top_psms_per_precursor
├─ fit_irt_model ← FIT #2, #3, #4... (per iteration)
└─ update_best_attempt! (stores fitted model)

process_convergence
├─ filter_top_psms_per_precursor
├─ fit_irt_model ← FIT #N (final)
└─ set_rt_to_irt_model!
```

### After (Single Fit)
```
run_single_phase_attempt
├─ collect PSMs (initial)
├─ filter_top_psms_per_precursor
└─ update_best_attempt! (stores PSMs)

run_single_phase (iteration loop)
├─ collect PSMs (adjusted)
├─ filter_top_psms_per_precursor
└─ update_best_attempt! (stores PSMs)

process_convergence
├─ Determine best PSMs (converged or best attempt)
├─ fit_irt_model ← FIT ONCE
└─ set_rt_to_irt_model!
```

## Timeline

- **Phase 1-2**: 30 minutes (data structure changes)
- **Phase 3**: 30 minutes (remove redundant fitting)
- **Phase 4**: 45 minutes (centralized fitting logic)
- **Phase 5**: 15 minutes (fallback handling)
- **Phase 6**: 1 hour (testing)

**Total estimated time**: 3 hours

## Conclusion

This refactoring provides measurable performance benefits with minimal risk. The changes are well-isolated to the ParameterTuningSearch module and don't affect downstream code. The clearer separation of concerns also makes the code easier to maintain and debug.
