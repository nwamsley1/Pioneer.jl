# Detailed Plan: Tolerance Expansion Post-Convergence Optimization

## Status: DETAILED IMPLEMENTATION PLAN

## Executive Summary

After achieving convergence with the current mass tolerance used for PSM collection, test whether expanding that collection tolerance by 50% yields more PSMs. If it does, refit the mass error model using the expanded PSM set. This ensures we don't miss valid PSMs due to conservative tolerance estimates during the collection phase.

## Detailed Problem Analysis

### Current Convergence Behavior

The current implementation has a subtle but important limitation:

1. **Collection Phase**: We use tolerance T₁ (e.g., 20 ppm) to collect PSMs
2. **Fitting Phase**: We fit a model from collected PSMs, yielding tolerance T₂ (e.g., 18 ppm)
3. **Problem**: T₂ is constrained by what was discovered with T₁
4. **Consequence**: We might miss valid PSMs that fall in the range T₁ to 1.5×T₁

### Key Insight

The fitted tolerance (T₂) tells us about the PSMs we found, but the collection tolerance (T₁) determines what PSMs we could find. By testing an expanded collection tolerance after convergence, we can verify whether our collection window was optimal.

## Detailed Implementation Plan

### 1. Track Collection Tolerance

First, we need to track the actual tolerance used for PSM collection, not just the fitted model.

#### Add to IterationState (types.jl)
```julia
mutable struct IterationState
    current_phase::Int64
    current_iteration_in_phase::Int64
    total_iterations::Int64
    phase_bias_shifts::Vector{Float32}
    converged::Bool
    collection_tolerance::Float32  # NEW: Track tolerance used for collection
    
    function IterationState()
        new(1, 0, 0, Float32[], false, 0.0f0)
    end
end
```

#### Update tolerance tracking in execute_strategy
```julia
function execute_strategy(strategy_num::Int, filtered_spectra, spectra, 
                         search_context, params, ms_file_idx, iteration_state)
    # Track the actual collection tolerance
    current_model = getMassErrorModel(search_context, ms_file_idx)
    iteration_state.collection_tolerance = (getLeftTol(current_model) + getRightTol(current_model)) / 2.0f0
    
    # ... rest of strategy execution
end
```

### 2. Core Expansion Function

Create a new function to handle the tolerance expansion test:

```julia
"""
    test_tolerance_expansion!(results, search_context, params, ms_file_idx,
                             current_psms, current_model, current_ppm_errs,
                             collection_tolerance, filtered_spectra, spectra)

After convergence, test if expanding the collection tolerance yields more PSMs.
If successful, updates results with the expanded PSM set and refitted model.

# Arguments
- `collection_tolerance`: The tolerance that was used to collect current_psms
- Other arguments as in standard search functions

# Returns
- `(final_psms, final_model, final_ppm_errs, was_expanded::Bool)`
"""
function test_tolerance_expansion!(
    results::ParameterTuningSearchResults,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    current_psms::DataFrame,
    current_model::MassErrorModel,
    current_ppm_errs::Vector{Float32},
    collection_tolerance::Float32,
    filtered_spectra::FilteredMassSpecData,
    spectra::MassSpecData
)
    # Configuration parameters
    expansion_factor = get_expansion_factor(params)  # Default 1.5
    
    current_psm_count = size(current_psms, 1)
    
    # Calculate expanded tolerance
    expanded_tolerance = collection_tolerance * expansion_factor
    
    # Create expanded model for collection
    # Keep the same bias, just expand the window
    expanded_model = MassErrorModel(
        getMassOffset(current_model),
        (expanded_tolerance, expanded_tolerance)
    )
    
    @info "Testing tolerance expansion for PSM collection:" *
          "\n  Original collection tolerance: ±$(round(collection_tolerance, digits=1)) ppm" *
          "\n  Expanded collection tolerance: ±$(round(expanded_tolerance, digits=1)) ppm" *
          "\n  Current PSM count: $current_psm_count"
    
    # Store original model
    original_model = getMassErrorModel(search_context, ms_file_idx)
    
    # Set expanded model for collection
    setMassErrorModel!(search_context, ms_file_idx, expanded_model)
    
    # Collect PSMs with expanded tolerance
    expanded_psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    expanded_psm_count = size(expanded_psms, 1)
    
    # Calculate improvement
    psm_increase = expanded_psm_count - current_psm_count
    improvement_ratio = current_psm_count > 0 ? psm_increase / current_psm_count : 0.0
    
    @info "Expansion results:" *
          "\n  Expanded PSM count: $expanded_psm_count" *
          "\n  PSM increase: $psm_increase ($(round(100*improvement_ratio, digits=1))%)"
    
    # Check if expansion was beneficial (any improvement)
    if expanded_psm_count <= current_psm_count
        # No improvement, restore original and return
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        @info "No improvement found, keeping original results"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Significant improvement found - refit model with expanded PSM set
    @info "Significant improvement found, refitting model with expanded PSM set"
    
    # Get matched fragments for the expanded PSM set
    fragments = get_matched_fragments(spectra, expanded_psms, search_context, params, ms_file_idx)
    
    if length(fragments) == 0
        # Failed to get fragments, restore original
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        @warn "Failed to extract fragments from expanded PSM set, keeping original"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Fit new mass error model from expanded PSM set
    refitted_model, refitted_ppm_errs = fit_mass_err_model(params, fragments)
    
    if refitted_model === nothing
        # Failed to fit model, restore original
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        @warn "Failed to fit mass error model from expanded PSM set, keeping original"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Success! Use the expanded results
    @info "Successfully expanded tolerance:" *
          "\n  Original fitted: ±$(round((getLeftTol(current_model) + getRightTol(current_model))/2, digits=1)) ppm" *
          "\n  Expanded fitted: ±$(round((getLeftTol(refitted_model) + getRightTol(refitted_model))/2, digits=1)) ppm" *
          "\n  PSM improvement: $psm_increase PSMs ($(round(100*improvement_ratio, digits=1))%)"
    
    # Update the model in search context
    setMassErrorModel!(search_context, ms_file_idx, refitted_model)
    
    return expanded_psms, refitted_model, refitted_ppm_errs, true
end
```

### 3. Integration with Convergence Check

Modify the convergence checking and storage process:

```julia
function check_and_store_convergence!(results, search_context, params, ms_file_idx, 
                                      psms, mass_err_model, ppm_errs, strategy_name::String,
                                      iteration_state::IterationState,  # NEW: Need iteration state
                                      filtered_spectra, spectra)  # NEW: Need for expansion
    if mass_err_model === nothing || psms === nothing
        return false
    end
    
    current_model = getMassErrorModel(search_context, ms_file_idx)
    
    if !check_convergence(psms, mass_err_model, current_model, ppm_errs, getMinPsms(params))
        return false
    end
    
    @info "Converged after $strategy_name"
    
    # NEW: Test tolerance expansion if enabled
    final_psms = psms
    final_model = mass_err_model
    final_ppm_errs = ppm_errs
    was_expanded = false
    
    if should_test_expansion(params)  # Check if expansion is enabled
        collection_tol = iteration_state.collection_tolerance
        
        if collection_tol > 0.0f0  # Ensure we have a valid collection tolerance
            final_psms, final_model, final_ppm_errs, was_expanded = test_tolerance_expansion!(
                results, search_context, params, ms_file_idx,
                psms, mass_err_model, ppm_errs,
                collection_tol, filtered_spectra, spectra
            )
            
            if was_expanded
                @info "Using expanded tolerance results for final model"
            end
        else
            @warn "Invalid collection tolerance, skipping expansion test"
        end
    end
    
    # Store final results (original or expanded)
    store_convergence_results!(
        results, search_context, ms_file_idx,
        final_psms, final_model, final_ppm_errs
    )
    
    return true
end
```

### 4. Parameter Configuration

#### Update ParameterTuningSearchParameters struct (types.jl)
```julia
struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # ... existing fields ...
    
    # Tolerance expansion settings
    enable_tolerance_expansion::Bool
    expansion_factor::Float32
    max_expanded_tolerance_factor::Float32  # Maximum allowed expansion relative to max_tol
    
    # ... rest of struct ...
end
```

#### Update constructor to extract new parameters
```julia
function ParameterTuningSearchParameters(params::PioneerParameters)
    # ... existing extraction ...
    
    # Extract tolerance expansion settings
    enable_expansion = if hasproperty(search_params, :tolerance_expansion) && 
                         hasproperty(search_params.tolerance_expansion, :enabled)
        Bool(search_params.tolerance_expansion.enabled)
    else
        true  # Default: enabled
    end
    
    expansion_factor = if hasproperty(search_params, :tolerance_expansion) && 
                         hasproperty(search_params.tolerance_expansion, :expansion_factor)
        Float32(search_params.tolerance_expansion.expansion_factor)
    else
        1.5f0  # Default: 50% expansion
    end
    
    # ... continue construction ...
end
```

#### Helper functions for parameter access
```julia
# Add to types.jl
should_test_expansion(params::ParameterTuningSearchParameters) = params.enable_tolerance_expansion
get_expansion_factor(params::ParameterTuningSearchParameters) = params.expansion_factor
```

### 5. JSON Configuration

Add to parameter JSON files:

```json
{
    "parameter_tuning": {
        "search_settings": {
            "tolerance_expansion": {
                "enabled": true,
                "expansion_factor": 1.5
            }
        }
    }
}
```

### 6. Update Main Process Flow

Modify `process_file!` to pass iteration_state through the pipeline:

```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    # ... initialization ...
    
    iteration_state = IterationState()
    
    # ... main iteration loop ...
    
    while continue_iteration
        # ... existing iteration logic ...
        
        # Execute strategy with iteration state
        psms, mass_err_model, ppm_errs = execute_strategy(
            current_strategy, filtered_spectra, spectra, 
            search_context, params, ms_file_idx, iteration_state  # Pass iteration state
        )
        
        # Check convergence with expansion capability
        if check_and_store_convergence!(
            results, search_context, params, ms_file_idx,
            psms, mass_err_model, ppm_errs, "iteration $total_iterations",
            iteration_state, filtered_spectra, spectra  # Pass additional parameters
        )
            converged = true
            break
        end
        
        # ... rest of iteration ...
    end
    
    # ... finalization ...
end
```

## Implementation Phases

### Phase 1: Core Infrastructure (Day 1)
1. Add `collection_tolerance` to `IterationState`
2. Implement tolerance tracking in `execute_strategy`
3. Create `test_tolerance_expansion!` function
4. Add helper functions for parameter access

### Phase 2: Integration (Day 2)
1. Modify `check_and_store_convergence!` to call expansion test
2. Update all callers to pass required parameters
3. Update `process_file!` to thread iteration_state through

### Phase 3: Configuration (Day 3)
1. Add new fields to `ParameterTuningSearchParameters`
2. Update parameter extraction in constructor
3. Add validation in `paramsChecks.jl`
4. Update default JSON configurations

### Phase 4: Testing & Refinement (Day 4-5)
1. Unit tests for `test_tolerance_expansion!`
2. Integration tests with known datasets
3. Performance benchmarking
4. Edge case handling

## Testing Strategy

### Unit Tests
```julia
@testset "Tolerance Expansion" begin
    # Test expansion with improvement
    @test test_expansion_with_improvement()
    
    # Test expansion without improvement
    @test test_expansion_no_improvement()
    
    # Test expansion with invalid model refit
    @test test_expansion_failed_refit()
    
    # Test disabled expansion
    @test test_expansion_disabled()
end
```

### Integration Tests
1. Run on E. coli test dataset
2. Verify PSM count improvements
3. Check downstream impact on other search methods
4. Validate final protein quantification

### Performance Tests
1. Measure runtime overhead (target: <10%)
2. Memory usage impact
3. Convergence rate changes

## Risk Analysis & Mitigation

### Risk 1: Infinite Expansion Loop
**Mitigation**: 
- Only test expansion once after convergence
- Cap maximum tolerance with `max_tolerance_factor`
- No recursive expansion

### Risk 2: Quality Degradation
**Mitigation**:
- Require minimum improvement ratio
- Validate refitted model parameters
- Maintain all existing quality filters

### Risk 3: Breaking Changes
**Mitigation**:
- Default to enabled but make configurable
- Backward compatible parameter structure
- Graceful fallback on any failure

## Success Criteria

1. **PSM Yield**: Any increase in valid PSMs for suitable datasets
2. **Quality**: FDR remains controlled at configured level
3. **Performance**: <10% increase in parameter tuning time
4. **Robustness**: No decrease in convergence success rate
5. **Compatibility**: All existing tests pass

## Rollback Plan

If issues arise:
1. Set `"enabled": false` in JSON configuration
2. Code paths automatically skip expansion
3. Behavior reverts to current implementation
4. No code changes required for rollback

## Documentation Updates

1. Update PARAMETER_TUNING_IMPROVEMENTS_SUMMARY.md
2. Add expansion parameters to CLAUDE.md
3. Update example JSON files with comments
4. Add diagnostic logging for troubleshooting

## Summary

This detailed plan implements a post-convergence tolerance expansion test that:
- Uses the collection tolerance (not fitted tolerance) as the basis for expansion
- Tests a 50% expansion to discover potentially missed PSMs
- Refits the model if significant improvement is found
- Is fully configurable and can be disabled if needed
- Maintains all quality controls and safety checks
- Integrates cleanly with the existing convergence framework

The implementation is designed to be safe, backward-compatible, and provide meaningful improvements in PSM discovery without sacrificing quality.