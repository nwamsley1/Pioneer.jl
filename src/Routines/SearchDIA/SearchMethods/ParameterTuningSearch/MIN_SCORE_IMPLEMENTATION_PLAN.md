# Multi-Score Search Implementation Plan for ParameterTuningSearch

## Executive Summary

This document outlines the implementation plan for enhancing ParameterTuningSearch with a multi-score search capability. Currently, the algorithm uses a single `min_score` parameter (default: 22) to filter fragment index matches. This enhancement will allow trying multiple score thresholds in sequence, improving the robustness of parameter tuning for challenging datasets.

## Problem Statement

The current implementation uses a fixed `min_score` threshold for filtering precursor matches in the fragment index search. When this threshold is too high (e.g., 22), some datasets may fail to produce enough PSMs for convergence. The proposed solution is to try multiple score thresholds systematically, starting with stricter thresholds and relaxing them if needed.

## Proposed Solution

### Configuration Change

**Current JSON structure:**
```json
{
  "parameter_tuning": {
    "fragment_index_params": {
      "min_score": 22
    }
  }
}
```

**Proposed JSON structure:**
```json
{
  "parameter_tuning": {
    "fragment_index_params": {
      "min_score": [22, 17, 15]  // Now accepts array or single value
    }
  }
}
```

### Algorithm Enhancement

Add a score threshold loop INSIDE each phase of the 3-phase loop. Each phase will try all score thresholds before moving to the next phase. This allows each bias setting (Phase 1: zero, Phase 2: positive, Phase 3: negative) to explore different score thresholds independently.

## Detailed Implementation Plan

### 1. Type System Updates

**File: `types.jl`**

```julia
# Modify ParameterTuningSearchParameters struct
struct ParameterTuningSearchParameters <: SearchParameters
    # ... existing fields ...
    min_index_search_scores::Vector{UInt8}  # Changed from min_index_search_score::UInt8
    # ... rest of fields ...
end

# Update constructor to handle both single value and array
function ParameterTuningSearchParameters(params::Dict{String, Any})
    # ... existing code ...
    
    # Handle min_score as either single value or array
    min_score_raw = params["parameter_tuning"]["fragment_index_params"]["min_score"]
    min_scores = if min_score_raw isa Vector
        UInt8.(min_score_raw)
    else
        [UInt8(min_score_raw)]  # Convert single value to array
    end
    
    # ... rest of constructor ...
end
```

### 2. Main Algorithm Modification

**Modified `run_single_phase` function:**

```julia
function run_single_phase(
    phase::Int64,
    filtered_spectra::FilteredMassSpecData,
    iteration_state::IterationState,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    @info "="^40
    @info "Starting Phase $phase (Bias: $(calculate_phase_bias_shift(phase, params)) ppm)"
    @info "="^40
    
    # Get all score thresholds to try
    score_thresholds = getMinIndexSearchScores(params)  # Returns Vector{UInt8}
    
    # NEW: Score threshold loop INSIDE each phase
    for (score_idx, min_score) in enumerate(score_thresholds)
        @info "  Phase $phase - Trying min_score = $min_score (threshold $score_idx of $(length(score_thresholds)))"
        
        # Update current score threshold
        setCurrentMinScore!(params, min_score)
        
        # Reset for new phase with appropriate bias shift
        reset_for_new_phase!(search_context, ms_file_idx, params, phase, iteration_state)
        
        # Step 0: Try at initial tolerance with phase bias
        psms_initial, psm_count = collect_and_log_psms(
            filtered_spectra, spectra, search_context,
            params, ms_file_idx, "at initial tolerance with Phase $phase bias and min_score=$min_score"
        )
        
        # Try to fit models and check convergence
        mass_err_model, ppm_errs, _ = fit_models_from_psms(
            psms_initial, spectra, search_context, params, ms_file_idx
        )
        
        if check_and_store_convergence!(
            results, search_context, params, ms_file_idx,
            psms_initial, mass_err_model, ppm_errs,
            "Phase $phase initial attempt with min_score=$min_score",
            iteration_state, filtered_spectra, spectra
        )
            @info "  ✓ Converged in Phase $phase with min_score=$min_score (initial tolerance)"
            return true
        end
        
        # Iteration loop with tolerance expansion
        iterations_per_phase = getIterationsPerPhase(getIterationSettings(params))
        
        for iter in 1:iterations_per_phase
            @info "  Phase $phase, Score $min_score, Iteration $iter"
            iteration_state.total_iterations += 1
            
            # Step 1: EXPAND mass tolerance
            expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                  getMassToleranceScaleFactor(getIterationSettings(params)))
            
            # Step 2: COLLECT PSMs with expanded tolerance
            psms_expanded, _ = collect_and_log_psms(
                filtered_spectra, spectra, search_context,
                params, ms_file_idx, "with expanded tolerance and min_score=$min_score"
            )
            
            # Step 3: FIT MODEL to determine optimal bias
            mass_err_model_for_bias, _, _ = fit_models_from_psms(
                psms_expanded, spectra, search_context, params, ms_file_idx
            )
            
            # Step 4: ADJUST BIAS based on fitted model
            if mass_err_model_for_bias !== nothing
                new_bias = getMassOffset(mass_err_model_for_bias)
                current_model = getMassErrorModel(search_context, ms_file_idx)
                old_bias = getMassOffset(current_model)
                
                # Update the bias while keeping the expanded tolerance
                adjusted_model = MassErrorModel(
                    new_bias,
                    (getLeftTol(current_model), getRightTol(current_model))
                )
                setMassErrorModel!(search_context, ms_file_idx, adjusted_model)
                
                @info "    Adjusted bias: $(round(old_bias, digits=2)) → $(round(new_bias, digits=2)) ppm"
            end
            
            # Step 5: COLLECT PSMs again with adjusted bias
            psms_adjusted, _ = collect_and_log_psms(
                filtered_spectra, spectra, search_context,
                params, ms_file_idx, "with adjusted bias and min_score=$min_score"
            )
            
            # Step 6: FIT FINAL MODELS
            mass_err_model, ppm_errs, psm_count = fit_models_from_psms(
                psms_adjusted, spectra, search_context, params, ms_file_idx
            )
            
            # Step 7: CHECK CONVERGENCE
            if check_and_store_convergence!(
                results, search_context, params, ms_file_idx,
                psms_adjusted, mass_err_model, ppm_errs,
                "Phase $phase, Iteration $iter with min_score=$min_score",
                iteration_state, filtered_spectra, spectra
            )
                @info "  ✓ Converged in Phase $phase with min_score=$min_score (iteration $iter)"
                return true
            end
        end
        
        @info "  Phase $phase with min_score=$min_score did not converge"
    end
    
    @info "Phase $phase completed without convergence (tried all score thresholds)"
    return false
end
```

### 3. Supporting Functions

```julia
# New getter for multiple scores
getMinIndexSearchScores(params::ParameterTuningSearchParameters) = params.min_index_search_scores

# Modify existing getter to return current active score
getMinIndexSearchScore(params::ParameterTuningSearchParameters) = params.current_min_score

# New setter for current score
function setCurrentMinScore!(params::ParameterTuningSearchParameters, score::UInt8)
    params.current_min_score = score
end
```

## Process Flow Schematic

```
╔════════════════════════════════════════════════════════════════════╗
║                    ENHANCED PARAMETERTUNINGSEARCH                  ║
║                 Multi-Score Mass Error & RT Calibration            ║
╚════════════════════════════════════════════════════════════════════╝
                                  │
                                  ▼
                      ┌──────────────────────┐
                      │  Initialize Models   │
                      │  - Mass Error: ±5ppm │
                      │  - Quad: Gaussian    │
                      └──────────────────────┘
                                  │
                                  ▼
        ┌────────────────────────────────────────────┐
        │      SCAN SCALING LOOP (Outer Loop)       │
        │  Initial: 500 → Scale 10x → Max: 80k      │
        └────────────────────────────────────────────┘
                                  │
                                  ▼
            ┌────────────────────────────────┐
            │  3-PHASE LOOP (Middle Loop)    │
            │  Phase 1: Bias = 0 ppm         │
            │  Phase 2: Bias = +max_tol     │
            │  Phase 3: Bias = -max_tol     │
            └────────────────────────────────┘
                                  │
                                  ▼
     ╔═══════════════════════════════════════════════════════╗
     ║    NEW: SCORE THRESHOLD LOOP (Inside Each Phase)      ║
     ║     Try each min_score in [22, 17, 15] sequence      ║
     ╚═══════════════════════════════════════════════════════╝
                                  │
                    ┌─────────────┼─────────────┐
                    ▼             ▼             ▼
              [Score: 22]   [Score: 17]   [Score: 15]
                    │             │             │
                    └─────────────┼─────────────┘
                                  ▼
                  [If not converged with any score,
                   move to next phase]
                                  │
                                  ▼
          ┌──────────────────────────────────┐
          │  ITERATION LOOP (Innermost)      │
          │  Up to 3 iterations per score    │
          │  - Expand tolerance              │
          │  - Detect & adjust bias          │
          │  - Check convergence             │
          └──────────────────────────────────┘
                                  │
                                  ▼
                        [Mass Error Estimation
                         Process (unchanged)]
```

## Detailed Execution Flow

```
SCAN ATTEMPT 1 (500 scans)
│
├─ PHASE 1 (Bias = 0 ppm)
│  ├─ Score: 22
│  │  ├─ Initial attempt → No convergence
│  │  ├─ Iteration 1: Expand tol → Collect → Adjust bias → No convergence
│  │  ├─ Iteration 2: Expand tol → Collect → Adjust bias → No convergence
│  │  └─ Iteration 3: Expand tol → Collect → Adjust bias → No convergence
│  ├─ Score: 17
│  │  ├─ Initial attempt → No convergence
│  │  ├─ Iteration 1: Expand tol → Collect → Adjust bias → No convergence
│  │  ├─ Iteration 2: Expand tol → Collect → Adjust bias → CONVERGED! ✓
│  │  └─ [Iteration 3 skipped - already converged]
│  └─ [Score: 15 skipped - already converged]
│
├─ [PHASE 2 skipped - already converged]
└─ [PHASE 3 skipped - already converged]

RESULT: Converged in Phase 1, Score 17, Iteration 2, with 500 scans
```

## Execution Flow Example

For `min_score: [22, 17, 15]`:

```
START
├─ Scan Count: 500
│  ├─ Phase 1 (bias=0)
│  │  ├─ Score 22: Try iterations 1-3 → No convergence
│  │  ├─ Score 17: Try iterations 1-3 → No convergence  
│  │  └─ Score 15: Try iterations 1-3 → No convergence
│  ├─ Phase 2 (bias=+max_tol)
│  │  ├─ Score 22: Try iterations 1-3 → No convergence
│  │  ├─ Score 17: Initial attempt → CONVERGED! ✓
│  │  └─ [Score 15 skipped - already converged]
│  └─ [Phase 3 skipped - already converged]
│
└─ [Further scan scaling skipped - already converged]

RESULT: Converged with Phase 2, min_score=17, 500 scans
```

## Benefits

1. **Robustness**: Each phase independently explores score thresholds with its specific bias setting
2. **Efficiency**: Tries stricter thresholds first within each phase, only relaxing if necessary
3. **Backwards Compatible**: Single value min_score still works (converted to single-element array)
4. **Transparency**: Clear logging shows which phase and score threshold achieved convergence
5. **Optimal Quality**: Prefers higher scores when possible, ensuring better PSM quality
6. **Phase-Specific Optimization**: Each bias setting (zero, positive, negative) gets its own chance to find the best score threshold

## Implementation Steps

1. **Phase 1**: Update type definitions
   - Modify `ParameterTuningSearchParameters` struct to store `min_index_search_scores::Vector{UInt8}`
   - Add field for `current_min_score::UInt8` to track active score
   - Update constructor to handle both array and single value inputs
   - Add getter/setter functions for score management

2. **Phase 2**: Modify `run_single_phase` function
   - Add score threshold loop inside each phase (after phase setup, before iterations)
   - Update `reset_for_new_phase!` calls to occur for each score threshold
   - Ensure proper logging to show phase, score, and iteration clearly
   - Return immediately upon convergence with any score

3. **Phase 3**: Update parameter passing
   - Ensure current `min_score` is properly passed to `filterPrecursorMatches!` 
   - Modify `library_search` to use `getMinIndexSearchScore(params)` for current score
   - Update all relevant function signatures to support dynamic score

4. **Phase 4**: Testing
   - Test with single value (backwards compatibility)
   - Test with array of values `[22, 17, 15]`
   - Verify convergence happens at expected phase/score combinations
   - Check logging output for clarity and debugging ease

## Risk Assessment

- **Low Risk**: Changes are localized to ParameterTuningSearch module
- **Medium Complexity**: Adds score loop inside phase loop, increasing total iterations
- **Performance Impact**: Maximum iterations = scan_attempts × 3 phases × N scores × 3 iterations
- **Testing Required**: Need to verify behavior with challenging datasets

## Alternative Approaches Considered

1. **Score Loop Outside Phase Loop** (Original proposal)
   - Rejected: Would try all scores with Phase 1 bias before trying Phase 2, less flexible

2. **Parallel Score Search**: Try all scores in parallel
   - Rejected: Would use N× memory and complicate result merging

3. **Adaptive Score Reduction**: Dynamically reduce score based on PSM count
   - Rejected: Less predictable, harder to debug and reproduce

4. **Score in Iteration Loop**: Put score loop inside iteration loop
   - Rejected: Would reset tolerance expansion for each score, inefficient

## Conclusion

This implementation provides a systematic way to handle varying data quality by trying multiple score thresholds within each phase. The approach:
- Gives each phase (bias setting) equal opportunity to find convergence with different scores
- Maintains clear separation of concerns (bias → score → tolerance expansion)
- Is backwards compatible and minimally invasive
- Provides better convergence for challenging datasets while maintaining high quality for good data

**Approval Request**: Please review this revised implementation plan where the score loop is placed inside each phase, and provide feedback or approval to proceed.