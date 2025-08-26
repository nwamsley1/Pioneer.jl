# Plan to Fix Excessive Progress Bars During Model Selection

## Problem Statement
During ScoringSearch model selection phase, multiple progress bars appear even though the code attempts to suppress them:
- Model comparison tests 4 models (SimpleXGBoost, AdvancedXGBoost, ProbitRegression, SuperSimplified)
- Each model training shows ~8 progress bars (1 for training + 7 for file sorting)
- Total: ~32 progress bars during model selection
- Only the final model training should show progress bars

## Root Causes

### 1. Incomplete Output Suppression
**Current Code** (`score_psms.jl` line 160-167):
```julia
redirect_stdout(devnull) do
    score_precursor_isotope_traces_in_memory(
        psms_copy, file_paths, precursors, config,
        match_between_runs, max_q_value_xgboost_rescore,
        max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore,
        false  # show_progress = false during comparison
    )
end
```

**Issue**: Only redirects stdout, but ProgressBars write to stderr

### 2. File Sorting Always Shows Progress
**Current Code** (`ArrowOperations.jl` line 121, 128):
```julia
Threads.@threads for ref in ProgressBar(refs)  # Always shows
for ref in ProgressBar(refs)                   # Always shows
```

**Issue**: No `show_progress` parameter to control visibility

### 3. Progress Parameter Not Fully Propagated
- `show_progress=false` passed to training function
- But `sort_file_by_keys!` doesn't accept/use this parameter
- Called 7+ times during each model's post-processing

## Proposed Solution

### Phase 1: Fix Output Redirection During Model Selection

#### A. Update `score_psms.jl` (lines 159-167)
**Change from:**
```julia
redirect_stdout(devnull) do
    score_precursor_isotope_traces_in_memory(...)
end
```

**Change to:**
```julia
# Redirect both stdout AND stderr during model comparison
original_stdout = stdout
original_stderr = stderr
try
    redirect_stdout(devnull)
    redirect_stderr(devnull)
    
    score_precursor_isotope_traces_in_memory(
        psms_copy, file_paths, precursors, config,
        match_between_runs, max_q_value_xgboost_rescore,
        max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore,
        false  # show_progress = false
    )
finally
    redirect_stdout(original_stdout)
    redirect_stderr(original_stderr)
end
```

### Phase 2: Add Progress Control to File Operations

#### A. Update `ArrowOperations.jl` Function Signatures

1. **Modify `sort_file_by_keys!` for single file** (line ~60):
```julia
function sort_file_by_keys!(ref::FileReference, keys::Symbol...; 
                           reverse::Union{Bool, Vector{Bool}}=false,
                           show_progress::Bool=true)  # ADD THIS
```

2. **Modify `sort_file_by_keys!` for multiple files** (line ~113):
```julia
function sort_file_by_keys!(refs::Vector{<:FileReference}, keys::Symbol...; 
                           reverse::Union{Bool, Vector{Bool}}=false, 
                           parallel::Bool=true,
                           show_progress::Bool=true)  # ADD THIS
```

#### B. Update Progress Bar Logic (lines 121-130)
**Change from:**
```julia
if parallel && length(refs) > 1
    Threads.@threads for ref in ProgressBar(refs)
        if exists(ref)
            sort_file_by_keys!(ref, keys...; reverse=reverse)
        end
    end
else
    for ref in ProgressBar(refs)
        if exists(ref)
            sort_file_by_keys!(ref, keys...; reverse=reverse)
```

**Change to:**
```julia
if parallel && length(refs) > 1
    if show_progress
        Threads.@threads for ref in ProgressBar(refs)
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
            end
        end
    else
        Threads.@threads for ref in refs  # No ProgressBar
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
            end
        end
    end
else
    if show_progress
        for ref in ProgressBar(refs)
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
            end
        end
    else
        for ref in refs  # No ProgressBar
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
            end
        end
    end
end
```

### Phase 3: Thread Progress Through Call Chain

#### A. Update `score_precursor_isotope_traces_in_memory` signature
In `score_psms.jl` (line ~273):
```julia
function score_precursor_isotope_traces_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    show_progress::Bool = true
)
```

#### B. Pass `show_progress` to Child Functions
Update calls to pass `show_progress`:
- `train_xgboost_model_in_memory` 
- `train_probit_model_in_memory`
- Any `sort_file_by_keys!` calls within these functions

#### C. Update `ScoringSearch.jl` Sorting Calls
For all 6 calls to `sort_file_by_keys!` in `ScoringSearch.jl`, determine context:
- During model selection: Add `show_progress=false`
- During final processing: Keep default `show_progress=true`

Since we can't easily determine context in ScoringSearch.jl, we'll need to:
1. Add a `show_progress` field to `ScoringSearchResults` or
2. Use a thread-local variable to track if we're in model selection

**Simpler approach**: Add parameter to functions that call sorting:
```julia
# In summarize_results! and any function that calls sort_file_by_keys!
function process_scoring_step!(refs, ...; show_progress::Bool=true)
    sort_file_by_keys!(refs, :prob, :target; 
                      reverse=[true, true], 
                      show_progress=show_progress)
end
```

### Phase 4: Handle Progress in ML Training

#### A. Check `percolatorSortOf.jl`
Ensure that `show_progress` parameter is properly handled:
- Line 52: `pbar = show_progress ? ProgressBar(...) : nothing`
- Line 88: `show_progress && update(pbar)`

This appears to be already correct, but verify it works when `show_progress=false`.

### Phase 5: Testing Strategy

1. **Test Model Selection with Suppressed Output**:
   ```julia
   # Should show NO progress bars during selection
   # Only show progress for final model
   ```

2. **Test Direct Training (no selection)**:
   ```julia
   # Should show normal progress bars
   ```

3. **Test File Operations**:
   ```julia
   # Test sort_file_by_keys! with show_progress=false
   # Test sort_file_by_keys! with show_progress=true
   ```

## Implementation Order

1. **Fix stderr redirection** (Phase 1) - Immediate improvement
2. **Add show_progress to ArrowOperations** (Phase 2) - Core functionality
3. **Thread parameter through scoring functions** (Phase 3) - Complete integration
4. **Verify ML training respects parameter** (Phase 4) - Ensure consistency
5. **Test all scenarios** (Phase 5) - Validate solution

## Expected Results

### Before Fix:
```
[ Info: Model comparison: Testing 4 models on 50000 PSMs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
[32 progress bars appear here during model selection]
  SimpleXGBoost: 1234 IDs at q ≤ 0.01
  AdvancedXGBoost: 1456 IDs at q ≤ 0.01
  ProbitRegression: 1100 IDs at q ≤ 0.01  
  SuperSimplified: 900 IDs at q ≤ 0.01
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✓ Selected: AdvancedXGBoost (1456 IDs)
[ Info: Training final model: AdvancedXGBoost
[8 progress bars for final training]
```

### After Fix:
```
[ Info: Model comparison: Testing 4 models on 50000 PSMs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  SimpleXGBoost: 1234 IDs at q ≤ 0.01
  AdvancedXGBoost: 1456 IDs at q ≤ 0.01
  ProbitRegression: 1100 IDs at q ≤ 0.01
  SuperSimplified: 900 IDs at q ≤ 0.01
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✓ Selected: AdvancedXGBoost (1456 IDs)
[ Info: Training final model: AdvancedXGBoost
[8 progress bars for final training - these remain visible]
```

## Benefits

1. **Cleaner Output**: 75% reduction in progress bars during model selection
2. **Faster Perception**: Model comparison appears faster without visual clutter
3. **Clear Distinction**: User sees comparison results, then final training progress
4. **Backward Compatible**: Default behavior unchanged (progress shown)
5. **Flexible Control**: Can disable all progress if needed

## Risks and Mitigation

### Risk 1: Breaking Existing Code
- **Mitigation**: Use default parameter values (show_progress=true)
- All existing code continues to work unchanged

### Risk 2: Stderr Redirection Side Effects  
- **Mitigation**: Use try/finally to ensure streams restored
- Only redirect during model comparison phase

### Risk 3: Thread Safety with Redirects
- **Mitigation**: Redirection happens in single-threaded context
- Only during model selection, not during parallel operations

## Files to Modify

1. `/src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl` - Fix stderr redirect
2. `/src/utils/FileOperations/io/ArrowOperations.jl` - Add show_progress parameter
3. `/src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl` - Pass parameter through

## Estimated Time: 45 minutes

- Phase 1: 10 minutes (stderr redirect)
- Phase 2: 15 minutes (ArrowOperations update)
- Phase 3: 10 minutes (parameter threading)
- Phase 4: 5 minutes (verification)
- Phase 5: 5 minutes (testing)