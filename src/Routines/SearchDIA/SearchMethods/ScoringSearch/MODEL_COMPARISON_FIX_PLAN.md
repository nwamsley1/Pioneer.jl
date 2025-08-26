# Model Comparison File Modification Fix Plan

## Problem Statement

The current model comparison implementation in `select_psm_scoring_model` has critical issues:

1. **File Modification During Testing**: Each candidate model writes its scores back to Arrow files during evaluation
2. **File Restoration Overhead**: After testing each model, we restore original files from memory, causing:
   - Bus errors due to excessive memory operations
   - Memory pressure from storing all original DataFrames
   - Inefficient I/O with multiple read/write cycles
3. **Race Conditions**: Potential file corruption from rapid write/restore cycles

### Current Problematic Flow
```
For each model candidate:
  1. Train model (writes to Arrow files)
  2. Read files to count passing targets
  3. Restore original files from memory (causes bus error)
  4. Repeat for next model
```

## Root Cause Analysis

The file writing happens at the end of the scoring functions:
- `sort_of_percolator_in_memory!` (lines 123-127 in percolatorSortOf.jl): Drops columns then writes to Arrow files
- `train_probit_model_in_memory` (lines 390-394): Similar pattern of writing results to files

The key insight is that the file writing is already separated from the scoring logic - it happens at the very end. We can refactor this more cleanly by **separating concerns**: one function scores, another function writes.

## Solution Architecture

### Design Principles
1. **Separation of Concerns**: Separate model training/scoring from file I/O completely
2. **Single Write Principle**: Only write files once with the final selected model
3. **Minimal Code Duplication**: Refactor existing functions rather than creating copies
4. **Clean Abstractions**: Clear distinction between scoring and I/O operations

### Proposed Flow
```
Model Selection Phase (no file writes):
  1. Load PSMs into memory once
  2. For each model candidate:
     - Train model on in-memory data
     - Score PSMs in-memory
     - Evaluate performance on DataFrame
  3. Select best model

Application Phase (single file write):
  4. Train selected model (or reuse if already trained)
  5. Write results to Arrow files once using new utility function
```

## Implementation Plan

### Phase 1: Refactor Core ML Functions to Separate Scoring from Writing

#### 1.1 Refactor `sort_of_percolator_in_memory!` in `percolatorSortOf.jl`

**Current structure:**
```julia
function sort_of_percolator_in_memory!(psms, file_paths, ...)
    # ... scoring logic ...
    dropVectorColumns!(psms)
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    return models
end
```

**Refactored structure:**
```julia
# Modified original function - no longer writes files, no longer needs file_paths
function sort_of_percolator_in_memory!(psms::DataFrame, features::Vector{Symbol}, ...)
    # ... all scoring logic remains unchanged ...
    # Remove lines 123-127 (dropVectorColumns and writing loop)
    # Remove file_paths parameter from signature
    return models
end

# New utility function for writing scored PSMs (export from percolatorSortOf.jl)
function write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})
    dropVectorColumns!(psms)
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
end
```

#### 1.2 Update `train_xgboost_model_in_memory` in `score_psms.jl`

**Current structure:**
```julia
function train_xgboost_model_in_memory(best_psms, file_paths, ...)
    # ... setup ...
    return sort_of_percolator_in_memory!(best_psms, file_paths, features, ...)
end
```

**Refactored structure:**
```julia
function train_xgboost_model_in_memory(best_psms, file_paths, ...)
    # ... setup ...
    # Note: file_paths is kept in signature for compatibility but not used
    models = sort_of_percolator_in_memory!(best_psms, features, ...)  # No file_paths passed
    # NO FILE WRITING HERE - will be done at higher level
    return models
end
```

#### 1.3 Update `train_probit_model_in_memory` in `score_psms.jl`

**Current structure:**
```julia
function train_probit_model_in_memory(best_psms, file_paths, ...)
    # ... probit regression ...
    dropVectorColumns!(best_psms)
    for (ms_file_idx, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    return nothing
end
```

**Refactored structure:**
```julia
function train_probit_model_in_memory(best_psms, file_paths, ...)
    # ... probit regression logic remains unchanged ...
    # REMOVE lines 390-394 (dropVectorColumns and writing loop)
    # NO FILE WRITING HERE - will be done at higher level
    return nothing
end
```

#### 1.4 NO CHANGES to `score_precursor_isotope_traces_in_memory`

The function remains as-is. File writing will be controlled at a higher level in the call chain.

### Phase 2: Modify Model Comparison Logic

#### 2.1 Update `select_psm_scoring_model`

**Remove these sections:**
```julia
# Lines 156-158: File backup logic - DELETE ENTIRELY
original_file_contents = Dict{String, DataFrame}()
for fpath in file_paths
    if endswith(fpath, ".arrow")
        original_file_contents[fpath] = DataFrame(Arrow.Table(fpath))
    end
end

# Lines 187-189: File restoration logic - DELETE ENTIRELY
for (fpath, original_df) in original_file_contents
    writeArrow(fpath, original_df)
end
```

**Keep the model training call unchanged** (lines 167-171):
```julia
# This stays the same - it already doesn't write files after our Phase 1 changes
score_precursor_isotope_traces_in_memory(
    psms_copy, file_paths, precursors, config,
    match_between_runs, max_q_value_xgboost_rescore,
    max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
)
```

**Replace target counting** (line 174):
```julia
# Old (reads from files):
target_count = count_passing_targets(file_paths, qvalue_threshold)

# New (reads from DataFrame):
target_count = count_passing_targets(psms_copy, qvalue_threshold)
```

#### 2.2 Modify existing `count_passing_targets` function

Change the function to accept a DataFrame instead of file paths:

```julia
# Old signature:
function count_passing_targets(file_paths::Vector{String}, qvalue_threshold::Float64)
    # ... reads from files ...
end

# New signature (replaces the old one):
function count_passing_targets(scored_psms::DataFrame, qvalue_threshold::Float64)
    # Use existing q-value calculation logic
    if :prob in propertynames(scored_psms) && :target in propertynames(scored_psms)
        # Calculate q-values from probabilities
        qvals = Vector{Float32}(undef, nrow(scored_psms))
        get_qvalues!(scored_psms.prob, scored_psms.target, qvals)
        
        # Count targets passing threshold
        return sum((scored_psms.target .== true) .& (qvals .<= qvalue_threshold))
    else
        error("DataFrame must contain :prob and :target columns")
    end
end
```

### Phase 3: Handle File Writing After Model Selection

#### 3.1 Modify `score_precursor_isotope_traces` to write files AFTER model selection

**Current flow in Case 3 (lines 99-117):**
```julia
else
    # Case 3: In-memory with automatic model comparison (<100K)
    @user_info "Running automatic model comparison for $psms_count PSMs (< 100K)"
    model_config = select_psm_scoring_model(
        best_psms, file_paths, precursors, match_between_runs,
        max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
        min_PEP_neg_threshold_xgboost_rescore
    )
end

@user_info "Selected PSM scoring model: $(model_config.name)"

# Execute selected model using unified function
models = score_precursor_isotope_traces_in_memory(
    best_psms, file_paths, precursors, model_config,
    match_between_runs, max_q_value_xgboost_rescore,
    max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
)
```

**Updated flow - add explicit file writing:**
```julia
else
    # Case 3: In-memory with automatic model comparison (<100K)
    @user_info "Running automatic model comparison for $psms_count PSMs (< 100K)"
    model_config = select_psm_scoring_model(
        best_psms, file_paths, precursors, match_between_runs,
        max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
        min_PEP_neg_threshold_xgboost_rescore
    )
end

@user_info "Selected PSM scoring model: $(model_config.name)"

# Execute selected model using unified function (no file writing due to Phase 1 changes)
models = score_precursor_isotope_traces_in_memory(
    best_psms, file_paths, precursors, model_config,
    match_between_runs, max_q_value_xgboost_rescore,
    max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
)

# NEW: Explicitly write the scored PSMs to files
write_scored_psms_to_files!(best_psms, file_paths)
```

#### 3.2 Also update Case 2 for consistency

**Case 2 (lines 94-108) also needs the explicit write:**
```julia
if psms_count >= MAX_FOR_MODEL_SELECTION  # 100K
    # Case 2: In-memory with default/advanced XGBoost (no comparison)
    @user_info "Using in-memory advanced XGBoost for $psms_count PSMs (< $max_psms_in_memory but ≥ 100K)"
    model_config = create_default_advanced_xgboost_config()
end
# ... existing code ...

# Execute selected model (no file writing due to Phase 1 changes)
models = score_precursor_isotope_traces_in_memory(
    best_psms, file_paths, precursors, model_config,
    match_between_runs, max_q_value_xgboost_rescore,
    max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
)

# NEW: Explicitly write the scored PSMs to files
write_scored_psms_to_files!(best_psms, file_paths)
```

### Phase 4: Out-of-Memory Case - No Changes Needed

#### 4.1 Leave Case 1 (out-of-memory) unchanged

The out-of-memory path uses `score_precursor_isotope_traces_out_of_memory!` which calls `sort_of_percolator_out_of_memory!`. This path:
- Does NOT do model comparison
- Can continue to write files directly
- No refactoring needed for this case

**Rationale**: The out-of-memory path processes files incrementally and needs to write them as it goes. Since there's no model comparison in this path, there's no need to separate scoring from writing.

### Phase 5: Function Access

No imports or exports needed. The `write_scored_psms_to_files!` function will be available through Pioneer's existing `importScripts()` mechanism which loads all functions.

## Summary of Complete Solution

### Files to Modify

1. **`src/utils/ML/percolatorSortOf.jl`**:
   - Remove lines 123-127 from `sort_of_percolator_in_memory!`
   - Remove `file_paths` parameter from function signature
   - Add new `write_scored_psms_to_files!` function

2. **`src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`**:
   - Update `train_xgboost_model_in_memory` to not pass `file_paths` to `sort_of_percolator_in_memory!`
   - Remove lines 390-394 from `train_probit_model_in_memory`
   - Update `select_psm_scoring_model`:
     - Remove file backup logic (lines 156-158)
     - Remove file restore logic (lines 187-189)
     - Update `count_passing_targets` call to pass DataFrame instead of file_paths
   - Modify existing `count_passing_targets` function to accept DataFrame instead of file paths
   - Update `score_precursor_isotope_traces` to call `write_scored_psms_to_files!` after model training

3. **No changes to**:
   - `score_precursor_isotope_traces_out_of_memory!` (Case 1)
   - `sort_of_percolator_out_of_memory!`
   - Any import/export statements

### Key Benefits

1. **Eliminates Bus Errors**: No more storing/restoring entire file contents in memory
2. **Faster Model Comparison**: No file I/O during model testing
3. **Cleaner Code**: Clear separation between scoring logic and file I/O
4. **Single Write**: Files only written once with the selected model's scores
5. **Minimal Changes**: Surgical refactoring that doesn't change the overall architecture

## Implementation Order

1. First, modify `percolatorSortOf.jl` to separate scoring from writing
2. Then update `score_psms.jl` to use the new structure
3. Test with SearchDIA on test data to verify correctness

## Risk Assessment

### Potential Issues

1. **DataFrame Size**: Keeping scored DataFrames in memory
   - Mitigation: Only during comparison phase, cleared after
   
2. **Column Compatibility**: Ensuring all required columns present
   - Mitigation: Validation checks, comprehensive testing

3. **Score Consistency**: Ensuring no-write variants produce same scores
   - Mitigation: Unit tests comparing outputs

### Rollback Plan

If issues arise:
1. Keep original functions intact during development
2. Add feature flag to toggle new behavior
3. Gradual rollout with monitoring

## Implementation Timeline

1. **Phase 1** (2-3 hours): Create non-writing functions
2. **Phase 2** (1-2 hours): Update model comparison logic
3. **Phase 3** (30 min): Validate main entry point
4. **Phase 4** (1 hour): Add safety checks and logging
5. **Testing** (2-3 hours): Comprehensive test suite
6. **Documentation** (30 min): Update code documentation

Total estimated time: 7-10 hours

## Success Criteria

1. ✅ No file modifications during model comparison
2. ✅ Single file write with selected model
3. ✅ No bus errors or memory issues
4. ✅ Faster model comparison execution
5. ✅ Same final results as current implementation
6. ✅ Clean, maintainable code structure

## Future Improvements

1. **Parallel Model Evaluation**: Test models concurrently since no file conflicts
2. **Caching**: Cache model training for repeated experiments
3. **Streaming Evaluation**: Process files in chunks for very large datasets
4. **Model Persistence**: Save trained models for reuse

## Code Organization

New files to create:
- `score_psms_no_write.jl` - Non-writing training functions
- `test/UnitTests/test_model_comparison_no_write.jl` - Test suite

Modified files:
- `score_psms.jl` - Update `select_psm_scoring_model`
- `percolatorSortOf.jl` - Add no-write variant

## Notes

- The no-write functions should be internal implementation details
- Public API remains unchanged
- Consider deprecating file backup/restore pattern elsewhere in codebase
- This pattern (separate evaluation from application) could benefit other search methods