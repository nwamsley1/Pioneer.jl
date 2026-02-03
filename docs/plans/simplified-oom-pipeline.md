# Plan: Simplified OOM Percolator Pipeline

## Current Problem

The current OOM implementation (`sort_of_percolator_out_of_memory!`) duplicates most of the logic from `sort_of_percolator_in_memory!`, leading to:
- ~400 lines of duplicate code
- Multiple subtle bugs (pairing, MBR_num_runs, transfer candidates)
- Difficult maintenance

## Proposed Solution

**Key Insight**: The only reason we need OOM is that we can't load ALL PSMs into memory. But we CAN:
1. Load a SAMPLE of PSMs for training
2. Apply the trained models to files one-at-a-time

**Approach**:
- Call the exact same `sort_of_percolator_in_memory!` on sampled PSMs
- Add ONE new function to apply models to all files

---

## Architecture Overview

```
Phase 0: Assign pair_ids globally (existing code - no change)
    ↓
Phase 1: Sample complete pairs → load into memory
    ↓
Phase 2: Call sort_of_percolator_in_memory!(sampled_psms) → returns models
    ↓
Phase 3: Apply models to ALL files (new function)
    - 3 passes through files (one per iteration)
    - MBR tracker built during iteration 2
    - MBR features applied during iteration 3
```

---

## Detailed Design

### Phase 0: Pair ID Assignment (No Change)

Existing `assign_pair_ids_oom!` function handles this:
- Pass 1: Collect unique precursors from all files
- Pass 2: Assign pair_ids using same algorithm as in-memory
- Pass 3: Write pair_ids to all files

### Phase 1: Sample PSMs for Training (Minor Change)

Existing `sample_complete_pairs_for_training` function:
- Collects pair counts from all files
- Selects complete pairs up to memory limit
- Loads selected PSMs into memory

**No changes needed** - returns a DataFrame ready for in-memory processing.

### Phase 2: Train on Sample (New - Reuse In-Memory)

```julia
# Simply call the in-memory function on sampled data!
models = sort_of_percolator_in_memory!(
    sampled_psms,
    features,
    match_between_runs;
    # ... all the same parameters
)
```

The in-memory function:
- Assigns pair_ids (already done, but harmless to redo on sample)
- Trains 6 models (2 folds × 3 iterations)
- Returns `Dict{UInt8, LightGBMModelVector}` where each fold has 3 models

**Output structure**:
```
models[fold_1] = [model_iter1, model_iter2, model_iter3]
models[fold_2] = [model_iter1, model_iter2, model_iter3]
```

### Phase 3: Apply Models to All Files (New Function)

This is the only new code needed. Here's the detailed algorithm:

```julia
function apply_models_to_files!(
    models::Dict{UInt8, LightGBMModelVector},
    file_paths::Vector{String},
    features::Vector{Symbol},
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    iter_scheme::Vector{Int} = [100, 200, 200]
)
    n_iterations = length(iter_scheme)
    mbr_start_iter = n_iterations  # MBR features used in last iteration
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    # MBR tracker - accumulates best PSMs across all files
    mbr_tracker = Dict{Tuple{UInt32, Tuple{Int8,Int8}}, MBRTrackerValue}()

    # Process each iteration
    for iter in 1:n_iterations
        current_features = iter < mbr_start_iter ? non_mbr_features : features
        is_last_iter = (iter == n_iterations)
        is_mbr_accumulation_iter = (iter == mbr_start_iter - 1)  # Usually iter 2

        # Pass through all files for this iteration
        for file_path in file_paths
            # Load file (copycols=true for writability)
            df = DataFrame(Arrow.Table(file_path); copycols=true)

            # Apply quantile binning if needed
            apply_quantile_bins!(df, bin_edges)

            # Predict using fold-specific models
            trace_probs = predict_by_cv_fold(models, df, current_features, iter)

            if is_mbr_accumulation_iter && match_between_runs
                # Accumulate MBR statistics for iteration 2
                update_mbr_tracker!(mbr_tracker, df, trace_probs, max_q_value_lightgbm_rescore)
            end

            if is_last_iter && match_between_runs
                # Apply MBR features from tracker before final prediction
                apply_mbr_features!(df, mbr_tracker)
                # Re-predict with MBR features
                trace_probs = predict_by_cv_fold(models, df, features, iter)
            end

            # Write predictions to file
            df.trace_prob = trace_probs

            if is_last_iter
                # Final iteration: compute final outputs
                write_final_predictions(file_path, df, trace_probs, ...)
            else
                # Intermediate: just save trace_prob for MBR computation
                Arrow.write(file_path, df)
            end
        end

        @info "Iteration $iter complete"
    end
end

function predict_by_cv_fold(models, df, features, iter)
    trace_probs = zeros(Float32, nrow(df))
    for (fold_idx, fold_models) in pairs(models)
        fold_rows = findall(==(fold_idx), df.cv_fold)
        if !isempty(fold_rows)
            model = fold_models[iter]  # Get model for this iteration
            trace_probs[fold_rows] = predict(model, df[fold_rows, :])
        end
    end
    return trace_probs
end
```

---

## MBR Feature Handling

The tricky part is MBR features, which require knowing the "best PSM in other runs" for each pair.

### Current Approach (Complex)
- Interleaved with training
- Separate tracker logic in OOM
- Different from in-memory

### New Approach (Simple)
- **Iteration 1**: Predict all files, no MBR
- **Iteration 2**: Predict all files, accumulate MBR tracker
- **Iteration 3**: Apply MBR features from tracker, then predict

The MBR tracker from iteration 2 captures:
- Best trace_prob per (pair_id, isotopes) across all files
- Best PSM's features (weights, irts, etc.) for comparison
- Which runs have passing q-values

This is then used to compute MBR features for iteration 3:
- `MBR_num_runs`: Count of other passing runs
- `MBR_max_pair_prob`: Best prob from other run
- `MBR_rv_coefficient`: Correlation with best other run
- etc.

---

## Code Changes Summary

### Files to Modify

| File | Change |
|------|--------|
| `score_psms.jl` | Replace OOM logic with: sample → in-memory train → apply to files |
| `percolatorSortOf.jl` | Add `apply_models_to_files!` function (~150 lines) |
| `percolatorSortOf.jl` | DELETE `sort_of_percolator_out_of_memory!` (~400 lines) |

### Net Code Change
- **Add**: ~150 lines (apply_models_to_files!)
- **Remove**: ~400 lines (sort_of_percolator_out_of_memory!)
- **Net**: ~250 lines removed

---

## Questions to Consider

1. **Quantile bin edges**: Currently computed on the sample. Should we compute on ALL files first?
   - Current: Compute on sample, apply to all
   - Alternative: Quick pass through all files just for bin edge computation

2. **MBR feature iteration**: The in-memory version updates MBR features after EACH prediction within a fold. The simplified OOM would update after ALL files are predicted for an iteration.
   - This is a slight difference in timing
   - Should be equivalent because MBR features are only used in the final iteration

3. **Memory for MBR tracker**: The tracker stores best PSMs per (pair_id, isotopes). With ~100K pairs, this is manageable (~10-50MB).

4. **Number of file passes**:
   - Current OOM: 3 passes (one per iteration, but each pass does both read and write)
   - New OOM: 3 passes (same)
   - No change in I/O

---

## Implementation Steps

1. **Create `apply_models_to_files!` function** in percolatorSortOf.jl
   - Start simple: just apply models without MBR
   - Add MBR tracker accumulation for iteration 2
   - Add MBR feature application for iteration 3

2. **Modify `score_psms.jl`** OOM branch to:
   - Call `sort_of_percolator_in_memory!` on sample
   - Call `apply_models_to_files!` with returned models

3. **Test** with 275mb threshold
   - Compare precursor counts to in-memory
   - Should now match exactly (same training code!)

4. **Delete** `sort_of_percolator_out_of_memory!` once verified

---

## Expected Benefits

1. **Correctness**: By reusing in-memory code exactly, we eliminate all the subtle bugs
2. **Maintainability**: One training implementation to maintain, not two
3. **Simplicity**: ~250 fewer lines of code
4. **Consistency**: OOM and in-memory guaranteed to produce same training behavior
