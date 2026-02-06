# Out-of-Memory (OOM) LightGBM Scoring Pipeline Implementation Plan

## Git Commit Strategy

**IMPORTANT:** Make commits frequently throughout implementation:
- **Before starting**: Commit any uncommitted work on the branch
- **After each phase**: Commit completed infrastructure, streaming functions, etc.
- **After each major function**: Commit working incremental progress
- **After integration**: Commit the assembled OOM entry point
- **After testing**: Commit with test results noted

Commit messages should be descriptive (e.g., "Add streaming pair_id assignment for OOM scoring", "Implement global q-value calculation across files").

---

## Overview

Implement an out-of-memory scoring pipeline that allows training LightGBM models on datasets too large to fit in memory by:
1. Sampling PSMs by complete pairs (preserving pair_id integrity)
2. Training models on the sampled data
3. Applying models to ALL files via streaming
4. Computing GLOBAL q-values across all files (critical for correct iteration filtering)
5. Re-sampling eligible PSMs for each iteration based on global q-values
6. Computing MBR features across files without loading all data

## Critical Design Insight: Global Q-Values

**Problem**: The in-memory implementation calculates q-values on ALL data to determine:
- Which PSMs pass the q-value threshold for iterations 2 and 3
- Which weak targets get negative-mined (PEP threshold)

If we only calculate q-values on sampled data, we get biased estimates because:
- Q-value = (# decoys at score) / (# targets at score)
- Sampling changes this ratio, leading to incorrect filtering

**Solution**: After each iteration:
1. Apply model to ALL files (streaming)
2. Calculate GLOBAL q-values across ALL files
3. Use global q-values to determine training eligibility for next iteration
4. Re-sample fresh from eligible PSMs (different PSMs each iteration is OK)

## Configuration Changes

### New Parameter
Add to `optimization.machine_learning` in `defaultSearchParams.json`:

```json
"max_psm_memory_mb": 0  // 0 = use max_psms_in_memory; >0 = memory-based threshold in MB
```

### Files to Modify
- `assets/example_config/defaultSearchParams.json` - Add new parameter
- `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl` - Add default value
- `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl` - Add validation

---

## Implementation Steps

### Step 1: Memory Estimation Function
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Calculate exact memory based on Arrow schema (column types + row count):
```julia
function estimate_arrow_dataframe_memory_mb(file_paths::Vector{String})::Float64
    # Read schema from first file header (no data loading)
    schema = Arrow.Table(file_paths[1]).schema

    # Calculate bytes per PSM based on column types
    bytes_per_psm = sum(sizeof_column_type(col.type) for col in schema)

    # Count total rows across all files (lazy - just reads metadata)
    total_rows = sum(length(Arrow.Table(fp)) for fp in file_paths)

    return (bytes_per_psm * total_rows) / (1024 * 1024)
end

function calculate_max_sampleable_psms(file_paths::Vector{String}, max_mb::Int)::Int
    schema = Arrow.Table(file_paths[1]).schema
    bytes_per_psm = sum(sizeof_column_type(col.type) for col in schema)
    return floor(Int, (max_mb * 1024 * 1024) / bytes_per_psm)
end
```

This gives exact memory estimates without loading data, and tells us exactly how many PSMs we can sample within the memory budget.

### Step 2: OOM Decision Logic
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Modify `score_precursor_isotope_traces` (line 76-77) to use memory-based threshold:
```julia
use_oom = if max_psm_memory_mb > 0
    estimate_arrow_dataframe_memory_mb(file_paths) > max_psm_memory_mb
else
    psms_count >= max_psms_in_memory
end
```

### Step 2b: Global Quantile Feature Preprocessing (DISABLED FOR NOW)
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

**STATUS: Disabled for initial OOM implementation**

For now, the OOM version will NOT use quantile-binned features (`_qbin` columns). The feature set will exclude these features.

**TODO for future optimization:**
- Use a streaming statistics package (e.g., OnlineStats.jl) to compute quantile boundaries without loading all values into memory
- Then stream through files to add `_qbin` columns
- This avoids the O(N) memory requirement for collecting all feature values

```julia
# PLACEHOLDER - not implemented in initial version
function add_global_quantile_features_streaming!(file_paths, features_to_bin, n_bins)
    @warn "Quantile features disabled for OOM version - using non-binned features only"
    # Future: Use OnlineStats.jl for streaming quantile estimation
    # Future: Stream through files to apply boundaries
    return nothing
end
```

**For initial implementation:** Simply exclude `_qbin` features from the OOM feature set. The model will use raw feature values instead.

### Step 3: Streaming pair_id Assignment
**File:** `src/utils/ML/percolatorSortOf.jl`

**Arrow.jl Memory Behavior (Important):**
- `Arrow.Table(file)` uses memory mapping (mmap) - OS maps file bytes into virtual address space
- Columns are `ArrowVector` views into memory-mapped bytes, NOT copies
- `Tables.columntable(Arrow.Table(...))` MATERIALIZES data into Julia arrays (copies to RAM)
- Memory mapping ≠ lazy loading: files are mapped upfront, OS handles paging

**Conservative OOM approach:** Iterate through files ONE AT A TIME for guaranteed memory control:

```julia
function collect_unique_precursors_streaming(file_paths::Vector{String})
    precursor_info = Dict{UInt32, @NamedTuple{
        irt_pred::Float32,
        cv_fold::UInt8,
        isotopes_captured::Tuple{Int8,Int8}
    }}()

    # Process files one at a time for predictable memory usage
    for file_path in file_paths
        table = Arrow.Table(file_path)  # Memory-mapped, columns are views

        # Iterate through this file's rows
        for i in 1:length(table.precursor_idx)
            pid = table.precursor_idx[i]
            if !haskey(precursor_info, pid)
                precursor_info[pid] = (
                    irt_pred = table.irt_pred[i],
                    cv_fold = table.cv_fold[i],
                    isotopes_captured = table.isotopes_captured[i]
                )
            end
        end
        # table goes out of scope, OS can reclaim mapped memory
    end
    return precursor_info
end
```

**Note:** `Arrow.Table(file_paths::Vector)` also works and concatenates files virtually, but processing files individually gives us more control over memory and allows progress tracking per file.

1. `collect_unique_precursors_streaming(file_paths)` - Iterate files one-by-one, collect unique precursor info
2. `assign_pair_ids_streaming!(file_paths)` - Assign pair_ids:
   - Collect unique precursors via file-by-file streaming
   - Calculate iRT bins using existing `getIrtBins()`
   - Group by (irt_bin, cv_fold, isotopes_captured)
   - Shuffle and pair consecutive precursors (seed 1844)
   - Return `Dict{UInt32, UInt32}` mapping precursor_idx → pair_id
3. `write_pair_ids_to_files!(file_paths, precursor_to_pair_id)` - Add pair_id column to all files

### Step 4: Pair-Preserving Sampling
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Create `sample_psms_by_pairs(file_paths, precursor_to_pair_id, target_pair_count)`:
- Get unique pair_ids from dictionary
- Sample complete pairs (not individual PSMs)
- Load only PSMs belonging to sampled pairs
- Maintains pair integrity for MBR feature computation

### Step 5: OOM Training Entry Point (Corrected Iterative Workflow)
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Create `score_precursor_isotope_traces_oom!()` with proper global q-value calculation:

```julia
function score_precursor_isotope_traces_oom!(
    second_pass_folder, file_paths, precursors, model_config,
    match_between_runs, max_q_value_lightgbm_rescore,
    max_q_value_mbr_itr, min_PEP_neg_threshold_itr, n_quantile_bins;
    target_pair_count = 500_000,
    iter_scheme = [100, 200, 200]
)
    # Step 0: Quantile features DISABLED for OOM version
    # TODO: Future - implement streaming quantile computation with OnlineStats.jl
    # For now, OOM uses non-binned features only (exclude _qbin from feature set)

    # Step 1: Assign pair_ids to all files (streaming, one-time)
    @info "OOM Step 1: Assigning pair_ids..."
    precursor_to_pair_id = assign_pair_ids_streaming!(file_paths)
    write_pair_ids_to_files!(file_paths, precursor_to_pair_id)

    models = Dict{UInt8, Vector{LightGBM.Booster}}()  # Store models per CV fold

    for (itr, num_rounds) in enumerate(iter_scheme)
        @info "OOM Iteration $itr: Training with $num_rounds rounds"

        # Step 2: Sample PSMs for this iteration
        # Note: Quantile features already in files from Step 0
        if itr == 1
            # Iteration 1: Sample from ALL PSMs (no q-value filtering)
            sampled_psms = sample_psms_by_pairs(file_paths, precursor_to_pair_id,
                                                 target_pair_count)
        else
            # Iterations 2+: Sample from training-eligible PSMs only
            # (decoys + targets passing GLOBAL q-value threshold)
            sampled_psms = sample_psms_by_pairs_with_qvalue_filter(
                file_paths, precursor_to_pair_id, target_pair_count,
                max_q_value_lightgbm_rescore
            )

            # Iteration 3: Also apply negative mining using global PEP
            if itr == 3
                apply_negative_mining!(sampled_psms, min_PEP_neg_threshold_itr)
            end
        end

        # Quantile features already present from global preprocessing (Step 0)

        # Step 3: Train single iteration model on sampled data
        # Note: We train ONE iteration at a time, not all 3 together
        features = get_features_for_iteration(itr, match_between_runs)
        itr_models = train_single_iteration_cv!(sampled_psms, features, num_rounds)

        # Store models for this iteration
        for (fold, model) in itr_models
            if !haskey(models, fold)
                models[fold] = Vector{LightGBM.Booster}()
            end
            push!(models[fold], model)
        end

        # Step 4: Apply model to ALL files (streaming)
        apply_models_streaming!(file_paths, itr_models, features)

        # Step 5: Calculate GLOBAL q-values across ALL files
        # This is critical - q-values must be based on ALL data, not just sampled
        calculate_global_qvalues_streaming!(file_paths)

        # Step 6: Compute MBR features for next iteration (if needed)
        if match_between_runs && itr < length(iter_scheme)
            compute_mbr_features_streaming!(file_paths, precursor_to_pair_id,
                                            max_q_value_lightgbm_rescore)
        end

        sampled_psms = nothing
        GC.gc()
    end

    # Final MBR computation
    if match_between_runs
        compute_mbr_features_streaming!(file_paths, precursor_to_pair_id,
                                        max_q_value_lightgbm_rescore)
    end

    return models
end
```

**Key Differences from In-Memory Version:**
1. Each iteration trains a SINGLE model (not all 3 together)
2. Models are applied to ALL files after each iteration
3. Global q-values are calculated from ALL files, not just sampled data
4. Fresh sampling each iteration based on global q-values
5. Different PSMs may be sampled in each iteration (this is correct behavior)

### Step 6: Streaming Model Application
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Create `apply_models_streaming!(file_paths, models, features)`:
- Load each file one at a time
- Apply appropriate CV fold model
- Write back `trace_prob` column
- Use existing `writeArrow()` for atomic file updates

**File editing approach (initial implementation):**
```julia
function apply_models_streaming!(file_paths, models, features)
    for file_path in file_paths
        # SIMPLE APPROACH: Load full DataFrame, edit, write back
        # TODO: Future optimization - edit single column without loading all columns
        # Could potentially use paired files (editable columns vs static columns)
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Apply model based on CV fold
        for (fold, model) in models
            mask = df.cv_fold .== fold
            if any(mask)
                df.trace_prob[mask] = predict(model, df[mask, features])
            end
        end

        writeArrow(file_path, df)
    end
end
```

**Future optimization ideas:**
- Edit single Arrow column without loading all columns (if Arrow.jl supports this)
- Use paired files: one with frequently-edited columns (trace_prob, q_value, MBR features), one with static features
- Memory-mapped editing if possible

### Step 6b: Global Q-Value Calculation (NEW - Critical)
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Create `calculate_global_qvalues_streaming!(file_paths)`:

**Pass 1:** Collect global score distribution
```julia
function collect_score_distribution_streaming(file_paths)
    # SIMPLE APPROACH: Collect all scores into memory
    # TODO: Future optimization - use streaming/histogram approach to avoid O(N) memory
    # Could use OnlineStats.jl or histogram-based approximation

    target_scores = Float32[]
    decoy_scores = Float32[]

    # Process files one at a time
    for file_path in file_paths
        table = Arrow.Table(file_path)  # Memory-mapped views

        for i in 1:length(table.trace_prob)
            if table.target[i]
                push!(target_scores, table.trace_prob[i])
            else
                push!(decoy_scores, table.trace_prob[i])
            end
        end
        # table goes out of scope after each file
    end

    sort!(target_scores, rev=true)
    sort!(decoy_scores, rev=true)
    return target_scores, decoy_scores
end
```

**Pass 2:** Write q-values to all files
```julia
function write_global_qvalues_streaming!(file_paths, target_scores, decoy_scores)
    for file_path in file_paths
        # Load full DataFrame for editing (see Step 6 note on future optimization)
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Calculate q-value for each PSM using global score distribution
        # MUST use same algorithm as in-memory version (get_qvalues!)
        df[!, :q_value] = calculate_qvalues_from_global_distribution(
            df.trace_prob, df.target, target_scores, decoy_scores
        )

        writeArrow(file_path, df)
    end
end
```

**CRITICAL:** The q-value calculation must use the exact same algorithm as the in-memory `get_qvalues!` function to ensure identical filtering behavior.

### Step 6c: Q-Value Filtered Sampling
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Create `sample_psms_by_pairs_with_qvalue_filter(...)`:

**CRITICAL:** Filtering rules MUST match in-memory `get_training_data_for_iteration!`:
- **Iteration 2+:** Include decoys + targets passing q-value threshold
- **Iteration 3:** Additionally apply negative mining (relabel weak targets as decoys based on PEP)

```julia
function sample_psms_by_pairs_with_qvalue_filter(
    file_paths, precursor_to_pair_id, target_pair_count, max_q_value
)
    # Pass 1: Find pair_ids that have training-eligible PSMs
    # Training-eligible = decoys OR targets passing q-value threshold
    eligible_pairs = Set{UInt32}()

    for file_path in file_paths
        table = Arrow.Table(file_path)  # Memory-mapped views

        for i in 1:length(table.precursor_idx)
            is_training_eligible = !table.target[i] ||  # All decoys
                                   table.q_value[i] <= max_q_value  # Passing targets
            if is_training_eligible
                push!(eligible_pairs, precursor_to_pair_id[table.precursor_idx[i]])
            end
        end
        # table goes out of scope after each file
    end

    # Sample from eligible pairs (sample pair_ids, not individual PSMs)
    sampled_pairs = sample_pairs(eligible_pairs, target_pair_count)

    # Pass 2: Load only training-eligible PSMs from sampled pairs
    return load_training_eligible_psms_for_pairs(
        file_paths, precursor_to_pair_id, sampled_pairs, max_q_value
    )
end

function load_training_eligible_psms_for_pairs(
    file_paths, precursor_to_pair_id, sampled_pairs, max_q_value
)
    result_dfs = DataFrame[]
    for file_path in file_paths
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Filter to: sampled pairs AND training-eligible
        mask = [
            precursor_to_pair_id[df.precursor_idx[i]] in sampled_pairs &&
            (!df.target[i] || df.q_value[i] <= max_q_value)
            for i in 1:nrow(df)
        ]

        if any(mask)
            push!(result_dfs, df[mask, :])
        end
    end
    return vcat(result_dfs...)
end
```

### Step 7: Streaming MBR Feature Computation
**File:** `src/utils/ML/percolatorSortOf.jl`

**Key insight:** Can't do groupby across files, so use dictionary keyed by (pair_id, isotopes_captured) to accumulate best-per-run data, then compute MBR features in second pass.

**Data structures:**
```julia
# Per-run best PSM data for MBR computation
struct RunBestPSM
    trace_prob::Float32
    weights::Vector{Float32}      # For rv_coefficient calculation
    irts::Vector{Float32}         # For rv_coefficient calculation
    weight::Float32               # Scalar for log2_weight_ratio
    log2_intensity_explained::Float32
    irt_residual::Float32         # irt_pred - irt_obs
    is_decoy::Bool
    passes_qvalue::Bool
end

# All MBR data for one (pair_id, isotopes) combination
struct PairMBRData
    # Sparse: only runs with PSMs for this pair
    run_best::Dict{UInt32, RunBestPSM}  # ms_file_idx -> best PSM data
end

# Main dictionary: O(unique_pairs × isotope_combinations) entries
const PairMBRDict = Dict{
    @NamedTuple{pair_id::UInt32, isotopes::Tuple{Int8,Int8}},
    PairMBRData
}
```

**Pass 1:** `collect_pair_mbr_data_streaming(file_paths, precursor_to_pair_id, max_q_value)`
```julia
function collect_pair_mbr_data_streaming(file_paths, precursor_to_pair_id, max_q_value)
    pair_mbr = PairMBRDict()

    # Process files one at a time for predictable memory usage
    for file_path in file_paths
        table = Arrow.Table(file_path)  # Memory-mapped views

        for i in 1:length(table.precursor_idx)
            pair_id = precursor_to_pair_id[table.precursor_idx[i]]
            key = (pair_id = pair_id, isotopes = table.isotopes_captured[i])
            ms_file_idx = table.ms_file_idx[i]
            trace_prob = table.trace_prob[i]

            # Initialize pair entry if needed
            if !haskey(pair_mbr, key)
                pair_mbr[key] = PairMBRData(Dict{UInt32, RunBestPSM}())
            end

            run_data = pair_mbr[key].run_best

            # Update if this is the best PSM for this run
            if !haskey(run_data, ms_file_idx) || trace_prob > run_data[ms_file_idx].trace_prob
                run_data[ms_file_idx] = RunBestPSM(
                    trace_prob,
                    copy(table.weights[i]),      # Small vectors, OK to copy
                    copy(table.irts[i]),
                    table.weight[i],
                    table.log2_intensity_explained[i],
                    table.irt_pred[i] - table.irt_obs[i],  # irt_residual
                    table.decoy[i],
                    table.q_value[i] <= max_q_value
                )
            end
        end
        # table goes out of scope after each file
    end

    return pair_mbr
end
```

**Pass 2:** `apply_mbr_features_streaming!(file_paths, pair_mbr, precursor_to_pair_id)`
```julia
function apply_mbr_features_streaming!(file_paths, pair_mbr, precursor_to_pair_id)
    for file_path in file_paths
        # Must load DataFrame for editing (see Step 6 note on future optimization)
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Initialize MBR columns
        initialize_mbr_columns!(df)

        for i in 1:nrow(df)
            pair_id = precursor_to_pair_id[df.precursor_idx[i]]
            key = (pair_id = pair_id, isotopes = df.isotopes_captured[i])

            if !haskey(pair_mbr, key)
                mark_mbr_missing!(df, i)
                continue
            end

            run_data = pair_mbr[key].run_best
            current_run = df.ms_file_idx[i]

            # Find best PSM from a DIFFERENT run (same logic as in-memory version)
            best_other = find_best_other_run(run_data, current_run)

            if isnothing(best_other)
                mark_mbr_missing!(df, i)
                continue
            end

            # Compute MBR features (same formulas as summarize_precursors!)
            compute_mbr_features!(df, i, best_other, run_data)
        end

        writeArrow(file_path, df)
    end
end
```

**Memory:** Dictionary size = O(unique_pairs × isotope_combinations × avg_runs_per_pair)
- Doesn't grow with total PSM count
- Each RunBestPSM is ~100-200 bytes (including small vectors)
- For 1M pairs × 2 isotope combos × 5 runs = ~1-2 GB worst case, typically much less

**Reuse opportunity:** The MBR feature formulas in `compute_mbr_features!` should match `summarize_precursors!` exactly. Could potentially extract shared helper functions.

### Step 8: Single-Iteration Training Function
**File:** `src/utils/ML/percolatorSortOf.jl`

Create `train_single_iteration_cv!(psms, features, num_rounds)`:
```julia
function train_single_iteration_cv!(
    psms::DataFrame,
    features::Vector{Symbol},
    num_rounds::Int;
    feature_fraction = 0.5,
    learning_rate = 0.15,
    # ... other LightGBM params
)
    # Train one model per CV fold (similar to existing code but for single iteration)
    models = Dict{UInt8, LightGBM.Booster}()
    unique_folds = unique(psms.cv_fold)

    for test_fold in unique_folds
        train_mask = psms.cv_fold .!= test_fold
        train_data = psms[train_mask, :]

        bst = train_booster(train_data, features, num_rounds; ...)
        models[test_fold] = bst
    end

    return models
end
```

This extracts the per-iteration training logic from `sort_of_percolator_in_memory!` so it can be called separately for OOM processing.

---

## Code Reuse Strategy

### Functions to Reuse Directly
- `getIrtBins()` - iRT bin calculation
- `assign_pair_ids()` - Core pairing logic within bins
- `update_pair_statistics()` - Running statistics for MBR
- `writeArrow()` - Atomic file writing
- `train_booster()` - LightGBM training wrapper
- Model configurations from `model_config.jl`

### Functions to Extract/Refactor
- **Extract from `sort_of_percolator_in_memory!`:**
  - Single-iteration CV training → `train_single_iteration_cv!()`
  - Feature selection per iteration → `get_features_for_iteration()`
  - Negative mining logic → `apply_negative_mining!()`

### New Streaming Functions
| Function | Purpose |
|----------|---------|
| `assign_pair_ids_streaming!()` | Assign pair_ids without loading all data |
| `apply_models_streaming!()` | Apply CV models to all files |
| `calculate_global_qvalues_streaming!()` | Compute q-values from all files |
| `sample_psms_by_pairs_with_qvalue_filter()` | Sample eligible PSMs by pair |
| `collect_pair_statistics_streaming()` | Gather MBR stats across files |
| `apply_mbr_features_streaming!()` | Write MBR features to files |

### Key Shared Data Structures
- `PairStatistics` named tuple for MBR tracking (already defined)
- LightGBM model vectors and training infrastructure
- Feature sets from `model_config.jl`

### Relationship to In-Memory Version
The in-memory `sort_of_percolator_in_memory!` trains all 3 iterations in one function because it can calculate q-values instantly on the full dataset. The OOM version must:
1. Train one iteration at a time
2. Apply to ALL files after each iteration
3. Calculate GLOBAL q-values before next iteration
4. Re-sample based on global q-values

Both versions share the core training logic (`train_booster`, feature sets, hyperparameters) but differ in the iteration orchestration.

---

## Critical Files to Modify

| File | Changes |
|------|---------|
| `score_psms.jl` | Add OOM entry point, memory estimation, sampling by pairs, streaming model application |
| `percolatorSortOf.jl` | Add streaming pair_id assignment, streaming MBR computation |
| `defaultSearchParams.json` | Add `max_psm_memory_mb` parameter |
| `paramDefaults.jl` | Add default value for new parameter |
| `paramsChecks.jl` | Add validation for new parameter |

---

## Verification Plan

### Unit Tests
1. Test pair_id assignment consistency between streaming and in-memory versions
2. Test that sampled pairs are never split
3. Test MBR feature computation matches in-memory on small datasets

### Integration Tests
1. Run on small dataset with artificially low `max_psm_memory_mb` to force OOM path
2. Compare results (FDR, protein groups) with in-memory processing
3. Verify memory usage stays within threshold

### Test Commands
```julia
# Force OOM processing on small test dataset
params = load_params("test_params.json")
params.optimization.machine_learning.max_psm_memory_mb = 10  # Force OOM
SearchDIA(params)

# Compare results
in_memory_results = load_results("results_in_memory/")
oom_results = load_results("results_oom/")
compare_fdr_curves(in_memory_results, oom_results)
```

---

## Memory Usage Analysis

### In-Memory (Current)
- All PSMs loaded: O(N) where N = total PSMs
- Can be 10+ GB for large experiments

### OOM (Proposed)
- Sampled PSMs: O(S) where S = target_pair_count × avg_psms_per_pair
- Pair statistics dictionary: O(P) where P = unique pairs × ~100 bytes
- **Global score arrays**: O(N × 4 bytes) for collecting all trace_prob values
- File processing: O(F) where F = largest single file
- Total: O(S + P + N×4 + F)

### Memory Tradeoff Note
The global q-value calculation requires storing all `trace_prob` values (4 bytes each) to compute the score distribution. For 10M PSMs, this is ~40 MB - much smaller than loading all PSM features (~10 GB). This is an acceptable tradeoff for correct q-value calculation.

If even this is too much, we could use:
- Histogram-based approximation (fixed memory)
- Streaming quantile estimation (fixed memory, approximate)
- File-based sorting (disk I/O tradeoff)

---

## Edge Cases to Handle

1. **Single-file experiments**: Skip cross-file MBR, process normally
2. **Very few PSMs**: Fall back to in-memory processing
3. **Missing pair members**: Handle orphan PSMs in MBR computation
4. **Interrupted processing**: Ensure atomic file writes prevent corruption
5. **Empty pair statistics**: Handle pairs with no passing PSMs gracefully

---

## Implementation Order

### Phase 1: Infrastructure
1. Configuration changes - Add `max_psm_memory_mb` parameter
2. Memory estimation function (schema-based, exact calculation)
3. OOM decision logic in `score_precursor_isotope_traces`
4. OOM feature set definition (exclude `_qbin` features for now)

### Phase 2: Core Streaming Operations
5. Streaming pair_id assignment (`assign_pair_ids_streaming!`) - uses lazy `Arrow.Table(file_paths)`
6. Basic pair-preserving sampling (`sample_psms_by_pairs`)
7. Streaming model application (`apply_models_streaming!`) - simple load/edit/write approach

### Phase 3: Global Q-Value Pipeline (Critical)
8. Score distribution collection (`collect_score_distribution_streaming`) - simple O(N) memory approach
9. Global q-value writing (`write_global_qvalues_streaming!`) - must match in-memory algorithm
10. Q-value filtered sampling (`sample_psms_by_pairs_with_qvalue_filter`) - must match in-memory filtering rules

### Phase 4: Single-Iteration Training
11. Extract `train_single_iteration_cv!` from existing code
12. Extract `get_features_for_iteration` logic (OOM version excludes _qbin features)
13. Extract `apply_negative_mining!` logic - must match in-memory rules exactly

### Phase 5: MBR Streaming
14. Pair statistics collection (`collect_pair_statistics_streaming`)
15. MBR feature application (`apply_mbr_features_streaming!`)

### Phase 6: Integration
16. Assemble `score_precursor_isotope_traces_oom!` entry point
17. Integration testing against in-memory version
18. Memory profiling and optimization

### Future Optimizations (TODO)
- Quantile features via streaming stats (OnlineStats.jl)
- Single-column Arrow editing without loading full DataFrame
- Paired files (editable vs static columns)
- Histogram-based q-value approximation for O(1) memory

---

## Current Codebase Reference

### Existing OOM Infrastructure (Currently Disabled)
- `score_psms.jl` line 76-77: OOM path disabled with `if false`
- `percolatorSortOf.jl` lines 437-805: Commented out OOM version with patterns for streaming

### Key Existing Functions to Reuse
- `getIrtBins()` (percolatorSortOf.jl:94-107): iRT binning - reuse directly
- `assign_pair_ids()` (percolatorSortOf.jl:140-168): Core pairing logic - reuse directly
- `update_pair_statistics()` (percolatorSortOf.jl:57-92): Running MBR statistics - reuse directly
- `train_booster()`: LightGBM training wrapper - reuse directly
- `writeArrow()`: Atomic file writing - reuse directly

### Functions to Refactor
- `sort_of_percolator_in_memory!()` (percolatorSortOf.jl:242-432):
  - **Cannot reuse directly** for OOM because it calculates q-values on sampled data only
  - Need to extract: single-iteration training, feature selection, negative mining
  - The iteration loop must be controlled by OOM entry point to interleave global q-value calculation

### Functions to Replace
- `sample_psms_for_lightgbm()` (score_psms.jl): Samples by PSM, not by pair
  - Replace with `sample_psms_by_pairs()` for pair integrity
- `load_psms_for_lightgbm()` (score_psms.jl): Loads all PSMs
  - Replace with streaming operations for OOM
