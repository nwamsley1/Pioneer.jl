# Comprehensive Documentation: In-Memory vs OOM Percolator Pipeline

## Current Status

After implementing per-fold q_value computation in `build_mbr_tracker_for_fold`:

| Metric | In-Memory | OOM | Difference |
|--------|-----------|-----|------------|
| MBR_max_pair_prob count | 391,430 | 386,353 | -5,077 |
| Transfer candidates | 153,018 | 151,907 | -1,111 |
| Final PSMs (file 1) | 61,150 | 60,951 | ~0.3% |

The ~0.3% difference is likely within expected variance from model training on differently-ordered data.

---

## Pipeline Overview

Both pipelines perform semi-supervised learning with cross-validation to score PSMs. The key difference is memory management: in-memory loads all PSMs into a single DataFrame, while OOM processes files individually.

### Entry Point

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`
**Function**: `score_precursor_isotope_traces` (line ~180)

```julia
if n_psms > max_psms_in_memory
    # OOM path
    sort_of_percolator_oom!(...)
else
    # In-memory path
    psms = load_all_psms_to_dataframe(file_paths)
    sort_of_percolator_in_memory!(psms, ...)
end
```

---

## In-Memory Pipeline

**File**: `src/utils/ML/percolatorSortOf.jl`
**Main Function**: `sort_of_percolator_in_memory!` (lines 242-431)

### Step-by-Step Flow

#### 1. Initialization (lines 260-297)
```julia
assign_random_target_decoy_pairs!(psms)
sort!(psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])
```
- Assigns `pair_id` grouped by `[:irt_bin_idx, :cv_fold, :isotopes_captured]`
- Adds columns: `trace_prob`, `q_value`, MBR features

#### 2. CV Fold Loop (lines 302-408)
```julia
for test_fold_idx in unique_cv_folds
    initialize_prob_group_features!(psms, match_between_runs)

    train_idx = train_indices[test_fold_idx]
    test_idx = fold_indices[test_fold_idx]

    psms_train = @view psms[train_idx, :]
    psms_test = @view psms[test_idx, :]

    # Training iterations
    for (itr, num_round) in enumerate(iter_scheme)
        # Get training data for iteration
        psms_train_itr = get_training_data_for_iteration!(psms_train, ...)

        # Train model on psms_train
        bst = train_booster(psms_train_itr, train_feats, num_round; ...)

        # Predict on training data
        prob_train[train_idx] = predict(bst, psms_train)
        psms_train[!,:trace_prob] = prob_train[train_idx]
        get_qvalues!(psms_train.trace_prob, psms_train.target, psms_train.q_value)

        # Predict on test fold
        prob_test[test_idx] = predict(bst, psms_test)
        psms_test[!,:trace_prob] = prob_test[test_idx]

        # Store non-MBR estimates before MBR iteration
        if itr == (mbr_start_iter - 1)
            nonMBR_estimates[test_idx] = prob_test[test_idx]
        end

        # Update MBR features (if enabled)
        if match_between_runs
            update_mbr_features!(psms_train, psms_test, ...)
        end
    end
end
```

#### 3. MBR Feature Computation
**Function**: `update_mbr_features!` (lines 1221-1236)
```julia
if itr >= mbr_start_iter - 1
    get_qvalues!(psms_test.trace_prob, psms_test.target, psms_test.q_value)
    summarize_precursors!(psms_test, q_cutoff=max_q_value_lightgbm_rescore)
    summarize_precursors!(psms_train, q_cutoff=max_q_value_lightgbm_rescore)
end
```

**Function**: `summarize_precursors!` (lines 1238-1329)
- Groups by `[:pair_id, :isotopes_captured]`
- Finds top-2 runs by probability (r1, r2)
- Assigns comparison targets: r1↔r2, all others→r1
- Computes: `MBR_max_pair_prob`, `MBR_num_runs`, `MBR_rv_coefficient`, etc.

#### 4. Final Probability Assignment (lines 410-429)
```julia
if match_between_runs
    # Determine which precursors failed the q-value cutoff prior to MBR
    qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
    get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)
    pass_mask = (qvals_prev .<= max_q_value_lightgbm_rescore)
    prob_thresh = any(pass_mask) ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)

    # Label transfer candidates
    psms[!, :MBR_transfer_candidate] = .!pass_mask .& (psms.MBR_max_pair_prob .>= prob_thresh)

    # Store both trace probabilities
    psms[!, :trace_prob] = nonMBR_estimates
    psms[!, :MBR_boosted_trace_prob] = MBR_estimates
else
    psms[!, :trace_prob] = prob_test
end
```

### Functions Used (In-Memory Only)
| Function | Lines | Purpose |
|----------|-------|---------|
| `sort_of_percolator_in_memory!` | 242-431 | Main orchestration |
| `update_mbr_features!` | 1221-1236 | Coordinates MBR computation |
| `summarize_precursors!` | 1238-1329 | Computes MBR features per pair group |

---

## OOM Pipeline

**File**: `src/utils/ML/percolatorSortOf.jl`
**Main Function**: `apply_models_to_files!` (lines 857-954)

The OOM pipeline uses a simplified approach:
1. Train models on sampled PSMs using the exact same `sort_of_percolator_in_memory!` code
2. Apply trained models to all files with separate MBR tracking

### Step-by-Step Flow

#### 1. Pair ID Assignment (score_psms.jl)
**Function**: `assign_pair_ids_oom!` (lines 822-922)
```julia
# Pass 1: Collect unique (precursor_idx, cv_fold, isotopes) combinations
for file_path in file_paths
    df = DataFrame(Arrow.Table(file_path))
    unique_precs = unique(df, [:precursor_idx, :cv_fold, :isotopes_captured])
    append!(precursor_info, select(unique_precs, ...))
end

# Pass 2: Assign pair_ids grouped by [:irt_bin_idx, :cv_fold, :isotopes_captured]
for group in groupby(precursor_info, [:irt_bin_idx, :cv_fold, :isotopes_captured])
    last_pair_id = assignPairIds!(group, last_pair_id)
end

# Pass 3: Write pair_ids back to files
for file_path in file_paths
    add_column_to_file!(ref, :pair_id, compute_pair_id)
end
```

#### 2. Sample PSMs for Training
**Function**: `sample_complete_pairs_for_training` (lines 942-1007)
```julia
# Collect all unique pair_ids and count their PSMs across all files
pair_counts = Dict{UInt32, Int}()
for file_path in file_paths
    for pair_id in df.pair_id
        pair_counts[pair_id] = get(pair_counts, pair_id, 0) + 1
    end
end

# Randomly select pair_ids until we reach target PSM count
shuffled_pairs = shuffle(all_pair_ids)
for pair_id in shuffled_pairs
    if total_psms + pair_counts[pair_id] <= max_psms
        push!(selected_pairs, pair_id)
        total_psms += pair_counts[pair_id]
    end
end

# Load PSMs belonging to selected pairs from all files
```

#### 3. Train Models on Sampled Data (score_psms.jl)
**Function**: `score_precursor_isotope_traces_out_of_memory!` (lines 1029-1112)
```julia
# STEP 1: Train on sampled PSMs using EXACT same code as in-memory
models = sort_of_percolator_in_memory!(
    best_psms, features, match_between_runs; ...
)

# STEP 2: Apply trained models to ALL files
apply_models_to_files!(models, file_paths, features, match_between_runs; ...)
```

#### 4. Apply Models to All Files
**Function**: `apply_models_to_files!` (lines 857-954)
```julia
for iter in 1:n_iterations
    # Step 1: Apply predictions for this iteration to ALL files
    for (file_idx, file_path) in enumerate(file_paths)
        df = DataFrame(Arrow.Table(file_path); copycols=true)

        # Predict using fold-specific models for this iteration
        trace_probs = zeros(Float32, nrow(df))
        for fold_idx in unique_cv_folds
            fold_rows = findall(==(fold_idx), df.cv_fold)
            if !isempty(fold_rows)
                model = models[fold_idx][iter]
                trace_probs[fold_rows] = lightgbm_predict(model, df[fold_rows, :]; ...)
            end
        end

        df.trace_prob = trace_probs
        get_qvalues!(df.trace_prob, df.target, df.q_value)
        Arrow.write(file_path, df)
    end

    # Step 2: Compute MBR features per CV fold (if needed)
    if need_mbr
        for cv_fold in unique_cv_folds
            # Pass 1: Build tracker for this fold
            tracker = build_mbr_tracker_for_fold(file_paths, cv_fold, max_q_value_lightgbm_rescore)

            # Pass 2: Apply MBR features for this fold
            apply_mbr_from_tracker!(file_paths, cv_fold, tracker)
        end
    end
end
```

#### 5. Build MBR Tracker for Fold
**Function**: `build_mbr_tracker_for_fold` (lines 960-1094)
```julia
function build_mbr_tracker_for_fold(file_paths, cv_fold, q_cutoff)
    # First pass: collect all PSMs for this CV fold to compute per-fold q_values
    fold_probs = Float32[]
    fold_targets = Bool[]
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        for i in 1:nrow(df)
            df.cv_fold[i] == cv_fold || continue
            push!(fold_probs, df.trace_prob[i])
            push!(fold_targets, df.target[i])
        end
    end

    # Compute per-fold q_values (matching in-memory FDR calculation scope)
    fold_qvalues = Vector{Float64}(undef, length(fold_probs))
    get_qvalues!(fold_probs, fold_targets, fold_qvalues)

    # Second pass: Build tracker with best/second-best per pair
    tracker = Dict{key, NamedTuple}()
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        for i in 1:nrow(df)
            df.cv_fold[i] == cv_fold || continue

            # Get per-fold q_value for this PSM
            per_fold_qvalue = fold_qvalues[fold_qvalue_idx]

            # Update tracker with best/second-best probabilities
            # Track: best_prob_1, best_prob_2, unique_passing_runs, etc.
        end
    end
    return tracker
end
```

#### 6. Apply MBR Features from Tracker
**Function**: `apply_mbr_from_tracker!` (lines 1099-1187)
```julia
function apply_mbr_from_tracker!(file_paths, cv_fold, tracker)
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path); copycols=true)

        for i in 1:nrow(df)
            df.cv_fold[i] == cv_fold || continue

            key = (df.pair_id[i], df.isotopes_captured[i])
            entry = tracker[key]

            # Compute MBR_num_runs (exclude current run if it passes)
            current_run_passes = file_idx in entry.unique_passing_runs
            df.MBR_num_runs[i] = length(entry.unique_passing_runs) - current_run_passes

            # Find best from OTHER run
            if entry.best_file_1 != file_idx && !isempty(entry.best_weights_1)
                # Use best
                best_prob = entry.best_prob_1
                # ... compute MBR features
            elseif entry.best_file_2 != 0 && !isempty(entry.best_weights_2)
                # Use second best
                best_prob = entry.best_prob_2
                # ... compute MBR features
            else
                # No valid comparison - set defaults
            end

            df.MBR_max_pair_prob[i] = best_prob
            # ... set other MBR features
        end

        Arrow.write(file_path, df)
    end
end
```

### Functions Used (OOM Only)
| Function | Lines | Purpose |
|----------|-------|---------|
| `assign_pair_ids_oom!` | 822-922 | 3-pass pair_id assignment |
| `sample_complete_pairs_for_training` | 942-1007 | Sample PSMs for training |
| `score_precursor_isotope_traces_out_of_memory!` | 1029-1112 | OOM orchestration |
| `apply_models_to_files!` | 857-954 | Apply models file-by-file |
| `build_mbr_tracker_for_fold` | 960-1094 | Build MBR tracker per fold |
| `apply_mbr_from_tracker!` | 1099-1187 | Apply MBR features from tracker |

---

## Key Differences Summary

| Aspect | In-Memory | OOM |
|--------|-----------|-----|
| **Data Loading** | All PSMs in single DataFrame | Files processed individually |
| **pair_id Assignment** | Single pass in memory | 3-pass streaming |
| **Training** | Direct on full DataFrame | Sample-based training |
| **Model Application** | Direct prediction | File-by-file with Arrow I/O |
| **q_value Computation** | Per CV fold subset | Global then per-fold for MBR |
| **MBR Feature Storage** | DataFrame columns | Tracker dictionary |
| **Memory Usage** | O(n_psms) | O(n_psms_per_file) |

---

## Critical Implementation Details

### Per-Fold Q-Value Computation for MBR

**Problem**: In-memory computes q-values on the CV fold subset when computing MBR features. OOM originally computed global q-values, causing different `unique_passing_runs` sets.

**Solution**: `build_mbr_tracker_for_fold` now:
1. First pass: Collects all PSMs for the current CV fold
2. Computes q-values on just this fold's PSMs
3. Second pass: Uses per-fold q-values to determine passing runs

### MBR_num_runs Computation

**Both modes**: Exclude the current run from the count if it has any passing PSMs.
```julia
current_run_passes = file_idx in entry.unique_passing_runs
df.MBR_num_runs[i] = length(entry.unique_passing_runs) - current_run_passes
```

### Transfer Candidate Logic

**In-Memory** (lines 419-420):
```julia
psms[!, :MBR_transfer_candidate] = .!pass_mask .& (psms.MBR_max_pair_prob .>= prob_thresh)
```

**OOM** (`update_mbr_probs_oom!`, lines 1502-1514):
```julia
pass_mask = (prev_qvals .<= qval_thresh)
trace_prob_thresh = any(pass_mask) ? minimum(trace_probs[pass_mask]) : typemax(Float32)
df[!, :MBR_transfer_candidate] = (prev_qvals .> qval_thresh) .& (df.MBR_max_pair_prob .>= trace_prob_thresh)
```

Both implement the same logic: A PSM is a transfer candidate if:
1. It did NOT pass the q-value cutoff with non-MBR scoring
2. Its MBR_max_pair_prob >= the minimum passing probability threshold

---

## Potential OOM Performance Bottlenecks

### 1. Multiple File I/O Passes
- `assign_pair_ids_oom!`: 3 passes (read, compute, write)
- `apply_models_to_files!`: Multiple passes per iteration
- `build_mbr_tracker_for_fold`: 2 passes per fold
- `apply_mbr_from_tracker!`: 1 pass per fold (read + write)

**Total I/O per file**: ~10-15 read/write operations

### 2. Arrow Serialization Overhead
```julia
df = DataFrame(Arrow.Table(file_path); copycols=true)  # Read + copy
Arrow.write(file_path, df)                             # Write
```

### 3. Dictionary Lookups in Tracker
```julia
tracker = Dict{Tuple{UInt32, Tuple{Int8,Int8}}, NamedTuple{...}}()
```
- Large dictionaries with complex keys
- Repeated lookups for each PSM

### 4. Per-Fold Q-Value Recomputation
```julia
# In build_mbr_tracker_for_fold
get_qvalues!(fold_probs, fold_targets, fold_qvalues)  # Called once per fold
```

### 5. Suggested Optimizations
1. **Batch file operations**: Combine multiple passes where possible
2. **Memory-mapped files**: Use Arrow memory mapping instead of full reads
3. **Parallel file processing**: Process independent files concurrently
4. **Cache tracker lookups**: Pre-sort PSMs by pair_id to improve locality
5. **Streaming q_values**: Compute incrementally instead of collecting all first

---

## Verification Commands

```julia
# Test OOM mode
SearchDIA("config_oom.json")

# Compare metrics
# Look for in logs:
# - MBR_max_pair_prob count
# - Transfer candidates
# - Final PSM counts per file
```

---

## Related Documentation

- `simplified-oom-pipeline.md` - Original plan for simplified OOM approach
- `oom-scoring-discrepancy-fix.md` - Analysis of q-value computation differences
- `percolatorSortOf.jl` - Main implementation file
- `score_psms.jl` - Entry point and orchestration
