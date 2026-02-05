# Out-of-Memory Percolator Scoring with Trait Abstraction

## Goal

Add out-of-memory (OOM) processing capability to percolator scoring using a `MemoryStrategy` trait. The algorithm remains identical - only data access changes. Keep implementations side-by-side for maintainability.

**IMPORTANT**: This commit serves as a checkpoint. If the implementation doesn't go well, revert to this commit and rework the plan.

---

## Design Overview

### New Trait: `MemoryStrategy`

```julia
abstract type MemoryStrategy end

struct InMemoryProcessing <: MemoryStrategy end

struct OutOfMemoryProcessing <: MemoryStrategy
    max_training_psms::Int
    fold_file_paths::Dict{UInt8, Vector{String}}  # CV fold -> file paths
end
```

### Key Design Decisions

1. **Sample by pair_id**: Sample pair_ids (not individual PSMs) to keep pairs together for MBR
2. **Sample once upfront**: Use same sampled pairs for all training iterations
3. **Q-values during training**: Compute on sampled training data only
4. **Q-values during prediction**: Two-pass approach - predict all files, then compute global q-values across all files in CV fold
5. **Implementations side-by-side**: Trait-dispatched helpers defined together for easy maintenance

---

## Workflow Comparison

### In-Memory (Current)

```
1. Load all PSMs into memory
2. For each CV fold:
   a. Training: Train on training split (view)
   b. Prediction: Predict on test split (view)
3. All q-values computed on in-memory views
```

### Out-of-Memory (New)

```
Phase 0: Setup
  1. Scan all files to collect pair_id counts per file
  2. Sample pair_ids proportionally (respecting max_training_psms budget)
  3. Load ONLY sampled PSMs into memory for training

Phase 1: Training (same algorithm, sampled data)
  For each CV fold:
    - Use sampled training data (already in memory)
    - Train models exactly as in-memory
    - Q-values computed on sampled data only

Phase 2: Prediction (file-by-file with two passes)
  For each CV fold:
    Pass 1 - Predict:
      For each file in fold:
        - Load file
        - Apply model for each iteration
        - Write predictions back to file

    Pass 2 - Global Q-values:
      - Read all scores from fold files
      - Compute global q-values
      - Write q-values back to files

    MBR Features (if enabled):
      - Collect pair statistics across files (streaming)
      - Apply to each file
```

---

## Implementation Plan

### Step 1: Add MemoryStrategy Trait (scoring_traits.jl)

```julia
#############################################################################
# 8. MemoryStrategy - In-Memory vs Out-of-Memory Processing
#############################################################################

abstract type MemoryStrategy end

"""
    InMemoryProcessing <: MemoryStrategy

All PSMs loaded into memory. Current default behavior.
"""
struct InMemoryProcessing <: MemoryStrategy end

"""
    OutOfMemoryProcessing <: MemoryStrategy

PSMs processed file-by-file. Only sampled training data held in memory.
"""
struct OutOfMemoryProcessing <: MemoryStrategy
    max_training_psms::Int
    fold_file_paths::Dict{UInt8, Vector{String}}
end
```

### Step 2: Add OOM Helper Functions (percolator_generic.jl)

Keep helpers together with in-memory versions for maintainability:

```julia
#############################################################################
# Memory Strategy Helpers (In-Memory and OOM side-by-side)
#############################################################################

#= Pair ID Sampling - OOM only =#

"""
Scan files to get pair_id -> (file_path, count) mapping.
"""
function scan_pair_ids(file_paths::Vector{String})
    pair_counts = Dict{UInt32, Int}()
    pair_to_files = Dict{UInt32, Vector{String}}()

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        pair_ids = tbl.pair_id
        for pid in pair_ids
            pair_counts[pid] = get(pair_counts, pid, 0) + 1
            files = get!(pair_to_files, pid, String[])
            if isempty(files) || files[end] != fpath
                push!(files, fpath)
            end
        end
    end

    return pair_counts, pair_to_files
end

"""
Sample pair_ids proportionally to achieve target PSM count.
Returns Set of sampled pair_ids.
"""
function sample_pair_ids(pair_counts::Dict{UInt32, Int}, max_psms::Int)
    total_psms = sum(values(pair_counts))
    if total_psms <= max_psms
        return Set(keys(pair_counts))  # Use all
    end

    # Sample proportionally
    sample_rate = max_psms / total_psms
    sampled = Set{UInt32}()
    rng = MersenneTwister(1776)

    for (pid, count) in pair_counts
        if rand(rng) < sample_rate
            push!(sampled, pid)
        end
    end

    return sampled
end

"""
Load only PSMs with sampled pair_ids from files.
"""
function load_sampled_psms(file_paths::Vector{String}, sampled_pairs::Set{UInt32})
    dfs = DataFrame[]
    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        df = DataFrame(tbl)
        mask = [pid in sampled_pairs for pid in df.pair_id]
        if any(mask)
            push!(dfs, df[mask, :])
        end
    end
    return vcat(dfs...)
end

#= Training Data Access - dispatched by MemoryStrategy =#

function get_training_psms(::InMemoryProcessing, psms, train_indices, fold)
    return get_view(psms, train_indices)
end

function get_training_psms(memory::OutOfMemoryProcessing, _, _, fold)
    # Get files for training (all folds except test fold)
    train_files = String[]
    for (f, paths) in memory.fold_file_paths
        if f != fold
            append!(train_files, paths)
        end
    end

    # Scan and sample pair_ids
    pair_counts, _ = scan_pair_ids(train_files)
    sampled_pairs = sample_pair_ids(pair_counts, memory.max_training_psms)

    # Load sampled PSMs
    df = load_sampled_psms(train_files, sampled_pairs)
    return DataFramePSMContainer(df, Val(:unsafe))
end

#= Prediction Processing - dispatched by MemoryStrategy =#

function apply_predictions!(::InMemoryProcessing, psms_test, test_indices,
                           fold_models, iteration_rounds, config, ...)
    # Current in-memory prediction logic
    process_fold_iterations!(PredictionPhase(), psms_test, ...)
end

function apply_predictions!(memory::OutOfMemoryProcessing, _, _, fold,
                           fold_models, iteration_rounds, config, ...)
    test_files = memory.fold_file_paths[fold]

    # Pass 1: Predict on each file
    for fpath in test_files
        df = DataFrame(Arrow.Table(fpath))
        psms = DataFramePSMContainer(df, Val(:unsafe))

        for (itr, _) in enumerate(iteration_rounds)
            model = fold_models[itr]
            probs = predict_scores(model, psms)
            set_column!(psms, :trace_prob, probs)
        end

        # Write predictions back
        Arrow.write(fpath, to_dataframe(psms))
    end

    # Pass 2: Compute global q-values across all test files
    compute_global_qvalues!(test_files)

    # MBR: Collect pair statistics, then apply
    if has_mbr_support(config.mbr_update)
        pair_stats = collect_pair_statistics(test_files)
        apply_mbr_features!(test_files, pair_stats)
    end
end

#= Global Q-value Computation (OOM) =#

function compute_global_qvalues!(file_paths::Vector{String})
    # Collect all scores and targets
    all_scores = Float32[]
    all_targets = Bool[]
    file_ranges = Tuple{Int,Int}[]  # (start, end) for each file

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        start_idx = length(all_scores) + 1
        append!(all_scores, tbl.trace_prob)
        append!(all_targets, tbl.target)
        push!(file_ranges, (start_idx, length(all_scores)))
    end

    # Compute global q-values
    q_values = zeros(Float64, length(all_scores))
    get_qvalues!(all_scores, all_targets, q_values)

    # Write back to each file
    for (i, fpath) in enumerate(file_paths)
        df = DataFrame(Arrow.Table(fpath))
        start_idx, end_idx = file_ranges[i]
        df.q_value = q_values[start_idx:end_idx]
        Arrow.write(fpath, df)
    end
end

#= MBR Statistics Collection (OOM) =#

function collect_pair_statistics(file_paths::Vector{String})
    # Collect running statistics per pair_id
    pair_stats = Dict{UInt32, PairStatistics}()

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        for i in 1:length(tbl.pair_id)
            pid = tbl.pair_id[i]
            stats = get!(pair_stats, pid, PairStatistics())
            update_stats!(stats, tbl, i)
        end
    end

    return pair_stats
end

function apply_mbr_features!(file_paths::Vector{String}, pair_stats)
    for fpath in file_paths
        df = DataFrame(Arrow.Table(fpath))
        for i in 1:nrow(df)
            pid = df.pair_id[i]
            stats = pair_stats[pid]
            apply_mbr_from_stats!(df, i, stats)
        end
        Arrow.write(fpath, df)
    end
end
```

### Step 3: Modify Main Entry Point (percolator_generic.jl)

```julia
"""
Main entry point - dispatches based on memory strategy.
"""
function percolator_scoring!(
    data,  # AbstractPSMContainer or file paths
    config::ScoringConfig;
    memory::MemoryStrategy = InMemoryProcessing(),
    show_progress::Bool = true,
    verbose::Bool = false
)
    return percolator_scoring_impl!(data, config, memory; show_progress, verbose)
end

# In-memory implementation (current code, minimal changes)
function percolator_scoring_impl!(
    psms::AbstractPSMContainer,
    config::ScoringConfig,
    ::InMemoryProcessing;
    show_progress::Bool = true,
    verbose::Bool = false
)
    # ... existing implementation ...
end

# OOM implementation
function percolator_scoring_impl!(
    _,  # psms not used - we load from files
    config::ScoringConfig,
    memory::OutOfMemoryProcessing;
    show_progress::Bool = true,
    verbose::Bool = false
)
    unique_cv_folds = collect(keys(memory.fold_file_paths))
    models = Dict{UInt8, Vector{Any}}()

    # Phase 1: Training (one fold at a time)
    for test_fold in unique_cv_folds
        # Load sampled training data
        training_psms = get_training_psms(memory, nothing, nothing, test_fold)

        # Apply pairing to sampled data
        assign_pairs!(training_psms, config.pairing)
        initialize_mbr_columns!(training_psms, config.mbr_update)

        # Train exactly as in-memory
        fold_models = Vector{Any}(undef, length(get_iteration_rounds(config.iteration_scheme)))
        # ... training loop using existing process_fold_iterations! ...

        models[test_fold] = fold_models
        training_psms = nothing  # Free memory
        GC.gc()
    end

    # Phase 2: Prediction (file-by-file)
    for test_fold in unique_cv_folds
        fold_models = models[test_fold]
        apply_predictions!(memory, nothing, nothing, test_fold, fold_models, ...)
    end

    return models
end
```

---

## Files to Modify

1. **scoring_traits.jl**
   - Add `MemoryStrategy`, `InMemoryProcessing`, `OutOfMemoryProcessing`

2. **percolator_generic.jl**
   - Add OOM helper functions (pair sampling, file iteration, global q-values)
   - Add `percolator_scoring_impl!` with memory dispatch
   - Keep in-memory and OOM implementations together

3. **mbr_update.jl** (optional)
   - Add streaming MBR statistics collection
   - Add `PairStatistics` struct for aggregation

---

## Key Invariants

1. **Same algorithm**: Training iterations, model selection, feature usage all identical
2. **Pair integrity**: All PSMs with same pair_id are either all sampled or all excluded
3. **Q-value correctness**:
   - Training: Local q-values on sampled data
   - Prediction: Global q-values across all files in CV fold
4. **MBR correctness**: Pair statistics aggregated across files before applying

---

## Verification

```julia
# Test with small dataset (should use in-memory)
SearchDIA("small_test_params.json")

# Test with large dataset (should use OOM)
SearchDIA("large_test_params.json")

# Compare results at same q-value thresholds
# - PSM counts should be similar (not identical due to sampling)
# - Protein group counts should be similar
```

---

## Performance Notes

- **Memory**: Only sampled training PSMs + one file at a time during prediction
- **I/O**: Two passes over test files (predict, then global q-values)
- **MBR**: Additional pass for pair statistics if MBR enabled
