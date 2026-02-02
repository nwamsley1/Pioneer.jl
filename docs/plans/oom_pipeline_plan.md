# Plan: Out-of-Memory (OOM) DIA Search Pipeline

## Overview
Restore and modernize OOM capabilities so Pioneer.jl can process datasets too large to fit in memory by training ML models on sampled PSMs and applying them file-by-file.

---

## File 1: ScoringSearchParameters (Add `max_psm_memory_mb`)

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

### Change 1.1: Add parameter to struct (line ~51)
```julia
struct ScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    # LightGBM parameters
    max_psms_in_memory::Int64
    max_psm_memory_mb::Int64  # NEW: User-specified memory limit for PSMs
    min_best_trace_prob::Float32
    # ... rest unchanged
```

### Change 1.2: Add to constructor (line ~93)
```julia
new{typeof(isotope_trace_type)}(
    Int64(ml_params.max_psms_in_memory),
    Int64(get(ml_params, :max_psm_memory_mb, 8000)),  # NEW: Default 8GB
    Float32(ml_params.min_trace_prob),
    # ... rest unchanged
```

---

## File 2: score_psms.jl (Enable OOM Path)

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

### Change 2.1: Add memory-based threshold calculation (after line 74)
```julia
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    max_psms_in_memory::Int64,
    max_psm_memory_mb::Int64,  # NEW parameter
    n_quantile_bins::Int64,
    q_value_threshold::Float32 = 0.01f0,
    ms1_scoring::Bool = true
)
    # Step 1: Count PSMs and determine processing approach
    psms_count = get_psms_count(file_paths)

    # NEW: Calculate memory-based threshold
    psm_size_bytes = 500  # Estimated bytes per PSM with features
    memory_based_threshold = floor(Int64, max_psm_memory_mb * 1e6 / psm_size_bytes)
    effective_threshold = min(max_psms_in_memory, memory_based_threshold)
```

### Change 2.2: Enable OOM path (replace line 77)
```julia
    # FROM:
    # if false  # psms_count >= max_psms_in_memory

    # TO:
    if psms_count >= effective_threshold
        # OOM MODE ACTIVATED
        @user_info "Using OUT-OF-MEMORY mode: $psms_count PSMs exceeds $effective_threshold threshold"
        @debug_l1 "\n[OOM] Out-of-memory percolator activated"
        @debug_l1 "\n[OOM]   PSM count: $psms_count"
        @debug_l1 "\n[OOM]   Memory threshold: $max_psm_memory_mb MB → $effective_threshold PSMs"
        @debug_l1 "\n[OOM]   Files: $(length(file_paths))"

        # Phase 0: Assign pair_ids across all files
        phase0_timing = @timed assign_pair_ids_oom!(file_paths)
        @debug_l1 "\n[OOM] Phase 0 (pair assignment): $(round(phase0_timing.time, digits=2))s"

        # Sample complete pair groups for training
        best_psms = sample_complete_pairs_for_training(
            file_paths,
            effective_threshold,
            precursors
        )

        # Add quantile-binned features
        features_to_bin = [:prec_mz, :irt_pred, :weight, :tic]
        add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)

        # Use AdvancedLightGBM config for OOM
        model_config = create_default_advanced_lightgbm_config(ms1_scoring)

        # Train and apply OOM percolator
        models = score_precursor_isotope_traces_out_of_memory!(
            best_psms,
            file_paths,
            precursors,
            model_config,
            match_between_runs,
            max_q_value_lightgbm_rescore,
            max_q_value_mbr_itr,
            min_PEP_neg_threshold_itr
        )
    else
        # In-memory processing (unchanged)
        # ...
```

### Change 2.3: Add `assign_pair_ids_oom!` function
```julia
"""
    assign_pair_ids_oom!(file_paths::Vector{String})

Assign pair_ids to ALL PSMs across ALL files before sampling.
This ensures pair groups are never split between training sample and OOM data.
"""
function assign_pair_ids_oom!(file_paths::Vector{String})
    @debug_l1 "\n[OOM] Collecting precursor info from $(length(file_paths)) files..."

    # Pass 1: Collect all unique (precursor_idx, irt_pred, cv_fold, isotopes_captured, target)
    precursor_info = DataFrame(
        precursor_idx = UInt32[],
        irt_pred = Float32[],
        cv_fold = UInt8[],
        isotopes_captured = Tuple{Int8,Int8}[],
        target = Bool[]
    )

    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        # Get unique precursors from this file
        unique_precs = unique(df, [:precursor_idx, :cv_fold, :isotopes_captured])
        append!(precursor_info, select(unique_precs,
            :precursor_idx, :irt_pred, :cv_fold, :isotopes_captured, :target))
    end

    # Deduplicate across files
    unique!(precursor_info, [:precursor_idx, :cv_fold, :isotopes_captured])

    @debug_l1 "\n[OOM] Found $(nrow(precursor_info)) unique precursor-isotope combinations"

    # Pass 2: Assign pair_ids using same algorithm as in-memory
    # Sort for deterministic ordering regardless of file processing order
    sort!(precursor_info, [:cv_fold, :isotopes_captured, :irt_pred, :precursor_idx])

    # Add irt_bin_idx
    precursor_info[!, :irt_bin_idx] = getIrtBins(precursor_info.irt_pred)

    # Assign pair_ids within groups
    last_pair_id = zero(UInt32)
    precursor_info[!, :pair_id] = zeros(UInt32, nrow(precursor_info))

    for group in groupby(precursor_info, [:irt_bin_idx, :cv_fold, :isotopes_captured])
        last_pair_id = assignPairIds!(group, last_pair_id)
    end

    # Create lookup: (precursor_idx, cv_fold, isotopes) -> pair_id
    pair_lookup = Dict{Tuple{UInt32, UInt8, Tuple{Int8,Int8}}, UInt32}()
    for row in eachrow(precursor_info)
        key = (row.precursor_idx, row.cv_fold, row.isotopes_captured)
        pair_lookup[key] = row.pair_id
    end

    @debug_l1 "\n[OOM] Assigned $(last_pair_id) pair_ids"

    # Pass 3: Write pair_ids back to all files
    @debug_l1 "\n[OOM] Writing pair_ids to files..."
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))

        # Add pair_id column
        df[!, :pair_id] = map(eachrow(df)) do row
            key = (row.precursor_idx, row.cv_fold, row.isotopes_captured)
            get(pair_lookup, key, zero(UInt32))
        end

        # Write back
        Arrow.write(file_path, df)
    end

    return pair_lookup
end
```

### Change 2.4: Add `sample_complete_pairs_for_training` function
```julia
"""
    sample_complete_pairs_for_training(file_paths, max_psms, precursors)

Sample COMPLETE pair groups for training. Never splits a pair between
training sample and OOM data.
"""
function sample_complete_pairs_for_training(
    file_paths::Vector{String},
    max_psms::Int64,
    precursors::LibraryPrecursors
)
    @debug_l1 "\n[OOM] Sampling complete pairs for training (target: $max_psms PSMs)..."

    # Collect all unique pair_ids and their PSM counts
    pair_counts = Dict{UInt32, Int}()
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        for pair_id in df.pair_id
            pair_counts[pair_id] = get(pair_counts, pair_id, 0) + 1
        end
    end

    all_pair_ids = collect(keys(pair_counts))
    @debug_l1 "\n[OOM] Found $(length(all_pair_ids)) unique pair_ids"

    # Shuffle and select pairs until we reach target PSM count
    Random.seed!(1776)
    shuffled_pairs = shuffle(all_pair_ids)

    selected_pairs = Set{UInt32}()
    total_psms = 0
    for pair_id in shuffled_pairs
        if total_psms + pair_counts[pair_id] <= max_psms
            push!(selected_pairs, pair_id)
            total_psms += pair_counts[pair_id]
        end
        if total_psms >= max_psms * 0.9  # Stop at 90% to avoid overshoot
            break
        end
    end

    @debug_l1 "\n[OOM] Selected $(length(selected_pairs)) pairs with $total_psms PSMs"

    # Load PSMs belonging to selected pairs
    sampled_dfs = DataFrame[]
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        mask = in.(df.pair_id, Ref(selected_pairs))
        if any(mask)
            push!(sampled_dfs, df[mask, :])
        end
    end

    sampled_psms = vcat(sampled_dfs...)
    @debug_l1 "\n[OOM] Loaded $(nrow(sampled_psms)) PSMs for training"

    return sampled_psms
end
```

---

## File 3: percolatorSortOf.jl (OOM Percolator)

**File**: `src/utils/ML/percolatorSortOf.jl`

### Change 3.1: Uncomment and update `sort_of_percolator_out_of_memory!`

The existing commented code needs these key updates:

```julia
"""
    sort_of_percolator_out_of_memory!(sampled_psms, file_paths, features, ...)

Out-of-memory percolator that:
1. Trains on sampled PSMs (complete pair groups)
2. Applies models file-by-file
3. Computes MBR features using tracker dictionary
"""
function sort_of_percolator_out_of_memory!(
    sampled_psms::DataFrame,
    file_paths::Vector{String},
    features::Vector{Symbol},
    match_between_runs::Bool = true;
    max_q_value_lightgbm_rescore::Float32 = 0.01f0,
    max_q_value_mbr_itr::Float32 = 0.20f0,
    min_PEP_neg_threshold_itr::Float32 = 0.90f0,
    iter_scheme::Vector{Int} = [100, 200, 200],
    # ... other params
)
    @debug_l1 "\n[OOM] Starting OOM percolator training on $(nrow(sampled_psms)) sampled PSMs"

    # Initialize MBR tracker dictionary
    mbr_tracker = Dictionary{
        @NamedTuple{pair_id::UInt32, isotopes::Tuple{Int8,Int8}},
        MBRTrackerEntry
    }()

    #=== PHASE A: Train on sampled PSMs ===#
    training_timing = @timed begin
        # This uses the existing in-memory logic but on the sampled subset
        models = train_cv_models_on_sample!(
            sampled_psms,
            features,
            match_between_runs,
            max_q_value_lightgbm_rescore,
            max_q_value_mbr_itr,
            min_PEP_neg_threshold_itr,
            iter_scheme
        )
    end
    @debug_l1 "\n[OOM] Phase A (training): $(round(training_timing.time, digits=2))s"

    #=== PHASE B: Apply to all files ===#
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]
    mbr_start_iter = length(iter_scheme)

    for (itr, _) in enumerate(iter_scheme)
        iter_timing = @timed begin
            train_feats = itr < mbr_start_iter ? non_mbr_features : features

            # Pass 1: Predict all files, update MBR tracker
            pass1_timing = @timed begin
                for file_path in file_paths
                    predict_and_update_tracker!(
                        file_path, models, train_feats, itr,
                        mbr_tracker, max_q_value_lightgbm_rescore,
                        match_between_runs && itr >= mbr_start_iter - 1
                    )
                end
            end
            @debug_l1 "\n[OOM]   Iteration $itr Pass 1 (predict): $(round(pass1_timing.time, digits=2))s"

            # Pass 2: Apply MBR features from tracker (only at mbr_start_iter - 1)
            if match_between_runs && itr == mbr_start_iter - 1
                pass2_timing = @timed begin
                    for file_path in file_paths
                        apply_mbr_features_from_tracker!(file_path, mbr_tracker)
                    end
                end
                @debug_l1 "\n[OOM]   Iteration $itr Pass 2 (MBR features): $(round(pass2_timing.time, digits=2))s"
            end
        end
        @debug_l1 "\n[OOM] Iteration $itr total: $(round(iter_timing.time, digits=2))s"
    end

    #=== PHASE C: Final scoring and transfer candidates ===#
    final_timing = @timed begin
        for file_path in file_paths
            finalize_oom_scores!(
                file_path, models, features, non_mbr_features,
                mbr_start_iter, match_between_runs,
                max_q_value_lightgbm_rescore
            )
        end
    end
    @debug_l1 "\n[OOM] Phase C (finalize): $(round(final_timing.time, digits=2))s"

    return models
end
```

### Change 3.2: Add MBR tracker struct
```julia
"""
MBR tracker entry for OOM processing.
Stores best/second-best scores and chromatogram data per (pair_id, isotopes) group.
"""
struct MBRTrackerEntry
    # Best scoring PSM
    best_prob_1::Float32
    best_ms_file_idx_1::UInt32
    best_log2_weights_1::Vector{Float32}
    best_irts_1::Vector{Float32}
    best_irt_residual_1::Float32
    best_weight_1::Float32
    best_log2_intensity_explained_1::Float32
    is_best_decoy_1::Bool

    # Second best (for when current run IS the best)
    best_prob_2::Float32
    best_ms_file_idx_2::UInt32
    best_log2_weights_2::Vector{Float32}
    best_irts_2::Vector{Float32}
    best_irt_residual_2::Float32
    best_weight_2::Float32
    best_log2_intensity_explained_2::Float32
    is_best_decoy_2::Bool

    # Runs passing q-value threshold
    unique_passing_runs::Set{UInt16}
end
```

### Change 3.3: Add `predict_and_update_tracker!` function
```julia
"""
Predict on a single file and update the MBR tracker.
"""
function predict_and_update_tracker!(
    file_path::String,
    models::Dict{UInt8, LightGBMModelVector},
    features::Vector{Symbol},
    iteration::Int,
    tracker::Dictionary,
    q_cutoff::Float32,
    update_mbr::Bool
)
    df = DataFrame(Arrow.Table(file_path))

    # Predict using CV-fold-appropriate models
    trace_probs = zeros(Float32, nrow(df))
    for (fold_idx, fold_models) in pairs(models)
        fold_mask = df.cv_fold .== fold_idx
        if any(fold_mask)
            model = fold_models[iteration]
            trace_probs[fold_mask] = predict(model, df[fold_mask, features])
        end
    end

    # Update trace_prob in file
    df[!, :trace_prob] = trace_probs

    # Update MBR tracker if needed
    if update_mbr
        # Compute q-values for this file
        qvals = zeros(Float32, nrow(df))
        get_qvalues!(trace_probs, df.target, qvals)

        for i in 1:nrow(df)
            key = (pair_id = df.pair_id[i], isotopes = df.isotopes_captured[i])
            prob = trace_probs[i]

            if haskey(tracker, key)
                entry = tracker[key]
                # Update tracker with new score
                tracker[key] = update_mbr_tracker_entry(
                    entry, prob, df, i, qvals[i] <= q_cutoff
                )
            else
                # Create new entry
                tracker[key] = create_mbr_tracker_entry(df, i, prob, qvals[i] <= q_cutoff)
            end
        end
    end

    # Write updated file
    Arrow.write(file_path, df)
end
```

### Change 3.4: Add `apply_mbr_features_from_tracker!` function
```julia
"""
Apply MBR features to a file using the tracker dictionary.
"""
function apply_mbr_features_from_tracker!(
    file_path::String,
    tracker::Dictionary
)
    df = DataFrame(Arrow.Table(file_path))

    # Initialize MBR columns
    n = nrow(df)
    df[!, :MBR_max_pair_prob] = zeros(Float32, n)
    df[!, :MBR_best_irt_diff] = zeros(Float32, n)
    df[!, :MBR_log2_weight_ratio] = zeros(Float32, n)
    df[!, :MBR_log2_explained_ratio] = zeros(Float32, n)
    df[!, :MBR_rv_coefficient] = zeros(Float32, n)
    df[!, :MBR_is_best_decoy] = trues(n)
    df[!, :MBR_num_runs] = zeros(Int32, n)
    df[!, :MBR_is_missing] = falses(n)

    for i in 1:nrow(df)
        key = (pair_id = df.pair_id[i], isotopes = df.isotopes_captured[i])

        if !haskey(tracker, key)
            df.MBR_is_missing[i] = true
            continue
        end

        entry = tracker[key]
        current_file = df.ms_file_idx[i]

        # Use second-best if current file is the best
        if entry.best_ms_file_idx_1 == current_file && entry.best_prob_2 > 0
            # Use second-best data
            best_log2_weights = entry.best_log2_weights_2
            best_irts = entry.best_irts_2
            best_irt_residual = entry.best_irt_residual_2
            best_weight = entry.best_weight_2
            best_log2_ie = entry.best_log2_intensity_explained_2
            is_best_decoy = entry.is_best_decoy_2
            max_pair_prob = entry.best_prob_2
        elseif entry.best_prob_1 > 0
            # Use best data
            best_log2_weights = entry.best_log2_weights_1
            best_irts = entry.best_irts_1
            best_irt_residual = entry.best_irt_residual_1
            best_weight = entry.best_weight_1
            best_log2_ie = entry.best_log2_intensity_explained_1
            is_best_decoy = entry.is_best_decoy_1
            max_pair_prob = entry.best_prob_1
        else
            df.MBR_is_missing[i] = true
            continue
        end

        # Compute MBR features (same logic as in-memory summarize_precursors!)
        current_residual = Float32(df.irt_pred[i] - df.irt_obs[i])

        df.MBR_max_pair_prob[i] = max_pair_prob
        df.MBR_best_irt_diff[i] = abs(best_irt_residual - current_residual)
        df.MBR_is_best_decoy[i] = is_best_decoy
        df.MBR_log2_weight_ratio[i] = log2(df.weight[i] / best_weight)
        df.MBR_log2_explained_ratio[i] = df.log2_intensity_explained[i] - best_log2_ie

        # Compute RV coefficient
        current_log2_weights = log2.(df.weights[i])
        current_irts = df.irts[i]
        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, current_log2_weights)
        best_irts_padded, irts_padded = pad_rt_equal_length(best_irts, current_irts)
        df.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_irts_padded, weights_padded, irts_padded)

        # MBR_num_runs: count of passing runs excluding current
        n_passing = length(entry.unique_passing_runs)
        current_in_passing = current_file in entry.unique_passing_runs
        df.MBR_num_runs[i] = n_passing - (current_in_passing ? 1 : 0)
    end

    Arrow.write(file_path, df)
end
```

---

## CRITICAL: Potential Divergence Points

### Acceptable Differences (Inherent to Sampling)
1. **Model weights** - Different training data → different models (OK)
2. **Q-value distribution** - Sample differs from full (OK)
3. **PEP-based relabeling** - Different score distribution (OK)

### Must Handle Correctly
4. **MBR Feature Timing** - Compute for ALL files at `mbr_start_iter - 1`
5. **Iteration State** - Update trace_prob in all files between iterations
6. **MBR Tracker** - Handle best/second-best when same file has both
7. **`MBR_num_runs`** - Requires seeing all files per pair
8. **CV Fold Handling** - Apply correct model per fold in each file

---

## Timing & Logging

### OOM Mode Indicator
```julia
@user_info "Using OUT-OF-MEMORY mode: $psms_count PSMs exceeds $threshold"
@debug_l1 "\n[OOM] Out-of-memory percolator activated"
@debug_l1 "\n[OOM]   PSM count: $psms_count"
@debug_l1 "\n[OOM]   Threshold: $threshold"
@debug_l1 "\n[OOM]   Files: $(length(file_paths))"
```

### Phase Timing Pattern
```julia
phase_timing = @timed begin
    # ... phase code ...
end
@debug_l1 "\n[OOM] Phase X: $(round(phase_timing.time, digits=2))s"
```

---

## Verification

1. **Unit tests**: Pair assignment, sampling, tracker updates
2. **Accuracy test**: OOM vs in-memory on same dataset (< 5% difference at 1% FDR)
3. **Memory test**: Peak memory stays below `max_psm_memory_mb`
4. **Integration test**: Downstream pipelines work unchanged

---

## Implementation Order

1. Add `max_psm_memory_mb` parameter to `ScoringSearchParameters`
2. Add `assign_pair_ids_oom!` function
3. Add `sample_complete_pairs_for_training` function
4. Add MBR tracker struct and helper functions
5. Implement `sort_of_percolator_out_of_memory!`
6. Enable OOM path in `score_precursor_isotope_traces`
7. Test and validate
