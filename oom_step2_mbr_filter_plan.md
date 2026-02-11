# OOM Step 2: Streaming MBR Filter + Per-File Precursor Aggregation

## Context

After `percolator_scoring!` (fully OOM-safe), **ScoringSearch.jl Step 2** (lines 461-505) loads ALL PSMs from every file into a single `merged_df` DataFrame for:
1. MBR filter training + application (`apply_mbr_filter!`)
2. Precursor probability aggregation (`transform!(groupby(...))`)

This is the **last remaining O(N) allocation** in the scoring pipeline. For a 10M PSM dataset with ~30 columns, `merged_df` costs ~1.2 GB plus copies during filtering.

**Three key insights make this fixable with O(1-file) memory:**

1. **Precursor aggregation groups by `(precursor_idx, ms_file_idx)`** — since `ms_file_idx` is a group key, each group lives within a single file. The merge was never needed.
2. **MBR filter models are low-capacity** (Probit = linear, LightGBM = depth-3 × 100 rounds). Training on a **sampled subset** gives nearly identical models. Same pattern as `_sample_by_pair_id` in scoring_workspace.jl:314.
3. **FTR threshold computation** can use the **streaming sort-merge** pattern from `_compute_prob_threshold_from_merged` (scoring_workspace.jl:705) — per-file sorted temp Arrow → merge → single-pass stream.

**Constraint**: At no point should more than one file's data be in memory simultaneously (beyond a bounded training sample).

## Files to Modify

| File | Changes |
|------|---------|
| `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl` | Replace Step 2 (lines 459-505) with single function call |
| `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl` | Add 3 new functions |

## Existing Functions Reused (unchanged)

| Function | Location | Used in |
|----------|----------|---------|
| `train_probit_model_df` / `predict_probit_model_df` | scoring_interface.jl:324, 366 | Phase 2, 3, 4 |
| `train_lightgbm_model_df` / `lightgbm_predict` | scoring_interface.jl:382, lightgbm_utils.jl | Phase 2, 3, 4 |
| `select_mbr_features` | scoring_interface.jl:454 | Phase 2 |
| `stream_sorted_merge` | MergeOperations.jl:419 | Phase 3 |
| `PSMFileReference` / `mark_sorted!` / `file_path` | FileReferences.jl, SortStateManagement.jl | Phase 3 |
| `write_arrow_file` | ArrowOperations.jl:161 | Phase 4 |
| `writeArrow` | writeArrow.jl:38 | Phase 3 (temp files) |

## Temp File Cleanup Pattern

Use the exact same pattern as `_compute_prob_threshold_from_merged` (scoring_workspace.jl:726-728):

```julia
tbl = nothing              # release Arrow mmap reference
GC.gc(false)               # minor GC to release file handles
rm(merged_path, force=true) # delete temp file
```

And for per-file temp refs, use the pattern from scoring_workspace.jl:867:
```julia
for ref in temp_refs
    rm(file_path(ref), force=true)
end
```

These are established patterns in the codebase, including Windows-safe handling.

---

## Architecture: 4-Phase Streaming

```
Phase 1: SAMPLE — Stream files one-at-a-time, reservoir-sample candidates for training
Phase 2: TRAIN  — Train per-fold Probit + LightGBM models on sample (in-memory, bounded)
Phase 3: EVAL   — Stream files to predict ALL candidates → sort-merge-stream FTR → select method
Phase 4: APPLY  — Stream files to re-predict + apply threshold + aggregate precursor probs
```

Each phase loads **at most one file** at a time. Phases 1, 3, 4 each make one pass through all files.

---

## Step 1: Add `_aggregate_trace_to_precursor_probs!`

**File**: `scoring_interface.jl` (near line 505, before `get_quant_necessary_columns`)

This extracts the inline `transform!(groupby(...))` logic from ScoringSearch.jl into a reusable per-file helper.

```julia
"""
    _aggregate_trace_to_precursor_probs!(df, has_mbr)

Per-file Bayesian aggregation of trace-level → precursor-level probabilities.
Groups by (precursor_idx, ms_file_idx). Since ms_file_idx is constant within
a single file, this is effectively grouping by precursor_idx alone.
"""
function _aggregate_trace_to_precursor_probs!(df::DataFrame, has_mbr::Bool)
    prob_agg = p -> begin
        trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
        clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
    end
    if has_mbr && hasproperty(df, :MBR_boosted_trace_prob)
        transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
                   :MBR_boosted_trace_prob => prob_agg => :MBR_boosted_prec_prob)
    end
    transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => prob_agg => :prec_prob)
end
```

---

## Step 2: Add `_compute_ftr_threshold_streaming`

**File**: `scoring_interface.jl`

Same pattern as `_compute_prob_threshold_from_merged` (scoring_workspace.jl:705-730), but computing FTR instead of FDR. Uses the same `tbl = nothing; GC.gc(false); rm(...)` cleanup pattern.

```julia
"""
    _compute_ftr_threshold_streaming(merged_path, target_ftr) -> (threshold, n_passing)

Stream through a merged Arrow file (sorted by score desc) computing FTR.
Returns the minimum score where FTR ≤ target_ftr and the count of candidates
passing at that threshold. Cleans up the merged file after reading.

Same pattern as _compute_prob_threshold_from_merged (scoring_workspace.jl:705).
"""
function _compute_ftr_threshold_streaming(merged_path::String, target_ftr::Float64)
    tbl = Arrow.Table(merged_path)
    scores = tbl.score
    bad_flags = tbl.is_bad_transfer
    n = length(scores)

    n_total = 0
    n_bad = 0
    threshold = typemax(Float64)
    n_passing = 0

    @inbounds for i in 1:n
        n_total += 1
        if bad_flags[i]
            n_bad += 1
        end
        if n_total > 0 && n_bad / n_total <= target_ftr
            threshold = Float64(scores[i])
            n_passing = n_total
        end
    end

    # Release mmap before deleting (same pattern as scoring_workspace.jl:726-728)
    tbl = nothing
    GC.gc(false)
    rm(merged_path, force=true)
    return threshold, n_passing
end
```

---

## Step 3: Add `apply_mbr_filter_and_aggregate_per_file!`

**File**: `scoring_interface.jl` (after existing `apply_mbr_filter!` at ~line 151)

This is the main orchestrator. Only the training columns are ever collected into memory (bounded by `MAX_MBR_TRAINING_CANDIDATES`). All other work is per-file.

### Phase 1: SAMPLE

Stream through files. For each file, identify candidates (`MBR_transfer_candidate == true`), extract only the ~10 training columns. Reservoir-sample into a bounded DataFrame.

**Training columns** (only what `select_mbr_features` + model training actually use):
```julia
training_columns = [
    # Features used by select_mbr_features()
    :trace_prob, :irt_error, :ms1_ms2_rt_diff,
    :MBR_max_pair_prob, :MBR_best_irt_diff, :MBR_rv_coefficient,
    :MBR_log2_weight_ratio, :MBR_log2_explained_ratio,
    # For ThresholdFilter scoring
    :MBR_boosted_trace_prob,
    # For CV fold splitting
    :cv_fold,
    # For bad transfer labeling
    :target, :decoy, :MBR_is_best_decoy
]
```

~40 bytes/candidate × 100K max sample = **~4 MB** training data.

```julia
# Phase 1: Stream files, reservoir-sample candidates
sample_df = DataFrame()
n_candidates_total = 0
rng = MersenneTwister(1844)

for ref in refs
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
    n_cand = sum(mask)
    n_cand == 0 && continue
    n_candidates_total += n_cand

    cols_available = filter(c -> hasproperty(df, c), training_columns)
    file_candidates = df[mask, cols_available]

    if nrow(sample_df) == 0
        sample_df = file_candidates
    else
        append!(sample_df, file_candidates)
    end
    # Downsample if over budget
    if nrow(sample_df) > MAX_MBR_TRAINING_CANDIDATES
        keep = randperm(rng, nrow(sample_df))[1:MAX_MBR_TRAINING_CANDIDATES]
        sample_df = sample_df[sort(keep), :]
    end
end
```

### Phase 2: TRAIN

Train per-fold models on the bounded sample. Each model is keyed by **test fold** (i.e., model for fold X was trained on data where `cv_fold != X`):

```julia
# Compute bad transfer labels on sample
sample_labels = (
    (sample_df.target .& coalesce.(sample_df.MBR_is_best_decoy, false)) .|
    (sample_df.decoy .& .!coalesce.(sample_df.MBR_is_best_decoy, false))
)

feature_cols = select_mbr_features(sample_df)
feature_data = sample_df[:, feature_cols]
folds = sort(unique(sample_df.cv_fold))

# Train Probit: test_fold => (model, valid_cols)
probit_models = Dict{eltype(folds), Any}()
probit_valid_cols = Dict{eltype(folds), Vector{Symbol}}()
probit_ok = true
try
    for fold in folds
        train_mask = sample_df.cv_fold .!= fold
        model, vcols = train_probit_model_df(feature_data[train_mask, :], sample_labels[train_mask], params)
        probit_models[fold] = model
        probit_valid_cols[fold] = vcols
    end
catch e
    @user_warn "Probit training failed: $(typeof(e)) — $(e)"
    probit_ok = false
end

# Train LightGBM: test_fold => booster
lgbm_models = Dict{eltype(folds), Any}()
lgbm_ok = true
try
    for fold in folds
        train_mask = sample_df.cv_fold .!= fold
        lgbm_models[fold] = train_lightgbm_model_df(feature_data[train_mask, :], sample_labels[train_mask], params)
    end
catch e
    @user_warn "LightGBM training failed: $(typeof(e)) — $(e)"
    lgbm_ok = false
end

sample_df = nothing; feature_data = nothing  # free training data
```

### Phase 3: EVAL

For each file: predict ALL candidates with all methods, write per-file sorted temp Arrow files. Then for each method: `stream_sorted_merge` → `_compute_ftr_threshold_streaming` → threshold + n_passing. Select best method.

```julia
# Build list of active methods
method_names = String["Threshold"]
probit_ok && push!(method_names, "Probit")
lgbm_ok && push!(method_names, "LightGBM")
method_temp_refs = Dict(m => PSMFileReference[] for m in method_names)

# Single pass through all files — predict candidates with all methods
for ref in refs
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
    n_cand = sum(mask)
    n_cand == 0 && continue

    cand_indices = findall(mask)
    cand_features = df[cand_indices, feature_cols]

    is_bad = (
        (df.target[cand_indices] .& coalesce.(df.MBR_is_best_decoy[cand_indices], false)) .|
        (df.decoy[cand_indices] .& .!coalesce.(df.MBR_is_best_decoy[cand_indices], false))
    )

    for method_name in method_names
        scores = if method_name == "Threshold"
            Float64.(df.MBR_boosted_trace_prob[cand_indices])
        elseif method_name == "Probit"
            out = zeros(Float64, n_cand)
            for fold in folds
                fm = df.cv_fold[cand_indices] .== fold
                any(fm) || continue
                out[fm] = predict_probit_model_df(probit_models[fold], cand_features[fm, :], probit_valid_cols[fold])
            end
            out
        else  # LightGBM
            out = zeros(Float64, n_cand)
            for fold in folds
                fm = df.cv_fold[cand_indices] .== fold
                any(fm) || continue
                out[fm] = lightgbm_predict(lgbm_models[fold], cand_features[fm, :])
            end
            out
        end

        # Sort by score desc, write temp Arrow
        perm = sortperm(scores; rev=true)
        temp_path = tempname() * "_mbr_ftr_$(method_name).arrow"
        writeArrow(temp_path, DataFrame(score=scores[perm], is_bad_transfer=is_bad[perm]))
        tref = PSMFileReference(temp_path)
        mark_sorted!(tref, :score)
        push!(method_temp_refs[method_name], tref)
    end
end

# For each method: merge → stream FTR → threshold + n_passing
method_results = Dict{String, Tuple{Float64, Int}}()  # method => (threshold, n_passing)
for method_name in method_names
    trefs = method_temp_refs[method_name]
    if isempty(trefs)
        method_results[method_name] = (typemax(Float64), 0)
        continue
    end
    merged_path = tempname() * "_mbr_ftr_merged_$(method_name).arrow"
    stream_sorted_merge(trefs, merged_path, :score; reverse=true)
    for tref in trefs
        rm(file_path(tref), force=true)  # same pattern as scoring_workspace.jl:867
    end
    threshold, n_passing = _compute_ftr_threshold_streaming(
        merged_path, Float64(params.max_MBR_false_transfer_rate))
    method_results[method_name] = (threshold, n_passing)
end

# Select best method (most passing candidates)
best_method = method_names[argmax([method_results[m][2] for m in method_names])]
best_threshold, best_n_passing = method_results[best_method]

# Log (same format as current scoring_interface.jl:120-124)
@user_info "MBR Method Selection:"
for method_name in method_names
    threshold, n_pass = method_results[method_name]
    marker = method_name == best_method ? " ✓" : ""
    @user_info "  $(method_name): $(n_pass)/$(n_candidates_total) pass ($(round(100*n_pass/n_candidates_total, digits=1))%)$marker"
end
```

### Phase 4: APPLY

For each file: re-predict best method, apply threshold, zero out bad transfers, aggregate precursor probs.

```julia
for ref in refs
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    cand_mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
    n_cand = sum(cand_mask)

    if n_cand > 0
        cand_indices = findall(cand_mask)

        # Re-predict with best method (cheap: Probit=dot product, LightGBM=tree eval)
        if best_method == "Threshold"
            for row in cand_indices
                if Float64(df.MBR_boosted_trace_prob[row]) < best_threshold
                    df.MBR_boosted_trace_prob[row] = 0.0f0
                    df.trace_prob[row] = 0.0f0
                end
            end
        else
            cand_features = df[cand_indices, feature_cols]
            scores = zeros(Float64, n_cand)
            for fold in folds
                fm = df.cv_fold[cand_indices] .== fold
                any(fm) || continue
                if best_method == "Probit"
                    scores[fm] = predict_probit_model_df(
                        probit_models[fold], cand_features[fm, :], probit_valid_cols[fold])
                else
                    scores[fm] = lightgbm_predict(lgbm_models[fold], cand_features[fm, :])
                end
            end
            for (i, row) in enumerate(cand_indices)
                if scores[i] < best_threshold
                    df.MBR_boosted_trace_prob[row] = 0.0f0
                    df.trace_prob[row] = 0.0f0
                end
            end
        end

        # Zero out bad transfers regardless of score
        for row in cand_indices
            is_bad = (df.target[row] && coalesce(df.MBR_is_best_decoy[row], false)) ||
                     (df.decoy[row] && !coalesce(df.MBR_is_best_decoy[row], false))
            if is_bad
                df.MBR_boosted_trace_prob[row] = 0.0f0
                df.trace_prob[row] = 0.0f0
            end
        end

        # Diagnostic columns
        df[!, :MBR_candidate] = df.MBR_transfer_candidate
        df[!, :MBR_transfer_q_value] = Vector{Union{Missing, Float32}}(missing, nrow(df))
    else
        df[!, :MBR_candidate] = falses(nrow(df))
        df[!, :MBR_transfer_q_value] = Vector{Union{Missing, Float32}}(missing, nrow(df))
    end

    # Precursor aggregation (per-file, no merge needed)
    _aggregate_trace_to_precursor_probs!(df, true)
    write_arrow_file(ref, df)
end
```

---

## Step 4: Replace Step 2 in ScoringSearch.jl

### BEFORE (lines 459-505):

```julia
        # Step 2: Apply MBR filtering and calculate precursor probabilities
        step2_time = @elapsed begin
            merged_scores_path = joinpath(temp_folder, "merged_trace_scores.arrow")
            sort_file_by_keys!(second_pass_refs, :trace_prob, :target; reverse=[true, true])
            stream_sorted_merge(second_pass_refs, merged_scores_path, :trace_prob, :target;
                               reverse=[true, true])

            merged_df = DataFrame(Arrow.Table(merged_scores_path))
            sqrt_n_runs = floor(Int64, sqrt(length(getFilePaths(getMSData(search_context)))))

            if hasproperty(merged_df, :MBR_boosted_trace_prob)
                apply_mbr_filter!(merged_df, params)

                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :MBR_boosted_trace_prob => (p -> begin
                               trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               trace_prob = clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(trace_prob)
                           end) => :MBR_boosted_prec_prob)

                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :trace_prob => (p -> begin
                               trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               trace_prob = clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(trace_prob)
                           end) => :prec_prob)
            else
                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :trace_prob => (p -> begin
                               trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               trace_prob = clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(trace_prob)
                           end) => :prec_prob)
            end

            for (file_idx, ref) in zip(valid_file_indices, second_pass_refs)
                sub_df = merged_df[merged_df.ms_file_idx .== file_idx, :]
                write_arrow_file(ref, sub_df)
            end
        end
```

### AFTER:

```julia
        # Step 2: Apply MBR filtering and calculate precursor probabilities (per-file OOM)
        step2_time = @elapsed begin
            apply_mbr_filter_and_aggregate_per_file!(second_pass_refs, valid_file_indices, params)
        end
```

**Removed**: `sort_file_by_keys!`, `stream_sorted_merge` into `merged_scores_path`, `merged_df = DataFrame(...)`, `sqrt_n_runs`, inline aggregation, per-file write-back loop.

---

## Memory Comparison

| Component | Current | OOM |
|-----------|---------|-----|
| Step 2 merged_df | O(N_total × 30 cols) | **eliminated** |
| Training data | O(N_candidates × 30 cols) | O(min(100K, N_candidates) × 10 cols) ≈ 4 MB |
| Peak per-file | part of merged_df | O(largest_file) |
| Temp FTR Arrow files | N/A | O(N_candidates × 5 bytes × 3 methods) |
| **Total Step 2** | **~120 bytes × N_total** | **~4 MB sample + max(file) + 15 × N_candidates bytes temp disk** |

---

## MBR_transfer_q_value

Per-candidate FTR q-values require globally sorted scores. In this implementation, `MBR_transfer_q_value` is set to `missing` for all rows. This is safe:
- The column is diagnostic, not used for filtering (filtering is done by zeroing `trace_prob`)
- The current code wraps the computation in try-catch (scoring_interface.jl:127-147)
- `get_quant_necessary_columns` expects the column to exist (it does, as `Union{Missing, Float32}`)

Can add streaming q-value computation in a follow-up.

## What Stays the Same

- All existing `apply_mbr_filter!` functions (kept for in-memory path/testing)
- All `train_and_evaluate`, `run_cv_training` functions
- `select_mbr_features`, `calibrate_ml_threshold`, `get_ftr_threshold`
- All downstream steps (3-23) — read from per-file refs with same columns

## Verification

1. `SearchDIA("./data/ecoli_test/ecoli_test_params.json")` — exercises OOM path
2. Verify "MBR Method Selection:" log output appears
3. Verify downstream columns: `MBR_boosted_prec_prob`, `prec_prob`, `MBR_candidate`, `MBR_transfer_q_value`
4. Verify no regression in precursor/protein identification counts
5. Check ecoli test's "No MBR transfer candidates found" warning still appears
