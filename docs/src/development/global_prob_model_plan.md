# Global Precursor Probability Model — Design and Implementation Plan

## Goal

Replace the current heuristic `:global_prob` aggregation (log-odds of per-run `:prec_prob`) with a dedicated ML model trained at the precursor level. The model predicts a single experiment‑wide probability per `precursor_idx` using summary features aggregated across runs. The predicted probability becomes the new `:global_prob` and is assigned back to all PSM rows for that `precursor_idx` in `merged_df` and downstream artifacts.

## Rationale

- The present `logodds(p, sqrt_n_runs)` is simple and robust but cannot exploit richer cross‑run signals (e.g., score dispersion, coverage, consistency).
- A small tree‑based model (LightGBM, same backend used elsewhere) can learn non‑linearities and interactions from a compact, per‑precursor feature vector, improving separation of targets and decoys at the precursor level.

## Scope

- In scope: New aggregation machinery, training pipeline, cross‑validated fit, prediction to all `precursor_idx`, and integration into ScoringSearch’s pipeline to populate `:global_prob`.
- Out of scope: Changes to per‑run `:prec_prob` estimation, downstream q‑value/PEP logic beyond switching `:global_prob` source, model persistence across projects.

## Data and Features

### Source tables

- `merged_df` after Step 2/3 of ScoringSearch: has one row per `(precursor_idx, ms_file_idx)` with columns at minimum `:precursor_idx`, `:prec_prob`, `:target`, possibly `:trace_prob`, etc.

### Aggregation object (custom type)

Define a compact summary per `precursor_idx`. We will construct a Dictionary keyed by `precursor_idx` whose values are instances of a small immutable struct:

```julia
struct GlobalPrecFeatures{N}
    top_scores::NTuple{N, Float32}      # descending top-N :prec_prob across runs (pad with -1.0f0)
    n_scores_available::Int32            # actual number of runs with scores (≤ N)
    stats::NTuple{5, Float32}           # (mean, max, min, std, skewness) - use -1.0f0 for undefined
    counts::NTuple{3, Int32}            # (n_runs_with_score, n_runs_total, n_above_thresh)
    deltas::NTuple{2, Float32}          # (top1 - top2, top2 - top3) - use -1.0f0 if undefined
    logodds_baseline::Float32           # current baseline logodds value (can be used as feature or standalone score)
end
```

Notes:
- `N` (default 5) controls the number of highest scores retained; pad missing with -1.0f0 sentinel value.
- `n_scores_available` tracks the actual number of runs where this precursor was observed.
- `n_runs_total` equals the number of MS files considered for the experiment.
- `n_above_thresh` uses a configurable `prec_prob` threshold (e.g., 0.95) to capture run coverage quality.
- `skewness` computed using streaming statistics (no array storage required); set to -1.0f0 if undefined.
- Delta features set to -1.0f0 when fewer scores available than needed for comparison.
- `logodds_baseline` equals the existing `logodds(p, sqrt_n_runs)` result, which returns a probability in [0,1] after applying sigmoid to average log-odds. This can serve both as a feature in the model and as the baseline score for comparison.

### Label assignment

- Precursor label `target::Bool` taken from library metadata (same across runs) via existing lookup on `precursor_idx`.

### Feature engineering details

- Sorting and padding for `top_scores` guarantee fixed‑length inputs for the model.
- Use -1.0f0 sentinel value consistently for all undefined/missing values (not 0.0 or NaN).
- Memory-efficient streaming calculations: mean, std, skewness computed without storing arrays.
- Standardization is not required for tree‑based models, but we will clamp/clip probabilities into `[eps, 1-eps]` when derived metrics use logs (e.g., baseline log‑odds).

## Training Pipeline

### Dataset construction

1. Group `merged_df` by `:precursor_idx` and collect `:prec_prob` values (one per run where observed).
   - Note: `merged_df` is memory-mapped, allowing efficient full scan without loading entire table into memory.
   - **Streaming approach**: For each group, collect only the top-N scores (sorted descending) and compute streaming statistics (mean, std, skewness) on-the-fly without materializing a full in-memory array of all scores across all runs. This is critical for scalability when datasets have many runs.
2. Build `GlobalPrecFeatures{N}` per group using streaming statistics and record `target`.
3. Convert the dictionary to a unique DataFrame with one row per `precursor_idx` and columns:
   - `precursor_idx::UInt32`, feature columns expanded from the struct, and `target::Bool`.

**Important**: Training happens once on the aggregated global features table, not separately for each run.

### Cross-validation (match MBR/regular scoring)

- Use Pioneer's existing CV folds supplied by the spectral library: `getCvFold(precursors, precursor_idx)` (from `StandardLibraryPrecursors.pid_to_cv_fold`).
- This mapping ensures all products of a given peptide (and all precursors belonging to the same protein group) end up in the same CV fold, consistent with MBR and regular scoring.
- Typical folds are two values (0 and 1). Perform CV by holding out one fold and training on the remaining fold(s), rotating across all unique fold values.
- No stratification required - keep the library-defined folds as-is.

### Model choice

- LightGBM only (binary logistic objective with probability output). No XGBoost.

### Hyperparameters (match MBR scoring; not user-configured)

- num_leaves: 63
- max_depth: 10
- learning_rate: 0.15
- min_data_in_leaf: 1
- feature_fraction: 0.5
- bagging_fraction: 0.5 (with `bagging_freq = 1`)
- min_gain_to_split: 0.0
- num_iterations: 200 (aligned with the final LightGBM round used in MBR)
- objective: binary (logistic), outputs calibrated probabilities
- metrics tracked: AUC, binary logloss (diagnostics only)

### Calibration

- Tree outputs are calibrated probabilities under logistic loss; optionally add isotonic calibration on out‑of‑fold predictions if needed. Start without extra calibration for simplicity.

### Evaluation and Fallback Logic

- Report CV AUC and logloss.
- Compare ROC/AUC of model vs. baseline `logodds_baseline`.
- **Q-value computation for method selection**:
  - Both model predictions and baseline logodds scores are probabilities in [0,1]
  - Assign chosen score to `:global_prob` column in `merged_df`
  - Use existing q-value pipeline: `stream_sorted_merge` (sorts by `:global_prob` and `:target`) followed by `get_precursor_global_qval_spline` (computes q-values from sorted scores)
- **Fallback decision**: After computing both model-based and baseline global q-values:
  - Count precursors passing q-value threshold with each method
  - Use whichever method passes more precursors (typically indicates better separation)
  - Log comparison with @info showing: model_passing, baseline_passing, method_selected
- **Minimum data requirements**:
  - At least 1 target AND 1 decoy per CV fold required for training
  - If requirements not met, automatically fallback to baseline logodds method
- Sanity check distributions on held‑out fold.

## Integration Points

Target file: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

1. After current computation of per‑run `:prec_prob` and before the current `:global_prob` calculation:
   - Build the per‑precursor dictionary of `GlobalPrecFeatures{N}` using streaming calculations.
   - Check minimum data requirements (≥1 target and ≥1 decoy per fold).
   - Train the model with CV and generate out‑of‑fold predictions for training diagnostics.
   - Fit final model on all data (single training on global features).
2. Generate per‑precursor predictions `global_prob_pred[precursor_idx]` using LightGBM.
3. Compute q-values for both model predictions and baseline logodds.
4. Compare performance and select the better method (more precursors passing threshold).
5. Assign chosen `:global_prob` to every row of `merged_df` by vectorized lookup.
6. Continue existing pipeline steps that consume `:global_prob` (e.g., global q‑value spline and filtering).

## Configuration

- Keep the feature behind a simple enable flag (e.g., `optimization.machine_learning.global_prob_model.enabled`; default false) and a small set of feature controls only (e.g., `top_n`, `prec_prob_threshold`).
- Do not expose model hyperparameters in JSON; they are fixed to the MBR settings listed above.

## File/Type Additions

- **New type**: `GlobalPrecFeatures{N}` in `src/structs/GlobalProb.jl` (will be loaded via `importScripts()` by adding to the structs list in `importScripts.jl`).
- **New package dependency**: `OnlineStats` (for streaming mean/std/skewness calculations) — add `using OnlineStats` in `src/Pioneer.jl`.
- New utils: `build_global_prec_features(merged_df; top_n, prob_thresh, n_runs)` to compute the dictionary and DataFrame using streaming statistics.
- New ML helper: `train_global_prob_model(df, folds, params)` returning fitted model and OOF predictions.
- New scorer: `predict_global_prob(model, feature_df) -> Dict{UInt32, Float32}` mapping `precursor_idx` to probability.
- New evaluator: `compare_global_prob_methods(model_scores, baseline_scores, merged_df, params, ctx)` to compute q-values for both methods and select the better one.

## Pseudocode Sketch

```julia
# After prec_prob computed per (precursor_idx, ms_file_idx)
merged_df = DataFrame(Arrow.Table(merged_scores_path))  # memory-mapped
n_runs = length(getFilePaths(getMSData(ctx)))
sqrt_n_runs = floor(Int, sqrt(n_runs))

dict = Dict{UInt32, GlobalPrecFeatures{TOPN}}()
labels = Dict{UInt32, Bool}()

# Process with streaming statistics (minimize in-memory arrays)
for g in groupby(merged_df, :precursor_idx)
    pid = first(g.precursor_idx)

    # Extract probabilities and determine top-N
    probs_view = g.prec_prob  # view, not a copy
    n_available = length(probs_view)

    # For top_scores: materialize only top-N by partial sorting
    top_n_vals = partialsort(probs_view, 1:min(n_available, TOPN); rev=true)
    top_scores = ntuple(i -> i <= length(top_n_vals) ? Float32(top_n_vals[i]) : -1.0f0, TOPN)

    # Streaming statistics using OnlineStats (single pass, no full array)
    stats_series = Series(Mean(), Variance(), Skewness())
    for p in probs_view
        fit!(stats_series, p)
    end

    if n_available > 0
        meanp = Float32(value(stats_series[1]))
        maxp = Float32(maximum(probs_view))
        minp = Float32(minimum(probs_view))
        stdp = n_available > 1 ? Float32(sqrt(value(stats_series[2]))) : -1.0f0
        skewp = n_available >= 3 ? Float32(value(stats_series[3])) : -1.0f0
    else
        meanp = maxp = minp = stdp = skewp = -1.0f0
    end

    nabove = count(>=(prob_thresh), probs_view)

    # Delta features with sentinel values
    del12 = (n_available >= 2) ? (top_scores[1] - top_scores[2]) : -1.0f0
    del23 = (n_available >= 3) ? (top_scores[2] - top_scores[3]) : -1.0f0

    # Baseline: logodds returns a probability (sigmoid of avg log-odds)
    lodds_baseline = n_available > 0 ? logodds(probs_view, sqrt_n_runs) : -1.0f0

    dict[pid] = GlobalPrecFeatures{TOPN}(
        top_scores,
        Int32(n_available),
        (meanp, maxp, minp, stdp, skewp),
        (Int32(n_available), Int32(n_runs), Int32(nabove)),
        (del12, del23),
        Float32(lodds_baseline),
    )
    labels[pid] = library_is_target(pid)
end

feat_df = features_to_dataframe(dict, labels)  # 1 row per pid
folds = [getCvFold(getPrecursors(getSpecLib(ctx)), pid) for pid in feat_df.precursor_idx]

# Check minimum data requirements
fold_counts = countmap(zip(folds, feat_df.target))
has_min_data = all(f ->
    get(fold_counts, (f, true), 0) >= 1 && get(fold_counts, (f, false), 0) >= 1,
    unique(folds)
)

if has_min_data
    model, oof_pred = train_global_prob_model(
        feat_df,
        folds;
        num_iterations=200,
        learning_rate=0.15,
        num_leaves=63,
        max_depth=10,
        feature_fraction=0.5,
        bagging_fraction=0.5,
        min_data_in_leaf=1,
        min_gain_to_split=0.0,
    )
    model_scores = predict_global_prob(model, feat_df)  # Dict{UInt32, Float32}
    baseline_scores = Dict(pid => feat.logodds_baseline for (pid, feat) in dict)

    # Compare methods: compute q-values using existing pipeline for each method
    # (stream_sorted_merge + get_precursor_global_qval_spline)
    use_model, model_passing, baseline_passing = compare_global_prob_methods(
        model_scores, baseline_scores, merged_df, params, ctx
    )

    global_prob_map = use_model ? model_scores : baseline_scores

    @info "Global prob method selection" model_passing baseline_passing use_model
else
    # Fallback to baseline if insufficient data
    global_prob_map = Dict(pid => feat.logodds_baseline for (pid, feat) in dict)
    @info "Insufficient data for model training, using baseline logodds"
end

# Assign chosen :global_prob to merged_df
merged_df.global_prob = [global_prob_map[pid] for pid in merged_df.precursor_idx]

# Continue with existing pipeline steps 5-6:
# Step 5: stream_sorted_merge on :global_prob and :target
# Step 6: get_precursor_global_qval_spline to compute final q-values
```

## Outputs and Artifacts

- Save training diagnostics to results folder: CV metrics, ROC curves, feature importances, and OOF prediction histogram.
- Optionally, persist the fitted model parameters (JSON) for reproducibility.

## Testing Plan

- Unit tests for feature builder: correct top‑N padding with -1.0f0 sentinel values, stats computation, counts, baseline.
- Test handling of edge cases: single-run precursors, missing values, undefined statistics.
- Small synthetic dataset: verify CV split stability and no leakage across groups.
- Integration test: run ScoringSearch on test data with the feature flag enabled; check that `:global_prob` is replaced and that q‑value distributions remain sensible.
- Regression: compare target/decoy separation AUC vs. baseline log‑odds on the same data.
- Verify method selection logic: ensure more passing precursors correctly trigger method selection.

## Risks and Mitigations

- Leakage across folds via related precursors: Use protein or peptide grouping to partition folds (already handled by library CV folds).
- Class imbalance: Use class weights or scale_pos_weight if needed; track `fdr_scale_factor` for consistency with q‑value estimation.
- Insufficient data: Require minimum 1 target + 1 decoy per fold; automatic fallback to baseline.
- Dependency sprawl: Use LightGBM only (already present); no new dependencies.
- Overfitting on small datasets: Keep feature set compact; compare against baseline performance.
- Memory usage with many files: Use streaming statistics and memory-mapped Arrow tables.

## Rollout

1. Implement behind `optimization.machine_learning.global_prob_model.enabled` (default false).
2. Land feature with comprehensive tests and docs; keep old log‑odds as fallback.
3. Evaluate on public and internal datasets; if improvements are consistent, flip default.
