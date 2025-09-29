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
    top_scores::NTuple{N, Float32}  # descending top-N :prec_prob across runs (pad with 0.0f0)
    stats::NTuple{5, Float32}       # (mean, max, min, std, iqr)
    counts::NTuple{3, Int32}        # (n_runs_with_score, n_runs_total, n_above_thresh)
    deltas::NTuple{2, Float32}      # (top1 - top2, top2 - top3)
    logodds_baseline::Float32       # current baseline (for ensembling/ablation)
end
```

Notes:
- `N` (default 5) controls the number of highest scores retained; pad missing with zeros.
- `n_runs_total` equals the number of MS files considered for the experiment.
- `n_above_thresh` uses a configurable `prec_prob` threshold (e.g., 0.95) to capture run coverage quality.
- `iqr` computed from observed per‑run scores for robustness; if <4 observations, fallback to 0.
- `logodds_baseline` equals the existing `logodds(p, sqrt_n_runs)` result for a comparable reference feature.

### Label assignment

- Precursor label `target::Bool` taken from library metadata (same across runs) via existing lookup on `precursor_idx`.

### Feature engineering details

- Sorting and padding for `top_scores` guarantee fixed‑length inputs for the model.
- Standardization is not required for tree‑based models, but we will clamp/clip probabilities into `[eps, 1-eps]` when derived metrics use logs (e.g., baseline log‑odds).

## Training Pipeline

### Dataset construction

1. Group `merged_df` by `:precursor_idx` and collect `:prec_prob` values (one per run where observed).
2. Build `GlobalPrecFeatures{N}` per group and record `target`.
3. Convert the dictionary to a unique DataFrame with one row per `precursor_idx` and columns:
   - `precursor_idx::UInt32`, feature columns expanded from the struct, and `target::Bool`.

### Cross-validation (match MBR/regular scoring)

- Use Pioneer’s existing CV folds supplied by the spectral library: `getCvFold(precursors, precursor_idx)` (from `StandardLibraryPrecursors.pid_to_cv_fold`).
- This mapping ensures all products of a given peptide (and all precursors belonging to the same protein group) end up in the same CV fold, consistent with MBR and regular scoring.
- Typical folds are two values (0 and 1). Perform CV by holding out one fold and training on the remaining fold(s), rotating across all unique fold values.

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

### Evaluation

- Report CV AUC and logloss.
- Compare ROC/AUC of model vs. baseline `logodds_baseline` to justify replacement.
- Sanity check distributions on held‑out fold.

## Integration Points

Target file: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

1. After current computation of per‑run `:prec_prob` and before the current `:global_prob` calculation:
   - Build the per‑precursor dictionary of `GlobalPrecFeatures{N}`.
   - Train the model with CV and generate out‑of‑fold predictions for training diagnostics.
   - Fit final model on all data.
2. Generate per‑precursor predictions `global_prob_pred[precursor_idx]` using LightGBM trained with the library‑defined folds.
3. Assign `:global_prob = global_prob_pred[precursor_idx]` to every row of `merged_df` by a single vectorized lookup or join.
4. Continue existing pipeline steps that consume `:global_prob` (e.g., global q‑value spline and filtering).

## Configuration

- Keep the feature behind a simple enable flag (e.g., `optimization.machine_learning.global_prob_model.enabled`; default false) and a small set of feature controls only (e.g., `top_n`, `prec_prob_threshold`).
- Do not expose model hyperparameters in JSON; they are fixed to the MBR settings listed above.

## File/Type Additions

- New type: `GlobalPrecFeatures{N}` in an appropriate module (e.g., `src/structs/GlobalProb.jl`).
- New utils: `build_global_prec_features(merged_df; top_n, prob_thresh, n_runs)` to compute the dictionary and DataFrame.
- New ML helper: `train_global_prob_model(df, folds, params)` returning fitted model and OOF predictions.
- New scorer: `predict_global_prob(model, feature_df) -> Dict{UInt32, Float32}` mapping `precursor_idx` to probability.

## Pseudocode Sketch

```julia
# After prec_prob computed per (precursor_idx, ms_file_idx)
merged_df = DataFrame(Arrow.Table(merged_scores_path))
n_runs = length(getFilePaths(getMSData(ctx)))

dict = Dict{UInt32, GlobalPrecFeatures{TOPN}}()
labels = Dict{UInt32, Bool}()
for g in groupby(merged_df, :precursor_idx)
    pid = first(g.precursor_idx)
    probs = sort!(Float32.(g.prec_prob); rev=true)
    top_scores = ntuple(i -> i <= length(probs) ? probs[i] : 0.0f0, TOPN)
    meanp = mean(probs); maxp = maximum(probs); minp = minimum(probs)
    stdp = length(probs) > 1 ? std(probs) : 0f0
    iqrp = length(probs) >= 4 ? (quantile(probs, 0.75f0) - quantile(probs, 0.25f0)) : 0f0
    nabove = count(>=(prob_thresh), probs)
    del12 = (top_scores[1] - top_scores[2])
    del23 = (top_scores[2] - top_scores[3])
    lodds = logodds(probs, floor(Int, sqrt(n_runs)))  # baseline
    dict[pid] = GlobalPrecFeatures{TOPN}(
        top_scores,
        (meanp, maxp, minp, stdp, iqrp),
        (length(probs), n_runs, nabove),
        (del12, del23),
        lodds,
    )
    labels[pid] = library_is_target(pid)
end

feat_df = features_to_dataframe(dict, labels)  # 1 row per pid
folds = [getCvFold(getPrecursors(getSpecLib(ctx)), pid) for pid in feat_df.precursor_idx]
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
global_prob_map = predict_global_prob(model, feat_df)

# Assign back to merged_df
merged_df.global_prob = [global_prob_map[pid] for pid in merged_df.precursor_idx]
```

## Outputs and Artifacts

- Save training diagnostics to results folder: CV metrics, ROC curves, feature importances, and OOF prediction histogram.
- Optionally, persist the fitted model parameters (JSON) for reproducibility.

## Testing Plan

- Unit tests for feature builder: correct top‑N padding, stats, counts, baseline.
- Small synthetic dataset: verify CV split stability and no leakage across groups.
- Integration test: run ScoringSearch on test data with the feature flag enabled; check that `:global_prob` is replaced and that q‑value distributions remain sensible.
- Regression: compare target/decoy separation AUC vs. baseline log‑odds on the same data.

## Risks and Mitigations

- Leakage across folds via related precursors: Use protein or peptide grouping to partition folds.
- Class imbalance: Use class weights or scale_pos_weight; also track `fdr_scale_factor` for consistency with q‑value estimation.
- Dependency sprawl: Use LightGBM only (already present); no new dependencies.
- Overfitting on small datasets: Use early stopping, CV monitoring, and keep feature set compact.

## Rollout

1. Implement behind `optimization.machine_learning.global_prob_model.enabled` (default false).
2. Land feature with comprehensive tests and docs; keep old log‑odds as fallback.
3. Evaluate on public and internal datasets; if improvements are consistent, flip default.
