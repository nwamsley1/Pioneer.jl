# ScoringSearch: Two‑Pass Probit Regression (Design Plan)

Author: Pioneer.jl coding assistant  
Status: Proposal (for review)  
Target: `src/Routines/SearchDIA/SearchMethods/ScoringSearch`

## Background

- FirstPassSearch currently uses an iterative probit scheme (`score_main_search_psms!` in `FirstPassSearch/utils.jl`).
  - Round 1: Train probit on all training data (initial q-values seeded from `:scribe`).
  - Later rounds: Re‑train probit on a subset: all decoys + high‑confidence targets (q ≤ threshold), then rescore and update q-values.

- ScoringSearch probit path (`probit_regression_scoring_cv!` in `score_psms.jl`) is single‑pass per CV fold:
  - For each `:cv_fold`, fit once with all train labels and predict on the held‑out fold.
  - Adds MBR compatibility columns but does not use MBR features.

Goal: Add a two‑pass (two‑round) probit flow to ScoringSearch that mirrors FirstPassSearch’s “subset of targets as true positives in the second round”.

## Scope and Non‑Goals

- In scope:
  - Add a second training pass per CV fold, based on q‑value filtering of the first pass predictions.
  - Keep intercept handling consistent (explicit `:intercept` feature and DataFrame column).
  - Preserve current MBR behavior for probit (compatibility columns only; no MBR features in training).
  - Reuse the existing ScoringSearch q‑value threshold configuration (see “Parameters”).

- Out of scope:
  - Changing XGBoost behavior or its iteration scheme.
  - Adding MBR features to probit training.
  - Changing how CV folds are assigned (still expect `[0,1]`).

## Parameters

- Threshold for selecting positive targets in pass 2:
  - Use a hardcoded q‑value threshold of 1% (0.01) for probit second pass selection of high‑confidence positives.
  - We deliberately do not reuse the XGBoost q‑value threshold here.

- FDR scale factor for q‑value computation:
  - Preferred: compute using the library target/decoy ratio (akin to `getLibraryFdrScaleFactor(search_context)`).
  - Practical in ScoringSearch: we don’t have `search_context` in the probit helpers. Plan: compute a fallback scale factor from the PSMs in scope (ratio of targets/decoys), with a default of `1.0f0`. If desired, we can plumb the library scale factor from the higher‑level call.

## Data and Features

- Intercept: already enforced in probit path; the feature set for probit models includes `:intercept`, and the column is inserted if missing. We will preserve this behavior.
- MBR: no MBR features used for training. After scoring, we continue to add `:MBR_prob` (copy of `:prob`) and `:MBR_is_best_decoy` (missing) for downstream compatibility.

## Algorithm: Two‑Pass Probit per CV Fold

For each `fold ∈ unique(psms.cv_fold)`:

1) Split data
   - `train_mask = (cv_fold != fold)`; `test_mask = (cv_fold == fold)`.
   - Validate minimum train data (≥100 targets and ≥100 decoys). If insufficient, set `prob[test_mask] = 0.5` and continue.

2) Pass 1 — Initial fit and predict
   - Ensure `:intercept` column exists; ensure `:intercept` is in features vector (dedup if needed).
   - Fit `β₁ = Pioneer.ProbitRegression` on `train_data[:, features]` vs `train_data.target`.
   - Predict `prob₁[test_mask] = Pioneer.ModelPredictProbs!(...)`.

3) Compute q‑values (for selection) and PEP (for negative mining)
   - Build temporary vectors on the training set: `qvals` and `pep`.
   - Compute `get_qvalues!(prob₁_train, target_train, qvals)` and `get_PEP!(prob₁_train, target_train, pep)`.
   - Selection masks: high‑confidence positives via `qvals ≤ 0.01`; mined negatives via `pep ≥ neg_mining_pep_threshold`.

4) Pass 2 — Subset re‑fit and predict
   - Build masks (on the training set):
     - `decoys_train = !target`
     - `confident_pos = target & (qvals ≤ 0.01)`
     - `mined_negs = target & (pep ≥ neg_mining_pep_threshold)`
   - `second_mask = decoys_train | confident_pos | mined_negs` and `second_labels = false` except `true` for `confident_pos` entries.
   - Fit `β₂` on `train_data[second_mask, features]` vs `second_labels`.
   - Predict `prob₂[test_mask]` and overwrite `prob[test_mask] = prob₂[test_mask]`.

5) Finalization
   - After looping all folds: clamp probs to `[1e-6, 0.9999]`, add MBR compatibility columns, and drop temporary columns (e.g., `:intercept`, `:cv_fold`).

Notes:
- This mirrors FirstPassSearch’s `score_main_search_psms!`: first broad fit, then restricted positives in the second fit.
- We do not iterate beyond 2 passes to keep runtime predictable and behavior aligned with the request.

## Integration Points and API Changes

1) `train_probit_model_in_memory`
   - Current: `train_probit_model_in_memory(best_psms, file_paths, precursors, model_config, match_between_runs)`
   - Change: add `rescore_q_threshold::Float32` and optional `fdr_scale_factor::Float32=1.0f0`
   - Caller: `score_precursor_isotope_traces_in_memory(...)` has `max_q_value_xgboost_rescore` available; pass it as `rescore_q_threshold`. For `fdr_scale_factor`, either:
     - compute from `best_psms` (ratio targets/decoys), or
     - plumb a library‑derived value from upstream in a follow‑up change.

2) `probit_regression_scoring_cv!`
   - Current: `(psms, file_paths, features, match_between_runs)`
   - Change: add keyword `neg_mining_pep_threshold::Float32` (default 0.90) and implement a two‑pass flow.
   - Second pass uses hardcoded q ≤ 0.01 for selecting positive targets and converts low‑confidence targets with `PEP ≥ neg_mining_pep_threshold` into negatives for training.

3) Model Config
   - No additional change required; probit models already include `:intercept` as a feature. XGBoost models do not.

## Edge Cases & Fallbacks

- Insufficient training data per fold (targets < 100 or decoys < 100): set `prob[test_mask] = 0.5` and skip to next fold.
- No/constant features: filter features by `hasproperty(psms, f)`; if the set becomes empty, warn and set 0.5.
- Numerical issues (NaN/Inf): detected and replaced with 0.5; clamp probability range to avoid downstream overflow.
- `:cv_fold` not in data: warn and treat all as a single fold (train/test split by random mask or use full data with K=1) — optional enhancement; initially require `:cv_fold`.

## Testing Plan

- Unit tests:
  - Intercept: probit models include `:intercept` in feature list; XGBoost models do not.
  - Two‑pass behavior: with a small synthetic dataset, show that pass‑2 rescoring changes some probabilities when threshold is tight.
  - Fallbacks: insufficient data per fold returns 0.5 on the fold; no NaN/Inf after clamping.

- Integration tests:
  - Run model comparison on a small test set; confirm both probit configs run and do not regress existing outputs.
  - Validate that downstream MBR filtering still operates with the compatibility columns.

## Acceptance Criteria

- Probit path runs two passes per CV fold and uses the subset of high‑confidence targets (q ≤ threshold) plus all decoys for pass 2.
- No change to XGBoost behavior or parameter structures.
- Intercept is present in probit features and added to the DataFrame when missing; it is not present for XGBoost.
- Tests covering model list size and intercept presence pass; new tests for two‑pass logic added and passing.

## Rollout

- Implement behind behavior change default (enabled). If issues arise, we can add a `:two_pass_probit` toggle in `ModelConfig.hyperparams` to disable it temporarily.

---

References:
- First pass iterative probit: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl` (`score_main_search_psms!`).
- Current ScoringSearch probit: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl` (`probit_regression_scoring_cv!`).
