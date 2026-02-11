# Percolator Scoring Workflow: Branch Comparison Report

**Branches compared:** `develop` vs `feature/oom_percolator`
**Date:** 2026-02-11
**Repository:** Pioneer.jl

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Background: What Percolator Scoring Does](#2-background)
3. [Develop Branch Architecture](#3-develop-branch-architecture)
4. [Feature/OOM Branch Architecture](#4-featureoom-branch-architecture)
5. [Detailed Comparison](#5-detailed-comparison)
6. [File Inventory](#6-file-inventory)
7. [Data Flow Diagrams](#7-data-flow-diagrams)
8. [Key Algorithms](#8-key-algorithms)
9. [Integration with ScoringSearch](#9-integration-with-scoringsearch)
10. [Trade-offs and Design Decisions](#10-trade-offs-and-design-decisions)

---

## 1. Executive Summary

The `feature/oom_percolator` branch is a major architectural refactoring of Pioneer's percolator-style PSM scoring system. The key changes are:

| Aspect | `develop` | `feature/oom_percolator` |
|--------|-----------|--------------------------|
| **Architecture** | Single 1,237-line monolithic function | Trait-based system across 13 focused files (2,697 lines total) |
| **ML files** | 11 files in `src/utils/ML/` | 23 files in `src/utils/ML/` |
| **OOM support** | Hardcoded `if false` (disabled) | Fully implemented via `ArrowFilePSMContainer` + sidecar pattern |
| **Extensibility** | Add code branches to monolithic function | Add new trait types; zero changes to core loop |
| **Entry point** | `sort_of_percolator_in_memory!()` (DataFrame-only) | `percolator_scoring!(psms::AbstractPSMContainer, config::ScoringConfig)` |
| **Core algorithm** | Identical semi-supervised CV training with iterative negative mining and MBR | Identical (preserved exactly) |

The refactoring **preserves the scoring algorithm exactly** while decomposing it into composable, independently testable trait implementations. The OOM path enables scoring datasets that exceed available RAM by keeping PSMs on disk as Arrow files and only loading subsets for training and per-file prediction.

---

## 2. Background

### What Percolator Scoring Does

Pioneer uses a semi-supervised machine learning approach (inspired by the Percolator algorithm) to rescore peptide-spectrum matches (PSMs). The workflow:

1. **Pair** target and decoy precursors randomly within retention-time bins
2. **Train** a classifier (LightGBM or probit) using cross-validation folds
3. **Iterate** with progressively harder negative mining (converting low-confidence targets to negatives)
4. **Optionally** compute match-between-runs (MBR) features from paired precursors across MS files
5. **Finalize** by computing q-values and identifying MBR transfer candidates

The output is a `trace_prob` score for each PSM, plus optional `MBR_boosted_trace_prob` and `MBR_transfer_candidate` columns for MBR-enabled searches.

---

## 3. Develop Branch Architecture

### Single-File Design

On develop, virtually all scoring logic lives in one file:

```
src/utils/ML/percolatorSortOf.jl   (1,237 lines)
```

Supporting files exist for LightGBM utilities, probit regression, and FDR computation, but the main scoring loop, pairing, MBR feature computation, training data selection, and finalization are all inline within `sort_of_percolator_in_memory!()`.

### Include Order (develop)

```julia
"fdrUtilities.jl",
"ftrUtilities.jl",
"lightgbm_utils.jl",
"percolatorSortOf.jl",       # <-- monolithic: pairing + training + MBR + finalization
"piecewiseLinearFunction.jl",
"probitRegression.jl",
"spectralLinearRegression.jl",
"uniformBasisCubicSpline.jl",
"wittakerHendersonSmoothing.jl",
"libraryBSpline.jl"
```

### Main Function Signature (develop)

```julia
function sort_of_percolator_in_memory!(
    psms::DataFrame,
    features::Vector{Symbol},
    match_between_runs::Bool = true;
    max_q_value_lightgbm_rescore::Float32 = 0.01f0,
    max_q_value_mbr_itr::Float32 = 0.20f0,
    min_PEP_neg_threshold_itr = 0.90f0,
    feature_fraction::Float64 = 0.5,
    learning_rate::Float64 = 0.15,
    min_data_in_leaf::Int = 1,
    bagging_fraction::Float64 = 0.5,
    min_gain_to_split::Float64 = 0.0,
    max_depth::Int = 10,
    num_leaves::Int = 63,
    iter_scheme::Vector{Int} = [100, 200, 200],
    show_progress::Bool = true,
    verbose_logging::Bool = false
)
```

All hyperparameters, iteration control, and behavioral switches are passed as keyword arguments. The function only accepts `DataFrame` input.

### Workflow (develop)

```
sort_of_percolator_in_memory!(psms, features, match_between_runs)
│
├── 1. assign_random_target_decoy_pairs!(psms)     [inline, ~120 lines]
│   ├── getIrtBins(irts)
│   └── assign_pair_ids()
│
├── 2. sort!(psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])
│
├── 3. Allocate output arrays: prob_test, prob_train, MBR_estimates, nonMBR_estimates
│
├── 4. Pre-compute fold_indices / train_indices via Dict comprehension
│
├── 5. For each CV fold:
│   ├── initialize_prob_group_features!(psms, match_between_runs)
│   ├── Split train/test as DataFrame views
│   │
│   └── For each iteration (itr, num_round) in iter_scheme:
│       ├── get_training_data_for_iteration!()     [inline, ~40 lines]
│       │   └── Itr 1: all data; Itr 2+: q-value filter + PEP negative mining
│       ├── train_booster()                         [inline helper, ~25 lines]
│       ├── predict on train → prob_train, compute q-values
│       ├── predict on test  → prob_test
│       ├── store nonMBR baseline at mbr_start_iter - 1
│       └── update_mbr_features!() → summarize_precursors!()  [~100 lines]
│
├── 6. Finalize:
│   ├── MBR path: compute transfer candidates, set trace_prob + MBR_boosted_trace_prob
│   └── Non-MBR path: set trace_prob = prob_test
│
└── Return models Dict{UInt8, LightGBMModelVector}
```

### OOM Support (develop)

The OOM path is **completely disabled** via a hardcoded `if false` guard in `score_psms.jl`:

```julia
if false  # psms_count >= max_psms_in_memory
    # Case 1: Out-of-memory processing with default LightGBM (DISABLED)
    ...
```

The develop branch had aspirational OOM code but it was never functional.

---

## 4. Feature/OOM Branch Architecture

### Trait-Based Decomposition

The monolithic function has been decomposed into 6 independent trait hierarchies combined via a `ScoringConfig` struct:

```
ScoringConfig{M, P, T, F, I, B}
├── model::M         ← PSMScoringModel      (LightGBMScorer | ProbitScorer)
├── pairing::P       ← PairingStrategy      (RandomPairing | NoPairing)
├── training_data::T ← TrainingDataStrategy  (QValueNegativeMining | AllDataSelection)
├── feature_selection::F ← FeatureSelectionStrategy (IterativeFeatureSelection | StaticFeatureSelection)
├── iteration_scheme::I  ← IterationScheme   (FixedIterationScheme | SinglePassScheme)
└── mbr_update::B    ← MBRUpdateStrategy    (PairBasedMBR | NoMBR)
```

Additionally, a 7th trait (`ScoringPhase`) separates training from prediction behavior:

```
ScoringPhase
├── TrainingPhase   — trains models, predicts on training data
└── PredictionPhase — applies pre-trained models to held-out test data
```

### Data Abstraction Layer

A new `AbstractPSMContainer` hierarchy decouples the scoring algorithm from data storage:

```
AbstractPSMContainer
├── DataFramePSMContainer       — wraps DataFrame (in-memory)
│   └── DataFramePSMContainerView — zero-copy view for train/test splits
└── ArrowFilePSMContainer       — file-backed (OOM)
    └── ArrowFileGroup          — one (ms_file_idx, cv_fold) Arrow file pair
```

### Workspace Abstraction

A new `AbstractScoringWorkspace` manages CV fold setup, output arrays, and model storage:

```
AbstractScoringWorkspace
├── InMemoryScoringWorkspace     — pre-computed fold indices + global arrays
└── ArrowFileScoringWorkspace    — wraps InMemoryScoringWorkspace for sampled training
                                   + per-file streaming prediction
```

### Include Order (feature/oom_percolator)

```julia
# Core ML utilities (unchanged)
"fdrUtilities.jl",
"ftrUtilities.jl",
"lightgbm_utils.jl",
"probitRegression.jl",
"piecewiseLinearFunction.jl",
"spectralLinearRegression.jl",
"uniformBasisCubicSpline.jl",
"wittakerHendersonSmoothing.jl",
"libraryBSpline.jl",

# Trait-based scoring system (NEW, in dependency order)
"psm_container.jl",         # AbstractPSMContainer abstraction
"arrow_psm_container.jl",   # ArrowFilePSMContainer (file-backed OOM)
"scoring_traits.jl",        # 6 trait types + ScoringPhase
"scoring_config.jl",        # ScoringConfig struct
"pairing.jl",               # PairingStrategy implementations
"model_training.jl",        # PSMScoringModel train/predict
"training_selection.jl",    # TrainingDataStrategy implementations
"feature_selection.jl",     # FeatureSelectionStrategy implementations
"iteration_scheme.jl",      # IterationScheme implementations
"mbr_update.jl",            # MBRUpdateStrategy implementations
"scoring_workspace.jl",     # Workspace + OOM sampling + OOM prediction
"percolator_generic.jl",    # Main percolator_scoring! function
"percolatorSortOf.jl"       # Legacy wrapper (now 350 lines, down from 1,237)
```

### Main Function Signature (feature/oom_percolator)

```julia
function percolator_scoring!(
    psms::AbstractPSMContainer,
    config::ScoringConfig;
    show_progress::Bool = true,
    verbose::Bool = false
) -> Dict{UInt8, Vector{Any}}
```

All algorithm configuration is encoded in the type-parameterized `ScoringConfig`. The function accepts any `AbstractPSMContainer`, enabling both in-memory and OOM execution via dispatch.

### Workflow (feature/oom_percolator)

```
percolator_scoring!(psms, config)
│
├── 1. prepare_training_data!(psms, config)          [dispatched on container type]
│   ├── DataFramePSMContainer:
│   │   ├── assign_pairs!(psms, config.pairing)
│   │   ├── sort_container!(psms, sort_cols)
│   │   └── initialize_mbr_columns!(psms, config.mbr_update)
│   │
│   └── ArrowFilePSMContainer:
│       ├── Phase 1: Read 4 lightweight columns globally, assign pairs (~19 bytes/PSM)
│       └── Phase 2: Per-file sort + sidecar creation (~32 bytes/PSM for mutable cols)
│
├── 2. setup_scoring_workspace(psms, config)          [dispatched on container type]
│   ├── DataFramePSMContainer → InMemoryScoringWorkspace
│   │   └── Pre-compute fold/train indices, allocate output arrays
│   │
│   └── ArrowFilePSMContainer → ArrowFileScoringWorkspace
│       └── 3-pass pair-sampled training:
│           1. Count PSMs per pair_id (read only pair_id column)
│           2. Shuffle pair_ids, greedily select within budget
│           3. Pre-allocate + fill DataFrame from matching rows
│
├── 3. PHASE 1 — Train all models
│   └── For each fold:
│       └── process_fold_iterations!(TrainingPhase(), workspace, fold, ...)
│           └── For each iteration:
│               ├── get_training_subset(TrainingPhase(), ...)  → filter by q-value/PEP
│               ├── get_or_train_model(TrainingPhase(), ...)   → train + store
│               ├── predict_scores(model, psms_view)
│               ├── compute_and_set_qvalues!(TrainingPhase(), ...) → for next iter
│               └── update_mbr!(TrainingPhase(), ...)
│
├── 4. PHASE 2 — Predict on held-out test data
│   └── For each fold:
│       ├── InMemoryScoringWorkspace:
│       │   └── process_fold_iterations!(PredictionPhase(), workspace, fold, ...)
│       │       └── Same loop, but retrieves stored models instead of training
│       │
│       └── ArrowFileScoringWorkspace:
│           └── process_fold_iterations!(PredictionPhase(), workspace, fold, ...)
│               └── SPECIALIZED: per-file streaming prediction
│                   ├── For each Arrow file:
│                   │   ├── Load data + scores sidecar
│                   │   ├── Filter to fold, apply model, write predictions to sidecar
│                   │   └── Store nonMBR baseline at pre-MBR iteration
│                   └── Cross-file MBR: _compute_mbr_for_fold!()
│                       └── Streaming 2-pass MBR with disk-based FDR (O(1) memory)
│
├── 5. store_final_predictions!(workspace, fold, mbr_strategy)
│
└── 6. finalize_scoring!(workspace, config.mbr_update)  [dispatched on workspace type]
    ├── InMemoryScoringWorkspace:
    │   ├── PairBasedMBR: compute q-values on baseline → find threshold → mark transfer candidates
    │   └── NoMBR: set trace_prob = prob_test
    │
    └── ArrowFileScoringWorkspace:
        ├── PairBasedMBR: 2-pass over all files
        │   ├── Pass 1: Collect baselines + targets → global q-values → probability threshold
        │   └── Pass 2: Merge trace_prob, MBR_boosted_trace_prob, MBR features into data files
        └── NoMBR: Copy trace_prob from sidecars to data files
```

---

## 5. Detailed Comparison

### 5.1 Pairing

| | `develop` | `feature/oom_percolator` |
|-|-----------|--------------------------|
| **Location** | Inline in `percolatorSortOf.jl` (~120 lines) | `pairing.jl` (136 lines) |
| **Input** | `DataFrame` only | `AbstractPSMContainer` |
| **Algorithm** | Identical: random pairs within iRT bins (seed=1844, bin_size=1000) | Identical |
| **Extensibility** | Requires modifying inline code | Add new `PairingStrategy` subtype |
| **OOM support** | N/A | Global pairing reads only 4 lightweight columns (~19 bytes/PSM) |

Both branches use the same deterministic pairing algorithm:
1. Compute iRT bins of size 1000 via `getIrtBins()`
2. Group by `(irt_bin, cv_fold, isotopes_captured)`
3. Within each group, shuffle `precursor_idx` with seed 1844
4. Pair consecutive elements: `(1,2), (3,4), ...`

### 5.2 Cross-Validation Setup

| | `develop` | `feature/oom_percolator` |
|-|-----------|--------------------------|
| **Fold storage** | `Dict{fold => indices}` in local variables | `InMemoryScoringWorkspace` struct |
| **Train/test split** | `@view psms[train_idx, :]` / `@view psms[test_idx, :]` | `DataFramePSMContainerView(parent, indices)` |
| **Output arrays** | Local `Float32` vectors | Workspace fields: `prob_test`, `prob_train`, `nonMBR_estimates`, `MBR_estimates` |
| **Model storage** | `Dict{UInt8, LightGBMModelVector}` | `Dict{UInt8, Vector{Any}}` in workspace |

The OOM branch adds the concept of workspace accessors dispatched on `ScoringPhase`:
```julia
get_phase_view(ws, TrainingPhase(), fold)   → training set view
get_phase_view(ws, PredictionPhase(), fold) → test set view
get_phase_output(ws, TrainingPhase())       → prob_train array
get_phase_output(ws, PredictionPhase())     → prob_test array
```

### 5.3 Training Loop

| | `develop` | `feature/oom_percolator` |
|-|-----------|--------------------------|
| **Loop structure** | Single nested loop: `for fold` → `for itr` | Two-phase: Phase 1 trains all folds, then Phase 2 predicts all folds |
| **Training data selection** | `get_training_data_for_iteration!()` inline | `select_training_data()` dispatched on `TrainingDataStrategy` |
| **Feature selection** | Manual: `itr < mbr_start_iter ? non_mbr_features : features` | `get_features()` dispatched on `FeatureSelectionStrategy` |
| **Model training** | `train_booster()` inline helper | `train_model()` dispatched on `PSMScoringModel` |
| **Prediction** | `predict(bst, psms_test)` | `predict_scores(model, psms_view)` |

**Critical difference in loop structure:** On develop, training and prediction happen in the same loop iteration — the model is trained on the training set, then immediately applied to the test set. On the OOM branch, Phase 1 trains all folds first (storing models), then Phase 2 applies stored models to test data. This separation is required for OOM because:
- Training happens on a sampled in-memory subset
- Prediction must stream through all files on disk

The scoring algorithm itself (negative mining, MBR iteration, q-value computation) is preserved identically.

### 5.4 Training Data Selection

Both branches implement the same iterative negative mining strategy:

- **Iteration 1:** Use all training data
- **Iterations 2+:**
  1. Compute PEP from current `trace_prob` scores
  2. Convert worst-scoring targets (PEP ≥ `min_PEP_neg_threshold_itr`) to negatives (flip `target` to `false`)
  3. Keep all decoys + targets with `q_value ≤ max_q_value_lightgbm_rescore`

On develop this is inline in `get_training_data_for_iteration!()`. On the OOM branch it's factored into `training_selection.jl` with dispatch on `QValueNegativeMining` vs `AllDataSelection`.

### 5.5 MBR Feature Computation

| | `develop` | `feature/oom_percolator` |
|-|-----------|--------------------------|
| **Core function** | `summarize_precursors!()` (inline, ~100 lines) | Same `summarize_precursors!()` (moved but identical) |
| **When called** | `update_mbr_features!()` handles both train + test | Split: `update_mbr_features_train_only!()` / `update_mbr_features_test_only!()` |
| **OOM MBR** | N/A | `_compute_mbr_for_fold!()`: loads fold data from ALL files, runs `summarize_precursors!()`, distributes results back |

The MBR algorithm computes 9 features per PSM within each `(pair_id, isotopes_captured)` group:

| Feature | Description |
|---------|-------------|
| `MBR_max_pair_prob` | Best probability from paired PSM in another run |
| `MBR_best_irt_diff` | Absolute iRT residual difference to best paired PSM |
| `MBR_rv_coefficient` | RV coefficient (multivariate correlation of weights × RTs) |
| `MBR_log2_weight_ratio` | log2(current weight / best weight) |
| `MBR_log2_explained_ratio` | log2(current explained intensity / best explained intensity) |
| `MBR_is_best_decoy` | Whether the best paired PSM is a decoy |
| `MBR_num_runs` | Number of other runs with passing matches in this pair |
| `MBR_transfer_candidate` | Computed during finalization: failed baseline but best pair prob ≥ threshold |
| `MBR_is_missing` | No valid paired PSM found |

The OOM MBR computation is the most complex part of the OOM path because pairs span multiple MS files. The solution uses a fully streaming 2-pass approach:

**Pass 1:** For each file, read `trace_prob` + `decoy` columns, build the pair dictionary (best 2 runs per pair), and write a per-file sorted temp Arrow file with `(trace_prob, target)`.

**Between passes:** Use `stream_sorted_merge` (existing k-way merge infrastructure) to globally sort all temp files by `trace_prob` descending, then stream through the merged file doing a single forward FDR pass to find the probability threshold — O(1) heap memory.

**Pass 2:** For each file, compute MBR features using the pair dictionary and write updated scores to sidecars.

**Key optimization:** Per-PSM q-values are never written to sidecars during MBR — they were unused (the ML model doesn't use `q_value` as a feature, and finalization recomputes q-values independently from `nonMBR_trace_prob`). This eliminates the `all_qvalues` vector (8 bytes/PSM) entirely. Combined with the disk-based sort, total memory for FDR computation dropped from N × 13 bytes to ~5-10 MB regardless of dataset size.

### 5.6 Finalization

Both branches use the same finalization logic:

**With MBR:**
1. Compute q-values on the non-MBR baseline (`nonMBR_estimates`)
2. Find the probability threshold where targets pass at q ≤ 0.01
3. Mark transfer candidates: PSMs that failed baseline but whose best pair prob ≥ threshold
4. Set `trace_prob = nonMBR_estimates` and `MBR_boosted_trace_prob = MBR_estimates`

**Without MBR:**
1. Set `trace_prob = prob_test`

The OOM branch adds a 2-pass finalization for Arrow files:
- **Pass 1:** Collect baselines + targets from all files, compute global q-values
- **Pass 2:** Merge final columns back into data files, copy MBR feature columns from sidecars

### 5.7 OOM Mechanics (feature/oom_percolator only)

#### File Layout

After `prepare_training_data!`, each (ms_file_idx, cv_fold) pair has two files:

```
{base}_fold0.arrow           ← immutable: all feature columns + pair_id + irt_bin_idx
{base}_fold0_scores.arrow    ← mutable sidecar: trace_prob, MBR_* (~32 bytes/PSM)
```

The sidecar pattern means scoring iterations only rewrite ~32 bytes/PSM instead of the full ~300-500 bytes/PSM of feature data.

#### 3-Pass Pair-Sampled Training

To train models, the OOM path samples complete target-decoy pairs from disk within a PSM budget:

1. **Pass 1 — Count:** Read only the `pair_id` column from each file. Build `Dict{pair_id => count}`.
2. **Pass 2 — Select:** Shuffle pair_ids deterministically (seed=1776). Greedily select pairs until budget exhausted.
3. **Pass 3 — Fill:** Pre-allocate DataFrame. Copy matching rows from each file using BitVector mask on `pair_id`.

This produces a representative `DataFramePSMContainer` that feeds into a normal `InMemoryScoringWorkspace`.

#### Per-File Streaming Prediction

During Phase 2, the `ArrowFileScoringWorkspace` overrides `process_fold_iterations!` for `PredictionPhase`:

```
For each iteration:
  For each Arrow file:
    Load data + scores sidecar
    Filter to current fold
    Build temporary DataFramePSMContainer
    Apply pre-trained model → predictions
    Write predictions to sidecar
  If MBR iteration:
    _compute_mbr_for_fold!() → cross-file MBR
```

This ensures only one file's data is in memory at a time.

---

## 6. File Inventory

### Files unique to `develop`

| File | Lines | Role |
|------|-------|------|
| `percolatorSortOf.jl` | 1,237 | Monolithic: pairing, training, MBR, finalization |

### Files unique to `feature/oom_percolator` (12 new files)

| File | Lines | Role |
|------|-------|------|
| `scoring_traits.jl` | 273 | 6 trait types + `ScoringPhase` + `ProbitModel` wrapper |
| `scoring_config.jl` | 109 | `ScoringConfig{M,P,T,F,I,B}` + factory functions |
| `psm_container.jl` | 265 | `AbstractPSMContainer`, `DataFramePSMContainer`, `DataFramePSMContainerView` |
| `arrow_psm_container.jl` | 162 | `ArrowFilePSMContainer`, `ArrowFileGroup` |
| `percolator_generic.jl` | 330 | `percolator_scoring!()` entry point + phase-dispatched helpers |
| `scoring_workspace.jl` | 667 | `InMemoryScoringWorkspace`, `ArrowFileScoringWorkspace`, OOM sampling + prediction |
| `pairing.jl` | 136 | `assign_pairs!()` for `RandomPairing` / `NoPairing` |
| `model_training.jl` | 101 | `train_model()` / `predict_scores()` for LightGBM / Probit |
| `training_selection.jl` | 73 | `select_training_data()` for `QValueNegativeMining` / `AllDataSelection` |
| `feature_selection.jl` | 52 | `get_features()` for `IterativeFeatureSelection` / `StaticFeatureSelection` |
| `iteration_scheme.jl` | 41 | `get_iteration_rounds()` for `FixedIterationScheme` / `SinglePassScheme` |
| `mbr_update.jl` | 139 | `initialize_mbr_columns!()`, `update_mbr_features_train_only!()` / `test_only!()` |

### Files modified on `feature/oom_percolator`

| File | develop | OOM branch | Change |
|------|---------|------------|--------|
| `percolatorSortOf.jl` | 1,237 lines | 350 lines | **-887 lines** — now a thin wrapper around `percolator_scoring!()` |
| `score_psms.jl` | OOM path disabled (`if false`) | OOM path enabled with `ArrowFilePSMContainer` | OOM is functional + `force_oom` flag added |
| `scoring_workspace.jl` | N/A | 667 lines | New file: largest single file, contains both workspace types + OOM logic |

### Files unchanged between branches

| File | Lines | Role |
|------|-------|------|
| `fdrUtilities.jl` | — | `get_qvalues!()`, `get_PEP!()` |
| `ftrUtilities.jl` | — | FTR computation |
| `lightgbm_utils.jl` | — | LightGBM API wrapper |
| `probitRegression.jl` | — | Probit regression implementation |
| `piecewiseLinearFunction.jl` | — | Piecewise linear functions |
| `spectralLinearRegression.jl` | — | Huber-loss spectral regression |
| `uniformBasisCubicSpline.jl` | — | Cubic spline implementation |
| `wittakerHendersonSmoothing.jl` | — | Smoothing |
| `libraryBSpline.jl` | — | B-spline models |

---

## 7. Data Flow Diagrams

### Develop Branch: In-Memory Only

```
                        ┌──────────────────┐
                        │  DataFrame (RAM)  │
                        │  All PSMs loaded  │
                        └────────┬─────────┘
                                 │
                    ┌────────────▼────────────┐
                    │  assign_random_target_   │
                    │  decoy_pairs!()          │
                    │  + sort by pair_id       │
                    └────────────┬─────────────┘
                                 │
              ┌──────────────────▼──────────────────┐
              │  For each CV fold:                   │
              │  ┌─────────────────────────────────┐ │
              │  │ For each iteration:             │ │
              │  │  Train LightGBM on train set    │ │
              │  │  Predict on train + test sets   │ │
              │  │  Compute q-values (train)       │ │
              │  │  Update MBR features            │ │
              │  └─────────────────────────────────┘ │
              └──────────────────┬──────────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │  Finalize: set trace_prob│
                    │  + MBR columns           │
                    └────────────┬─────────────┘
                                 │
                        ┌────────▼────────┐
                        │  DataFrame (RAM) │
                        │  with scores     │
                        └─────────────────┘
```

### OOM Branch: File-Backed Path

```
  ┌──────────────────────────────────────────┐
  │  Arrow files on disk (per ms_file × fold)│
  │  {base}_fold0.arrow, _fold1.arrow, ...   │
  └───────────────────┬──────────────────────┘
                      │
      ┌───────────────▼───────────────────┐
      │  prepare_training_data!()          │
      │  Phase 1: Read 4 cols globally     │
      │           (~19 bytes/PSM)          │
      │           Assign pairs globally    │
      │  Phase 2: Per-file sort + sidecar  │
      │           creation (~32 bytes/PSM) │
      └───────────────┬───────────────────┘
                      │
      ┌───────────────▼───────────────────┐
      │  setup_scoring_workspace()         │
      │  3-pass pair-sampled training:     │
      │  1. Count PSMs per pair_id         │
      │  2. Shuffle + select within budget │
      │  3. Fill DataFrame from matching   │
      │     rows → InMemoryScoringWorkspace│
      └───────────────┬───────────────────┘
                      │
    ┌─────────────────▼─────────────────────┐
    │  PHASE 1: Train on sampled subset     │
    │  (identical to in-memory training)    │
    │  Store models per fold                │
    └─────────────────┬─────────────────────┘
                      │
    ┌─────────────────▼─────────────────────┐
    │  PHASE 2: Stream prediction per file  │
    │  For each fold:                       │
    │    For each Arrow file:               │
    │      Load data + sidecar              │
    │      Filter to fold                   │
    │      Apply stored model → predictions │
    │      Write to sidecar                 │
    │    Cross-file MBR per fold:           │
    │      Pass 1: per-file sorted temps    │
    │      Merge → stream FDR threshold     │
    │      Pass 2: apply MBR features       │
    └─────────────────┬─────────────────────┘
                      │
    ┌─────────────────▼─────────────────────┐
    │  finalize_scoring!()                  │
    │  Pass 1: Collect baselines → global   │
    │          q-values → threshold          │
    │  Pass 2: Merge final cols into data   │
    │          files (trace_prob, MBR_*)    │
    └─────────────────┬─────────────────────┘
                      │
  ┌───────────────────▼──────────────────────┐
  │  Arrow files on disk (scored)            │
  │  Now contain trace_prob, MBR_* columns   │
  └──────────────────────────────────────────┘
```

---

## 8. Key Algorithms

### 8.1 Target-Decoy Pairing (Both Branches)

```
Input: PSMs with irt_pred, cv_fold, isotopes_captured, precursor_idx

1. bin_idx = getIrtBins(irt_pred)          # bin_size=1000
2. groups = groupby(PSMs, [bin_idx, cv_fold, isotopes_captured])
3. For each group:
     precursors = unique(precursor_idx)
     shuffle!(MersenneTwister(1844), precursors)
     pair_id = 1
     For i in 1:2:length(precursors):
       assign pair_id to all PSMs of precursors[i]
       if i+1 ≤ length: assign same pair_id to precursors[i+1]
       pair_id += 1

Output: pair_id column added to PSMs
```

### 8.2 Iterative Training (Both Branches)

```
iter_scheme = [100, 200, 200]   # boosting rounds per iteration
mbr_start_iter = 3              # MBR features used in last iteration only

Iteration 1 (100 rounds):
  Training data = ALL PSMs (no filtering)
  Features = base features only (no MBR_*)
  → Train model, predict, compute q-values

Iteration 2 (200 rounds):
  Training data = targets with q ≤ 0.01 + all decoys
                  + targets with PEP ≥ 0.90 converted to negatives
  Features = base features only
  → Train model, predict, compute q-values
  → Store nonMBR baseline (this is mbr_start_iter - 1)
  → Compute MBR features via summarize_precursors!()

Iteration 3 (200 rounds):
  Training data = same filtering as iteration 2
  Features = base features + MBR features
  → Train model, predict
  → Update MBR features again
```

### 8.3 OOM Pair Sampling (feature/oom_percolator only)

```
Budget: max_training_psms (default 50M, or all if force_oom)

Pass 1: Count
  pair_counts = Dict{UInt32, Int}()
  For each Arrow file:
    For pid in file.pair_id:
      pair_counts[pid] += 1

Pass 2: Select
  pids = shuffle(MersenneTwister(1776), keys(pair_counts))
  selected = []
  total = 0
  For pid in pids:
    if total + pair_counts[pid] ≤ budget:
      push!(selected, pid)
      total += pair_counts[pid]

Pass 3: Fill
  sampled_df = pre_allocate(total rows, schema from first file)
  For each Arrow file:
    mask = [pid ∈ selected for pid in file.pair_id]
    copy matching rows to sampled_df

Return: DataFramePSMContainer(sampled_df)
```

---

## 9. Integration with ScoringSearch

### Develop Branch

```julia
# score_psms.jl — OOM disabled
if false  # psms_count >= max_psms_in_memory
    # DISABLED
else
    best_psms = load_psms_for_lightgbm(second_pass_folder)
    # ... model selection logic ...
    models = score_precursor_isotope_traces_in_memory(best_psms, ...)
    write_scored_psms_to_files!(best_psms, file_paths)
end
```

All PSMs are loaded into a single DataFrame regardless of dataset size.

### Feature/OOM Branch

```julia
# score_psms.jl — OOM enabled
if force_oom || psms_count >= max_psms_in_memory
    # OOM path: construct ArrowFilePSMContainer from existing fold files
    container = ArrowFilePSMContainer(file_paths; max_training_psms=max_training)
    config = build_scoring_config(model_config, match_between_runs, ...)
    models = percolator_scoring!(container, config; show_progress=true, verbose=true)
    # Finalization already wrote scores back to data files — no write_scored_psms_to_files! needed
else
    # In-memory path: unchanged
    best_psms = load_psms_for_lightgbm(second_pass_folder)
    # ... same model selection logic ...
    models = score_precursor_isotope_traces_in_memory(best_psms, ...)
    write_scored_psms_to_files!(best_psms, file_paths)
end
```

The OOM path constructs an `ArrowFilePSMContainer` directly from the existing per-fold Arrow files on disk, builds a `ScoringConfig`, and calls the generic `percolator_scoring!()`. No intermediate loading step is required. The `force_oom` flag allows testing the OOM path even on small datasets.

New parameters added to `score_precursor_isotope_traces`:
- `force_oom::Bool = false` — force OOM processing regardless of PSM count
- `max_training_psms::Int64 = 50_000_000` — PSM budget for training sampling

---

## 10. Trade-offs and Design Decisions

### Benefits of the Refactoring

1. **OOM capability:** Datasets exceeding RAM can now be scored by keeping PSMs on disk and streaming through files. Training uses a representative sample; prediction operates per-file.

2. **Composability:** New scoring strategies can be added by implementing a single trait method without touching the core loop. For example, adding a neural network scorer requires only a new `PSMScoringModel` subtype + `train_model()` / `predict_scores()` methods.

3. **Testability:** Each trait implementation can be unit-tested in isolation (pairing logic, training selection, feature selection, etc.).

4. **Phase separation:** Decoupling training from prediction enables the OOM streaming prediction pattern and makes the algorithm flow clearer.

5. **Code reduction:** `percolatorSortOf.jl` went from 1,237 to 350 lines. The total code increased (2,697 lines across 13 new files), but each file has a single clear responsibility.

### Costs and Risks

1. **Complexity:** 13 new files with trait dispatch can be harder to follow than a single function for developers unfamiliar with the pattern.

2. **I/O overhead (OOM path):** The sidecar pattern writes scores back to Arrow files after each iteration/fold. For small datasets, this is slower than in-memory arrays. The sidecar optimization (~32 bytes/PSM vs ~300-500 bytes/PSM) mitigates this.

3. **Training representativeness (OOM path):** Training on a sampled subset (vs. all data) may reduce model quality for extremely large datasets where the sample is a small fraction. The pair-aware sampling ensures target-decoy pairs stay together, which is critical for the pairing-based algorithm.

4. **MBR cross-file cost (OOM path):** `_compute_mbr_for_fold!()` uses a fully streaming approach — files are processed one at a time with per-file sorted temp Arrow files and a k-way merge for global FDR computation. The pair dictionary (mapping pair_id → best 2 runs) remains in memory, but the q-value/FDR computation uses O(1) heap memory via memory-mapped Arrow files. Temp disk usage is ~5 bytes/PSM (automatically cleaned up).

5. **Known bug preserved:** `ProbitRegression` in `model_training.jl:65` expects a DataFrame but receives a `Matrix{Float64}`. This pre-existing bug exists on both branches.

### Design Decision: Two-Phase Training/Prediction

On develop, each CV fold trains and predicts in a single pass through the iterations. On the OOM branch, all folds train first (Phase 1), then all folds predict (Phase 2). This is necessary for OOM because:

- Training happens on sampled in-memory data (fast, no I/O)
- Prediction must stream through disk files (slow, I/O-bound)
- Interleaving training and streaming would thrash the disk

For in-memory data, the two-phase approach produces identical results because training and prediction are independent across folds.

### Design Decision: Sidecar Pattern

Rather than rewriting entire Arrow files (~300-500 bytes/PSM) on each scoring iteration, mutable columns are stored in a separate "sidecar" file (~32 bytes/PSM). This reduces I/O by ~10x per iteration while keeping feature data immutable and safely memory-mappable.

---

*Report generated from source code analysis of Pioneer.jl branches `develop` and `feature/oom_percolator`.*
