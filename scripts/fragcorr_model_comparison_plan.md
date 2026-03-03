# Plan: FRAGCORR Model Comparison — Do eigengap features improve ML scoring?

## Context

From the `fragcorr_vs_features.jl` results, `eigengap_raw` and `lambda1_frac_raw` are strong unsupervised FRAGCORR scores (~107K targets at 1% FDR). Since they capture multi-scan fragment coherence (fundamentally different from single-scan spectral match features like scribe), they may be complementary. This script tests whether adding FRAGCORR features to LightGBM and probit models improves target-decoy discrimination.

## Data

Same 6 file pairs in `/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_test/temp_data/`:
- **PSM files** (`first_pass_psms/*.arrow`): 1 row/precursor, columns include `ms_file_idx`, `precursor_idx`, `scribe`, `spectral_contrast`, `city_block`, `entropy_score`, `y_count`, `err_norm`, `target`, `score`, `prob`
- **CORR files** (`first_pass_corr/*.arrow`): many rows/precursor, `precursor_idx, intensity_1..5, rt, target`

## File: `scripts/fragcorr_model_comparison.jl` (~250 lines)

### Step 1: Load & Join Data

Include `fragcorr_unsupervised_scores.jl` for FRAGCORR scoring + q-value utilities. Load PSM+CORR pairs, compute FRAGCORR scores, inner join on `(precursor_idx, target)`, concatenate all 6 files. Same pattern as `fragcorr_vs_features.jl` `load_and_join`.

### Step 2: Define Feature Sets

```julia
BASELINE_FEATURES = [:scribe, :spectral_contrast, :city_block, :entropy_score, :y_count, :err_norm]
FRAGCORR_FEATURES = [:eigengap_raw, :lambda1_frac_raw, :ev1_raw, :median_corr_raw, :best_frag_corr_mean]
AUGMENTED_FEATURES = vcat(BASELINE_FEATURES, FRAGCORR_FEATURES)
```

Note: drop `poisson` (1 target at all thresholds — non-discriminative in this data).

### Step 3: Cross-Validated Model Training

**CV fold assignment**: Assign random 3-fold CV based on `precursor_idx` (hash-based for reproducibility — `precursor_idx % 3 + 1`). This ensures the same precursor always falls in the same fold.

**4 experimental conditions** (2 models x 2 feature sets):
1. LightGBM + baseline features
2. LightGBM + augmented features
3. Probit + baseline features
4. Probit + augmented features

**LightGBM training** — use `LightGBM.jl` directly:
- `LGBMClassification` with hyperparams similar to Pioneer's SimpleLightGBM config (max_depth=4, num_leaves=15, learning_rate=0.1, 200 iterations)
- Train on folds != k, predict held-out fold k
- Collect all held-out predictions

**Probit training** — standalone implementation (~40 lines):
- Port the core IRLS loop from `src/utils/ML/probitRegression.jl` without `@turbo` dependency (use plain loops + `SpecialFunctions.erf`)
- Same CV structure as LightGBM
- Key functions: `probit_train(X, y; max_iter=30)` -> beta, `probit_predict(X, beta)` -> probabilities

### Step 4: Evaluate & Compare

For each of the 4 conditions:
- Compute q-values from held-out predictions using `compute_qvalues`
- Count targets at thresholds [0.1%, 0.5%, 1%, 2%, 5%, 10%]

Also include reference lines:
- Existing `score`/`prob` from the arrow files (already-trained pipeline model)
- Best individual FRAGCORR score (`eigengap_raw`) as unsupervised baseline

### Step 5: Output

**Table**: Targets at q-value thresholds for all conditions
**Plot**: Q-value curves (6 lines: 4 conditions + `prob` reference + eigengap unsupervised)
**LightGBM feature importances**: Print gain-based importances for augmented model to see how much weight the FRAGCORR features get

## Key Reuse

| Function | Source | Purpose |
|---|---|---|
| `load_fragcorr_data()` | `fragcorr_unsupervised_scores.jl` | Load CORR arrow, group by precursor |
| `score_all_precursors()` | `fragcorr_unsupervised_scores.jl` | Compute 51 FRAGCORR scores |
| `compute_qvalues()` | `fragcorr_unsupervised_scores.jl` | Target-decoy q-value calc |
| `targets_at_qvalue()` | `fragcorr_unsupervised_scores.jl` | Count targets at threshold |
| LightGBM API pattern | `src/utils/ML/lightgbm_utils.jl` | `LGBMClassification`, `fit!`, `predict` |
| Probit IRLS algorithm | `src/utils/ML/probitRegression.jl` | Core algorithm (simplified port) |

## Dependencies

```julia
using Arrow, DataFrames, Tables           # Already in base script
using Statistics, LinearAlgebra           # Already in base script
using StatsBase                           # Already in base script
using Plots                               # Already in base script
using LightGBM                            # Pioneer dependency
using SpecialFunctions                    # For probit erf()
```

## Verification

1. `include("scripts/fragcorr_model_comparison.jl")` — runs full pipeline
2. Check that augmented LightGBM >= baseline LightGBM at all thresholds (if FRAGCORR helps)
3. Check feature importances — FRAGCORR features should have nonzero gain
4. Probit comparison provides a linear-model sanity check
5. Compare all models against existing `prob` reference from the pipeline
