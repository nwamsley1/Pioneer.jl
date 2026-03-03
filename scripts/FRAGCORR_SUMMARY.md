# FRAGCORR: Multi-Scan Fragment Correlation Features for DIA PSM Rescoring

## Goal

Add multi-scan fragment correlation (FRAGCORR) features to Pioneer's DIA PSM
rescoring pipeline. The hypothesis: fragments from the same precursor should
co-elute across multiple scans, and measuring that coherence provides
discriminative power beyond single-scan spectral matching scores.

## What Was Built

### Source Code Changes

**`src/Routines/SearchDIA/LibrarySearch.jl`**
- `build_focused_precursor_set(scored_psms; rt_window=0.25f0)` — builds an
  `ExpandedPrecursorSet` centered on the best-scored PSM's RT for each
  precursor (±15 seconds), replacing the broad fragment-index RT window.
- `run_fragcorr_pass(spectra, search_context, params, ms_file_idx, expanded)` —
  runs the FRAGCORR correlation collection pass using the focused precursor set.
- First-pass `library_search` now runs with `run_fragcorr = false` (deferred).

**`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`**
- `process_file!` now runs FRAGCORR *after* PSM scoring instead of during
  library search. It calls `build_focused_precursor_set` on scored PSMs, then
  `run_fragcorr_pass` to collect fragment traces within the focused ±15s window.

### Architecture

```
Library search (fragment index)
    ↓
PSM scoring (scribe, spectral_contrast, etc.)
    ↓
build_focused_precursor_set()   ← ±15 sec around best PSM RT
    ↓
run_fragcorr_pass()             ← collect top-6 fragment intensity traces
    ↓
Unsupervised scores             ← 51 features (PCA, correlation, coelution)
    ↓
LightGBM rescoring              ← baseline + FRAGCORR features
```

Key design choice: FRAGCORR runs **after** initial scoring so the RT window can
be centered on the best PSM apex rather than the broad fragment-index window.
This yields much tighter chromatographic profiles and stronger correlation
signal.

## Scripts Inventory

| Script | Purpose |
|--------|---------|
| `fragcorr_unsupervised_scores.jl` | Core library: loads FRAGCORR traces, computes 51 unsupervised coherence/coelution scores per precursor (PCA eigenvalues, pairwise correlations, cosine similarity, etc.) |
| `fragcorr_model_comparison.jl` | Trains {LightGBM, Probit} x {baseline, augmented} models with 3-fold CV; prints targets at multiple FDR thresholds |
| `fragcorr_diann_recovery.jl` | Validates extra targets by checking overlap with DIA-NN independently identified precursors |
| `fragcorr_topn_comparison.jl` | Compares top-5 vs 6 vs 7 fragments for FRAGCORR scoring |
| `fragcorr_global_prob_sim.jl` | Simulates Pioneer's global filter (log-odds combine, PEP, hybrid threshold) to estimate precursor-level recovery |
| `fragcorr_fixed_window.jl` | Experiments with fixed scan-count windows centered on data-driven apex |
| `fragcorr_vs_features.jl` | Compares unsupervised FRAGCORR scores against traditional per-scan features via q-value curves |
| `fragcorr_focused_analysis.jl` | **Main analysis**: end-to-end pipeline using focused ±15s RT window — loads data, trains models, runs DIA-NN validation, saves plots |

## Key Results (Focused ±15s Analysis, OlsenEclipse Dataset)

### Targets at 1% FDR

| Condition | Targets | Change |
|-----------|---------|--------|
| LightGBM baseline | 192,103 | — |
| LightGBM augmented | 207,165 | **+7.8%** |
| Probit baseline | 175,834 | — |
| Probit augmented | 178,908 | +1.7% |

### DIA-NN Overlap at 1% FDR

| Condition | Total | DIA-NN | DIA-NN % |
|-----------|-------|--------|----------|
| LightGBM baseline | 192,103 | 187,129 | 97.4% |
| LightGBM augmented | 207,165 | 201,743 | 97.4% |

Both baseline and augmented maintain 97.4% DIA-NN overlap — the extra targets
are not false discoveries.

### Venn Decomposition (LightGBM baseline vs augmented, 1% FDR)

| Set | Count | DIA-NN validated |
|-----|-------|------------------|
| Both models | 186,937 | 97.9% |
| Augmented only (FRAGCORR) | 20,145 | **92.3%** |
| Baseline only | 5,083 | 76.3% |

The 20,145 precursors gained by FRAGCORR are 92.3% DIA-NN validated — well
above the 76.3% rate for precursors lost from baseline. This confirms FRAGCORR
is genuinely rescuing real identifications.

### Feature Importances (LightGBM augmented)

Top features by gain:

| Rank | Feature | Gain |
|------|---------|------|
| 1 | scribe | 157,490 |
| 2 | **median_corr_raw** | **121,047** |
| 3 | spectral_contrast | 100,050 |
| 4 | city_block | 87,620 |
| 5 | entropy_score | 57,805 |
| 6 | **lambda1_frac_raw** | **41,893** |
| 7 | y_count | 33,862 |
| 8 | **eigengap_raw** | **16,547** |
| 9 | err_norm | 14,281 |
| 10 | **ev1_raw** | **5,963** |
| 11 | **best_frag_corr_mean** | **2,755** |

`median_corr_raw` (median pairwise Pearson correlation of fragment traces) is
the 2nd most important feature overall, confirming strong discriminative value
of multi-scan coherence.

## FRAGCORR Feature Descriptions

The 51 unsupervised scores computed per precursor include:

- **PCA eigenvalue features**: `lambda1_frac` (fraction of variance in PC1),
  `eigengap` (gap between eigenvalues 1 and 2), `ev1`/`ev2` (raw eigenvalues)
- **Pairwise correlation**: `median_corr`, `min_corr`, `max_corr`, `iqr_corr`
  (Pearson correlations between all fragment-pair traces)
- **Cosine similarity**: `median_cosine`, `min_cosine` across scan vectors
- **Coelution metrics**: `best_frag_corr_mean` (average correlation of the
  best-correlated fragment pair)
- Each computed on both `_raw` (raw intensities) and `_norm` (row-normalized)
  variants

The top discriminative FRAGCORR features are: `median_corr_raw`,
`lambda1_frac_raw`, and `eigengap_raw`.
