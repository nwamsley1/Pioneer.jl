# Peak-Shape Feature A/B — 2026-05-12

Catalog of 23 candidate features for chromatographic peak shape, edge
proximity, and DIA-NN-style profile-based scoring. Tested one batch at a
time on 2-file Olsen Exploris (`/tmp/my_tune.json`). For each: per-file
LGBM gain from the importance dump after `[after paircomp]` log.

**Baseline (before any of these features, m1frag branch):** f1 q≤.01 ~ 41,213,
f2 q≤.01 ~ 41,183, final IDs ~ 88,875 (per the morning report).

LGBM noise on 2-file is ±400-500 IDs; per-feature gain values are reliable
since they come from same-run trees.

## Batch A — Edge-distance features (your idea + extensions)

For each PSM at position `k` within its precursor's scan-sorted PSM list
of size `N`.

| # | Feature | Definition | Gain (file 1, file 2 average) | Status | Notes |
|---|---|---|---:|---|---|
| A1 | `n_scans_left_capped3` | `min(k-1, 3)` | 0 | ❌ | No splits — redundant with `relative_position` |
| A2 | `n_scans_right_capped3` | `min(N-k, 3)` | 0 | ❌ | No splits |
| A3 | `min_edge_distance` | `min(A1, A2)` | 0 | ❌ | No splits |
| A4 | `is_at_edge` | `A1==0 \|\| A2==0` (bool) | 0 | ❌ | No splits |
| A5 | `relative_position` | `(k-1) / max(N-1, 1)` | 225 | ✓ | Modest signal |
| A6 | `dist_to_relative_center` | `\|0.5 - A5\|` | 302 | ✓ | Strongest of batch |

## Batch B — Weight-chromatogram peak-shape stats

Operates on the (irt, weight) chromatogram per precursor.

| # | Feature | Definition | Gain | Status | Notes |
|---|---|---|---:|---|---|
| B1 | `weight_apex_relative_pos` | `argmax(w_chrom) / max(N-1, 1)` | 1,176 | ✓ | Medium |
| B2 | `weight_chrom_skewness` | std skewness `E[((x-μ)/σ)^3]` of weight | 2,628 | ✓ | **Top Batch-B** |
| B3 | `weight_chrom_kurtosis` | excess kurtosis `E[((x-μ)/σ)^4] - 3` | 2,281 | ✓ | Strong |
| B4 | `apex_to_edge_weight_ratio` | `max(w) / max(w_first, w_last)` | 2,243 | ✓ | Strong |
| B5 | `n_above_hm_left_of_apex` | scans left of apex with w > 0.5·max(w) | 0 | ❌ | Captured by skewness |
| B6 | `n_above_hm_right_of_apex` | symmetric | 392 | ✓ | Weak |
| B7 | `hm_asymmetry` | `\|B5 - B6\| / max(B5+B6, 1)` | 0 | ❌ | Redundant with skewness |

## Batch C — DIA-NN-style Gaussian fit

Fit `log(w) = a + b·irt + c·irt²` per precursor (linear regression).

| # | Feature | Definition | Gain | Status | Notes |
|---|---|---|---:|---|---|
| C1 | `weight_chrom_gaussian_r2` | R² of log-quadratic fit on weight chrom | 1,721 | ✓ | Medium |
| C2 | `weight_chrom_gaussian_sigma` | fitted σ (peak width); requires c<0 | 1,445 | ✓ | Medium |
| C3 | `weight_chrom_gaussian_apex_irt` | fitted apex iRT (peak center) | 1,530 | ✓ | Medium |

## Batch D — DIA-NN direct ports (`pShape`, `pBestCorrDelta`, `pTotCorrSum`)

Direct ports from `diann.cpp` lines 7700-7717.

| # | Feature | Definition | Gain | Status | Notes |
|---|---|---|---:|---|---|
| D1 | `shape_m2` | normalized weight at scan-position k-2 (vs current PSM scan) | 578 | ✓ | Weak |
| D2 | `shape_m1` | normalized weight at k-1 | 455 | ✓ | Weak |
| D3 | `shape_0` | normalized weight at k | 795 | ✓ | Weak-medium |
| D4 | `shape_p1` | normalized weight at k+1 | 597 | ✓ | Weak |
| D5 | `shape_p2` | normalized weight at k+2 | 583 | ✓ | Weak |
| D6 | `gof_minus_max_gof_precursor` | `gof - max(gof over precursor's scans)` (≤0) | 915 | ✓ | Medium |
| D7 | `gof_fraction_of_total_precursor` | `log(gof / (sum(gof) + 1))` | 0 | ❌ | Captured by D6 |

## Per-batch run results

(Filled in as each batch finishes.)

| Batch | Final IDs | f1 q≤.01 | f2 q≤.01 | PGs q≤.01 | Δ vs baseline |
|---|---:|---:|---:|---:|---:|
| Baseline (m1frag, no peak-shape features) | 88,875 | 41,213 | 41,183 | 12,386 | — |
| + Batch A | 88,280 | 41,304 | 40,768 | 12,334 | −595 (noise) |
| + Batch B | 88,612 | 41,140 | 40,913 | 12,311 | −263 (noise) |
| + Batch C | 88,217 | 41,101 | 40,825 | 12,260 | −658 (noise) |
| + Batch D | 88,798 | 40,951 | 40,800 | **12,506** | −77 IDs; **+120 PGs** |

## Synthesis (after all 4 batches)

23 features tested across 4 batches. Cumulative on 2-file Olsen at end of
Batch D (all A+B+C+D active): **88,798 precursors / 12,506 PGs** vs baseline
**88,875 / 12,386**. Final IDs essentially flat (LGBM noise ±400) but **PG
count went up +120** — first PG-positive result on Olsen for any feature
batch we've added.

### Strong keepers (gain > 1500)

| Feature | Batch | Gain | Why |
|---|---|---:|---|
| `weight_chrom_skewness` | B | 2,628 | Best single new feature |
| `apex_to_edge_weight_ratio` | B | 2,243 | Clean apex isolation signal |
| `weight_chrom_kurtosis` | B | 2,281 | Spike-vs-broad signal |
| `weight_chrom_gaussian_r2` | C | 1,721 | DIA-NN-style fit quality |
| `weight_chrom_gaussian_apex_irt` | C | 1,530 | |

### Medium keepers (400-1500)

`weight_chrom_gaussian_sigma` (C, 1,445), `weight_apex_relative_pos` (B, 1,176),
`gof_minus_max_gof_precursor` (D, 915), `shape_0` (D, 795), `shape_p1` (D, 597),
`shape_p2` (D, 583), `shape_m2` (D, 578), `shape_m1` (D, 455),
`dist_to_relative_center` (A, 302), `relative_position` (A, 225).

### To drop (no LGBM splits — 6 features)

- A1-A4 `n_scans_left_capped3`, `n_scans_right_capped3`, `min_edge_distance`, `is_at_edge` — redundant with `relative_position`
- B5 `n_above_hm_left_of_apex`, B7 `hm_asymmetry` — captured by `weight_chrom_skewness`
- D7 `gof_fraction_of_total_precursor` — captured by `gof_minus_max_gof_precursor`

### Net pattern

The peak-shape features (B + C + D) are doing real work even if they don't
move 2-file q≤.01 IDs much — they're already overlapping with `worst_max_residual_11scan`
(the dominant feature from the M+1 work) which itself captures a lot of
the same chromatographic-quality signal. The PG-count improvement (+120) is
the most interesting result and may be a signal that these features improve
protein-level scoring more than PSM-level.

Expect bigger lifts on 23-file Olsen and MTAC where LGBM noise is smaller
relative to the lift.

## Batch E — alphaDIA inspiration

From `/Users/nathanwamsley/Downloads/alphadia-main/alphadia/search/scoring/features/{fragment_features,profile_features}.py`. Candidates for future testing.

| # | Feature | alphaDIA `feature_array[]` | Definition | Why it's novel for us |
|---|---|---|---|---|
| E1 | `pred_obs_area_corr` | 18 | Pearson(observed_fragment_area, library_predicted_intensity) | We have spectral_contrast (cosine); Pearson is a different statistic |
| E2 | `pred_obs_height_corr` | 19 | Pearson(observed_fragment_height, library_predicted_intensity) | Mirror but uses peak height instead of area |
| E3 | `frag_area_explained` | 22 | Σ(predicted_intensity_norm[i] for observed_area[i]>0) | Fraction of predicted intensity mass that was observed |
| E4 | `log_b_intensity` | 25 | `log(Σ b_ion_intensities + 1)` | We have b_count but not log-summed intensity |
| E5 | `log_y_intensity` | 26 | `log(Σ y_ion_intensities + 1)` | Same for y |
| E6 | `log_b_y_intensity_ratio` | 27 | E4 − E5 | Real peptides have characteristic b/y ratios |
| E7 | `top3_ms2_mass_error_mean` | 41 | mean ppm error of top-3 highest-intensity fragments | We have `err_norm` but it's all fragments; top-3 mass errors are cleaner |
| E8 | `mean_ms2_mass_error` | 42 | mean ppm error across all fragments | |
| E9 | `n_b_y_overlapping` | 43 | # b/y ions at conflicting sequence positions | Co-occurrence of b and y at conflicting positions = contamination |
| E10 | `mean_fragment_to_median_corr` | 31 | mean Pearson(each fragment, median profile across fragments) | Different reference: median rather than consensus-best |
| E11 | `top3_b_ion_correlation` | 34 | mean correlation of top-3 b-ions to median profile | Per-ion-type breakdown |
| E12 | `top3_y_ion_correlation` | 36 | mean correlation of top-3 y-ions to median profile | |
| E13 | `weighted_fragment_fwhm` | 38 | Σ (per-fragment FWHM × predicted_intensity) | Per-fragment FWHM weighted by intensity (we have per-precursor only) |
| E14 | `delta_frame_peak_center` | 40 | (median fragment apex position − search-window center) × observation_importance | RT shift signal relative to search-window center |
| E15 | `template_frame_correlation` | 33 | Weighted correlation of each fragment to a TEMPLATE (predicted) chromatogram shape | Predicted-vs-observed shape similarity |

### Notes on alphaDIA's approach

- They feed the LightGBM ~46 features (vs our 70+); many fewer but more targeted.
- They use **median profile** as the chromatographic reference, not consensus-best. We could test this as a variant of `frag_corr_best_*`.
- They have explicit **b vs y ion separation** in profile correlations. We have b_count / y_count / longest_y but no per-ion-type chromatogram correlations.
- The **b/y overlapping count** (E9) is interesting: when a b-ion at position 5 and a y-ion at position 7 of a 10-residue peptide both score, the y-ion's position 4 from C-term contradicts b-ion presence at N-term residue 5+1=6. Real peptides shouldn't have many such contradictions; chimeric hits often do.
- Their `delta_frame_peak` (E14) measures whether the fragments' apex is near the search-window center. We've handled this via `irt_dist_to_weight_apex` but theirs is per-fragment-aggregated.
