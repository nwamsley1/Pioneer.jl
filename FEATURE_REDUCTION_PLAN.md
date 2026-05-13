# Feature Reduction Plan — Beyond the 66/71 Smart Composite

Date: 2026-05-13
Branch: `feature/mbr-batch-f` (current commit: `143d6b87`)

## Current state

After the smart-composite bake-in:

- **MainSearch PRESCORE_FEATURES**: 66 features
- **ScoringSearch ADVANCED_FEATURE_SET**: 71 features (= PRESCORE + 5 MBR-only)
- The two LGBMs now train on the same feature columns apart from 5 MBR
  signals that depend on cross-run pooling.

Verified gains on 2-file Olsen / MTAC vs full-feature baseline (with MBR):

| Dataset | Per-file Δ | Final IDs Δ | PGs Δ |
|---------|-----------:|------------:|------:|
| Olsen   |       +364 |        +240 |  +112 |
| MTAC    |       −797 |        +672 |  +359 |

## Goal

Push further reductions to:
1. Save compute (struct field captures, per-precursor correlations, MS1
   bookkeeping).
2. Simplify the feature engineering surface area.
3. Eventually enable LGBM column subsampling for additional speed.

We stop when removing a feature loses more IDs than we tolerate
(roughly: > 1% drop sustained on either dataset's final IDs).

## Removal candidates by tier

### Tier 1 — Expensive AND already showed positive in earlier ablation

| Feature | Earlier single-drop Δ (Olsen per-file) | Compute cost saved |
|---|---:|---|
| `frag_corr_mean_pairwise`        | +557 | Spearman — 15 rank-sorts per precursor. Most expensive single feature. |
| `frag_corr_m0_m1_top1`           | (only `m0_m1_corr` member kept) | Removing this feature lets us delete the entire M+1 fragment-intensity capture machinery: 6 fields on `MainUnscoredPSM` / `MainSearchScoredPSM`, the per-scan capture in `apply_main_scoring!`, and the M+1 chromatogram correlation pass in `_add_fragment_chromatogram_features!`. |
| `poisson`                        | +429 | Statistical artifact from `total_ions`. Likely redundant. |
| `fitted_hellinger`               | +357 | Similar information to `fitted_manhattan_distance`. |
| `weight_ratio_at_scan` + `weight_rank_at_scan` | +548 / +505 | Both require cross-precursor scan lookups. |
| `irt_dist_to_weight_apex`        | +366 | Requires per-precursor apex lookup. |

### Tier 2 — Likely redundant pairs

| Feature | Redundant with | Reason |
|---|---|---|
| `rt_fwhm` | `irt_fwhm` | Same chromatogram, different time units; earlier `fwhm_keep_irt_fwhm` was +282 (winner). |
| `num_scans` | `n_scans` | Both per-precursor scan counts — likely identical. |
| `irt_pred` (vs `irt_diff`) | — | `irt_diff = irt_pred − irt_obs`; the LGBM has both. Drop `irt_pred`. |
| `best_rank_iso` (vs `topn_iso`) | — | Both encode "isotope rank coverage." |
| `total_ions_iso` (vs `topn_iso`) | — | Similar isotope-count signal. |

### Tier 3 — Per-family further reduction

| Family | Current size | Candidate reduction |
|---|---:|---|
| `win5` (5-scan window features) | 6 (`best_gof_5scan`, `best_manhattan_5scan`, `best_max_residual_5scan`, three `irt_dist_best_*_5scan`) | Keep 1–2 winners. Earlier `keep_best_manhattan_5scan` was +549. |
| `frag1..6_int` (per-rank M0 fragment intensities) | 6 | Try keeping `frag1_int` + a sum/mean only — would let us delete 5 struct fields if it holds. |
| MS1 envelope: `ms1_m1_to_m0_ratio`, `ms1_m1_to_m0_pred`, `ms1_envelope_dev_log2` | 3 | Test keep only `ms1_envelope_dev_log2` (which is the log-difference of obs vs pred). |

## Execution plan

### Phase 1 — Tier 1 single-drops (one feature at a time)

Each ablation runs the full pipeline on 2-file Olsen Exploris + 2-file
MTAC 3P (MBR on, FTR-PEP α=0.005), records per-file MainSearch q≤.01
sum, final IDs, PGs. Comparison anchor: the baked 66/71 baseline.

1. Drop `frag_corr_mean_pairwise`
2. Drop `frag_corr_m0_m1_top1`
3. Drop `poisson`
4. Drop `fitted_hellinger`
5. Drop `weight_ratio_at_scan` + `weight_rank_at_scan`
6. Drop `irt_dist_to_weight_apex`

### Phase 2 — Tier 2 redundant pairs

7. Drop `rt_fwhm`
8. Drop `irt_pred` (keep `irt_diff`)
9. Drop `num_scans` (keep `n_scans`)
10. Drop `best_rank_iso` (keep `topn_iso`)
11. Drop `total_ions_iso` (keep `topn_iso`)

### Phase 3 — Tier 3 family reductions

12. `win5_keep_best_manhattan_5scan` (drop 5 of 6)
13. `frag*_int` reduction (test keep top-3 only first)
14. MS1 envelope reduction

### Phase 4 — Composite minimal

Compose all the Tier 1/2/3 winners (those that don't cost > 1% IDs)
into a single reduced blacklist, verify on both datasets. Bake into
code if the composite holds.

### Phase 5 (downstream) — Compute cleanup

If `frag_corr_m0_m1_top1` survives removal:
- Delete `frag1..6_int_m1` fields from `MainUnscoredPSM` and `MainSearchScoredPSM`.
- Delete the M+1 branch in `apply_main_scoring!`.
- Delete the M+1 correlation pass in `_add_fragment_chromatogram_features!`.

If `frag_corr_mean_pairwise` survives removal:
- Delete the Spearman branch in `_add_fragment_chromatogram_features!`.

Then re-test that nothing else relies on the dropped computations.

### Phase 6 (LGBM speed) — column subsampling

Once features are minimal, try LightGBM `feature_fraction_bynode` /
`feature_fraction` < 1 to randomly subsample features per tree node.
Reduces training time per tree.

## Tolerance

Treat as "drop OK" when ALL of:
- Olsen final IDs Δ ≥ −300 (~ noise floor for 2-file)
- MTAC final IDs Δ ≥ −600 (~ 0.3%)
- Olsen / MTAC per-file Δ ≥ −500

Otherwise, keep the feature.
