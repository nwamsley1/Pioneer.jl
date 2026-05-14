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

## Tier-1 Results (2026-05-13)

Tested 6 single-feature drops on 2-file Olsen Exploris + 2-file MTAC 3P
(MBR on, FTR-PEP α=0.005). Baselines: Olsen 90,593 IDs / 12,445 PGs;
MTAC 179,461 IDs / 20,354 PGs.

| Feature | Olsen Δ IDs | Olsen Δ PGs | MTAC Δ IDs | MTAC Δ PGs | Verdict |
|---|---:|---:|---:|---:|---|
| `frag_corr_mean_pairwise` | **+471** | **+140** | **+649** | **+8** | **DROP** |
| `frag_corr_m0_m1_top1` | −133 | +6 | **−771** | **−157** | KEEP (expensive — revisit) |
| `poisson` | −120 | +12 | **−800** | **−134** | KEEP |
| `fitted_hellinger` | +221 | +68 | **−1,437** | **−197** | KEEP |
| `weight_ratio_at_scan` + `weight_rank_at_scan` | +241 | +50 | +30 | −89 | borderline KEEP |
| `irt_dist_to_weight_apex` | +344 | +21 | −294 | −172 | KEEP |

**Only `frag_corr_mean_pairwise` is cross-dataset safe.** Dropped from
both `PRESCORE_FEATURES` and `ADVANCED_FEATURE_SET`. This eliminates the
Spearman 15-rank-sort-per-precursor computation — meaningful runtime win.

**Note on `frag_corr_m0_m1_top1`**: passed the cross-dataset "KEEP" test
(MTAC loses 771 IDs / 157 PGs without it). BUT this single feature is
the only consumer of the entire M+1 fragment-intensity capture machinery
— 6 fields on `MainUnscoredPSM` + `MainSearchScoredPSM`, the `iso_idx==1`
branch in `apply_main_scoring!`, and the M+1 chromatogram reshape +
correlation pass in `_add_fragment_chromatogram_features!`. The
ID/compute trade-off is real: a 771-ID MTAC loss may be acceptable in
exchange for substantially less per-scan capture cost and simpler code.
Revisit this decision once Tier-2 and Tier-3 reductions are done — at
that point we'll know how many IDs we have to play with overall.

Cross-dataset pattern: features that look droppable on Olsen alone
(positive single-drop Δ) tend to be load-bearing on MTAC. The minimal
cross-validated reduction is smaller than the per-dataset minimum.

New sizes:
- PRESCORE_FEATURES: 66 → 65
- ADVANCED_FEATURE_SET: 71 → 70

Next: Tier 2 (redundant pairs) and Tier 3 (per-family further reduction)
on the new 65/70 baseline.

## Tier-2 Results (2026-05-13)

Tested 5 single-feature drops on 2-file Olsen + MTAC using the 65/70
baseline (post-Tier-1). Baselines: Olsen 91,064 IDs / 12,585 PGs;
MTAC 180,110 IDs / 20,362 PGs.

| Drop              | Olsen Δ IDs | Olsen Δ PGs | Olsen pf Δ | MTAC Δ IDs | MTAC Δ PGs | MTAC pf Δ | Verdict |
|-------------------|------------:|------------:|-----------:|-----------:|-----------:|----------:|---------|
| `rt_fwhm`         |        −189 |         −34 |          0 |       −331 |        −22 |         0 | KEEP    |
| `num_scans`       |        −196 |         −76 |          0 |       −176 |        −29 |         0 | KEEP    |
| `irt_pred`        |        −404 |        −165 |       +340 |     −1,143 |       −146 |       −75 | KEEP    |
| `best_rank_iso`   |        −399 |        −139 |       +265 |       −375 |        +50 |        +4 | KEEP    |
| `total_ions_iso`  |        −378 |         −12 |       +190 |       −723 |       −110 |      −551 | KEEP    |

**0 of 5 candidates safe to drop end-to-end.** All cost ≥ 176 final IDs
on at least one dataset.

Key pattern:
- `rt_fwhm` and `num_scans` have per-file Δ = 0 on BOTH datasets — truly
  MainSearch-redundant. The end-pipeline penalty comes from ScoringSearch.
  Could be candidates for "drop from PRESCORE only, keep in ADVANCED" —
  saves per-file LGBM training time (small) without losing IDs.
- `irt_pred`, `best_rank_iso`, `total_ions_iso` slightly *help* per-file
  LGBM when dropped (Olsen pf +190 to +340) — they may add confusion at
  per-file level — but ScoringSearch needs them for cross-file
  differentiation.

No code changes from Tier 2. Sizes remain 65 PRESCORE / 70 ADVANCED.
