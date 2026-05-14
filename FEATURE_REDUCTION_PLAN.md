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

## Tier-2 Drop-All-5 (2026-05-13)

Tested all 5 Tier-2 features dropped simultaneously
(`rt_fwhm,num_scans,irt_pred,best_rank_iso,total_ions_iso`).

| Dataset | Baseline (65/70) | Drop-all-5 | Δ IDs | Δ PGs | pf Δ |
|---|---:|---:|---:|---:|---:|
| Olsen | 91,064 / 12,585 | **91,202 / 12,577** | **+138** | −8 | **+463** |
| MTAC  | 180,110 / 20,362 | 179,705 / 20,372 | −405 | **+10** | −1,062 |

**Surprise: the combined drop is better than any single drop**, and well
inside tolerance on both datasets. The 5 features were interacting
destructively — removing any one made things worse (each cost 176–1,143
final IDs individually) but removing all 5 together is net-neutral on
Olsen IDs (+138 IDs) and roughly neutral on MTAC (−0.2%).

**Verdict: DROP ALL 5.** New sizes: 60 PRESCORE / 65 ADVANCED.

Apply this drop next: edit `PRESCORE_FEATURES` and `ADVANCED_FEATURE_SET`
to remove `rt_fwhm, num_scans, irt_pred, best_rank_iso, total_ions_iso`.

## Pre-filter threshold diagnostic (2026-05-13)

Goal: identify cheap row-level filters that remove PSMs essentially
nothing will score above PEP 0.9, to shrink the data flowing into LGBM.

Setup: Olsen 2-file run with `PIONEER_MAIN_PEP_FILTER_THR=1.0` (PEP
filter effectively off) and `delete_temp=false`. Counted how many rows
in the main-search and second-pass PSM tables each candidate rule would
drop.

**main_search_psms (per-fold, 4 files):**

| Rule | File 1 fold0 | File 1 fold1 | File 2 fold0 | File 2 fold1 |
|---|---:|---:|---:|---:|
| `topn==0 && topn_iso==0` | 2.0% | 1.9% | 1.4% | 1.3% |
| `total_ions <= 2`        | **8.5%** | **8.3%** | **20.4%** | **21.2%** |
| `b_count + y_count == 0` | 0.0% | 0.0% | 0.1% | 0.1% |
| `best_rank > 6`          | 0.3% | 0.2% | 0.3% | 0.3% |
| `weight < 1e-4`          | 0.0% | 0.0% | 0.0% | 0.0% |
| **UNION of all 5**       | 9.6% | 9.4% | 21.0% | 21.9% |

**second_pass_psms (2 files):**

| Rule | File 1 (842,829) | File 2 (983,142) |
|---|---:|---:|
| `topn==0 && topn_iso==0` | 1.9% | 1.3% |
| `total_ions <= 2`        | 8.4% | **20.8%** |
| `b_count + y_count == 0` | 0.0% | 0.1% |
| `best_rank > 6`          | 0.3% | 0.3% |
| `weight < 1e-4`          | 0.0% | 0.0% |
| **UNION of all 5**       | 9.5% | 21.4% |
| Union AND `prob >= 0.5`  | 5.0% | 11.6% |

**Findings:**
- `total_ions <= 2` is the only rule that meaningfully shrinks data
  (8–21% of rows). Everything else is < 0.5%.
- `weight < 1e-4` removes nothing — the per-precursor deconv weight is
  almost always ≥ 1e-4 by the time these tables are written.
- **Caveat: 5–11.6% of the union-dropped rows have main-search
  `prob ≥ 0.5`.** That's not the same as "would be kept at q ≤ 0.01"
  but it's a clear signal that `total_ions <= 2` is NOT safe — many
  short-peptide PSMs (especially in File 2) score reasonably well and
  would be erased by this rule.
- The cheap, safe filters (`topn==0 && topn_iso==0`, `b_count+y_count==0`,
  `weight < 1e-4`, `best_rank > 6`) collectively remove < 2.5%. Not
  worth a new code path on their own.

**Recommendation:** Do not add a pre-LGBM `total_ions <= 2` filter
without an ablation showing IDs are preserved. The cheap rules are too
small to bother with. The PEP-0.9 filter is already doing the bulk of
the work.

## BitVec min_excess_rate 0.03 → 0.02 (2026-05-13)

Tested lowering the bitvec-LUT excess-rate threshold from 0.03 to 0.02
on the 60/65 feature set. Made `BITVEC_MIN_EXCESS_RATE` env-var
configurable via `PIONEER_BITVEC_MIN_EXCESS_RATE` (default unchanged at 0.03).

| Dataset | bitvec=0.03 | bitvec=0.02 | Δ IDs | Δ PGs |
|---|---:|---:|---:|---:|
| Olsen | 91,202 / 12,577 | **92,397 / 12,595** | **+1,195** (+1.31%) | +18 |
| MTAC  | 179,954 / 20,395 | **180,510 / 20,406** | **+556** (+0.31%) | +11 |

Both datasets positive — Olsen substantially so. **Not changing default** —
recorded for future consideration. Worth re-confirming on YeastMBR /
HelaOnly before baking in.

## Tier-3 baseline (60/65 set, bitvec=0.03)

| Dataset | IDs | PGs | pf q≤.01 sum |
|---|---:|---:|---:|
| Olsen | 91,202 | 12,577 | 83,198 |
| MTAC  | 179,954 | 20,395 | (from log) |


## Tier-3 Results (2026-05-13)

Tested 9 single-feature / family drops on 60/65 baseline.
Olsen 2-file baseline: 91,202 IDs / 12,577 PGs / pf=83,198.
MTAC 2-file baseline: 179,896 IDs / 20,481 PGs / pf=157,678.

| Drop | Olsen ΔIDs | Olsen ΔPGs | Olsen Δpf | MTAC ΔIDs | MTAC ΔPGs | MTAC Δpf | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| `drop_win11` (worst_*_11scan ×2)   | +276 | −142 | **−790** | **+619** | **+188** | **+2373** | borderline (Olsen pf bad) |
| `drop_win5_keep_manhattan` (×5)    | +287 | +39  | −184 | +190 | −107 | −385 | borderline |
| `drop_frag_int_keep_frag1` (×5)    | **−223** | −60 | −3  | +177 | −118 | +358 | KEEP (Olsen IDs lose) |
| `drop_ms1_env_keep_dev` (m0_ratio, m0_pred) | +187 | −23 | **−779** | +95 | +8 | +345 | borderline (Olsen pf bad) |
| `drop_weight_at_scan_pair`         | +291 | −50 | +103 | +14 | −47 | −431 | borderline |
| `drop_best_max_residual_3scan`     | +271 | +44 | −19 | **−1468** | **−158** | +425 | **KEEP** (MTAC disaster) |
| `drop_ms1_m1_intensity`            | **+448** | −41 | −8 | +127 | −135 | −78 | **DROP** (cleanest) |
| `drop_top3_ms2_mass_error_mean`    | +156 | −113 | −253 | +402 | −42 | +392 | borderline |
| `drop_log_by_ratio_m0`             | +27 | −136 | −101 | +395 | −41 | +507 | borderline (PGs) |

**Clean wins (positive IDs both datasets, pf within tolerance):**
- `drop_ms1_m1_intensity` is the only single-drop with positive IDs on both
  datasets AND no per-file MainSearch impact.

**MBR-rescue caveat:** `drop_win11` and `drop_ms1_env_keep_dev` each lose
~780 per-file Olsen IDs but recover positive at final ID level via MBR.
That is end-pipeline improvement but a feature-quality regression at the
MainSearch LGBM level.

**Hard KEEPs (Tier-3):**
- `frag2..6_int` — Olsen loses 223 IDs.
- `best_max_residual_3scan` — MTAC loses 1,468 IDs.

Next: test combined drop of promising candidates (replicating Tier-2
"drop-all-5" finding that individually-borderline drops can be net
positive together).

## drop_clean10 8-file Olsen verification (2026-05-14)

Tested drop_clean10 on the 8-file Olsen Exploris MBR validation set
(rep1+rep2 of E5/E20/E30/E45 conditions) to check whether 2-file
results scale.

Blacklist (10 features): `ms1_m1_intensity, best_gof_5scan,
best_max_residual_5scan, irt_dist_best_gof_5scan,
irt_dist_best_manhattan_5scan, irt_dist_best_max_residual_5scan,
weight_ratio_at_scan, weight_rank_at_scan, top3_ms2_mass_error_mean,
log_by_ratio_m0`

|  | Baseline | drop_clean10 | Δ | % |
|---|---:|---:|---:|---:|
| Final IDs | 424,639 | 421,611 | −3,028 | **−0.71%** |
| Final PGs | 51,873  | 51,787  | −86    | −0.17% |
| pf sum    | 355,948 | 353,167 | −2,781 | −0.78% |

Per-file MainSearch q≤.01:

| File | Baseline | drop_clean10 | Δ | % |
|---|---:|---:|---:|---:|
| 1 | 45,836 | 45,298 | −538 | −1.17% |
| 2 | 44,447 | 43,900 | −547 | −1.23% |
| 3 | 45,501 | 45,421 |  −80 | −0.18% |
| 4 | 45,675 | 45,094 | −581 | −1.27% |
| 5 | 47,846 | 47,455 | −391 | −0.82% |
| 6 | 47,193 | 46,875 | −318 | −0.67% |
| 7 | 40,248 | 39,939 | −309 | −0.77% |
| 8 | 39,202 | 39,185 |  −17 | −0.04% |

Loss is fairly uniform (0.18–1.27% per file). Not hiding a catastrophic
failure on a single file. 3 files (1, 2, 4) lose 1.2% each — real
per-file LGBM regression, not pure MBR artifact.

The 2-file experiment was misleadingly optimistic: 2-file Olsen showed
−0.26% final / −0.62% pf, but the 8-file shows −0.71% final / −0.78%
pf — about 3× worse on final IDs.

**Decision (user): record but don't apply yet. <1% drop is acceptable
if we accumulate more drops; revisit after wider feature sweep.**
