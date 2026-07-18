# Shadow + global_max — improvements progress

Branch `feat/shadow-global-max-improvements`, off `feat/global-max-shadow-scoring-on-develop`.
Goal: improve the strongest two-round variant so far (shadow-decoy regularized `global_max`)
on **data completeness** (precursors quantified in more runs) while keeping cross-run false
transfer near zero.

## Where we are (baseline validation)

`shadow + global_max` (their branch defaults: `KNN_SCOPE=:all_runs`, `KNN_AGG=:max`,
`SHADOW_DECOY_MODE=:symmetric`), run on the KEAP1 KO data (6 files: 3 H292 WT + 3 KEAP1-KO,
`Pioneer_Human_canon_std.poin`, MBR off, `TWO_ROUND=1`):

| config | rows @ q<.01 | unique | KEAP1→KO leak |
|---|---|---|---|
| predicted baseline (develop, MBR off) | 660,369 | 137,859 | 0 |
| shadow + global_max | 701,304 (+6.2%) | 139,563 (+1.2%) | 1 |

- Gain is **breadth, not discovery** (+6.2% rows vs +1.2% unique) — precursors identified in
  more runs. That is the win we care about.
- The **shadow regularization controls `global_max`'s leakage**: 1 KEAP1→KO leak, vs 4 for the
  unregularized global_max in our offline harness (naive irt_knn_max, KEAP1 6 runs). Matches
  our offline shadow/twin result (+1.4% rows, 0 leak).
- Round-2 top gain is `smoothed_2d_shadow_hellinger` (17.5%), not the cross-run feature — the
  1:1-marginal shadow regularization is holding; the cross-run feature is used jointly, not as a
  standalone lever.
- Their own FTR-27 (27 files, 3-species) numbers: shadow+global_max +16.6% correct IDs at <15%
  feature reliance — the raw-count payoff needs many-run cross-condition structure; KEAP1's 6
  files are the *leakage* stress test, not the count contest.

Note: `global_max` is written at the **fold-file** level ("11 other runs / 12 files"), but this
is NOT a leak — precursors are uniquely fold-assigned, so a precursor appears once per real run
and `max over other fold-files` == `max over its other-run instances`.

## Improvement plan (priority order)

1. **Cross-run shadow-spectrum agreement feature** (in progress). Alongside `global_max`, add the
   **donor↔acceptor top-8 smoothed-fragment Hellinger**: for each (precursor, run), compare this
   run's shadow spectrum to the best-donor run's (the run supplying the global_max). A *true*
   transfer (same real peptide) → high spectral agreement; a *false* transfer (KEAP1 into a KO
   run) → low agreement (KO acceptor is noise). This is the signal `global_max` lacks — it is
   *why* global_max leaks. As a second shadow-regularized mbr feature it lets round-2 learn
   "high global_max AND high agreement → real; high global_max but low agreement → decoy",
   buying breadth without leak. Reuses the top-8 `frag*_smoothed_intensity` columns + Hellinger
   from our empirical-library work.
2. **Condition-aware scope** — `global_max` restricted to the run's abundance-correlation cluster
   (z-corr→SVD→modularity, see reference:run_grouping_method) so one cross-condition run cannot
   seed the transfer. Complements #1.
3. **Multi-feature round-2** — `global_max` + `shadow_hellinger` + `cluster_max`, all
   shadow-regularized; let the LGBM weigh them.

## Context we must not lose (from the broader investigation)

- **False-transfer axis is the differentiator.** naive cross-run features leak (global-max is the
  MBR trap); the shadow-decoy (symmetric twin) is the regularizer that caps leakage.
- **Shadow-decoy == our offline `offline_score_twin.jl`**: per real target → shadow_decoy
  (nearest-iRT decoy's columns + target's mbr_score, labeled decoy); per real decoy →
  shadow_target. Makes the mbr_score marginal exactly 1:1 → forces joint use with base features.
- **DIA-NN's MBR** = build empirical library from *confident* IDs, refine iRT from best run,
  re-search the *predicted* library (2nd pass). It transfers a *template + iRT prior*, requiring
  run-B evidence — NOT a score. Our cross-run *score* features can transfer without evidence,
  which is why they leak and need the shadow regularization.
- **Empirical-library thread** (separate branch `feat/empirical-library-isotopes`): searching an
  observed-intensity library recovers ~86% of predicted-lib unique precursors with a two-library
  swap + raw-intensity read; mutation decoys break its FDR (empirical targets re-match trivially,
  mutation decoys can't), observed-derived decoys are the only self-consistent null. Parked in
  favor of this stronger shadow-scoring line.
