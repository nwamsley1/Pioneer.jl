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

## Results log

### #1 cross-run shadow-spectrum agreement (`<KNN_COL>_shadow_hel`) — KEAP1 KO

| config | rows @ q<.01 | unique | KEAP1→KO | round-2 gain |
|---|---|---|---|---|
| predicted baseline | 660,369 | 137,859 | 0 | — |
| shadow + global_max | 701,304 (+6.2%) | 139,563 | 1 | global_max 2.8% |
| + shadow-spectrum agreement | 710,222 (+7.6%) | 140,075 | 2 | gm 2.8%, **shadow_hel 1.5%** |

- **Completeness +1.3% rows** over plain shadow+global_max (the win we want); feature IS used (1.5% gain).
- Did **NOT** reduce leakage (1→2, both noise-level @1% FDR over ~119k). The LGBM used the agreement
  feature as another *borrow* signal, not a *reject* guard, as hoped.
- Modest because KEAP1 has only 6 runs (few donors, thin cross-run spectral structure). The real
  test is a many-run dataset (FTR-27 / 60-file yeast KO) where the agreement signal is rich.
- FDR-A (real-only) 722,596 vs FDR-B (shadow-inclusive) 698,093.

### #shadow_hel on many-run (3P) + MULTI_FEATURE

3P (OlsenExploris 3-proteome, 23 files, mixture → low leakage risk):

| variant | rows | vs base | unique |
|---|---|---|---|
| baseline | 1,335,695 | — | 79,168 |
| shadow + global_max | 1,418,289 | +6.2% | 80,331 |
| + shadow_hel | 1,437,599 | +7.6% | 80,712 |
| MULTI_FEATURE (all 5) | 1,445,999 | +8.3% | 81,057 |
| cluster_max (unregularized, earlier) | 1,474,289 | +10.4% | 81,412 |

KEAP1 MULTI_FEATURE: 710,943 rows (+7.7%), leak 2 — ≈ +shadow_hel; cluster features inert (6 runs → cluster≈global).

Findings:
- shadow_hel adds a consistent ~+1.3–1.4% breadth on BOTH KEAP1 and 3P.
- MULTI_FEATURE is SAFE (KEAP1 leak held at 2 with 5 cross-run features; no feature dominated) and
  AUTO-SELECTS (3P gains: global_mean 2.5%, cluster_mean 1.3%, global_max 1.3%, shadow_hel 1.2%; the
  model preferred mean/cluster on the mixture, global on KEAP1). It beat pure global_max (+6.2%→+8.3%).
- BUT the uniform 1:1 shadow graft OVER-DAMPS: +8.3% still trails unregularized cluster_max (+10.4%).
  The regularization that buys leak-safety caps the completeness ceiling on low-leakage data.
- On 3P (mixture) cluster_max is the completeness leader; on KEAP1 (true absence) the shadow-regularized
  global variants are the safe choice. Dataset-dependent.

NEXT (selective regularization): graft the shadow only onto the LEAKAGE-PRONE global features
(global_max/global_mean), NOT the condition-scoped cluster features (cluster_max/cluster_mean) —
those are already leak-resistant by scope. Should let the model lean fully on cluster_* (recover
toward +10.4%) while keeping global_* regularized (safe). i.e. regularize by leakage risk, not uniformly.

TODO also: (b) make shadow_hel a *guard* not a booster (gate global_max by agreement: global_max×(1−hel));
(c) condition-aware scope is now IN via cluster_* — tune KNN_K / try hard z-corr→modularity clusters.

## Improvement plan (priority order)

1. **Cross-run shadow-spectrum agreement feature** (DONE, KEAP1 — see results log; modest, needs
   many-run test). Alongside `global_max`, add the
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

### EWZ (yeast/human two-proteome false-transfer test) — DECISIVE

GO113 = Hela-only (no yeast); yeast passing there = FALSE TRANSFER (baseline floor 5,792).
GO114 = Hela+Yeast (real yeast). Boost = total/GO114-yeast; Safety = GO113-yeast leak.

| variant | total | yeast real (GO114) | yeast LEAK (GO113) | leak Δ |
|---|---|---|---|---|
| baseline (MBR off) | 1,836,800 | 98,877 | 5,792 | — |
| shadow + global_max | 1,923,346 (+4.7%) | 105,395 | 6,489 | +697 |
| MULTI uniform graft | 1,971,421 (+7.3%) | 107,621 | 8,286 | +2,494 |
| MULTI selective graft (GRAFT_SCOPE=global) | 2,028,424 (+10.4%) | 112,033 | 6,974 | +1,182 |
| cluster_max (earlier) | 2,023,330 (+10.2%) | 110,516 | 6,869 | +1,077 |
| global/irt_knn NO shadow (MBR trap) | 2,025,997 (+10.3%) | 110,598 | 10,366 | +4,574 |

WINNER: **MULTI selective graft** — highest boost (+10.4%, beats cluster_max) at controlled leak
(+1,182, less than half the unregularized trap's +4,574). Validates "regularize by leakage risk,
not uniformly": grafting only the leakage-prone global_* (leaving condition-scoped cluster_*
ungrafted) lets the model boost from the leak-SAFE cluster features, so it needn't over-rely on
the grafted-but-leaky global features. Selective OUT-SAFES uniform (+1,182 vs +2,494) AND boosts
more (+10.4% vs +7.3%). shadow+global_max confirms the shadow guard works (MBR trap leak 10,366 -> 6,489).

RECOMMENDATION: MULTI_FEATURE + GRAFT_SCOPE=global is the strong default — dataset-adaptive
(cluster on mixtures, global where safe), maximum completeness at controlled false transfer.
