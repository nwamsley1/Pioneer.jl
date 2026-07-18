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

### #gated_hel guard (option 1) — KEAP1 KO

`global_max_gated = global_max * (1 - shadow_hel)` as a SINGLE round-2 feature (ENV `GATED_HEL=1`);
disagreement can only withhold a transfer, never add one (shadow_hel_guard.pdf, option 1).

| config | rows @ q<.01 | uniq | KEAP1→KO leak | round-2 gain |
|---|---|---|---|---|
| baseline (MBR off) | 663,312 | 139,236 | 0 | — |
| shadow + global_max | 705,821 (+6.4%) | 140,957 | 1 | global_max 2.8% |
| shadow + shadow_hel | 714,718 (+7.7%) | 141,474 | 2 | gm 2.8%, shadow_hel 1.5% |
| **GATED_HEL guard** | 698,774 (+5.3%) | 140,167 | 1 | **global_max_gated 2.7%** |
| MULTI + gate-global | 744,208 (+12.2%) | 142,340 | 6 | global_mean_gated + RAW cluster |

- Feature IS used (2.7%, ≈ plain global_max's 2.8%); the single gated column replaced the pair.
- **MULTI + gate-global on KEAP1**: biggest breadth (+12.2%) but leak rises to 6. On only 6 runs
  cluster ≈ global (not enough runs to be condition-scoped), so leaving `cluster_*` RAW reintroduces
  the leak the gate removed from the global channel — the raw cluster channel now carries it. The
  variant is safe *where cluster is genuinely condition-scoped* (many-run EWZ, below); on few-run
  data raw cluster is itself a leak source. (Still 6/142k = noise-level @1% FDR, but relatively up.)
- On KEAP1 the gate is **NOT a win**: it costs breadth (+5.3% vs +6.4%) for **no safety gain** — the
  leak was already at the floor (1). A subtract-only guard can only lose where there's no leak to
  remove. Some true WT-run transfers with imperfect spectral agreement (A<1) got damped too.
- Expected. KEAP1 (6 runs, leak≈floor) is the wrong test for a guard. The decisive test is **EWZ**,
  where shadow+global_max leaks +697 and MULTI leaks +1,182–2,494 — there the gate has real leak to
  cut toward the floor while (ideally) keeping the boost.

### #twin_score alongside global_max (Avenue 1, both grafted)

`TWIN_SCORE=1` → features `[global_max, twin_score, delta_irt]`, graft BOTH global_max and
twin_score (keeps the 1:1 marginal, unlike raw cluster). twin_score = round-1 s1 in the single
most-cosine-similar run. Hypothesis: it's a leak discriminator (a false transfer's twin is a
same-condition run where the precursor is absent → twin_score≈0).

| dataset | variant | IDs | vs base | leak | round-2 twin_score |
|---|---|---|---|---|---|
| KEAP1 | shadow+global_max | 705,821 | +6.4% | 1 | — |
| KEAP1 | **TWIN_SCORE** | 707,251 | **+6.6%** | **1** | 0.9% |
| EWZ | shadow+global_max | 1,935,183 | +4.9% | +931 | — |
| EWZ | **TWIN_SCORE** | 1,951,516 | **+5.8%** | **+1,149** | 2.0% |

- **Safe everywhere, small breadth win**: KEAP1 +6.6% vs +6.4% at the SAME leak floor (1); EWZ +5.8%
  vs +4.9%. Feature IS used (0.9% KEAP1, 2.0% EWZ). Grafting kept it leak-controlled — leak 1 on
  KEAP1, NOT the 6 the raw-cluster variants leaked. The "regularize what you use" design holds.
- **BUT it did NOT reduce leak — hypothesis partly wrong.** EWZ leak rose +931→+1,149 (scales with
  the extra IDs; leak/1000-new-IDs 10.3→10.7, ~flat). twin_score only discriminates CONDITION-CROSSING
  transfer. Neither dataset's residual leak is that kind: KEAP1 has no run structure (twin is a
  coin-flip, can be a WT run); EWZ's residual yeast leak is SYSTEMATIC-REPRODUCIBLE (~480 yeast
  precursors matching human runs in 14–17/20 HO runs — a false transfer's HO-run twin ALSO carries
  the false match → twin_score high, not 0). So the twin can't veto them.
- **Verdict: twin_score is a safe incremental booster (≈ +0.2–0.9 pp over shadow+global_max, leak
  controlled), not a leak-reducing lever on these datasets.** At matched leak it's dominated on EWZ
  by MULTI+gate-global (+9.3% at +1,168 vs twin's +5.8% at +1,149) — but that variant is unsafe on
  few-run data (KEAP1 leak 6), whereas twin_score stays safe across run counts. Reasonable as a small,
  universally-safe upgrade to the shadow+global_max default; not the breakthrough we hoped.

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

### EWZ — gated_hel guard (option 1)

Re-ran the reference set fresh (fresh snapshot; totals shifted vs the doc table above, ordering
identical). Leak = pure-YEAST (`species=="YEAST"`) precursors passing qval<0.01 in a GO113
(Hela-only) run; `analyze_ewz.jl`. The gate is a SINGLE-feature variant (`global_max_gated`
replaces the global_max+shadow_hel pair), so its fair head-to-head is **shadow+global_max**.

| variant | total | Δ% | yeast real (GO114) | yeast LEAK (GO113) | Δleak | leak/real-gain |
|---|---|---|---|---|---|---|
| baseline (MBR off) | 1,844,621 | — | 95,345 | 2,355 | — | — |
| shadow + global_max | 1,935,183 | +4.9% | 102,184 | 3,286 | +931 | 13.6% |
| MULTI uniform graft | 1,983,736 | +7.5% | 104,354 | 5,031 | +2,676 | 29.7% |
| MULTI selective graft | 2,041,759 | +10.7% | 108,766 | 3,729 | +1,374 | 10.2% |
| cluster_max | 2,035,693 | +10.4% | 107,105 | 3,637 | +1,282 | 10.9% |
| **GATED_HEL guard** | 1,898,391 | **+2.9%** | 98,587 | 2,683 | **+328** | **10.1%** |

- **The gate works as designed.** vs shadow+global_max: leak **+931 → +328** (65% cut, toward the
  2,355 floor); boost **+4.9% → +2.9%**. `global_max_gated` round-2 importance 2.5% (feature used).
- **Best transfer PRECISION of any variant**: leak per new-true-yeast 13.6% → **10.1%** — the gate
  preferentially removes FALSE transfers (high global_max, low agreement), exactly the PDF claim.
  Ties selective-graft (10.2%), crushes uniform graft (29.7%).
- **BUT as a single feature it caps breadth** (+2.9% is the lowest boost) — the subtract-only gate
  also damps true transfers with imperfect spectral agreement (real yeast 102,184 → 98,587).
- **Verdict: the gate is a pure safety/precision lever, not a breadth lever.** On its own it is the
  SAFEST cross-run variant but the least complete. The high-value use is **gating the leaky global
  channel WITHIN the multi-feature model** (gate global_max/global_mean, leave cluster_* raw) — that
  should keep selective-graft's +10.7% breadth from the leak-safe cluster features while applying the
  precision gate only where leakage originates. NEXT: MULTI + gated-global variant.

### EWZ — MULTI + gate-global (ENV `MULTI_FEATURE=1 GATE_GLOBAL=1`)

Gate ONLY the global channel inside the multi-feature model: feed `global_{max,mean}_gated =
global_*·(1-shadow_hel)` + RAW `cluster_{max,mean}` + `delta_irt`; graft the gated global cols only.

| variant | total | Δ% | yeast real | yeast LEAK | Δleak | leak/real-gain |
|---|---|---|---|---|---|---|
| baseline (MBR off) | 1,844,621 | — | 95,345 | 2,355 | — | — |
| GATED single feature | 1,898,391 | +2.9% | 98,587 | 2,683 | +328 | 10.1% |
| MULTI selective graft | 2,041,759 | +10.7% | 108,766 | 3,729 | +1,374 | 10.2% |
| cluster_max | 2,035,693 | +10.4% | 107,105 | 3,637 | +1,282 | 10.9% |
| **MULTI + gate-global** | 2,015,772 | **+9.3%** | 107,472 | 3,523 | **+1,168** | **9.6%** |

- **Fixes the single gate's breadth collapse**: +2.9% → **+9.3%** (the RAW cluster features carry the
  breadth; the gate only precision-filters the global channel). Round-2 leans on `global_mean_gated`
  (1.9%) + cluster features; `global_max_gated` small (0.1%).
- **Safer than multi_sel AND better precision**: leak +1,168 vs +1,374 (−15%), leak/real-gain 9.6%
  vs 10.2% — the best precision of the high-breadth variants.
- **The tradeoff is real but modest**: gives up ~1.4 pp breadth (26k rows) vs multi_sel for the
  safety/precision gain. It sits exactly where predicted — between multi_sel (max breadth) and the
  pure gate (max safety). Good default when false transfer is the priority; multi_sel when raw
  completeness is.
