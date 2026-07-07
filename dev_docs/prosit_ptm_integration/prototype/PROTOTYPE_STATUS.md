# Localization-FLR scoring — prototype status (2026-07-07, branch feature/localization-flr-scoring)

Built and validated the core of the localization-scoring plan
(`../localization_flr_model_plan.pdf`) while unattended. **Nothing here is wired into the Pioneer
package** — it is a validated prototype (analysis-only) ready to be promoted.

## What's implemented + verified

1. **Core scoring algorithm** — `loc_score_core.jl` + `test_loc_score_core.jl` (**15 tests, all pass**).
   - General "distinguishing = cumulative mod content differs" definition (handles competitors that
     differ at >1 site).
   - `s(A,B)` competitive ratio on presence/absence; `S_score` = min over competitors.
   - Competitor sets per the settled design: **target vs its siblings; decoy vs ALL targets incl its
     parent (kept).**
   - Per-site accumulation with the joint-only flag (2b(ii)).
   - **The tests reproduce the plan's `S A S A S A K` worked example exactly**: target `A`→1.0, decoy
     `A'` (parent kept)→0.0, decoy `A'` (parent dropped)→1.0 escape, mislocalized target→0.0.

2. **Real-data feature extractor** — `~/prosit_phospho_test/extract_sds_features.jl`.
   - Runs on the combined 12-file search (`dilution_curve2/combined`, lib `yps_decoy2.poin`).
   - Computes b/y fragment m/z from sequence+mods (validated: b2/y2 match hand calc), matches to the
     observed spectrum at each instance's PSM scan (`scan_idx` verified 200/200 -> correct MS2 scan
     with precursor in isolation window), 20 ppm.
   - Emits per `C>1` instance: `S_vs_targets`, `self_frac`, `frac`, `iso_group_size`, `C`, label.

3. **LGBM ablation** — `~/prosit_phospho_test/ablation_lgbm.jl`.

## Results (10k balanced sample: 5,561 targets / 4,438 decoys)

**S_vs_targets separates targets from decoys, strongest for localizable (low-C) strata** (as designed):

| C | target median S | decoy median S |
|---|---|---|
| 2 | 0.80 | 0.15 |
| 3 | 0.59 | 0.17 |
| 4-6 | 0.36 | 0.15 |
| 7+ | 0.18 | 0.14 |

**Ablation (70/30 split; "FLR" = sample D/R, ~1:1, so *comparison* is meaningful, not the absolute):**

| features | AUC | test IDs @10% | @5% |
|---|---|---|---|
| frac | 0.643 | 4 | 4 |
| frac+C | 0.647 | 10 | 1 |
| frac+C+S | 0.68 | **286** | **191** |

=> the site-determining S evidence is the key axis; frac+C alone cannot control FLR, frac+C+S recovers
1-2 orders of magnitude more IDs at matched decoy-FLR.

## Caveats / honest limitations

- **The "FLR" here is a balanced-sample `D/R`, NOT the calibrated per-ID FLR.** Real calibration needs
  the true target:decoy ratio and, crucially, the **standards ground-truth check** (predicted FLR vs
  empirical error) — the Phase-2 decision gate — which is **not yet done**.
- **Leakage caught + removed:** `n_comp` (competitor count) leaks the label (`C` for a decoy, `C-1`
  for its target); with it in, AUC was a bogus 0.999. Dropped it. `S_vs_targets` itself is computed
  with the asymmetric competitor set (decoy killed by parent) — that is the *intended* design (decoy
  models a mis-localized target), and the modest 0.68 AUC (heavy target/decoy overlap in high-C)
  confirms it is a real, imperfect localization signal, not a pure label-encoder.
- **Per-site features (2b(ii)) are in the core but NOT yet in the extractor/LGBM** — only the aggregate
  `S_vs_targets` (min) is used. Wiring `min_site_support`, `n_sites_joint_only` in is a next step.
- Single mod type (phospho) on `yps_decoy2.poin`. Oxidation (`yps_std.poin`) needs the multi-type m/z
  (mod-mass per type) in the extractor.
- Sampled (10k of ~95k); scoring is fast, scale-up is just runtime.

## Next steps (Phase 2 gate)

1. **Standards calibration** — score the Coon standards' `C>1` peptidoforms, compare per-ID FLR to the
   known-site empirical error (the definitive test the sample-D/R cannot give).
2. Add per-site S features to the extractor + LGBM (expect a further lift on multi-site strata).
3. Proper per-ID FLR from the model score at the real target:decoy ratio; MBR handling; report.
4. Promote to package code once calibrated.

## How to run

```
julia dev_docs/prosit_ptm_integration/prototype/test_loc_score_core.jl        # core tests
# from ~/prosit_phospho_test:
julia --project=Pioneer.jl extract_sds_features.jl dilution_curve2/combined/results yps_decoy2.poin 10000 out.arrow
julia --project=Pioneer.jl ablation_lgbm.jl out.arrow
```
