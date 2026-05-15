# feature/mbr-batch-f — campaign results summary (2026-05-15)

End-to-end summary of the three layered campaigns landed on `feature/mbr-batch-f`:

1. **Feature reduction** (Tier 2–5): 66/71 → 47/52 LGBM features
2. **MainSearch + Precursor Scoring perf optimization**: −70% on Main Search wall time
3. **MBR rtv1–rtv4** (companion Pioneer_MBR_CodeChanges_2026-05-14.md): new MBR features + slim 12-feature FTR LGBM

## Headline: YeastMBR rtv4 validation on top of perf-optimized branch

20-file (Combined) + 10-file (HelaOnly) using `altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted_entr1.poin`, `PIONEER_FTR_MODE=pep`, `PIONEER_FTR_ALPHA=0.005`.

| Metric | Doc baseline (9749589f) | Doc rtv4 | **Our HEAD (2150a2f6)** |
|---|---:|---:|---:|
| Precursor any-YEAST Dennis FTR | 3.18% | 2.87% | **1.63%** |
| Precursor strict-YEAST FTR | 2.54% | 2.20% | **1.17%** |
| ΔTotal q-passing (B−A, precursor) | 6,133 | 7,149 | **8,523** |
| Combined wall time (20-file + HelaOnly) | ~56 min | ~56 min | **16.4 min** |

**Net effect: dramatically lower species-mismatch FTR + more legitimate transfers + 3.4× faster.**

## Per-step wall-time deltas (8-file Olsen, MBR=on)

| Step | Pre-trim baseline | HEAD | Δ |
|---|---:|---:|---:|
| Main Search | 442.04s | **134.36s** | **−70%** |
| Precursor Scoring | 50.14s | **29.37s** | **−42%** |
| Total Runtime | 639s | **316s** | **−51%** |

vs develop branch (no MBR Batch F + no feature additions) on the same 8-file MBR=on:
- Develop: 114s total (lightweight, no MBR machinery active on small dataset)
- HEAD: 235s total — 2× slower, but **+18.8% IDs / +5.8% PGs** delivered

## Feature counts

| Set | Pre-campaign | HEAD |
|---|---:|---:|
| MainSearch PRESCORE_FEATURES | 66 | **47** |
| ScoringSearch ADVANCED_FEATURE_SET | 71 | **54** (52 + 2 new rtv MBR) |
| FTR_FEATURES_F_TRUE (MBR-FTR LGBM) | 51 | **12** |

## Threading touched

Six per-precursor / per-PSM passes in MainSearch + 1 in PrecursorScoring now use `Threads.@threads :static`:
- `_add_neighborhood_features_fused!` (3/5/11-scan windows fused into one pass)
- `_add_fragment_chromatogram_features!`
- `_add_ms1_chromatogram_features!`
- `add_ms1_features!` per-PSM lookup (via `Threads.@spawn` chunks)
- `compute_mbr_features_dual!` per-PSM loop
- `_add_neighborhood_features_window!` (now superseded by fused variant)

Each writes to disjoint row indices; bit-identical results vs single-threaded.

## Code changes summary

| File | Lines net | What |
|---|---:|---|
| `MainSearch/scoring.jl` | minor | SCORING_LGBM_HP 600/0.033 → 200/0.10; `compute_infold` plumbing |
| `MainSearch/features.jl` | −500 | M+1 machinery deletion + trim 3 feature builders + threading + sortperm/segmented walks + Float16 outputs + fusion |
| `CommonSearchUtils/fusedMatch.jl` | −30 | M+1 capture in `apply_main_scoring!` deleted |
| `PSMs/UnscoredPSMs.jl`, `PSMs/ScoredPSMs.jl` | −30 | M+1 fields + `log_by_ratio_m01` removed |
| `PrecursorScoringSearch/mbr_features.jl` | +60 | rtv1+rtv2 features + threading |
| `PrecursorScoringSearch/mbr_ftr.jl` | −10 | FTR feature list 51 → 12 (rtv3+rtv4) |
| `PrecursorScoringSearch/model_config.jl` | −15 | ADVANCED_FEATURE_SET 71 → 54 |
| `PrecursorScoringSearch/score_psms.jl` | +5 | `compute_infold = match_between_runs` + persist `:trace_prob_infold` |
| `BitVecCalibration/BitVecCalibration.jl` | +5 | env-var override `PIONEER_BITVEC_MIN_EXCESS_RATE` (default unchanged) |

## Open items for follow-up

1. **Per-file species-prior MBR feature** — from `BATCH_F_PLAN.md` open questions; the structural fix for the ~2% strict-YEAST Dennis FTR floor.
2. **Tighter `MBR_DONOR_Q_THRESHOLD`** (0.01 → 0.001) — only confident donors qualify.
3. **MTAC MBR=on validation** confirmed 18.8% more IDs vs develop at 2× wall time on 6-file run; full 23-file pre-merge validation would solidify the IDs/wall-time tradeoff.
4. **PrecursorScoringSearch `apply_mbr_filter_paired!`** — still single-threaded but mostly LGBM-internal (~3-5s on 8-file). Marginal gain available.

## Commits

Final commits on top of pre-campaign tip `780d32b3` (merge-base with develop) — push 1 (perf trims), push 2 (feature reduction), push 3 (rtv1-4 MBR work). Final HEAD: **`2150a2f6`**.
