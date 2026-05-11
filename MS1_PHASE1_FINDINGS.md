# MS1 Phase 1 — Findings Log

Branch: `feature/ms1-phase1`. Updated 2026-05-11.

Tracks experiments and decisions from the MS1 phase-1 + feature-engineering +
deconv-tuning work on the non-entrap Olsen Exploris test set
(`/Users/nathanwamsley/Data/RegressionTestsLite/Olsen_3P_Exploris_one`,
2 files, library `altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted.poin`).

Baseline reference: `nonentrap_pmm_paircomp` = **62,451 IDs @ q≤.001 / 86,047 @ q≤.01 / 12,199 protein groups @ q≤.01**.

---

## 1. Feature work (committed)

### MS1 point-lookup features (`d4ba4b72`, `8e3651c5`, `883d0de1`)
Per-PSM lookup of the nearest MS1 scan for M0 and M+1 (charge-aware m/z, ±ppm window).
Adds `ms1_m0_intensity`, `ms1_m1_intensity`, `ms1_m0_mass_err_ppm`. M+2 / M+3 were trialed
and pruned — kept the top 7 most-informative MS1 fields.

### MS1 chromatogram correlation features (`a04fa661`)
Per-precursor across-scan correlations: `ms1_corr_weight_m0`, `ms1_corr_m0_m1`,
`ms1_apex_offset_irt`, `ms1_weight_apex_to_m0_apex_irt`.
**These rank top 5 in feature importance** (gain 7,900–8,500 each).

### MS2 fragment chromatogram features (`3a0e9129`)
Per-rank top-6 fragment intensities (`frag1_int..frag6_int`, captured inline by
`apply_main_scoring!` during the fused scan loop), plus per-precursor frag-chromatogram
correlations (`frag_corr_top1_top2`, `frag_corr_top1_top3`, `frag_corr_top1_weight`,
`frag_corr_mean_pairwise`, `frag_corr_min_pairwise`, `frag_corr_top3_weight`,
`frag_apex_dispersion_irt`, `n_correlated_fragments`). Validated 2026-05-10 to add
~+2,088 IDs at q≤.01 vs MS1-only baseline (entrap1, paired EFDR ~0.0107).
**`n_correlated_fragments` is the #2 most-important feature overall** (gain 13,188).

### Smoothness feature (`7fe47402`)
Σ(second-derivative)² of the weight chromatogram per precursor; large = bumpy
(decoy/noise), small = real peak. Used by PrecursorScoringSearch only (gain 849).

### Feature-set pruning (this branch, post-importance analysis)
Dropped low/zero-gain features:
- From `PRESCORE_FEATURES`: `:rank1_matched`, `:top3_matched`, `:top5_matched` (all gain ≤ 187).
- From `ADVANCED_FEATURE_SET`: `:charge` (119), `:b_count` (94), `:percent_theoretical_ignored` (0).
- **Kept `:Mox`** (gain 91) per explicit request — useful even at low gain.

### Features tried and rejected
- **MS1↔fragment direct cross-correlations** (`ms1_corr_m0_fragsum`,
  `ms1_corr_m0_frag1`, `ms1_corr_m0_best`, `ms1_apex_offset_fragsum_irt`):
  −3,325 IDs at q≤.001 vs baseline, +477 at q≤.01. Mixed; reverted (`b3de15fe`).
  Mechanism likely too correlated with existing chromatogram features to add
  high-confidence signal.

---

## 2. Deconvolution + scoring infrastructure

### Pair competition (`00870d93`)
Drop the lower-scoring partner of each target-decoy pair before computing q-values
and PEP. Default-on (`PIONEER_PAIR_COMPETITION` defaults to `"1"`; set `"0"` to disable).
Lives in `MainSearch/MainSearch.jl:290–...`. Ratio dropped at q≤.001 is ~3:4 target:decoy
(decoys lose more, confirming target advantage). Validated to win q≤.001 by a big margin
and lose q≤.01 by ~1% vs the no-paircomp baseline.

### LASSO deconvolution solver (`4037fa0b`)
ENV-gated alternative to default PMM. Set `PIONEER_DECONV_SOLVER=lasso` and
optionally `PIONEER_LASSO_LAMBDA=0.01`. Validated:
- Best LASSO + keepzero loses 1,267–2,737 IDs vs PMM at headline metric.
- Not recommended as default; available for experimentation.

### 2-pass refinement (`7fe47402`)
Re-run library_search after PEP/qval filter, retrain LGBM. Off by default
(`PIONEER_TWO_PASS_QVAL` / `PIONEER_TWO_PASS_PEP`). Validated: marginal lift not
worth the complexity vs paircomp alone — paircomp won q≤.001 by a big margin,
lost q≤.01 by only 1%, so paircomp made default.

### MS2 ppm floor (reverted, `b3de15fe`)
`PIONEER_MS2_MIN_PPM` was an env-gated floor on the intensity-aware mass tolerance.
Bit-identical to baseline in one test (62,451 / 86,047 / 12,199), suggesting calibrated
tolerances never went below 10 ppm. Reverted with the xcorr features.

### Hyperparameter tuning (`04c68246`)
`SHARED_LGBM_HP` reduced to production-realistic 200 iter / lr=0.10.

---

## 3. Per-method PSM type refactor (`017ffe94`..`e12fd00d`)

**Goal:** stop tuning paths from carrying MainSearch-only fragment-chromatogram
captures (`frag1..6_int`, `matched_rank_mask`, `rank1/top3/top5_matched`).

**Audit finding:** the per-scan write cost was already zero in tuning paths thanks
to the `record_match!(::FusedQuadEst/RTIndexed, ...)` no-op overloads. The remaining
waste was in `Score!`: tuning still produced `MainSearchScoredPSM` rows with the
MainSearch-only fields zero-initialized and copied, plus the wider struct allocated
per-thread.

**Refactor:** introduced parallel slim type families:
- `MainUnscoredPSM` / `MainSearchScoredPSM` — full (frag captures, rank-mask, rank1/3/5_matched).
- `TuningUnscoredPSM` / `TuningScoredPSM` — slim, drop those fields.

`record_match!` dispatches on `(FusedKind, eltype(unscored_psms))`:
- `FusedStandard` + `MainUnscoredPSM` → `apply_main_scoring!`
- `FusedStandard` + `TuningUnscoredPSM` → `apply_tuning_scoring!` (no frag/mask stores)
- `FusedQuadEst` / `FusedRTIndexed` over any UnscoredPSM → no-op (unchanged)

`SimpleLibrarySearch` per-thread scratch now carries both families. Each search method
picks via `get_scored_psms` / `get_unscored_psms` dispatch on params type.
ParameterTuning / QuadTuning / IntegrateChromatograms all route to the tuning family.

**Verification:** ecoli integration test bit-identical (2,371 PSMs / 1,341 protein groups).
Olsen Exploris with the refactor + the xcorr features produced bit-identical results to
the same code without the refactor — refactor is a true behavioral no-op.

**Memory cost:** ~400 KB extra per thread for the parallel buffer set; one-time
allocation, not per-scan.

---

## 4. Currently active defaults (post-refactor, post-pruning)

| Knob | Default | Override |
|---|---|---|
| Deconv solver | **PMM** | `PIONEER_DECONV_SOLVER=lasso\|alasso\|ols` |
| Pair competition | **ON** | `PIONEER_PAIR_COMPETITION=0` to disable |
| 2-pass refinement | OFF | `PIONEER_TWO_PASS_QVAL=<frac>` / `PIONEER_TWO_PASS_PEP=<frac>` |
| Keep zero weight | OFF | `PIONEER_KEEP_ZERO_WEIGHT=1` |
| MS2 ppm floor | OFF (no-op in practice) | (env var removed; was `PIONEER_MS2_MIN_PPM`) |
| iRT σ-tolerance multiplier | 3 | `PIONEER_IRT_SIGMA_TOL=4` (validated; tried but didn't lift) |
| Max frag rank | 255 (no cap) | `PIONEER_MAX_FRAG_RANK=6` (top-6 only) |

---

## 5. Feature importance (mean gain, Olsen Exploris, post-refactor)

Top 20 over both stages (MainSearch per-file LGBM = 4 calls; PrecursorScoring
experiment-wide LGBM = 2 calls; shared features get 6):

| Feature | mean gain | stage |
|---|---:|---|
| `best_gof_3scan` | 56,417 | MainSearch |
| `n_correlated_fragments` | 13,188 | MainSearch |
| `ms1_corr_weight_m0` | 8,425 | MainSearch |
| `ms1_corr_m0_m1` | 7,901 | MainSearch |
| `max_matched_residual` | 7,616 | both |
| `frag_corr_mean_pairwise` | 6,615 | MainSearch |
| `frag_corr_top1_top3` | 3,902 | MainSearch |
| `frag_apex_dispersion_irt` | 3,780 | MainSearch |
| `ms1_weight_apex_to_m0_apex_irt` | 3,244 | MainSearch |
| `irt_error` | 3,163 | both |
| `best_max_residual_3scan` | 3,131 | MainSearch |
| `frag_corr_top3_weight` | 3,113 | MainSearch |
| `frag_corr_top1_top2` | 2,967 | MainSearch |
| `frag_corr_top1_weight` | 2,931 | MainSearch |
| `spectrum_peak_count` | 2,860 | MainSearch |
| `frag_corr_min_pairwise` | 2,765 | MainSearch |
| `best_manhattan_3scan` | 2,630 | MainSearch |
| `ms1_m0_mass_err_ppm` | 2,618 | MainSearch |
| `frag1_int` | 2,407 | MainSearch |
| `fitted_manhattan_distance` | 2,331 | MainSearch |

**Caveat:** LightGBM is multithreaded with `bagging_fraction=0.25` and
`feature_fraction=0.5`, so run-to-run variance at q≤.001 can be on the order of
a few thousand IDs. Headline metrics (q≤.01, protein groups) are more stable.

---

## 6. Running the pipeline from the terminal

Copy `dev/tuning_config_template.json` to a working dir and edit `paths`.
Then:

```bash
# Defaults — PMM solver, paircomp on
julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/path/to/my_config.json")'

# Try LASSO at λ=0.01
PIONEER_DECONV_SOLVER=lasso PIONEER_LASSO_LAMBDA=0.01 \
  julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/path/to/my_config.json")'

# Disable paircomp (e.g. for baseline comparison)
PIONEER_PAIR_COMPETITION=0 \
  julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/path/to/my_config.json")'

# 2-pass refinement at PEP filter 0.75
PIONEER_TWO_PASS_PEP=0.75 \
  julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/path/to/my_config.json")'
```

Reading out the IDs:

```julia
using Arrow, DataFrames
df = DataFrame(Arrow.Table("/path/to/results/precursors_long.arrow"))
count((df.global_qval .<= 0.001) .& df.target)    # q≤.001
count((df.global_qval .<= 0.01)  .& df.target)    # q≤.01

pg = DataFrame(Arrow.Table("/path/to/results/protein_groups_long.arrow"))
count((pg.global_qval .<= 0.01)  .& pg.target)    # PGs @ q≤.01
```
