# CLAUDE.md - ParameterTuningSearch

Guidance for Claude Code when working with `ParameterTuningSearch` — the first search method in the SearchDIA pipeline, responsible for fitting per-file fragment mass-error and RT-to-iRT alignment models before the rest of the pipeline runs.

## Algorithm: scout → collection → MS1 fit

The historical "3-phase bias-shift × multi-score × scan-scaling" iteration was replaced by a simpler scout-then-collect flow, followed by an MS1 mass-error fit using the MS2-accepted PSMs:

```
┌─────────────────────────────────────────────────────────┐
│ initialize_models!                                      │
│  - Mass error model: bias=0, tol=±WIDE_SCOUT_TOL_PPM   │
│    (= 100 ppm; placeholder, overwritten by phase 1)    │
│  - Quad transmission: SquareQuadModel(0.0)             │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ Phase 1: Wide Scout (MS2)                               │
│  Mass model: ±WIDE_SCOUT_TOL_PPM (100 ppm hardcoded)   │
│  Target PSMs: WIDE_SCOUT_TARGET_PSMS                    │
│  Initial scans: WIDE_SCOUT_INITIAL_SCANS                │
│  Top-N peaks: TUNING_TOPN_PEAKS                         │
│  → Fit bias + Gaussian σ from wide-scout fragments      │
│  → Replaces the placeholder mass error model            │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ Phase 2: Collection (MS2)                               │
│  Mass model: from Phase 1 fit                           │
│  Target PSMs: TUNING_MIN_SAMPLES (1200)                 │
│  Initial scans: TUNING_MIN_COLLECT_SCANS (5000)         │
│  → Refine model with PSMs collected at the fitted tol  │
│  → Fit final RT alignment spline                        │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ Phase 3: MS1 Mass-Error Fit                             │
│  Input: MS2-accepted PSMs from Phase 2                  │
│  For each PSM: find nearest MS1 scan by RT, match       │
│  M+0/M+1/M+2 within ±MS1_DIAG_WINDOW_PPM (100 ppm) of   │
│  the predicted precursor m/z; record ppm residuals.     │
│  Fit MassErrorModel(median, ±MS1_DIAG_TOL_K_MAD·MAD)    │
│  (k = 5). Install via setMs1MassErrorModel! and write   │
│  qc_plots/ms1_mass_error_plots/{file}.pdf.              │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
              Store models in SearchContext
```

All tolerance/scan/score thresholds are constants in `constants.jl` (`WIDE_SCOUT_TOL_PPM`, `WIDE_SCOUT_TARGET_PSMS`, `WIDE_SCOUT_INITIAL_SCANS`, `TUNING_TOPN_PEAKS`, `TUNING_MIN_SAMPLES`, `TUNING_MIN_COLLECT_SCANS`, `TUNING_SCORE_TIERS`). None are user-tunable — the historical `parameter_tuning.*` JSON section was removed; `init_mass_tol_ppm_guess` was always shadowed by the wide scout's hardcoded 100 ppm anyway.

## Files

| File | Purpose |
|---|---|
| `ParameterTuningSearch.jl` | Main implementation: `initialize_models!`, scout/collect phases, `store_final_results!`, summarize_results! |
| `types.jl` | `ParameterTuningSearchParameters` (no JSON-fed fields), `ParameterTuningSearchResults`, `IterationState`, `ParameterTuningDiagnostics`, `ParameterHistory` |
| `utils.jl` | Per-file helpers: `library_search`, `add_tuning_search_columns!`, `filter_and_score_psms!`, `fit_mass_err_model`, `fit_irt_model`, `fit_intensity_mass_error.jl`. Operates on a `Vector{MassErrSample}` (12-byte struct) rather than the deleted `Vector{FragmentMatch}`. |
| `constants.jl` | All hardcoded tuning constants (replaces the deleted `parameter_tuning` JSON section). |
| `fit_intensity_mass_error.jl` | Optional `IntensityMassErrorModel` fitting pipeline. Defines `EMPK_QUANTILE = 0.95` and `EMPK_CLAMP_HI = 4.5f0` for the coverage-k fit. |
| `ms1_diagnostic.jl` | MS1 mass-error fit: `collect_ms1_residuals` (per-PSM ±100 ppm isotopologue match), `fit_ms1_model_from_residuals` (median + 5·MAD), `generate_ms1_residual_histogram`. Constants `MS1_DIAG_WINDOW_PPM = 100.0f0`, `MS1_DIAG_TOL_K_MAD = 5.0f0`. |

## Mass error fitting

`fit_mass_err_model` consumes a `Vector{MassErrSample}` written by `run_fused_masserr!` (the ParameterTuning variant of `run_fused!` in `CommonSearchUtils/fusedMatch.jl`). Each sample is 12 bytes:
```julia
struct MassErrSample
    theoretical_mz::Float32
    observed_mz::Float32
    intensity::Float32
end
```
This replaced the 32-byte `FragmentMatch` per-thread buffer (~2.7× smaller; 320 KB → 120 KB per thread on Olsen Astral).

## Convergence

A phase is considered converged when:
1. PSM count ≥ `min_psms` (constant)
2. |fitted bias| < `init_tolerance / 4`
3. Fitted tolerance is within sane bounds

If Phase 1 (wide scout) fails to converge after exhausting its scan budget, the run falls back to:
- A SquareQuad transmission model (already set by `initialize_models!`)
- Conservative defaults (±50 ppm) for the mass error model
and emits a `@user_warn` so downstream methods know calibration was approximate. `get_fallback_parameters` can also borrow from earlier successfully-tuned files in a multi-file run.

## MS1 mass-error fit

After the MS2 fit completes per file, `collect_ms1_residuals` reuses the same MS2-accepted PSMs to fit an MS1 model:

1. Build a `scan_id → nearest MS1 scan` map (mirrors `MainSearch/features.jl:add_ms1_lookup_features!`).
2. For each PSM, look at the precursor's predicted M+0/M+1/M+2 m/z. Search the nearest MS1 spectrum within ±`MS1_DIAG_WINDOW_PPM` (100 ppm) of each target.
3. Record signed `(observed − theoretical)/theoretical × 1e6` ppm residuals (one per matched isotopologue).
4. `fit_ms1_model_from_residuals` returns `MassErrorModel(median, (k·MAD, k·MAD))` with `k = MS1_DIAG_TOL_K_MAD = 5`. MAD is normalized so 1·MAD ≈ σ for Gaussian residuals.
5. Install via `setMs1MassErrorModel!(search_context, ms_file_idx, model)`.

The MainSearch `add_ms1_lookup_features!` and `IntegrateChromatogramsSearch` MS1 path both read this model: target m/z is shifted by the fitted bias before peak matching, and the tolerance window is ±5·MAD.

Diagnostic histogram per file saved to `qc_plots/ms1_mass_error_plots/{file}.pdf` (median = red vline; ±k·MAD = black dashed).

Fallback: if fewer than 10 residuals are collected the fit returns `nothing` and the per-file model dict is left unset; `getMs1MassErrorModel` then returns the default ±30 ppm placeholder and emits a `@user_warn`.

## Outputs

Per file, ParameterTuningSearch:
- Updates `SearchContext.mass_error_model[ms_file_idx]` (plot saved to `qc_plots/mass_err_plots/`)
- Updates `SearchContext.ms1_mass_error_model[ms_file_idx]` (plot saved to `qc_plots/ms1_mass_error_plots/`)
- Updates `SearchContext.irt_rt_map[ms_file_idx]` and `rt_irt_map` (plot saved to `qc_plots/rt_alignment_plots/`)
- Updates `SearchContext.irt_errors[ms_file_idx]`

These models are consumed by every subsequent search method.

## Things that have been removed

The historical doc described:
- 3-phase bias-shift loop (zero / +max_tol / −max_tol) — gone; replaced by single wide-scout phase whose 100 ppm window is wide enough to absorb any bias
- Multi-score iteration (`min_score: [22, 17, 15]`) inside each phase — gone; the scout uses a single threshold
- `mass_tolerance_scale_factor` expansion within iterations — gone; the wide scout's 100 ppm is the only setting
- `IterationSettings` struct + `init_mass_tol_ppm_guess` JSON knob — deleted
- 5 ppm "initial models" — replaced by a 100 ppm placeholder that the wide scout immediately overwrites

If you find references to those in older notes/docs/PRs, they describe a previous design.

## Common pitfalls

- **`getFragTolPpm(::ParameterTuningSearchParameters)`** no longer exists. The fallback was `getInitMassTolPpm(iteration_settings)`; both removed. If you need a placeholder tolerance, use `WIDE_SCOUT_TOL_PPM` directly.
- **The mass error model set by `initialize_models!` is a placeholder** that the wide-scout phase overwrites within the same function. Don't rely on it for downstream code.
- **`MassErrSample` carries no rank/charge/ion-type fields.** Filtering of mass-error samples by rank/charge/y-ion happens **inside** `run_fused_masserr!` before the sample is written; downstream consumers (`fit_mass_err_model`, `fit_intensity_mass_error_model`, `extract_fragment_plot_data`, `generate_wide_scout_plot`) only need theoretical/observed/intensity.
