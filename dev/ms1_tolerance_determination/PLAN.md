# Plan: Online MS1 ppm tolerance determination

## Motivation

`PIONEER_MS1_PPM_TOL` is currently a static knob (default 100 ppm) read directly in `MainSearch/features.jl:367`. The 5-seed SCP_250pg_3ms bench (2026-05-21) showed tightening it to 3 ppm gives +5.4% IDs. Static defaults can't reach this without per-instrument hand-tuning. Solution: fit MS1 ppm online in `ParameterTuningSearch`, the same way MS2 ppm is fit.

## Existing infrastructure (already in place — half-wired)

This is mostly plumbing — the data structures and one consumer already exist:

- `SearchContext.ms1_mass_error_model::Dict{Int64, AbstractMassErrorModel}` — defined in `SearchTypes.jl:258`.
- `getMs1MassErrorModel(s, idx)` / `setMs1MassErrorModel!(s, idx, m)` — defined in `SearchTypes.jl:529-539`. Getter falls back to `MassErrorModel(0.0, (30.0, 30.0))` with a `@user_warn` when missing.
- `qc_plots/ms1_mass_error_plots/` — created in `ParameterTuningSearch.jl:110-111`.
- **`IntegrateChromatogramsSearch/utils.jl:702`** already calls `getMs1MassErrorModel(...)` — currently hits the warn-and-default path on every run because nothing populates the dict.
- `MainSearch/features.jl:367` still reads `PIONEER_MS1_PPM_TOL` env var (default 100) — not yet routed through the model.

So **nothing populates `ms1_mass_error_model`** anywhere. That's the missing piece.

## Implementation

### Step 1 — Collect MS1 residuals during ParameterTuningSearch

In `utils.jl` (next to `fit_mass_err_model`), add `fit_ms1_mass_err_model!()`. Its job per file:

For each top-scoring PSM accepted by the MS2 calibration step (i.e. the PSMs already used to fit `mass_err_model`):

1. Identify the apex MS1 scan(s) for the PSM (RT-bracketed; reuse RT-window logic already in chromatogram extraction).
2. Compute theoretical M+0 m/z = `prec_mz`. Optionally also M+1 = `prec_mz + NEUTRON/charge`, M+2 = `prec_mz + 2·NEUTRON/charge`.
3. Search the MS1 spectrum within ±50 ppm (wide) for each isotopologue.
4. Record `MassErrSample(theoretical_mz, observed_mz, intensity)` per matched isotopologue. This is the same 12-byte struct used for fragments — no new type needed.
5. Require ≥ 2 matched isotopologues to include the PSM (more robust than M+0 alone, especially on SCP).

Feed the resulting `Vector{MassErrSample}` straight into the existing `fit_mass_err_model(...)` — the fit is m/z-agnostic, so the same code works for MS1.

### Step 2 — Wire into the per-file flow

In `ParameterTuningSearch.jl`, after the MS2 mass-error fit completes per file (around `store_final_results!`):

```julia
ms1_model = fit_ms1_mass_err_model!(spectra, search_context, ms_file_idx, accepted_psms)
if ms1_model !== nothing
    setMs1MassErrorModel!(search_context, ms_file_idx, ms1_model)
else
    # Fallback: ±10 ppm offset 0 (tighter than the current 100 default but safe).
    setMs1MassErrorModel!(search_context, ms_file_idx, MassErrorModel(0.0f0, (10.0f0, 10.0f0)))
    @user_warn "MS1 mass error fit failed for $(parsed_fname); using ±10 ppm fallback"
end
```

The fallback should be tighter than today's 100 ppm to actually help — 10 ppm is a defensible Orbitrap/Astral default.

### Step 3 — Route MainSearch through the model

`MainSearch/features.jl:367`:

```julia
ms1_ppm_tol::Float32 = Float32(parse(Float64, get(ENV, "PIONEER_MS1_PPM_TOL", "100.0")))
```

becomes:

```julia
ms1_ppm_tol::Float32 = if haskey(ENV, "PIONEER_MS1_PPM_TOL")
    Float32(parse(Float64, ENV["PIONEER_MS1_PPM_TOL"]))   # explicit override wins
else
    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    max(getLeftTol(mem), getRightTol(mem))                # use fitted half-width
end
```

Env var keeps working as an override for reproducing historical runs / debugging.

### Step 4 — QC

Generate the same Gaussian-fit residual plot as MS2, drop it in `qc_plots/ms1_mass_error_plots/`. The plotting machinery in `fit_intensity_mass_error.jl` / mass-err plot helpers should drop in with the m/z-agnostic `MassErrSample` input.

### Step 5 — Safety clamp

Inside `fit_ms1_mass_err_model!`, clamp the fitted tolerance to **[1.0, 50.0] ppm** before returning. Prevents:
- Pathologically narrow fits on contaminated samples.
- Pathologically wide fits when the PSM pool is thin.

### Step 6 — Diagnostic logging

Log the per-file fitted (μ, σ) to the existing `ParameterTuningDiagnostics` so users see what was chosen vs the static 100 default.

## Validation

3 seeds × 4 datasets (SCP_250pg_3ms, Exploris, MTAC, Astral8) at the new fitted MS1 tolerance vs the static `PIONEER_MS1_PPM_TOL=100` baseline. Targets:

- SCP_250pg_3ms: ≥ 6,033 IDs (the MS1=3 ceiling we got with lr=0.20).
- Exploris, MTAC, Astral8: ≥ static-default IDs (no regression).
- Fitted ppm per dataset should be sensible (Orbitrap ~3-10 ppm; Astral ~2-5 ppm; Exploris ~5-15 ppm).

## Risk

The trickiest piece is reliable monoisotopic peak extraction at low MS1 intensity (SCP). Mitigations:
- ≥ 2 matched isotopologues required (not M+0 alone).
- Pull from many PSMs (the same pool MS2 calibration uses, typically thousands).
- ±10 ppm fallback when the fit fails.

## Effort

~1 day. The fit reuses existing machinery; the per-file pass over MS1 scans is the main new code (~150 lines).

## Files touched

- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl` — new `fit_ms1_mass_err_model!()`.
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl` — invoke + `setMs1MassErrorModel!`.
- `src/Routines/SearchDIA/SearchMethods/MainSearch/features.jl:367` — read from model with env override.
- Possibly `SearchTypes.jl` — update the `getMs1MassErrorModel` fallback warn to be quieter (since the model will be populated by default).

## Open questions

- Per-file fit (matches MS2) or experiment-wide? Recommendation: **per-file** (mirrors MS2; handles instrument drift across runs).
- Include MS1 offset correction? Yes — it's free (the fit already produces μ) and matters on miscalibrated instruments.
- Keep `PIONEER_MS1_PPM_TOL` env var as override? Yes — useful for debugging.
- Should we similarly switch `IntegrateChromatogramsSearch` from `getMs1MassErrorModel` (already wired) to also accept the env override for consistency? Yes — same one-liner.
