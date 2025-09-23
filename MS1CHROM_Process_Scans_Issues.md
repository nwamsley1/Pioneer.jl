# MS1CHROM `process_scans!` Review (SecondPassSearch)

This note tracks identified issues and quick fixes for the MS1 path:

- Wrong mass-error model: `ScoreFragmentMatches!` in the MS1 path passes `getMassErrorModel(...)` (MS2) instead of the MS1 model. Should pass `getMs1MassErrorModel(...)` (or the local `mem`).
- Misordered/duplicated sizing + init: Array sizing and weight initialization previously occurred before building the `id_to_col` map, and targeted MS2 arrays (`getSpectralScores`/`getUnscoredPsms`) in the MS1 path.
- Redundant MS2 arrays in MS1 path: Resize only the MS1 arrays (`getMs1SpectralScores`, `getMs1UnscoredPsms`) for MS1; avoid touching MS2 arrays here.
- Capacity risk for matches: `ion_matches`/`ion_misses` are allocated with a fixed size (10,000). If `matchPeaks!` returns more, risk overflow. Consider dynamic growth like `ion_templates`.
- Dead or unused scratch: `precs_temp`/`prec_temp_size` are maintained but not used downstream. Remove or wire up if intended for later logic.
- Boolean style: Prefer `||` to `|` in `((scan_idx < 1) | (scan_idx > length(spectra)))` for idiomatic short-circuit boolean logic.
- Diagnostics scope: Current grouping stats warn per-thread (scan chunk). If a per-file summary is desired, aggregate across tasks and emit once after `vcat` in the caller.

Status:
- [x] Fix MS1 mass error model usage (use MS1 model in `ScoreFragmentMatches!`).
- [x] Build `id_to_col` from m/z grouping before any use; size MS1 arrays; zero-init group weights.
- [ ] Consider dynamic growth of `ion_matches`/`ion_misses`.
- [ ] Remove or justify `precs_temp` usage.
- [ ] Switch bitwise `|` to logical `||` for bounds check.
- [ ] Consider moving diagnostics aggregation to caller for per-file summary.

