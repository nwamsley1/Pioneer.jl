# Koina Empirical Findings (verified 2026-07-02)

Server: https://koina.wilhelmlab.org  — reachable from this env via curl/urllib. `/v2/health/ready` = 200.
Triton inference server. Endpoints:
- `GET /v2/models/{model}/config` — I/O contract
- `POST /v2/repository/index` — full model list
- `POST /v2/models/{model}/infer` — inference

## Relevant models (verified present)
### Fragment intensity
- `Prosit_2019_intensity` — classic, no PTM
- `Prosit_2020_intensity_HCD` / `_CID` / `_TMT`
- `Prosit_2024_intensity_PTMs_gl` — **PTM-capable** (general)
- `Prosit_2024_intensity_cit` — citrullination
- `Prosit_2025_intensity_22PTM`, `_40PTM`, `_ptm2` — newer PTM sets
- `Altimeter_2024_*` (current Pioneer model; spline-based)
- AlphaPeptDeep_ms2_generic, ms2pip variants, UniSpec

### Retention time / iRT
- `Chronologer_RT` (current Pioneer RT model)
- `Prosit_2019_irt`, `Prosit_2020_irt_TMT`
- `Prosit_2024_irt_PTMs_gl` — **matched iRT for the PTM intensity model**
- `Prosit_2025_irt_22PTM`, `_40PTM`, `_ptm2`

## I/O contracts (verified)
`Prosit_2020_intensity_HCD` inputs:
  peptide_sequences STRING[-1], precursor_charges INT32[1], collision_energies FP32[1]
  outputs: intensities FP32[174], mz FP32[174], annotation STRING[174]   (fixed 174 slots)

`Prosit_2024_intensity_PTMs_gl` inputs (NOTE extra input):
  peptide_sequences, precursor_charges, collision_energies, **fragmentation_types STRING[1]** ("HCD"/"CID")
  outputs: intensities/mz/annotation FP32/STRING[174]

`Prosit_2024_irt_PTMs_gl` inputs: peptide_sequences only -> outputs: irt FP32[1]
`Chronologer_RT` inputs: peptide_sequences only -> outputs: rt FP32[1]

max_batch_size = 1000 for Prosit models.

## Verified inference results
Phosphopeptide round-trip on `Prosit_2024_intensity_PTMs_gl` (charge 2, CE 27, HCD):
- PEPTIDEK: 28 nonzero peaks; y6+1 dominant
- S[UNIMOD:21]PEPTIDEK: 32 peaks; b2 shifts 227.103 -> 265.058 (+79.966 phospho on S1) — CORRECT
- ELVIS[UNIMOD:21]LIVESK: 40 peaks
Matched iRT (`Prosit_2024_irt_PTMs_gl`): [4.06, 17.59, 108.12] for the three.

**Chronologer DOES support PTMs** (verified accepts + shifts RT):
- PEPTIDEK -> 4.20; S[UNIMOD:21]PEPTIDEK -> 5.56 (phospho); M[UNIMOD:35] -> 4.63 (ox); C[UNIMOD:4] -> 4.72 (carbamidomethyl)

## Modification notation
ProForma / UNIMOD inline: `S[UNIMOD:21]` (phospho), `[UNIMOD:35]` ox-M, `[UNIMOD:4]` carbamidomethyl-C.

## Key design implications
1. Prosit gives RAW intensities normalized to base peak = 1.0 (no NCE calibration; CE is a raw input, not calibrated per-run like Altimeter).
2. Prosit predicts a FIXED 174-slot vector (b/y ions, charges 1-3, up to ~30 residues) — must filter zero/negative (-1 = not applicable) entries.
3. PTM intensity model needs `fragmentation_types` — classic 2020 models do not.
4. Isotopes: Prosit predicts monoisotopic fragment intensity only (not total abundance) — Pioneer must add isotope envelope itself (contrast Altimeter which handles isotopes/reisotope internally via Altimeter_2024_isotopes/reisotope).
5. RT for PTMs: TWO options — reuse Chronologer (already integrated, supports phospho) OR Prosit matched iRT. Recommend Chronologer for parity in unmodified/phospho case first.
