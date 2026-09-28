# Precursor scoring maintenance

See [README.md](README.md) for the training, prediction, and memory lifecycle.

- `PrecursorScoringSearch.jl` orchestrates fold merging, global scoring, calibration,
  filtering, and run similarity.
- `score_psms.jl` provides the unified MBR/non-MBR training entry point and attaches
  aligned predictions during fold merging.
- `pass1_oom.jl` selects the fixed representative training pool and streams predictions.
- `scoring_interface.jl` and `utils.jl` provide file-based calibration helpers.
- Protein inference and protein scoring live in `ProteinScoringSearch`.
- Post-integration donor features and transfer rescoring live in `MBR`.

Preserve CV fold assignments and verify sidecar row identity before attaching
predictions. Keep feature loading, sorting, and prediction workspaces bounded.
Tests cover training pools, fold merging, grouped score calibration, and MBR
sidecar alignment under `test/UnitTests` and `test/Routines/SearchDIA`.
