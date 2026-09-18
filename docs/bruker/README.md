# Bruker timsTOF `.d` → Pioneer Arrow (prototype converter)

Pure-Julia reader and converters for timsTOF diaPASEF bundles (no Bruker SDK). See
`docs/bruker_timstof_progress.md` for the format, the decisions, and the results behind these.

Environment: `Project.toml` here (Arrow, CodecZstd, DataFrames, SQLite). `julia --project=docs/bruker`.

- `tdf_reader.jl` — SQLite + zstd block decoder (bit-identical to `opentimstdf`), linear m/z calibration from
  `CalibrationInfo`, IM calibration.
- `tdf_to_arrow.jl` — packet conversion (one row per frame × IM scan inside a window); options `--im-bin N`,
  `--merge-tol K`, `--zstd-only`, `--no-zstd`. Includes `tdf_reader.jl`.
- `tdf_centroid_to_arrow.jl` — IM-smoothed, m/z-centroided slices (the format the searches use). Includes
  `tdf_to_arrow.jl`. Options: `--im-sigma 5 --mz-sigma 1 --stride 8 --cull-q 0 --centroid wmean --sum-scale
  [--ms1-stride K] [--ms1-cull-q Q] [--split-cull] [--min-scans N] [--zstd]`. Output name encodes the settings.
- `convert_many.jl` — several conversions of one `.d` in a single Julia session:
  `julia -t 12 --project=docs/bruker docs/bruker/convert_many.jl <run.d> <out_dir> "5 1 8 0 wmean sum" "5 1 8 0.05 wmean sum ms1q=0"`.
- `test_centroid_equiv.jl` — checks the sparse per-slice accumulation against a dense-matrix reference.
- `count_passing.jl` — precursors at 1% / 0.1% FDR from a run's `temp_data/passing_psms`
  (`julia --project=<Pioneer> count_passing.jl <run_dir>...`).

Best settings so far: `--im-sigma 5 --mz-sigma 1 --stride 8 --centroid wmean --sum-scale`, with `--cull-q 0.05`
at 50 ng and no cull at 250 pg. Files are large (2.6× the raw peak count at 250 pg); see the progress document §7.
