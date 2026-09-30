# The `.scxs` container

Centroided SCIEX SWATH scans, one zstd block per scan. It's the `.tdfs` layout of TimsSlices.jl with a
per-scan m/z calibration and no ion mobility. `SciexWiff.convert` writes it; Pioneer reads it with
`ScxsMassSpecData` (branch `feat/sciex-scxs`).

```
<name>.scxs/
  meta.json    format_version (1), source files, windows, scan ranges, centroid parameters,
               bin_scale, int_scale, zstd_level, converter, n_scans, n_peaks, blocks_bytes
  scans.arrow  one row per scan (= one Pioneer scan) in acquisition order; empty Idx records are
               not written, so rows line up one to one with msConvert's spectra
  blocks.bin   one zstd block per scan, in scans.arrow order, no headers
```

## blocks.bin

The format is **byte-identical to TimsSlices' per-slice `.tdfs` blocks**, so
`TimsSlices.decode_slice!` decodes them. A block is one zstd frame that decompresses to
`8 · n_peaks` bytes. That's four byte planes of `2 · n_peaks` bytes each (plane p holds byte p
of every little-endian UInt32 word). Reassembled, the word stream holds, per peak:

```
bin_delta   stored bin minus the previous peak's stored bin; the accumulator starts at 0xFFFFFFFF,
            so the first delta is bin + 1
intensity   UInt32
```

Peaks are sorted by stored bin and unique; centroids that round to the same bin are merged by
summing. An empty scan has `block_size` 0 and `n_peaks` 0.

- stored bin `q = round(bin_scale · position)`, where position is the centroid in raw TOF units
  (one decoded bin step is 8 raw units)
- **m/z = (cal_a · (q / bin_scale / 40 − cal_b))²**, with `cal_a`, `cal_b` taken from the scan's
  row. SCIEX recalibrates during the run, so this is per scan.
- intensity = stored / `int_scale`. Before scaling it is the centroid area (100 · Σ I·Δm/z), the
  scale the vendor centroider reports.

The defaults are `bin_scale = 4` (0.31 ppm steps at m/z 400) and `int_scale = 10`. Coarser grids
barely shrink the file (paper run: 3.8 B/peak at 1 / 1 against 4.5 at 4 / 10).

## scans.arrow

| column | type | meaning |
|---|---|---|
| record | Int32 | 1-based Idx record in the `.wiff` |
| cycle | Int32 | acquisition cycle (1-based); Pioneer's `cycle_idx` |
| experiment | Int16 | experiment in the cycle (1 = MS1, k+1 = SWATH window k) |
| ms_order | UInt8 | 1 or 2 |
| retention_time | Float32 | minutes |
| low_mz, high_mz | Float32 | the experiment's acquired m/z range |
| tic | Float32 | the instrument's TIC (Idx), identical to msConvert's |
| center_mz, isolation_width | Float32 | the SWATH window (NaN for MS1) |
| cal_a, cal_b | Float64 | the scan's m/z calibration |
| n_peaks | Int32 | centroids in the scan |
| block_offset, block_size | Int64, Int32 | the scan's block in blocks.bin |

## The Arrow output

`convert(...; params = ConvertParams(format = :both))` also writes Pioneer's standard scan Arrow
(the `convertMzML` columns plus `cycle_idx`). It holds the same quantised peaks as Float32 values,
so a search from either output sees identical spectra. This is checked in two places:
`scripts/validate_outputs.jl`, and on the Pioneer side in
`test/UnitTests/test_scxs_mass_spec_data.jl`.
