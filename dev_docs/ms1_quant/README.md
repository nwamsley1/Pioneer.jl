# MS1 quantification (DDA and DIA)

## Acquisition mode
Each file's acquisition mode is decided at search start (`src/Routines/SearchDIA/acquisition_mode.jl`): the
converter's `acquisition_mode` label when present, otherwise a rule on the data (coverage of the MS2 isolation windows
per MS1 cycle: DIA tiles the precursor range every cycle, DDA samples it). `search.acquisition_mode` (`auto` | `dda` | `dia`)
overrides it. A search that mixes DDA and DIA files is refused.

## DDA
A DDA precursor is fragmented in about one MS2 scan, so it is quantified from MS1
(`IntegrateChromatogramsSearch/ms1_quant.jl`):
- per-file intensity-aware MS1 mass-error model fitted on the most confident precursors;
- per MS1 scan, the candidates' M0-M2 isotope envelopes are fitted to the peaks (Huber with a small ridge);
- per-file extraction window and smoothing scale calibrated from the MS1 peaks around the identifying scans, then a
  second extraction re-centred on each precursor's MS1 apex (`PIONEER_MS1_TWO_PASS=0` disables the second pass);
- RT x ion-mobility integration on timsTOF data.
MBR and PSM scoring leave out the features that assume an MS2 elution profile.

## DIA (optional)
`maxLFQ.ms1_quantification = true` adds MS1 quantification next to MS2 for DIA files. After MS2 integration, each
precursor's MS1 trace is extracted and integrated exactly like its MS2 trace: candidacy centred on the MS2 apex, the
MS2 chromatogram window, a narrow mobility band, window-relative smoothing. MS2 results are unchanged. Outputs:
`peak_area_ms1`, `peak_area_ms1_normalized`, `quant_withheld_ms1` (precursors); `abundance_ms1`,
`total_peak_area_ms1`, `n_precursors_quantified_ms1` (protein groups); `precursors_wide_ms1` and
`protein_groups_wide_ms1`; MS1 counts in `run_summary.tsv`.

On 250 pg three-proteome data MS1 is as accurate as MS2 on the Astral with about half the CV, and about equal on
timsTOF. It is not reliable at high load (Astral 50 ng): low-abundance precursors and crowded MS1 scans compress the
ratios. Off by default.

## Segmented and BoxCar MS1
When a file's MS1 scans measure different m/z ranges, a precursor's MS1 trace only uses the scans that measure its
M0 or M1 (`CommonSearchUtils/ms1_coverage.jl`, `Ms1Coverage`), and the per-PSM MS1 features read the nearest MS1 scan
that measures the isolation centre. A scan measures [lowMz, highMz], or for a multiplexed (BoxCar / msx) scan the
boxes listed in its filter string (`Full msx ms [500.0000-510.0000, 530.0000-540.0000, ...]`). Files whose MS1 scans
all measure the same range are unchanged. On simulated segmented Astral data (MS1 split in two) the MS1 areas equal
those of the original data; on BoxCar (PXD006109 HeLa/E. coli) the E. coli ratio is 2.64 log2 (expected 2.58), with
lower CV than the standard runs of the same samples.

## Known limitations
- convertRaw does not yet write the boxes of a BoxCar scan into its scan header. Until it does, each BoxCar scan
  counts as measuring its whole range and BoxCar MS1 areas are wrong.
- The 2D (timsTOF) integrator has no withheld-area rule, so `quant_withheld_ms1` is always false there.
