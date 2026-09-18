# Bruker timsTOF diaPASEF support — progress and decisions (branch `feat/bruker-timstof`)

Written 2026-09-18 after four days of work (2026-09-15 to 09-18). This is the memory of the branch: what the
raw data looks like, what was built, every decision with its measurement, and what is open. The working notes
that this distils live outside the repo in `~/BrukerTims/NOTES.md` (Nathan's machine, ~1,000 lines).

## 1. Goal and test data

Read Bruker timsTOF diaPASEF `.d` bundles natively, convert them to Pioneer's Arrow schema, and search them with
the ion-mobility dimension used rather than collapsed. Test files (PRIDE PXD070049, the LFQ benchmark
"Generation Beta"; plus a small HeLa file from PXD027359):

| file | instrument | gradient | load | raw `.d` | raw peaks MS1 / MS2 | use |
|---|---|---|---|---|---|---|
| `LFQ_Ultra2_diaPASEF_5min_50ng_Ecoli_01` | Ultra 2 | 5 min | 50 ng | 2.7 GB | 587 M / 503 M | development file |
| `LFQ_Ultra_diaPASEF_15min_250pg_Human_01` | Ultra | 15 min | 250 pg | 2.5 GB | 722 M / 209 M | single-cell-scale file |
| `LFQ_Ultra2_diaPASEF_15min_50ng_Human_01` | Ultra 2 | 15 min | 50 ng | 6.9 GB | 1.18 G / 1.60 G | typical input (first converted 09-17) |
| `..._HeLa_50ng_5_6min_DIA_..._25186` | Pro | 5.6 min | 50 ng | 0.35 GB | | tiny regression file |

All are stepped diaPASEF (not synchro-PASEF): 8 window groups × 3 fixed 25 Da isolation windows, each parked over
a distinct IM scan range, 24 windows with 24 distinct m/z centres over 400–1000 m/z, ~940 IM scans per frame,
0.97 s cycle. MS1 frames are unfiltered (quad pass-through, 100–1700 m/z). Reference: DIA-NN 2.6.0 library-free
single-file: E. coli 14,849 precursors, human 250 pg 9,932.

## 2. The raw format (validated reader, bit-identical to `opentimstdf`)

`analysis.tdf` is SQLite (`Frames`, `DiaFrameMsMsInfo`, `DiaFrameMsMsWindows`, `CalibrationInfo`, `MzCalibration`,
`TimsCalibration`, `GlobalMetadata`). `analysis.tdf_bin` holds one block per frame at `Frames.TimsId`:

- `u32 block_size`, `u32 scan_count`, then zstd-compressed payload ("TimsCompressionType" 2).
- Decompressed length `4 · (NumScans + 2 · NumPeaks)`, **byte-transposed**: four byte-planes of that length;
  `u32[i] = b0[i] | b1[i] << 8 | b2[i] << 16 | b3[i] << 24`.
- `u32[0] = NumScans`; `u32[s+1] / 2` = peak count of scan `s` (last scan takes the remainder); then an
  interleaved `(tof_delta, intensity)` stream. **TOF is delta-coded per scan** (accumulator starts at
  0xFFFFFFFF, so delta 1 = bin 0); intensity is a plain u32.
- Net: **2.7 bytes per peak** on disk. Each scan is a list of integer TOF-bin centroids (one entry per ion per
  scan, entries ≥ 5 bins apart); the same ion's bin jitters ±1–2 between adjacent IM scans.
- m/z: `√(m/z) = a + b · bin`, fit by least squares on `CalibrationInfo` (`√ReferencePeakMasses` vs
  `(MeasuredTimesOfFlight − DigitizerDelay)/DigitizerTimebase`); exact (0.000 ppm) on the Ultra 2 file's 3
  reference peaks, ≤ 2 ppm on older timsControl files. Do NOT use the boundary model from
  `MzAcqRange*/DigitizerNumSamples` (−20…+7 ppm). Bin width in ppm = 7.0 at m/z 200, 4.4 at 500, 3.1 at 1000
  (631,387 bins over 100–1700 on the Ultra; 638,506 on the Ultra 2).
- 1/K0: linear in scan from `OneOverK0AcqRange*` and `NumScans`, ~4% slope error on the HeLa file; the search
  fits its own per-file scan-vs-library-1/K0 line instead (§5).
- Collision energy is an eV ramp in scan (window table value = ramp at the window's mid-scan).

Decode speed of the pure-Julia reader (SQLite.jl + CodecZstd.jl): ~35–43 M peaks/s single-threaded. Bruker's
SDK is Windows/Linux only; not needed.

## 3. Conversion: packets → IM-smoothed, m/z-centroided slices

Decided 2026-09-15: Arrow rows = raw packets (one row per frame × IM scan inside a quad window), Pioneer's
existing schema plus `frameId`, `imScan`, `windowGroup`, `cycle_idx`, `collisionEnergyEvField`. Packet files:
E. coli 8.5 GB, human 250 pg 8.2 GB. Searchable, but tuning never converged on packets (one precursor = ~7
adjacent-scan PSMs) and adjacent-scan centroid jitter split ions on exact-bin merging.

**Centroiding converter** (`docs/bruker/tdf_centroid_to_arrow.jl`, driver `convert_many.jl`): per frame and DIA
window (MS1: whole frame), for slices every `stride` scans: (1) IM Gaussian (σ scans) accumulated per TOF bin,
(2) m/z Gaussian (σ bins) over dense runs of nearby bins, (3) local maxima with a footprint walk (≤ 4 bins),
m/z = intensity-weighted mean bin (sub-bin), intensity = footprint sum, (4) optional culls: raw-intensity
quantile (`--cull-q`, per level with `--ms1-cull-q`/`--split-cull`), persistence (`--min-scans`). `--sum-scale`
scales the IM kernel to the stride so a centroid carries the ion's summed intensity over the slice's scans
(without it a raw-quantile cull is ~8× too strict). Multi-record-batch output above 1 G peaks (Arrow list
offsets are Int32). Verified equal to a dense-matrix implementation.

Sweep results (E. coli, 1% FDR): m/z kernel essential (−24% without), σ = 5 scans optimum, stride 8 the knee,
weighted mean > Gaussian apex. **Best settings**: σ5, m/z σ1, stride 8, wmean, sum-scaled; cull q0.05 at 50 ng
(free), **no cull at 250 pg** (any cull costs IDs: singletons are the signal there).

## 4. Search results (single file, no MBR, targets at 1% / 0.1% FDR)

| stage | E. coli 50 ng | human 250 pg |
|---|---|---|
| packets, first search | 9.1 k | |
| + im_error feature, square quad | 12,112 | 2,834 |
| centroided slices (best cull) | 14,026 / — | 7,610 / — |
| + MS1 lookup at the precursor's mobility | 14,622 / 12,675 | 9,349 / 5,700 (tuning window ×2) |
| + default tuning window (×1) | 14,751 / 12,979 | 9,670 / 6,164 |
| + n_scans_in_window + weight_frac_in_cycle | ~14,800 | 9,914 / 5,274 |
| + new tuning (this branch head) | 14,956 / 12,886 (MS1 uncalled) | 9,646 / 6,270 |
| DIA-NN 2.6.0 | 14,849 | 9,932 |

Runs are deterministic to the digit (checked on single-file Bruker and 3-file SCP Astral runs).

## 5. What changed in Pioneer (by commit)

- **Library** (`1171f0e5d`): `library_params.im_model` (alphapept_ccs | im2deep via Koina) → `ccs`,
  `inv_ion_mobility` precursor columns; `prec_partition_width`. Published libraries have no mobility columns;
  `add_im_to_lib.jl` (workspace) annotates one in place.
- **Dev hooks** (`23778684a`, `dd37b4e3a`): `PIONEER_STOP_AFTER`, `PIONEER_TUNING_DUMP_PSMS`, `PIONEER_INDEX_DUMP_DIR`.
- **Tuning on packets** (`dc79f0955`): stratified (RT × window × IM-bin) scan priority order for files with an
  `imScan` column (`get_ms2_scan_priority_order_im`); unique-precursor convergence unit (superseded, see below).
- **NCE keyed on the scan's collision energy** (`2f62b7193`): `CeBinnedNceModel`; the eV ramp gives per-scan NCE.
- **IM calibration + `im_error`** (`4cff15ab3`, `0c3f6f379`, `6f1440f10`): `add_im_error!` fits per-charge lines of
  library 1/K0 vs IM scan on prob > 0.9 targets; `im_error` = |residual|/σ is a scoring feature (rank ~18/72);
  models stored on the SearchContext (`getImModel`); QC plots at `qc_plots/ion_mobility_model/`. Fitted σ ≈
  0.014 1/K0 ≈ 15 packets for z=2 (library-prediction error, not the ~14-packet instrument peak width).
- **IntegrateChromatograms on packets** (`5f055d7f9`, `532233bfa`): IM gate (±3σ of the fitted line,
  `PIONEER_CHROM_IM_SIGMA`), `PIONEER_CHROM_RT_TOL`, fitted-weight dump (`PIONEER_CHROM_DUMP_DIR`, needs
  `match_between_runs=false`), `PIONEER_CHROM_ALT_FILES`. 2D chromatograms are clean blobs ~14 scans × ~2 cycles.
  **1D integration is not valid on slices** (WH smoothing divides by RT spacing); a 2D design is still open.
- **Quad model** (`c3896a0e6`, `a92265804`): on files with `imScan`, QuadTuning uses `SquareQuadModel` at the
  reported 25 Da window (single-packet isotope ratios are noise). +2.6%.
- **MS1 lookup at the precursor's mobility** (`1232e783b`, `build_scan_to_ms1`): nearest MS1 frame by RT, then
  the MS1 row at the nearest `imScan`. Before: MS1 features populated for 2.8% of PSMs. +20% human, +3.5% E. coli.
  Explainer with quoted code: `~/BrukerTims/MS1_features_explained.{md,pdf}`.
- **Mobility features** (`fa79b3554`, `22ddfd176`, `fac0643fd`): `n_scans_in_window` (the precursor's PSMs in the
  same cycle = slices at that RT) and `weight_frac_in_cycle` (weight / sum of the precursor's weights in the
  cycle). Together +2.5% at 1% FDR on human 250 pg; each alone less. Tried and dropped: mobility apex offsets
  (`831d4fc4f`, reverted `f222a2ba7`), per-slice and per-cycle fragment-correlation effective-n, weight/max ratio.
- **Parameter tuning** (`231ee3b8d`): second-order stopping (back off a score tier when the decaying marginal
  yield cannot reach the target), every tier grows from `initial_scans`, scout counts PSM rows, collection counts
  rows after a best-3-per-precursor cap. Human tuning 182 s → 24 s at −0.8% IDs; SCP Astral 250 pg conditions
  unchanged or +1%, tuning 2–3× faster. Rejected on the way: tuning-only slice stride (−7%), fraction cap (Nathan),
  scout min-batches (a batch can be one scan), never-back-off-on-last-tier (−5%). Per-checkpoint and per-phase
  summary logs. Flags read at call time: `PIONEER_TUNING_SECOND_ORDER`, `PIONEER_TUNING_MAX_PER_PREC`.
- **Window multiplier**: `PIONEER_TUNING_IRT_TOL_MULT` (`42a21b26c`) exists; ×2 was best before the MS1 fix,
  ×1 (default) is best after it (+5.7% human). Leave unset.

## 6. Things learned the hard way

- **Never gate feature lists or behaviour with a load-time `const` read from ENV.** It is evaluated at
  precompile and cached; every later process inherits whichever value compiled it. Two experiments (feature
  subsets, tuning flag) were silently invalid this way. All dev flags now read ENV inside the function.
- `LibrarySearch.jl` unwraps `TaskFailedException` and rethrows, which discards the worker's backtrace; a
  bounds error inside `process_scans_fused!` shows only the `fetch` frame. Print `e.task.backtrace` when hunting.
- **UInt16 column limit** (open bug): the fused design matrix's column counter and `AbstractPrecursorMap{UInt16}`
  cap a scan at 65,535 matched precursors. The 50 ng human slice file (MS1 slices of 34 K peaks at the scout's
  ±100 ppm) exceeds it → `BoundsError ... 65536-element Vector{UInt32} at index [0]` in `finalize_column!`.
  Fix: widen to UInt32 through `SparseArrayFused`, `PrecursorMap`, `run_fused!`, `finalize_column!`,
  `deconvolutionArrayUtils.jl`; or cap matched precursors per scan in the scout.
- The scan priority order front-loads richness (TIC top of every RT bin first), so the marginal yield of batch 2
  is always far below batch 1 on every file; a decay rule must not read that first drop as the file's decay.
- The per-precursor 3-PSM cap + unique-precursor counting (`dc79f0955`) applied to all file types and roughly
  tripled the effective tuning targets on Astral; the per-phase units in `231ee3b8d` replace it.

## 7. The open problem: file size

The IM smoothing that makes the search work multiplies the data. Measured (σ5, stride 8):

| file | raw `.d` | slice Arrow | peaks vs raw | notes |
|---|---|---|---|---|
| human 250 pg | 2.5 GB | 19.3 GB | 2.6× | no cull (any cull costs IDs) |
| E. coli 50 ng | 2.7 GB | 5.7 GB | 0.6× | MS2+MS1 q0.05 cull |
| human 50 ng | 6.9 GB | 25.2 GB | 1.1× | MS2 q0.05, MS1 uncalled |

Causes: 8.25 B/peak (Float32 m/z + Float32 intensity + validity bitmaps) vs Bruker's 2.7; and every MS1 ion
yields one centroid in every slice it persists in (the σ5 kernel reaches ±15 scans; even σ2 leaves 1.3× raw).
**Narrowing the kernel does not help** (E1, 2026-09-17): IDs fall in proportion to bytes — human 250 pg
σ3 −7.6%, σ2.5 −10.7%, σ2 −15.2%; E. coli −2.7 / −4.3 / −5.7%. The kernel is integrating each ion over its
mobility peak; that integration is the signal at low input.

Plan (`~/BrukerTims/PLAN_2026-09-17_file_size.md`): (a) memory-mapped fixed-width integer encoding — UInt32
fixed-point TOF bin (bin·256, keeps the sub-bin centroid) + UInt16 √-encoded intensity (6 B/peak), then UInt16
per-scan deltas with an escape (4 B/peak, decoded into scratch at the scan visit); (b) a different treatment of
MS1 (stride, or smoothing at load time from compact packets); (c) load-aware per-level cull. zstd Arrow is
rejected: Arrow.jl decompresses whole buffers into RAM at open, so it trades disk for RSS.

## 8. Workspace (outside the repo, Nathan's machine)

`~/BrukerTims/`: `pride/` (bundles), `arrow/` (converted files), `proto/` (converter, reader, ~40 analysis
scripts; env `proto/Project.toml`), `search/` (configs, `run_many.jl` one-session driver, batch scripts,
`out_*` results), `lib/` (human library with mobility added), `scp/` (SCP Astral conditions copied from RIS for the
tuning validation), `NOTES.md`, `REPORT_2026-09-16_bruker_centroiding.md`, `MS1_features_explained.{md,pdf}`,
`PLAN_2026-09-17_file_size.md`. Copies of the reports on RIS `NTW/`. DIA-NN reference runs on RIS
`NTW/BrukerTims_DIANN/`. Best Arrow files to keep: human 250 pg `*_cen_s5_m1_k8_q0_wmean_sum.arrow`,
E. coli `*_cen_s5_m1_k8_q0.05_wmean_sum{,_ms1q0}.arrow`, human 50 ng `*_cen_s5_m1_k8_q0.05_wmean_sum_ms1q0.arrow`.
Launch long jobs with `nohup` (no `setsid` on macOS); run several searches in one Julia session (JIT ≈ 110 s each).
