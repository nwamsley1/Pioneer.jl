# `.wiff` + `.wiff.scan` format: what we have verified

Two ZenoTOF 7600 SWATH runs so far:

- **paper**: JPST002949 `BenchSample_B_nswath4_25ng` (173 windows, 917 cycles).
  Oracle: the msConvert 3.0.21229 vendor-centroided mzML.
- **pride**: PXD036786 `..._K562_0.98ng_1` (60 windows, 1699 cycles). No oracle yet, so
  only oracle-free checks against the Idx.

Status labels: **VERIFIED** means checked on every scan, against the oracle where one exists.
**INTERNAL** means checked on every scan, but only against other fields in the file.
**OPEN** means not resolved.

Scripts in `scripts/` reproduce every number here (`validate_decode.jl`, `internal_check.jl`,
`details.jl`, `raw_vs_centroid.jl`).

## 1. `.wiff` container

A CFB v3 container with 512-byte sectors (`src/CFB.jl`). The paper file has 3,528 streams
and the PRIDE file 1,329. The parts that matter:

| Stream | Content |
|---|---|
| `SampleSubtree/Sample1/Idx` | scan index (§2) |
| `MethodSubtree/Method1/DeviceMethod0/SWATHMethod` | 40-byte preamble + 20-byte `{f64 lo, f64 hi, u32}` records, one per window (paper 173, PRIDE 60) |
| `.../Period0/Experiment{N}/...` | one Experiment per cycle slot: Experiment0 is MS1, then one per window (174 and 61) |
| `SampleSubtree/Sample1/TOFCalibrationData`, `TDCStatistics`, `Itc` | present, but **not needed**; calibration is carried per block in `.wiff.scan` (§4) |

Both files hold a single sample (`Sample1`). The paper file's mzML names show that some
`.wiff` files from this study hold several samples.

## 2. `Idx`: VERIFIED against the oracle on the paper run

A 32-byte header (`00 00 00 00 04 00 00 00` + 24 zero bytes), then **54-byte records, one per
(cycle, experiment) slot, in acquisition order**: record `r` is cycle `r ÷ E + 1`, experiment
`r % E + 1`, where `E` is the number of SWATH windows plus 1. Both files divide exactly.

| Offset | Type | Meaning | Evidence |
|---|---|---|---|
| 0x00 | u32 | byte offset of the block's `ffffffff`, **relative to byte 44** of `.wiff.scan` | every block found |
| 0x04 | u32 | block size in bytes (0 means an empty scan) | blocks tile `.wiff.scan` exactly to EOF |

The u32 offsets wrap once a `.wiff.scan` passes 4 GiB (seen on a 7.65 GB ZT run). Blocks are contiguous and in
order, so `ScanIndex` unwraps them: each block starts after the previous one ends.
| 0x08 | f64 | **retention time in milliseconds** | = mzML RT × 60000, max error 3e-13 min, 151,423/151,423 |
| 0x10 | u16 | `6` in every record (**not** the MS level, contrary to OpenSXRaw) | OPEN |
| 0x12 | f64 | TIC | = mzML TIC, max relative error 5e-13 |
| 0x1A | f64 | base-peak intensity | = mzML BPI exactly; = max decoded intensity, 100% of scans in both files |
| 0x22 | f64 | unknown; equals the minimum decoded intensity in 59% of scans | OPEN |
| 0x2A | f64 | base-peak **bin** (§4) | = bin of the decoded maximum in 99.98% (paper) and 99.997% (PRIDE) |
| 0x32 | u32 | always 0 | |

- The MS level comes from the experiment index: experiment 1 is MS1, the rest are MS2 (151,423/151,423).
- **Empty scans**: 8,135 paper records have size 0 and TIC 0. msConvert omits exactly these
  records and no others.
- The 23 base-peak bin mismatches are all very intense peaks (74k–157k counts), off by 1–9
  bins. OPEN, but they don't affect decoding.

## 3. `.wiff.scan` layout: INTERNAL on both files

```
[44-byte file header]
repeat per non-empty Idx record, in Idx order:
  [varint len][protobuf metadata, len bytes]
  [ff ff ff ff][u32 start bin][00][token stream][0xff × 1..4 padding]
```

- **File header (44 B)**: `u32 2200`, `u64 file size`, 8 zero bytes, `11 11 11 11`, `u32 2200`,
  `u32 1`, `u64 (file size − 44)`, 4 zero bytes.
- **Blocks are contiguous**: each block's metadata starts where the previous block ended.
  Idx `off + 44` points at the `ffffffff`, and `off + 44 + size` is where the next block's
  metadata starts.
- **Padding**: `0xff` bytes pad the block to a multiple of 4. The paper run has 1, 2, 3 and 4
  padding bytes in roughly equal numbers. Intensity bytes are never `0xff` (§4), so trailing
  `0xff` bytes can be stripped without ambiguity.
- **Metadata protobuf**, whose varint length prefix can be longer than one byte:
  - field 1 (len-delimited) holds `{1: f64 a, 2: f64 b}`, the per-block calibration.
  - field 2 (repeated, len-delimited) holds `{1: varint start_bin, 2: varint x}`. The first entry
    is always `(start bin, 8)`, and its start equals the `u32` after `ffffffff`. Later entries
    appear about every 512 bytes of token data, with `x = 512k + (8…13)`. The start bin is an
    actual peak's bin (99.9%), and that peak's token begins 0–3 bytes before offset `x`. This
    looks like a **seek index** for random access by m/z. It is **not needed for sequential
    decoding**; the exact meaning of `x` is OPEN.
- The paper run has 1–219 segments per block (median about 4).

## 4. Token stream and calibration

**Grammar.** It decodes all 151,423 paper blocks and all 95,552 PRIDE blocks, with no leftover
bytes. Each token is a delta prefix followed by an intensity field. Let `m` be the running sum
of the deltas, **including the first token's delta**:

| Delta prefix | Δ |
|---|---|
| `00–7b` (no prefix byte; this byte starts the intensity) | 1 |
| `80–fb` | `byte − 0x7f` (1–124) |
| `fc v` | `v + 1` |
| `fd lo hi` | `u16 + 1` |
| `fe b0 b1 b2 b3` | `u32 + 1` (**4 bytes**; sciexwiff treats `fe` as Δ = 127, which is wrong here) |

| Intensity field | Value |
|---|---|
| `00–7b` | literal |
| `7c v` | u8 |
| `7d lo hi` | u16 |
| `7e b0 b1 b2 b3` | u32 (**4 bytes**; see below) |

**Correction (2026-09-25):** `7e` carries a 4-byte intensity. Reading it as 3 bytes (as first written here)
leaves its high byte, 0x00, to be parsed as a spurious "next bin, intensity 0" token, shifting every later peak of
the block by one step (~10 ppm). In SWATH this only affects intensities > 65,535 (the 23 base-peak mismatches
below were exactly these); in ZT MS2 (100 per ion) it hits every peak above ~655 ions. With the fix, the decoded
base-peak bin equals the Idx base-peak bin in 151,423/151,423 (paper) and 95,552/95,552 (PRIDE) blocks, and no
zero-intensity points remain in any file.

Neither file has intensity bytes ≥ `0x7f` and neither has an `ff` delta prefix. The raw
128–255 intensity form in sciexwiff never occurs here.

**Bin position**: `bin = start_bin + 8·m`. This matches Idx `0x2A` exactly, as noted above.

**Calibration: VERIFIED against the oracle on the paper run.**

```
m/z = ( a · (bin / 40 − b) )²
```

`a` and `b` come from that block's metadata. Tested against the mzML base-peak m/z of every
scan: **median error 0.00 ppm, maximum 1.52 ppm, n = 151,423**. Equivalently, `t = bin/40` is
a flight time, `b` is t₀ and `a` is the √(m/z) slope. One stored bin step (8 units) is about
10 ppm at m/z 400 and about 6 ppm at m/z 1000. So the vendor's base-peak m/z is reported at
the bin itself, not interpolated. Neither reference reader's formula works here: OpenSXRaw's linear
`slope·bin + intercept` and sciexwiff's `(a/5·n + b)²` both fail. The paper run has only 4
distinct `(a, b)` pairs, so the calibration is effectively constant. The PRIDE run has 3,078,
so it is updated every few scans.

## 5. Isolation windows: VERIFIED against the oracle on the paper run

`SWATHMethod`: a 40-byte preamble whose `u16` at byte 38 is the window count (`0xad` = 173),
then 20-byte records `{f64 lo, f64 hi, u32 = 0}`. **Experiment k+1 is window k**. All 173
match the mzML exactly: `[target − lower offset, target + upper offset] = [lo, hi]`, so
Pioneer's `centerMz = (lo+hi)/2` and `isolationWidthMz = hi − lo`. On this method the windows
are 2.9 Da wide with no overlap, from 400.0 to 901.7.

## 6. What the raw data look like

- Profile but **sparse**: zero bins are dropped. In the paper run the median is 776 stored
  bins per scan and the maximum 90,594. 24.6% of neighbouring stored bins are adjacent
  (Δ = 1), and 3.8% are gaps of ≥ 1000.
- Isolated noise events have intensity about 40–44, so intensity isn't simply 1 per ion.
- One MS1 scan (record 52200) compared with the vendor centroids: 36,323 raw bins become 5,782
  centroids. Real peptide peaks span roughly 9–37 adjacent bins. The vendor centroid m/z is
  within about 0–15 ppm of a naive intensity-weighted mean of each run of adjacent bins. The
  naive runs merge overlapping isotope peaks, so this is an upper bound on the difference.
  Vendor centroid intensity is about 0.4× the run's summed intensity and about 1.4× its
  apex. The vendor centroider does more than sum the run: OPEN, Phase 4.
- Raw intensity sum vs Idx TIC: the TIC is 0.085% lower (median). The TIC is not the plain
  sum of stored intensities: OPEN.

## 7. Not yet verified

- Scan window (mzML: MS1 400–1250, MS2 100–1500; probably `MassRangeEx`), collision energy,
  and the multi-sample `.wiff` layout.
- Whether older (TripleTOF/Analyst) files use the same block format. OpenSXRaw's 56-byte block
  header and `82 05 00 00` file header suggest they don't.
- The profile-vs-oracle decode check (brief §7 Phase 3). It needs a profile msConvert run; the
  existing mzML is centroided.

## 8. Compression, compared with Bruker tdf_bin (2026-09-24)

`scripts/compression_cmp.jl`. SCIEX is the paper run (192 M stored bins). Bruker is PXD053462
HeLa 1 ng on a timsTOF Pro (133 M stored peaks, codec 2). Bytes per stored peak:

| | SCIEX 7600 | Bruker timsTOF Pro |
|---|---|---|
| native, whole file | **2.40** | **3.34** |
| native peak payload | 2.33 (byte tokens, no entropy coder) | 3.12 (u32 words → byte planes → zstd) |
| the other vendor's codec on the same peaks | 1.72–1.85 (Bruker codec 2, per scan or per cycle) | 4.22 (SCIEX tokens); 3.35 with zstd added |
| zstd-3 on the native SCIEX tokens | 1.73–1.83 | |
| within-scan bin delta, median (q90) | 17 (459); 24.8% adjacent | 1372 (15,565); 0% adjacent |
| intensity, median (q99) | 40 (312) | 78 (197) |
| bin width at m/z 400 | ~9.8 ppm (8 raw units) | ~7.8 ppm (394,428 TOF bins, 100–1700) |

- **Both formats store integer (TOF-bin delta, intensity) pairs** and calibrate with √(m/z) linear
  in the bin. SCIEX's byte code is a real compression (3.4× smaller than Float32 pairs), but it
  has no entropy stage. zstd would save another ~25%, the same as Bruker's shuffle + zstd.
- **Bruker spends more bytes per peak because its peaks are sparser.** A single IM scan has
  almost no adjacent bins, and the median gap is 1,372 bins. SCIEX scans sum all ion mobilities,
  so peaks are profile runs (a median gap of 17 bins). "Peak" means a stored non-zero bin in
  both cases, not an ion.

## 9. Centroiding and mass accuracy (2026-09-24)

`src/Centroid.jl`. Evaluated with `scripts/mass_accuracy2.jl`, `fragment_accuracy.jl`,
`noise_filter.jl` and `vendor_centroid_study.jl` on the paper run.

**What the vendor centroider does (from its output only):**
- Intensity is `100 × Σ I·Δm/z`, the peak area in count·Da. On isolated peaks the median
  factor is 100.6 (IQR 97.7–105.7), with log-log exponents 1.01 on the summed intensity and
  0.52 on m/z. So vendor intensities carry a √(m/z) weight compared with a plain sum.
- **Every peak built from a single non-zero bin is dropped**, in MS1 and MS2 (0.0–0.1% of them
  are kept). 94–98% of peaks with ≥ 2 non-zero bins are kept.
- Its m/z sits within about ±1–1.4 ppm (MAD) of an intensity-weighted mean, pulled toward the
  most intense part of noisy peaks. We have not reproduced it exactly.

**Peak shape**: σ ≈ 1 step (8 raw units) at every m/z, and FWHM ≈ 2.3 steps (about 30 ppm at
m/z 200, about 15 ppm at 1000). Long runs of adjacent bins are tails, overlaps and noise.

**Our centroider**: gaps count as zeros. The profile is smoothed with a Gaussian of σ = 2 steps,
and every smoothed maximum is a peak. From the raw bins within ±3 steps (stopping at valleys) it
reports the intensity-weighted mean position and the area. A peak needs ≥ 2 non-zero raw bins.

**Mass accuracy against theoretical m/z** (DIA-NN 1% FDR precursors, 23,403). The error is the
most intense centroid within ±15 ppm, in the scan nearest the apex.

| | Vendor (msConvert) | Ours (default) |
|---|---|---|
| MS1 precursor MAD, low / mid / high intensity third | 3.61 / 3.70 / 2.70 ppm | 3.42 / 3.41 / 2.63 ppm |
| MS1 precursors found (of the vendor's) | 100% | 96.1% |
| MS2 fragments matched | 109,025 | 106,360 |
| MS2 fragment median / MAD | −8.09 / 4.15 ppm | −8.17 / 4.08 ppm |
| MS2 paired (same fragment) MAD, all / vendor-I > median | 4.01 / 3.27 | 4.03 / 3.31 (median difference 0.00) |
| MS2 centroids per scan | 442 | 417 |
| log-intensity correlation with vendor (paired fragments) | | 0.943, median ratio 0.94 |

Mass errors carry a systematic offset: −3 ppm in MS1 and −8 ppm in MS2. It appears the same in
vendor and ours, so it's the instrument calibration in this run, not our decoding. Pioneer's
mass-error model absorbs it.

Things that did worse: a weighted mean over the whole run (tails pull it, +0.5 ppm MAD);
splitting at gaps between stored bins (it fragments sparse MS2 peaks); and smoothing with
σ ≤ 1.5 steps for MS2.

## 10. Speed (Apple M-series, 14 cores; paper run, 461 MB `.scan`, 159,558 records)

| | 1 thread | 8 threads |
|---|---|---|
| open (CFB + Idx + windows) | 0.45 s | 0.43 s |
| decode all blocks (192 M bins) | 0.72 s (268 M bins/s) | 0.12 s |
| decode + centroid (25.5 M centroids) | 3.7 s | 0.62 s |

`.wiff.scan` is memory-mapped; each thread reuses one `ScanBuffer` and one `CentroidBuffer`
(no per-scan allocation). The centroider skips clusters with fewer stored bins than `min_bins`
and smooths by scattering each stored bin's kernel instead of convolving a mostly-zero array
(3.4× faster).
