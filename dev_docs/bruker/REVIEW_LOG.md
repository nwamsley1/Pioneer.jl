# timsTOF branch review log

Decisions from the interactive review of `feat/tdfs-reader` against `develop`. Nothing here is applied
until it is marked **applied** with a commit.

## Product decisions (2026-09-23)

- TimsSlices.jl goes public and is registered in Julia General; AGPL-3.0-or-later (same as Pioneer),
  license header in every source file.
- Keep the `.tdfs` name for now (= "TDF slices"); define it in `docs/format.md`. Revisit later.
- Library build: a Bruker/timsTOF yes/no option; when yes, ion mobility is predicted with `alphapept_ccs`.
  **Applied (GUI toggle -> `library_params.im_model`) `5bd43f747`.**
- Searching `.tdfs` data with a library that has no ion-mobility predictions: hard error with an
  interpretable message, raised before any search work. **Applied `3f14465dc`.**
- Mixed vendors in one search (`.arrow` + `.tdfs`): hard error. **Applied `3f14465dc`.**
- The GUI must be able to build Bruker libraries, convert `.d`, and search `.tdfs`.
  **Applied `8ff7fc1a3` (convertBruker executable, `pioneer convert-bruker`) + `5bd43f747` (GUI).**
  Type-checked and Rust-tested only; not yet clicked through in the app. Open: the library calibration-file
  picker does not take `.tdfs`; no precompile target for convertBruker / a `.tdfs` search (needs a fixture).
- Before merging: verify parameter tuning is not broken or degraded on non-Bruker datasets
  (regression tests plus local analyses).
- Rebase onto `develop` before opening the PR.
- Remove all development-only logging, profiling and data dumps.
- After the PR: collect a varied set of Bruker datasets (PRIDE and other repositories; different ramps,
  gradients, loads, instruments) and check that the parameters generalise.

## Parameter overrides (`PIONEER_*` env hooks)

| hook | decision | status |
|---|---|---|
| `CHROM_IM_SIGMA` | Delete. Dead: default 0 leaves `im_lines` empty, so the library-line gate in `collect_rt_window_precursors!` never runs. Removes `CHROM_IM_TOL_SIGMA`, `chrom_im_tol_sigma`, `im_lines_by_charge`, the `im_scan/im_lib/im_lines/im_tol_sigma` args, and their tests in `test_chrom_im_gate.jl`. | **applied** `723dee6ee` |
| `CHROM_RT_TOL` | Delete. Dead: default 0 means the override branch never runs. | **applied** `723dee6ee` |
| `TUNING_IRT_TOL_MULT` | Delete (default 1.0 is a no-op multiply). | **applied** `723dee6ee` |
| `TUNING_SECOND_ORDER` | Delete the hook and the "off" branch; second-order stays on. | **applied** `723dee6ee` |
| `TUNING_MAX_PER_PREC` | Keep `TUNING_MAX_PSMS_PER_PRECURSOR = 3` as a constant; drop the ENV wrapper. Covered by the tuning regression check. | **applied** `723dee6ee` |
| `IRT_TOL_MULT` | Revert to develop's plain `4.0f0` default (the value was never changed). | **applied** `723dee6ee` |
| `CHROM_IM_SCANS` | Replace with `CHROM_IM_WINDOW_K0 = 0.055` (half-width, 1/K0); per file, scans = `ceil(0.055 / abs(slope))` using the file's calibration slope. See below. | **applied** `723dee6ee` |
| `IM_GATE_SIGMA` | **Keep 4σ** (A/B done 2026-09-23: neutral, nothing against it; see below). Drop the ENV wrapper. Re-examine on the post-PR dataset collection. | **applied** `723dee6ee` |
| `CHROM_IM_BAND_K0` | Keep 0.021. Convert band AND window to scans with the instrument slope from `.tdfs` meta, rounding up (option a, 2026-09-23). Band becomes 25 scans on current data (was 23). | **applied** `723dee6ee` (`getImSlope`, `im_half_width_scans`) |
| dump/debug hooks (`CHROM_DUMP_DIR`, `CHROM_ALT_FILES`, `CHROM_DUMP_ONLY`, `INDEX_DUMP_DIR`, `INDEX_DUMP_STRIDE`, `TUNING_DUMP_PSMS`, `STOP_AFTER`) | Delete with the code they guard. | **applied** `723dee6ee` |

### `CHROM_IM_SCANS` / `CHROM_IM_WINDOW_SCANS = 64`

What it does: in chromatogram extraction, a precursor is extracted from a slice only if the slice's IM
scan lies within ±64 IM scans of the IM scan of the precursor's best PSM (`collect_rt_window_precursors!`,
`utils.jl`). Precursors with no entry in `precursor_im_map` are not restricted. It is the mobility analogue of
the per-precursor RT window.

Findings:
- It is a half-width in raw IM scans, not in 1/K0. Every dataset measured so far used the same TIMS ramp
  (1/K0 0.64–1.45, 936–953 scans, 0.00085–0.000865 1/K0 per scan), so ±64 scans has only ever meant
  ±0.054–0.055 1/K0. A method with a narrower mobility range packs more scans per 1/K0 unit: at 0.85–1.30
  over ~930 scans, ±64 scans is about ±0.031, less than 2× the measured median FWHM (0.018). Untested and
  could clip peaks.
- The integration band (`CHROM_IM_BAND_K0`) is already in 1/K0 and converted per file. The window can do
  the same. The exact per-file slope is available from the instrument calibration (`meta.json`
  `im_slope_1overK0_per_scan`), which does not depend on a fit.
- The comment "64 scans (8 slices)" reads as the total width, but it is the half-width (±8 slices at stride 8).
- The comment also says ±3 slices holds 97–98% of weight, but ±4 slices (32 scans) clipped the apex for 4.5%
  of precursors. Checked 2026-09-23 (below): the centre is right; the far signal is mostly not the precursor's peak.

**Is the window centred on the mobility apex?** Measured with the ±96-scan dump (`pride_hye/dump_imwide`,
h50 A_REP1 / B_REP1, ~84k passing precursors each; script `proto/chrom2d/im_centre.jl`). Apex = argmax of the precursor's summed deconvolution weight per IM scan.

| | A_REP1 | B_REP1 |
|---|---|---|
| apex exactly at the best PSM's IM scan | 74.5% | 75.2% |
| apex within 1 slice (8 scans) | 93.8% | 94.7% |
| apex more than 32 scans away | 4.0% | 3.4% |
| apex more than 64 scans away | 1.89% | 1.58% |
| mean signed offset (scans) | −0.53 | −0.44 |
| weight within ±32 / ±48 / ±64 scans of centre (median, share of the ±96 total) | 87.9 / 93.7 / 96.9% | 89.1 / 94.4 / 97.3% |

Reading: the best PSM is centred on the apex with no bias. The ~2% whose "apex" is more than 64 scans away cannot be
the same mobility peak (median FWHM 0.018 1/K0 is about 21 scans), so it is most likely co-isolated interference
picked up by the deconvolution, not a clipped peak. That agrees with the ±96 re-search changing nothing. The
~3% of weight beyond ±64 is outside the ±0.021 integration band anyway. So ±0.055 1/K0 is well supported.

**Decision (2026-09-23):** `CHROM_IM_WINDOW_K0 = 0.055`, rounded UP to whole IM scans per file. On current data
that gives 64 scans (slope 0.000865) or 65 (slope 0.00085), so behaviour is essentially unchanged. Add a test
for the conversion and the rounding.

### `IM_GATE_SIGMA` / `IM_GATE_TOL_SIGMA = 4.0`, `IM_GATE_SIGMA_MULT = (z2 1.0, z3 1.8, other 2.2)`

Added in `fc2598d8a` (2026-09-21).

What it does: in the fragment-index search, a candidate precursor is emitted for a slice only if its library 1/K0
is within `4 × mult(z) × σ` of the z2 calibration line evaluated at the slice's IM scan. Hot-path cost is one
load and one FMA per candidate (`passes_im_gate`, `PartitionedFragmentIndex/search.jl`), dispatched at compile
time (`EmitToBuffer{F,G}`, with `G = Nothing` when off). No efficiency concern.

Where the line comes from:
- Built per `library_search` call by `build_im_gate` (`LibrarySearch.jl`), from `getImModel(ctx, file)[2]`.
  Returns `nothing` (no gate) unless: the file has IM scans, the library has 1/K0, and a z2 line exists.
- Parameter Tuning fits the line (`ParameterTuningSearch.jl` ~L748) on target z2 PSMs, needing
  `TUNING_IM_MIN_CALIB = 50` of them. MainSearch refits all charges later (`add_im_error!`), but that is after
  its own `library_search`. So in practice every gated search uses tuning's z2 line.
- σ = 1.4826·MAD of residuals, per file, so the gate is in data units and scales with the ramp (unlike the
  old `CHROM_IM_WINDOW_SCANS`). At z2 σ ≈ 0.0127 1/K0 (h50 A_REP1) the gate is about ±0.051 1/K0 for z2 and
  ±0.091 for z3.

Evidence on record (commit message and code comment): at 4σ the gate removes 6–9% of candidates and 0.2–0.6%
of passing PSMs (3σ: 12–19% and 0.9–1.9%). z3/z4 sit on the z2 line with no offset but 1.5–1.9× (z3) and
2.2× (z4) its scatter. The MainSearch per-charge fit on h50 A_REP1 agrees: z3 1.72×, z4 1.83×, z1 2.9×.

Concerns:
1. **Only the cost side was measured.** Nothing records what the gate buys in final IDs, FDR, quant or
   runtime. A 6–9% cut in candidates is modest, and the gate loses 0.2–0.6% of true PSMs. `im_error` is also a
   scoring feature, so wrong-mobility candidates are already penalised downstream. Needs a gate-on/off A/B.
2. **Tuning feeds itself.** Once tuning fits the line, later tuning iterations are gated by it. The line is
   fitted by ordinary least squares (outlier-sensitive; only σ is robust) from as few as 50 PSMs. A poor early
   fit could prune true candidates in the rest of tuning. Relevant to the "tuning not degraded" check.
3. **Multipliers are tied to the predictor.** They are alphapept_ccs's error structure (every library so far:
   HYE via `lib/add_im_to_lib.jl`, E. coli via BuildSpecLib; both alphapept_ccs). Fine since alphapept_ccs is
   now the fixed choice, but say so in the comment.
4. `other = 2.2` covers z1 too, where measured scatter is 2.9×. Moot at default `min_charge = 2`.
5. Stale docstring on `ImGate` (`search.jl`): says "per-charge line ... pooled line where a charge has none",
   but every line is the z2 line scaled.
6. No unit test for `build_im_gate` / `passes_im_gate`. The only test touching the multipliers
   (`test_chrom_im_gate.jl`) goes with the dead chromatogram gate.

Side finding for Phase 4: the HYE library used in the comparison was NOT built with `im_model`. IM was added
afterwards with `~/BrukerTims/lib/add_im_to_lib.jl`. The E. coli library was built end to end (`im_model:
alphapept_ccs`). The full BuildSpecLib path has therefore run once, on E. coli only.

### `CHROM_IM_BAND_K0 = 0.021` and the scan-conversion slope

The band is converted to IM scans in `IntegrateChromatogramsSearch.jl` (~L446) as
`max(1, round(Int, band_k0 / |b|))`, where `b` is the slope of MainSearch's *pooled fitted line* (library 1/K0
vs IM scan). That fitted slope is not the instrument's: on h50 A_REP1 the fit gives −0.000926 1/K0/scan and
the `.tdfs` calibration (`meta.json`) gives −0.000865, about 7% apart. Currently ±0.021 → round(22.7) = 23 scans. With
the instrument slope and ceil (the rule chosen for the window) it becomes ceil(24.3) = 25 scans.

For the PR, both the band and the new `CHROM_IM_WINDOW_K0` should use the same slope and the same rounding.
Options: (a) instrument slope + ceil for both. Physically right, but moves the tested band by about 2 scans.
(b) Keep the fitted slope for both. Matches what was tested, but depends on the fit and on the predictor.
`TdfsMassSpecData` does not expose the meta slope yet (it reads only `mz_lo/mz_hi` from meta).

## Review 1: `src/structs/MassSpecData/TdfsMassSpecData.jl`

Correction: format versions ARE checked. `TimsSlices.open_tdfs` rejects a mismatched `format_version`
(`codec/tdfs.jl:154`). Only its message could be friendlier ("re-convert with TimsSlices >= x").

### 1-A. Peak element type (`Union{Missing,Float32}`): decided 2026-09-23 — **applied `c0247a6f0`**
(outputs identical on timsTOF h50 and Astral; the `.tdfs` precompile workload item is still open)
The union is inherited, not needed: Thermo/mzML Arrow files declare peak columns nullable, and the hot kernels
declare exactly `AbstractArray{Union{Missing,Float32}}` (`fusedMatch.jl:488-489, 795, 1003-1004`;
`PartitionedFragmentIndex/search.jl:216-217`), which, arrays being invariant, rejects `Vector{Float32}`.
**Decision:** widen those signatures to `AbstractArray{<:Union{Missing,Float32}}` and make the `.tdfs` reader
use plain `Vector{Float32}`. Own commit. Verify: Thermo regression identical, `.tdfs` search identical.
Add a `.tdfs` search to the packaging precompile workload (the kernels get a second specialisation).
**Long term:** PioneerConverter should write plain Float32 peak columns (no nullable), removing the union at
the source. Pre-existing on develop, not ours: `FilteredMassSpecData.getMzArray` copies into a Union vector on
every call for the same reason.

### 1-B. Per-thread decode buffer: audit 2026-09-23
The buffer is picked by `Threads.threadid()`, and callers receive views into it. Two ways it can break:
(1) a task yields while holding views, another task runs on that thread and decodes into the same buffer;
(2) the same task fetches a second scan while still using the first one's views.

Audit of every `getMzArray`/`getIntensityArray` call site in `src/` (process_scans_fused, fragment-index
scorer, ParameterTuning/QuadTuning/HuberTuning mass-error and fused loops, IntegrateChromatograms MS2 and MS1
builders, MainSearch MS1 features, PrecursorScoring wide-window features, ms1_diagnostic, QuadTuning probe):
- No site fetches a second scan while the first scan's views are live. The MS1 feature paths copy into their own
  cache immediately.
- No yield-capable call (logging, locks, IO, channels, wait/fetch) sits between a fetch and the last use. The only
  logging in those files is outside the scan loops.
So it is correct today. But nothing enforces either rule; a single `@debug_l1` or lock added inside a scan
loop would silently corrupt peaks under contention.

Key point: in Julia the unit of ownership should be the **task**, not the thread. A task can move between
threads at a yield, and one thread can interleave several tasks. Pioneer already gives each spawned task its
own `search_data[first(thread_task)]`. The robust fix puts the decode buffer there.

**Decision (user, 2026-09-23):** one decode buffer per task; whatever thread runs the task uses that task's buffer.
**Implemented 2026-09-23 (uncommitted, working tree of `feat/tdfs-reader`):**
- `PeakDecodeBuffer` (was `TdfsSliceScratch`) + `getPeaks!(buf, spectra, scan) -> (mz, intensity)`. `.tdfs` decodes into
  `buf`; Arrow/Filtered data ignore it; `IndexedMassSpecData` forwards. Cache key is (file uid, scan), since task
  buffers are reused across files. The per-thread buffer vector is gone; `getMzArray`/`getIntensityArray` on
  `TdfsMassSpecData` now throw, so a missed call site fails loudly.
- `SimpleLibrarySearch.decode_buf` (+ `getDecodeBuffer`), one per task's search_data: main scoring loop, the three
  tuning searches, both chromatogram builders, `score_psms!` TIC sum.
- Own tasks, own buffers: fragment index (one per worker), MainSearch MS1 lookup (one per `parallel_foreach!`
  chunk), both wide-window loops (now chunked `parallel_foreach!`, buffer + peak scratch per chunk, which also
  removes develop's `threadid()` scratch there).
- Sequential sites with a local buffer: `FilteredMassSpecData` constructor and `append!`, QuadTuning probe,
  MS1 diagnostic.
- No `getMzArray`/`getIntensityArray` call remains in `src/` outside `structs/MassSpecData/`.

**Evidence:**
- `test_tdfs_mass_spec_data.jl`: 31/31 pass (-t 4). New: task stress (4× tasks per thread, `yield()` between
  fetch and check), (file, scan) cache key across two files, `getMzArray` on `.tdfs` throws, Arrow ignores buffer.
- Same stress against the old design (`faddc7040`): without yield 0/32,000 corrupted; **with one `yield()`,
  26,061/32,000 fetches (81%) returned another scan's peaks.** New design: 0/32,000.
- Thermo equivalence, `test/integration/search_ecoli.json`, old vs new code: precursors/protein groups, long and
  wide, all cell-for-cell IDENTICAL.
- `.tdfs` equivalence, h250 condA (3 files), old vs new code: precursors (23,113 × 144) / protein groups, long
  and wide, all cell-for-cell IDENTICAL. Total allocation 102.64 vs 102.63 GB, every step within 0.01 GB.
  Wall times (1171 vs 1319 s) are not comparable: both ran beside the gate A/B, and untouched steps moved as much
  as changed ones (library loading +30%, chromatogram integration −20%). Rerun on a quiet machine, 12 threads.
- Decode micro-benchmark (`proto/decode_bench.jl`, 1 thread, every slice of h250_A_REP1, 1,597,190 slices):
  old 8.08 µs/scan, new 8.14 µs/scan (noise; machine busy), identical checksum. End-to-end timing: pending a
  quiet machine.
- End-to-end, quiet machine, 12 threads, h250 condA (3 files) (`pride_hye/perf_ab.sh`): old 546.7 s, new 541.8 s.
  Every step within ~1 s (Main Search 174.9 vs 174.5, Parameter Tuning 147.1 vs 147.2); allocation 102.0 vs
  102.4 GB. Outputs again cell-for-cell IDENTICAL. **No performance cost.**
- Warm-start ABBA A/B (2026-09-23, `~/BrukerTims/perf_ab/`: old/new/new/old, one process each, warm-up searches
  first, MBR on, 12 threads), 2 runs per version: timsTOF h50 2 files 579.4 -> 577.7 s, 177.24 -> 177.89 GiB
  allocated; Astral Olsen HYE 200 ng 2 files 278.2 -> 276.6 s, 110.22 -> 110.55 GiB; GC time and peak RSS
  (~49 GiB) unchanged; outputs cell-for-cell IDENTICAL on both. **Committed `ea7b0117a`.**
- Found in passing (pre-existing): results depend on thread count. The same 3 files give 23,113 precursor rows at
  4 threads and 23,669 at 12 (old and new code agree at each count). Worth understanding before the PR; it
  affects reproducibility claims.
- Found in passing (pre-existing, TimsSlices): `read_slice!` allocates ~235 bytes per slice (inside
  `decode_slice!`: zstd / untranspose path), about 375 MB of garbage per full pass over a 1.6M-slice file. Fix in the
  Phase 0 TimsSlices work.

Other notes:
- Pre-existing on develop: `wide_window_features.jl:620` indexes scratch by `threadid()` inside a dynamic
  `@threads` loop (same pattern, same "no yield" reliance). `:static` loops (IntegrateChromatograms L544,
  MainSearch features L1058/L1427) are sticky and therefore safe.
- Performance aside: the fragment index is partition-major, so a `.tdfs` scan is decoded once per partition it
  is relevant to (the buffer caches one scan). With 25 Da windows that is about 2 decodes per scan. Measure
  before acting.

### IM gate A/B result (2026-09-23)

HYE, 12 files each, MBR off, branch HEAD `faddc7040`; `pride_hye/gate_ab.jl`, `proto/chrom2d/gate_ab_compare.jl`
(1% run and global q-value; quant needs ≥3 replicates per condition).

| | prec/run | PG/run | prec CV | E. coli prec ±0.5 | E. coli prot ±0.5 |
|---|---|---|---|---|---|
| 50 ng on | 83,648 | 9,862 | 11.0% | 74.3% | 85.1% |
| 50 ng off | 83,704 | 9,828 | 11.0% | 74.9% | 85.7% |
| 250 pg on | 6,962 | 1,898 | 19.3% | 63.4% | 77.6% |
| 250 pg off | 6,894 | 1,896 | 19.1% | 60.4% | 80.9% |

50 ng: no measurable difference. 250 pg: +1% precursors with the gate, accuracy differences in both directions
on ~3.4k quantified precursors (noise). Runtimes not comparable (50 ng arms overlapped other jobs; 250 pg: on 708 s,
off 660 s). Decision per user rule: keep.
