# Fragment index: precursor partition width (2026-09-29)

The fragment index splits precursors into m/z partitions of `library_params.prec_partition_width` Da (a BuildSpecLib
parameter, fixed into the library; at the time of the sweep: default 5 Da, 10 Da when `im_model` is set). Each scan scores the precursors of
every partition overlapping its isolation window: partitions much wider than the window score out-of-window
precursors that are discarded later, partitions much narrower than it cost a fixed overhead per partition visited.
This sweep measures that trade-off on five datasets. Companion to [UINT32_LOCAL_IDS.md](UINT32_LOCAL_IDS.md).

## Method

- **Libraries differ only in the index.** Each dataset's library was rebuilt at every width by an index-only rebuild
  (`~/BrukerTims/pwsweep/rebuild_index.jl`): fragments and precursors are shared with the source library, each
  precursor's fragments are put back in rank order (the order `buildPionLib` indexes; `getRank` is stored in the
  compact fragments), and `build_partitioned_index_from_lib` runs with the BuildSpecLib constants. Rebuilding at the
  source's own settings reproduces its index byte for byte (HYE timsTOF 10 Da / UInt32) or field for field (the
  older prebuilt `Pioneer_Human_canon_std`, whose file differs only in serialisation).
- **Local IDs `"auto"`:** UInt16 when every partition fits 65,535 precursors (1 and 2.5 Da here), else UInt32, so
  every width is its nominal width.
- **Timing:** all searches of a session in one Julia process (`--threads 12 --gcthreads 8,1`) after warm-up searches
  covering both local-ID types on the data type; every (dataset, width) searched twice, the second pass in reverse
  order. Tables give the mean of the two passes (both shown); IDs were identical between passes. MBR on.
- Code: PR #502 at 36c40b88f. Sources: `results*.tsv` in `~/BrukerTims/pwsweep`, tables by `make_tables.py`.

## Summary

| data | isolation window | fastest width | default vs fastest |
|---|---|---|---|
| Bruker HYE 50 ng | 25 Da | 10 Da | 10 Da: fastest (2.5 Da +21.5%, 25 Da +13.3%) |
| Olsen Exploris 500 ng | 14.7 Th | flat (10 Da nominal) | 5 Da: +1.2%, within pass-to-pass spread |
| SCP 250 pg | 4.4 Th | flat (2.5 Da nominal) | 5 Da: +2.2%, within pass-to-pass spread |
| SCIEX ZT three-proteome | 2.9 Da | 2.5 ≈ 5 Da | 5 Da: +0.4% |
| Olsen Astral 200 ng | 2 Th | 2.5 ≈ 5 Da | 5 Da: +1.0% (session 2) |

- **The fastest width tracks the isolation window**, and it only matters where the fragment index is a large share
  of the run time (timsTOF, Astral, SCIEX). Partitions narrower than the window (1 Da) are always slower
  (+8-14%); very wide ones hurt narrow-window data (25 Da: +14-20%).
- **IDs have no systematic width dependence** (at most about 1% on precursors, no trend). Each library gives identical
  counts on repeat, so differences are the candidate order changing the downstream result, not noise.
- **Exception, low-load data:** on SCP 250 pg the IDs swing 5% at 1% FDR and 44% at 0.1% FDR (3,846-5,535) across
  widths with no trend: small index changes move the scoring result a lot there (see the tuning-fragility item).
- **Recommendation: keep the defaults** (5 Da; 10 Da for timsTOF), within about 1% of the fastest width everywhere tested.
  2.5 Da is never worse on narrow-window data but the gain is too small for a new default.

## Adopted rule

BuildSpecLib takes the approximate acquisition isolation window width, `library_params.isolation_window_width` (m/z,
GUI field "Isolation window width"), and uses `prec_partition_width = clamp(isolation_window_width, 2.5, 10)` Da.
Unset, it is 5 (5 Da partitions), which avoids the worst case at either end (within 4% of the fastest width on every
dataset above). An explicit `prec_partition_width` still overrides. The GUI prefills 25 for a timsTOF library.

## Results

### Bruker timsTOF Ultra 2, HYE 50 ng 15 min, 2 files (diaPASEF, 25 Da windows)

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.5 | UInt16 | 270 | 502.9 (514.1 / 491.7) | +21.5% | 359.3 | 223.1 | 172,893 | 139,567 | 19,713 | 18,298 |
| 5 | UInt32 | 135 | 430.4 (431.2 / 429.5) | +4.0% | 291.1 | 154.2 | 173,478 | 139,923 | 19,798 | 18,271 |
| 10 | UInt32 | 68 | 413.9 (414.7 / 413.2) | +0.0% | 270.7 | 128.2 | 174,733 | 139,765 | 19,917 | 18,381 |
| 13 | UInt32 | 52 | 432.6 (433.7 / 431.6) | +4.5% | 291.8 | 122.8 | 174,149 | 139,149 | 19,822 | 18,205 |
| 25 | UInt32 | 27 | 468.8 (483.8 / 453.8) | +13.3% | 322.3 | 138.7 | 174,205 | 141,280 | 19,901 | 18,478 |

### Thermo Astral (Olsen), HYE 200 ng 30 min, 2 files (2 Th windows), session 1

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.5 | UInt16 | 270 | 268.5 (272.7 / 264.3) | +2.0% | 130.5 | 22.0 | 414,360 | 319,213 | 26,954 | 25,874 |
| 5 | UInt32 | 135 | 263.2 (262.2 / 264.2) | +0.0% | 127.4 | 22.9 | 409,745 | 316,324 | 27,033 | 25,945 |
| 10 | UInt32 | 68 | 279.4 (283.2 / 275.5) | +6.1% | 136.8 | 29.1 | 412,272 | 316,356 | 27,043 | 26,018 |
| 13 | UInt32 | 52 | 280.2 (280.4 / 280.0) | +6.5% | 138.1 | 33.2 | 412,105 | 318,281 | 27,003 | 25,992 |
| 25 | UInt32 | 27 | 315.0 (318.5 / 311.5) | +19.7% | 169.9 | 53.5 | 409,449 | 317,421 | 26,923 | 25,688 |

### Thermo Astral (Olsen), same files, session 2 (adds 1 Da)

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | UInt16 | 675 | 288.5 (296.6 / 280.4) | +11.6% | 145.1 | 26.5 | 411,559 | 314,291 | 26,987 | 25,890 |
| 2.5 | UInt16 | 270 | 258.6 (258.2 / 259.0) | +0.0% | 126.1 | 19.4 | 414,360 | 319,213 | 26,954 | 25,874 |
| 5 | UInt32 | 135 | 261.3 (264.1 / 258.5) | +1.0% | 128.2 | 23.9 | 409,745 | 316,324 | 27,033 | 25,945 |

### Thermo Exploris (Olsen), HYE 500 ng 30 SPD, 2 files (14.7 Th windows)

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.5 | UInt16 | 270 | 106.1 (113.1 / 99.1) | +4.3% | 26.3 | 4.6 | 120,660 | 94,617 | 14,552 | 13,347 |
| 5 | UInt32 | 135 | 103.0 (103.5 / 102.4) | +1.2% | 26.8 | 4.8 | 120,429 | 94,385 | 14,551 | 13,315 |
| 10 | UInt32 | 68 | 101.7 (101.8 / 101.6) | +0.0% | 26.0 | 4.6 | 120,214 | 93,684 | 14,614 | 13,456 |
| 13 | UInt32 | 52 | 103.7 (104.1 / 103.3) | +2.0% | 27.2 | 4.7 | 120,397 | 93,428 | 14,625 | 13,305 |
| 25 | UInt32 | 27 | 104.3 (104.6 / 104.0) | +2.6% | 27.4 | 5.3 | 120,725 | 94,327 | 14,731 | 13,304 |

### SCIEX ZenoTOF 7600, three-proteome nswath4, 3 files (2.9 Da windows)

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | UInt16 | 675 | 322.5 (328.6 / 316.3) | +7.7% | 165.9 | 23.6 | 298,929 | 218,777 | 27,055 | 25,332 |
| 2.5 | UInt16 | 270 | 299.4 (299.5 / 299.4) | +0.0% | 155.6 | 17.9 | 298,830 | 218,756 | 27,268 | 25,636 |
| 5 | UInt32 | 135 | 300.5 (299.4 / 301.7) | +0.4% | 158.3 | 21.4 | 298,458 | 221,328 | 27,207 | 25,329 |
| 10 | UInt32 | 68 | 315.7 (316.0 / 315.4) | +5.4% | 170.3 | 24.4 | 298,485 | 220,871 | 27,152 | 25,757 |
| 25 | UInt32 | 27 | 342.5 (343.1 / 341.8) | +14.4% | 187.4 | 43.1 | 298,070 | 219,393 | 27,008 | 25,223 |

### Thermo single-cell (SCP), human 250 pg, 3 files (4.4 Th windows, FAIMS)

| width (Da) | local IDs | partitions | total s (pass 1 / 2) | vs fastest | Main Search s | fragment index s | precursors 1% | precursors 0.1% | protein groups 1% | protein groups 0.1% |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | UInt16 | 675 | 74.0 (80.6 / 67.4) | +14.0% | 13.1 | 0.5 | 11,187 | 4,197 | 3,836 | 2,370 |
| 2.5 | UInt16 | 270 | 64.9 (64.6 / 65.2) | +0.0% | 12.6 | 0.3 | 11,318 | 5,535 | 4,032 | 2,575 |
| 5 | UInt32 | 135 | 66.3 (68.8 / 63.9) | +2.2% | 12.1 | 0.3 | 11,030 | 3,846 | 3,734 | 2,581 |
| 10 | UInt32 | 68 | 65.3 (65.1 / 65.5) | +0.6% | 13.3 | 0.5 | 11,234 | 5,038 | 4,112 | 3,010 |
| 25 | UInt32 | 27 | 70.8 (72.6 / 68.9) | +9.0% | 14.0 | 1.1 | 11,588 | 4,393 | 4,113 | 3,087 |

