# Plan: fold TimsSlices.jl and SciexWiff.jl into Pioneer.jl

Status: in progress (2026-09-28). Decisions agreed with Nathan are marked **decided**.

## Why

Pioneer is the only consumer of both packages. Folding them in removes the unregistered dependency
(`[sources]`), the General registration review (TimsSlices #169479, comment posted asking to close it; SciexWiff
#169493 had already merged as v0.1.0, no further versions planned there), and lets CI and the precompile workload
cover `.d` -> `.tdfs` / `.wiff` -> `.scxs` conversion and search, which neither covers today.

## Where

- **decided** `src/vendor/TimsSlices/`, `src/vendor/SciexWiff/`: each keeps its own `module ... end` (both export
  clashing names: `ConvertParams`, `convert`, `bin_to_mz`), included from `src/Pioneer.jl`; call sites stay
  qualified (`TimsSlices.open_tdfs`). `convert` renamed `convert_run` in both.
- TimsSlices lands on the timsTOF branch (PR #502, `feat/tdfs-reader-clean`); SciexWiff on the SCIEX branch (PR #503).
- Files copied with the source commits recorded (TimsSlices `b3e6d12`, SciexWiff `e4d9097`); not copied:
  `bin/`, `bench/`, per-package Project/Manifest. `docs/format.md` (`.tdfs` spec) and the `.scxs` notes move here.

## Dependencies

Drop TimsSlices / SciexWiff from `[deps]`, `[compat]`, `[sources]`; add direct deps SQLite, JSON3, Zstd_jll, Mmap
(already in Pioneer's Manifest: no new packages or binaries).

## Tests

- Package tests move to `test/UnitTests/formats/` (synthetic: codec, smoothing, params, CFB/container parsing).
- Real-data equivalence tests stay opt-in (`PIONEER_FORMATS_TEST_DATA`, `PIONEER_FORMATS_TEST_BIG`).
- Fixture tests in CI: convert each fixture, compare to stored checksums, read with `TdfsMassSpecData` /
  `ScxsMassSpecData`, run `convertBruker` / `convertSciex` end to end. Skip with a message when the Zenodo
  fixtures are absent (`temp/zenodo`).

## Precompile (`src/build/snoop.jl`)

New targets: `convertBruker`, `SearchDIA_tdfs` (committed IM E. coli library), `convertSciex`, `SearchDIA_scxs`
(committed E. coli library; the SCIEX fixture is a three-proteome run). Verify with `--trace-compile` and a
first-search timing in a locally built app.

## Fixtures

| # | fixture | where | how |
|---|---|---|---|
| 1 | truncated E. coli `.d` (analysis.tdf + analysis.tdf_bin) | **decided** Zenodo | PXD070049 E. coli 50 ng 5 min diaPASEF; limited RT range; MS2 frames keep 1-2 windows (frames re-encoded, needs a TDF frame encoder); MS1 trimmed to the kept windows |
| 2 | checksums of converting #1 | repo | fixture script |
| 3 | `.tdfs` from #1 | generated | test / precompile time |
| 4 | ion-mobility E. coli test library | **decided** repo | committed E. coli test library + alphapept_ccs 1/K0 |
| 5 | patched `.wiff` | **decided** Zenodo | JPST002949 `BenchSample_B_nswath4_25ng` (three-proteome, ZenoTOF 7600; local copy, SHA-256 matches SciexWiff's manifest). `Idx` stream (8.6 MB, 159,558 records = 917 cycles x 174 slots) edited in place: kept records re-pointed, dropped records size 0 (read as empty). Stream length unchanged, so the compound file layout is untouched; file stays 17.5 MB |
| 6 | truncated `.wiff.scan` | **decided** Zenodo | 44-byte header + kept blocks (limited RT, 1-2 windows per cycle) |
| 7 | checksums of converting #5 + #6 | repo | fixture script |
| 8 | `.scxs` from #5 + #6 | generated | test / precompile time |
| 9 | precompile configs | repo | JSON |

Fixture scripts: `test/fixtures/tools/make_d_fixture.jl`, `make_wiff_fixture.jl`. Zenodo: a new version of record
21812589 (existing files + `ecoli_tims_d_fixture.zip` + `sciex_wiff_fixture.zip`), uploaded by Nathan; then update
`ZENODO_RECORD`, add the downloads, bump the cache key in `.github/actions/precompile-data`, and use the action in
the unit-test workflow too.

## Checks before done

1. `.tdfs` / `.scxs` from the full source runs byte-identical (old packages vs vendored), converter metadata aside.
2. Search results identical before/after on an E. coli `.tdfs` and a `.scxs` run.
3. All test parts pass locally and in CI; fixture tests pass.
4. Precompile trace shows the new methods; no first-search compile pause in a built app.
