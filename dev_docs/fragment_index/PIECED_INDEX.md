# Pieced fragment index

Status: implemented and verified on a 305M-precursor synthetic library (2026-10-07). Not yet wired into
BuildSpecLib or SearchDIA library loading: only the experiment harness (`PIONEER_SYNTH_*`, LibrarySearch.jl) builds
and searches pieced indexes so far.

## What it is

A partitioned fragment index stored as several **pieces**. Each piece is an ordinary
`LocalPartitionedFragmentIndex(32)` over a consecutive run of the initial precursor-m/z bins, so the pieces
together have exactly the partitions of the single index over all bins. The search loads one piece at a time,
searches every scan whose isolation window overlaps it, and frees it before loading the next. A scan's candidates
are the union over pieces; pieces hold disjoint precursors, so nothing is double counted.

## Code

- `src/structs/SpectralLibrary/PartitionedFragmentIndex/build.jl`
  - `initial_partitions(prec_mzs, partition_width)`: Step 1 of the builder (global ids per initial m/z bin), split
    out so a piece built from a consecutive run of bins gets the same partitions as the single index.
  - `build_partitioned_index_from_selection(...; initial_partition_pids)`: builds from a given run of bins.
  - `_build_local_partition`: `sizehint!(...; shrink = true)` on the SoA bin arrays. `resize!` kept one slot per
    fragment of capacity, so a freshly built index held ~2.5x its real size (60 GB vs 23 GB at 305M precursors);
    serialized/reloaded indexes were never affected.
- `.../PartitionedFragmentIndex/pieces.jl` (new)
  - `IndexPiece`, `PiecedFragmentIndex`, manifest `index_pieces.json`.
  - `build_index_pieces(sel, dir; partition_width, ..., id_type_request, max_piece_bytes)`: groups consecutive
    bins by an upper-bound size estimate (17 B per fragment + 4 B per precursor; conservative, real is ~9 B per
    fragment, so pieces come out at ~half the cap), rebuilds by halving if a piece is still too large, resolves the
    local ID type per piece, holds one piece in memory at a time.
  - `write_index_piece` / `read_index_piece`: raw arrays (header + per-partition vectors), loads at disk speed.
  - `index_bytes(pfi)`: exact array bytes of an index.
  - `for_each_index_piece(f, pfi, scan_prec_min, scan_prec_max)`: loads only pieces some scan window overlaps.
- `.../PartitionedFragmentIndex/search.jl`
  - `searchFragmentIndexPartitionMajorHinted` accepts a `PiecedFragmentIndex` as well as a single index.
    Per-thread result buffers and write positions (`thread_counts`) persist across pieces; `_search_index!` is
    the per-index function barrier (scan mapping, counters sized/typed per piece, threaded `_run_thread`, which now
    takes its starting write position). A single index is searched as its only piece (same code path).
- `src/importScripts.jl`: `pieces.jl` included between `build.jl` and `search.jl`.

## Verification

- Unit level (`experiments/immuno_synth/test_pieces.jl`, 3P library): 17 UInt16 pieces and 23 UInt32 pieces
  reproduce every partition of the single index exactly (fragments, bins, rt bins, skip hints, bounds) after a
  write/read round trip.
- Existing tests: partitionedFragmentIndex, buildPartitionedIndex, partitionedFragmentIndex32,
  test_build_fragment_index_exact, test_counter, test_build_determinism: 1079/1079 pass.
- Search level, RIS, Olsen Astral E40H50Y10_01, 305.5M-precursor synthetic library (real 3P + 298.7M fake
  non-enzymatic human 8-12mers), 5 Da partitions, UInt32, 32 threads:

  | | single index | 12 pieces (1.8-2.3 GB) |
  |---|---|---|
  | exact-candidate capture | md5 0cf7611c... | identical md5 |
  | 23 count fields (calibration, candidates, unique) | | identical |
  | index load | 94 s (compressed .jls) | 7.5 s (manifest) |
  | fragment index | 178 s | 125 s |
  | extra full pass (exact capture) | 88 s | 122 s |
  | bitvec calibration (1 batch) | 2.5 s | 20.6 s (every piece reloaded per call) |
  | peak RSS | 40.4 GB | 25.0 GB |

  One run each, on different nodes: the fragment-index speedup is suggestive, not established.

## Open questions / next steps

- **Do small pieces speed up ordinary searches?** Hypothesis: per-thread working set (counter + current
  partition's arrays) and memory traffic shrink, and piece-level scheduling is more cache friendly. Test on normal
  libraries (human / 3P, 5-10 M precursors) with pieces in the MB range, held IN MEMORY (no reload): the reload
  cost above is I/O, not search. A cheap first test: search a single in-memory index piece-by-piece (group its
  partitions) vs all at once, comparing fragment-index time at 5/10 Da on Astral + timsTOF data.
- Every call reloads every overlapping piece: fine for Main Search (one call per file), bad for tuning steps that
  call `library_search` many times on small scan subsets (parameter / quad / NCE tuning, bitvec batches). Either
  keep pieces resident when they fit, or tune on a small calibration library.
- Tighten the size estimate (~9 B per fragment) if fewer, larger pieces are wanted.
- Wire into BuildSpecLib (`write_fragment_indexes`) and library loading (`loadSpectralLibrary`,
  `getPartitionedIndex`) behind a size threshold; the presearch index too.
