# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with `CommonSearchUtils` — the algorithmic core shared across SearchDIA's search methods.

## Overview

`CommonSearchUtils` houses the per-scan inner loop, RT/work distribution helpers, and isotope/spline helpers that every SearchDIA method invokes. The historical pipeline (`selectTransitions! → matchPeaks! → sort → buildDesignMatrix! → sortSparse! → ScoreFragmentMatches!`) was replaced by a single fused per-precursor scan pass; classic helpers and types (`SparseArray`, `FragmentMatch`, `UnmatchedIon`, `PrecursorMatch`, etc.) are gone.

## Files

| File | Purpose |
|---|---|
| `fusedMatch.jl` | `run_fused!` — the fused per-precursor match+score+build loop. Replaces the entire classic pipeline (transition selection, peak matching, design-matrix construction, scoring) in one pass that writes `SparseArrayFused` in CSC order. |
| `fusedScan.jl` | Low-level primitives: `FusedScratch` (per-thread scratch buffers), `SparseArrayFused` accessors, `finalize_column!`, `bsearch_hybrid`, `verify_mz_sorted`. |
| `getFragIsotopes.jl` | `getFragIsotopes!` for `PartialPrecCapture` / `FullPrecCapture` strategies (originally lived in the deleted `selectTransitions/fillTransitionList.jl`). |
| `isotope_utils.jl` | Shared isotope helpers (called by `MainSearch` and `IntegrateChromatograms`). |
| `buildRTIndex.jl` | RT index construction for chromatogram extraction. |
| `partitionThreadTasks.jl` | Scan-aware thread work distribution. |
| `deconvolutionArrayUtils.jl` | `initResiduals!`, sparse-matrix solver helpers (`AbstractSparseDesignMatrix` interface). |
| `rt_alignment_utils.jl` | Spline-based RT alignment helpers. |

## Fused match+score+build (`fusedMatch.jl`)

### Algorithm

`run_fused!(kind, Hs, unscored_psms, id_to_col, scratch, ...)` performs in a single pass per scan:

1. **Per-precursor pre-match filters** (extracted helpers):
   - `passes_irt_filter(prec_irt, scan_irt, irt_tol)`
   - `quad_window_with_iso_bounds(qfunc, prec_charge, iso_err_bounds)` → window
   - `passes_prec_mz_filter(prec_mz, low, high)`
2. **Iso-pass loop** (1 pass for `FusedStandard`, 3 for `FusedQuadEst`):
   - `setup_transmission!` populates `prec_trans_buf`
   - `frag_isotope_scale` applies per-pass scaling
3. **Per-fragment loop**: rank filter → `getFragIsotopes!` → conservative half-width
4. **Per-isotope loop**: `iso_mz_for` → `in_frag_mz_window` → `match_window` → `bsearch_hybrid` → `scan_for_nearest_in_window`
5. **Match recording**: `push_match!` / `push_miss!` into `FusedScratch`; `record_match!` (trait-dispatched) updates `unscored_psms`
6. **Column finalize**: `finalize_column!` flushes the precursor's scratch into `Hs` in CSC order with per-(col,row) deduplication

### Trait-based dispatch (`FusedSearchKind`)

| Kind | Used by | Iso passes | `do_prec_check` | Match recording |
|---|---|---|---|---|
| `FusedStandard` | `MainSearch`, NCE tuning | 1 | true | `apply_main_scoring!` updates `MainUnscoredPSM` |
| `FusedRTIndexed` | `IntegrateChromatogramsSearch` | 1 | false (filters applied upstream) | no-op (chromatogram points read post-deconv) |
| `FusedQuadEst` | `QuadTuningSearch` | 3 (one per precursor isotope) | false | no-op |

The compiler emits one specialization of `run_fused!` per concrete `K`, folding the trait-dispatched constants away.

### Pre-condition: m/z-sorted detailed_fragments

`run_fused!` requires per-precursor fragments to be **m/z-sorted within each precursor**. Production libraries built with current `BuildSpecLib` satisfy this automatically (see `sort_detailed_fragments_by_mz!` in `BuildSpecLib/build/build_poin_lib.jl`). Older libraries that fail the precondition must be rebuilt from FASTA — there is no in-place rescue path. `verify_mz_sorted` (in `fusedScan.jl`) is called once at MainSearch setup and throws with a clear message if the precondition fails.

## Reusable per-thread buffers (`FragIndexScratch`)

`searchFragmentIndexPartitionMajorHinted` (the partitioned fragment-index search) used to reallocate `Counter`s + Int32/UInt32 buffers per call (~50–100 MB per call × multiple calls per file). These are now lifted to `SearchContext.frag_index_scratch::FragIndexScratch` (`src/structs/Counter.jl`) — grown lazily via `prepare!` and reused across BitVecCalibration's adaptive batches and per-file MainSearch calls. On 5G heap-size-hint runs this dropped MainSearch GC time from ~12% → ~7% and total runtime ~9%.

## Performance patterns

- **`@inline` helpers** for the run_fused! inner loop (preserves SIMD + folding through trait dispatch)
- **Pre-allocated buffers** on `SimpleLibrarySearch` reused across scans within a file
- **`@inbounds @fastmath`** in hot loops after careful validation
- **Sparse storage** via `SparseArrayFused` (5-array CSC layout with packed UInt8 meta)
- **Counter-based scoring** in fragment-index search (branchless `inc!` / `or!`)
- **Counting-sort merges** in `_collect_frag_index_results!`

## Testing

Per-helper unit tests live in `test/UnitTests/`:

- `test_fused_prec_filters.jl` — 179 tests for the 10 extracted helpers
  (`passes_irt_filter`, `quad_window_with_iso_bounds`, `iso_mz_for`,
  `in_frag_mz_window`, `conservative_half_width`, `match_window`,
  `compute_ppm_err`, `push_match!`, `push_miss!`, `passes_prec_mz_filter`)
- `test_run_fused.jl` — 149 orchestration tests for `run_fused!` itself,
  built on a `make_fused_fixture` helper that constructs a complete
  synthetic `SparseArrayFused`/`FusedScratch`/`DensePrecMap`/peaks/fragments
  call site. Covers: empty range, all-match/all-miss, do_prec_check
  branches, FusedRTIndexed bypass, FusedQuadEst (3 columns/precursor),
  rank filter, frag_mz_bounds clipping, scratch growth, multi-iso,
  iso_anchor advancement, intensity-dependent tolerance stress, M+1/M+0
  iso overlap, finalize_column! per-(col,row) deduplication, missing
  peak intensity, empty fragment lists, full unscored_psms field
  contract (b/y/p ions, longest_y, best_rank, error accumulation)
- `test_sort_detailed_fragments.jl` — 22 tests for the BuildSpecLib m/z-sort helper
- `test_counter.jl` — 468 tests for `Counter`, `CountFilter`, `LUTFilter`,
  `PatternAccumulator`, `filter_counter!`

Total: **818 unit tests** on the fused match pipeline + Counter, all in <2s.

## Common pitfalls

- **Iteration order matters in the partitioned-index builder.** `build_partitioned_index_from_lib` (`src/structs/SpectralLibrary/PartitionedFragmentIndex/build.jl:109-134`) assigns bitmask bits via `rank += 1` based on iteration position, NOT `getRank(frag)`. Sorting `detailed_fragments` by m/z **before** building the index would silently re-bind which fragment owns which bit. `BuildSpecLib` builds the index first, then sorts.
- **Counter `reset!` vs reuse.** `Counter.reset!` is O(n_active) — only zeros touched IDs. Reused via `FragIndexScratch` across files; `_run_thread` calls `reset!(lc)` between scans within the same call.
- **`@view` over `SparseArrayFused.rowval[1:n_vals]`** in solver code — the trailing slots are stale from prior scans. Always slice to `n_vals`.

## Files removed (historical reference)

The following classic-pipeline helpers and types are gone — if you see them referenced in older docs/tests, they no longer exist:

- `matchPeaks.jl`, `buildDesignMatrix.jl`, `selectTransitions/`
- `SparseArray` (replaced by `SparseArrayFused`)
- `sortSparse!` (no longer needed; `run_fused!` writes CSC directly)
- `FragmentMatch`, `UnmatchedIon`, `PrecursorMatch` (replaced by direct `SparseArrayFused` writes for matches; `MassErrSample` for ParameterTuning)
- `ScoreFragmentMatches!`, `ModifyFeatures!` (replaced by `apply_main_scoring!` called inline from `run_fused!`)
- `ArrayDict` (replaced by `AbstractPrecursorMap{V}` with `DensePrecMap` / `SparsePrecMap`)
