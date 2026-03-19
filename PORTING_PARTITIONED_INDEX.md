# Porting the Partitioned Fragment Index to `feature/bypass-first-pass`

This document describes how to bring the partitioned fragment index search (from `feature/fast-fragment-index-search`) onto `feature/bypass-first-pass`. It is written as instructions for a Claude working on the target branch.

---

## What This Changes

The monolithic fragment index search in `LibrarySearch.jl` is replaced with a **partitioned fragment index** that achieves a **5-7x speedup** on the fragment index search phase. The partitioned index is built at search time from the `SpectralLibrary` (~16s), not stored in `.poin` format, so it requires zero changes to library loading or serialization.

### Why It's Faster

1. **5-Da precursor m/z partitions**: Fragments are grouped by their precursor m/z into ~142 partitions (5 Da each). Within a partition, ALL fragments unconditionally match the precursor m/z criterion, eliminating the per-fragment binary search on precursor m/z that dominates the monolithic search.

2. **UInt16 local IDs (4-byte fragments)**: Each partition has at most 65,535 precursors, so fragment entries use `UInt16` local IDs + `UInt8` score = 4 bytes (vs 12 bytes for `IndexFragment{Float32}`). A per-partition `local_to_global::Vector{UInt32}` maps back to global IDs after scoring.

3. **L1/L2-cache-friendly counters**: The scoring counter is `LocalCounter{UInt16, UInt8}` (~192 KB), fitting entirely in L1/L2 cache. The monolithic `Counter{UInt32, UInt8}` is indexed by global precursor ID (potentially millions of entries).

4. **SoA layout + SIMD search**: Fragment bin metadata uses struct-of-arrays (`SoAFragBins{T}`) with 4 contiguous arrays instead of `Vector{FragIndexBin{T}}`. The `highs` array is scanned with 8-wide SIMD (`F32x8` via `llvmcall`) for finding the first matching fragment bin. A hybrid binary→SIMD search (threshold=128) handles both small and large ranges efficiently.

5. **5-Da skip hints**: Each fragment bin stores a `UInt16` hint: "how many bins ahead is +5 Da in the `lows` array?" This enables O(1) lower-bound advancement between consecutive peaks (which are sorted by m/z). Only ~10% of peaks need exponential-doubling fallback.

6. **Partition-major threading**: The outer loop is over partitions, the inner loop fans scans across threads. All threads process the same partition simultaneously, sharing the partition's data in L2/L3 cache.

### What Stays the Same

- `.poin` library format and loading (`loadSpectralLibrary.jl`) — unchanged
- `getPSMS()` function — unchanged (receives same `scan_to_prec_idx` + `precursors_passed_scoring` format)
- `selectTransitions!`, `matchPeaks!`, `buildDesignMatrix!` — unchanged
- `queryFragmentIndex.jl` (monolithic search) — kept as fallback, not deleted
- `LibraryFragmentIndex.jl` types — kept for Arrow loading

---

## Files to Create (3 source + 2 test)

### 1. `src/structs/PartitionedFragmentIndex.jl`

**Source**: `feature/fast-fragment-index-search:src/structs/PartitionedFragmentIndex.jl`

Contains the production types:

```julia
# SoA layout for fragment bins
struct SoAFragBins{T<:AbstractFloat}
    lows::Vector{T}
    highs::Vector{T}          # padded with 7 Inf sentinels for SIMD safety
    first_bins::Vector{UInt32}
    last_bins::Vector{UInt32}
end

# 4-byte fragment: UInt16 local precursor ID + UInt8 score
struct LocalFragment
    local_id::UInt16
    score::UInt8
end

# Per-partition index
struct LocalPartition{T<:AbstractFloat}
    fragment_bins::SoAFragBins{T}
    rt_bins::Vector{FragIndexBin{T}}
    fragments::Vector{LocalFragment}
    local_to_global::Vector{UInt32}  # local_id → global prec_id
    n_local_precs::UInt16
    skip_hints::Vector{UInt16}       # per-frag-bin: bins to +5 Da in lows
end

# Top-level partitioned index
struct LocalPartitionedFragmentIndex{T<:AbstractFloat}
    partitions::Vector{LocalPartition{T}}
    partition_bounds::Vector{Tuple{T, T}}  # (prec_mz_min, prec_mz_max) per partition
    n_partitions::Int
end

# Type-correct counter for UInt16 keys (Pioneer's Counter hardcodes zero(T) which
# doesn't generalize cleanly)
mutable struct LocalCounter{I, C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
end

const HINT_LINEAR_THRESHOLD = UInt32(128)
const MAX_LOCAL_PRECS = 65535
```

Also defines: accessor functions (`getFragBins`, `getRTBins`, `getFragments`, `getSkipHints`, `getPartitions`, `getPartition`, `getNPartitions`), `get_partition_range()` (binary search for partition overlap), `inc!`/`reset!` for `LocalCounter`.

**Note on `getPrecID`/`getScore` for `LocalFragment`**: These are defined directly on the struct in this file (not inheriting from `Ion`). They return `UInt16` and `UInt8` respectively, matching what `searchFragmentBinUnconditional!` expects.

### 2. `src/Routines/SearchDIA/CommonSearchUtils/buildPartitionedIndex.jl`

**Source**: `feature/fast-fragment-index-search:src/Routines/SearchDIA/CommonSearchUtils/buildPartitionedIndex.jl`

**Main function**: `build_partitioned_index_from_lib(spec_lib::SpectralLibrary; partition_width=5.0f0, frag_bin_tol_ppm=2.5f0, rt_bin_tol=3.0f0, rank_to_score=UInt8[8,4,4,2,2,1,1], y_start_index=UInt8(4), b_start_index=UInt8(3))`

The build algorithm:
1. Scan all precursor m/z values → compute partition count
2. Assign precursors to partitions by `floor((prec_mz - min_prec_mz) / partition_width)`
3. Split any partition exceeding 65,535 precursors (UInt16 limit)
4. For each partition: iterate precursors → extract `DetailedFrag` entries → apply ion-type/isotope filters → assign rank-based scores → build `SimpleFrag{Float32}` entries with local UInt16 IDs
5. Build `LocalPartition`: sort by iRT → bin by RT tolerance → sort by frag m/z → bin by PPM tolerance → build SoA arrays → pad `highs` with 7 Inf sentinels → compute skip hints

**Internal helpers** (also in this file):
- `_build_local_partition()` — builds one `LocalPartition` from `SimpleFrag[]`
- `_build_local_frag_bins!()` — builds frag m/z bins with SoA layout
- `_compute_skip_hints()` — binary search to find "+5 Da in lows" per bin

**Key dependency**: Uses `SimpleFrag{Float32}` (defined in `LibraryIon.jl`), `FragIndexBin{Float32}` (from `LibraryFragmentIndex.jl`), `DetailedFrag` (from `LibraryIon.jl`), `getPrecFragRange` (from `LibraryIon.jl`).

**Accessor usage**: `getMz(dfrag)` for `DetailedFrag` (lowercase z, defined on `AltimeterFragment`), `getMZ(x)` for `SimpleFrag` (uppercase Z, defined on `Ion`), `getPrecMZ(x)` for `SimpleFrag` (defined on `LibraryFragmentIon`), `getPrecID(x)` for `SimpleFrag` (defined on `LibraryFragmentIon`).

### 3. `src/Routines/SearchDIA/CommonSearchUtils/searchPartitionedIndex.jl`

**Source**: `feature/fast-fragment-index-search:src/Routines/SearchDIA/CommonSearchUtils/searchPartitionedIndex.jl`

Contains all search logic:

**SIMD primitives** (at top of file):
```julia
const F32x8 = NTuple{8, Core.VecElement{Float32}}
_vbroadcast8(x)     # broadcast Float32 to 8-wide vector
_vload8(arr, i)     # unsafe load 8 Float32s from array
_vcmpge_mask(a, b)  # SIMD >= comparison, returns UInt8 bitmask via llvmcall
_find_first_ge(highs, start, stop, threshold)  # SIMD scan + scalar tail
_findFirstFragBin_hybrid(highs, lb, ub, frag_min, simd_cutoff)  # binary→SIMD
```

**Scoring**:
- `searchFragmentBinUnconditional!(counter, fragments, range)` — scores all fragments in range, no precursor m/z binary search. This is the key performance change vs the monolithic `searchFragmentBin!`.

**Hinted search**:
- `queryFragmentHinted!(counter, max_idx, lb, ub, soa_frag_bins, fragments, frag_min, frag_max, hints, prev_mz, threshold)` — the core per-peak search function. Uses hint-based lb advancement, hint-based ub guess, exponential doubling fallback, hybrid binary→SIMD find.
- `_score_partition_hinted!(local_counter, partition, irt_lo, irt_hi, masses, mass_err_model)` — scores one partition across all peaks in a scan.
- `_find_rt_bin_start(rt_bins, irt_low)` — binary search for first RT bin with high >= irt_low.

**Orchestrator**:
- `searchFragmentIndexPartitionMajorHinted(scan_to_prec_idx, pfi, spectra, all_scan_idxs, n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)` — the partition-major entry point. Pre-computes per-scan properties, builds partition→scan mapping, spawns n_threads workers. Each worker: for each partition, score interleaved scan slices, filter by min_score and prec_mz bounds, push `(scan_idx, global_prec_id)` tuples. After all workers finish, collects results into `scan_to_prec_idx` and returns flat `precursors_passed` vector.

### 4. `test/UnitTests/partitionedFragmentIndex.jl`

**Source**: `feature/fast-fragment-index-search:test/UnitTests/partitionedFragmentIndex.jl`

Tests for: type construction, SIMD `_find_first_ge`, `_findFirstFragBin_hybrid`, `searchFragmentBinUnconditional!`, `queryFragmentHinted!` vs brute-force reference, threshold sweep, `_score_partition_hinted!` consistency, `_find_rt_bin_start`, `get_partition_range`, stress test (2000 bins, 20 peaks).

The test file defines helper functions (`make_soa_frag_bins`, `make_soa_hints`, `make_test_partition`, `brute_force_query!`, `extract_local_scores`) used to build synthetic test data.

**Imports needed** (add to top of test file or runtests.jl):
```julia
using Pioneer: SoAFragBins, LocalFragment, LocalCounter, LocalPartition,
    LocalPartitionedFragmentIndex, FragIndexBin, MassErrorModel,
    HINT_LINEAR_THRESHOLD, MAX_LOCAL_PRECS,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getPrecID, getScore, getLow, getHigh, getSubBinRange,
    F32x8, _vbroadcast8, _vload8, _vcmpge_mask, _find_first_ge,
    _findFirstFragBin_hybrid, searchFragmentBinUnconditional!,
    queryFragmentHinted!, _score_partition_hinted!, _find_rt_bin_start
```

### 5. `test/UnitTests/buildPartitionedIndex.jl`

**Source**: `feature/fast-fragment-index-search:test/UnitTests/buildPartitionedIndex.jl`

Tests for: `_build_local_partition` from synthetic `SimpleFrag[]`, `_compute_skip_hints` correctness, empty partition handling, `LocalPartitionedFragmentIndex` construction + `get_partition_range`.

**Imports needed**:
```julia
using Pioneer: SoAFragBins, LocalFragment, LocalPartition, LocalPartitionedFragmentIndex,
    FragIndexBin, MAX_LOCAL_PRECS, build_partitioned_index_from_lib,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getLow, getHigh, getSubBinRange,
    _build_local_partition, _compute_skip_hints, SimpleFrag
```

---

## Files to Modify (3)

### 1. `src/importScripts.jl`

Add `"PartitionedFragmentIndex.jl"` to the struct includes, immediately after `"LibraryFragmentIndex.jl"`:

```julia
include_files!(
    joinpath(package_root, "src","structs"),
    [
        ...
        "LibraryFragmentIndex.jl",
        "PartitionedFragmentIndex.jl",    # ← ADD THIS LINE
        "IsotopeTraceType.jl",
        ...
    ]
)
```

The two new files in `CommonSearchUtils/` (`buildPartitionedIndex.jl`, `searchPartitionedIndex.jl`) are picked up automatically by the `safe_include_directory!` call on the `CommonSearchUtils` directory — no explicit include needed.

### 2. `src/Routines/SearchDIA/LibrarySearch.jl`

This is the critical integration point. The target branch has **4 functions** that call the monolithic search and need updating:

#### a. `LibrarySearch()` (main search function)

**Before** (target branch): Calls `searchFragmentIndex()` per thread → `fetch.(tasks)` → per-thread `precursors_passed_scoring[thread_id]` → `getPSMS()`.

**After**: Build partitioned index → `searchFragmentIndexPartitionMajorHinted()` (single call, handles threading internally) → flat `precursors_passed_scoring` → `getPSMS()`.

Replace the body of `LibrarySearch()` with:

```julia
function LibrarySearch(
    spectra::MassSpecData,
    ms_file_idx::UInt32,
    fragment_index::FragmentIndex{Float32},  # kept in signature for API compat
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    nce_model::NceModel{Float32},
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel,
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures,
        P<:FragmentIndexSearchParameters}

    thread_tasks = partition_scans(spectra, Threads.nthreads())
    n_threads = length(thread_tasks)
    precursor_mzs = getMz(getPrecursors(spec_lib))

    # Build partitioned index from library
    # For presearch (irt_tol >= typemax(Float32)), use infinite RT bin tolerance
    rt_bin_tol = irt_tol >= typemax(Float32) ? typemax(Float32) : 3.0f0
    partitioned_index = build_partitioned_index_from_lib(spec_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=2.5f0, rt_bin_tol=rt_bin_tol)

    # Collect all valid MS2 scan indices
    all_scan_idxs = Int[]
    for tt in thread_tasks
        append!(all_scan_idxs, last(tt))
    end
    filter!(si -> si > 0 && si <= length(spectra) && getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)

    # Partition-major hinted search (single call, handles threading internally)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    precursors_passed_scoring = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    # getPSMS phase: fan out across threads as before
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            return getPSMS(
                ms_file_idx, spectra, last(thread_task),
                getPrecursors(spec_lib), getFragmentLookupTable(spec_lib),
                nce_model, scan_to_prec_idx, precursors_passed_scoring,
                search_data[thread_id], params, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end

    return fetch.(tasks)
end
```

**Critical change in getPSMS call**: The old code passed `precursors_passed_scoring[thread_id]` (per-thread vector). The new code passes the single flat `precursors_passed_scoring` vector. This works because `scan_to_prec_idx[scan_idx]` contains ranges that index directly into the flat vector, and each thread only accesses scans in its own `thread_task`.

#### b. `LibrarySearchNceTuning()`

Same pattern: replace `searchFragmentIndex()` per-thread calls + `materialize()` with one `build_partitioned_index_from_lib()` + one `searchFragmentIndexPartitionMajorHinted()`. The NCE grid loop over `getPSMS()` stays the same, but passes the flat `precursors_passed_scoring` instead of `precursors_passed_scoring[thread_id]`.

#### c. `fragment_index_search_only()` — TARGET BRANCH ONLY

This function exists only on `feature/bypass-first-pass`. It runs Phase 1 (fragment index search) without getPSMS scoring. **Update it** to use the partitioned index:

```julia
function fragment_index_search_only(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:FragmentIndexSearchParameters}
    spec_lib = getSpecLib(search_context)
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt_spline = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = haskey(getIrtErrors(search_context), ms_file_idx) ? getIrtErrors(search_context)[ms_file_idx] : Float32(Inf)
    precursor_mzs = getMz(getPrecursors(spec_lib))

    thread_tasks = partition_scans(spectra, Threads.nthreads())
    n_threads = length(thread_tasks)

    # Build partitioned index
    rt_bin_tol = irt_tol >= typemax(Float32) ? typemax(Float32) : 3.0f0
    partitioned_index = build_partitioned_index_from_lib(spec_lib;
        partition_width=5.0f0, frag_bin_tol_ppm=2.5f0, rt_bin_tol=rt_bin_tol)

    # Collect valid scan indices
    all_scan_idxs = Int[]
    for tt in thread_tasks
        append!(all_scan_idxs, last(tt))
    end
    filter!(si -> si > 0 && si <= length(spectra) && getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    precursors_passed = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    return scan_to_prec_idx, precursors_passed
end
```

This is simpler than the old version because `searchFragmentIndexPartitionMajorHinted` already returns a single flat vector with correctly-indexed `scan_to_prec_idx` ranges — no per-thread merge needed.

#### d. `searchFragmentIndex()` and `getRTWindow()`

**Keep `searchFragmentIndex()` unchanged** — it serves as fallback and is tested in existing tests.

**Keep `getRTWindow()` unchanged.**

**Keep `write_fragment_index_matches()` and `load_fragment_index_matches()` unchanged** — they work on `scan_to_prec_idx` + `precursors_passed` which have the same format.

### 3. `test/runtests.jl`

Add includes for the two new test files after the existing unit test includes:

```julia
    include("./UnitTests/RazoQuadModel.jl")
    include("./UnitTests/partitionedFragmentIndex.jl")  # ← ADD
    include("./UnitTests/buildPartitionedIndex.jl")      # ← ADD
```

### 4. `src/structs/MassErrorModel.jl` — Float32 type-conversion fix

The target branch still uses `1e6` (Float64 literal) in `getCorrectedMz` and `getMzBoundsReverse`. This causes Float32→Float64→Float32 promotion on every call in the hot path. Fix:

```julia
# getCorrectedMz: change 1e6 → 1f6
function getCorrectedMz(mem::MassErrorModel, mz::Float32)
    return mz - getMassOffset(mem)*(mz/1f6)
end

# getMzBoundsReverse: change 1e6 → 1f6
function getMzBoundsReverse(mem::MassErrorModel, mass::Float32)
    ppm = mass/1f6
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return mass - r_tol, mass + l_tol
end

# getMzBounds: change 1e6 → 1f6
function getMzBounds(mem::MassErrorModel, mass::Float32)
    ppm = mass/1f6
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return mass - l_tol, mass + r_tol
end
```

Also remove the commented-out `(mem::MassErrorModel)(mass::Float32)` callable if it still has `Float32(1e6)` — change to `1f6` or delete the dead code.

---

## What NOT to Port

- **`experiments/PartitionedFragIndex/`** — All experiment files. These are development-only benchmarking/profiling scripts. They stay on the experiment branch.
- **`NativeFragmentIndex` and `materialize()`** — The target branch doesn't have these. The partitioned index is built directly from the `SpectralLibrary`, not from a materialized fragment index. No need to port them.
- **Experimental types from `partitioned_types.jl`**: `PartitionedFragmentIndex`, `CompactFragment`, `CompactPartition`, `CompactPartitionedFragmentIndex`, `BitmaskCounter` — these were intermediate experiments, only the `Local*` types are production.
- **Non-hinted search functions**: `searchScanPartitioned!`, `searchFragmentIndexPartitioned`, `searchFragmentIndexPartitionMajor`, `searchFragmentIndexPartitionMajorBitmask`, `queryFragmentPartitioned!`, `_score_partition!`, `_score_partition_bitmask!`, all bitmask variants.

---

## Verification Steps

1. **Module loads**: `julia --project=. -e 'using Pioneer'` — should not error
2. **Unit tests**: `julia --project=. -e 'using Pkg; Pkg.test()'` — all existing + new tests pass
3. **Integration test**: `julia --threads=auto --project=. -e 'using Pioneer; SearchDIA("./data/ecoli_test/ecoli_test_params.json")'` — should complete without error

---

## Potential Merge Conflicts

1. **`src/importScripts.jl`**: Both branches modify the struct include list. Resolution: just ensure `"PartitionedFragmentIndex.jl"` appears after `"LibraryFragmentIndex.jl"`.

2. **`src/Routines/SearchDIA/LibrarySearch.jl`**: Major divergence. The target branch has extra functions (`fragment_index_search_only`, `write_fragment_index_matches`, `load_fragment_index_matches`) and a modified `getPSMS`. Resolution: keep all target branch additions; replace only `LibrarySearch()`, `LibrarySearchNceTuning()`, and `fragment_index_search_only()` bodies as described above.

3. **`test/runtests.jl`**: Both branches add test includes. Resolution: add both sets of includes.

4. **`src/structs/MassErrorModel.jl`**: If the target branch hasn't changed this file, just apply the `1e6` → `1f6` fix. If it has changes, merge carefully — the fix is just replacing the literal, not restructuring.

---

## Architecture Diagram

```
LibrarySearch() / LibrarySearchNceTuning() / fragment_index_search_only()
    │
    ├── build_partitioned_index_from_lib(spec_lib)        [buildPartitionedIndex.jl]
    │       │
    │       ├── Assign precursors to partitions by prec_mz
    │       ├── For each partition: extract DetailedFrags → SimpleFrag → LocalFragment
    │       ├── Build SoAFragBins (SoA layout, SIMD-padded highs)
    │       ├── Compute skip hints (5-Da in lows)
    │       └── Return LocalPartitionedFragmentIndex
    │
    ├── searchFragmentIndexPartitionMajorHinted(...)       [searchPartitionedIndex.jl]
    │       │
    │       ├── Pre-compute per-scan: irt_lo/hi, prec_min/max
    │       ├── Build partition→scan mapping
    │       ├── Spawn n_threads workers:
    │       │       For each partition k:
    │       │           For interleaved scan slices:
    │       │               _score_partition_hinted!(local_counter, partition, ...)
    │       │                   │
    │       │                   ├── _find_rt_bin_start() — binary search for RT bin
    │       │                   ├── For each mass in sorted peaks:
    │       │                   │       queryFragmentHinted!(...)
    │       │                   │           ├── Hint-based lb advancement
    │       │                   │           ├── Hint-based ub guess
    │       │                   │           ├── Exponential doubling fallback (10%)
    │       │                   │           ├── _findFirstFragBin_hybrid() / _find_first_ge()
    │       │                   │           └── searchFragmentBinUnconditional!() — score all
    │       │                   │
    │       │               Filter by min_score + prec_mz bounds
    │       │               Map local_id → global_id via local_to_global
    │       │               Push (scan_idx, global_pid) to thread results
    │       │
    │       └── Collect thread results → scan_to_prec_idx + precursors_passed
    │
    └── getPSMS(scan_to_prec_idx, precursors_passed, ...)  [unchanged]
```

---

## Parameter Defaults (benchmark-validated)

| Parameter | Value | Reason |
|-----------|-------|--------|
| `partition_width` | `5.0f0` Da | Safe for wide acquisition windows; ~142 partitions for typical libraries |
| `frag_bin_tol_ppm` | `2.5f0` ppm | Tight selectivity, essentially free (no speed cost vs 10 ppm) |
| `rt_bin_tol` | `3.0f0` iRT | 30% faster than 1.0 iRT; downstream `selectTransitions!` does exact iRT filtering anyway |
| `HINT_LINEAR_THRESHOLD` | `UInt32(128)` | Optimal on Astral data after Float32 type fix; pure SIMD up to 128 elements |
| `rank_to_score` | `UInt8[8,4,4,2,2,1,1]` | Same as monolithic index: 7 fragments per precursor, weighted by rank |
