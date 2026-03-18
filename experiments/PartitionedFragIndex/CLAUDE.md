# Partitioned Fragment Index Experiment

## Goal

Speed up the fragment index search phase of Pioneer's DIA analysis by partitioning
the monolithic `NativeFragmentIndex` by precursor m/z. This eliminates the per-fragment
binary search on precursor m/z and enables cache-friendly scoring with small counters.

## Current Best Result

**5-Da hint + SoA + hybrid binary→SIMD + Float32 type-conversion fix: 6.53x speedup**
over the baseline monolithic index on Astral data (70.1s → 10.7s, 16 threads).
Best threshold shifted from 32 to 128 after eliminating Float64 promotion overhead.

## Architecture

The partitioned index splits fragments by precursor m/z into ~142 partitions (5 Da width).
Each partition uses:
- **LocalFragment** (4 bytes): `UInt16` local precursor ID + `UInt8` score
- **LocalCounter** (192 KB per thread): type-correct `UInt16`/`UInt8` counter that fits in L1 cache
- **Partition-major threading**: all threads process the same partition simultaneously for shared L2/L3 cache utilization
- **SoAFragBins**: struct-of-arrays layout for fragment bins (4 contiguous arrays instead of array-of-structs), enabling SIMD scans on the `highs` array
- **Hybrid binary→SIMD search**: binary search narrows the range, SIMD `_find_first_ge` finishes the last ≤32 elements

## Key Files

| File | Purpose |
|------|---------|
| `partitioned_types.jl` | Type definitions: `SoAFragBins`, `LocalFragment`, `LocalCounter`, `LocalPartition`, etc. |
| `build_partitioned_index.jl` | Build `LocalPartitionedFragmentIndex` from spectral library (SoA construction) |
| `search_partitioned_index.jl` | SIMD primitives, hybrid search, partition-major scoring |
| `benchmark_step_and_hints.jl` | Benchmark: baseline vs hybrid binary→SIMD with threshold sweep + correctness validation |
| `benchmark.jl` | Full benchmark comparing baseline vs partitioned (7-frag, 3-frag, bitmask) |
| `test_synthetic.jl` | 10 synthetic correctness tests |
| `test_hinted_search.jl` | Correctness tests for hinted search (threshold extremes, PPM partitions, etc.) |
| `profile_partitioned.jl` | PProf profiling of the partitioned search |
| `diagnose.jl` | Diagnostic tool for score mismatch analysis |
| `diagnose_hint_jumps.jl` | Diagnostic v2: validates 5-Da hint + advancing lb approach |

## Search Optimization Experiments (SearchOpt*)

Profiling showed **36% of search time** is spent in `findFirstFragmentBin` (binary search)
and **13%** in `exponentialFragmentBinSearch`. These are the main optimization targets.

### SearchOptA: m/z Lookup Table
- Pre-build coarse lookup table mapping m/z → approximate frag bin index
- One reciprocal multiply (no division) to get bucket index
- Linear scan from approximate position
- **Status**: implemented, not yet benchmarked standalone

### SearchOptB: Delta-Based Skip
- Use m/z delta between consecutive peaks to estimate frag bin jump distance
- Multiply delta by precomputed `inv_avg_bin_width` (one multiply, no division)
- Small deltas → linear scan; large deltas → jump + refine
- **Status**: implemented, not yet benchmarked standalone
- **Issue**: `inv_avg_bin_width` is a poor estimate because variable-width bins have wildly different widths

### SearchOptC: Fixed-Width Bins with O(1) Lookup
- Replace variable-width frag bins with fixed 0.005 Da bins
- CSR-style `UInt16` pointer array per RT bin maps bin index → fragment range
- Given peak m/z: `bin = floor(Int32, (mz - mz_min) * INV_BIN_WIDTH) + 1` → one multiply + one load
- **No exponential search, no binary search**
- **Result: 4.21x speedup** (17.2s vs 72.4s baseline, 22.8s variable-bin partitioned)
- Memory: 2.04 GB (1.85 GB frag_ptrs + 167 MB fragments + 24 MB local-to-global)

### SearchOptD: 5-Da Direct Hint + Advancing LB
- **Idea**: Two synergistic changes to the variable-bin hinted search:
  1. **New hint**: `hint[j] = k` where `getLow(frag_bins[j+k]) - getLow(frag_bins[j]) >= 5.0 Da`
     (direct measurement, not density average). Stored as UInt16. `est_step = hint[lb] * (delta_mz / 5.0)`
  2. **lb = first_matching_bin**: Return `(first_match, ub)` instead of `(saved_lb, ub)`.
     Provably safe: `findFirstFragmentBin` returns the first bin with `getHigh >= frag_mz_min`;
     all prior bins have `getHigh < frag_mz_min`, so next peak (higher m/z) can never match them.
  3. **Hint-based lb advancement**: Before searching, advance lb using the hint.
     `delta_mz > 5 Da`: advance by hint (provably safe — only 5 Da in getLow).
     `delta_mz <= 5 Da`: advance by `hint * delta_mz / 5 * factor`.
  4. **Hint-based UB guess**: `ub_guess = new_lb + est_step * overshoot_factor`.
     Exponential doubling only when guess is insufficient.
- **Diagnostic v2 results** (100 Astral scans):
  - Per-peak jump: median 13 bins (vs ~1409 with old non-advancing lb)
  - Hint accuracy (est/actual): median 0.86 — well-calibrated
  - UB guess sufficient: 90% (only 10% need exponential fallback)
  - **LB overshoot: 0.0% at ALL factors (0.5–1.0)** — linearity within 5 Da validated
  - Range sizes: median 17, 49% ≤16 (linear scan), 82% ≤64
  - Bin widths: tiny (median ~0.003 Da)
  - Recommended: `lb_factor=1.0`, `ub_overshoot_factor=1.5`
- **Result: 5.68x speedup** (12.6s, AoS layout)

### SearchOptE: SoA Layout + SIMD Linear Scan
- **SoA layout** (`SoAFragBins{T}`): Split `FragIndexBin` array-of-structs into 4 parallel
  arrays (`lows`, `highs`, `first_bins`, `last_bins`). Each `getHigh` scan now touches
  only the contiguous `highs` array (4 bytes/element) instead of striding through 16-byte
  AoS records. 4x better cache utilization for field-specific scans.
- **SIMD primitives**: `_find_first_ge(highs, start, stop, threshold)` uses `F32x8`
  (`NTuple{8, Core.VecElement{Float32}}`) with `llvmcall` for `fcmp oge` + `bitcast` to
  scan 8 elements per iteration. Scalar tail handles remaining elements.
- **Hybrid binary→SIMD** (`_findFirstFragBin_hybrid`): For ranges > threshold, branchless
  binary search narrows to ≤threshold elements, then `_find_first_ge` finishes with SIMD.
  For ranges ≤ threshold, direct SIMD scan (no binary overhead).
- **SIMD padding**: `highs` array gets 7 extra `Inf` sentinel elements for safe `_vload8`.
- **Result (before type-fix)**: 5.74x at threshold=32, 5.86x at threshold=1M (pure SIMD)
- Memory: ~585 MB (same as SearchOptD — SoA is same total size as AoS)

### SearchOptF: Float32 Type-Conversion Fix ← CURRENT BEST
- **`1e6` → `1f6` in MassErrorModel**: `getCorrectedMz` and `getMzBoundsReverse` used
  Float64 literal `1e6`, causing Float32→Float64 promotion then Float64→Float32 cast on
  every call. Replaced with `1f6` (Float32 literal, exactly representable). Removed
  unnecessary `Float32()` wrapping on return values.
- **`unsafe_trunc(UInt32, ...)` in hint computation**: Replaced `floor(Int, x)` →
  `UInt32(max(..., 1))` chain (two checked conversions) with `max(unsafe_trunc(UInt32, x),
  one(UInt32))`. Safe because `est_step_f` is always a small positive float.
- **`HINT_LINEAR_THRESHOLD` as UInt32**: Eliminates per-peak `UInt32(linear_threshold)`
  conversion in the inner loop.
- **Result: 6.52x at threshold=32, 6.53x at threshold=128 (new best)**
- Type conversion overhead reduced from 6.1% → 4.0% of active time
- Remaining conversions are in `_find_first_ge` (UInt32↔Int at entry/return) and
  UInt32 array indexing — inherent to Int-indexed SIMD design.

#### Hybrid Threshold Sweep Results (after type-conversion fix)

| Threshold | Time (s) | Speedup | Strategy |
|-----------|----------|---------|----------|
| 1 (pure binary) | 12.58 | 5.57x | Binary all the way |
| 8 | 11.72 | 5.98x | ≤8 SIMD, >8 binary→SIMD-8 |
| 16 | 11.25 | 6.24x | ≤16 SIMD, >16 binary→SIMD-16 |
| 32 | 11.06 | 6.34x | ≤32 SIMD, >32 binary→SIMD-32 |
| 64 | 10.97 | 6.39x | ≤64 SIMD, >64 binary→SIMD-64 |
| **128** | **10.73** | **6.53x** | **≤128 SIMD, >128 binary→SIMD-128** |
| 256 | 11.09 | 6.32x | ≤256 SIMD, >256 binary→SIMD-256 |
| 1M (pure SIMD) | 10.90 | 6.44x | SIMD the whole range |

## Performance Progression

| Version | Speedup | Search Time | Memory |
|---------|---------|-------------|--------|
| Baseline (monolithic, 7 frags) | 1.0x | 72–74s | 525 MB |
| Full frag bin copy (reference) | 1.46x | 49s | 3,367 MB |
| Variable bins, CompactFragment | 3.13x | 23s | 680 MB |
| Variable bins, LocalFragment + local counter | 3.12x | 23s | 541 MB |
| Partition-major threading | 3.17x | 23s | 541 MB |
| Fixed-width bins, O(1) lookup (SearchOptC) | 4.21x | 17s | 2,037 MB |
| 5-Da hint + advancing lb, AoS (SearchOptD) | 5.68x | 12.6s | ~585 MB |
| 5-Da hint + SoA + hybrid binary→SIMD (SearchOptE) | 5.74x | 13.0s | ~585 MB |
| **+ Float32 type-conversion fix (SearchOptF)** | **6.52x** | **10.8s** | **~585 MB** |

## Profile Breakdown (SearchOptF, hinted SoA+SIMD + type fix)

Based on 45,348 active samples in `_score_partition_hinted!`:

| Component | % of Search | Notes |
|-----------|-------------|-------|
| `getCorrectedMz` | 6.0% | Pure Float32 ops (no conversions after `1f6` fix) |
| `getMzBoundsReverse` | 7.3% | Pure Float32 ops (no conversions after `1f6` fix) |
| `queryFragmentHinted!` total | 81.0% | |
| — Hint computation (lb advance + ub guess) | 20.6% | `unsafe_trunc` + UInt32 arithmetic |
| — Exponential doubling fallback | 8.9% | ~10% of peaks need this |
| — SIMD/binary search | 13.0% | `_find_first_ge` + `_findFirstFragBin_hybrid` |
| — `searchFragmentBinUnconditional!` (scoring) | 18.6% | Fragment matching + `inc!` counter |
| — Fragment bin loop overhead | 13.6% | Bounds checks, `getSubBinRange` |
| Iteration + getMzArray | 6.3% | Mass array iteration + RT bin loop |
| Remaining type conversions | 4.0% | `_find_first_ge` UInt32↔Int, array indexing |

### Previous Profile (variable-bin partitioned, before SearchOptC)

| Component | % of time |
|-----------|-----------|
| `findFirstFragmentBin` (binary search) | 35.8% |
| `getCorrectedMz` + `getMzBoundsReverse` | 16.4% |
| `queryFragmentPartitioned!` overhead | 13.2% |
| `searchFragmentBinUnconditional!` (scoring) | 10.4% |
| Mass iteration loop | 8.7% |
| Other (RT bins, partition setup) | 15.5% |

## Correctness

The partitioned index produces slightly different results than the monolithic baseline
due to independent RT/frag bin construction per partition. On 500 Astral scans, ~496 show
ID set differences (a few dozen precursors per scan near bin boundaries). This is the
expected cost of independent binning, not a bug. End-to-end proteomics validation (PSMs
passing FDR) is the appropriate correctness metric.

The SIMD and binary search paths produce **identical results**: threshold=1 (pure binary)
vs threshold=1000000 (pure SIMD) validated on 500 scans with 92,359 precursor IDs,
0 mismatches in both ID sets and scores.

## Checkpoints (git tags)

- `checkpoint/compact-fragment-3x-speedup` — CompactFragment, 3.13x, safe fallback
- `checkpoint/local-counter-3x-speedup` — LocalFragment + local counter, 3.12x
- `checkpoint/partition-major-3.17x` — partition-major threading, 3.17x
- `checkpoint/pre-search-optimization` — before SearchOpt* experiments

## Running

### Prerequisites

1. **Pioneer.jl** checked out and instantiated (`julia --project=. -e 'using Pkg; Pkg.instantiate()'`)
2. **Astral DIA data** converted to Arrow format (`.arrow` files from `convertMzML`)
3. **Spectral library** in `.poin` format (built by Pioneer)

### Setup: Create a params.json

Create a JSON config file pointing to your data. The benchmark scripts only use the
calibration pipeline (ParameterTuning, NceTuning, QuadTuning) and fragment index search —
they do not run the full quantification pipeline. Example:

```json
{
    "paths": {
        "library": "/path/to/your_library.poin",
        "ms_data": "/path/to/arrow_files_directory",
        "results": "/path/to/output_results"
    },
    "acquisition": {
        "nce": 26,
        "quad_transmission": {
            "fit_from_data": true,
            "overhang": 0.25,
            "smoothness": 5.0
        }
    }
}
```

All other parameters use Pioneer defaults. The `ms_data` directory should contain the
`.arrow` files produced by `Pioneer.convertMzML()`. The `library` path should point to
a `.poin` spectral library directory. The `results` directory will be created if needed.

### Running benchmarks

All commands run from the Pioneer.jl project root:

```bash
# Synthetic correctness tests (no data needed)
julia --project=. experiments/PartitionedFragIndex/test_synthetic.jl

# Main benchmark: baseline vs hinted SoA+SIMD with threshold sweep + correctness validation
# This is the primary benchmark — reports speedup, threshold sweep, and binary-vs-SIMD validation.
julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_step_and_hints.jl /path/to/params.json

# Profile with PProf (generates .pb.gz flame graph)
julia --threads=auto --project=. experiments/PartitionedFragIndex/profile_partitioned.jl /path/to/params.json

# Full benchmark comparing baseline vs partitioned (7-frag, 3-frag, bitmask)
julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark.jl /path/to/params.json

# Fixed-bin benchmark (SearchOptC)
julia --threads=auto --project=. experiments/PartitionedFragIndex/SearchOptC/benchmark_opt_c.jl /path/to/params.json
```

### What the benchmark does

`benchmark_step_and_hints.jl`:
1. Loads data and runs calibration (~2 min)
2. Builds the partitioned index with SoA layout and 5-Da hints (~26s)
3. Runs baseline monolithic search (~70s)
4. Runs hinted SoA+SIMD search at default threshold (~11s)
5. Sweeps thresholds [1, 8, 16, 32, 64, 128, 256, 1M] to find optimal
6. Validates correctness: threshold=1 (pure binary) vs threshold=1M (pure SIMD) on 500 scans — must produce 0 mismatches

### Expected output

On a 16-thread machine with Astral DIA data (~404K MS2 scans):
- Baseline: ~70s
- Best partitioned: ~10.7s (6.5x speedup)
- Correctness: PASS (0 mismatches on 92K+ precursor IDs)
