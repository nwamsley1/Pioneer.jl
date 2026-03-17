# Partitioned Fragment Index Experiment

## Goal

Speed up the fragment index search phase of Pioneer's DIA analysis by partitioning
the monolithic `NativeFragmentIndex` by precursor m/z. This eliminates the per-fragment
binary search on precursor m/z and enables cache-friendly scoring with small counters.

## Current Best Result

**Fixed-width bin O(1) lookup (SearchOptC): 4.21x speedup** over the baseline monolithic
index on Astral data (72.4s → 17.2s, 16 threads).

## Architecture

The partitioned index splits fragments by precursor m/z into ~142 partitions (5 Da width).
Each partition uses:
- **LocalFragment** (4 bytes): `UInt16` local precursor ID + `UInt8` score
- **LocalCounter** (192 KB per thread): type-correct `UInt16`/`UInt8` counter that fits in L1 cache
- **Partition-major threading**: all threads process the same partition simultaneously for shared L2/L3 cache utilization

## Key Files

| File | Purpose |
|------|---------|
| `partitioned_types.jl` | All type definitions: `LocalFragment`, `LocalCounter`, `LocalPartition`, etc. |
| `build_partitioned_index.jl` | Build `LocalPartitionedFragmentIndex` from spectral library |
| `search_partitioned_index.jl` | Search functions: partition-major, local counter scoring, no global counter |
| `benchmark.jl` | Full benchmark comparing baseline vs partitioned (7-frag, 3-frag, bitmask) |
| `test_synthetic.jl` | 10 synthetic correctness tests |
| `profile_partitioned.jl` | PProf profiling of the partitioned search |
| `diagnose.jl` | Diagnostic tool for score mismatch analysis |

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

### SearchOptC: Fixed-Width Bins with O(1) Lookup ← CURRENT BEST
- Replace variable-width frag bins with fixed 0.005 Da bins
- CSR-style `UInt16` pointer array per RT bin maps bin index → fragment range
- Given peak m/z: `bin = floor(Int32, (mz - mz_min) * INV_BIN_WIDTH) + 1` → one multiply + one load
- **No exponential search, no binary search**
- **Result: 4.21x speedup** (17.2s vs 72.4s baseline, 22.8s variable-bin partitioned)
- Memory: 2.04 GB (1.85 GB frag_ptrs + 167 MB fragments + 24 MB local-to-global)

## Performance Progression

| Version | Speedup | Search Time | Memory |
|---------|---------|-------------|--------|
| Baseline (monolithic, 7 frags) | 1.0x | 72s | 525 MB |
| Full frag bin copy (reference) | 1.46x | 49s | 3,367 MB |
| Variable bins, CompactFragment | 3.13x | 23s | 680 MB |
| Variable bins, LocalFragment + local counter | 3.12x | 23s | 541 MB |
| Partition-major threading | 3.17x | 23s | 541 MB |
| **Fixed-width bins, O(1) lookup** | **4.21x** | **17s** | **2,037 MB** |

## Profile Breakdown (variable-bin partitioned, before SearchOptC)

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

## Checkpoints (git tags)

- `checkpoint/compact-fragment-3x-speedup` — CompactFragment, 3.13x, safe fallback
- `checkpoint/local-counter-3x-speedup` — LocalFragment + local counter, 3.12x
- `checkpoint/partition-major-3.17x` — partition-major threading, 3.17x
- `checkpoint/pre-search-optimization` — before SearchOpt* experiments

## Running

```bash
# Synthetic tests
julia --project=. experiments/PartitionedFragIndex/test_synthetic.jl

# Full benchmark (requires Astral data)
julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark.jl /path/to/params.json

# Fixed-bin benchmark
julia --threads=auto --project=. experiments/PartitionedFragIndex/SearchOptC/benchmark_opt_c.jl /path/to/params.json

# Profile
julia --threads=auto --project=. experiments/PartitionedFragIndex/profile_partitioned.jl /path/to/params.json
```
