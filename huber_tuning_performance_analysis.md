# HuberTuning Performance Analysis

## Problem Statement

KYSE70_TEST File 1 HuberTuning is taking 5-7 minutes when it should take ~30 seconds based on scan count scaling.

## Performance Comparison

### Yeast Dataset (Baseline)
- **File 1**: 11,128 scans, 13,436 PSMs
- **Time**: 14 seconds total
- **Per-scan rate**: 0.00126 seconds/scan

### KYSE70 Dataset (Problem Case)
- **File 1**: 22,702 scans, 32,263 PSMs
- **Expected time** (2x scans): 22,702 × 0.00126s = **29 seconds**
- **Actual time**: 5-7 minutes (300-420 seconds)
- **Slowdown**: **10-15x slower than expected**

This is NOT explained by 2x scan count or 2.4x PSM count.

## What Are RT Bins?

### Concept

RT bins are a spatial indexing data structure used to quickly find precursors within a retention time window. Think of them as "time buckets" that organize precursors by when they elute.

### Structure

```julia
# From CommonSearchUtils/buildRTIndex.jl
struct RetentionTimeIndex
    rt_bins::Vector{RTBin}              # Array of time bins
    precursor_ranges::Vector{UnitRange}  # Precursor indices per bin
    sorted_precursors::Vector{UInt32}    # Precursors sorted by m/z within bins
end

# Each bin covers a time range (default 0.1 min)
struct RTBin
    lb::Float32  # Lower bound (minutes)
    ub::Float32  # Upper bound (minutes)
end
```

### Example

For a 90-minute gradient with 0.1 min bins:
```
Bin 1:  [0.0  - 0.1  min] → Precursors 1-50
Bin 2:  [0.1  - 0.2  min] → Precursors 51-120
Bin 3:  [0.2  - 0.3  min] → Precursors 121-180
...
Bin 900: [89.9 - 90.0 min] → Precursors 45000-45100
```

### How They're Used in HuberTuning

**Location**: `HuberTuningSearch/utils.jl` lines 230-233

```julia
# Calculate RT window for current scan
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...)
irt_stop = searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...)
```

**Process**:
1. Convert scan's observed RT to iRT using RT model
2. Add/subtract RT tolerance to get window: `[irt - tol, irt + tol]`
3. Binary search to find bin range covering the window
4. Select all precursors in those bins for deconvolution

### Why More Bins = More Overhead

#### Dataset Comparison

**Yeast** (5 min gradient):
- Total RT range: 1.7 - 6.6 min = 4.9 min
- Number of bins: 4.9 / 0.1 = **~50 bins**
- RT tolerance: 0.85 min = **~9 bins per window**

**KYSE70** (90 min gradient):
- Total RT range: 7.5 - 89 min = 81.5 min
- Number of bins: 81.5 / 0.1 = **~815 bins**
- RT tolerance: 0.12 min = **~1.2 bins per window**

#### Overhead Sources

**1. Binary Search Depth**
```
log₂(50) = 5.6 comparisons  (yeast)
log₂(815) = 9.7 comparisons (KYSE70)
```
- 1.7x more comparisons per scan
- Minor contributor (~10% overhead)

**2. Cache Locality**
- **Yeast**: 50 bins fit in L1/L2 cache (few KB)
- **KYSE70**: 815 bins may spill to L3 or RAM
- Cache misses cause 10-100x slowdown per memory access
- Scanning through `precursor_ranges` becomes expensive

**3. Memory Access Patterns**
```julia
# For each scan, must access:
for bin_idx in irt_start:irt_stop
    prec_range = rt_index.precursor_ranges[bin_idx]  # Array lookup
    for prec_id in rt_index.sorted_precursors[prec_range]  # Another array access
        # ... select precursor
    end
end
```

With 815 bins:
- More scattered memory accesses
- Worse prefetcher performance
- Higher latency per precursor lookup

**4. Transition Selection Overhead**
Each scan calls `select_transitions_for_huber!` which:
1. Iterates through RT bins
2. Accesses precursor data for each bin
3. Copies fragment data into search structures

More bins = more iterations = more overhead, even if final precursor count is similar.

## Theories for 10-15x Slowdown

### Theory 1: Design Matrix Size Explosion ⚠️

**Hypothesis**: Longer gradients cause more precursors to co-elute, creating much larger deconvolution problems.

**Mechanism**:
```julia
# In process_delta_values! around line 398-410
buildDesignMatrix!(Hs, ion_matches, ion_misses, nmatches, nmisses, id_to_col)
```

The design matrix `Hs` is:
- **Rows**: Unique observed peaks
- **Columns**: Precursors present in the scan
- **Density**: Depends on how many fragments each precursor has

**Yeast (5 min gradient)**:
- Narrow peaks (good chromatography)
- Less co-elution
- Typical matrix: 200 rows × 10-20 cols = 2,000-4,000 entries

**KYSE70 (90 min gradient)**:
- Potentially broader peaks
- More complex samples (more proteins)
- More chance of co-elution
- Typical matrix: 300 rows × 30-50 cols = 9,000-15,000 entries

**Impact on solveHuber!**:
- Huber solver complexity scales with matrix size
- Larger matrices = more iterations to converge
- Could easily cause 5-10x slowdown per scan

**Evidence needed**: Log matrix dimensions in `process_delta_values!`

### Theory 2: Convergence Issues in solveHuber! ⚠️

**Hypothesis**: The Huber deconvolution solver takes many more iterations to converge on KYSE70 data.

**Code location**: `utils.jl` lines 425-443
```julia
solveHuber!(
    Hs, residuals, weights, δ,
    params.lambda,
    params.max_iter_newton,      # Default: 100
    params.max_iter_bisection,   # Default: 100
    params.max_iter_outer,       # Default: 100
    params.accuracy_newton,      # Default: 1e-5
    params.accuracy_bisection,   # Default: 1e-4
    params.max_diff,             # Default: 1e-5
    NoNorm()
)
```

**Potential causes**:
1. **Ill-conditioned matrices**: More co-eluting precursors = harder to resolve
2. **Poor initial guesses**: Weights initialized from previous delta values
3. **Numerical instability**: Large range of intensity values
4. **Delta value sensitivity**: Some δ values cause slow convergence

**Impact**:
- If average iterations doubles: 2x slower
- If some scans hit max iterations (100): 10x+ slower
- Testing 26 delta values × slow convergence = multiplicative effect

**Evidence needed**:
- Log actual iteration counts from solveHuber!
- Check if hitting max_iter limits
- Profile which delta values are slow

### Theory 3: RT Index Overhead (Explained Above) 🔍

**Hypothesis**: The 16x more RT bins cause significant overhead in precursor lookup and selection.

**Direct overhead**: ~1.7x from binary search depth (minor)

**Indirect overhead**: Cache misses and memory access patterns (potentially major)

**Evidence from code**:
- Lines 230-247: RT bin lookup happens once per scan
- Lines 244-247: Transition selection iterates through bin ranges
- If this section takes 10ms instead of 1ms → 9ms × 22,702 scans = 204s overhead

**Likelihood**: Moderate contributor (2-3x), not 10-15x

### Theory 4: More Precursors per Scan 🎯

**Hypothesis**: Despite tighter RT tolerance (0.12 vs 0.85 min), KYSE70 might still have more precursors per scan.

**Why?**
1. **Denser library**: KYSE70 might have more precursors in the library
2. **Better RT model**: Tighter tolerance but more precursors pass filters
3. **Sample complexity**: More proteins expressed → more co-elution

**Calculation**:
- Yeast: 11,128 scans, 13,436 PSMs → 1.2 PSMs/scan
- KYSE70: 22,702 scans, 32,263 PSMs → 1.4 PSMs/scan

**But**: The PSMs counted are AFTER q-value filtering. During HuberTuning, we consider ALL precursors in the RT window, not just high-confidence ones.

**Evidence needed**: Log actual precursor count from `scan_to_prec` dict size per scan

## Proposed Solutions

### Solution 1: Limit to Top N PSMs (Immediate Mitigation)

**Implementation**: `HuberTuningSearch/HuberTuningSearch.jl` around line 187

```julia
# Add to HuberTuningSearchParameters
max_psms_for_huber::Int  # Default: 10000

# After line 186-187
best_psms = get_best_psms(search_context, params.q_value_threshold)
file_psms = filter(row -> row.ms_file_idx == ms_file_idx, best_psms)

# NEW: Limit PSMs if too many
if nrow(file_psms) > params.max_psms_for_huber
    sort!(file_psms, :best_prob, rev=true)  # Sort by probability
    n_original = nrow(file_psms)
    file_psms = first(file_psms, params.max_psms_for_huber)
    @user_info "HuberTuning: Limited to top $(params.max_psms_for_huber) PSMs (from $n_original) for File $ms_file_idx"
end
```

**Benefits**:
- Quick to implement (5 lines of code)
- Limits worst-case performance
- Top PSMs are most informative for delta optimization
- User can adjust threshold via parameters

**Tradeoffs**:
- Might miss optimal delta if excluded PSMs have different characteristics
- Less comprehensive grid search

**Recommended limit**: 10,000-15,000 PSMs

### Solution 2: Add Detailed Profiling (Diagnostic)

**Purpose**: Identify exact bottleneck to guide optimization

**Locations to instrument**:

**A. Per-scan timing** (`utils.jl` around line 223):
```julia
scan_time = @elapsed begin
    # ... existing scan processing code ...
end
if scan_time > 0.05  # Log slow scans (>50ms)
    @user_info "Slow scan $scan_idx: $(round(scan_time, digits=3))s"
end
```

**B. Transition selection timing** (around line 244):
```julia
transition_time = @elapsed begin
    ion_idx, _ = select_transitions_for_huber!(...)
end
```

**C. Peak matching timing** (around line 245):
```julia
matching_time = @elapsed begin
    nmatches, nmisses = match_peaks_for_huber!(...)
end
```

**D. Delta processing timing** (around line 266):
```julia
delta_processing_time = @elapsed begin
    process_delta_values!(...)
end
```

**E. Inside process_delta_values!, per-delta timing**:
```julia
for δ in delta_grid
    huber_time = @elapsed begin
        solveHuber!(...)
    end
    # Track slow deltas
end
```

**Output**: Create timing breakdown showing where time is spent per scan

### Solution 3: Optimize RT Index Access (If Needed)

If profiling shows RT bin overhead is significant:

**Option A**: Cache bin lookups for consecutive scans in similar RT ranges

**Option B**: Use wider bins (0.2 or 0.5 min) for long gradients to reduce bin count

**Option C**: Switch to spatial hash instead of sorted array for very long gradients

### Solution 4: Huber Solver Optimization (If Needed)

If profiling shows convergence issues:

**Option A**: Tighten accuracy requirements after initial delta screen:
- Loose tolerance (1e-3) for broad search
- Tight tolerance (1e-5) for final optimization

**Option B**: Better initial guesses using previous delta results

**Option C**: Early termination if weight changes are minimal

## Recommended Approach

### Phase 1: Quick Fix (Now)
1. **Add PSM limit** (10K default) to prevent runaway cases
2. **Keep running** KYSE70 with limit to see if it helps

### Phase 2: Root Cause Analysis (Next)
1. **Add profiling** to identify bottleneck
2. **Compare timings** between yeast (fast) and KYSE70 (slow)
3. **Analyze profiles** to see if it's:
   - RT index overhead
   - Design matrix size
   - Solver convergence
   - Something else

### Phase 3: Targeted Optimization (Future)
Based on profiling results, implement specific optimizations for the identified bottleneck.

## Open Questions

1. What's the actual time for KYSE70 File 1 HuberTuning currently?
2. How does this compare to develop branch?
3. Is the slowdown consistent across all 8 files or just File 1?
4. What's the typical design matrix size (rows × cols)?
5. What's the average iteration count for solveHuber!?
6. Are we hitting max_iter limits frequently?
