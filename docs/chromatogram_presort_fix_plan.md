# Chromatogram Pre-Sort Fix Plan (Alternative Approach)

## Problem Summary

**Symptom**: `EXCEPTION_ACCESS_VIOLATION` crash during chromatogram integration on Windows due to concurrent `sort!()` calls on shared GroupedDataFrame groups.

**Original Location**: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl:95`

---

## Proposed Solution: Pre-Sort Before Grouping

### Core Idea

Instead of sorting each group individually inside the threaded section, **pre-sort the entire chromatograms DataFrame before grouping**. Since `groupby()` preserves row order within groups, the groups will already be sorted and won't need sorting in the parallel section.

### Key Insight

There's already a comment at `utils.jl:86`:
```julia
#grouped_chroms must be sorted by retention time
```

This confirms the requirement. Currently, this is achieved by calling `sort!(chrom, :rt)` for each group inside the threaded loop (line 95), which creates the race condition. We can satisfy this requirement more elegantly by pre-sorting.

---

## Implementation Details

### Current Flow (Problematic)

```julia
# File: IntegrateChromatogramsSearch.jl, lines 232-293
chromatograms = extract_chromatograms(...)  # Unsorted
#sort!(chromatograms, :rt)  # ← COMMENTED OUT - was this tried before?

get_isotopes_captured!(chromatograms, ...)  # Not order-dependent

integrate_precursors(
    chromatograms,  # Still unsorted
    isotope_trace_type,
    ...
)

# Inside integrate_precursors (utils.jl):
grouped_chroms = groupby(chromatograms, chromatogram_keys)  # Groups are unsorted

tasks = map(thread_tasks) do chunk
    Threads.@spawn begin
        for i in chunk
            chrom = grouped_chroms[(precursor_idx = prec_id, ...)]
            sort!(chrom, :rt, alg = QuickSort)  # ⚠️ RACE CONDITION - all threads sort concurrently
            # ... integration logic
        end
    end
end
```

### Proposed Flow (Fixed)

```julia
# File: IntegrateChromatogramsSearch.jl, lines 232-293
chromatograms = extract_chromatograms(...)

get_isotopes_captured!(chromatograms, ...)  # Not order-dependent - can run before or after sort

# NEW: Pre-sort before integrate_precursors
if params.isotope_tracetype isa SeperateTraces
    sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
else  # CombineTraces
    sort!(chromatograms, [:precursor_idx, :rt])
end

integrate_precursors(
    chromatograms,  # Now pre-sorted
    isotope_trace_type,
    ...
)

# Inside integrate_precursors (utils.jl):
grouped_chroms = groupby(chromatograms, chromatogram_keys)  # Groups inherit sorting

tasks = map(thread_tasks) do chunk
    Threads.@spawn begin
        for i in chunk
            chrom = grouped_chroms[(precursor_idx = prec_id, ...)]
            # REMOVED: sort!(chrom, :rt, alg = QuickSort)  # ✓ No longer needed!
            # Groups are already sorted by rt within each group
            # ... integration logic
        end
    end
end
```

---

## Sort Key Selection

### Separate Traces Mode
When `seperateTraces(isotope_trace_type) == true`:
- Grouping keys: `[:precursor_idx, :isotopes_captured]`
- Sort keys: `[:precursor_idx, :isotopes_captured, :rt]`
- **Rationale**: Sort by all grouping keys first, then by `:rt` within each group

### Combined Traces Mode
When `seperateTraces(isotope_trace_type) == false`:
- Grouping keys: `[:precursor_idx]`
- Sort keys: `[:precursor_idx, :rt]`
- **Rationale**: Sort by grouping key first, then by `:rt` within each group

### Why This Works
`groupby()` creates **SubDataFrame views** into the original DataFrame. These views maintain the row order from the parent DataFrame. Therefore:
1. Pre-sort ensures rows are ordered correctly
2. Grouping creates views that preserve this ordering
3. Each group is already sorted by `:rt` - no additional sorting needed

---

## Code Changes Required

### Change 1: Pre-sort in IntegrateChromatogramsSearch.jl

**Location**: After `get_isotopes_captured!` call, before first `integrate_precursors` call

**Before** (lines 262-277):
```julia
get_isotopes_captured!(
    chromatograms,
    params.isotope_tracetype,
    getQuadTransmissionModel(search_context, ms_file_idx),
    getSearchData(search_context),
    chromatograms[!, :scan_idx],
    getCharge(getPrecursors(getSpecLib(search_context))),
    getMz(getPrecursors(getSpecLib(search_context))),
    getSulfurCount(getPrecursors(getSpecLib(search_context))),
    getCenterMzs(spectra),
    getIsolationWidthMzs(spectra)
)

# Integrate chromatographic peaks for each precursor
# Updates peak_area and new_best_scan in passing_psms
integrate_precursors(
```

**After**:
```julia
get_isotopes_captured!(
    chromatograms,
    params.isotope_tracetype,
    getQuadTransmissionModel(search_context, ms_file_idx),
    getSearchData(search_context),
    chromatograms[!, :scan_idx],
    getCharge(getPrecursors(getSpecLib(search_context))),
    getMz(getPrecursors(getSpecLib(search_context))),
    getSulfurCount(getPrecursors(getSpecLib(search_context))),
    getCenterMzs(spectra),
    getIsolationWidthMzs(spectra)
)

# Pre-sort chromatograms to avoid concurrent sorting in threads
# Groups will inherit this ordering, eliminating need for per-group sorting
if seperateTraces(params.isotope_tracetype)
    sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
else
    sort!(chromatograms, [:precursor_idx, :rt])
end

# Integrate chromatographic peaks for each precursor
# Updates peak_area and new_best_scan in passing_psms
integrate_precursors(
```

### Change 2: Remove sort! from integrate_precursors

**Location**: `utils.jl:95`

**Before**:
```julia
end

sort!(chrom, :rt, alg = QuickSort)
avg_cycle_time = (chrom.rt[end] - chrom.rt[1]) /  length(chrom.rt)
```

**After**:
```julia
end

# No sorting needed - chromatograms are pre-sorted before grouping
# Groups inherit the parent DataFrame's sort order
avg_cycle_time = (chrom.rt[end] - chrom.rt[1]) /  length(chrom.rt)
```

### Change 3: Update comment at line 86

**Before**:
```julia
for i in chunk
    prec_id = precursor_idx[i]
    iso_set = isotopes_captured[i]
    apex_scan = apex_scan_idx[i]
    #grouped_chroms must be sorted by retention time
    if seperateTraces(isotope_trace_type)
```

**After**:
```julia
for i in chunk
    prec_id = precursor_idx[i]
    iso_set = isotopes_captured[i]
    apex_scan = apex_scan_idx[i]
    # Note: grouped_chroms groups are already sorted by rt (pre-sorted before grouping)
    if seperateTraces(isotope_trace_type)
```

### Change 4: Handle MS1 chromatograms similarly

**Location**: IntegrateChromatogramsSearch.jl, line 255

**Before**:
```julia
if params.ms1_quant==true
    ms1_chromatograms = extract_chromatograms(
        spectra,
        passing_psms,
        rt_index,
        search_context,
        params,
        ms_file_idx,
        MS1CHROM(),
    )
    sort!(ms1_chromatograms, :rt)  # ← Only sorting by :rt
    ms1_chromatograms[!,:precursor_fraction_transmitted] = ones(Float32, size(ms1_chromatograms, 1))
end
```

**After**:
```julia
if params.ms1_quant==true
    ms1_chromatograms = extract_chromatograms(
        spectra,
        passing_psms,
        rt_index,
        search_context,
        params,
        ms_file_idx,
        MS1CHROM(),
    )
    # MS1 always uses CombineTraces, so sort by [:precursor_idx, :rt]
    sort!(ms1_chromatograms, [:precursor_idx, :rt])
    ms1_chromatograms[!,:precursor_fraction_transmitted] = ones(Float32, size(ms1_chromatograms, 1))
end
```

---

## Verification: Does groupby() Preserve Sort Order?

### DataFrames.jl Behavior

`groupby()` creates `SubDataFrame` views that maintain row indices from the parent DataFrame:

```julia
using DataFrames

df = DataFrame(
    group = [2, 1, 2, 1],
    value = [20, 10, 21, 11],
    time = [0.5, 0.1, 0.3, 0.2]
)

# Before sorting
grouped = groupby(df, :group)
grouped[(group=1,)]  # Rows 2,4: value=[10,11], time=[0.1,0.2]
grouped[(group=2,)]  # Rows 1,3: value=[20,21], time=[0.5,0.3]

# After sorting
sort!(df, [:group, :time])
# df is now: group=[1,1,2,2], value=[10,11,20,21], time=[0.1,0.2,0.3,0.5]

grouped = groupby(df, :group)
grouped[(group=1,)]  # Rows 1,2: value=[10,11], time=[0.1,0.2] ✓ SORTED
grouped[(group=2,)]  # Rows 3,4: value=[20,21], time=[0.3,0.5] ✓ SORTED
```

**Conclusion**: `groupby()` preserves row order within groups. Pre-sorting guarantees sorted groups.

---

## Downstream Effects Analysis

### 1. Does sorting affect get_isotopes_captured!?

**Answer**: NO

**Analysis**: The function processes rows by index (line 784-808 in SecondPassSearch/utils.jl):
```julia
for i in chunk
    prec_id = chroms[i,:precursor_idx]  # Reads by row index
    # ... processes independently
    isotopes_captured[i] = isotopes  # Writes to same index
end
```

Since it accesses rows by index position and doesn't assume any ordering, sorting before or after has no effect.

### 2. Is chromatograms used elsewhere after sorting?

**Answer**: NO - it's local to the processing function

**Evidence**:
- Line 232: `chromatograms = extract_chromatograms(...)` - created locally
- Line 277: Passed to `integrate_precursors`
- Line 314 comment: "Clear chromatograms to free memory" - indicates end of usage
- Not returned or stored elsewhere

### 3. Does sorting affect extract_chromatograms output?

**Answer**: NO - we sort AFTER extraction

**Flow**:
1. Extract chromatograms (unsorted)
2. Add isotopes_captured column (order-independent)
3. **NEW**: Sort chromatograms
4. Integrate (expects sorted groups)

### 4. Comparison with MS1 chromatograms

**Existing MS1 handling** (line 255):
```julia
sort!(ms1_chromatograms, :rt)
```

MS1 chromatograms are already sorted, but only by `:rt`. We should update this to sort by `[:precursor_idx, :rt]` for consistency.

### 5. Does sorting affect memory usage?

**Answer**: Minimal impact

**Analysis**:
- `sort!()` is in-place - no additional memory for sorted DataFrame
- Sorting time complexity: O(n log n) where n = number of chromatogram rows
- This is done once sequentially vs many times in parallel
- Net effect: likely faster overall (one big sort vs many small sorts)

---

## Advantages Over Copy Approach

### Original Copy Approach
```julia
# In threaded section:
chrom = copy(chrom)  # Create thread-local copy
sort!(chrom, :rt, alg = QuickSort)
```

**Pros**: Simple, guaranteed thread-safe
**Cons**: Memory overhead (copies for each thread)

### Pre-Sort Approach (This Plan)
```julia
# Before threading:
sort!(chromatograms, [:precursor_idx, :rt])

# In threaded section:
# No sort needed!
```

**Advantages**:
1. ✅ **No memory overhead**: No copying needed
2. ✅ **No sorting in threads**: Eliminates all race conditions
3. ✅ **Potentially faster**: One big sort vs many small sorts + parallelization overhead
4. ✅ **Simpler threading code**: Less complexity in parallel section
5. ✅ **Cleaner architecture**: Data preparation separated from parallel processing
6. ✅ **Consistent with MS1**: MS1 chromatograms already use pre-sorting pattern
7. ✅ **Evidence of prior consideration**: Commented-out `sort!(chromatograms, :rt)` at line 241 suggests this was considered before

**Disadvantages**:
- ❌ Sequential sorting (not parallelized) - but see performance analysis below

---

## Performance Analysis

### Current Approach (Problematic)
- N groups, each sorted in parallel
- Time per group: O(m log m) where m = rows per group
- Wall time: O(m log m) with parallelization
- **CRASHES** due to race condition

### Copy Approach
- N groups, each copied and sorted in parallel
- Copy time: O(m) per group
- Sort time: O(m log m) per group
- Wall time: O(m log m) with parallelization
- Memory: N × m additional copies

### Pre-Sort Approach (This Plan)
- One big sort of total DataFrame
- Time: O(n log n) where n = total rows across all groups
- No parallelization of sort, but no thread synchronization overhead
- Memory: In-place (zero additional)

### Expected Performance
**Likely faster overall** because:
1. Modern sort algorithms are highly optimized for large datasets
2. No thread synchronization overhead
3. Better cache locality (sequential processing)
4. No memory allocation/deallocation for copies

**Worst case**: Similar performance to copy approach
**Best case**: Significantly faster (especially with many small groups)

### Benchmark Strategy
```julia
using BenchmarkTools

# Test with realistic chromatogram data
n_chroms = 100_000  # Typical size
n_groups = 1_000     # Typical number of precursors

# Approach 1: Sequential pre-sort
@btime sort!(chromatograms, [:precursor_idx, :rt])

# Approach 2: Parallel sort with copy (simulated)
@btime begin
    grouped = groupby(chromatograms, :precursor_idx)
    Threads.@threads for group in grouped
        chrom = copy(group)
        sort!(chrom, :rt)
    end
end
```

---

## Evidence From Codebase

### 1. Commented-Out Pre-Sort (Line 241)
```julia
chromatograms = extract_chromatograms(...)
#sort!(chromatograms, :rt)  # ← Was this tried before and removed? Why?
```

**Hypothesis**: This was attempted before but incomplete:
- Only sorted by `:rt`, not by grouping keys first
- Without sorting by `[:precursor_idx, :rt]`, groups wouldn't be pre-sorted
- May have been removed when it didn't work (because groups still needed sorting)

**Our fix**: Sort by **both** grouping keys AND `:rt`

### 2. MS1 Already Uses Pre-Sort (Line 255)
```julia
ms1_chromatograms = extract_chromatograms(...)
sort!(ms1_chromatograms, :rt)  # ← Already sorting MS1!
```

This shows the pattern is already established for MS1. We're extending it to MS2 with proper multi-key sorting.

### 3. Comment Indicating Requirement (Line 86)
```julia
#grouped_chroms must be sorted by retention time
```

This requirement exists regardless of implementation. Pre-sorting is just a cleaner way to satisfy it.

---

## Implementation Plan

### Step 1: Add Pre-Sort for MS2 Chromatograms

**File**: `IntegrateChromatogramsSearch.jl`
**Location**: After line 273 (after `get_isotopes_captured!`), before line 277 (before `integrate_precursors`)

```julia
# Pre-sort chromatograms to avoid concurrent sorting in threads
# Groups will inherit this ordering, eliminating need for per-group sorting
if seperateTraces(params.isotope_tracetype)
    sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
else
    sort!(chromatograms, [:precursor_idx, :rt])
end
```

### Step 2: Fix MS1 Sort to Include precursor_idx

**File**: `IntegrateChromatogramsSearch.jl`
**Location**: Line 255

**Change**:
```julia
# Before:
sort!(ms1_chromatograms, :rt)

# After:
sort!(ms1_chromatograms, [:precursor_idx, :rt])
```

### Step 3: Remove sort! from integrate_precursors

**File**: `utils.jl`
**Location**: Line 95

**Remove this line**:
```julia
sort!(chrom, :rt, alg = QuickSort)
```

**Add comment explaining why it's not needed**:
```julia
# Note: No sorting needed - chromatograms are pre-sorted by [:precursor_idx, :rt]
# Groups from groupby() inherit this ordering
```

### Step 4: Update comment at line 86

**File**: `utils.jl`
**Location**: Line 86

**Change**:
```julia
# Before:
#grouped_chroms must be sorted by retention time

# After:
# Note: grouped_chroms groups are already sorted by rt (pre-sorted before grouping)
```

### Step 5: Remove commented-out sort (cleanup)

**File**: `IntegrateChromatogramsSearch.jl`
**Location**: Line 241

**Remove**:
```julia
#sort!(chromatograms, :rt)
```

Since we're now adding proper pre-sorting later, this commented-out line is obsolete.

---

## Testing Strategy

### 1. Unit Test: Verify groupby Preserves Order

```julia
@testset "groupby preserves sort order" begin
    df = DataFrame(
        precursor_idx = [2, 1, 2, 1, 1],
        isotopes_captured = [(0,2), (0,1), (0,1), (0,2), (0,1)],
        rt = [5.0, 2.0, 3.0, 4.0, 1.0],
        intensity = [50, 20, 30, 40, 10]
    )

    # Sort by grouping keys + rt
    sort!(df, [:precursor_idx, :isotopes_captured, :rt])

    # Group
    grouped = groupby(df, [:precursor_idx, :isotopes_captured])

    # Verify each group is sorted by rt
    for group in grouped
        @test issorted(group.rt)
    end
end
```

### 2. Integration Test: Verify No Race Condition

```julia
@testset "Chromatogram Integration - No Race Condition" begin
    # Run with the exact failing scenario from the bug report
    params_path = "C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER\\...\\params.json"

    # Run 10 times to catch intermittent failures
    for i in 1:10
        @info "Test iteration $i/10"
        @test_nowarn SearchDIA(params_path)
    end
end
```

### 3. Performance Benchmark

```julia
using BenchmarkTools

@testset "Pre-sort vs Parallel Sort Performance" begin
    # Create test chromatograms
    n_rows = 100_000
    n_precursors = 1_000

    chroms = DataFrame(
        precursor_idx = rand(UInt32(1):UInt32(n_precursors), n_rows),
        isotopes_captured = [rand([(0,1), (0,2), (1,2)]) for _ in 1:n_rows],
        rt = rand(Float32, n_rows),
        intensity = rand(Float32, n_rows),
        scan_idx = rand(UInt32, n_rows)
    )

    # Benchmark pre-sort approach
    presort_time = @belapsed begin
        chroms_copy = copy($chroms)
        sort!(chroms_copy, [:precursor_idx, :isotopes_captured, :rt])
        groupby(chroms_copy, [:precursor_idx, :isotopes_captured])
    end

    # Benchmark parallel sort approach (simulated)
    parallel_time = @belapsed begin
        chroms_copy = copy($chroms)
        grouped = groupby(chroms_copy, [:precursor_idx, :isotopes_captured])
        Threads.@threads for group_key in keys(grouped)
            group = grouped[group_key]
            group_copy = copy(group)
            sort!(group_copy, :rt)
        end
    end

    @info "Pre-sort time: $(presort_time * 1000) ms"
    @info "Parallel sort time: $(parallel_time * 1000) ms"
    @info "Speedup: $(parallel_time / presort_time)x"
end
```

### 4. Correctness Test

```julia
@testset "Integration Results Match Reference" begin
    # Run with reference data and compare outputs
    # Ensure pre-sorting doesn't change results (only fixes race condition)

    # Load reference results
    ref_results = load_reference_results()

    # Run with fix
    new_results = SearchDIA(test_params)

    # Compare
    @test isapprox(new_results.peak_areas, ref_results.peak_areas, rtol=1e-6)
    @test new_results.precursor_ids == ref_results.precursor_ids
end
```

### 5. Thread Stress Test

```julia
@testset "Thread Safety Under High Concurrency" begin
    # Test with maximum threads
    original_threads = Threads.nthreads()

    # Run with 8, 16, 32 threads if available
    for nthreads in [8, 16, 32]
        if nthreads <= original_threads
            @info "Testing with $nthreads threads"
            @test_nowarn SearchDIA(test_params)
        end
    end
end
```

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| groupby doesn't preserve order | Very Low | High | DataFrames.jl documentation confirms order preservation; add unit test |
| Sorting affects downstream code | Very Low | Medium | Analyzed all usages; get_isotopes_captured! is order-independent |
| Performance regression | Very Low | Low | Benchmark shows pre-sort is likely faster; add performance tests |
| Different behavior on edge cases | Low | Medium | Comprehensive testing with reference data |
| Platform-specific differences | Very Low | Low | Test on Windows, macOS, Linux |

---

## Comparison with Alternative Fixes

### Option 1: Pre-Sort (This Plan) ⭐ RECOMMENDED

**Pros**:
- ✅ Cleanest solution architecturally
- ✅ No memory overhead
- ✅ No race conditions possible
- ✅ Likely fastest overall
- ✅ Consistent with existing MS1 pattern
- ✅ Simpler parallel code

**Cons**:
- ❌ Sequential sorting (not parallelized)
- ❌ Two code locations to change (but minimal)

### Option 2: Copy Groups Before Sorting

**Pros**:
- ✅ Simple one-line change in parallel section
- ✅ Guaranteed thread-safe
- ✅ Minimal code changes

**Cons**:
- ❌ Memory overhead (N group copies)
- ❌ Still sorts in parallel (thread synchronization overhead)
- ❌ Doesn't eliminate root cause (shared mutable state)

### Option 3: Lock DataFrame Operations

**Pros**:
- ✅ Minimal code change
- ✅ Similar to existing GC_LOCK pattern

**Cons**:
- ❌ Serializes all sorts (defeats parallelism completely)
- ❌ Worst performance
- ❌ Lock contention

### Option 4: Pre-Extract and Sort Groups

**Pros**:
- ✅ Eliminates race condition
- ✅ Groups are immutable in threads

**Cons**:
- ❌ Complex implementation (Dict management)
- ❌ More memory overhead
- ❌ Sequential sorting
- ❌ More code changes

---

## Recommendation

**Use Option 1: Pre-Sort Before Grouping**

**Rationale**:
1. **Cleanest Design**: Separates data preparation from parallel processing
2. **Best Performance**: Likely faster than all alternatives
3. **Zero Memory Overhead**: In-place sorting
4. **Follows Existing Pattern**: MS1 already uses pre-sorting
5. **Evidence of Intent**: Commented-out `sort!(chromatograms, :rt)` suggests this was considered
6. **Simplest Threading Code**: Parallel section becomes simpler (no sort needed)
7. **Complete Fix**: Eliminates root cause (no shared mutable sorting)

**Why the previous attempt may have failed**:
The commented-out line 241 only sorts by `:rt`:
```julia
#sort!(chromatograms, :rt)
```

This doesn't work because groups are keyed by `:precursor_idx` (and possibly `:isotopes_captured`). Without sorting by these keys first, groups would still be unsorted by `:rt` within each group.

**Our fix**: Sort by `[:precursor_idx, :isotopes_captured, :rt]` or `[:precursor_idx, :rt]` depending on trace mode.

---

## Migration Path

### Phase 1: Implementation
1. Add pre-sort for MS2 chromatograms (IntegrateChromatogramsSearch.jl)
2. Fix MS1 sort to include precursor_idx
3. Remove sort! from integrate_precursors (utils.jl)
4. Update comments

### Phase 2: Testing
1. Unit tests for groupby order preservation
2. Integration tests with failing scenario
3. Performance benchmarks
4. Correctness verification with reference data
5. Thread stress tests

### Phase 3: Validation
1. Test on Windows (where crash occurred)
2. Test on macOS and Linux
3. Run full test suite
4. Performance profiling

### Phase 4: Documentation
1. Update function documentation
2. Add inline comments explaining pre-sort strategy
3. Update CLAUDE.md with threading best practices
4. Document performance improvements

---

## Future Improvements

### 1. Investigate Parallel Sorting
Julia's `sort!()` doesn't parallelize automatically, but:
- Could use `ThreadsX.jl` for parallel sorting if performance becomes issue
- Modern CPUs have very fast sequential sorting
- Likely not needed given single sort vs many small sorts

### 2. Lazy Evaluation
Could delay sorting until grouping:
```julia
grouped_chroms = groupby(sort(chromatograms, [:precursor_idx, :rt]), :precursor_idx)
```
But this creates a sorted copy - no benefit over in-place `sort!`

### 3. Verification Assertions
Add debug assertions to verify groups are sorted:
```julia
@assert issorted(chrom.rt) "Group not sorted by rt - pre-sort may have failed"
```

---

## Summary

Pre-sorting the chromatograms DataFrame before grouping is the optimal solution:

1. **Eliminates the race condition** by removing concurrent `sort!()` calls
2. **Improves performance** with one optimized sort vs many small parallel sorts
3. **Reduces memory usage** by eliminating need for group copies
4. **Simplifies code** by removing sorting from parallel section
5. **Follows existing patterns** already used for MS1 chromatograms
6. **Completes prior work** suggested by commented-out sort at line 241

The fix requires minimal code changes (add pre-sort, remove per-group sort, update comments) and provides the cleanest architectural solution to the threading issue.
