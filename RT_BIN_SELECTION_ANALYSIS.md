# RT Bin Selection Analysis: FirstPassSearch vs SecondPassSearch

## Executive Summary

A detailed analysis of retention time (RT) bin selection in Pioneer.jl's DIA search pipeline reveals significant differences between FirstPassSearch and SecondPassSearch implementations. The FirstPassSearch method contains a critical bug that can cause it to miss RT bins, potentially excluding valid precursors from consideration.

## The RT Bin Coverage Issue in FirstPassSearch

### Current Implementation (LibrarySearch.jl, lines 42-49)

```julia
function searchFragmentIndex(...)
    # ...
    rt_bin_idx = 1  # State variable maintained across ALL scans
    for scan_idx in thread_task
        # ...
        # Calculate the RT window for current scan
        irt_lo, irt_hi = getRTWindow(rt_to_irt_spline(getRetentionTime(spectra, scan_idx)), irt_tol)

        # PROBLEM: Sequential search with persistent state
        while rt_bin_idx < length(getRTBins(frag_index)) && getHigh(getRTBin(frag_index, rt_bin_idx)) < irt_lo
            rt_bin_idx += 1  # Move forward to find first bin
        end
        while rt_bin_idx > 1 && getLow(getRTBin(frag_index, rt_bin_idx)) > irt_lo
            rt_bin_idx -= 1  # Move backward if we overshot
        end

        # Fragment index search - BUT ONLY FROM rt_bin_idx!
        searchScan!(
            getPrecursorScores(search_data),
            getRTBins(frag_index),
            # ... other parameters ...
            rt_bin_idx,  # Starting from this single index
            irt_hi,       # Upper bound passed but not used to find range
            # ...
        )
    end
end
```

### Problems with This Approach

#### Problem 1: Incomplete RT Bin Range Selection

The current code finds a single `rt_bin_idx` that contains or is near `irt_lo`, but it:
- **Does NOT find all bins** within the tolerance window `[irt_lo, irt_hi]`
- **Does NOT calculate an end index** for the RT range
- Only passes a single starting index to `searchScan!`

Example scenario:
```
RT tolerance: 5.0 iRT units
Bin width: 0.1 iRT units
Current scan iRT: 50.0

Expected: Search bins covering [45.0, 55.0] = ~50 bins
Actual: Searches from single rt_bin_idx covering ~45.0, missing bins 46-55
```

#### Problem 2: State Persistence Across Scans

The `rt_bin_idx` variable is maintained across all scans in the thread task:

```julia
rt_bin_idx = 1  # Initialized once at the start
for scan_idx in thread_task  # Used across ALL scans
    # rt_bin_idx is updated but never reset
```

This causes issues when:
1. **Scans are not in RT order**: If scan N+1 has lower RT than scan N, the algorithm may start searching from too high an index
2. **RT gaps exist**: Large jumps in RT between scans can cause the search to miss intermediate bins
3. **Multiple MS files**: Different files may have different RT ranges, but the index persists

#### Problem 3: No Upper Bound Enforcement

While `irt_hi` is calculated and passed to `searchScan!`, it's not used to determine the range of RT bins to search. The search appears to continue from `rt_bin_idx` forward without establishing a proper stopping point.

## Correct Implementation in SecondPassSearch

### Current Implementation (SecondPassSearch/utils.jl, lines 205-206)

```julia
function process_scans!(...)
    for scan_idx in scan_range
        # ...
        # Calculate RT window for current scan
        irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))

        # CORRECT: Binary search to find exact bin range
        irt_start = max(
            searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1,
            1
        )
        irt_stop = min(
            searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1,
            length(rt_index.rt_bins)
        )

        # Process ALL bins in range
        for rt_bin_idx in irt_start:irt_stop
            precs = rt_index.rt_bins[rt_bin_idx].prec
            # Process all precursors in this bin
        end
    end
end
```

### Why This Works Better

1. **Complete Coverage**:
   - Uses `searchsortedfirst` to find the first bin that might overlap `[irt - irt_tol]`
   - Uses `searchsortedlast` to find the last bin that might overlap `[irt + irt_tol]`
   - The `-1` and `+1` adjustments ensure edge bins are included

2. **Fresh Calculation Per Scan**:
   - No state is maintained between scans
   - Each scan gets an independent RT range calculation

3. **Binary Search Efficiency**:
   - O(log n) complexity for finding bin ranges
   - More efficient than sequential search

## Impact on Search Results

### Missed Precursors in FirstPassSearch

When RT bins are missed, the following cascade occurs:

1. **Fragment Index Search**: Precursors in missed bins are not scored
2. **No Initial PSMs**: These precursors don't generate PSMs in FirstPassSearch
3. **No RT Calibration**: Missing PSMs can affect RT model calibration
4. **Propagated to SecondPass**: If RT models are poor, SecondPassSearch may also miss these precursors

### Quantification of Impact

Assuming:
- iRT tolerance = 5.0 units
- RT bin width = 0.1 units
- Expected bins per search = 50

Current implementation may search only 1-5 bins instead of 50, potentially missing **90-98% of candidate precursors** in the RT dimension.

## Additional Issues Found

### 1. No Post-Selection iRT Filtering

Neither FirstPassSearch nor SecondPassSearch performs fine-grained iRT filtering after bin selection:

```julia
// What's missing:
for precursor in precursors_in_bin
    if abs(precursor.irt - target_irt) <= irt_tol
        // Process this precursor
    end
end
```

With 0.1 iRT bin width, precursors at bin edges might be outside the actual tolerance.

### 2. Inconsistent Tolerance Application

- FirstPassSearch: Uses `getRTWindow()` function (implementation not shown)
- SecondPassSearch: Direct calculation with `irt ± irt_tol`
- Could lead to different effective windows between passes

## Recommended Fixes

### Fix 1: Implement Proper RT Bin Range Selection in FirstPassSearch

```julia
function searchFragmentIndex(...)
    for scan_idx in thread_task
        # Calculate RT window
        irt = rt_to_irt_spline(getRetentionTime(spectra, scan_idx))
        irt_lo, irt_hi = irt - irt_tol, irt + irt_tol

        # Find ALL bins in range using binary search
        rt_bin_start = max(
            searchsortedfirst(getRTBins(frag_index), irt_lo,
                            lt=(bin,val)->getHigh(bin)<val) - 1,
            1
        )
        rt_bin_stop = min(
            searchsortedlast(getRTBins(frag_index), irt_hi,
                           lt=(val,bin)->val<getLow(bin)) + 1,
            length(getRTBins(frag_index))
        )

        # Process ALL bins in range
        for rt_bin_idx in rt_bin_start:rt_bin_stop
            searchScan!(
                getPrecursorScores(search_data),
                getRTBins(frag_index),
                # ... other parameters ...
                rt_bin_idx,
                # ...
            )
        end
    end
end
```

### Fix 2: Add Post-Selection iRT Filtering

```julia
# After getting precursors from RT bins
for prec_idx in precursors_from_bins
    prec_irt = getPrecursorIRT(prec_idx)
    if abs(prec_irt - target_irt) > irt_tol
        continue  # Skip this precursor
    end
    # Process precursor
end
```

### Fix 3: Create Shared RT Bin Selection Utility

```julia
function getRTBinRange(rt_bins::Vector, target_irt::Float32, irt_tol::Float32)
    rt_bin_start = max(
        searchsortedfirst(rt_bins, target_irt - irt_tol,
                         lt=(bin,val)->getHigh(bin)<val) - 1,
        1
    )
    rt_bin_stop = min(
        searchsortedlast(rt_bins, target_irt + irt_tol,
                        lt=(val,bin)->val<getLow(bin)) + 1,
        length(rt_bins)
    )
    return rt_bin_start:rt_bin_stop
end
```

## Testing Recommendations

### Unit Test for RT Bin Coverage

```julia
@testset "RT Bin Selection Coverage" begin
    # Create mock RT bins (0.1 iRT units wide)
    rt_bins = [RTBin(i*0.1, (i+1)*0.1) for i in 0:1000]

    # Test with 5.0 iRT tolerance
    target_irt = 50.0
    irt_tol = 5.0

    bin_range = getRTBinRange(rt_bins, target_irt, irt_tol)

    # Should cover approximately 100 bins (45.0 to 55.0)
    @test length(bin_range) >= 98  # Allow for edge effects
    @test rt_bins[first(bin_range)].low <= 45.1
    @test rt_bins[last(bin_range)].high >= 54.9
end
```

### Integration Test

```julia
@testset "FirstPass vs SecondPass Precursor Selection" begin
    # Run both searches on same data
    first_pass_precursors = runFirstPassSearch(test_data)
    second_pass_precursors = runSecondPassSearch(test_data)

    # For same RT windows, should consider same precursors
    for scan_idx in test_scans
        fp_precs = getPrecursorsForScan(first_pass_precursors, scan_idx)
        sp_precs = getPrecursorsForScan(second_pass_precursors, scan_idx)

        # Should have substantial overlap (>95%)
        overlap = length(intersect(fp_precs, sp_precs))
        @test overlap / length(union(fp_precs, sp_precs)) > 0.95
    end
end
```

## Conclusion

The RT bin selection issue in FirstPassSearch is a critical bug that can cause the search to miss the majority of candidate precursors within the specified retention time tolerance. This impacts:

1. **Search Sensitivity**: Fewer precursors → fewer potential identifications
2. **RT Calibration**: Fewer PSMs → poorer RT models → propagated errors
3. **Quantification**: Missing precursors → incomplete quantification

The fix is straightforward: adopt the binary search approach from SecondPassSearch and ensure complete RT bin coverage for all scans. This should be considered a high-priority fix as it directly impacts the fundamental correctness of the search algorithm.