# Fix for EXCEPTION_ACCESS_VIOLATION in Chromatogram Integration

## Problem Summary

Multiple `EXCEPTION_ACCESS_VIOLATION` crashes occurring during parallel chromatogram integration at ~17% progress (line 1 of error log shows "16.7%"). Analysis of the crash dumps from `/Users/n.t.wamsley/Desktop/Nov-11-bug.txt` reveals the issue.

## Key Evidence

**From crash dumps:**
```
Exception: EXCEPTION_ACCESS_VIOLATION at 0x236d6ea98fc
jl_gc_alloc_ at C:/workdir/src\gc-stock.c:797
jl_gc_small_alloc_inner at C:/workdir/src\gc-stock.c:735

Allocations: 998818349335 (Pool: 998817920628; Big: 428707); GC: 34190
```

**Critical observations:**
1. ✅ Crash happens **during GC allocation** (`jl_gc_alloc_`), not during Arrow array access
2. ✅ **~1 trillion allocations** before crash (998,818,349,335)
3. ✅ **34,190 GC runs** - garbage collector running almost constantly
4. ✅ Stack traces show crashes in multiple threads simultaneously
5. ✅ Crashes at `chainedvector.jl:100`, `getMzArray`, `getIntensityArray` - but always during GC

## Root Cause: Allocation Bomb 💣

**NOT an Arrow thread-safety issue** - it's a memory allocation catastrophe!

### Allocation Hotspots Found

**Location:** `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`

#### Hotspot 1: Chromatogram Vector Growth (Line 362-364)
```julia
# CURRENT CODE (BAD):
for scan_idx in scan_range
    # ... process scan ...
    if nmatches > 2
        for j in 1:prec_temp_size
            rt_idx += 1
            if rt_idx + 1 > length(chromatograms)
                append!(chromatograms, Vector{MS2ChromObject}(undef, 500000))  # ❌ DISASTER!
            end
        end
    end
end
```

**Why this is catastrophic:**
- Allocates **500,000 new elements** every time it needs to grow
- Julia's `append!` creates a NEW vector to hold old + new data (copy semantics)
- With 24 threads, this happens constantly across all threads
- Old vectors become garbage immediately
- Example: If this triggers 100 times per thread → 24 × 100 × 500k = **1.2 billion objects allocated**

#### Hotspot 2: List Comprehension in Hot Loop (Line 332)
```julia
# CURRENT CODE (BAD):
if getIdToCol(search_data).size > length(weights)
    new_entries = getIdToCol(search_data).size - length(weights) + 1000
    resize!(weights, length(weights) + new_entries)
    resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
    append!(getUnscoredPsms(search_data),
        [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])  # ❌ Allocates intermediate array!
end
```

**Why this is bad:**
- List comprehension `[... for _ in 1:new_entries]` allocates temporary array
- Happens potentially every scan in the hot loop
- Multiplied across 24 threads

#### Hotspot 3: Multiple Resizes Per Scan (Lines 328-333)
```julia
# CURRENT CODE:
if getIdToCol(search_data).size > length(weights)
    new_entries = getIdToCol(search_data).size - length(weights) + 1000
    resize!(weights, length(weights) + new_entries)                      # Resize 1
    resize!(getSpectralScores(search_data), ...)                         # Resize 2
    append!(getUnscoredPsms(search_data), ...)                           # Resize 3
end
```

**Impact:**
- 3 separate allocations per resize event
- Triggered frequently in the scan loop
- Across 24 threads = allocation storm

### The Crash Sequence

1. **24 threads** all running `build_chromatograms` simultaneously
2. Each thread allocating **massive vectors** (500k elements) repeatedly
3. **GC runs constantly** (34,190 times!) trying to free memory
4. During GC, Julia needs to allocate memory for GC bookkeeping
5. **GC allocation fails** because system is out of memory / address space
6. Crash occurs in `jl_gc_alloc_` with `EXCEPTION_ACCESS_VIOLATION`

**Note:** The Arrow/ChainedVector references in the stack trace are **red herrings** - they just happen to be what was executing when GC tried (and failed) to allocate.

## Recommended Solution

### Phase 1: Fix Allocation Hotspots ⭐ (CRITICAL)

#### Fix 1: Exponential Growth for Chromatograms (Line 362-364)

```julia
# BEFORE (BAD):
chromatograms = Vector{MS2ChromObject}(undef, 500000)  # Line 232
# ... later in loop ...
if rt_idx + 1 > length(chromatograms)
    append!(chromatograms, Vector{MS2ChromObject}(undef, 500000))  # ❌
end

# AFTER (GOOD):
# Pre-allocate with better estimate
estimated_points = length(scan_range) * 100  # ~100 points per scan average
chromatograms = Vector{MS2ChromObject}(undef, max(estimated_points, 10000))
# ... later in loop ...
if rt_idx + 1 > length(chromatograms)
    resize!(chromatograms, length(chromatograms) * 2)  # ✅ Exponential growth, not fixed chunks
end
```

**Expected impact:** Reduces allocations from O(N²) to O(log N) where N = total points

#### Fix 2: Remove List Comprehension (Line 332)

```julia
# BEFORE (BAD):
append!(getUnscoredPsms(search_data),
    [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])  # ❌

# AFTER (GOOD):
psms = getUnscoredPsms(search_data)
old_length = length(psms)
resize!(psms, old_length + new_entries)  # ✅ Direct resize, no intermediate allocation
for i in (old_length + 1):length(psms)
    psms[i] = eltype(psms)()
end
```

**Expected impact:** Eliminates temporary array allocation in hot loop

#### Fix 3: Pre-allocate with Generous Capacity (Line 232, 328-333)

```julia
# BEFORE (BAD):
chromatograms = Vector{MS2ChromObject}(undef, 500000)
# weights and other vectors: allocated on-demand

# AFTER (GOOD):
# Pre-allocate all working vectors with generous initial capacity
estimated_points = length(scan_range) * 100
chromatograms = Vector{MS2ChromObject}(undef, max(estimated_points, 10000))

# Also pre-allocate other frequently-resized vectors in search_data
# (weights, spectral_scores, unscored_psms) with capacity of 100k+
```

### Phase 2: Verify Fix

**Test metrics to monitor:**
```julia
# Add at start of build_chromatograms:
alloc_start = Base.gc_num().allocd

# Add at end:
alloc_end = Base.gc_num().allocd
alloc_diff = alloc_end - alloc_start
println("Thread $(thread_id) allocated: $(alloc_diff / 1e9) GB")
```

**Success criteria:**
- ✅ Allocations drop from ~1 trillion to < 10 billion (100x reduction minimum)
- ✅ GC runs drop from 34k to < 500
- ✅ No crashes during integration
- ✅ Similar or better performance

### Phase 3: Additional Optimizations (If Needed)

If Phase 1 fixes aren't sufficient, consider:

#### Option A: Disable GC During Critical Section
```julia
function build_chromatograms(...)
    GC.enable(false)
    try
        # ... main processing loop ...
    finally
        GC.enable(true)
    end
end
```

**Pros:** Prevents GC from running during allocation-heavy section
**Cons:** Only safe if Phase 1 allocations are fixed first; risk of OOM

#### Option B: Add Arrow Access Locking (Unlikely to be needed)
```julia
# In src/structs/MassSpecData.jl
const ARROW_ACCESS_LOCK = ReentrantLock()

function getMzArray(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer) where T
    lock(ARROW_ACCESS_LOCK) do
        getMzArrays(ms_data)[scan_idx]
    end
end
# ... similar for other accessors
```

**Note:** This is unlikely to be the actual issue, but included as a fallback if GC fixes don't fully resolve crashes.

## Files to Modify

1. **Primary fix:**
   - `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`
     - Lines 232 (initial allocation)
     - Lines 328-333 (resize logic)
     - Lines 362-364 (chromatogram growth)

2. **Secondary (if needed):**
   - `src/structs/MassSpecData.jl` (Arrow locking, only if Phase 1 insufficient)

## Testing Plan

### Minimal Test
```julia
# Run on the failing dataset with allocation tracking
julia --threads=24 --gcthreads=1,1

using Pioneer
# Enable allocation tracking
Profile.Allocs.@profile sample_rate=0.01 begin
    SearchDIA("path/to/failing/params.json")
end
```

### Validation Checks
1. ✅ Monitor allocation count (should be < 10B vs ~1T before)
2. ✅ Monitor GC runs (should be < 500 vs 34k before)
3. ✅ Verify no crashes during chromatogram integration
4. ✅ Check memory usage stays reasonable (< system RAM)
5. ✅ Verify output correctness (same PSMs and quantification as expected)

## Expected Outcome

**Before fixes:**
- Allocations: ~1,000,000,000,000 (1 trillion)
- GC runs: 34,190
- Crash: Yes, at ~17% progress

**After fixes:**
- Allocations: ~10,000,000,000 (10 billion) - **100x reduction**
- GC runs: < 500 - **70x reduction**
- Crash: No
- Performance: Similar or improved (less GC overhead)

## Why Other Theories Were Considered (But Wrong)

### ❌ Theory 1: Arrow Thread-Safety
**Seemed plausible because:**
- Stack traces showed `chainedvector.jl:100`, `getMzArray`, `getIntensityArray`
- Multiple threads accessing shared `spectra` object
- Arrow's ChainedVector could have mutable internal state

**But evidence against:**
- Scan ranges don't overlap (round-robin distribution)
- Crashes happen in `jl_gc_alloc_`, not during Arrow indexing
- 1 trillion allocations before crash points to allocation issue, not thread safety

**Verdict:** Red herring. Arrow access just happened to be on the call stack when GC failed.

### ❌ Theory 3: Thread ID Mismatch
**Considered:** Maybe threads are sharing search_data incorrectly

**Evidence against:**
- Each thread gets unique `thread_id` from partition
- `search_data = getSearchData(search_context)[thread_id]` is thread-local
- No indication of data corruption, just crashes

### ❌ Theory 4: Buffer Overflow
**Considered:** Maybe writing past end of `chromatograms` array

**Evidence against:**
- Bounds checking would catch this (different error)
- Crash is in GC allocation, not array write
- Allocation count is the smoking gun

## Confidence Level

**Phase 1 fixes (allocation hotspots): 95% confidence** ✅

This is almost certainly the root cause based on:
1. Crash location (`jl_gc_alloc_`)
2. Allocation count (1 trillion!)
3. GC pressure (34k runs)
4. Specific allocation hotspots found in code
5. Timing (crashes at 17% when allocations accumulate)

**Phase 2 fixes (Arrow locking): 5% confidence**

Only include this if Phase 1 doesn't fully resolve the issue. It's a belt-and-suspenders approach.

## References

- Error log: `/Users/n.t.wamsley/Desktop/Nov-11-bug.txt`
- Main file: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`
- Existing GC locking pattern: `src/utils/writeArrow.jl:20` (for reference if needed)
- Thread partitioning: `src/Routines/SearchDIA/CommonSearchUtils/partitionThreadTasks.jl`
