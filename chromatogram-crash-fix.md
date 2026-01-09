# Fix for Chromatogram Integration Crash

## Problem

`EXCEPTION_ACCESS_VIOLATION` crashes at ~17% progress during parallel chromatogram integration.

**Crash signature:**
```
Exception: EXCEPTION_ACCESS_VIOLATION at 0x236d6ea98fc
jl_gc_alloc_ at C:/workdir/src\gc-stock.c:797
Allocations: 998818349335 (Pool: 998817920628; Big: 428707); GC: 34190
```

## Root Cause

**Memory allocation catastrophe** causing GC to fail.

**Evidence:**
- ~1 trillion allocations before crash
- 34,190 GC runs (running almost constantly)
- Crash during GC allocation, not during data access

## The Fix

Three allocation hotspots in `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`:

### Fix 1: Exponential Growth for Chromatograms (Lines 232, 362-364)

**Current code (BAD):**
```julia
# Line 232
chromatograms = Vector{MS2ChromObject}(undef, 500000)

# Lines 362-364
if rt_idx + 1 > length(chromatograms)
    append!(chromatograms, Vector{MS2ChromObject}(undef, 500000))  # ❌ Allocates 500k every time!
end
```

**Fixed code:**
```julia
# Line 232 - Pre-allocate with better estimate
estimated_points = length(scan_range) * 100  # ~100 points per scan average
chromatograms = Vector{MS2ChromObject}(undef, max(estimated_points, 10000))

# Lines 362-364 - Exponential growth instead of fixed chunks
if rt_idx + 1 > length(chromatograms)
    resize!(chromatograms, length(chromatograms) * 2)  # ✅ Grows 2x, not +500k
end
```

### Fix 2: Remove List Comprehension Allocation (Line 332)

**Current code (BAD):**
```julia
append!(getUnscoredPsms(search_data),
    [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])  # ❌ Temp array!
```

**Fixed code:**
```julia
psms = getUnscoredPsms(search_data)
old_length = length(psms)
resize!(psms, old_length + new_entries)  # ✅ Direct resize
for i in (old_length + 1):length(psms)
    psms[i] = eltype(psms)()
end
```

### Fix 3: Pre-allocate Working Vectors (Lines 228-231)

**Current code:**
```julia
Hs = getHs(search_data)
weights = getTempWeights(search_data)
precursor_weights = getPrecursorWeights(search_data)
residuals = getResiduals(search_data)
```

**Ensure these are pre-allocated with sufficient capacity in SearchDataStructures initialization**
- Check their initial sizes are adequate (100k+ elements)
- This avoids frequent resizes in the hot loop at lines 328-333

## Testing

### Run with allocation tracking:
```julia
julia --threads=24 --gcthreads=1,1

# At start of build_chromatograms, add:
alloc_start = Base.gc_num().allocd

# At end, add:
alloc_end = Base.gc_num().allocd
@info "Thread allocations (GB): $((alloc_end - alloc_start) / 1e9)"
```

### Success criteria:
- ✅ Allocations: < 10 billion (vs ~1 trillion before)
- ✅ GC runs: < 500 (vs 34k before)
- ✅ No crashes during integration
- ✅ Memory usage: < system RAM

## Expected Impact

**Before:**
- Allocations: 1,000,000,000,000 (1 trillion)
- GC runs: 34,190
- Result: Crash at 17%

**After:**
- Allocations: ~10,000,000,000 (10 billion) → **100x reduction**
- GC runs: < 500 → **70x reduction**
- Result: Completes successfully

**Performance:** Similar or better (less GC overhead)

## Implementation Checklist

- [ ] Fix 1: Line 232 - Better initial allocation
- [ ] Fix 1: Lines 362-364 - Change `append!` to `resize!` with 2x growth
- [ ] Fix 2: Line 332 - Remove list comprehension, use loop
- [ ] Fix 3: Verify SearchDataStructures pre-allocations
- [ ] Add allocation tracking (temporary, for testing)
- [ ] Test on failing dataset
- [ ] Verify allocation count drops to < 10B
- [ ] Verify no crashes
- [ ] Remove allocation tracking
- [ ] Commit fix
