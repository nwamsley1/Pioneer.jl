# Thread Safety Investigation: Segfault Theories

## Problem Summary

Intermittent `EXCEPTION_ACCESS_VIOLATION` crashes during chromatogram integration on Windows:
- **Error 1** (windows_err.txt): Crashes at 67.9% completion in `rtIndexTransitionSelection.jl:69` - **BEFORE previous fix**
- **Error 2** (windows_err2.txt): Crashes at 45.2% completion in garbage collector during spline evaluation - **AFTER previous fix**
- **IMPORTANT**: The previous fix (thread-local Sets) **SOLVED Error 1**. Error 2 is a **NEW, DIFFERENT issue**.
- Crashes occur at seemingly random points during execution
- Some searches with slightly different libraries complete successfully

## Current Focus: Error 2 (NEW ISSUE)

**Error 1 has been FIXED by the previous patch (thread-local Sets). We are now dealing with Error 2, which is a different problem.**

## Error Details

### Error 1: Set Membership Check Crash (FIXED ✅)
```
EXCEPTION_ACCESS_VIOLATION at 0x20ab7dde088
#_select_transitions_impl!#382 at rtIndexTransitionSelection.jl:69
```

**Line 69 in rtIndexTransitionSelection.jl:**
```julia
if !isnothing(precursors_passing) && prec_idx ∉ precursors_passing
```

This was a Set membership check crash caused by concurrent DataFrame access. **RESOLVED** by creating thread-local Sets.

### Error 2: Garbage Collector Crash (CURRENT ISSUE ⚠️)
```
EXCEPTION_ACCESS_VIOLATION at 0x7fffea530640
ijl_gc_mark_queue_objarray at libjulia-internal.dll (unknown line)
jl_parallel_gc_threadfun (multiple GC threads)
ijl_box_float32 (boxing Float32 values)
splevl at libraryBSpline.jl:61 [inlined]
getIntensity at LibraryIon.jl:362 [inlined]
getFragIsotopes! at fillTransitionList.jl:189
```

**Critical Stack Trace Details:**
- **Line 31 in error**: `ijl_box_float32` - Julia is boxing Float32 values (allocation!)
- **Line 35**: `splevl` - Spline evaluation (math operations causing boxing)
- **Line 36**: `getIntensity` - Reading from shared SplineFragmentLookup
- **Lines 9-22**: Multiple GC threads running `ijl_gc_mark_queue_objarray`
- **Lines 45-50**: GC marking functions (`gc_setmark_big`, `gc_mark_outrefs`)

**What's Happening:**
1. Multiple threads simultaneously access shared `SplineFragmentLookup` data
2. Spline evaluation involves Float32 math that gets boxed (triggers allocation)
3. Allocation triggers GC safepoint
4. GC runs in parallel and tries to mark shared library objects
5. GC sees corrupted pointers because multiple threads are reading the same data without synchronization
6. **SEGFAULT** in GC marking code

## Theory Ranking for Error 2 (Current Issue)

**For addressing Error 2 specifically:**
1. ⭐⭐⭐⭐⭐ Theory 2: GC + Shared Library Data Access (PRIMARY FOCUS)
2. ⭐⭐⭐⭐ Theory 3: Shared Spectral Library Structures Without Protection
3. ⭐⭐⭐ Theory 4: Windows Memory Model Issues
4. ⭐ Theory 1: Sets (No longer relevant - fixed)

---

## Theory 1: Julia Set Implementation Not Thread-Safe for Concurrent Reads (FIXED ✅)

### Status: RESOLVED - This theory explained Error 1, which has been fixed

### Confidence: VERY HIGH (for Error 1)

### Root Cause
Julia's hash-based `Set` implementation is **not guaranteed to be thread-safe for concurrent reads**, especially on Windows with weaker memory ordering.

**This issue has been resolved** by the previous patch that created thread-local Sets.

### Evidence
1. **Direct crash location**: Error 1 crashes exactly at `prec_idx ∉ precursors_passing` (Set membership check)
2. **Non-deterministic behavior**: Crashes at different progress points (45% vs 67%)
3. **Library-dependent**: Different library completes successfully (different Set sizes/hash patterns)
4. **Platform-specific**: Only observed on Windows (weaker memory model than Linux)
5. **Multiple concurrent readers**: Many threads checking Set membership simultaneously

### Technical Details

**Set Implementation Issues:**
- Hash table uses open addressing or chaining with pointer manipulation
- Cache line sharing between threads reading same/nearby hash buckets
- No memory barriers for read operations
- Hash function itself may involve shared state

**Why Concurrent Reads Can Fail:**
```julia
# Thread A reads Set during GC safepoint
if prec_idx ∉ precursors_passing  # Hash lookup starts

# Thread B or GC modifies internal state
# (even if logically read-only, GC can relocate objects)

# Thread A continues with stale pointer -> SEGFAULT
```

**Windows Memory Model:**
- x86-64 Windows has relaxed memory ordering vs Linux
- Compiler optimizations may reorder loads
- Cache coherency delays more visible under load

### Current Implementation (Lines 176-189 in utils.jl)

```julia
# Line 176: Create thread-local Sets (CORRECT)
precursor_sets = [Set(passing_psms[!, :precursor_idx]) for _ in 1:Threads.nthreads()]

# Line 180-188: Pass thread-specific Set to each task (CORRECT)
tasks = map(thread_tasks) do thread_task
    Threads.@spawn begin
        thread_id = first(thread_task)
        return build_chromatograms(
            spectra,
            last(thread_task),
            precursor_sets[thread_id],  # ← Each thread gets its own Set
            ...
        )
    end
end
```

**However, this doesn't solve the problem because:**
1. Each Set is still being read concurrently from the same thread multiple times
2. GC can run while threads are reading Sets
3. Julia Sets don't have internal synchronization even for reads

### Why Previous Fix Didn't Work

The previous fix created thread-local Sets to avoid concurrent **DataFrame access**, but:
- Sets themselves are not thread-safe for concurrent reads
- GC interactions with Set internals cause issues
- No explicit memory barriers or atomic operations

## Theory 2: Garbage Collector + Concurrent Shared Library Access ⭐⭐⭐⭐⭐

### Status: PRIMARY THEORY FOR ERROR 2

### Confidence: VERY HIGH

### Root Cause
GC runs asynchronously and tries to mark shared spectral library objects while multiple threads are actively reading the same structures. The GC sees inconsistent pointer states because the object graph is being traversed concurrently without synchronization.

### Evidence
1. **Error 2 explicit GC crash**: `ijl_gc_mark_queue_objarray`
2. **During spline evaluation**: Accessing shared fragment/spline coefficient arrays
3. **Object graph traversal**: GC marking objects that threads are simultaneously reading
4. **Multiple allocation-heavy operations**: Spline evaluation, array accesses

### Technical Details

**GC Interaction Pattern:**
```julia
# Thread performs spline evaluation
getIntensity(frag, spline_data)  # Accesses:
  → frag.intensity (NTuple - boxed floats)
  → spline_data.knots (NTuple - boxed floats)
  → lookup.nce_model[] (Ref - mutable box)

# GC safepoint triggered (any allocation)
# GC thread tries to mark these objects
# → Sees inconsistent pointers if main thread mid-access
# → EXCEPTION_ACCESS_VIOLATION
```

**Why This Happens:**
- Spline coefficients stored as `NTuple{N, Float32}` (can be boxed)
- Multiple threads reading same tuples simultaneously
- GC doesn't know these are being concurrently accessed
- No GC.@preserve blocks protecting reads

### Why This Is The Primary Issue for Error 2

Unlike Error 1 (which was about concurrent DataFrame/Set access), Error 2 is specifically about:
1. **Allocation during critical section**: `ijl_box_float32` shows Float32 boxing causing allocations
2. **GC triggered mid-operation**: Allocation triggers GC while accessing shared library
3. **Parallel GC threads**: Multiple GC threads try to mark the same shared objects
4. **No GC protection**: No `GC.@preserve` blocks protecting shared library access
5. **Stack trace clearly shows**: GC marking code crashing while user thread does spline math

## Theory 3: Shared Spectral Library Data Structures ⭐⭐⭐

### Confidence: MEDIUM-HIGH

### Root Cause
Multiple threads concurrently reading shared spectral library structures, especially mutable references.

### Evidence
1. **Mutable Ref in library**: `SplineFragmentLookup` contains `nce_model::Base.Ref{<:NceModel}`
2. **All threads read same library**: `getFragmentLookupTable(getSpecLib(search_context))`
3. **Spline data access**: Fragment vectors, coefficient arrays, precursor metadata
4. **Error 2 during library access**: Crashes during `getIntensity()` which reads library data

### Technical Details

**Shared Library Structures (from LibraryIon.jl:550-556):**
```julia
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}        # Shared across threads
    prec_frag_ranges::Vector{UInt64}              # Shared across threads
    knots::NTuple{M, T}                           # Shared across threads
    nce_model::Base.Ref{<:NceModel{T}}           # MUTABLE reference!
    degree::Int64                                 # Shared across threads
end
```

**Concurrent Access Pattern:**
```julia
# All threads call (line 274 in utils.jl):
getFragmentLookupTable(getSpecLib(search_context))

# Which returns THE SAME SplineFragmentLookup instance
# All threads then:
getFragments(lookup)              # Read frags vector
getPrecFragRange(lookup, idx)     # Read prec_frag_ranges
getSplineData(lookup, charge, mz) # Read knots, nce_model.Ref
```

**The nce_model Ref Issue:**
```julia
# Line 565-566: Setter modifies Ref
function setNceModel!(lookup::SplineFragmentLookup{N,M,T}, new_nce_model::NceModel{T})
    lookup.nce_model[] = new_nce_model  # Modifies shared mutable box
end

# Line 578: Getter dereferences Ref
getNCE(lfp, prec_charge, prec_mz) → reads lfp.nce_model[]
```

**Why This Causes Issues:**
1. `setNceModel!` called per MS file (line 188-190 in IntegrateChromatogramsSearch.jl)
2. Multiple threads reading `nce_model[]` concurrently
3. GC can relocate Ref box during concurrent access
4. No memory barriers protecting Ref access

### Shared Precursor Data

**Lines 274-278 in utils.jl:**
```julia
getFragmentLookupTable(getSpecLib(search_context)),
precs_temp,
getMz(getPrecursors(getSpecLib(search_context))),          # Shared vector
getCharge(getPrecursors(getSpecLib(search_context))),      # Shared vector
getSulfurCount(getPrecursors(getSpecLib(search_context))), # Shared vector
```

All threads read same precursor property vectors simultaneously.

## Theory 4: Windows Memory Model and Cache Coherency ⭐⭐

### Confidence: MEDIUM

### Root Cause
Windows has weaker memory ordering guarantees than Linux, and without explicit memory barriers, threads can see stale or inconsistent data.

### Evidence
1. **Platform-specific crashes**: Only observed on Windows
2. **Non-deterministic timing**: Cache coherency delays
3. **Hash table operations**: Involve cache line sharing
4. **No explicit barriers**: Code assumes sequential consistency

### Technical Details

**Memory Ordering Differences:**
```
Linux (x86-64):    Strong memory model, fewer reorderings
Windows (x86-64):  Weaker guarantees, more compiler/CPU reorderings
```

**Hash Table Cache Line Sharing:**
```
Thread 1: Reads hash bucket at address 0x1000
Thread 2: Reads hash bucket at address 0x1040
→ Both in same 64-byte cache line
→ Cache line ping-pongs between cores
→ Performance degradation + potential stale reads
```

**No Memory Barriers:**
Julia's Set implementation doesn't use atomic operations or memory barriers for read-only operations, assuming sequential consistency.

## Theory 5: Large Set Size Causing Hash Table Pathologies ⭐⭐

### Confidence: MEDIUM-LOW

### Root Cause
With 312 files and thousands of precursors per file, Sets contain 10,000+ elements. Large hash tables under concurrent access can exhibit pathological behavior.

### Evidence
1. **Large dataset**: 312 files processed
2. **Many precursors**: Thousands per thread's Set
3. **Hash collisions**: More likely with large Sets
4. **Memory pressure**: Large Sets + spline data + fragment arrays

### Technical Details

**Set Size Estimation:**
```
- ~3,350,136 passing PSMs (from error log)
- Distributed across precursors
- Each thread's Set: Potentially 5,000-15,000 unique precursor IDs
- Hash table size: ~2x elements = 10,000-30,000 buckets
- Memory per Set: ~100-300 KB
```

**Large Hash Table Issues:**
- More cache misses during lookups
- Higher collision probability
- More memory allocations during resizing
- GC pressure from large data structures

## Cross-Theory Connections

```
Set Concurrent Access (Theory 1)
    ↓ corrupts internal pointers
GC Inconsistent State (Theory 2)
    ↓ tries to mark corrupted objects
Shared Library Access (Theory 3)
    ↓ amplifies with multiple readers
Windows Memory Model (Theory 4)
    ↓ makes timing issues worse
Large Set Size (Theory 5)
```

## Diagnostic Evidence Summary

| Observation | Theory 1 | Theory 2 | Theory 3 | Theory 4 | Theory 5 |
|-------------|----------|----------|----------|----------|----------|
| Crash at Set membership check | ⭐⭐⭐⭐⭐ | ⭐ | ⭐ | ⭐⭐ | ⭐⭐ |
| GC crash during spline eval | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐ |
| Non-deterministic crashes | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐ |
| Windows-only crashes | ⭐⭐⭐ | ⭐⭐ | ⭐ | ⭐⭐⭐⭐ | ⭐ |
| Library-dependent behavior | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ | ⭐ | ⭐⭐⭐ |
| Random crash locations | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ | ⭐⭐ |

## Recommended Fixes for Error 2 (Priority Order)

### Priority 1: Add GC Protection for Critical Sections ⭐⭐⭐⭐⭐

**Problem:** GC running during concurrent reads of shared library data during spline evaluation
**Solution:** Use `GC.@preserve` blocks to prevent GC from interrupting critical operations

```julia
# Option A: Protect the spline evaluation (fillTransitionList.jl around line 189)
total_fragment_intensity = GC.@preserve frag spline_data begin
    getIntensity(frag, spline_data)
end

# Option B: Protect the entire fragment processing loop (rtIndexTransitionSelection.jl around line 113-141)
for frag_idx in precursor_fragment_range
    frag = fragment_ions[frag_idx]
    frag.rank > max_frag_rank && continue

    GC.@preserve frag spline_data isotopes precursor_transmission begin
        getFragIsotopes!(
            prec_estimation_type,
            isotopes,
            precursor_transmission,
            # ... rest of arguments
        )

        transition_idx = addTransitionIsotopes!(
            transitions,
            transition_idx,
            frag,
            isotopes,
            # ... rest of arguments
        )
    end
end
```

**Why This Works:**
- Prevents GC from running while threads access shared library data
- Protects the entire spline evaluation pipeline from GC interruption
- Minimal performance impact (GC still runs, just not during critical sections)
- Thread-safe because each thread has its own execution context

### Priority 2: Disable GC During Processing (Testing Only) ⭐⭐⭐⭐

**Problem:** Need to confirm GC is the issue
**Solution:** Temporarily disable GC during chromatogram integration to validate theory

```julia
# In extract_chromatograms function, wrap parallel tasks:
GC.enable(false)
try
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            # ... processing code
        end
    end
    result = vcat(fetch.(tasks)...)
finally
    GC.enable(true)
    GC.gc()  # Run GC after processing completes
end
```

**Important:** This is for **diagnostic purposes only**. If disabling GC fixes the crash, it confirms Theory 2.

### Priority 3: Copy Shared Library Data to Thread-Local Storage ⭐⭐⭐⭐

**Problem:** Multiple threads reading shared spectral library structures
**Solution:** Copy needed data to thread-local arrays

```julia
# Before parallel processing, per thread:
thread_local_data = [
    ThreadLocalLibraryData(
        copy(getMz(getPrecursors(spec_lib))),
        copy(getCharge(getPrecursors(spec_lib))),
        copy(getSulfurCount(getPrecursors(spec_lib)))
    ) for _ in 1:Threads.nthreads()
]
```

### Priority 4: Make nce_model Ref Access Atomic ⭐⭐⭐

**Problem:** Concurrent reads of mutable Ref
**Solution:** Use atomic reference or copy per thread

```julia
# Option A: Atomic Ref
nce_model::Threads.Atomic{NceModel{T}}

# Option B: Thread-local copy (preferred)
# Set NCE model once per thread at start of build_chromatograms
thread_local_nce = lookup.nce_model[]  # Copy once
# Use thread_local_nce instead of dereferencing Ref repeatedly
```

### Priority 5: Add Memory Barriers (Windows-Specific) ⭐⭐

**Problem:** Weak memory ordering on Windows
**Solution:** Add explicit memory barriers

```julia
# Before critical reads:
Base.Threads.atomic_fence()
value = shared_data_structure
Base.Threads.atomic_fence()
```

### Priority 6: Reduce Thread Count for Testing ⭐

**Problem:** More threads = more concurrency = more crashes
**Solution:** Test with JULIA_NUM_THREADS=1,2,4 to isolate issue

## Testing Strategy

### Phase 1: Validate Theory 1 (Set Safety)
1. Replace Sets with sorted vectors
2. Run same workload that crashed
3. If crashes stop → Theory 1 confirmed

### Phase 2: Validate Theory 2 (GC)
1. Add GC.@preserve blocks
2. Disable GC during critical section (for testing only)
3. Monitor GC stats

### Phase 3: Validate Theory 3 (Shared Library)
1. Copy library data to thread-local storage
2. Profile memory access patterns
3. Use ThreadSanitizer if available

### Phase 4: Stress Testing
1. Run with different thread counts (1, 2, 4, 8, 16)
2. Test with larger/smaller libraries
3. Monitor Windows Event Viewer for additional crash info

## Implementation Priority for Error 2

**Immediate (This PR) - Focus on GC Protection:**
1. Add GC.@preserve blocks around spline evaluation (Priority 1)
2. Test with GC disabled to confirm theory (Priority 2 - diagnostic only)

**If GC protection insufficient (Follow-up PR):**
3. Copy shared library data to thread-local storage (Priority 3)
4. Make nce_model Ref atomic or thread-local (Priority 4)

**Platform-Specific (If still seeing Windows issues):**
5. Add memory barriers for Windows (Priority 5)
6. Performance profiling and optimization

## Success Criteria

1. ✅ No EXCEPTION_ACCESS_VIOLATION crashes on Windows
2. ✅ Consistent completion across multiple runs
3. ✅ Works with various library sizes
4. ✅ Scalable to higher thread counts
5. ✅ No performance regression (< 5% slowdown acceptable)

## Additional Notes

### Why Thread-Local Sets Solved Error 1 But Not Error 2

The previous fix (thread-local Sets) **successfully resolved Error 1** by avoiding concurrent DataFrame access. However, Error 2 is a completely different issue:

**Error 1 (FIXED):**
- Concurrent access to DataFrame → Set membership check crash
- Solution: Thread-local Sets ✅

**Error 2 (CURRENT):**
- Concurrent access to shared spectral library during GC
- Different root cause: GC marking + shared library access
- Requires different solution: GC protection or thread-local library copies

### Why Some Runs Succeed

Non-deterministic crashes depend on:
- Exact thread scheduling (varies by CPU load)
- GC timing (depends on allocation patterns)
- Hash function output (varies by data content)
- Cache coherency timing (varies by CPU model)

A slightly different library might:
- Have different hash collision patterns
- Trigger GC at different times
- Distribute precursors differently across Sets
- Change memory layout just enough to avoid crashes

This is classic **heisenbug** behavior - timing-dependent race conditions.
