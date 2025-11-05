# NCE Model Design Analysis: Root Cause and Solutions

## The Core Problem

### Current Design Issue

The `SplineFragmentLookup` struct contains a **mutable reference** to an NCE model:

```julia
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    nce_model::Base.Ref{<:NceModel{T}}  # ← MUTABLE REFERENCE!
    degree::Int64
end
```

### How It's Used (The Anti-Pattern)

**Setting NCE Model (once per MS file, BEFORE parallel processing):**
```julia
# IntegrateChromatogramsSearch.jl:188-191
setNceModel!(
    getFragmentLookupTable(getSpecLib(search_context)),  # ← Shared library
    getNceModelModel(search_context, ms_file_idx)        # ← File-specific NCE
)

# LibraryIon.jl:565-567
function setNceModel!(lookup::SplineFragmentLookup{N,M,T}, new_nce_model::NceModel{T})
    lookup.nce_model[] = new_nce_model  # ← Modifies shared mutable box
end
```

**Reading NCE Model (during parallel processing, many threads simultaneously):**
```julia
# rtIndexTransitionSelection.jl:121
getSplineData(lookup, prec_charge, prec_mz)

# LibraryIon.jl:575-580
function getSplineData(lfp::SplineFragmentLookup{N,M,T}, prec_charge::UInt8, prec_mz::T)
    return SplineType(
        getKnots(lfp),
        getNCE(lfp, prec_charge, prec_mz),  # ← Calls getNCE
        getDegree(lfp)
    )
end

# LibraryIon.jl:591-593
function getNCE(lfp::SplineFragmentLookup, prec_charge::UInt8, prec_mz::T)
    return lfp.nce_model[](prec_mz, prec_charge)  # ← DEREFERENCES MUTABLE REF!
end
```

**The Critical Path to Crash:**
```julia
# fillTransitionList.jl:189
getFragIsotopes!(...)
  → getIntensity(frag, spline_data)
    → splevl(getNCE(spline_data), ...)  # Reads NCE from SplineType
      → Float32 math operations
        → Boxing allocations  # ijl_box_float32
          → GC triggered
            → GC tries to mark lfp.nce_model Ref
              → Sees inconsistent pointer state
                → SEGFAULT in ijl_gc_mark_queue_objarray
```

## Why This Is Wrong

### 1. Mutable Shared State

**The Design Smell:**
```julia
Global Shared Resource: SplineFragmentLookup (one instance for entire library)
    ↓
Mutable Field: nce_model::Base.Ref
    ↓
Modified Per MS File: setNceModel!(lookup, file_specific_model)
    ↓
Read By All Threads: Multiple threads call getNCE(lookup, ...)
    ↓
NO SYNCHRONIZATION: No locks, atomics, or GC protection
```

**Why `Base.Ref` Is Problematic:**
- `Base.Ref{T}` is a mutable heap-allocated box
- Dereferencing `ref[]` involves pointer dereference
- GC can relocate the box or its contents
- No atomic guarantees for concurrent access
- GC marking threads see the Ref while user threads read it

### 2. Semantic Mismatch

**NCE Model Varies By:**
- MS File (different instruments/settings)
- Potentially search method (FirstPass vs SecondPass)

**SplineFragmentLookup Represents:**
- Global spectral library (universal across all files)
- Immutable fragment data
- Spline coefficients
- Knot points

**The Mismatch:**
```
File-Specific Data (NCE Model)
    stored in
Global Invariant Structure (Spectral Library)
```

This violates the **Single Responsibility Principle** and creates a **temporal coupling** where the library's behavior depends on when you set the NCE model.

### 3. Concurrency Anti-Pattern

**Classic Race Condition:**
```
Time  Thread 1                    Thread 2                    GC Thread
----  -------------------------  -------------------------  ----------------------
t0    Process file A             Process file A
t1    Read nce_model[]           Read nce_model[]
t2    Call NCE model             Call NCE model
t3    splevl() allocates         splevl() allocates
t4    GC triggered               GC triggered               Mark nce_model Ref
t5                                                           See corrupted pointer
t6                                                           SEGFAULT
```

**No Protection:**
- No mutex/lock around Ref access
- No atomic operations
- No GC.@preserve blocks
- No memory barriers
- Assumes benign concurrent reads (WRONG!)

## Why It Seemed to Work Before

### False Sense of Security

1. **Low Probability Race:**
   - GC triggers are timing-dependent
   - Specific allocation patterns required
   - Depends on thread scheduling
   - May work 90%+ of the time

2. **Platform Differences:**
   - Linux: Stronger memory ordering, less visible
   - Windows: Weaker memory model, more crashes
   - macOS: Intermediate behavior

3. **Workload Dependent:**
   - Small libraries: Fewer allocations, less GC
   - Fewer threads: Less contention
   - Different libraries: Different hash patterns, different crash probability

## Proposed Solutions (Ranked)

### Solution 1: Pass NCE Model as Explicit Parameter ⭐⭐⭐⭐⭐

**The Clean Architecture Fix**

**Rationale:** NCE model is file-specific context, not library invariant data. It should be passed as a parameter, not stored in shared state.

**Changes Required:**

#### Step 1: Add NCE parameter to getSplineData
```julia
# OLD (LibraryIon.jl:575-580):
function getSplineData(lfp::SplineFragmentLookup{N,M,T}, prec_charge::UInt8, prec_mz::T)
    return SplineType(
        getKnots(lfp),
        getNCE(lfp, prec_charge, prec_mz),  # ← Gets from mutable Ref
        getDegree(lfp)
    )
end

# NEW:
function getSplineData(
    lfp::SplineFragmentLookup{N,M,T},
    nce_model::NceModel{T},  # ← Explicit parameter!
    prec_charge::UInt8,
    prec_mz::T
) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        nce_model(prec_mz, prec_charge),  # ← Direct call, no Ref deref
        getDegree(lfp)
    )
end
```

#### Step 2: Remove mutable Ref from SplineFragmentLookup
```julia
# OLD:
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    nce_model::Base.Ref{<:NceModel{T}}  # ← REMOVE THIS
    degree::Int64
end

# NEW:
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    # No nce_model field!
    degree::Int64
end
```

#### Step 3: Thread transition selection to pass NCE model
```julia
# OLD (rtIndexTransitionSelection.jl:116-121):
transition_idx = fillTransitionListPrecomputed!(
    transitions,
    prec_estimation_type,
    getPrecFragRange(lookup, prec_idx),
    getFragments(lookup),
    getSplineData(lookup, prec_charge, prec_mz),  # ← Gets NCE from lookup
    ...
)

# NEW:
transition_idx = fillTransitionListPrecomputed!(
    transitions,
    prec_estimation_type,
    getPrecFragRange(lookup, prec_idx),
    getFragments(lookup),
    getSplineData(lookup, nce_model, prec_charge, prec_mz),  # ← Explicit NCE
    ...
)
```

#### Step 4: Store NCE model in SearchData (thread-local)
```julia
# SearchData already has thread-local storage
# Add nce_model field to SearchDataStructures or pass as parameter

# In IntegrateChromatogramsSearch.jl:
nce_model = getNceModelModel(search_context, ms_file_idx)  # Get once per file

# Pass to build_chromatograms:
build_chromatograms(
    spectra,
    scan_range,
    precursors_passing,
    rt_index,
    search_context,
    search_data,
    params,
    ms_file_idx,
    nce_model,  # ← NEW parameter
    chrom_type
)
```

#### Step 5: Thread NCE through selectTransitions! call chain
```julia
# selectTransitions! signature change:
function selectTransitions!(
    transitions::Vector{DetailedFrag{Float32}},
    strategy::RTIndexedTransitionSelection,
    prec_estimation_type::PrecEstimation,
    lookup::LibraryFragmentLookup,
    nce_model::NceModel,  # ← NEW parameter
    precs_temp::Vector{UInt32},
    # ... rest of parameters
)
```

**Advantages:**
- ✅ Eliminates mutable shared state
- ✅ Makes dependencies explicit
- ✅ Thread-safe by design
- ✅ No performance overhead
- ✅ Clearer semantics (NCE is context, not library property)
- ✅ No GC issues (no shared mutable Ref)
- ✅ Type-stable
- ✅ Testable (easy to pass different NCE models)

**Disadvantages:**
- ⚠️ Requires threading parameter through multiple call sites
- ⚠️ Signature changes propagate through codebase
- ⚠️ Need to update all selectTransitions! implementations
- ⚠️ Need to update all getSplineData callers

**Estimated Changes:**
- ~15-20 function signatures updated
- ~30-40 call sites updated
- ~200-300 lines of code changes
- Low risk (compile-time type checking)

### Solution 2: Thread-Local NCE Model Cache ⭐⭐⭐⭐

**The Pragmatic Middle Ground**

**Rationale:** Keep NCE model in library for convenience, but cache it in thread-local storage to avoid repeated Ref dereferencing.

**Implementation:**

#### Add NCE model field to SearchDataStructures
```julia
# In SearchTypes.jl or similar:
mutable struct SearchDataStructures{T<:AbstractFloat}
    # ... existing fields ...
    nce_model::Union{Nothing, NceModel{T}}  # ← Thread-local cache
end
```

#### Cache NCE model at start of file processing
```julia
# In build_chromatograms (utils.jl):
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    # Cache NCE model once per thread at start
    if search_data.nce_model === nothing
        search_data.nce_model = lookup.nce_model[]  # One Ref deref at start
    end

    # Use cached model throughout processing
    cached_nce = search_data.nce_model

    # ... processing loop ...
    for scan_idx in scan_range
        # Pass cached_nce instead of dereferencing Ref repeatedly
        spline_data = SplineType(
            getKnots(lookup),
            cached_nce(prec_mz, prec_charge),  # No Ref deref!
            getDegree(lookup)
        )
    end
end
```

**Advantages:**
- ✅ Minimal code changes
- ✅ Reduces Ref dereferences (better performance)
- ✅ Thread-local cache (no sharing)
- ✅ GC only sees cached copy (less risk)

**Disadvantages:**
- ⚠️ Still has mutable Ref in library (design smell remains)
- ⚠️ Doesn't eliminate root cause
- ⚠️ Cache invalidation complexity
- ⚠️ Initial Ref deref still has small risk window

**Estimated Changes:**
- ~5-10 function changes
- ~50-100 lines of code
- Medium risk (cache invalidation bugs)

### Solution 3: Immutable Wrapper Struct ⭐⭐⭐

**The Composition Pattern**

**Rationale:** Create immutable wrapper that combines lookup + NCE model, pass wrapper instead of lookup.

**Implementation:**

```julia
# New immutable wrapper
struct SplineFragmentLookupWithContext{N,M,T<:AbstractFloat}
    lookup::SplineFragmentLookup{N,M,T}
    nce_model::NceModel{T}  # ← Immutable field!
end

# Delegation methods
getFragments(wrapper::SplineFragmentLookupWithContext) = getFragments(wrapper.lookup)
getPrecFragRange(wrapper::SplineFragmentLookupWithContext, idx) = getPrecFragRange(wrapper.lookup, idx)
getKnots(wrapper::SplineFragmentLookupWithContext) = getKnots(wrapper.lookup)
getDegree(wrapper::SplineFragmentLookupWithContext) = getDegree(wrapper.lookup)

# getSplineData now accesses immutable field
function getSplineData(wrapper::SplineFragmentLookupWithContext{N,M,T}, prec_charge, prec_mz)
    return SplineType(
        getKnots(wrapper),
        wrapper.nce_model(prec_mz, prec_charge),  # ← Immutable access!
        getDegree(wrapper)
    )
end

# Create wrapper per MS file
function process_file(search_context, ms_file_idx)
    base_lookup = getFragmentLookupTable(getSpecLib(search_context))
    nce_model = getNceModelModel(search_context, ms_file_idx)

    # Create immutable wrapper
    lookup_with_context = SplineFragmentLookupWithContext(
        base_lookup,
        nce_model
    )

    # Use wrapper in processing
    build_chromatograms(..., lookup_with_context, ...)
end
```

**Advantages:**
- ✅ Immutable by design (thread-safe)
- ✅ Clear composition pattern
- ✅ No mutable shared state
- ✅ Type-safe

**Disadvantages:**
- ⚠️ Need to change type parameters throughout call chain
- ⚠️ Delegation boilerplate
- ⚠️ More complex type system
- ⚠️ Still keeps Ref in base lookup (confusing)

**Estimated Changes:**
- ~20-30 type signature changes
- ~100-150 lines of code
- Medium risk (type parameter propagation)

### Solution 4: GC.@preserve Workaround ⭐⭐

**The Band-Aid**

**Rationale:** Just protect the Ref dereference from GC, don't fix the design.

```julia
function getNCE(lfp::SplineFragmentLookup, prec_charge::UInt8, prec_mz::T)
    GC.@preserve lfp begin
        return lfp.nce_model[](prec_mz, prec_charge)
    end
end
```

**Advantages:**
- ✅ Minimal code change (1-2 lines)
- ✅ Quick fix

**Disadvantages:**
- ❌ Doesn't fix root cause
- ❌ Design smell remains
- ❌ Performance overhead (GC barriers)
- ❌ May not fully solve the issue
- ❌ Doesn't address conceptual problem

## Recommendation

### Immediate (This PR): Solution 4 + Solution 2

**Short Term:** Add GC protection and thread-local caching to stop the crashes:
```julia
# Quick fix - protect Ref deref
function getNCE(lfp::SplineFragmentLookup, prec_charge::UInt8, prec_mz::T)
    GC.@preserve lfp begin
        return lfp.nce_model[](prec_mz, prec_charge)
    end
end

# Better fix - cache in thread-local storage
# Add nce_model cache to SearchDataStructures
# Copy once at start of file processing
```

### Follow-Up PR: Solution 1

**Long Term:** Refactor to pass NCE model as explicit parameter:
- Remove `nce_model` from SplineFragmentLookup
- Add `nce_model` parameter to `getSplineData`
- Thread through call chain
- Store in SearchDataStructures or pass as parameter

**Why This Order:**
1. **Safety First**: GC protection stops crashes immediately
2. **Incremental Improvement**: Thread-local cache improves performance and reduces risk
3. **Proper Fix**: Parameter passing fixes the design properly
4. **Testable Steps**: Each step can be tested independently

## Testing Strategy

### Validate Each Solution

**Test 1: Concurrent Access**
```julia
@testset "Concurrent NCE Access" begin
    # Spawn multiple threads reading same lookup
    # Verify no crashes
    # Monitor GC behavior
end
```

**Test 2: Different NCE Models Per File**
```julia
@testset "File-Specific NCE Models" begin
    # Process multiple files with different NCE models
    # Verify correct NCE used for each file
    # No cross-contamination
end
```

**Test 3: Stress Test**
```julia
@testset "High Concurrency Stress" begin
    # Run with maximum threads
    # Large library
    # Force frequent GC
    # Verify stability
end
```

## Success Criteria

1. ✅ No EXCEPTION_ACCESS_VIOLATION crashes
2. ✅ Correct NCE model used for each MS file
3. ✅ Thread-safe concurrent access
4. ✅ No performance regression (<5% acceptable)
5. ✅ Clean separation of concerns (NCE as context, not library state)

## Lessons Learned

### Design Principles Violated

1. **Mutable Shared State**: Global mutable state is concurrency anti-pattern
2. **Hidden Dependencies**: NCE model dependency not visible in function signatures
3. **Temporal Coupling**: Library behavior depends on when setNceModel! was called
4. **Single Responsibility**: Library manages both fragment data AND execution context
5. **Lack of Synchronization**: Concurrent access without protection

### Good Practices

1. **Immutability**: Prefer immutable data structures
2. **Explicit Dependencies**: Pass context as parameters
3. **Thread-Local Storage**: Per-thread data avoids sharing
4. **Composition Over Mutation**: Wrap with context instead of mutating
5. **Type Safety**: Compiler catches missing parameters

## Conclusion

The root cause is **mutable shared state** (the `nce_model::Base.Ref`) in a global structure (`SplineFragmentLookup`) being accessed concurrently without synchronization, triggering GC corruption.

The proper fix is **Solution 1**: Remove NCE model from library, pass as explicit parameter. This:
- Eliminates mutable shared state
- Makes dependencies explicit
- Is thread-safe by design
- Follows functional programming principles
- Improves testability

Short-term workarounds (GC protection + caching) can provide immediate relief while we implement the proper architectural fix.
