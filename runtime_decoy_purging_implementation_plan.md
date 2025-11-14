# Runtime Decoy Purging Implementation Plan

**Date:** 2025-11-14
**Branch:** `runtime-decoy-purging` (from `develop`)
**Author:** Claude Code

---

## 1. Overview

### Objective
Implement runtime decoy purging after FirstPassSearch, allowing users to reduce the fraction of decoys used in downstream searches without rebuilding the spectral library. This provides maximum flexibility for testing different decoy fractions while keeping library files consistent.

### Key Innovation
Instead of building libraries with fewer decoys, we purge a user-defined percentage of decoys from the `precursor_dict` after FirstPassSearch, then update the FDR scale factor to maintain accurate q-value calculations throughout the remaining search pipeline.

### Design Rationale
- **Flexibility:** One library file, test multiple decoy fractions
- **Experimentation-friendly:** Change fractions without rebuilding libraries
- **Stage-specific:** Could apply different fractions at different search stages
- **User-friendly:** Runtime parameter, not build-time commitment

---

## 2. Technical Background

### Current Architecture

#### FirstPassSearch Flow
```
FirstPassSearch
  ├─> performSearch!()
  │     └─> Process all MS files in parallel
  │           └─> Generate PSMs with q-values
  │
  └─> summarizeResults!()
        ├─> map_retention_times!()
        ├─> get_best_precursors_accross_runs!()  ← **PURGE HERE**
        │     ├─> Reads all FirstPass PSMs
        │     ├─> Filters by q-value (default: 0.01)
        │     └─> Returns precursor_dict with best PSM per precursor
        │
        ├─> setPrecursorDict!()  ← Store in SearchContext
        └─> create_rt_indices!()
```

#### Precursor Dict Structure
```julia
Dictionary{UInt32, @NamedTuple{
    best_prob::Float32,         # Best probability score
    best_ms_file_idx::UInt32,   # File with best PSM
    best_scan_idx::UInt32,      # Scan with best PSM
    best_irt::Float32,          # iRT at best identification
    mean_irt::Union{Missing, Float32},  # Mean iRT across identifications
    var_irt::Union{Missing, Float32},   # Variance of iRT
    n::Union{Missing, UInt16},  # Number of identifications
    mz::Float32                 # Precursor m/z
}}
```

**Key:** UInt32 precursor ID (index into library precursors)

#### FDR Scale Factor
Located in `SearchContext` (mutable struct):
```julia
n_library_targets::Int64
n_library_decoys::Int64
library_fdr_scale_factor::Float32  # = n_targets / max(n_decoys, 1)
```

Initialized in `SearchMethods.jl:358-366` by counting library precursors.

---

## 3. Implementation Plan

### 3.1 Add Configuration Parameter

#### Files to Modify
- `assets/example_config/defaultSearchParams.json`
- Any other search parameter templates

#### Parameter Addition
Add to the search parameters section (not library building):
```json
{
  "search": {
    "runtime_decoy_fraction": 1.0
  }
}
```

**Parameter Specification:**
- `runtime_decoy_fraction` (Float): Fraction of identified decoys to retain after FirstPassSearch
  - Range: (0.0, 1.0]
  - Default: 1.0 (no purging, maintains current behavior)
  - Example: 0.1 means keep only 10% of decoys randomly

**Important:** Value must be > 0.0. Setting to 0.0 would eliminate all decoys and break FDR calculations.

---

### 3.2 Parameter Validation

#### File
`src/Routines/SearchDIA/searchRAW.jl` (or parameter validation module)

#### Location
In parameter checking/loading function

#### Implementation
```julia
# Validate runtime_decoy_fraction
if haskey(params, "runtime_decoy_fraction")
    frac = params["runtime_decoy_fraction"]
    if !(frac isa Real)
        error("runtime_decoy_fraction must be a number, got: $(typeof(frac))")
    end
    if frac <= 0.0 || frac > 1.0
        error("runtime_decoy_fraction must be in range (0.0, 1.0], got: $frac")
    end
    params["runtime_decoy_fraction"] = Float64(frac)
else
    params["runtime_decoy_fraction"] = 1.0  # Default: no purging
end
```

---

### 3.3 Core Purging Function

#### File
`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

#### Location
After `get_best_precursors_accross_runs!()` call (line 563), before `setPrecursorDict!()` (line 599)

#### New Function: `purge_decoys_from_precursor_dict!`

```julia
"""
    purge_decoys_from_precursor_dict!(
        precursor_dict::Dictionary{UInt32, NamedTuple},
        spec_lib::SpectralLibrary,
        decoy_fraction::Float64,
        random_seed::Union{Int, Nothing} = nothing
    ) -> (n_targets_retained::Int64, n_decoys_retained::Int64)

Randomly removes a fraction of decoys from the precursor dictionary.

# Arguments
- `precursor_dict`: Dictionary mapping precursor IDs to best PSM info
- `spec_lib`: Spectral library containing precursor metadata (is_decoy flags)
- `decoy_fraction`: Fraction of decoys to RETAIN (0.0 < fraction ≤ 1.0)
- `random_seed`: Optional seed for reproducible decoy selection

# Returns
- Tuple of (n_targets_retained, n_decoys_retained) after purging

# Algorithm
1. Separate precursor IDs into targets and decoys using library is_decoy flags
2. Calculate number of decoys to keep: ceil(n_decoys * decoy_fraction)
3. Randomly select decoys to keep
4. Remove non-selected decoys from dictionary in-place
5. Return final counts for FDR scale factor update

# Notes
- Modifies precursor_dict in-place
- All targets are always retained
- Selection is random within decoys
- Reproducible if random_seed provided
"""
function purge_decoys_from_precursor_dict!(
    precursor_dict::Dictionary{UInt32, <:NamedTuple},
    spec_lib::SpectralLibrary,
    decoy_fraction::Float64,
    random_seed::Union{Int, Nothing} = nothing
)
    # Early return if no purging needed
    if decoy_fraction >= 1.0
        # Count targets and decoys for return values
        is_decoy = getIsDecoy(getPrecursors(spec_lib))
        n_targets = sum(!is_decoy[pid] for pid in keys(precursor_dict))
        n_decoys = sum(is_decoy[pid] for pid in keys(precursor_dict))
        return (n_targets, n_decoys)
    end

    # Get is_decoy flags from library
    is_decoy = getIsDecoy(getPrecursors(spec_lib))

    # Separate precursor IDs into targets and decoys
    target_pids = UInt32[]
    decoy_pids = UInt32[]

    for pid in keys(precursor_dict)
        if is_decoy[pid]
            push!(decoy_pids, pid)
        else
            push!(target_pids, pid)
        end
    end

    n_targets = length(target_pids)
    n_decoys_original = length(decoy_pids)

    # Calculate how many decoys to keep
    n_decoys_to_keep = ceil(Int, n_decoys_original * decoy_fraction)

    @user_info "Purging decoys: keeping $n_decoys_to_keep out of $n_decoys_original decoys ($(round(decoy_fraction * 100, digits=1))%)"

    # If keeping all decoys, return early
    if n_decoys_to_keep >= n_decoys_original
        return (n_targets, n_decoys_original)
    end

    # Create RNG for reproducible selection
    rng = isnothing(random_seed) ? Random.GLOBAL_RNG : Random.MersenneTwister(random_seed)

    # Randomly select decoys to keep
    shuffle!(rng, decoy_pids)
    decoys_to_keep = Set(decoy_pids[1:n_decoys_to_keep])
    decoys_to_remove = decoy_pids[(n_decoys_to_keep+1):end]

    # Remove non-selected decoys from dictionary
    for pid in decoys_to_remove
        delete!(precursor_dict, pid)
    end

    n_decoys_removed = n_decoys_original - n_decoys_to_keep
    @user_info "Removed $n_decoys_removed decoys from precursor dictionary"

    return (n_targets, n_decoys_to_keep)
end
```

---

### 3.4 Update FDR Scale Factor

#### File
`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

#### Location
After `purge_decoys_from_precursor_dict!()` call, before `setPrecursorDict!()`

#### Implementation
```julia
# In summarizeResults! function, after line 563

# Get runtime decoy fraction parameter
runtime_decoy_fraction = get(getParams(search_context), "runtime_decoy_fraction", 1.0)

# Purge decoys if fraction < 1.0
if runtime_decoy_fraction < 1.0
    random_seed = get(getParams(search_context), "runtime_decoy_random_seed", nothing)

    n_targets_retained, n_decoys_retained = purge_decoys_from_precursor_dict!(
        precursor_dict,
        getSpecLib(search_context),
        runtime_decoy_fraction,
        random_seed
    )

    # Update FDR scale factor in SearchContext
    new_fdr_scale_factor = Float32(n_targets_retained / max(n_decoys_retained, 1))
    search_context.library_fdr_scale_factor = new_fdr_scale_factor
    search_context.n_library_targets = n_targets_retained
    search_context.n_library_decoys = n_decoys_retained

    @user_info "Updated FDR scale factor from $(getLibraryFdrScaleFactor(search_context)) to $new_fdr_scale_factor"
    @user_info "Effective library composition: $n_targets_retained targets, $n_decoys_retained decoys"
end

setPrecursorDict!(search_context, precursor_dict)
```

---

### 3.5 Optional: Add Random Seed Parameter

For reproducibility, add optional seed parameter:

#### Configuration
```json
{
  "search": {
    "runtime_decoy_fraction": 0.1,
    "runtime_decoy_random_seed": 42
  }
}
```

This ensures the same decoys are purged across repeated runs with identical parameters.

---

## 4. Downstream Impact Analysis

### What Changes
✅ **SecondPassSearch:** Uses reduced precursor_dict, fewer precursors to search
✅ **RT Indices:** Built from reduced precursor_dict
✅ **ScoringSearch:** Uses updated FDR scale factor automatically
✅ **IntegrateChromatogramsSearch:** Works with reduced precursor set
✅ **MaxLFQSearch:** Quantifies reduced precursor set

### What Stays the Same
✅ **Library files:** Unchanged on disk
✅ **FirstPassSearch:** Searches full library (no speedup here)
✅ **Q-value calculations:** Automatically correct with updated scale factor
✅ **All FDR formulas:** Already handle arbitrary ratios correctly

### Performance Characteristics
- **Memory:** Reduced after FirstPassSearch (smaller precursor_dict)
- **FirstPassSearch time:** Unchanged (still searches full library)
- **SecondPassSearch time:** Reduced (fewer precursors to search)
- **ScoringSearch time:** Reduced (fewer PSMs to score)
- **Integration time:** Reduced (fewer chromatograms to extract)

---

## 5. Testing Strategy

### 5.1 Unit Tests

**Test Purging Logic:**
```julia
@testset "purge_decoys_from_precursor_dict!" begin
    # Create mock precursor_dict
    # Test fraction = 1.0 (no purging)
    # Test fraction = 0.5 (half decoys)
    # Test fraction = 0.1 (90% removed)
    # Test reproducibility with seed
    # Test edge cases (no decoys, all decoys)
end
```

**Test FDR Scale Factor Update:**
```julia
@testset "FDR scale factor after purging" begin
    # Verify correct calculation
    # Test with various fractions
    # Ensure scale factor > 1.0 after purging
end
```

### 5.2 Integration Tests

**Full Pipeline Test:**
1. Run SearchDIA with `runtime_decoy_fraction = 1.0` (baseline)
2. Run SearchDIA with `runtime_decoy_fraction = 0.1` (purged)
3. Compare:
   - Number of precursors in SecondPassSearch
   - Number of high-confidence PSMs
   - Q-value distributions
   - Final protein groups

**Expected Outcomes:**
- Fewer precursors searched in SecondPass (proportional to fraction)
- Slightly different q-values (due to different null distribution)
- Similar high-confidence identifications (robust targets survive)
- Computational savings in downstream stages

### 5.3 Validation Tests

**FDR Control Validation:**
- Plot empirical FDR vs target FDR at multiple thresholds
- Verify FDR control is maintained with purged decoys
- Check that target-decoy competition still works

**Reproducibility:**
- Run twice with same seed, verify identical results
- Run twice without seed, verify different decoy selection
- Verify q-values are deterministic given same decoy set

---

## 6. Files Modified Summary

### New Files
None (all modifications to existing files)

### Modified Files (3 total)

1. **`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`**
   - Add `purge_decoys_from_precursor_dict!()` function (~80 lines)
   - Modify `summarizeResults!()` to call purging and update FDR (~20 lines)
   - **Total: ~100 lines added/modified**

2. **`src/Routines/SearchDIA/searchRAW.jl`** (or parameter validation module)
   - Add parameter validation for `runtime_decoy_fraction` (~15 lines)
   - **Total: ~15 lines**

3. **`assets/example_config/defaultSearchParams.json`**
   - Add `runtime_decoy_fraction` parameter (2 lines)
   - Add optional `runtime_decoy_random_seed` parameter (1 line)
   - **Total: ~3 lines**

### Total Implementation Size
- **~120 lines of code**
- **3 files modified**
- **No changes to FDR calculation logic** (already correct)

---

## 7. Comparison: Runtime vs Library-Time Approach

| Aspect | Runtime Purging (This Plan) | Library-Time Reduction (Other Branch) |
|--------|----------------------------|--------------------------------------|
| **Flexibility** | ✅ High - change anytime | ❌ Low - rebuild required |
| **Library size** | ❌ Full size | ✅ Reduced size |
| **FirstPass speed** | ❌ No improvement | ✅ Faster (fewer precursors) |
| **Downstream speed** | ✅ Improved | ✅ Improved |
| **Memory usage** | ⚠️ Full until FirstPass | ✅ Reduced throughout |
| **User workflow** | ✅ Simple (one library) | ⚠️ More complex (multiple libraries) |
| **Experimentation** | ✅ Excellent | ❌ Cumbersome |
| **Production use** | ⚠️ Good | ✅ Optimal |

### Recommendation
- **Development/Research:** Use runtime purging for flexibility
- **Production pipelines:** Use library-time reduction for performance
- **Ideal:** Implement both, let users choose

---

## 8. Implementation Checklist

### Pre-Implementation
- [x] Understand precursor_dict structure
- [x] Locate FDR scale factor storage
- [x] Identify purging location in pipeline
- [x] Confirm downstream FDR usage
- [x] Create detailed implementation plan

### Parameter Configuration
- [ ] Add `runtime_decoy_fraction` to search params
- [ ] Add optional `runtime_decoy_random_seed`
- [ ] Add parameter validation

### Core Implementation
- [ ] Implement `purge_decoys_from_precursor_dict!()` function
- [ ] Add purging call in `summarizeResults!()`
- [ ] Update FDR scale factor after purging
- [ ] Add logging for purging statistics

### Testing
- [ ] Unit test purging function
- [ ] Unit test FDR scale factor update
- [ ] Integration test with full pipeline
- [ ] Validate FDR control maintained
- [ ] Test reproducibility with seeds

### Documentation
- [ ] Update user documentation
- [ ] Add example configurations
- [ ] Document performance characteristics
- [ ] Note differences from library-time approach

### Finalization
- [ ] Code review
- [ ] Performance benchmarking
- [ ] Commit changes
- [ ] Create pull request

---

## 9. Usage Examples

### Example 1: Default Behavior (No Purging)
```json
{
  "search": {
    "runtime_decoy_fraction": 1.0
  }
}
```
**Result:** All decoys retained, identical to current behavior

### Example 2: Aggressive Purging
```json
{
  "search": {
    "runtime_decoy_fraction": 0.1
  }
}
```
**Result:** Only 10% of decoys retained after FirstPassSearch, significant downstream speedup

### Example 3: Reproducible Purging
```json
{
  "search": {
    "runtime_decoy_fraction": 0.2,
    "runtime_decoy_random_seed": 42
  }
}
```
**Result:** 20% of decoys retained, same decoys selected every run

---

## 10. Future Enhancements

### Potential Extensions

1. **Adaptive Purging:**
   - Purge more aggressively for high-confidence decoys
   - Keep decoys with PSMs close to threshold

2. **Stratified Sampling:**
   - Ensure decoys span full m/z range
   - Maintain representation across RT bins

3. **Score-Based Purging:**
   - Prioritize removal of low-scoring decoys
   - Keep challenging decoys for robust FDR

4. **Stage-Specific Fractions:**
   - Different fractions for SecondPass vs Scoring
   - Progressive refinement of decoy set

5. **Quality-Weighted Purging:**
   - Consider PSM quality, not just random selection
   - Intelligent decoy selection based on competition

### Not in Current Scope
These enhancements require additional research and validation. Current implementation focuses on simple random purging with correct FDR calculations.

---

## 11. Risk Assessment

### Low Risk
✅ Parameter addition (well-defined, validated)
✅ Dictionary manipulation (straightforward delete operations)
✅ FDR scale factor update (mutable struct, direct access)

### Medium Risk
⚠️ **Downstream impacts:** Must verify all methods use updated scale factor
⚠️ **FDR validation:** Requires empirical testing of FDR control

### Mitigation Strategies
1. Comprehensive testing with multiple fraction values
2. Validation that fraction=1.0 produces identical results
3. FDR calibration plots at multiple thresholds
4. Code review focusing on FDR calculation usage

### Critical Invariants to Maintain
- `n_targets + n_decoys = length(precursor_dict)` after purging
- `library_fdr_scale_factor = n_targets / n_decoys` exactly
- All targets always retained (never purge targets)
- Decoy selection is random (no bias)

---

## 12. Performance Estimates

### Computational Savings (with 0.1 fraction)

**FirstPassSearch:** 0% (searches full library)
**SecondPassSearch:** ~40-50% (fewer precursors, but still processes all scans)
**ScoringSearch:** ~20-30% (fewer PSMs to score)
**IntegrateChromatograms:** ~50-60% (fewer chromatograms to extract)
**MaxLFQ:** ~50% (smaller quantification matrix)

### Memory Savings

**precursor_dict:** ~90% reduction (proportional to fraction)
**RT indices:** ~90% reduction
**PSM storage:** ~50% reduction (fewer decoys pass filters)
**Quantification matrices:** ~50% reduction

### Disk Savings

❌ **None** - Library files remain full size

---

## 13. Implementation Timeline

### Phase 1: Core Implementation (4-6 hours)
- Implement purging function
- Add parameter handling
- Update FDR scale factor
- Basic testing

### Phase 2: Testing & Validation (4-6 hours)
- Unit tests
- Integration tests
- FDR validation
- Performance benchmarking

### Phase 3: Documentation (2-3 hours)
- Code documentation
- User documentation
- Example configurations

### Total Estimated Time: 10-15 hours

---

## 14. Key Takeaways

### What This Approach Offers
1. **Maximum flexibility** - test different fractions without rebuilding
2. **Minimal code changes** - ~120 lines, 3 files
3. **Correct FDR calculations** - automatic with scale factor update
4. **Reproducible** - optional seeding for deterministic results

### What It Doesn't Offer
1. FirstPassSearch speedup (still searches full library)
2. Disk space savings (library file unchanged)
3. Early-stage memory savings (full library loaded initially)

### When to Use This vs Library-Time Reduction
- **Use runtime purging:** Research, experimentation, parameter tuning
- **Use library-time:** Production pipelines, resource-constrained environments
- **Use both:** Best of both worlds, user can choose

---

**End of Implementation Plan**
