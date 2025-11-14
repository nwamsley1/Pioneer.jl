# Fractional Decoy Generation Implementation Plan

**Date:** 2025-11-14
**Branch:** `fractional-decoys` (to be created from `develop`)
**Author:** Claude Code

---

## 1. Overview

### Objective
Modify Pioneer.jl to support fractional decoy generation, allowing users to specify what fraction of target sequences should have decoys generated. This reduces library size while maintaining statistically correct FDR calculations.

### Motivation
- Current system generates 1 decoy per target sequence (1:1 ratio)
- Large proteomes can result in very large libraries
- Many applications don't require 1:1 ratio for accurate FDR control
- Fractional decoys (e.g., 0.1 = 10% of targets) can significantly reduce library size

### Key Design Decisions
Based on user feedback:
1. **Application Level:** Base sequence level (before charge variants)
2. **Reproducibility:** Use existing or new random seed parameter
3. **Default Value:** 1.0 (100% decoys) for backward compatibility
4. **Minimum Decoys:** No minimum enforced (strict fraction application)

---

## 2. Technical Background

### Current Implementation

#### Decoy Generation Flow
```
prepare_chronologer_input() [chronologer_prep.jl:202-205]
  └─> add_decoy_sequences_grouped() [fasta_utils.jl:861-994]
        ├─> Groups entries by base sequence
        ├─> For EACH base sequence group:
        │     ├─> Generate one decoy sequence (shuffle or reverse)
        │     └─> Create decoy variants for all modifications/charges
        └─> Returns target + decoy entries combined
```

#### Current FDR Calculation
Located in `SearchMethods.jl:350-365`:
```julia
# Count targets and decoys in library
for i in 1:length(is_decoy_array)
    if is_decoy_array[i]
        n_decoys += 1
    else
        n_targets += 1
    end
end

# Calculate FDR scale factor
library_fdr_scale_factor = n_targets / max(n_decoys, 1)
```

Q-value calculation in `fdrUtilities.jl:48-75`:
```julia
qvals[i] = (decoys * fdr_scale_factor) / targets
```

**Critical Insight:** The existing FDR formula already handles arbitrary target/decoy ratios correctly. If we have 100 targets and 10 decoys (0.1 fraction), the scale factor is 10, meaning each decoy represents 10 targets' worth of false discoveries. **No changes to FDR/q-value code are needed.**

---

## 3. Implementation Plan

### 3.1 Add New Parameters

#### Files to Modify
1. `assets/example_config/defaultBuildLibParams.json`
2. `assets/example_config/defaultBuildLibParamsSimplified.json`

#### Changes
Add to `fasta_digest_params` section (after `decoy_method`):
```json
"decoy_fraction": 1.0,
"decoy_random_seed": null
```

**Parameter Specifications:**
- `decoy_fraction` (Float): Fraction of target base sequences to generate decoys for
  - Range: [0.0, 1.0]
  - Default: 1.0 (maintains current 1:1 behavior)
  - Example: 0.1 means generate decoys for 10% of base sequences

- `decoy_random_seed` (Integer or null): Random seed for reproducible decoy selection
  - Range: Any positive integer, or null for random selection
  - Default: null
  - Purpose: Ensures same decoys are selected across repeated builds

---

### 3.2 Parameter Validation

#### File
`src/Routines/BuildSpecLib/utils/check_params.jl`

#### Location
After line 95 (after `entrapment_method` validation, before `nce_params`)

#### Implementation
```julia
# Check decoy_fraction with default value
if !haskey(fasta_digest_params, "decoy_fraction")
    fasta_digest_params["decoy_fraction"] = 1.0
else
    decoy_fraction = fasta_digest_params["decoy_fraction"]
    if !(decoy_fraction isa Real)
        error("decoy_fraction must be a number, got: $(typeof(decoy_fraction))")
    end
    if decoy_fraction < 0.0 || decoy_fraction > 1.0
        error("decoy_fraction must be between 0.0 and 1.0, got: $decoy_fraction")
    end
    # Store as Float64 for consistency
    fasta_digest_params["decoy_fraction"] = Float64(decoy_fraction)
end

# Check decoy_random_seed (optional parameter)
if haskey(fasta_digest_params, "decoy_random_seed")
    seed = fasta_digest_params["decoy_random_seed"]
    if !isnothing(seed)
        if !(seed isa Integer)
            error("decoy_random_seed must be an integer or null, got: $(typeof(seed))")
        end
        if seed <= 0
            error("decoy_random_seed must be positive, got: $seed")
        end
    end
end
```

---

### 3.3 Core Decoy Generation Function

#### File
`src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

#### Function
`add_decoy_sequences_grouped()` (currently lines 861-994)

#### Signature Change
```julia
function add_decoy_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry};
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    decoy_method::String = "shuffle",
    decoy_fraction::Float64 = 1.0,
    random_seed::Union{Int,Nothing} = nothing
)::Vector{FastaEntry}
```

#### Implementation Details

**Step 1: Add imports (if not already present)**
At the top of the file, ensure:
```julia
using Random
```

**Step 2: Create RNG if seed provided**
After the `shuffle_seq` initialization (around line 880), add:
```julia
# Create random number generator for reproducible decoy selection
rng = isnothing(random_seed) ? Random.GLOBAL_RNG : Random.MersenneTwister(random_seed)
```

**Step 3: Select base sequences for decoy generation**
After grouping by base sequence (around line 890), add:
```julia
# Determine which base sequences will get decoys
total_base_sequences = length(groups)
n_decoys_to_generate = ceil(Int, total_base_sequences * decoy_fraction)

# Get list of base sequences and randomly select subset
all_base_seqs = collect(keys(groups))
if decoy_fraction < 1.0
    # Shuffle and take first n_decoys_to_generate
    shuffle!(rng, all_base_seqs)
    selected_base_seqs = Set(all_base_seqs[1:n_decoys_to_generate])
else
    # Generate decoys for all sequences (current behavior)
    selected_base_seqs = Set(all_base_seqs)
end

@user_info "Generating decoys for $n_decoys_to_generate out of $total_base_sequences base sequences (fraction: $decoy_fraction)"
```

**Step 4: Modify main loop**
In the main loop `for (base_seq, idxs) in groups` (around line 902), add condition:
```julia
for (base_seq, idxs) in groups
    # Skip if this base sequence not selected for decoy generation
    if base_seq ∉ selected_base_seqs
        continue
    end

    # ... rest of existing code unchanged ...
```

**Step 5: Update final statistics**
Before the return statement (around line 992), update logging:
```julia
@user_info "Decoy generation complete: $(length(decoy_entries)) decoy entries from $n_decoys_to_generate base sequences"
```

---

### 3.4 Update Function Call

#### File
`src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

#### Location
Lines 202-205

#### Changes
Replace:
```julia
# Step 6: Add decoys (GROUPED by base sequence; all mods share same decoy)
if _params.fasta_digest_params["add_decoys"]
    decoy_method = get(_params.fasta_digest_params, "decoy_method", "shuffle")
    fasta_entries = add_decoy_sequences_grouped(fasta_entries; decoy_method=decoy_method)
end
```

With:
```julia
# Step 6: Add decoys (GROUPED by base sequence; all mods share same decoy)
if _params.fasta_digest_params["add_decoys"]
    decoy_method = get(_params.fasta_digest_params, "decoy_method", "shuffle")
    decoy_fraction = get(_params.fasta_digest_params, "decoy_fraction", 1.0)
    decoy_random_seed = get(_params.fasta_digest_params, "decoy_random_seed", nothing)
    fasta_entries = add_decoy_sequences_grouped(
        fasta_entries;
        decoy_method=decoy_method,
        decoy_fraction=decoy_fraction,
        random_seed=decoy_random_seed
    )
end
```

---

## 4. Files NOT Modified

### FDR/Q-value Calculation (No Changes Required)

The following files require **zero modifications** because the existing FDR formula already handles arbitrary target/decoy ratios:

1. **`src/utils/ML/fdrUtilities.jl`**
   - `get_qvalues!()` function (lines 48-75)
   - `get_PEP!()` function (lines 96-131)
   - Formula: `qvals[i] = (decoys * fdr_scale_factor) / targets`
   - Already scales decoys appropriately

2. **`src/Routines/SearchDIA/SearchMethods/SearchMethods.jl`**
   - Scale factor calculation (lines 350-365)
   - Formula: `library_fdr_scale_factor = n_targets / max(n_decoys, 1)`
   - Automatically adjusts for any target/decoy ratio

3. **All search method implementations:**
   - FirstPassSearch
   - ParameterTuningSearch
   - ScoringSearch
   - All correctly use the scale factor from SearchContext

### Why No Changes Needed

**Mathematical Proof:**

With 1:1 ratio (current):
- 100 targets, 100 decoys
- Scale factor = 100/100 = 1
- Q-value = (decoys × 1) / targets

With 0.1 fraction (new):
- 100 targets, 10 decoys
- Scale factor = 100/10 = 10
- Q-value = (decoys × 10) / targets

Both formulas correctly estimate FDR. The scale factor compensates for the imbalanced library.

---

## 5. Testing Strategy

### Unit Tests
Create test cases for:

1. **Backward Compatibility**
   - Test with `decoy_fraction = 1.0`
   - Verify identical output to current implementation
   - Check that all base sequences get decoys

2. **Fractional Generation**
   - Test with `decoy_fraction = 0.1`
   - Verify approximately 10% of base sequences get decoys
   - Confirm exact count: `ceil(n_base_seqs * 0.1)`

3. **Reproducibility**
   - Run with same `decoy_random_seed` multiple times
   - Verify identical decoy selection
   - Test with different seeds produce different selections

4. **Edge Cases**
   - `decoy_fraction = 0.0`: No decoys generated
   - `decoy_fraction = 0.01`: Very few decoys (small library)
   - `decoy_fraction = 0.99`: Almost all sequences
   - `decoy_fraction = 1.0`: All sequences (current behavior)

5. **Parameter Validation**
   - Invalid fractions (< 0.0, > 1.0): Should error
   - Invalid seed (negative, non-integer): Should error
   - Missing parameters: Should use defaults

### Integration Tests

1. **Full Library Build**
   - Build library with 0.1 fraction
   - Verify library size reduction (~90% fewer decoys at precursor level)
   - Confirm library loads correctly

2. **FDR Calculation**
   - Run search on test data
   - Compare q-values from 1.0 vs 0.1 fraction libraries
   - Verify FDR control is maintained
   - Check that target PSMs have reasonable q-values

3. **Statistical Validity**
   - Build multiple libraries with different fractions (0.1, 0.3, 0.5, 1.0)
   - Run same search on all
   - Plot ROC curves to verify discrimination power
   - Ensure FDR control at various thresholds

---

## 6. Implementation Checklist

### Pre-Implementation
- [x] Research codebase structure
- [x] Identify all decoy generation sites
- [x] Identify all FDR calculation sites
- [x] Confirm FDR math handles fractional decoys
- [x] Get user input on design decisions
- [x] Create detailed implementation plan

### Parameter Configuration
- [ ] Update `defaultBuildLibParams.json`
- [ ] Update `defaultBuildLibParamsSimplified.json`
- [ ] Verify JSON syntax is valid

### Parameter Validation
- [ ] Add `decoy_fraction` validation to `check_params.jl`
- [ ] Add `decoy_random_seed` validation
- [ ] Test parameter validation with valid/invalid inputs

### Core Implementation
- [ ] Import `Random` module in `fasta_utils.jl`
- [ ] Update `add_decoy_sequences_grouped()` signature
- [ ] Add RNG initialization
- [ ] Add base sequence selection logic
- [ ] Modify main loop with selection check
- [ ] Update logging messages
- [ ] Update function call in `chronologer_prep.jl`

### Testing
- [ ] Test backward compatibility (fraction = 1.0)
- [ ] Test fractional generation (fraction = 0.1)
- [ ] Test reproducibility with seed
- [ ] Test edge cases (0.0, 0.01, 0.99, 1.0)
- [ ] Run full library build
- [ ] Verify FDR calculations
- [ ] Build and commit changes

---

## 7. Risk Assessment

### Low Risk
- **Parameter addition:** Simple, well-defined parameters
- **FDR calculations:** No changes needed (proven mathematically correct)
- **Backward compatibility:** Default value maintains current behavior

### Medium Risk
- **Random selection logic:** Need careful testing for edge cases
- **Reproducibility:** Must verify seed handling works correctly

### Mitigation Strategies
1. Extensive testing with various fraction values
2. Verification that fraction=1.0 produces identical results to current
3. Statistical validation on real data
4. Clear documentation of new parameters

---

## 8. Documentation Requirements

### Code Comments
- Document new parameters in function signatures
- Explain fraction-based selection logic
- Note that FDR calculations remain unchanged

### User Documentation
- Update parameter reference documentation
- Add example use cases for fractional decoys
- Explain impact on library size and FDR estimation
- Provide guidance on choosing appropriate fractions

### Example Configurations
Create example configs demonstrating:
- Standard use (fraction = 1.0)
- Reduced library (fraction = 0.1)
- Reproducible builds (with seed)

---

## 9. Performance Considerations

### Expected Performance Impact
- **Memory:** Reduced (fewer decoy entries stored)
- **Generation time:** Slightly faster (fewer decoys to generate)
- **Search time:** Faster (smaller library to search)
- **FDR calculation:** No change (same complexity)

### Optimizations
- Use `Set` for O(1) membership testing of selected base sequences
- Shuffle in-place with `shuffle!()` to avoid allocation
- Pre-calculate `n_decoys_to_generate` once

---

## 10. Future Enhancements

### Potential Extensions
1. **Minimum decoy count:** Optional parameter to ensure statistical validity
2. **Stratified sampling:** Select decoys proportionally across protein groups
3. **Targeted selection:** Prioritize decoys for high-abundance proteins
4. **Dynamic fractions:** Adjust fraction based on library size

### Not in Scope (Current Implementation)
- Precursor-level fractional decoys (too complex, inconsistent)
- Automatic fraction selection (requires domain expertise)
- Fraction validation based on statistical power (future research)

---

## 11. Summary

### Files Modified (5 total)
1. `assets/example_config/defaultBuildLibParams.json` - Add parameters
2. `assets/example_config/defaultBuildLibParamsSimplified.json` - Add parameters
3. `src/Routines/BuildSpecLib/utils/check_params.jl` - Validate parameters (~20 lines)
4. `src/Routines/BuildSpecLib/fasta/fasta_utils.jl` - Selection logic (~30 lines)
5. `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl` - Pass parameters (~5 lines)

### Lines of Code
- **Added:** ~55 lines
- **Modified:** ~10 lines
- **Total impact:** ~65 lines

### Key Benefits
- **Minimal changes:** Small, focused modification
- **No FDR changes:** Existing statistical framework already correct
- **Backward compatible:** Default behavior unchanged
- **Reproducible:** Optional seeding for consistent results
- **Flexible:** User controls exact fraction

### Implementation Time Estimate
- **Coding:** 1-2 hours
- **Testing:** 2-3 hours
- **Documentation:** 1 hour
- **Total:** 4-6 hours

---

## 12. References

### Key Code Locations

**Decoy Generation:**
- `src/Routines/BuildSpecLib/fasta/fasta_utils.jl:861-994`
- `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl:202-205`

**FDR Calculation:**
- `src/utils/ML/fdrUtilities.jl:48-75` (get_qvalues!)
- `src/utils/ML/fdrUtilities.jl:96-131` (get_PEP!)
- `src/Routines/SearchDIA/SearchMethods/SearchMethods.jl:350-365`

**Parameter Management:**
- `src/Routines/BuildSpecLib/utils/check_params.jl:59-96`
- `assets/example_config/defaultBuildLibParams.json`
- `assets/example_config/defaultBuildLibParamsSimplified.json`

---

**End of Implementation Plan**
