# Plan: Remove `use_for_quant` Dictionary from `InferenceResult`

**Date:** 2025-01-11
**Author:** Analysis by Claude Code
**Status:** Proposed - Awaiting Approval

---

## Executive Summary

The `use_for_quant` dictionary in `InferenceResult` is **completely redundant** after our recent changes that delete shared peptides from results. We can simplify the API by removing it and using presence in `peptide_to_protein` as the sole indicator of whether a peptide should be used for quantification.

**Key Insight:** After deleting shared peptides, if a peptide is in `peptide_to_protein`, it's always usable for quantification. If it's not in the dictionary, it shouldn't be used. This makes `use_for_quant` redundant.

---

## Current State Analysis

### The InferenceResult Struct
```julia
struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
    use_for_quant::Dictionary{PeptideKey, Bool}
end
```

### After Recent Changes (Commit: `refactor(protein_inference): simplify output by removing shared peptides`)

Shared peptides are deleted from BOTH dictionaries:
```julia
# In infer_proteins() - lines 478-484
for peptide_key in component_peptides
    if haskey(use_for_quant, peptide_key) && !use_for_quant[peptide_key]
        delete!(peptide_to_protein, peptide_key)
        delete!(use_for_quant, peptide_key)
    end
end
```

**Invariant:** For all `peptide_key`:
- If `peptide_key ∈ peptide_to_protein`, then `use_for_quant[peptide_key] == true`
- If `peptide_key ∉ peptide_to_protein`, then `peptide_key ∉ use_for_quant`

**Conclusion:** The `use_for_quant` dictionary is just a duplicate of the keys in `peptide_to_protein`.

---

## The Problem: Incorrect Fallback in `add_quantification_flag`

When a peptide exists in PSM data but NOT in inference results, the current code has **mismatched logic**:

### `add_inferred_protein_column` (lines 179-185) - ✅ CORRECT
```julia
if haskey(inference_result.peptide_to_protein, pep_key)
    inferred_proteins[i] = inference_result.peptide_to_protein[pep_key].name
else
    # Fallback to original protein - CORRECT behavior
    # Peptide wasn't in inference, so keep its original assignment
    inferred_proteins[i] = accession_numbers[i]
end
```

**This is appropriate:** Every PSM needs a protein assignment. If inference didn't include this peptide, fall back to the original protein from the library.

### `add_quantification_flag` (lines 219-224) - ❌ INCORRECT
```julia
if haskey(inference_result.use_for_quant, pep_key)
    use_for_quant[i] = inference_result.use_for_quant[pep_key]
else
    # Default to true if not in inference results
    use_for_quant[i] = true  # ← WRONG! Should be FALSE
end
```

**This is wrong:** If a peptide was deleted from inference results, it means it's a **shared peptide** that shouldn't be used for quantification. Defaulting to TRUE violates the principle of using only unique peptides.

### Why This is Problematic

If a shared peptide was deleted from inference results:
1. `inferred_protein_group` = original protein (e.g., "ProteinA" or "ProteinB") ✅ OK
2. `use_for_protein_quant` = TRUE (default) ❌ WRONG

**This means shared peptides are incorrectly marked as usable for quantification!**

The correct behavior should be:
- **In inference result** → Use inferred protein AND mark as usable (unique peptide)
- **NOT in inference result** → Use original protein BUT mark as NOT usable (shared peptide, deleted)

### When Can Peptides Be Missing from Inference Results?

**Case: Same peptide from different proteins (the normal shared peptide case)**
- PSM row 1: peptide "ABCD", protein "ProteinA"
- PSM row 2: peptide "ABCD", protein "ProteinB"
- PSM row 3: peptide "ABCD", protein "ProteinA" (duplicate row from same precursor)

The `unique(df, [:sequence, :accession_numbers, :is_decoy, :entrap_id])` call (line 126) creates TWO unique pairs:
- ("ABCD", "ProteinA")
- ("ABCD", "ProteinB")

After inference, if "ABCD" is shared between A and B → deleted from results.

When `add_quantification_flag` runs on ALL THREE PSM rows, it encounters a missing key and **incorrectly defaults to TRUE**.

---

## Proposed Solution: Remove `use_for_quant` Dictionary

### Rationale

1. **Eliminates redundancy** - Single source of truth
2. **Simpler API** - Smaller `InferenceResult` struct
3. **Less memory** - Half the dictionary storage
4. **Clearer semantics** - Presence in dict = usable for quantification
5. **Prevents inconsistency** - Can't have conflicting values
6. **Fixes the bug** - Forces consistent default behavior

### New Semantics

**Presence-based logic for quantification:**
- `haskey(peptide_to_protein, pep_key)` → **Use for quantification** (unique peptide assigned to minimal protein set)
- `!haskey(peptide_to_protein, pep_key)` → **Don't use for quantification** (shared peptide, deleted from inference)

**Protein assignment remains unchanged:**
- Peptides in inference → use inferred protein group
- Peptides not in inference → fall back to original protein (already correct)

### Changes Required

#### 1. Update `InferenceResult` Struct
**File:** `src/structs/protein_inference_types.jl:108-111`

```julia
# BEFORE
struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
    use_for_quant::Dictionary{PeptideKey, Bool}
end

# AFTER
struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
end
```

Update docstring to explain presence-based semantics.

#### 2. Update `infer_proteins()` Function
**File:** `src/utils/proteinInference.jl`

Remove all operations on `use_for_quant` dictionary:
- Remove dictionary initialization
- Remove insertions during peptide assignment
- Keep deletion of shared peptides from `peptide_to_protein` only
- Update return statement: `return InferenceResult(peptide_to_protein)`

#### 3. Update `add_quantification_flag()` Function
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/protein_inference_pipeline.jl:200-232`

```julia
# BEFORE
if haskey(inference_result.use_for_quant, pep_key)
    use_for_quant[i] = inference_result.use_for_quant[pep_key]
else
    # Default to true if not in inference results
    use_for_quant[i] = true
end

# AFTER
# Peptide is usable for quantification if it's in the inference results
# (i.e., it's a unique peptide assigned to a protein)
use_for_quant[i] = haskey(inference_result.peptide_to_protein, pep_key)
```

**Update function docstring** to explain the new logic.

#### 4. Update All Unit Tests
**File:** `test/UnitTests/test_protein_inference.jl`

Remove all assertions on `result.use_for_quant`:
- Remove `@test length(result.use_for_quant) == N` checks
- Remove `@test result.use_for_quant[peptides[i]] == true/false` checks
- Remove `@test !haskey(result.use_for_quant, peptides[i])` checks

Add comment explaining that quantification is determined by presence in `peptide_to_protein`.

#### 5. Update Empty Result Handling
**File:** `src/Routines/SearchDIA/SearchMethods/ScoringSearch/protein_inference_pipeline.jl:119-122`

```julia
# BEFORE
return InferenceResult(
    Dictionary{PeptideKey, ProteinKey}(),
    Dictionary{PeptideKey, Bool}()
)

# AFTER
return InferenceResult(
    Dictionary{PeptideKey, ProteinKey}()
)
```

#### 6. Update Documentation
**Files:**
- `src/Routines/SearchDIA/CLAUDE.md`
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/CLAUDE.md`
- `src/structs/protein_inference_types.jl` (docstrings)

Add clear explanation of presence-based semantics:
> **InferenceResult Semantics:**
> - Peptides present in `peptide_to_protein` are unique peptides assigned to the minimal protein set and should be used for quantification
> - Shared peptides are excluded from the result entirely (deleted after inference)
> - For PSMs with peptides not in the inference result: use original protein assignment but mark as `use_for_protein_quant = false`

---

## Impact Analysis

### Downstream Code That WILL Change

1. **`infer_proteins()` function** - Returns simpler struct
2. **`add_quantification_flag()` function** - Checks `peptide_to_protein` instead
3. **All unit tests** - Remove `use_for_quant` assertions

### Downstream Code That WON'T Change

✅ **`use_for_protein_quant` column in DataFrames** - Still populated correctly
✅ **`group_psms_by_protein()` function** - Still uses the column
✅ **MaxLFQ quantification** - Still filters by the column
✅ **Protein scoring** - Still uses the column
✅ **File formats and outputs** - No changes

**The change is purely internal to the inference algorithm.**

---

## Testing Strategy

1. **Run existing unit tests** - Verify Cases A-K still pass with updated assertions
2. **Run integration test** - `SearchDIA("./data/ecoli_test/ecoli_test_params.json")`
3. **Verify PSM outputs** - Check that `use_for_protein_quant` column is correctly populated
4. **Verify protein groups** - Ensure only unique peptides are used for scoring
5. **Verify MaxLFQ** - Ensure quantification uses correct peptides

---

## Risks and Mitigations

### Risk 1: Breaking Changes to Public API
**Mitigation:** `InferenceResult` is internal to protein inference module. No external code accesses it directly.

### Risk 2: Logic Errors in `add_quantification_flag`
**Mitigation:** The new logic is simpler (single check vs. fallback). Comprehensive tests will catch issues.

### Risk 3: Performance Impact
**Mitigation:** One dictionary instead of two → faster lookups, less memory. Performance improves.

### Risk 4: Documentation Debt
**Mitigation:** Update all relevant CLAUDE.md files and docstrings as part of implementation.

---

## Benefits Summary

1. **✅ Fixes the inconsistent default bug** - No more defaulting to TRUE incorrectly
2. **✅ Simpler code** - 50% fewer dictionary operations
3. **✅ Less memory** - One dictionary instead of two
4. **✅ Clearer semantics** - Presence = usable
5. **✅ Prevents future bugs** - Can't have conflicting values
6. **✅ Easier to understand** - Single source of truth

---

## Implementation Checklist

- [ ] Update `InferenceResult` struct definition and docstring
- [ ] Update `infer_proteins()` to remove `use_for_quant` operations
- [ ] Update `add_quantification_flag()` to check `peptide_to_protein`
- [ ] Update empty result handling in `apply_inference_to_dataframe()`
- [ ] Update all test cases to remove `use_for_quant` assertions
- [ ] Update CLAUDE.md documentation files
- [ ] Run full unit test suite (target: 129 passing tests)
- [ ] Run integration test (E. coli test dataset)
- [ ] Verify PSM output files have correct `use_for_protein_quant` values
- [ ] Commit with descriptive message

---

## Recommendation

**PROCEED** with this refactoring. The benefits clearly outweigh the risks, and it fixes a real bug in the fallback behavior while simplifying the codebase.
