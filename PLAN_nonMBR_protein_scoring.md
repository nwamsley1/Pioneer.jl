# Implementation Plan: Use nonMBR_prob for Protein Inference and Scoring

## Overview
Modify the protein inference pipeline to use `nonMBR_prob` instead of `prec_prob` for:
1. Initial protein group score calculation (pg_score)
2. CV fold assignment based on best-scoring peptide

This ensures protein-level scores are based on pre-MBR peptide probabilities, consistent with the Step 10 q-value recalculation.

## Background
Currently, protein group scores (`pg_score`) are calculated using `prec_prob`, which includes MBR enhancements. Since we now recalculate q-values in Step 10 using `nonMBR_prob`, we should also use `nonMBR_prob` for protein scoring to maintain consistency throughout the FDR calculation pipeline.

---

## Task 1: Modify Initial Protein Score Calculation

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/protein_inference_pipeline.jl`

**Change Location**: `group_psms_by_protein()` function (lines 263-276)

**Current Code** (lines 263-276):
```julia
# Calculate initial protein score (log-sum)
peptide_probs = gdf[gdf.use_for_protein_quant .== true, :prec_prob]
if isempty(peptide_probs)
    pg_score = 0.0f0
else
    # Use best probability per peptide
    unique_pep_probs = Float32[]
    for pep in quant_peptides
        pep_mask = (gdf.sequence .== pep) .& (gdf.use_for_protein_quant .== true)
        if any(pep_mask)
            push!(unique_pep_probs, maximum(gdf[pep_mask, :prec_prob]))
        end
    end
    pg_score = -sum(log.(1.0f0 .- unique_pep_probs))
end
```

**Proposed Change**:
```julia
# Calculate initial protein score (log-sum) using nonMBR probabilities
peptide_probs = gdf[gdf.use_for_protein_quant .== true, :nonMBR_prob]
if isempty(peptide_probs)
    pg_score = 0.0f0
else
    # Use best probability per peptide
    unique_pep_probs = Float32[]
    for pep in quant_peptides
        pep_mask = (gdf.sequence .== pep) .& (gdf.use_for_protein_quant .== true)
        if any(pep_mask)
            push!(unique_pep_probs, maximum(gdf[pep_mask, :nonMBR_prob]))
        end
    end
    pg_score = -sum(log.(1.0f0 .- unique_pep_probs))
end
```

**Changes**:
- Line 264: `:prec_prob` → `:nonMBR_prob`
- Line 273: `:prec_prob` → `:nonMBR_prob`

**Rationale**:
- `pg_score` is the initial protein score used for FDR calculation
- Should be based on the same probabilities used for final q-value calculation (nonMBR scores)
- The log-sum formula remains unchanged, only the input probabilities change

---

## Task 2: Modify CV Fold Assignment

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`

**Change Location**: `build_protein_cv_fold_mapping()` function (lines 1884, 1903-1904)

**Current Code** (lines 1884, 1903-1904):
```julia
# Line 1884 - required columns
required_columns = [:inferred_protein_group, :prec_prob, :precursor_idx]

# Lines 1903-1904 - finding best scoring PSM
best_idx = argmax(group.prec_prob)
best_score = group.prec_prob[best_idx]
```

**Proposed Changes**:
```julia
# Line 1884 - Update required columns
required_columns = [:inferred_protein_group, :nonMBR_prob, :precursor_idx]

# Lines 1903-1904 - Use nonMBR_prob for finding best scoring PSM
best_idx = argmax(group.nonMBR_prob)
best_score = group.nonMBR_prob[best_idx]
```

**Changes**:
- Line 1884: Add `:nonMBR_prob` to required columns (can keep `:prec_prob` for compatibility or replace)
- Line 1903: `group.prec_prob` → `group.nonMBR_prob`
- Line 1904: `group.prec_prob[best_idx]` → `group.nonMBR_prob[best_idx]`

**Rationale**:
- CV fold assignment is based on the highest-scoring peptide per protein
- Should use nonMBR scores to maintain consistency with protein score calculation
- The `best_score` field is used to track which peptide determined the CV fold

---

## Dependencies & Integration

**Upstream Dependencies**:
- ✅ `nonMBR_prob` column already added to PSM files in previous commit (Task 2 from previous plan)
- ✅ Column is written by `sort_of_percolator_in_memory!` at the end of Step 1

**Downstream Impact**:
- **Step 12** (Protein Inference): Will now use `nonMBR_prob` for initial `pg_score` calculation
- **Step 13** (CV Fold Mapping): Will now use `nonMBR_prob` to determine best peptide per protein
- **Step 14** (Probit Regression): Uses `pg_score` as input, so will be trained on nonMBR-based scores
- **Steps 15-22** (Protein Q-values): All downstream protein scoring uses `pg_score` and `global_pg_score`, which will now be based on nonMBR probabilities
- **Step 23** (PSM Updates): Backpropagates protein scores to PSMs

**Consistency Check**:
- ✅ Step 10: Uses `nonMBR_prob` for PSM q-value calculation
- ✅ Step 12: Will use `nonMBR_prob` for protein score calculation (after this change)
- ✅ Step 13: Will use `nonMBR_prob` for CV fold assignment (after this change)
- Result: Complete consistency - all FDR calculations use pre-MBR scores

---

## Testing Checklist

1. **Column Availability**: Verify `nonMBR_prob` exists in PSM files after Step 1
2. **Protein Inference**: Confirm Step 12 completes without missing column errors
3. **CV Fold Mapping**: Verify Step 13 successfully assigns CV folds
4. **Score Validation**: Check that `pg_score` values differ from previous runs (should be based on different probabilities)
5. **Downstream Steps**: Ensure Steps 14-23 complete successfully
6. **Final Output**: Verify protein groups have valid scores and q-values

---

## Alternative Approach Considered

**Option**: Keep using `prec_prob` for protein scoring, only use `nonMBR_prob` for Step 10 PSM q-values

**Rejected Because**:
- Creates inconsistency: PSM FDR based on nonMBR scores, protein FDR based on MBR scores
- Protein q-values would not be directly comparable to PSM q-values
- Harder to interpret results when different scoring methods used at different levels
- Current approach maintains full consistency: all FDR calculations use the same underlying scores

---

## Files Modified

1. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/protein_inference_pipeline.jl` (lines 264, 273)
2. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl` (lines 1884, 1903-1904)

## Expected Behavior

After this change:
- Protein groups will have slightly different `pg_score` values (based on pre-MBR peptide probabilities)
- CV fold assignments may differ for some proteins (if best peptide changes based on different scoring)
- Protein q-values will be calculated from nonMBR-based scores
- Complete consistency across PSM-level and protein-level FDR calculations
