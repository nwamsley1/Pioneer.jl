# Cross-Validation Fold Assignment Analysis and Plan

## Executive Summary

Pioneer.jl's current CV fold assignment creates data leakage when Match-Between-Runs (MBR) is enabled, leading to biased ML model performance and potentially explaining discrepancies in entrapment analysis. This document analyzes the current algorithm, identifies the core problem, and proposes a solution.

## Current CV Fold Assignment Algorithm

### Implementation Location
- **File**: `src/structs/LibraryIon.jl`
- **Functions**: `StandardLibraryPrecursors()` and `PlexedLibraryPrecursors()` constructors
- **Lines**: ~683-700

### Current Algorithm (Both Current Branch and Commit 02af558)

```julia
# Step 1: Group unique proteins
unique_proteins = unique(accession_numbers)

# Step 2: Random assignment by protein group
pg_to_cv_fold = Dictionary{String, UInt8}()
cv_folds = UInt8[0, 1]
Random.seed!(1776)  # Fixed seed for reproducibility
for pg in unique_proteins
    insert!(pg_to_cv_fold, pg, rand(cv_folds))  # Random 50/50 assignment
end

# Step 3: Assign all peptides from same protein to same fold
pid_to_cv_fold = Vector{UInt8}(undef, n)
for pid in range(1, n)
    pid_to_cv_fold[pid] = pg_to_cv_fold[accession_numbers[pid]]
end
```

### Key Characteristics
1. **Protein-based grouping**: All precursors from the same protein go to the same CV fold
2. **Random assignment**: Each protein group has 50% chance of being in fold 0 or 1
3. **Fixed seed**: Uses `Random.seed!(1776)` for reproducibility
4. **No balancing**: No consideration of target/decoy balance or pair relationships

## Target-Decoy Pairing Mechanism

### Implementation
- **Key function**: `getPartnerPrecursorIdx(precursors)[pid]`
- **Data structure**: `partner_precursor_idx` column in library
- **Creation**: Done during library building in `pair_decoys.jl`

### How Pairing Works
```julia
# During library building (pair_decoys.jl)
# 1. For each target peptide sequence
# 2. Create reversed decoy sequence  
# 3. Link target ↔ decoy via partner_precursor_idx
# 4. Assign unique pair_id to both members

# Example:
target_idx = 1234    # PEPTIDE
decoy_idx = 5678     # EDITPEP
partner_precursor_idx[1234] = 5678  # target points to decoy
partner_precursor_idx[5678] = 1234  # decoy points to target
pair_id[1234] = pair_id[5678] = 42  # same pair ID
```

### MBR Usage of Pairs
From FirstPassSearch.jl (commit 02af558):
```julia
if params.match_between_runs==true
    for (pid, val) in pairs(precursor_dict)
        partner_pid = getPartnerPrecursorIdx(precursors)[pid]
        if !ismissing(partner_pid)
            # If partner not identified, add it with same RT
            if !haskey(precursor_dict, partner_pid)
                insert!(precursor_dict, partner_pid, val)
                setPredIrt!(search_context, partner_pid, getIrt(precursors)[pid])
            end
        end
    end
end
```

**Key insight**: When MBR is enabled, if one member of a target-decoy pair is identified, its partner is automatically added and shares RT/scoring information.

## The Data Leakage Problem

### Scenario Analysis
Consider a target-decoy pair where:
- Target precursor (idx=1234, protein=P1) → CV fold 0
- Decoy precursor (idx=5678, protein=P2) → CV fold 1

### When Training Fold 0 (Testing Fold 1):
1. **Training set** includes target precursor 1234
2. **Test set** includes decoy precursor 5678
3. **MBR processing**: If decoy 5678 is identified in test set, target 1234 gets added to scoring
4. **Result**: Model trains on target while being evaluated on its paired decoy

### The Leakage Mechanism
- **Information flow**: Test fold decoy identification → training fold target addition
- **Shared features**: RT prediction, intensity patterns, scoring features
- **Biased evaluation**: Model sees "easy" targets in training that correspond to test decoys

## Evidence of Current Problems

### From Entrapment Analysis
- **Observation**: Entrapment sequences (fake targets) have different error rates than regular decoys
- **Possible cause**: ML model learns to distinguish entrapment-decoy pairs vs. regular target-decoy pairs
- **CV contribution**: Unbalanced folds may train model differently on different types of pairs

### From Recent Log Output Addition
The new logging shows fold composition:
```
[ Info: Fold 0 - Train: 3750 targets, 3750 decoys | Test: 1250 targets, 1250 decoys
[ Info: Fold 1 - Train: 3750 targets, 3750 decoys | Test: 1250 targets, 1250 decoys
```

**Question**: Are paired targets and decoys in the same fold? Current algorithm doesn't guarantee this.

## Mathematical Analysis

### Current Expected Imbalance
With random protein assignment:
- **Probability both members of pair in same fold**: P(same fold) = P(both fold 0) + P(both fold 1) = 0.5² + 0.5² = 0.5
- **Probability pair split across folds**: P(split) = 1 - 0.5 = 0.5

**Expected result**: ~50% of target-decoy pairs are split across CV folds when MBR is enabled.

### Impact on Training
For a dataset with N target-decoy pairs:
- **Leaky pairs**: ~N/2 pairs have information flow between train/test
- **Clean pairs**: ~N/2 pairs maintain proper CV separation
- **Bias magnitude**: Depends on how much MBR information transfer affects scoring

## Proposed Solution: Pair-Aware CV Assignment

### Algorithm Overview
Instead of assigning individual proteins to folds, assign target-decoy **pairs** to folds.

### Implementation Strategy

#### Option 1: Pair-First Assignment
```julia
# Step 1: Identify all target-decoy pairs
pair_to_proteins = Dictionary{UInt32, Vector{String}}()
for pid in 1:n
    pair_id = get_pair_id(precursor_table, pid)
    protein = accession_numbers[pid]
    if !haskey(pair_to_proteins, pair_id)
        pair_to_proteins[pair_id] = String[]
    end
    push!(pair_to_proteins[pair_id], protein)
end

# Step 2: Assign each pair to a fold
pair_to_cv_fold = Dictionary{UInt32, UInt8}()
Random.seed!(1776)
for pair_id in keys(pair_to_proteins)
    pair_to_cv_fold[pair_id] = rand(cv_folds)
end

# Step 3: Assign all proteins in each pair to same fold
pg_to_cv_fold = Dictionary{String, UInt8}()
for (pair_id, proteins) in pairs(pair_to_proteins)
    fold = pair_to_cv_fold[pair_id]
    for protein in proteins
        pg_to_cv_fold[protein] = fold
    end
end
```

#### Option 2: Protein-Pair Hybrid
```julia
# Step 1: Start with current protein assignment
# Step 2: For each protein, check if any of its precursors have partners
# Step 3: If partner protein in different fold, move it to same fold
# Step 4: Continue until no conflicts remain
```

### Validation
After assignment, verify:
1. **Pair integrity**: All members of each target-decoy pair are in same fold
2. **Balance preservation**: Fold sizes remain approximately equal
3. **Protein grouping**: All precursors from same protein still in same fold

## Implementation Plan

### Phase 1: Analysis
1. **Current state audit**: Count how many pairs are currently split
2. **Impact measurement**: Quantify performance difference with/without pair splits
3. **Data structure review**: Ensure `pair_id` and `partner_precursor_idx` are available during library construction

### Phase 2: Algorithm Development
1. **Prototype pair-aware assignment** in separate function
2. **Test on example dataset** to verify pair integrity
3. **Compare fold balance** with current method

### Phase 3: Integration
1. **Modify LibraryIon.jl constructors** to use pair-aware assignment
2. **Add validation checks** to ensure correctness
3. **Update logging** to show pair statistics per fold

### Phase 4: Validation
1. **Re-run entrapment analysis** with fixed CV assignment
2. **Compare ML performance** before/after fix
3. **Validate that target/decoy error rates** are now consistent

## Expected Outcomes

### Immediate Benefits
1. **Proper CV independence**: No information leakage between train/test folds
2. **Fair model evaluation**: More accurate cross-validation performance estimates
3. **Consistent entrapment analysis**: Entrapment and regular decoys treated equally

### Long-term Benefits
1. **More robust ML models**: Training on truly independent data
2. **Better generalization**: Models that perform consistently on new data
3. **Improved scientific validity**: Results that properly reflect model capabilities

## Risk Assessment

### Low Risk
- **Algorithm complexity**: Modification is straightforward
- **Performance impact**: Minimal computational overhead
- **Backward compatibility**: Results will be different but more correct

### Medium Risk
- **Fold imbalance**: Need to verify balanced fold sizes are maintained
- **Edge cases**: Handle unpaired precursors, missing partner information

### Mitigation Strategies
1. **Gradual rollout**: Test on small datasets first
2. **Comparison studies**: Run both algorithms in parallel initially
3. **Comprehensive validation**: Check multiple metrics and datasets

## Conclusion

The current CV fold assignment creates systematic data leakage when MBR is enabled, potentially explaining observed biases in entrapment analysis and ML model evaluation. The proposed pair-aware assignment ensures proper cross-validation independence while maintaining the protein-grouping principle that prevents overfitting.

This fix addresses a fundamental methodological issue that affects the scientific validity of Pioneer.jl's machine learning components and should be prioritized for implementation.