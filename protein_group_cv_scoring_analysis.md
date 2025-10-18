# Protein Group Cross-Validation Scoring Analysis

## Executive Summary

The `perform_probit_analysis_multifold` function implements protein-level probit regression with cross-validation (CV) to refine protein group scores. This document analyzes the complete workflow from CV fold assignment through model training and scoring.

**Key Finding**: The CV implementation appears **CORRECT** - proteins are assigned to CV folds based on their highest-scoring peptide, and probit models are properly trained on out-of-fold data and applied to held-out proteins.

---

## Table of Contents
1. [Overview](#overview)
2. [Complete Workflow](#complete-workflow)
3. [CV Fold Assignment Logic](#cv-fold-assignment-logic)
4. [Training and Scoring Process](#training-and-scoring-process)
5. [Data Flow Diagram](#data-flow-diagram)
6. [Critical Analysis](#critical-analysis)
7. [Potential Issues and Edge Cases](#potential-issues-and-edge-cases)
8. [Recommendations](#recommendations)

---

## Overview

### Purpose
Protein-level probit regression refines initial protein group scores (log-sum of peptide probabilities) using additional features like:
- `:pg_score` - Initial protein group score
- `:any_common_peps` - Whether protein has any fully-cleaved, unmodified peptides
- `:peptide_coverage` - Fraction of possible peptides observed (commented out in current version)
- `:n_possible_peptides` - Total peptides in library (commented out in current version)

### Key Files
- **Main function**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl:1888-2014`
- **Helper functions**: Lines 1654-1794 in same file
- **Library CV folds**: `src/structs/LibraryIon.jl:743`

---

## Complete Workflow

### Step-by-Step Process

#### **Phase 1: CV Fold Detection (Lines 1901-1903)**
```julia
unique_cv_folds = detect_unique_cv_folds(precursors)
```
- Extracts unique CV fold values from `precursors.pid_to_cv_fold`
- CV folds are assigned at **library creation time** (before any data analysis)
- Typically returns `[0, 1]` for 2-fold CV

**Source of CV Folds**:
- CV folds are stored in the spectral library and assigned to **precursors** (not proteins)
- Each precursor has a fixed CV fold assignment from the library

---

#### **Phase 2: Build Protein-to-CV-Fold Mapping (Lines 1906-1910)**

```julia
if protein_to_cv_fold === nothing
    psm_paths = [get_corresponding_psm_path(ref) for ref in pg_refs]
    protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, precursors)
end
```

##### Function: `build_protein_cv_fold_mapping` (Lines 1700-1756)

**Algorithm**:
1. Iterate through all PSM files
2. For each protein group:
   - Find all peptides belonging to the protein
   - Identify the **highest-scoring peptide** (by `prec_prob`)
   - Look up that peptide's `cv_fold` from the library
   - Assign the protein to that CV fold
3. If the same protein appears in multiple files:
   - Keep track of the globally highest-scoring peptide
   - Update protein's CV fold if a higher-scoring peptide is found

**Key Code (Lines 1734-1750)**:
```julia
for group in groupby(psms, :inferred_protein_group)
    protein_name = first(group.inferred_protein_group)

    # Find highest scoring PSM
    best_idx = argmax(group.prec_prob)
    best_score = group.prec_prob[best_idx]
    precursor_idx = group.precursor_idx[best_idx]

    # Get cv_fold from library (more reliable than PSM file)
    cv_fold = getCvFold(precursors, precursor_idx)

    # Update if this is the best score for this protein
    value = (best_score = best_score, cv_fold = cv_fold)
    if !haskey(protein_to_cv_fold, protein_name)
        insert!(protein_to_cv_fold, protein_name, value)
    elseif best_score > protein_to_cv_fold[protein_name].best_score
        protein_to_cv_fold[protein_name] = value
    end
end
```

**Result**:
- Dictionary mapping: `protein_name → (best_score, cv_fold)`
- Each protein assigned to exactly **ONE** CV fold
- Assignment based on highest-scoring peptide across all files

---

#### **Phase 3: Assign CV Folds to Protein Groups (Lines 1912-1913)**

```julia
assign_protein_group_cv_folds!(all_protein_groups, protein_to_cv_fold)
```

##### Function: `assign_protein_group_cv_folds!` (Lines 1771-1794)

**Algorithm**:
1. For each protein in `all_protein_groups` DataFrame:
   - Look up protein in the `protein_to_cv_fold` dictionary
   - Assign the CV fold to the protein
   - If protein not found in dictionary → assign to fold 0 (default)

**Key Code (Lines 1776-1786)**:
```julia
for (i, protein_name) in enumerate(all_protein_groups.protein_name)
    if haskey(protein_to_cv_fold, protein_name)
        cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
    else
        missing_count += 1
        cv_folds[i] = UInt8(0)  # Default to fold 0
    end
end
```

**Result**:
- `all_protein_groups` DataFrame now has `:cv_fold` column
- Each protein group instance (per file) assigned to the same global CV fold

---

#### **Phase 4: Feature Preparation (Lines 1921-1933)**

```julia
feature_names = [:pg_score, :any_common_peps]

# Apply feature filtering
adjust_any_common_peps!(feature_names, all_protein_groups)
remove_zero_variance_columns!(feature_names, all_protein_groups)
```

**Feature Filtering**:
1. **`adjust_any_common_peps!`** (Lines 100-108):
   - Checks if `:any_common_peps` is constant (all true or all false)
   - Removes it if constant to avoid singular matrix errors

2. **`remove_zero_variance_columns!`**:
   - Removes any features with zero variance

**Result**: Filtered feature list for probit regression

---

#### **Phase 5: Training Set Selection with FDR Control (Lines 1938-1944)**

```julia
# Determine positive training examples using 1% FDR on pg_score
qvals = Vector{Float32}(undef, n_proteins)
get_qvalues!(all_protein_groups.pg_score, all_protein_groups.target, qvals)
passing_mask = (qvals .<= 0.01f0) .& all_protein_groups.target
train_mask_FDR = passing_mask .| .!all_protein_groups.target
```

**Training Set Composition**:
- **Positive examples**: Target proteins passing 1% FDR threshold on initial `pg_score`
- **Negative examples**: ALL decoy proteins
- This is a **two-class** problem (high-confidence targets vs. decoys)

**Rationale**:
- Using only high-confidence targets (1% FDR) reduces label noise
- Similar to "semi-supervised" learning - learn from confident positives
- Prevents training on likely false positives

---

#### **Phase 6: Cross-Validation Model Training (Lines 1946-1965)**

```julia
for test_fold in unique_cv_folds
    # Get training data (all folds except test_fold)
    train_mask = (all_protein_groups.cv_fold .!= test_fold) .& train_mask_FDR

    # Check if we have sufficient data
    n_train_targets = sum(all_protein_groups[train_mask, :target])
    n_train_decoys = sum(.!all_protein_groups[train_mask, :target])

    if n_train_targets < 10 || n_train_decoys < 10
        @user_warn "Insufficient training data for fold $test_fold"
        continue
    end

    X_train = Matrix{Float64}(all_protein_groups[train_mask, feature_names])
    y_train = all_protein_groups[train_mask, :target]

    # Fit model
    β_fitted = fit_probit_model(X_train, y_train)
    models[test_fold] = β_fitted
end
```

**CV Training Process** (for 2-fold CV):

**Fold 0 as Test Set**:
- **Training data**: All proteins in fold 1 + FDR filtering
- **Model**: `models[0] = β_fitted`
- **Test data**: All proteins in fold 0 (scored later)

**Fold 1 as Test Set**:
- **Training data**: All proteins in fold 0 + FDR filtering
- **Model**: `models[1] = β_fitted`
- **Test data**: All proteins in fold 1 (scored later)

**Key Property**: Each model is trained **WITHOUT** seeing its corresponding test fold proteins

---

#### **Phase 7: Apply Models to Test Folds (Lines 1971-1987)**

```julia
all_protein_groups[!, :old_pg_score] = copy(all_protein_groups[!, :pg_score])

for test_fold in unique_cv_folds
    if !haskey(models, test_fold)
        @user_warn "No model available for fold $test_fold, keeping original scores"
        continue
    end

    test_mask = all_protein_groups.cv_fold .== test_fold
    n_test = sum(test_mask)

    if n_test > 0
        X_test = Matrix{Float64}(all_protein_groups[test_mask, feature_names])
        prob_scores = calculate_probit_scores(X_test, models[test_fold])
        all_protein_groups[test_mask, :pg_score] = Float32.(prob_scores)
    end
end
```

**Scoring Process**:
- For each CV fold:
  - Select proteins in that fold (`test_mask`)
  - Apply the model trained on **OTHER** folds
  - Update `:pg_score` column with probit-refined scores

**Result**:
- All proteins receive probit-refined scores
- Each protein scored by a model that **never** saw it during training

---

#### **Phase 8: Update Protein Group Files (Lines 2008-2010)**

```julia
if !isempty(pg_refs)
    apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names)
end
```

##### Function: `apply_probit_scores_multifold!` (Lines 1810-1852)

**Algorithm**:
1. For each protein group file:
   - Assign CV folds using the pre-built mapping
   - For each protein:
     - Determine its CV fold
     - Apply the corresponding model
     - Update `:pg_score` in the file

**Key Code (Lines 1820-1842)**:
```julia
for (i, protein_name) in enumerate(df.protein_name)
    if haskey(protein_to_cv_fold, protein_name)
        cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
    else
        cv_folds[i] = UInt8(0)
    end
end

for fold in unique(cv_folds)
    if !haskey(models, fold)
        continue
    end

    fold_mask = cv_folds .== fold
    if any(fold_mask)
        X = Matrix{Float64}(df[fold_mask, feature_names])
        scores = calculate_probit_scores(X, models[fold])
        df[fold_mask, :pg_score] = Float32.(scores)
    end
end
```

**Result**: Individual protein group files updated with probit scores

---

## CV Fold Assignment Logic

### Hierarchy of CV Fold Assignment

```
1. Library Precursors (Source of Truth)
   ├─ Each precursor has cv_fold assigned at library creation
   └─ Stored in: precursors.pid_to_cv_fold

2. PSMs (Inherit from Precursors)
   ├─ cv_fold copied from library during PSM creation
   └─ Column: :cv_fold

3. Proteins (Derived from Best Peptide)
   ├─ CV fold = CV fold of highest-scoring peptide
   ├─ Built using: build_protein_cv_fold_mapping()
   └─ Stored in: protein_to_cv_fold dictionary

4. Protein Groups (Assigned from Mapping)
   ├─ Each protein group assigned using protein_to_cv_fold
   └─ Column: :cv_fold
```

### Why This Design?

**Problem**: Proteins span multiple peptides, peptides may be in different CV folds

**Solution**: Assign protein to CV fold of its **most confident peptide**

**Rationale**:
1. Highest-scoring peptide is most likely to be correct
2. Protein's presence is most reliably indicated by this peptide
3. Consistent assignment across all files (global mapping)

**Example**:
```
Protein P1 has 3 peptides:
- Peptide A: prec_prob=0.95, cv_fold=0
- Peptide B: prec_prob=0.70, cv_fold=1
- Peptide C: prec_prob=0.60, cv_fold=0

→ Protein P1 assigned to cv_fold=0 (based on Peptide A)
```

---

## Training and Scoring Process

### Training Data Composition

For **Fold 0** as test set (Fold 1 trains the model):

```
Training Set = Proteins in Fold 1 meeting FDR criteria

Positive Examples:
- Proteins in Fold 1
- AND target = true
- AND q-value ≤ 0.01 on initial pg_score

Negative Examples:
- Proteins in Fold 1
- AND target = false (ALL decoys in Fold 1)
```

### Test Data

```
Test Set = ALL proteins in Fold 0
- Includes both targets and decoys
- No FDR filtering
- Model produces refined pg_score for all
```

### Model Application

Each protein scored by model trained on **opposite fold**:
- Fold 0 proteins → scored by model trained on Fold 1
- Fold 1 proteins → scored by model trained on Fold 0

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                     SPECTRAL LIBRARY                            │
│  (Created before any analysis)                                  │
│                                                                  │
│  precursors.pid_to_cv_fold: [0, 1, 0, 1, 0, ...]               │
│  ↓ Each precursor assigned to CV fold                          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                     PSM FILES (per MS run)                      │
│                                                                  │
│  PSMs inherit cv_fold from library via precursor_idx           │
│  Column added during second pass search                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              BUILD PROTEIN CV FOLD MAPPING                      │
│  (build_protein_cv_fold_mapping)                               │
│                                                                  │
│  For each protein:                                              │
│    1. Find all PSMs belonging to protein                        │
│    2. Identify highest-scoring PSM (by prec_prob)              │
│    3. Look up cv_fold of that PSM's precursor                  │
│    4. Assign protein to that cv_fold                           │
│                                                                  │
│  Result: protein_to_cv_fold dictionary                         │
│    "Protein1" → (best_score=0.95, cv_fold=0)                  │
│    "Protein2" → (best_score=0.88, cv_fold=1)                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│           PROTEIN GROUPS (aggregated across files)              │
│  (all_protein_groups DataFrame)                                │
│                                                                  │
│  Initial columns:                                               │
│    - protein_name, pg_score, n_peptides, etc.                  │
│                                                                  │
│  CV fold assignment:                                            │
│    - Look up protein in protein_to_cv_fold                     │
│    - Add :cv_fold column                                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│               TRAINING SET SELECTION                            │
│  (with 1% FDR filtering)                                       │
│                                                                  │
│  Positive Examples: Targets passing q ≤ 0.01                   │
│  Negative Examples: ALL decoys                                  │
│                                                                  │
│  Training mask per fold:                                        │
│    - Fold 0 model: Train on Fold 1 proteins + FDR filter      │
│    - Fold 1 model: Train on Fold 0 proteins + FDR filter      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  PROBIT MODEL TRAINING                          │
│  (fit_probit_model)                                            │
│                                                                  │
│  For each CV fold:                                              │
│    1. Extract training data (opposite fold + FDR)              │
│    2. Build feature matrix X                                    │
│    3. Fit probit regression: Φ^(-1)(P(target=1|X)) = Xβ       │
│    4. Store model coefficients β                               │
│                                                                  │
│  Result: models dictionary                                      │
│    models[0] = β_fitted (trained on Fold 1)                   │
│    models[1] = β_fitted (trained on Fold 0)                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  APPLY MODELS TO TEST FOLDS                     │
│  (calculate_probit_scores)                                     │
│                                                                  │
│  For each CV fold:                                              │
│    1. Select proteins in that fold (test set)                  │
│    2. Apply model trained on OTHER fold                        │
│    3. Compute: pg_score_new = Φ(Xβ)                           │
│    4. Update :pg_score column                                  │
│                                                                  │
│  Fold 0 proteins ← scored by models[0] (trained on Fold 1)    │
│  Fold 1 proteins ← scored by models[1] (trained on Fold 0)    │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│            UPDATE PROTEIN GROUP FILES                           │
│  (apply_probit_scores_multifold!)                             │
│                                                                  │
│  For each protein group file:                                   │
│    1. Assign CV folds using protein_to_cv_fold                │
│    2. Apply corresponding model to each protein                │
│    3. Write updated pg_score to file                           │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    DOWNSTREAM ANALYSIS                          │
│  - Calculate global q-values                                    │
│  - Filter by FDR threshold                                      │
│  - Update PSMs with protein scores                             │
└─────────────────────────────────────────────────────────────────┘
```

---

## Critical Analysis

### ✅ CORRECT Implementations

1. **No Data Leakage**:
   - Training data for fold i contains **ONLY** proteins from fold j (j ≠ i)
   - Test proteins never seen during model training
   - CV fold assignments are **pre-determined** from library, not data-driven

2. **Consistent Protein Assignment**:
   - Global mapping ensures same protein gets same CV fold across all files
   - Prevents same protein appearing in both train and test sets

3. **FDR-Based Training Selection**:
   - Using 1% FDR threshold creates high-confidence positive training set
   - Reduces label noise from false positive targets
   - Standard practice in semi-supervised learning

4. **Feature Filtering**:
   - Removes constant features to prevent singular matrices
   - Appropriate for probit regression

5. **Model Persistence**:
   - Models stored in dictionary for later application
   - Enables consistent scoring across protein group files

### 🔍 Design Decisions (Not Errors, But Worth Noting)

1. **Best Peptide Selection**:
   - **Current**: Protein CV fold = CV fold of highest-scoring peptide
   - **Alternative**: Could use majority vote of all peptides
   - **Justification**: Best peptide is most reliable indicator

2. **Feature Set**:
   - **Current**: Only `:pg_score` and `:any_common_peps` active
   - `:peptide_coverage` and `:n_possible_peptides` commented out
   - **Implication**: Limited feature set may reduce model complexity

3. **Training Set Size**:
   - Only targets passing 1% FDR used as positives
   - Could potentially use softer thresholds or all targets
   - **Trade-off**: Precision vs. training set size

---

## Potential Issues and Edge Cases

### 1. **Proteins Without Matching Peptides**

**Issue**: Line 1784-1786 assigns proteins to fold 0 if not found in mapping
```julia
if haskey(protein_to_cv_fold, protein_name)
    cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
else
    missing_count += 1
    cv_folds[i] = UInt8(0)  # Default to fold 0
end
```

**When This Happens**:
- Protein appears in protein group file but not in any PSM file
- Possible if protein filtering differs between stages

**Impact**:
- These proteins always assigned to fold 0
- Could create imbalanced fold sizes
- Minimal impact if rare

**Recommendation**: Log these cases for investigation

---

### 2. **Insufficient Training Data**

**Issue**: Lines 1954-1957 check for minimum 10 targets and 10 decoys
```julia
if n_train_targets < 10 || n_train_decoys < 10
    @user_warn "Insufficient training data for fold $test_fold"
    continue
end
```

**When This Happens**:
- Very small datasets
- Highly imbalanced CV folds
- Severe FDR filtering (almost no targets pass 1%)

**Impact**:
- Model not trained for that fold
- Proteins in that fold keep original scores (line 1974-1976)

**Recommendation**: Consider lower thresholds or merging folds

---

### 3. **Feature Availability**

**Issue**: Features like `:peptide_coverage` require protein catalog

**Current State**:
- `:peptide_coverage` and `:n_possible_peptides` commented out
- Only `:pg_score` and `:any_common_peps` active

**Impact**:
- Simpler model may have less discriminative power
- Faster computation, less potential overfitting

**Recommendation**:
- Document why features are disabled
- Consider A/B testing with/without these features

---

### 4. **Protein Appears in Multiple Runs with Different Scores**

**Scenario**:
```
Run 1: Protein P1, peptide A (score=0.95, cv_fold=0)
Run 2: Protein P1, peptide B (score=0.98, cv_fold=1)
```

**Current Behavior**:
- Protein P1 assigned to cv_fold=1 (higher score from peptide B)
- All instances of P1 across both runs get cv_fold=1

**Is This Correct?**:
- ✅ YES - ensures consistency
- Global assignment prevents P1 from being in both train and test
- Uses most confident detection to determine fold

---

### 5. **CV Fold Balance**

**Issue**: No guarantee that CV folds are balanced

**Factors Affecting Balance**:
- Library design determines precursor CV fold assignments
- Protein assignments based on best peptide scores
- Could create 80/20 split instead of 50/50

**Impact**:
- Imbalanced folds → smaller training sets
- Less representative models

**Recommendation**:
- Monitor fold distribution (currently done at line 1916-1919)
- Consider stratified splitting at library creation time

---

## Recommendations

### 1. **Add Diagnostics**

**Suggested Additions**:
```julia
# After line 1919 (fold distribution check)
n_total = nrow(all_protein_groups)
for fold in unique_cv_folds
    n_fold = get(fold_counts, fold, 0)
    pct = round(100.0 * n_fold / n_total, digits=1)
    @user_info "CV Fold $fold: $n_fold proteins ($pct%)"
end

# Report training set sizes
for test_fold in unique_cv_folds
    train_mask = (all_protein_groups.cv_fold .!= test_fold) .& train_mask_FDR
    n_train_targets = sum(all_protein_groups[train_mask, :target])
    n_train_decoys = sum(.!all_protein_groups[train_mask, :target])
    @user_info "Fold $test_fold training set: $n_train_targets targets, $n_train_decoys decoys"
end
```

### 2. **Validate Protein Mapping Coverage**

**Check**:
```julia
# After building protein_to_cv_fold mapping
n_proteins_in_groups = length(unique(all_protein_groups.protein_name))
n_proteins_in_mapping = length(protein_to_cv_fold)

if n_proteins_in_groups > n_proteins_in_mapping
    @user_warn "$(n_proteins_in_groups - n_proteins_in_mapping) proteins lack CV fold mapping"
end
```

### 3. **Document Feature Selection**

Add comments explaining why features are enabled/disabled:
```julia
# Current active features (2025-01)
feature_names = [
    :pg_score,          # Initial protein score (log-sum of peptide probs)
    :any_common_peps    # Has fully-cleaved, unmodified peptides
]

# Disabled features (require protein catalog, may cause overfitting)
# :peptide_coverage      # Fraction of possible peptides observed
# :n_possible_peptides   # Total peptides in library
```

### 4. **Consider Fallback for Small Datasets**

```julia
# If insufficient data for CV, fall back to single-fold training
if any([sum(all_protein_groups.cv_fold .== fold) < 20 for fold in unique_cv_folds])
    @user_warn "Small dataset detected, using pooled training instead of CV"
    # Train single model on all data with FDR filtering
    # Apply to all proteins
end
```

### 5. **Monitor Model Performance**

Track improvement metrics more granularly:
```julia
# Per-fold improvement
for fold in unique_cv_folds
    fold_mask = all_protein_groups.cv_fold .== fold
    old_passing = sum((old_qvalues[fold_mask] .<= 0.01f0) .& all_protein_groups.target[fold_mask])
    new_passing = sum((new_qvalues[fold_mask] .<= 0.01f0) .& all_protein_groups.target[fold_mask])
    @user_info "Fold $fold: $old_passing → $new_passing targets at 1% FDR"
end
```

---

## Conclusion

### Summary

The `perform_probit_analysis_multifold` implementation is **fundamentally correct**:

1. ✅ **No data leakage**: Proteins properly partitioned into train/test splits
2. ✅ **Consistent CV folds**: Global mapping prevents same protein in multiple folds
3. ✅ **Proper cross-validation**: Each protein scored by model that never saw it
4. ✅ **FDR-based training**: High-quality positive training examples
5. ✅ **Feature filtering**: Prevents singular matrix errors

### Key Insights

**CV Fold Assignment**:
- Determined by highest-scoring peptide's CV fold
- Consistent across all files (global mapping)
- Based on library CV folds (pre-assigned, not data-driven)

**Training Strategy**:
- Semi-supervised approach (confident positives + all negatives)
- 1% FDR threshold for positive selection
- Prevents training on likely false positives

**Model Application**:
- Each fold scored by model trained on other fold(s)
- Maintains proper train/test separation
- Updates both in-memory DataFrame and individual files

### Future Work

Consider testing:
1. Alternative protein CV fold assignment strategies (e.g., majority vote)
2. Different FDR thresholds for training set selection
3. Re-enabling `:peptide_coverage` and `:n_possible_peptides` features
4. Stratified CV fold assignment at library creation

---

**Document Created**: 2025-01-18
**Analysis Based On**: Pioneer.jl commit 437b097f
**Analyzed Function**: `perform_probit_analysis_multifold` (utils.jl:1888-2014)
