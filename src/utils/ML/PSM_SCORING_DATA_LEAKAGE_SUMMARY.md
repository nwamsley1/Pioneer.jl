# PSM Scoring and Data Leakage Investigation Summary

## Executive Summary

This document summarizes the investigation into percolator-style PSM (Peptide-Spectrum Match) scoring in Pioneer.jl and the discovery of potential data leakage issues when Match-Between-Runs (MBR) is enabled. The investigation revealed that the current cross-validation (CV) fold assignment strategy at the protein-group level can lead to biased machine learning model performance, particularly affecting entrapment analysis.

## Key Findings

### 1. Cross-Validation Data Leakage with MBR

**The Problem:**
When Match-Between-Runs (MBR) is enabled in FirstPassSearch, if one member of a target-decoy pair is identified, its partner is automatically added to the scoring dictionary with shared retention time (RT) and scoring information. This creates data leakage when target-decoy pairs are split across CV folds.

**Evidence from Code (FirstPassSearch.jl):**
```julia
if params.match_between_runs==true
    for (pid, val) in pairs(precursor_dict)
        partner_pid = getPartnerPrecursorIdx(precursors)[pid]
        if !ismissing(partner_pid)
            if !haskey(precursor_dict, partner_pid)
                insert!(precursor_dict, partner_pid, val)
                setPredIrt!(search_context, partner_pid, getIrt(precursors)[pid])
            end
        end
    end
end
```

**Impact:**
- Information flows from test fold to training fold through partner precursors
- ML models train on targets while being evaluated on their paired decoys
- Approximately 50% of target-decoy pairs are expected to be split across folds with random protein assignment

### 2. CV Fold Assignment Architecture

**Current Implementation (Protein-Group Based):**
- CV folds are assigned at library construction time in `LibraryIon.jl`
- All precursors from the same protein are assigned to the same CV fold
- Random 50/50 assignment per protein group with fixed seed (1776)
- No consideration of target-decoy pair relationships

**Key Insight:**
Each protein group contains both target and decoy versions of precursors. The current assignment ensures that target and decoy versions of the same protein sequence stay together, but doesn't guarantee that synthetic decoy proteins (with different accession numbers) stay with their target counterparts.

### 3. Feature Bias in MBR-Related Features

**Removed Features from ML Training:**
Based on the investigation, three potentially biased MBR features were identified and removed from training:
- `MBR_rv_coefficient` - May leak information about partner quality
- `MBR_best_irt_diff` - Could reveal partner RT predictions
- `MBR_is_best_decoy` - Direct indicator of target-decoy relationships

**Retained MBR Features:**
- `MBR_num_runs` - Number of runs where precursor was identified
- `MBR_max_pair_prob` - Maximum probability across runs
- `MBR_log2_weight_ratio` - Intensity ratio information
- `MBR_log2_explained_ratio` - Explained variance ratio

## Implemented Improvements

### 1. Enhanced Diagnostic Logging

Added comprehensive target/decoy/entrapment counting at multiple levels:

```julia
# Dataset-level counts
n_targets = sum(psms.target)
n_decoys = sum(.!psms.target)
n_entrapments = hasproperty(psms, :entrapment) ? sum(psms.entrapment) : 0

# Per-CV-fold counts
for fold in [0, 1]
    train_mask = psms.cv_fold .!= fold
    test_mask = psms.cv_fold .== fold
    # Count targets, decoys, entrapments in train/test sets
end
```

This helps diagnose:
- Overall dataset composition
- CV fold balance
- Potential bias in entrapment vs. regular sequences

### 2. Feature Importance Visualization

Fixed and enhanced feature importance printing using EvoTrees API:

```julia
# Correct API usage
importance_dict = EvoTrees.importance(bst)

# Display top features, 10 per line
sorted_features = sort(collect(importance_dict), by=x->x[2], rev=true)
for i in 1:10:length(sorted_features)
    # Format and display features
end
```

Now shows feature importance for:
- All iterations of training
- All CV folds
- Both gain and frequency metrics

### 3. Pair-Based CV Assignment (Investigated but Reverted)

**Attempted Solution:**
Implemented pair-based CV fold assignment using `pair_id` to ensure target-decoy pairs stay together:

```julia
# Group by pair_id instead of protein
unique_pairs = unique(pair_ids)
pair_to_cv_fold = Dictionary{UInt32, UInt8}()
for pair in unique_pairs
    insert!(pair_to_cv_fold, pair, rand(cv_folds))
end
```

**Outcome:**
- Successfully kept target-decoy pairs together
- However, reverted to protein-group assignment per user request
- The protein-group approach may be preferred for biological relevance

## Implications for Entrapment Analysis

### The Entrapment Problem

Entrapment sequences (synthetic "fake" targets) showed different error rates than regular decoys, potentially due to:

1. **ML Model Bias:** Models may learn to distinguish entrapment-decoy pairs from regular target-decoy pairs
2. **CV Fold Imbalance:** Unbalanced distribution of entrapment sequences across folds
3. **MBR Information Leakage:** Partner addition creating systematic biases

### Expected Impact of Fixes

The implemented improvements should:
- Provide better visibility into fold composition
- Remove potentially biased features from training
- Enable diagnosis of entrapment-specific issues
- Improve model generalization

## Mathematical Analysis

### Probability of Data Leakage

With random protein assignment and MBR enabled:
- P(pair in same fold) = 0.5² + 0.5² = 0.5
- P(pair split across folds) = 0.5
- Expected leaky pairs = N/2 (where N = total pairs)

### Impact on Model Training

For split pairs with MBR:
- Training set sees target with features
- Test set sees decoy with shared RT/scoring
- Model learns correlation between shared features

## Recommendations

### Short-term (Implemented)
1. ✅ Remove biased MBR features from training
2. ✅ Add comprehensive diagnostic logging
3. ✅ Fix feature importance display

### Medium-term (Proposed)
1. Consider pair-aware CV assignment for MBR scenarios
2. Add validation checks for fold balance
3. Implement stratified sampling for entrapment sequences

### Long-term (Future Work)
1. Develop MBR-aware CV strategy that prevents leakage
2. Investigate alternative protein grouping strategies
3. Consider separate models for MBR vs. non-MBR scenarios

## Technical Details

### Key Files Modified

1. **percolatorSortOf.jl**
   - Added target/decoy/entrapment counting
   - Fixed feature importance API usage
   - Enhanced diagnostic output

2. **score_psms.jl**
   - Removed biased MBR features from training
   - Maintained feature calculation for diagnostics

3. **LibraryIon.jl**
   - Investigated CV fold assignment logic
   - Temporarily implemented pair-based assignment
   - Reverted to protein-group assignment

### Commits

- Initial feature importance fix
- Comprehensive counting implementation
- MBR feature removal (before/after commits)
- Pair-based CV assignment (commit 5586c678)
- Reversion to protein-group (commit 1a1754a4)

## Conclusion

The investigation revealed systematic data leakage in PSM scoring when MBR is enabled, stemming from the interaction between CV fold assignment and partner precursor addition. While some improvements have been implemented (feature removal, enhanced diagnostics), the fundamental architectural issue of cross-fold information flow remains when target-decoy pairs are split across CV folds.

The current protein-group based assignment maintains biological coherence but doesn't prevent MBR-induced leakage. Future work should focus on developing MBR-aware cross-validation strategies that maintain both biological relevance and statistical independence.