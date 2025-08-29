# PercolatorSortOf.jl Changes Summary

This document summarizes the key changes made to `percolatorSortOf.jl` and explains why the results differ from the develop branch version.

## Overview of Changes

The primary change implements proper separation between training and test predictions in cross-validation to prevent data leakage. This follows the recommendations from ChatGPT discussion to ensure only out-of-fold predictions are used for final evaluation.

## Key Differences

### 1. Probability Storage Variables (Lines 38-40)

**Old Implementation (develop branch):**
```julia
prob_estimates = zeros(Float32, nrow(psms))
MBR_estimates  = zeros(Float32, nrow(psms))
```

**New Implementation (current branch):**
```julia
prob_test   = zeros(Float32, nrow(psms))  # final CV predictions
prob_train  = zeros(Float32, nrow(psms))  # temporary, used during training
MBR_estimates = zeros(Float32, nrow(psms)) # optional MBR layer
```

### 2. Prediction Assignment (Lines 88-102)

**Old Implementation:**
Used a single `predict_fold!` function that updated both training and test predictions in the `:prob` column directly:
```julia
predict_fold!(bst, psms_train, test_fold_psms, train_feats)
```

**New Implementation:**
Explicitly separates training and test predictions:
```julia
# **temporary predictions for training only**
prob_train[train_idx] = predict(bst, psms_train)
psms_train[!,:prob] = prob_train[train_idx]
get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)

# **predict held-out fold**
prob_test[test_idx] = predict(bst, psms_test)
psms_test[!,:prob] = prob_test[test_idx]
```

### 3. Final Probability Assignment (Lines 120-135)

**Old Implementation:**
Used `prob_estimates` for non-MBR cases:
```julia
if match_between_runs
    psms[!, :prob] = MBR_estimates
else
    psms[!, :prob] = prob_estimates
end
```

**New Implementation:**
Uses `prob_test` for non-MBR cases:
```julia
if match_between_runs
    # Determine which precursors failed the q-value cutoff prior to MBR
    qvals_prev = Vector{Float32}(undef, length(prob_test))
    get_qvalues!(prob_test, psms.target, qvals_prev)
    pass_mask = (qvals_prev .<= max_q_value_xgboost_rescore)
    prob_thresh = any(pass_mask) ? minimum(prob_test[pass_mask]) : typemax(Float32)
    # Label as transfer candidates only those failing the q-value cutoff but
    # whose best matched pair surpassed the passing probability threshold.
    psms[!, :MBR_transfer_candidate] .= (prob_test .< prob_thresh) .&
                                        (psms.MBR_max_pair_prob .>= prob_thresh)

    # Use the final MBR probabilities for all precursors
    psms[!, :prob] = MBR_estimates
else
    psms[!, :prob] = prob_test
end
```

## Impact on Results

### Why Results Are Different

1. **Data Leakage Prevention**: The old implementation mixed training and test predictions in the final `:prob` column. While `prob_estimates` was intended to store only out-of-fold predictions, the intermediate steps could contaminate this with in-sample predictions.

2. **Cleaner Separation**: The new implementation maintains strict separation:
   - `prob_train` stores temporary training predictions used only for q-value calculations and MBR feature computation within each fold
   - `prob_test` stores only held-out predictions that were never used for training

3. **More Accurate Evaluation**: By ensuring only true out-of-fold predictions are used for final evaluation, the new implementation provides more accurate estimates of model performance and eliminates potential overfitting artifacts.

### Expected Changes in Results

- **Slightly Different Probability Values**: Out-of-fold predictions are now guaranteed to be unbiased
- **More Conservative Performance Metrics**: True cross-validation without data leakage typically yields slightly more conservative (realistic) performance estimates
- **Improved Generalization**: Models trained with proper CV separation should generalize better to new data

## Technical Details

### Training vs Test Prediction Flow

**Each CV Fold Process:**
1. Train model on K-1 folds
2. Generate predictions for training folds → store in `prob_train[train_idx]`
3. Generate predictions for held-out fold → store in `prob_test[test_idx]`
4. Use training predictions temporarily for q-value calculations and MBR features
5. Training predictions are discarded after each fold; only test predictions persist

### MBR (Match Between Runs) Handling

The MBR logic uses the properly separated predictions:
- Pre-MBR q-values computed from `prob_test` (line 123)
- Transfer candidate determination based on `prob_test` vs `MBR_max_pair_prob`
- Final probabilities come from `MBR_estimates` when MBR is enabled

## Recommendations

1. **Use New Implementation**: The new version provides more rigorous cross-validation and should be preferred for all analyses.

2. **Expect Different Results**: When comparing to previous analyses, expect slightly different probability values and performance metrics due to the elimination of data leakage.

3. **Validation**: The new implementation better represents how the model will perform on truly unseen data, making it more suitable for production use.

4. **Backwards Compatibility**: While results differ, the interface and functionality remain the same. Existing code using this function should work without modification.

## References

This implementation follows best practices for cross-validation as discussed in the ChatGPT conversation, specifically addressing the concern about mixing training and test predictions that could lead to overly optimistic performance estimates.