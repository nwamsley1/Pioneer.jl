# Plan: Fix Probit Regression CV Fold Handling

## Problem Statement

The probit regression model in ScoringSearch incorrectly handles cross-validation folds, causing inconsistency with the XGBoost implementation and breaking protein-aware cross-validation.

### Current Issues:
1. **Arbitrary CV fold assignment**: Probit regression assigns CV folds based on row index (`psms[i, :cv_fold] = UInt8(((i - 1) % n_folds) + 1)`)
2. **Configurable n_folds**: Allows setting n_folds (currently 3) instead of using the library's 2-fold system
3. **Ignores protein grouping**: Overwrites protein-based CV assignments from the library
4. **Inconsistent with XGBoost**: XGBoost correctly uses `unique(psms[!, :cv_fold])` from library assignments

### How XGBoost Works (Correct Approach):
```julia
# In percolatorSortOf.jl:
unique_cv_folds = unique(psms[!, :cv_fold])  # Gets [0, 1] from library
cv_fold_col = psms[!, :cv_fold]              # Uses existing assignments
fold_indices = Dict(fold => findall(==(fold), cv_fold_col) for fold in unique_cv_folds)
train_indices = Dict(fold => findall(!=(fold), cv_fold_col) for fold in unique_cv_folds)
```

### How Library Assigns CV Folds:
- Library creates protein-based CV folds (values: 0 and 1)
- All peptides from the same protein get the same CV fold
- Assigned in `LibraryIon.jl`: `cv_folds = UInt8[0, 1]`
- PSMs receive these in SecondPassSearch: `psms[!,:cv_fold] = cv_fold`

## Implementation Plan

### 1. Fix `probit_regression_scoring_cv!` in `score_psms.jl`

#### Remove:
- The `n_folds` parameter entirely
- The arbitrary CV fold assignment loop
- The `:cv_fold` column creation

#### Add:
- Validation that `:cv_fold` column exists
- Detection of unique CV folds from data
- Use existing CV fold assignments

#### Modified function signature:
```julia
function probit_regression_scoring_cv!(
    psms::DataFrame,
    file_paths::Vector{String},  # Keep for compatibility, unused
    features::Vector{Symbol},
    match_between_runs::Bool
)
    # No n_folds parameter!
```

#### Implementation changes:
```julia
# Validate CV fold column exists
if !(:cv_fold in propertynames(psms))
    error("PSMs must have :cv_fold column from library assignments")
end

# Detect unique CV folds from data (should be [0, 1])
unique_cv_folds = sort(unique(psms[!, :cv_fold]))
n_folds = length(unique_cv_folds)

# Validate we have expected 2 folds
if unique_cv_folds != UInt8[0, 1]
    @user_warn "Unexpected CV folds: $unique_cv_folds (expected [0, 1])"
end

# Use existing CV folds for train/test splits
for fold_idx in unique_cv_folds
    test_mask = psms[!, :cv_fold] .== fold_idx
    train_mask = .!test_mask
    # ... rest of training logic
end
```

### 2. Update `train_probit_model` in `model_comparison.jl`

#### Changes needed:
- Remove `n_folds` from hyperparameters
- Ensure PSMs have `:cv_fold` before calling probit function
- Add CV fold from precursors if missing

```julia
function train_probit_model(config::ModelConfig, 
                          psms_train::DataFrame,
                          match_between_runs::Bool)
    
    # Create copy to avoid modifying original
    psms_copy = copy(psms_train)
    
    # Ensure CV fold column exists
    if !(:cv_fold in propertynames(psms_copy))
        error("PSMs must have :cv_fold column. This should come from SecondPassSearch.")
    end
    
    # Filter features
    available_features = filter(f -> f in propertynames(psms_copy), config.features)
    if match_between_runs
        mbr_features = [f for f in propertynames(psms_copy) if startswith(String(f), "MBR_")]
        append!(available_features, mbr_features)
    end
    
    # Call probit function WITHOUT n_folds parameter
    probit_regression_scoring_cv!(
        psms_copy,
        temp_file_paths,  # Dummy paths for compatibility
        available_features,
        match_between_runs
        # NO n_folds parameter!
    )
    
    return (model = nothing, trained_features = available_features)
end
```

### 3. Update ProbitRegression configuration in `model_comparison.jl`

Remove n_folds from hyperparameters:
```julia
ModelConfig(
    "ProbitRegression",
    :probit,
    REDUCED_FEATURE_SET,
    Dict(
        # Remove: :n_folds => 3,
        :max_iter => 30
    )
)
```

### 4. Update any other probit calls

Search for and update any other places that call probit_regression_scoring_cv! with n_folds parameter.

## Validation Steps

After implementation:
1. Verify probit uses same CV folds as XGBoost (0 and 1)
2. Check that peptides from same protein stay in same fold
3. Confirm model comparison works correctly with all three models
4. Test with datasets that have different protein distributions

## Benefits

1. **Consistency**: All models use same CV strategy
2. **Protein-aware CV**: Prevents data leakage between related peptides
3. **Fair comparison**: Models trained/tested on identical splits
4. **Simplified API**: No confusing n_folds parameter
5. **Correctness**: Respects library's protein-based fold assignments

## Testing Checklist

- [ ] Probit regression runs without n_folds parameter
- [ ] CV folds are [0, 1] from library
- [ ] Model comparison completes successfully
- [ ] Probit and XGBoost use identical train/test splits
- [ ] No arbitrary CV fold reassignment occurs
- [ ] Error thrown if cv_fold column missing
- [ ] Warning shown if unexpected CV fold values found