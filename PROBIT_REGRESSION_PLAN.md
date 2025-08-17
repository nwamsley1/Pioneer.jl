# Plan: Probit Regression Alternative for Small PSM Sets

## Overview
Implement probit regression as an alternative to XGBoost/EvoTrees for cases with under 100,000 PSMs. This will provide a simpler, more interpretable model for smaller datasets while maintaining the XGBoost approach for larger datasets.

## Motivation
- **Simplicity**: Probit regression is simpler and more interpretable than gradient boosting
- **Speed**: Faster training for small datasets without cross-validation overhead
- **Stability**: More stable with limited training data
- **Comparison**: Easy A/B testing by commenting/uncommenting different approaches

## Existing Infrastructure

Pioneer already has a robust probit regression implementation used in FirstPassSearch and ParameterTuningSearch:

### Core Functions (in `src/utils/ML/probitRegression.jl`)
- `ProbitRegression(β, X::DataFrame, y::Vector{Bool}, data_chunks, ...)` - Main training function
- `ModelPredict!(scores, psms::DataFrame, β, data_chunks)` - Predict Z-scores
- `ModelPredictProbs!(scores, psms::DataFrame, β, data_chunks)` - Predict probabilities (uses erf for CDF)

### Usage Examples in Codebase
1. **FirstPassSearch** (`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`):
   - Uses iterative training with best PSMs selection
   - Features: hyperscore, RT error, mass error, etc.
   
2. **ParameterTuningSearch** (`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`):
   - Single-pass training for scoring
   - Similar feature set

## Current State Analysis

### File: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

Current flow in `score_precursor_isotope_traces`:
```julia
if n_psms < 100_000
    # In-memory processing
    all_psms = load_all_psms()
    update_isotope_features!()
    sort_of_percolator_in_memory!(all_psms, ...)  # XGBoost/EvoTrees
else
    # Out-of-memory processing
    sort_of_percolator_out_of_memory!(...)
end
```

### Current `sort_of_percolator_in_memory!` Implementation
- Trains EvoTrees/XGBoost models with cross-validation
- Uses multiple features: spectral scores, RT errors, mass errors, etc.
- Implements iterative refinement with `max_prob` updates
- Handles both targets and decoys

## Proposed Implementation

### 1. Create New Function: `probit_regression_scoring_cv!`

Location: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

```julia
function probit_regression_scoring_cv!(
    all_psms::DataFrame,
    psm_scores::Vector{Symbol},
    ms_file_idx_col::Symbol,
    n_folds::Int64,
    spec_lib::SpectralLibrary
)
    # Use existing ProbitRegression from utils/ML/probitRegression.jl
end
```

### 2. Function Structure (Using Existing Infrastructure)

```julia
function probit_regression_scoring_cv!(
    all_psms::DataFrame,
    psm_scores::Vector{Symbol},
    ms_file_idx_col::Symbol,
    n_folds::Int64,
    spec_lib::SpectralLibrary
)
    # Step 1: Add CV fold column (same as XGBoost approach)
    all_psms[!, :cv_fold] = zeros(UInt8, size(all_psms, 1))
    for fold_idx in 1:n_folds
        fold_mask = (all_psms[!, ms_file_idx_col] .% n_folds) .== (fold_idx - 1)
        all_psms[fold_mask, :cv_fold] .= UInt8(fold_idx)
    end
    
    # Step 2: Initialize score columns
    all_psms[!, :prob] = zeros(Float32, size(all_psms, 1))
    
    # Step 3: Prepare features (same as XGBoost but simpler)
    # Use core scoring features
    feature_cols = [
        :hyperscore,
        :rt_irt_err,
        :ppm_err,
        :masked_ppm_err,
        :isotope_err,
        :missed_cleavages,
        :nce_err,
        :y_count,
        :b_count,
        :intercept  # Add intercept column (all ones) as in FirstPassSearch
    ]
    
    # Ensure intercept column exists
    if !(:intercept in propertynames(all_psms))
        all_psms[!, :intercept] = ones(Float32, size(all_psms, 1))
    end
    
    # Filter to available columns
    available_features = filter(col -> col in propertynames(all_psms), feature_cols)
    
    # Step 4: Train probit model per CV fold (similar to FirstPassSearch)
    tasks_per_thread = 10
    chunk_size = max(1, size(all_psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = Iterators.partition(1:size(all_psms, 1), chunk_size)
    
    for fold_idx in 1:n_folds
        # Get train/test masks
        test_mask = all_psms[!, :cv_fold] .== fold_idx
        train_mask = .!test_mask
        
        # Separate targets and decoys for training
        train_targets = train_mask .& all_psms[!, :target]
        train_decoys = train_mask .& .!all_psms[!, :target]
        
        # Skip if insufficient data
        if sum(train_targets) < 100 || sum(train_decoys) < 100
            @warn "Insufficient training data for fold $fold_idx"
            all_psms[test_mask, :prob] .= 0.5f0
            continue
        end
        
        # Train on this fold's training data
        train_data = all_psms[train_mask, :]
        train_targets_bool = train_data[!, :target]
        
        # Prepare chunks for training data
        train_chunk_size = max(1, size(train_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
        train_chunks = Iterators.partition(1:size(train_data, 1), train_chunk_size)
        
        # Initialize coefficients
        β = zeros(Float64, length(available_features))
        
        # Fit probit model using existing ProbitRegression
        β = Pioneer.ProbitRegression(
            β, 
            train_data[!, available_features], 
            train_targets_bool, 
            train_chunks, 
            max_iter = 30
        )
        
        # Predict probabilities on test set using ModelPredictProbs!
        test_data = all_psms[test_mask, :]
        test_probs = zeros(Float32, size(test_data, 1))
        
        test_chunk_size = max(1, size(test_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
        test_chunks = Iterators.partition(1:size(test_data, 1), test_chunk_size)
        
        Pioneer.ModelPredictProbs!(
            test_probs,
            test_data[!, available_features],
            β,
            test_chunks
        )
        
        all_psms[test_mask, :prob] = test_probs
    end
    
    # Step 5: Calculate q-values using existing function
    all_psms[!, :q_value] = Pioneer.getQvalues(
        all_psms[!, :prob],
        all_psms[!, :target],
        1.0f0  # prior
    )
    
    # Step 6: Select best PSM per precursor (same as XGBoost)
    all_psms[!, :best_psm] = zeros(Bool, size(all_psms, 1))
    gdf = DataFrames.groupby(all_psms, :precursor_idx)
    for (key, group) in pairs(gdf)
        best_idx = argmax(group[!, :prob])
        group[best_idx, :best_psm] = true
    end
    
    return nothing
end
```

### 3. Integration in `score_precursor_isotope_traces`

```julia
function score_precursor_isotope_traces(
    all_poss_isotopes::DataFrame,
    second_pass_psm_paths::Vector{String},
    spec_lib::SpectralLibrary,
    psm_scores::Vector{Symbol},
    max_best_rank::UInt8,
    q_value_threshold::Float32,
    prior::Float32,
    n_train_rounds::Int64,
    max_iter_per_round::Int64,
    n_folds::Int64,
    machine_learning_params::NamedTuple
)
    # ... existing code ...
    
    if n_psms < 100_000
        # In-memory processing
        all_psms = load_and_prepare_psms(...)
        
        # OPTION 1: Probit regression (NEW - SIMPLE, NO ITERATIVE REFINEMENT)
        # probit_regression_scoring_cv!(
        #     all_psms,
        #     psm_scores,
        #     ms_file_idx_col,
        #     n_folds,
        #     spec_lib
        # )
        
        # OPTION 2: XGBoost/EvoTrees (EXISTING - WITH ITERATIVE REFINEMENT)
        sort_of_percolator_in_memory!(
            all_psms,
            psm_scores,
            n_train_rounds,
            max_iter_per_round,
            prior,
            n_folds,
            ms_file_idx_col,
            machine_learning_params
        )
        
        # Rest of processing remains the same
        write_scored_psms(...)
    else
        # Out-of-memory processing (unchanged)
        sort_of_percolator_out_of_memory!(...)
    end
end
```

## Implementation Steps

### Phase 1: Create Probit Function
1. [ ] Add `probit_regression_scoring_cv!` function to `score_psms.jl`
2. [ ] Use existing `ProbitRegression` from `utils/ML/probitRegression.jl`
3. [ ] Use existing `ModelPredictProbs!` for probability predictions
4. [ ] Ensure intercept column is added (as in FirstPassSearch)

### Phase 2: Integration
1. [ ] Add function call in `score_precursor_isotope_traces`
2. [ ] Comment out by default (keep XGBoost active)
3. [ ] Ensure output columns match XGBoost approach (`:prob`, `:q_value`, `:best_psm`)
4. [ ] Test with small dataset

### Phase 3: Testing Strategy
1. [ ] Run with probit on ecoli_test dataset
2. [ ] Compare PSM counts with XGBoost
3. [ ] Compare q-value distributions
4. [ ] Compare protein group results
5. [ ] Benchmark speed differences

## Key Simplifications vs XGBoost

1. **No iterative refinement**: Single pass training per CV fold (like FirstPassSearch but with CV)
2. **No max_prob updates**: Direct probability calculation using erf/normal CDF
3. **No tree-based features**: Linear probit model only
4. **Simpler feature set**: Core spectral features without complex interactions
5. **No gradient boosting**: Single linear model per fold
6. **Existing infrastructure**: Reuses proven ProbitRegression implementation

## Expected Benefits

1. **Speed**: 5-10x faster for small datasets
2. **Interpretability**: Linear coefficients show feature importance
3. **Stability**: No hyperparameter tuning required
4. **Memory**: Lower memory footprint
5. **Simplicity**: Easier to debug and understand

## Risk Mitigation

1. **Backward compatibility**: Keep XGBoost as default
2. **Easy switching**: Simple comment/uncomment to switch
3. **Same output format**: Identical column structure
4. **Fallback logic**: Use XGBoost if probit fails

## Success Criteria

1. Probit regression runs without errors
2. Q-value calculation produces reasonable FDR
3. Protein identification within 10% of XGBoost
4. Speed improvement >3x for <20k PSMs
5. Memory usage reduced by >50%

## Code Location Summary

- Main implementation: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`
- Add new function: `probit_regression_scoring!`
- Modify function: `score_precursor_isotope_traces`
- Test with: `SearchDIA("./data/ecoli_test/ecoli_test_params.json")`

## Notes

- GLM.jl is already a dependency (used in FirstPassSearch)
- Can reuse q-value calculation from existing code
- Output must include: `:prob`, `:q_value`, `:best_psm` columns
- Consider adding a parameter flag later for automatic selection