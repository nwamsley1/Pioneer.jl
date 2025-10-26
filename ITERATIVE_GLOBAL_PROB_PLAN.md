# Iterative ML Training for Global Probability Prediction

## Background

### Current Implementation (Single-Pass Training)
The global probability ML system (`global_prob_ml.jl`) currently uses a single-pass training approach:

1. **Feature Calculation**: Extract 7 features from precursor probability scores across runs
   - logodds_score, n_observations, median_score, std_score, max_score, skewness, kurtosis
2. **Model Training**: Train two LightGBM models (one per CV fold) on ALL targets vs ALL decoys
3. **Prediction**: Generate out-of-fold predictions for global probability scores
4. **Model Selection**: Compare ML model vs logodds baseline using AUC

**Limitation**: Training on all data includes:
- Low-quality targets (ambiguous peptides with poor evidence)
- High-scoring decoys (false positives that look like real hits)
- This dilutes the signal and makes it harder for the model to learn discriminative patterns

### Percolator-Style Iterative Training (Reference Implementation)
The PSM-level scoring in `percolatorSortOf.jl` uses iterative refinement:

**Iteration 1**: Train on ALL targets vs ALL decoys
```julia
# From percolatorSortOf.jl:979-981
if itr == 1
    return copy(psms_train)  # Use all data
```

**Iteration 2+**: Refine training set
```julia
# From percolatorSortOf.jl:986-1004
# 1. Convert worst-scoring targets (high PEP) to negatives
PEPs = calculate_PEP(sorted_scores, sorted_targets)
idx_cutoff = findfirst(x -> x >= min_PEP_neg_threshold_itr, PEPs)
if !isnothing(idx_cutoff)
    worst_idxs = order[idx_cutoff:end]
    psms_train_itr.target[worst_idxs] .= false  # Relabel as negatives
end

# 2. Keep only:
#    - All decoys
#    - High-confidence targets (q_value <= threshold)
subset(psms_train_itr,
    [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value)))
```

**Benefits**:
- Focuses on clear positives vs clear negatives
- Reduces label noise from ambiguous targets
- Improves model discrimination by learning from high-confidence examples

## Proposed Implementation: Iterative Global Probability Training

### Algorithm Overview

**Step 1: Initial Training (Iteration 1)**
1. Calculate precursor-level features from prec_prob/MBR_boosted_prec_prob
2. Train LightGBM models on ALL targets vs ALL decoys
3. Generate out-of-fold predictions → global_prob_iter1
4. Calculate global q-values using FDR estimation

**Step 2: Refined Training (Iteration 2)**
1. **Identify high-confidence targets**:
   - Calculate PEP (Posterior Error Probability) from global_prob_iter1
   - Mark targets with PEP ≥ threshold as uncertain → relabel as negatives
   - Keep targets with global_qval ≤ threshold as high-confidence positives

2. **Create refined training set**:
   - Positives: Targets with low q-value AND low PEP
   - Negatives: All original decoys + relabeled uncertain targets

3. **Retrain models** on refined set
4. Generate final predictions → global_prob (or MBR_boosted_global_prob)

### Key Differences from PSM-Level Percolator

| Aspect | PSM-Level (percolatorSortOf) | Precursor-Level (global_prob_ml) |
|--------|------------------------------|-----------------------------------|
| **Training unit** | Individual PSM isotope traces | Unique precursors (aggregated across runs) |
| **Features** | ~50 spectral/RT/MS1 features | 7 cross-run distribution features |
| **Iterations** | Multiple (3-5 depending on MBR) | 2 (initial + refined) |
| **Data volume** | 10K-1M PSMs per file | 1K-100K unique precursors total |
| **CV strategy** | Fold assignment from library | Same (inherited from precursors) |
| **PEP calculation** | From trace_prob scores | From global_prob_iter1 scores |

## Implementation Details

### 1. New Function: `get_training_data_for_global_prob_iteration!`

Location: `global_prob_ml.jl` (add after `train_global_prob_models`)

```julia
"""
    get_training_data_for_global_prob_iteration!(
        precursor_features::DataFrame,
        itr::Int,
        global_prob_iter1::Union{Nothing, Vector{Float32}},
        max_global_qval_threshold::Float32,
        min_PEP_neg_threshold::Float32
    ) -> DataFrame

Filter and refine training data for iterative global probability training.

# Arguments
- `precursor_features`: DataFrame with precursor-level features and labels
- `itr`: Current iteration number (1 or 2)
- `global_prob_iter1`: Predictions from iteration 1 (nothing for iteration 1)
- `max_global_qval_threshold`: Max q-value for keeping targets (e.g., 0.01)
- `min_PEP_neg_threshold`: Min PEP for relabeling targets as negatives (e.g., 0.5)

# Returns
DataFrame with refined training set for current iteration

# Iteration 1
Returns copy of all data (no filtering)

# Iteration 2
1. Calculate PEP from global_prob_iter1
2. Relabel targets with PEP ≥ threshold as negatives (uncertain targets)
3. Filter to keep:
   - All original decoys
   - All relabeled uncertain targets (now negatives)
   - High-confidence targets (global_qval ≤ threshold)
"""
function get_training_data_for_global_prob_iteration!(
    precursor_features::DataFrame,
    itr::Int,
    global_prob_iter1::Union{Nothing, Vector{Float32}},
    max_global_qval_threshold::Float32,
    min_PEP_neg_threshold::Float32
)
    if itr == 1
        # Iteration 1: Train on all data
        return copy(precursor_features)
    else
        # Iteration 2: Refine training set
        @assert !isnothing(global_prob_iter1) "global_prob_iter1 required for iteration 2"

        # Create working copy to avoid modifying original
        features_refined = copy(precursor_features)

        # Add iteration 1 predictions to calculate PEP
        features_refined.global_prob_iter1 = global_prob_iter1

        # Calculate PEP from iteration 1 scores
        order = sortperm(features_refined.global_prob_iter1, rev=true)
        sorted_scores = features_refined.global_prob_iter1[order]
        sorted_targets = features_refined.target[order]

        PEPs = Vector{Float32}(undef, length(order))
        get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

        # Calculate global q-values from iteration 1 scores
        global_qvals = Vector{Float32}(undef, length(order))
        get_q_values!(sorted_scores, sorted_targets, global_qvals; doSort=false)

        # Map PEPs and q-values back to original order
        features_refined.PEP_iter1 = similar(PEPs)
        features_refined.global_qval_iter1 = similar(global_qvals)
        features_refined.PEP_iter1[order] = PEPs
        features_refined.global_qval_iter1[order] = global_qvals

        # Relabel uncertain targets (high PEP) as negatives
        uncertain_target_mask = features_refined.target .&&
                                (features_refined.PEP_iter1 .>= min_PEP_neg_threshold)
        features_refined.target[uncertain_target_mask] .= false

        @user_info "  Iteration 2: Relabeled $(sum(uncertain_target_mask)) uncertain targets as negatives"

        # Filter to keep:
        # - All negatives (original decoys + relabeled uncertain targets)
        # - High-confidence targets (q-value <= threshold)
        features_refined = subset(
            features_refined,
            [:target, :global_qval_iter1] => ByRow((t, q) ->
                (!t) || (t && q <= max_global_qval_threshold)
            )
        )

        @user_info "  Iteration 2: Training set size: $(nrow(features_refined)) precursors"
        @user_info "    - Targets: $(sum(features_refined.target))"
        @user_info "    - Negatives: $(sum(.!features_refined.target))"

        # Drop temporary columns
        select!(features_refined, Not([:global_prob_iter1, :PEP_iter1, :global_qval_iter1]))

        return features_refined
    end
end
```

### 2. Modified Function: `train_global_prob_models`

Update signature and logic to support iterative training:

```julia
function train_global_prob_models(
    precursor_features::DataFrame,
    params,  # ScoringSearchParameters
    n_iterations::Int = 2,  # NEW: Number of iterations (default 2)
    max_global_qval_threshold::Float32 = 0.01f0,  # NEW: Q-value threshold for iteration 2
    min_PEP_neg_threshold::Float32 = 0.5f0  # NEW: PEP threshold for relabeling
)::Tuple{Union{Nothing, Any}, Union{Nothing, Any}, Vector{Symbol}}

    # Feature columns (unchanged)
    feature_cols = Symbol[
        :logodds_score, :n_observations, :median_score, :std_score,
        :max_score, :skewness, :kurtosis
    ]

    @user_info "Training global probability models with $(n_iterations) iterations"
    @user_info "Features: $(join(string.(feature_cols), ", "))"

    # Check CV folds
    unique_folds = unique(precursor_features.cv_fold)
    if length(unique_folds) != 2
        @user_warn "Expected 2 CV folds, got $(length(unique_folds)). Using logodds only."
        return (nothing, nothing, feature_cols)
    end

    # Storage for iteration 1 predictions (needed for iteration 2 filtering)
    global_prob_iter1 = nothing

    # Iterative training loop
    for itr in 1:n_iterations
        @user_info "Global Probability Training - Iteration $itr"

        # Get training data for this iteration
        features_itr = get_training_data_for_global_prob_iteration!(
            precursor_features,
            itr,
            global_prob_iter1,
            max_global_qval_threshold,
            min_PEP_neg_threshold
        )

        # Train models on current iteration's data
        models = Vector{Any}(undef, 2)

        for (fold_idx, fold) in enumerate(unique_folds)
            test_mask = features_itr.cv_fold .== fold
            train_mask = .!test_mask

            train_data = features_itr[train_mask, feature_cols]
            train_labels = features_itr.target[train_mask]

            # Train LightGBM classifier
            classifier = build_lightgbm_classifier(
                num_iterations = 200,
                max_depth = 10,
                learning_rate = 0.1,
                feature_fraction = 1.0,
                bagging_fraction = 0.5,
                bagging_freq = 1,
                min_data_in_leaf = 100,
                min_gain_to_split = 1.0
            )

            models[fold_idx] = fit_lightgbm_model(
                classifier, train_data, train_labels;
                positive_label=true
            )

            @user_info "  Fold $fold (iter $itr): trained on $(sum(train_mask)) precursors, test $(sum(test_mask))"
        end

        # Generate predictions for this iteration (on FULL dataset, not filtered)
        # This ensures we have predictions for all precursors for iteration 2 filtering
        if itr < n_iterations
            global_prob_iter1 = predict_global_probs_cv(
                precursor_features,  # Use FULL dataset for predictions
                models[1],
                models[2],
                feature_cols
            )
        else
            # Final iteration - return the trained models
            # Print feature importances
            @user_info "Feature Importances (Final Iteration):"
            for (fold_idx, model) in enumerate(models)
                importances = lightgbm_feature_importances(model)
                if importances === nothing
                    @user_warn "  Fold $fold_idx: No feature importances available"
                else
                    feature_pairs = collect(zip(feature_cols, importances))
                    sort!(feature_pairs, by=x->x[2], rev=true)
                    @user_info "  Fold $fold_idx:"
                    for (feat, score) in feature_pairs
                        @user_info "    $feat: $(round(score, digits=3))"
                    end
                end
            end

            return models[1], models[2], feature_cols
        end
    end
end
```

### 3. Parameter Additions

Add to `defaultSearchParams.json`:

```json
"machine_learning": {
    ...existing params...
    "global_prob_n_iterations": 2,
    "global_prob_qvalue_threshold": 0.01,
    "global_prob_min_PEP_threshold": 0.5
}
```

Add to `ScoringSearchParameters` struct:

```julia
struct ScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    ...existing fields...

    # Global probability iterative training parameters
    global_prob_n_iterations::Int
    global_prob_qvalue_threshold::Float32
    global_prob_min_PEP_threshold::Float32

    function ScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        ...

        new{typeof(isotope_trace_type)}(
            ...existing fields...,

            # Global prob iteration parameters
            Int(get(ml_params, :global_prob_n_iterations, 2)),
            Float32(get(ml_params, :global_prob_qvalue_threshold, 0.01)),
            Float32(get(ml_params, :global_prob_min_PEP_threshold, 0.5))
        )
    end
end
```

### 4. Integration in ScoringSearch.jl Step 5

Update both MBR and non-MBR branches:

```julia
# Current (line ~507):
model_fold_0, model_fold_1, feature_cols = train_global_prob_models(
    precursor_features, params
)

# Updated:
model_fold_0, model_fold_1, feature_cols = train_global_prob_models(
    precursor_features,
    params,
    params.global_prob_n_iterations,
    params.global_prob_qvalue_threshold,
    params.global_prob_min_PEP_threshold
)
```

## Testing Strategy

### 1. Unit Tests
Create `test/UnitTests/GlobalProbML/test_iterative_training.jl`:

```julia
@testset "Iterative Global Probability Training" begin
    # Test iteration 1: all data
    features_iter1 = get_training_data_for_global_prob_iteration!(
        test_features, 1, nothing, 0.01f0, 0.5f0
    )
    @test nrow(features_iter1) == nrow(test_features)

    # Test iteration 2: refined data
    mock_predictions = randn(Float32, nrow(test_features))
    features_iter2 = get_training_data_for_global_prob_iteration!(
        test_features, 2, mock_predictions, 0.01f0, 0.5f0
    )
    @test nrow(features_iter2) < nrow(test_features)  # Should filter some
    @test sum(features_iter2.target) < sum(test_features.target)  # Relabeled some targets
end
```

### 2. Integration Test
Run on E. coli test dataset:

```julia
# Compare single-pass vs iterative
params_single_pass = copy(test_params)
params_single_pass.optimization.machine_learning.global_prob_n_iterations = 1

params_iterative = copy(test_params)
params_iterative.optimization.machine_learning.global_prob_n_iterations = 2

# Run both and compare:
# - Number of identified proteins at 1% FDR
# - Global probability AUC
# - Training set sizes in iteration 2
```

### 3. Validation Metrics
Track these metrics to validate improvement:

1. **AUC Improvement**:
   - Iteration 1 AUC vs Iteration 2 AUC
   - Should see modest improvement (0.01-0.05 typical)

2. **Training Set Refinement**:
   - Number of targets relabeled as negatives
   - Percentage of targets kept (should be ~50-80%)

3. **Downstream Impact**:
   - Number of protein groups at 1% FDR
   - Median protein q-values

4. **Computational Cost**:
   - Training time for iteration 2 (should be faster due to smaller dataset)
   - Total time overhead (target: <20% increase)

## Expected Results

### Improvements
1. **Better Discrimination**: Model learns from clearer positives/negatives
2. **Reduced Overfitting**: Fewer noisy labels in training set
3. **Higher Confidence**: Targets passing iteration 2 have stronger evidence

### Trade-offs
1. **Slightly Longer Training**: 2x model training (mitigated by smaller iteration 2 dataset)
2. **Parameter Sensitivity**: Need to tune q-value and PEP thresholds
3. **Diminishing Returns**: May not help datasets with very clear separation

## Implementation Phases

### Phase 1: Core Functions (Immediate)
- [ ] Implement `get_training_data_for_global_prob_iteration!`
- [ ] Modify `train_global_prob_models` for iterative training
- [ ] Add parameters to config files
- [ ] Update `ScoringSearchParameters` struct

### Phase 2: Integration (After Phase 1 works)
- [ ] Update ScoringSearch.jl Step 5 call sites
- [ ] Add parameter validation in paramsChecks.jl
- [ ] Test on small dataset (E. coli)

### Phase 3: Validation & Optimization (After Phase 2 works)
- [ ] Run on larger datasets (HeLa, yeast)
- [ ] Compare metrics: AUC, proteins identified, q-values
- [ ] Tune default thresholds based on results
- [ ] Document findings in this plan

### Phase 4: Advanced Features (Optional)
- [ ] Adaptive thresholds based on dataset size
- [ ] More than 2 iterations (if beneficial)
- [ ] Per-entrapment-group iteration (for entrapment analysis)

## References

- **Percolator Algorithm**: Käll et al. (2007) "Semi-supervised learning for peptide identification from shotgun proteomics datasets"
- **Pioneer PSM Scoring**: `src/utils/ML/percolatorSortOf.jl:969-1006` (get_training_data_for_iteration!)
- **Global Probability ML**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/global_prob_ml.jl`
