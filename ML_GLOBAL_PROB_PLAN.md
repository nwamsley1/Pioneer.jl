# Machine Learning Global Probability Prediction - Implementation Plan

## Overview
Replace the simple `logodds(prec_prob, sqrt_n_runs)` calculation in Step 5 of ScoringSearch with an adaptive ML-based approach. Train a LightGBM classifier to predict target/decoy labels at the precursor level using cross-run features, then use model predictions as `global_prob`. Compare ML vs logodds performance using AUC and automatically select the best approach.

## User Requirements Summary
- **Target variable**: Binary target/decoy classification (like percolator)
- **Top N features**: User-defined parameter (default=3), cannot exceed n_files
- **Model selection**: Compare AUC on first CV fold (ML vs logodds), use best
- **Training**: Single-pass (no iterative refinement), 2 models (one per CV fold)
- **Feature**: Include logodds score as a feature in ML model

---

## Implementation Phases

### Phase 1: Feature Calculation Helper Functions

**File**: Create new file `src/Routines/SearchDIA/SearchMethods/ScoringSearch/global_prob_ml.jl`

#### Feature Calculation Functions
Since we use `groupby` to get all prec_prob values for each precursor, we can directly compute features without needing accumulators.

```julia
"""
    calculate_top_n_features(prec_probs::AbstractVector{Float32}, n_top_scores::Int) -> Vector{Float32}

Extract top N scores from a vector of probabilities, padding with 0.0 if needed.
"""
function calculate_top_n_features(prec_probs::AbstractVector{Float32}, n_top_scores::Int)
    sorted_scores = sort(prec_probs, rev=true)
    top_n = Vector{Float32}(undef, n_top_scores)

    for i in 1:n_top_scores
        top_n[i] = i <= length(sorted_scores) ? sorted_scores[i] : 0.0f0
    end

    return top_n
end

"""
    calculate_precursor_features(prec_probs::AbstractVector{Float32},
                                 n_top_scores::Int,
                                 sqrt_n_runs::Int) -> NamedTuple

Calculate all features for a single precursor from its prec_prob values.
"""
function calculate_precursor_features(
    prec_probs::AbstractVector{Float32},
    n_top_scores::Int,
    sqrt_n_runs::Int
)
    # Basic statistics
    n_obs = length(prec_probs)
    mean_score = mean(prec_probs)
    median_score = median(prec_probs)
    std_score = std(prec_probs)
    min_score = minimum(prec_probs)
    max_score = maximum(prec_probs)

    # Logodds baseline
    logodds_score = logodds(prec_probs, sqrt_n_runs)

    # Top N scores
    top_n = calculate_top_n_features(prec_probs, n_top_scores)

    return (
        logodds_score = logodds_score,
        n_observations = n_obs,
        mean_score = mean_score,
        median_score = median_score,
        std_score = std_score,
        min_score = min_score,
        max_score = max_score,
        top_n_scores = top_n
    )
end
```

**Design choice**: Use DataFrames `groupby` to efficiently group all prec_prob values by precursor_idx, then compute features directly. This is simpler and more memory-efficient than accumulators since we need all values for logodds anyway.

---

### Phase 2: Feature Engineering Pipeline

#### Build Precursor-Level DataFrame
**Design**: Accept a `score_column` parameter to handle both MBR and non-MBR cases uniformly.

```julia
function build_precursor_features_df(
    merged_df::DataFrame,
    n_top_scores::Int,
    sqrt_n_runs::Int,
    score_column::Symbol = :prec_prob  # Can be :prec_prob or :MBR_boosted_prec_prob
)::DataFrame
    # Use groupby to get all score values for each precursor
    # This is similar to how current Step 5 works with transform!
    grouped = groupby(merged_df, :precursor_idx)

    # Pre-allocate result DataFrame
    n_precursors = length(grouped)
    precursor_features = DataFrame(
        precursor_idx = Vector{UInt32}(undef, n_precursors),
        target = Vector{Bool}(undef, n_precursors),
        cv_fold = Vector{UInt8}(undef, n_precursors),
        logodds_score = Vector{Float32}(undef, n_precursors),
        n_observations = Vector{Int}(undef, n_precursors),
        mean_score = Vector{Float32}(undef, n_precursors),
        median_score = Vector{Float32}(undef, n_precursors),
        std_score = Vector{Float32}(undef, n_precursors),
        min_score = Vector{Float32}(undef, n_precursors),
        max_score = Vector{Float32}(undef, n_precursors)
    )

    # Add top N score columns
    for i in 1:n_top_scores
        precursor_features[!, Symbol("top_$(i)_score")] = Vector{Float32}(undef, n_precursors)
    end

    # Process each group (precursor)
    for (idx, group) in enumerate(grouped)
        # Metadata (same for all rows in group)
        precursor_features.precursor_idx[idx] = first(group.precursor_idx)
        precursor_features.target[idx] = first(group.target)
        precursor_features.cv_fold[idx] = first(group.cv_fold)

        # Calculate features from score_column values (prec_prob or MBR_boosted_prec_prob)
        features = calculate_precursor_features(
            getproperty(group, score_column),
            n_top_scores,
            sqrt_n_runs
        )

        # Assign features to DataFrame
        precursor_features.logodds_score[idx] = features.logodds_score
        precursor_features.n_observations[idx] = features.n_observations
        precursor_features.mean_score[idx] = features.mean_score
        precursor_features.median_score[idx] = features.median_score
        precursor_features.std_score[idx] = features.std_score
        precursor_features.min_score[idx] = features.min_score
        precursor_features.max_score[idx] = features.max_score

        # Assign top N scores
        for i in 1:n_top_scores
            col_name = Symbol("top_$(i)_score")
            precursor_features[idx, col_name] = features.top_n_scores[i]
        end
    end

    return precursor_features
end
```

---

### Phase 3: Model Training with CV Folds

#### LightGBM Training Function
```julia
function train_global_prob_models(
    precursor_features::DataFrame,
    params  # ScoringSearchParameters
)::Tuple{Union{Nothing, Any}, Union{Nothing, Any}, Vector{Symbol}}
    # Feature columns (exclude metadata)
    feature_cols = Symbol[]
    for col in names(precursor_features)
        if col ∉ [:precursor_idx, :target, :cv_fold]
            push!(feature_cols, Symbol(col))
        end
    end

    @user_info "Training global probability models with $(length(feature_cols)) features"
    @user_info "Features: $(join(string.(feature_cols), ", "))"

    # Train one model per CV fold
    unique_folds = unique(precursor_features.cv_fold)
    @assert length(unique_folds) == 2 "Expected 2 CV folds, got $(length(unique_folds))"

    models = Vector{Any}(undef, 2)

    for (fold_idx, fold) in enumerate(unique_folds)
        test_mask = precursor_features.cv_fold .== fold
        train_mask = .!test_mask

        train_data = precursor_features[train_mask, feature_cols]
        train_labels = precursor_features.target[train_mask]

        # Train LightGBM classifier (use similar hyperparameters to PSM scoring SimpleLightGBM)
        classifier = build_lightgbm_classifier(
            num_iterations = 100,
            max_depth = 4,
            num_leaves = 15,
            learning_rate = 0.1,
            feature_fraction = 0.8,
            bagging_fraction = 0.8,
            bagging_freq = 1,
            min_data_in_leaf = 20,
            min_gain_to_split = 0.1
        )

        models[fold_idx] = fit_lightgbm_model(classifier, train_data, train_labels; positive_label=true)

        @user_info "  Fold $fold: trained on $(sum(train_mask)) precursors, test set $(sum(test_mask)) precursors"
    end

    return models[1], models[2], feature_cols
end
```

#### Prediction Function
```julia
function predict_global_probs_cv(
    precursor_features::DataFrame,
    model_fold_0::Any,
    model_fold_1::Any,
    feature_cols::Vector{Symbol}
)::Vector{Float32}
    predictions = Vector{Float32}(undef, nrow(precursor_features))

    # Predict out-of-fold
    for row_idx in 1:nrow(precursor_features)
        fold = precursor_features.cv_fold[row_idx]
        model = fold == 0 ? model_fold_1 : model_fold_0

        # Extract features for this row
        features_df = precursor_features[row_idx:row_idx, feature_cols]
        pred = lightgbm_predict(model, features_df)
        predictions[row_idx] = pred[1]
    end

    return predictions
end
```

---

### Phase 4: Model Selection via AUC Comparison

#### AUC Calculation Function
```julia
function calculate_auc(scores::Vector{Float32}, labels::Vector{Bool})::Float64
    # Sort by scores descending
    sorted_idx = sortperm(scores, rev=true)
    sorted_labels = labels[sorted_idx]

    n_pos = sum(sorted_labels)
    n_neg = length(sorted_labels) - n_pos

    if n_pos == 0 || n_neg == 0
        return 0.5  # Undefined AUC
    end

    # Calculate AUC using trapezoidal rule
    tpr = cumsum(sorted_labels) / n_pos
    fpr = cumsum(.!sorted_labels) / n_neg

    auc = 0.0
    for i in 2:length(fpr)
        auc += (fpr[i] - fpr[i-1]) * (tpr[i] + tpr[i-1]) / 2
    end

    return auc
end
```

#### Model Selection Function
```julia
function select_best_global_prob_method(
    precursor_features::DataFrame,
    ml_predictions::Vector{Float32}
)::Symbol
    # Use first CV fold for comparison (fold 0)
    fold_0_mask = precursor_features.cv_fold .== 0

    labels = precursor_features.target[fold_0_mask]
    logodds_scores = precursor_features.logodds_score[fold_0_mask]
    ml_scores = ml_predictions[fold_0_mask]

    auc_logodds = calculate_auc(logodds_scores, labels)
    auc_ml = calculate_auc(ml_scores, labels)

    @user_info "Global Probability Method Comparison (CV Fold 0):"
    @user_info "  Logodds AUC: $(round(auc_logodds, digits=4))"
    @user_info "  ML Model AUC: $(round(auc_ml, digits=4))"

    if auc_ml > auc_logodds
        @user_info "  ✓ Selected: ML Model (AUC improvement: +$(round(auc_ml - auc_logodds, digits=4)))"
        return :ml
    else
        @user_info "  ✓ Selected: Logodds (ML did not improve performance)"
        return :logodds
    end
end
```

---

### Phase 5: Integration into ScoringSearch.jl

#### Modify Step 5 in `performSearch!` function

**Location**: `ScoringSearch.jl` lines 477-507

**Design**: Separate paths for MBR on/off
- When MBR is ON: ML for `MBR_boosted_global_prob`, logodds for `global_prob`
- When MBR is OFF: ML for `global_prob`

**Replace**:
```julia
# Step 5: Calculate global probabilities from best traces
```

**With**:
```julia
# Step 5: Calculate global probabilities from best traces using ML
#@debug_l1 "Step 5: Calculating global probabilities from best traces..."
step5_time = @elapsed begin
    # Merge filtered (best-trace-only) PSMs for global probability calculations
    merged_best_traces_path = joinpath(temp_folder, "merged_best_traces.arrow")
    sort_file_by_keys!(filtered_refs, :precursor_idx)
    stream_sorted_merge(filtered_refs, merged_best_traces_path, :precursor_idx)

    merged_df = DataFrame(Arrow.Table(merged_best_traces_path))
    sqrt_n_runs = floor(Int64, sqrt(length(getFilePaths(getMSData(search_context)))))
    n_top_scores = min(params.n_top_scores_global_prob, length(valid_file_indices))

    if hasproperty(merged_df, :MBR_boosted_prec_prob)
        # MBR is ON: Use ML for MBR_boosted_global_prob, logodds for global_prob
        @user_info "Step 5: MBR enabled - calculating MBR_boosted_global_prob with ML"

        # Build features from MBR-boosted scores
        precursor_features = build_precursor_features_df(
            merged_df, n_top_scores, sqrt_n_runs, :MBR_boosted_prec_prob
        )

        # Train ML models
        model_fold_0, model_fold_1, feature_cols = train_global_prob_models(
            precursor_features, params
        )

        # Get predictions and select best method
        ml_predictions = Vector{Float32}(undef, 0)
        selected_method = :logodds  # Default fallback

        if !isnothing(model_fold_0) && !isnothing(model_fold_1)
            ml_predictions = predict_global_probs_cv(
                precursor_features, model_fold_0, model_fold_1, feature_cols
            )
            selected_method = select_best_global_prob_method(
                precursor_features, ml_predictions
            )
        else
            @user_warn "ML model training failed - falling back to logodds"
        end

        # Use selected method for MBR_boosted_global_prob
        if selected_method == :ml
            precursor_features.MBR_boosted_global_prob = ml_predictions
        else
            precursor_features.MBR_boosted_global_prob = precursor_features.logodds_score
        end

        # Create mapping and broadcast MBR_boosted_global_prob
        precursor_to_mbr_global_prob = Dict(
            row.precursor_idx => row.MBR_boosted_global_prob
            for row in eachrow(precursor_features)
        )
        merged_df.MBR_boosted_global_prob = [
            precursor_to_mbr_global_prob[idx] for idx in merged_df.precursor_idx
        ]

        # Calculate standard global_prob using simple logodds on prec_prob
        @user_info "  Calculating standard global_prob using logodds on prec_prob"
        transform!(groupby(merged_df, :precursor_idx),
                   :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
    else
        # MBR is OFF: Use ML for global_prob
        @user_info "Step 5: MBR disabled - calculating global_prob with ML"

        # Build features from regular prec_prob
        precursor_features = build_precursor_features_df(
            merged_df, n_top_scores, sqrt_n_runs, :prec_prob
        )

        # Train ML models
        model_fold_0, model_fold_1, feature_cols = train_global_prob_models(
            precursor_features, params
        )

        # Get predictions and select best method
        ml_predictions = Vector{Float32}(undef, 0)
        selected_method = :logodds  # Default fallback

        if !isnothing(model_fold_0) && !isnothing(model_fold_1)
            ml_predictions = predict_global_probs_cv(
                precursor_features, model_fold_0, model_fold_1, feature_cols
            )
            selected_method = select_best_global_prob_method(
                precursor_features, ml_predictions
            )
        else
            @user_warn "ML model training failed - falling back to logodds"
        end

        # Use selected method for global_prob
        if selected_method == :ml
            precursor_features.global_prob = ml_predictions
        else
            precursor_features.global_prob = precursor_features.logodds_score
        end

        # Create mapping and broadcast global_prob
        precursor_to_global_prob = Dict(
            row.precursor_idx => row.global_prob
            for row in eachrow(precursor_features)
        )
        merged_df.global_prob = [
            precursor_to_global_prob[idx] for idx in merged_df.precursor_idx
        ]
    end

    # Write updated data back to individual files
    for (file_idx, ref) in zip(valid_file_indices, filtered_refs)
        sub_df = merged_df[merged_df.ms_file_idx .== file_idx, :]
        write_arrow_file(ref, sub_df)
    end
end
```

---

### Phase 6: Parameter Addition

#### Add to `ScoringSearchParameters` struct

**File**: `ScoringSearch.jl` around line 64

**Add field**:
```julia
# Global probability ML parameters
n_top_scores_global_prob::Int
```

**Add to constructor** (around line 110):
```julia
# Global probability ML parameter (with default)
Int64(get(ml_params, :n_top_scores_global_prob, 3))
```

#### Add to `defaultSearchParams.json`

**File**: `assets/example_config/defaultSearchParams.json`

**Add to `machine_learning` section**:
```json
"n_top_scores_global_prob": 3
```

---

### Phase 7: Testing and Validation

#### Unit Tests
Create `test/UnitTests/ScoringSearch/test_global_prob_ml.jl`:
```julia
@testset "Global Probability ML" begin
    # Test feature calculation
    @testset "Feature Calculation" begin
        # Test calculate_top_n_features with n_obs < N, = N, > N
        # Test padding with 0.0 when n_obs < n_top_scores
        # Test that statistics match expected values
    end

    # Test feature DataFrame construction
    @testset "Feature Engineering" begin
        # Test build_precursor_features_df with mock merged_df
        # Test groupby processing
        # Test correct feature assignment
    end

    # Test AUC calculation
    @testset "AUC Calculation" begin
        # Test with known AUC values
        # Test edge cases (all targets, all decoys)
    end

    # Test model selection logic
    @testset "Model Selection" begin
        # Test AUC comparison logic
        # Test ML vs logodds selection
    end
end
```

#### Integration Test
Run full pipeline on E. coli test:
```julia
SearchDIA("test/test_config/ecoli_test_params.json")
```

**Validate**:
1. Models train without errors
2. AUC comparison executes correctly
3. global_prob values are in [0, 1] range
4. Downstream steps (protein inference) work correctly

---

## Implementation Checklist

### Files to Create
- [ ] `src/Routines/SearchDIA/SearchMethods/ScoringSearch/global_prob_ml.jl`
- [ ] `test/UnitTests/ScoringSearch/test_global_prob_ml.jl`

### Files to Modify
- [ ] `ScoringSearch.jl` - Step 5 replacement, parameter addition
- [ ] `ScoringSearch.jl` - Include new file in imports
- [ ] `defaultSearchParams.json` - Add `n_top_scores_global_prob` parameter
- [ ] `paramsChecks.jl` - Add validation for new parameter

### Implementation Order
1. **Phase 1**: Create structs and accumulation logic
2. **Phase 2**: Build feature engineering pipeline
3. **Phase 3**: Implement model training with CV
4. **Phase 4**: Add AUC calculation and model selection
5. **Phase 5**: Integrate into Step 5 of ScoringSearch
6. **Phase 6**: Add parameter to config files
7. **Phase 7**: Write tests and validate

---

## Design Decisions

### Using DataFrames groupby for Efficiency
**Decision**: Use `groupby(merged_df, :precursor_idx)` to process all prec_prob values per precursor
**Rationale**:
- **Idiomatic**: Matches existing pattern in Step 5 (`transform!(groupby(...))`)
- **Efficient**: DataFrames groupby is optimized for this use case
- **Simple**: No need for accumulators or running statistics since we get all values at once
- **Memory**: merged_df is already in memory, grouping just creates views
- **Logodds-friendly**: We need all values for logodds calculation anyway

**Why not OnlineStats?**
- We already have all prec_prob values in merged_df
- Logodds requires the full distribution anyway
- Simpler to compute mean/median/std directly on the vector
- No benefit from streaming when data is already loaded

### Why include logodds as a feature?
**Decision**: Add `logodds_score` as one of the ML model features
**Rationale**:
- ML model can learn when to trust logodds vs other features
- Provides a strong baseline feature
- Allows smooth degradation (if other features aren't informative, model uses logodds)

### Why compare on first fold only?
**Decision**: Use AUC on fold 0 for method selection
**Rationale**:
- Fast decision (no need to evaluate on all data)
- CV fold 0 should be representative
- Consistent with user requirement

### MBR Handling
**Decision**: Separate paths for MBR on/off with different column outputs

**When MBR is ON:**
1. Train ML model using `MBR_boosted_prec_prob` → output `MBR_boosted_global_prob`
2. Calculate `global_prob` from `prec_prob` using simple **logodds** (no ML)
3. **Result**: Both columns (`MBR_boosted_global_prob` and `global_prob`)

**When MBR is OFF:**
1. Train ML model using `prec_prob` → output `global_prob`
2. **Result**: Single column (`global_prob`)

**Rationale**:
- **Downstream compatibility**: Steps 6-10 expect `MBR_boosted_global_prob` when MBR is on
- **Protein inference**: Always needs standard `global_prob` from non-MBR scores
- **Correctness**: MBR-boosted scores get ML enhancement, while standard scores use simple logodds as baseline
- **Flexibility**: When MBR is off, full ML power is applied to standard prec_prob

---

## Expected Performance Impact

**Computational cost**:
- Feature engineering: ~0.5-1 second for 10K precursors (faster with OnlineStats)
- Model training: ~2-5 seconds for 2 models on 10K precursors
- AUC calculation: <1 second
- **Total overhead**: ~3-7 seconds for typical experiment

**Memory usage**:
- Precursor features DataFrame: ~400 KB for 10K precursors with 13 base features + top N
- No additional memory overhead (merged_df already loaded)
- groupby creates views, not copies

**Accuracy improvement**:
- Expected AUC improvement: +0.01 to +0.05 over logodds
- Benefit increases with larger experiments (more runs provide richer features)

---

## Future Enhancements (Out of Scope)

1. **Iterative training** (like PSM scoring) - add negative mining
2. **Additional features**:
   - RT variance across runs
   - Spectral angle variance
   - Mass error consistency
3. **Hyperparameter tuning** - optimize LightGBM parameters via grid search
