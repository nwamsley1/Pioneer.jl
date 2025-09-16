# ML-Based MBR Filtering Implementation Plan

## Overview

Replace the simple threshold-based MBR filtering with a machine learning model that predicts bad transfers from PSM features. This will provide more sophisticated filtering that can learn complex patterns in the data to identify risky MBR transfers.

## Current State Analysis

### Current Implementation (`get_ftr_threshold`)
- Uses simple quantile-based threshold on probability scores
- Binary decision: score >= threshold → keep, score < threshold → filter
- No consideration of feature-based patterns that indicate bad transfers

### Target Variable
- `is_bad_transfer`: Boolean indicating risky MBR transfers
- Current logic identifies bad transfers as:
  ```julia
  is_bad_transfer = candidate_mask .& (
      (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| # T->D
      (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false))   # D->T
  )
  ```

### Training Data
- **Include**: Only rows where `candidate_mask == true` (MBR transfer candidates)
- **Exclude**: Non-MBR candidates (already passed initial scoring)
- **Features**: Available PSM features in `merged_df`
- **Target**: `is_bad_transfer` binary classification

## Proposed Architecture

### 1. Feature Engineering
**Required Analysis**:
- Audit available features in `merged_df` at filtering time
- Identify MBR-specific features vs general PSM features
- Handle missing values and feature scaling

**Expected Feature Categories**:
- **Core PSM Features**: `prob`, `spectral_contrast`, `mass_error`, etc.
- **MBR-Specific Features**: `MBR_max_pair_prob`, `MBR_is_best_decoy`, `MBR_best_irt_diff`, etc.
- **RT/iRT Features**: Retention time predictions and errors
- **Intensity Features**: Signal quality and explained variance

### 2. Model Architecture
**Model Choice**: XGBoost/EvoTrees (consistent with existing codebase)
- Fast training for online filtering
- Handles mixed feature types well
- Interpretable feature importance
- Robust to missing values

**Model Configuration**:
```julia
config = EvoTreeClassifier(
    loss = :logloss,
    nrounds = 50,        # Fast training
    max_depth = 4,       # Prevent overfitting on small datasets
    eta = 0.1,           # Conservative learning rate
    subsample = 0.8,     # Bootstrap for robustness
    colsample = 0.8,     # Feature bagging
    gamma = 1.0          # Regularization
)
```

### 3. Implementation Strategy

#### Phase 1: Data Collection and Analysis
1. **Feature Audit**: Add logging to capture all available features
2. **Data Collection**: Save training examples across multiple runs
3. **EDA**: Analyze feature distributions and target balance
4. **Baseline**: Establish performance of current threshold method

#### Phase 2: Model Development
1. **Feature Selection**: Identify most predictive features
2. **Model Training**: Cross-validation with balanced sampling
3. **Validation**: Compare against current threshold method
4. **Hyperparameter Tuning**: Optimize for transfer candidate detection

#### Phase 3: Integration
1. **Online Training**: Train model within `apply_mbr_filter!`
2. **Scoring**: Replace threshold with model predictions
3. **Fallback**: Keep threshold method as backup
4. **Monitoring**: Log model performance metrics

## Detailed Implementation

### 1. Modified `apply_mbr_filter!` Function

```julia
function apply_mbr_filter!(
    merged_df::DataFrame,
    params,
    fdr_scale_factor::Float32,
)
    # Current logic for identifying candidates and bad transfers
    candidate_mask = merged_df.MBR_transfer_candidate
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| 
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false))
    )
    
    # NEW: ML-based filtering
    const MIN_MBR_CANDIDATES_FOR_ML = 20000  # Hardcoded threshold
    if sum(candidate_mask) >= MIN_MBR_CANDIDATES_FOR_ML
        filtered_prob_col = apply_ml_mbr_filter!(
            merged_df, candidate_mask, is_bad_transfer, params
        )
    else
        # Fallback to threshold method for small datasets
        @user_info "Using threshold-based MBR filtering ($(sum(candidate_mask)) candidates < $MIN_MBR_CANDIDATES_FOR_ML threshold)"
        filtered_prob_col = apply_threshold_mbr_filter!(
            merged_df, candidate_mask, is_bad_transfer, params
        )
    end
    
    merged_df._filtered_prob = filtered_prob_col
    return :_filtered_prob
end
```

### 2. ML Filter Implementation

```julia
function apply_ml_mbr_filter!(
    merged_df::DataFrame,
    candidate_mask::AbstractVector{Bool},
    is_bad_transfer::AbstractVector{Bool},
    params
)
    # Extract training data
    training_data = merged_df[candidate_mask, :]
    y_train = is_bad_transfer[candidate_mask]
    
    # Feature selection and preprocessing
    feature_cols = select_mbr_features(training_data)
    X_train = prepare_features(training_data[:, feature_cols])
    
    # Train model
    model = train_mbr_filter_model(X_train, y_train, params)
    
    # Generate predictions for all candidates
    bad_transfer_scores = predict(model, X_train)
    
    # Determine threshold based on desired FTR
    τ = calibrate_ml_threshold(
        bad_transfer_scores, y_train, params.max_MBR_false_transfer_rate
    )
    
    @user_info "ML MBR Filter: τ = $τ, $(sum(bad_transfer_scores .>= τ)) candidates filtered"
    
    # Apply filtering
    filtered_probs = copy(merged_df.prob)
    candidate_indices = findall(candidate_mask)
    
    for (i, idx) in enumerate(candidate_indices)
        if bad_transfer_scores[i] >= τ
            filtered_probs[idx] = 0.0f0
        end
    end
    
    return filtered_probs
end
```

### 3. Feature Engineering

```julia
function select_mbr_features(df::DataFrame)
    # Core features always available
    core_features = [:prob, :spectral_contrast, :mass_error, :rt_error]
    
    # MBR-specific features (if available)
    mbr_features = [
        :MBR_max_pair_prob, :MBR_best_irt_diff, :MBR_rv_coefficient,
        :MBR_log2_weight_ratio, :MBR_log2_explained_ratio, :MBR_num_runs
    ]
    
    # Additional features based on availability
    optional_features = [
        :intensity_explained, :fragment_coverage, :charge, :length,
        :missed_cleavages, :irt_error, :precursor_intensity
    ]
    
    # Filter to available columns
    available_features = Symbol[]
    for feature_set in [core_features, mbr_features, optional_features]
        for feature in feature_set
            if hasproperty(df, feature)
                push!(available_features, feature)
            end
        end
    end
    
    @user_info "Selected $(length(available_features)) features for MBR ML filtering"
    return available_features
end

function prepare_features(df::DataFrame)
    # Handle missing values
    for col in names(df)
        col_data = df[!, col]
        if eltype(col_data) <: Union{Missing, <:Number}
            # Replace missing with median for numeric
            median_val = median(skipmissing(col_data))
            df[!, col] = coalesce.(col_data, median_val)
        end
    end
    
    # Convert to Matrix for ML training
    return Matrix{Float64}(df)
end
```

### 4. Model Training and Calibration

```julia
function train_mbr_filter_model(X::Matrix, y::AbstractVector{Bool}, params)
    config = EvoTreeClassifier(
        loss = :logloss,
        nrounds = params.mbr_ml_nrounds,
        max_depth = params.mbr_ml_max_depth,
        eta = params.mbr_ml_eta,
        subsample = params.mbr_ml_subsample,
        colsample = params.mbr_ml_colsample,
        gamma = params.mbr_ml_gamma
    )
    
    # Handle class imbalance if needed
    if mean(y) < 0.1 || mean(y) > 0.9
        @user_warn "Imbalanced classes in MBR training: $(mean(y)*100)% positive"
    end
    
    model = fit(config, X, y; verbosity=0)
    
    # Feature importance logging
    importance = EvoTrees.importance(model)
    @user_info "Top MBR filter features: $(first(importance, 5))"
    
    return model
end

function calibrate_ml_threshold(scores::AbstractVector, y_true::AbstractVector{Bool}, target_ftr::Float64)
    # Find threshold that achieves target false transfer rate
    # Similar to get_qvalues! logic but for binary classification
    
    sorted_indices = sortperm(scores, rev=true)
    sorted_scores = scores[sorted_indices]
    sorted_labels = y_true[sorted_indices]
    
    cumulative_bad = cumsum(sorted_labels)
    cumulative_total = 1:length(sorted_labels)
    
    # False transfer rate = bad_transfers / total_transfers
    ftr = cumulative_bad ./ cumulative_total
    
    # Find first point where FTR <= target
    valid_idx = findfirst(x -> x <= target_ftr, ftr)
    
    if isnothing(valid_idx)
        @user_warn "Could not achieve target FTR $(target_ftr), using most restrictive threshold"
        return maximum(sorted_scores)
    end
    
    threshold = sorted_scores[valid_idx]
    achieved_ftr = ftr[valid_idx]
    
    @user_info "MBR ML threshold: $(threshold), achieved FTR: $(achieved_ftr)"
    return threshold
end
```

## Configuration Parameters

**Hardcoded Constants:**
```julia
const MIN_MBR_CANDIDATES_FOR_ML = 20000    # Minimum candidates to use ML (hardcoded)
```

**Future Parameters** (for later configuration):
```julia
# MBR ML Filtering Parameters (not implemented yet)
mbr_ml_nrounds::Int = 50                    # XGBoost rounds
mbr_ml_max_depth::Int = 4                   # Tree depth
mbr_ml_eta::Float64 = 0.1                   # Learning rate
mbr_ml_subsample::Float64 = 0.8             # Row sampling
mbr_ml_colsample::Float64 = 0.8             # Feature sampling
mbr_ml_gamma::Float64 = 1.0                 # Regularization
```

## Validation and Testing

### 1. Performance Metrics
- **Precision**: True bad transfers / Predicted bad transfers
- **Recall**: True bad transfers / Total bad transfers  
- **F1-Score**: Harmonic mean of precision and recall
- **FTR Control**: Achieved vs target false transfer rate
- **Downstream Impact**: Effect on final PSM/protein counts

### 2. Comparison Studies
- **Baseline**: Current threshold method performance
- **ML Method**: Proposed implementation
- **Hybrid**: ML + threshold fallback
- **Cross-Dataset**: Validate on different experimental setups

### 3. Edge Cases
- **Small Datasets**: < 100 candidates
- **Extreme Imbalance**: Very few bad transfers
- **Missing Features**: Handle gracefully
- **Computational Time**: Ensure acceptable performance

## Migration Strategy

### Phase 1: Development (Week 1-2)
1. Implement feature logging and data collection
2. Develop training pipeline offline
3. Validate model performance on historical data

### Phase 2: Integration (Week 3)
1. Integrate ML filter into `apply_mbr_filter!`
2. Add parameter configuration
3. Implement fallback mechanisms

### Phase 3: Validation (Week 4)
1. A/B testing with threshold method
2. Performance benchmarking
3. Parameter tuning and optimization

### Phase 4: Deployment
1. Documentation and user guidelines
2. Default parameter optimization
3. Monitoring and maintenance procedures

## Risks and Mitigation

### Technical Risks
- **Overfitting**: Use cross-validation and regularization
- **Feature Drift**: Monitor feature distributions over time
- **Computational Cost**: Profile and optimize training time
- **Integration Bugs**: Comprehensive testing with edge cases

### Scientific Risks
- **Bias Introduction**: Validate that ML doesn't introduce systematic bias
- **Generalization**: Test across different experimental conditions
- **Interpretability**: Ensure feature importance makes biological sense

### Mitigation Strategies
- **Gradual Rollout**: Optional parameter to enable/disable
- **Fallback System**: Always maintain threshold method as backup
- **Monitoring**: Log performance metrics for continuous validation
- **User Control**: Allow manual override of ML decisions

## Expected Benefits

1. **Improved Precision**: Better identification of truly bad transfers
2. **Adaptive Filtering**: Learns from data-specific patterns
3. **Feature Utilization**: Leverages all available PSM information
4. **Reduced False Positives**: More sophisticated decision boundary
5. **Interpretability**: Feature importance provides insights into transfer quality

## Success Metrics

- **Quantitative**: 10-20% improvement in precision at same recall
- **Qualitative**: Reduced manual review of filtered candidates
- **Robustness**: Consistent performance across different datasets
- **Efficiency**: Training time < 5 seconds for typical datasets
- **Adoption**: Positive user feedback and widespread usage