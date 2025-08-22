# Model Comparison Implementation Plan for ScoringSearch In-Memory Approach

## Overview

This document outlines the detailed implementation plan for adding model comparison capabilities to the ScoringSearch module when using the in-memory approach (< 100,000 PSMs). The system will test three models using an 80/20 train-validation split to select the best performer based on the number of targets passing q-value ≤ 0.01 threshold, then proceed with the selected model on the full dataset.

## Current Implementation Analysis

### Existing Models

1. **Complex XGBoost Model** (line 440-456 in score_psms.jl)
   - **Features**: 50+ features including spectral scores, RT errors, MS1 features, MBR features
   - **Hyperparameters**: colsample_bytree=0.5, eta=0.05, max_depth=10, iter_scheme=[100,200,200]
   - **Use case**: Large datasets (>100k PSMs) with comprehensive feature set

2. **Simple XGBoost Model** (line 609-626 in score_psms.jl)  
   - **Features**: Reduced feature set (~35 features, same as complex but different hyperparameters)
   - **Hyperparameters**: colsample_bytree=0.8, eta=0.1, max_depth=4, iter_scheme=[150,300,300]
   - **Use case**: Small datasets (<100k PSMs) with conservative regularization

3. **Probit Regression Model** (line 214-358 in score_psms.jl)
   - **Features**: Same feature set as XGBoost models but linear combination
   - **Algorithm**: Cross-validation probit regression with no iterative refinement
   - **Use case**: Very small datasets or when XGBoost fails

4. **Super Simplified Model** (line 537-545 in score_psms.jl)
   - **Features**: Only 5 core features: fitted_spectral_contrast, max_matched_residual, max_unmatched_residual, err_norm, log2_intensity_explained
   - **Algorithm**: Uses existing XGBoost or probit with minimal feature set
   - **Use case**: Very small datasets with limited feature availability

### Current Decision Logic

- **≥100k PSMs**: Use out-of-memory complex XGBoost model  
- **<100k PSMs**: Use in-memory approach with hard-coded model selection
  - Currently chooses simple XGBoost in most cases
  - Falls back to super simplified model for very small datasets

## Proposed Model Comparison Framework

### Model Configuration

#### Model 1: Simple XGBoost (Current Default)
```julia
ModelConfig(
    name = "SimpleXGBoost", 
    model_type = :xgboost,
    features = REDUCED_FEATURE_SET,  # ~35 key features (current implementation)
    hyperparams = Dict(
        :colsample_bytree => 0.8,
        :colsample_bynode => 0.8,
        :min_child_weight => 20,
        :gamma => 0.1,
        :subsample => 0.8,
        :max_depth => 4,
        :eta => 0.1,
        :iter_scheme => [150, 300, 300]
    )
)
```

#### Model 2: Probit Regression
```julia
ModelConfig(
    name = "ProbitRegression",
    model_type = :probit,
    features = REDUCED_FEATURE_SET,  # Same features as XGBoost for fair comparison
    hyperparams = Dict(
        :n_folds => 3,
        :max_iter => 30
    )
)
```

#### Model 3: Super Simplified Model
```julia
ModelConfig(
    name = "SuperSimplified",
    model_type = :xgboost,  # Use XGBoost with minimal features
    features = MINIMAL_FEATURE_SET,  # Only 5 core features
    hyperparams = Dict(
        :colsample_bytree => 0.8,
        :colsample_bynode => 0.8,
        :min_child_weight => 20,
        :gamma => 0.1,
        :subsample => 0.8,
        :max_depth => 4,
        :eta => 0.1,
        :iter_scheme => [150, 300, 300]
    )
)
```

### Feature Set Definitions

#### Reduced Feature Set (Simple XGBoost & Probit)
Full feature set from current in-memory implementation (~35 features):
```julia
REDUCED_FEATURE_SET = [
    # Core peptide properties
    :missed_cleavage, :Mox, :prec_mz, :sequence_length, :charge,
    # RT features
    :irt_pred, :irt_error, :irt_diff, :rt_diff, :ms1_irt_diff,
    # Spectral features
    :max_y_ions, :y_ions_sum, :longest_y, :y_count, :b_count, :isotope_count,
    :total_ions, :best_rank, :best_rank_iso, :topn, :topn_iso, :gof,
    # Quality metrics
    :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_residual, :max_unmatched_residual, :max_gof,
    :fitted_spectral_contrast, :spectral_contrast, :max_matched_ratio,
    :err_norm, :poisson, :weight, :log2_intensity_explained, :tic, :num_scans,
    :smoothness, :percent_theoretical_ignored, :scribe, :max_scribe,
    # MS1 features
    :weight_ms1, :gof_ms1, :max_matched_residual_ms1, :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1, :error_ms1, :m0_error_ms1, :n_iso_ms1,
    :big_iso_ms1, :ms1_features_missing
    # MBR features added automatically if match_between_runs=true
]
```

#### Minimal Feature Set (Super Simplified)
Only 5 core features from current super simplified implementation:
```julia
MINIMAL_FEATURE_SET = [
    :fitted_spectral_contrast,
    :max_matched_residual,
    :max_unmatched_residual,
    :err_norm,
    :log2_intensity_explained
]
```

## Implementation Plan

### Phase 1: Data Structures and Configuration

#### 1.1 Model Configuration Structure
```julia
struct ModelConfig
    name::String
    model_type::Symbol  # :xgboost or :probit
    features::Vector{Symbol}
    hyperparams::Dict{Symbol, Any}
end

struct ModelPerformance
    model_name::String
    n_targets_passing_qval::Int64  # Primary selection metric
    validation_auc::Float64
    validation_accuracy::Float64
    validation_sensitivity::Float64
    validation_specificity::Float64
    training_time::Float64
    n_features::Int
end
```

#### 1.2 Parameters Integration
Add to JSON parameters:
```julia
"machine_learning": {
    "enable_model_comparison": true,
    "validation_split_ratio": 0.2,
    "qvalue_threshold": 0.01,  # Q-value threshold for target counting
    "min_psms_for_comparison": 1000,
    "max_psms_for_comparison": 100000
}
```

### Phase 2: Data Splitting Framework

#### 2.1 Split Strategy
```julia
function create_train_validation_split(psms::DataFrame, 
                                     validation_ratio::Float64 = 0.2;
                                     seed::Int = 1776)
    """
    Creates simple 80/20 split with basic stratification by target/decoy status
    """
    
    Random.seed!(seed)
    
    # Simple stratification by target/decoy only
    target_indices = findall(psms.target)
    decoy_indices = findall(.!psms.target)
    
    # Split targets
    n_target_val = max(1, round(Int, length(target_indices) * validation_ratio))
    target_val = sample(target_indices, n_target_val, replace=false)
    target_train = setdiff(target_indices, target_val)
    
    # Split decoys  
    n_decoy_val = max(1, round(Int, length(decoy_indices) * validation_ratio))
    decoy_val = sample(decoy_indices, n_decoy_val, replace=false)
    decoy_train = setdiff(decoy_indices, decoy_val)
    
    train_indices = vcat(target_train, decoy_train)
    validation_indices = vcat(target_val, decoy_val)
    
    # Basic validation warnings
    if length(target_train) < 100
        @user_warn "Few target PSMs in training set: $(length(target_train))"
    end
    if length(decoy_train) < 100
        @user_warn "Few decoy PSMs in training set: $(length(decoy_train))"
    end
    if length(target_val) < 100
        @user_warn "Few target PSMs in validation set: $(length(target_val))"
    end
    if length(decoy_val) < 100
        @user_warn "Few decoy PSMs in validation set: $(length(decoy_val))"
    end
    
    return train_indices, validation_indices
end
```

#### 2.2 Split Summary (No Validation Required)
Basic logging of split characteristics is sufficient - no complex validation needed.

### Phase 3: Model Training Framework

#### 3.1 Unified Training Interface
```julia
function train_model(model_config::ModelConfig, 
                    psms_train::DataFrame,
                    match_between_runs::Bool) -> ModelResult
    """
    Trains a single model configuration and returns trained model + metadata
    """
    
    start_time = time()
    
    if model_config.model_type == :xgboost
        result = train_xgboost_model(model_config, psms_train, match_between_runs)
    elseif model_config.model_type == :probit
        result = train_probit_model(model_config, psms_train, match_between_runs)
    else
        error("Unknown model type: $(model_config.model_type)")
    end
    
    training_time = time() - start_time
    
    return ModelResult(
        model = result.model,
        model_config = model_config,
        training_time = training_time,
        n_features = length(model_config.features)
    )
end
```

#### 3.2 XGBoost Training Wrapper
```julia
function train_xgboost_model(config::ModelConfig, 
                           psms_train::DataFrame, 
                           match_between_runs::Bool) -> NamedTuple
    """
    Wrapper around existing sort_of_percolator_in_memory! for model comparison
    """
    
    # Filter features to those available and specified
    available_features = filter(f -> f in propertynames(psms_train), config.features)
    if match_between_runs
        mbr_features = [f for f in available_features if startswith(String(f), "MBR_")]
        append!(available_features, mbr_features)
    end
    
    # Create temporary file paths for this training (needed by existing function)
    temp_files = create_temporary_file_paths(psms_train)
    
    # Call existing training function with config hyperparameters
    models = sort_of_percolator_in_memory!(
        psms_train,
        temp_files,
        available_features,
        match_between_runs;
        config.hyperparams...
    )
    
    return (model = models, features = available_features)
end
```

#### 3.3 Probit Training Wrapper
```julia
function train_probit_model(config::ModelConfig, 
                          psms_train::DataFrame,
                          match_between_runs::Bool) -> NamedTuple
    """
    Wrapper around existing probit_regression_scoring_cv! for model comparison
    """
    
    # Create copy to avoid modifying original
    psms_copy = copy(psms_train)
    
    # Filter features
    available_features = filter(f -> f in propertynames(psms_copy), config.features)
    
    # Create temporary file paths
    temp_files = create_temporary_file_paths(psms_copy)
    
    # Call existing probit function
    probit_regression_scoring_cv!(
        psms_copy,
        temp_files, 
        available_features,
        match_between_runs;
        n_folds = config.hyperparams[:n_folds]
    )
    
    # Probit doesn't return models - extract coefficients if needed
    return (model = nothing, features = available_features, psms = psms_copy)
end
```

### Phase 4: Model Evaluation Framework

#### 4.1 Validation Performance Assessment
```julia
function evaluate_model_performance(model_result::ModelResult,
                                  psms_val::DataFrame,
                                  match_between_runs::Bool,
                                  qvalue_threshold::Float64 = 0.01) -> ModelPerformance
    """
    Evaluates trained model on validation set and computes performance metrics.
    Primary metric: number of targets passing q-value threshold
    """
    
    # Generate predictions on validation set
    if model_result.model_config.model_type == :xgboost
        val_predictions = predict_xgboost_validation(model_result, psms_val)
    elseif model_result.model_config.model_type == :probit
        val_predictions = predict_probit_validation(model_result, psms_val, match_between_runs)
    end
    
    # Compute q-values for validation predictions
    val_qvalues = zeros(Float64, length(val_predictions))
    get_qvalues!(val_predictions, psms_val.target, val_qvalues)
    
    # Primary selection metric: number of targets passing q-value threshold
    n_targets_passing = sum((psms_val.target) .& (val_qvalues .<= qvalue_threshold))
    
    # Secondary metrics for logging
    auc = compute_auc(val_predictions, psms_val.target)
    accuracy = compute_accuracy(val_predictions, psms_val.target)
    sensitivity, specificity = compute_sensitivity_specificity(val_predictions, psms_val.target)
    
    return ModelPerformance(
        model_result.model_config.name,
        n_targets_passing,  # Primary selection metric
        auc,
        accuracy, 
        sensitivity,
        specificity,
        model_result.training_time,
        model_result.n_features
    )
end
```

#### 4.2 Prediction Functions
```julia
function predict_xgboost_validation(model_result::ModelResult, 
                                  psms_val::DataFrame) -> Vector{Float32}
    """
    Generates XGBoost predictions on validation set using CV models
    """
    
    # Extract CV models from result
    cv_models = model_result.model
    features = model_result.trained_features
    
    # Predict using CV fold assignment
    predictions = zeros(Float32, nrow(psms_val))
    for (fold_idx, model) in cv_models
        fold_mask = psms_val.cv_fold .== fold_idx
        if any(fold_mask)
            predictions[fold_mask] = predict(model[end], psms_val[fold_mask, :])  # Use final iteration model
        end
    end
    
    return predictions
end

function predict_probit_validation(model_result::ModelResult,
                                 psms_val::DataFrame, 
                                 match_between_runs::Bool) -> Vector{Float32}
    """
    For probit, re-train on training set and predict validation
    (since probit doesn't store models, need to retrain)
    """
    
    # Alternative: Store probit coefficients during training phase
    # For now, extract predictions from probability column of trained psms
    # This requires modification to return trained coefficients
    
    error("Probit validation prediction not yet implemented - need to store coefficients")
end
```

#### 4.3 Performance Metrics
```julia
function compute_auc(predictions::Vector{Float32}, targets::Vector{Bool}) -> Float64
    """Uses existing ROC/AUC computation from Pioneer"""
    # Implementation using existing ROC utilities
end

function compute_accuracy(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0) -> Float64
    predicted_labels = predictions .> threshold
    return mean(predicted_labels .== targets)
end

function compute_sensitivity_specificity(predictions::Vector{Float32}, targets::Vector{Bool}, threshold::Float32 = 0.5f0) -> Tuple{Float64, Float64}
    predicted_labels = predictions .> threshold
    
    tp = sum(predicted_labels .& targets)
    fn = sum(.!predicted_labels .& targets) 
    tn = sum(.!predicted_labels .& .!targets)
    fp = sum(predicted_labels .& .!targets)
    
    sensitivity = tp / (tp + fn)  # True positive rate
    specificity = tn / (tn + fp)  # True negative rate
    
    return sensitivity, specificity
end
```

### Phase 5: Model Selection and Integration

#### 5.1 Model Selection Logic
```julia
function select_best_model(performances::Vector{ModelPerformance}) -> String
    """
    Selects best model based on number of targets passing q-value threshold.
    Tie-breaking: highest AUC, then fastest training time.
    """
    
    # Primary metric: number of targets passing q-value threshold
    max_targets = maximum([p.n_targets_passing_qval for p in performances])
    best_candidates = filter(p -> p.n_targets_passing_qval == max_targets, performances)
    
    if length(best_candidates) == 1
        best_model = best_candidates[1]
    else
        # Tie-breaking: highest AUC
        max_auc = maximum([p.validation_auc for p in best_candidates])
        auc_candidates = filter(p -> p.validation_auc == max_auc, best_candidates)
        
        if length(auc_candidates) == 1
            best_model = auc_candidates[1]
        else
            # Final tie-breaking: fastest training time
            best_model = auc_candidates[argmin([p.training_time for p in auc_candidates])]
        end
    end
    
    @user_info "Selected $(best_model.model_name) with $(best_model.n_targets_passing_qval) targets passing q-value threshold"
    
    return best_model.model_name
end
```

#### 5.2 Full Dataset Training
```julia
function train_selected_model_full_dataset(best_model_name::String,
                                         model_configs::Vector{ModelConfig},
                                         psms_full::DataFrame,
                                         file_paths::Vector{String},
                                         match_between_runs::Bool,
                                         max_q_value_xgboost_rescore::Float32,
                                         max_q_value_xgboost_mbr_rescore::Float32,
                                         min_PEP_neg_threshold_xgboost_rescore::Float32) -> Any
    """
    Trains the selected best model on the full dataset (100% of PSMs) using existing procedures
    """
    
    selected_config_idx = findfirst(c -> c.name == best_model_name, model_configs)
    config = model_configs[selected_config_idx]
    
    @user_info "Training selected model $(config.name) on full dataset ($(nrow(psms_full)) PSMs)"
    
    # Use existing training procedures exactly as implemented
    if config.model_type == :xgboost
        # Call existing sort_of_percolator_in_memory! with selected config
        return sort_of_percolator_in_memory!(
            psms_full,
            file_paths,
            config.features,
            match_between_runs;
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore,
            config.hyperparams...
        )
    elseif config.model_type == :probit
        # Call existing probit_regression_scoring_cv! procedure
        probit_regression_scoring_cv!(
            psms_full,
            file_paths,
            config.features,
            match_between_runs;
            n_folds = config.hyperparams[:n_folds]
        )
        return nothing  # Probit doesn't return models
    end
end
```

### Phase 6: Integration with Existing Pipeline

#### 6.1 Main Entry Point Modification
```julia
function score_precursor_isotope_traces_in_memory_with_comparison!(
    best_psms::DataFrame,
    file_paths::Vector{String}, 
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    enable_model_comparison::Bool = false,
    validation_split_ratio::Float64 = 0.2,
    qvalue_threshold::Float64 = 0.01
)
    """
    Enhanced version of score_precursor_isotope_traces_in_memory! with model comparison.
    Only applies to in-memory approach (<100k PSMs).
    """
    
    n_psms = size(best_psms, 1)
    
    # Check if model comparison should be enabled (in-memory only)
    if !enable_model_comparison || n_psms < 1000 || n_psms >= 100000
        @user_info "Model comparison disabled - falling back to existing logic"
        return score_precursor_isotope_traces_in_memory!(
            best_psms, file_paths, precursors, match_between_runs,
            max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore, 
            min_PEP_neg_threshold_xgboost_rescore
        )
    end
    
    @user_info "Starting model comparison with $(n_psms) PSMs (in-memory approach)"
    
    # Phase 1: Create train/validation split
    train_indices, val_indices = create_train_validation_split(best_psms, validation_split_ratio)
    psms_train = best_psms[train_indices, :]
    psms_val = best_psms[val_indices, :]
    
    @user_info "Split: $(nrow(psms_train)) training, $(nrow(psms_val)) validation PSMs"
    
    # Phase 2: Define model configurations  
    model_configs = create_model_configurations()
    
    # Phase 3: Train all models on training set
    model_results = ModelResult[]
    for config in model_configs
        @user_info "Training model: $(config.name)"
        result = train_model(config, psms_train, match_between_runs)
        push!(model_results, result)
    end
    
    # Phase 4: Evaluate all models on validation set
    performances = ModelPerformance[]
    for result in model_results
        @user_info "Evaluating model: $(result.model_config.name)"
        perf = evaluate_model_performance(result, psms_val, match_between_runs, qvalue_threshold)
        push!(performances, perf)
    end
    
    # Phase 5: Select best model and log results
    log_model_comparison_results(performances)
    best_model_name = select_best_model(performances)
    
    # Phase 6: Train selected model on full dataset using existing procedures
    @user_info "Training $(best_model_name) on full dataset"
    models = train_selected_model_full_dataset(
        best_model_name, model_configs, best_psms, file_paths, 
        match_between_runs, max_q_value_xgboost_rescore,
        max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
    )
    
    return models
end
```

#### 6.2 Results Logging and Diagnostics
```julia
function log_model_comparison_results(performances::Vector{ModelPerformance})
    """
    Logs detailed comparison results for all models
    """
    
    @user_info "Model Comparison Results:"
    @user_info "========================"
    
    # Sort by primary metric (targets passing q-value threshold)
    sorted_perfs = sort(performances, by=p->p.n_targets_passing_qval, rev=true)
    
    for (i, perf) in enumerate(sorted_perfs)
        @user_info "$(i). $(perf.model_name):"
        @user_info "   Targets Passing Q≤0.01: $(perf.n_targets_passing_qval)"
        @user_info "   AUC: $(round(perf.validation_auc, digits=4))"
        @user_info "   Accuracy: $(round(perf.validation_accuracy, digits=4))"  
        @user_info "   Sensitivity: $(round(perf.validation_sensitivity, digits=4))"
        @user_info "   Specificity: $(round(perf.validation_specificity, digits=4))"
        @user_info "   Training Time: $(round(perf.training_time, digits=2))s"
        @user_info "   Features: $(perf.n_features)"
        @user_info ""
    end
    
    # Also write to file for detailed analysis
    write_model_comparison_report(performances)
end

function write_model_comparison_report(performances::Vector{ModelPerformance}, 
                                     output_dir::String)
    """
    Writes detailed CSV report of model comparison results
    """
    
    report_path = joinpath(output_dir, "model_comparison_report.csv")
    
    df = DataFrame(
        model_name = [p.model_name for p in performances],
        n_targets_passing_qval = [p.n_targets_passing_qval for p in performances],
        validation_auc = [p.validation_auc for p in performances],
        validation_accuracy = [p.validation_accuracy for p in performances],
        validation_sensitivity = [p.validation_sensitivity for p in performances], 
        validation_specificity = [p.validation_specificity for p in performances],
        training_time_seconds = [p.training_time for p in performances],
        n_features = [p.n_features for p in performances]
    )
    
    CSV.write(report_path, df)
    @user_info "Model comparison report written to: $report_path"
end
```

## Implementation Timeline

### Sprint 1 (Week 1): Core Infrastructure
- [ ] Implement ModelConfig and ModelPerformance data structures
- [ ] Add parameters to JSON schema and validation
- [ ] Implement data splitting framework with stratification
- [ ] Create basic model training wrappers

### Sprint 2 (Week 2): Model Training and Evaluation  
- [ ] Implement XGBoost training wrapper with hyperparameter passing
- [ ] Implement Probit training wrapper with coefficient extraction
- [ ] Build validation prediction and performance evaluation framework
- [ ] Implement performance metrics (AUC, accuracy, sensitivity/specificity)

### Sprint 3 (Week 3): Integration and Selection
- [ ] Implement model selection logic with configurable metrics
- [ ] Build full dataset training pipeline
- [ ] Integrate with existing score_precursor_isotope_traces_in_memory!
- [ ] Add comprehensive logging and reporting

### Sprint 4 (Week 4): Testing and Optimization
- [ ] Unit tests for all new components
- [ ] Integration tests with real data
- [ ] Performance optimization and memory profiling
- [ ] Documentation and user guide

## Testing Strategy

### Unit Tests
- [ ] Data splitting maintains stratification
- [ ] Model configurations are valid
- [ ] Performance metrics are accurate
- [ ] Model selection logic works correctly

### Integration Tests  
- [ ] Full pipeline with ecoli test data
- [ ] Different dataset sizes (1k, 10k, 50k PSMs)
- [ ] With and without MBR enabled
- [ ] Edge cases (very unbalanced data, single file)

### Performance Tests
- [ ] Memory usage profiling
- [ ] Runtime comparison vs existing approach
- [ ] Scalability with different thread counts

## Configuration Options

### JSON Parameters
```json
{
  "machine_learning": {
    "enable_model_comparison": true,
    "validation_split_ratio": 0.2,
    "qvalue_threshold": 0.01,
    "min_psms_for_comparison": 1000,
    "max_psms_for_comparison": 100000,
    "comparison_models": {
      "simple_xgboost": {
        "enabled": true, 
        "hyperparams": {
          "colsample_bytree": 0.8,
          "eta": 0.1,
          "max_depth": 4
        }
      },
      "probit_regression": {
        "enabled": true,
        "hyperparams": {
          "n_folds": 3,
          "max_iter": 30
        }
      },
      "super_simplified": {
        "enabled": true,
        "hyperparams": {
          "colsample_bytree": 0.8,
          "eta": 0.1,
          "max_depth": 4
        }
      }
    }
  }
}
```

## Risk Mitigation

### Potential Issues and Solutions

1. **Insufficient Validation Data**
   - **Risk**: Small datasets may not have enough validation data for reliable assessment
   - **Solution**: Minimum PSM thresholds and fallback to existing logic

2. **Training Time Overhead**
   - **Risk**: Training 3 models instead of 1 increases runtime
   - **Solution**: Parallel training where possible, time limits, user configuration

3. **Memory Usage**
   - **Risk**: Storing multiple models and datasets simultaneously
   - **Solution**: Sequential training/evaluation with cleanup, memory monitoring

4. **Model Selection Instability**
   - **Risk**: Different random splits might select different models
   - **Solution**: Fixed random seeds, multiple validation runs option

5. **Feature Compatibility**
   - **Risk**: Some features may not be available in all datasets
   - **Solution**: Feature availability checking, graceful degradation

## Success Metrics

### Technical Metrics
- [ ] Implementation completed on schedule
- [ ] All tests passing
- [ ] Memory usage within 150% of baseline
- [ ] Runtime overhead < 2x for model comparison

### Scientific Metrics
- [ ] Model selection identifies approach with most targets passing q-value threshold
- [ ] Consistent selection of best model when re-run with different random seeds
- [ ] Improved target identification on validation sets compared to single model
- [ ] Training time overhead acceptable (target: <2x baseline)

## Future Enhancements

### Phase 2 Features (Post-MVP)
- [ ] Ensemble model combining multiple approaches
- [ ] Automated hyperparameter optimization 
- [ ] Cross-validation model selection
- [ ] Feature importance analysis and selection
- [ ] Integration with AutoML frameworks
- [ ] Real-time model performance monitoring

This implementation plan provides a comprehensive roadmap for adding robust model comparison capabilities to Pioneer's ScoringSearch module while maintaining compatibility with existing workflows and ensuring scientific rigor through proper validation procedures.