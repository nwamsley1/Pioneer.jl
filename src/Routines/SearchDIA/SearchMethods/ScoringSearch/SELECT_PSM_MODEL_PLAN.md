# Plan: Complete select_psm_scoring_model Function

## Problem Statement

The current `select_psm_scoring_model` function has several issues:
1. It calls `score_precursor_isotope_traces_in_memory_with_comparison!` which uses validation splits
2. It always returns a hardcoded SimpleXGBoost config regardless of actual comparison results
3. It doesn't properly integrate with the new unified architecture

## Solution Design

### Key Principles
- **No validation split**: Use full dataset for training each model
- **Training error for selection**: Count targets passing q-value threshold on training data
- **Clean separation**: Each model trains on a deepcopy to avoid side effects
- **Simple selection**: Choose model with highest number of passing targets

### Algorithm Flow

```
1. For PSMs ≥ 100K:
   → Return default advanced XGBoost config immediately

2. For PSMs < 100K:
   a. Get all model configurations (SimpleXGBoost, ProbitRegression, SuperSimplified)
   b. For each model:
      - Create deepcopy of PSM DataFrame
      - Train model using score_precursor_isotope_traces_in_memory
      - Count targets passing q-value threshold
      - Store performance metrics
   c. Select model with best performance
   d. Return winning ModelConfig
```

## Implementation Details

### 1. Updated `select_psm_scoring_model` Function

```julia
function select_psm_scoring_model(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    validation_split_ratio::Float64,  # Keep for API compatibility but unused
    qvalue_threshold::Float64,
    output_dir::String
)
    psms_count = size(best_psms, 1)
    
    if psms_count >= MAX_FOR_MODEL_SELECTION
        @user_info "Using default advanced XGBoost for $psms_count PSMs (≥ 100K)"
        return create_default_advanced_xgboost_config()
    else
        @user_info "Running model comparison for $psms_count PSMs (< 100K)"
        
        # Get model configurations
        model_configs = create_model_configurations()
        best_model_config = nothing
        best_target_count = 0
        
        for config in model_configs
            @user_info "Training $(config.name) model..."
            
            # Create deepcopy to avoid side effects
            psms_copy = deepcopy(best_psms)
            
            # Train model
            score_precursor_isotope_traces_in_memory(
                psms_copy, file_paths, precursors, config,
                match_between_runs, max_q_value_xgboost_rescore,
                max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
            )
            
            # Count passing targets
            target_count = count_passing_targets(file_paths, qvalue_threshold)
            @user_info "  $(config.name): $target_count targets at q-value ≤ $qvalue_threshold"
            
            # Update best model if this one is better
            if target_count > best_target_count
                best_target_count = target_count
                best_model_config = config
            end
        end
        
        @user_info "Selected best model: $(best_model_config.name) with $best_target_count passing targets"
        return best_model_config
    end
end
```

### 2. New Helper Function: `count_passing_targets`

```julia
"""
    count_passing_targets(file_paths::Vector{String}, qvalue_threshold::Float64) -> Int

Counts the number of target PSMs passing the q-value threshold across all files.

# Arguments
- `file_paths`: Vector of Arrow file paths containing scored PSMs
- `qvalue_threshold`: Q-value threshold for counting (typically 0.01)

# Returns
- Total count of target PSMs with q_value ≤ threshold
"""
function count_passing_targets(file_paths::Vector{String}, qvalue_threshold::Float64)
    total_count = 0
    
    for fpath in file_paths
        if !endswith(fpath, ".arrow")
            continue
        end
        
        # Read PSM file
        psms = DataFrame(Arrow.Table(fpath))
        
        # Count targets passing threshold
        # Use 'prob' column if q_value not yet calculated
        if :q_value in propertynames(psms)
            passing = psms[psms.target .& (psms.q_value .<= qvalue_threshold), :]
        elseif :prob in propertynames(psms)
            # If q-values not calculated, use probability threshold
            # This is a fallback - ideally q-values should be calculated
            passing = psms[psms.target .& (psms.prob .>= (1.0 - qvalue_threshold)), :]
        else
            @user_warn "No scoring column found in $fpath"
            continue
        end
        
        total_count += nrow(passing)
    end
    
    return total_count
end
```

## Benefits of This Approach

1. **Simplicity**: No complex validation splitting logic
2. **Fairness**: All models trained on same full dataset
3. **Direct optimization**: Selects model that actually produces most identifications
4. **No overfitting concerns**: Since we're selecting based on the same metric used in production
5. **Clean implementation**: Each model trains independently via deepcopy

## Potential Issues and Mitigations

### Issue 1: Memory usage from deepcopy
**Mitigation**: Only happens for <100K PSMs, so memory impact is limited

### Issue 2: Training time for 3 models
**Mitigation**: Only for small datasets (<100K), training is fast

### Issue 3: File I/O for counting targets
**Mitigation**: Arrow files are efficient, only happens 3 times total

## Testing Strategy

1. **Unit test**: Verify `count_passing_targets` correctly counts targets
2. **Integration test**: Ensure model selection chooses best performer
3. **Memory test**: Verify deepcopy doesn't cause memory issues
4. **File integrity**: Ensure original files aren't modified during comparison

## Migration Path

1. Implement new `select_psm_scoring_model` function
2. Add `count_passing_targets` helper
3. Test with small dataset
4. Eventually remove `score_precursor_isotope_traces_in_memory_with_comparison!`

## Code Cleanup

After implementation:
- Remove validation split logic from model_comparison.jl
- Keep ModelConfig structures and feature definitions
- Keep individual model training functions
- Remove complex comparison framework code