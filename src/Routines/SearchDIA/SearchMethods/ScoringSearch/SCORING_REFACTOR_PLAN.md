# Plan: Refactor PSM Scoring Architecture

## Implementation Status: NOT STARTED

### ❌ ACTUAL STATE (After Analysis):
1. **Probit Regression NOT Fixed**: Still has `n_folds=3` parameter and creates its own CV folds (lines 261-270)
2. **No Constants Added**: `MAX_FOR_MODEL_SELECTION` constant doesn't exist in code
3. **No Model Selection Function**: `select_psm_scoring_model()` function doesn't exist
4. **Model Comparison Framework**: EXISTS and is complete in `model_comparison.jl`

### 🔴 CRITICAL ISSUES:
- Main entry point still uses `enable_model_comparison` parameter (should be automatic)
- Legacy function `score_precursor_isotope_traces_in_memory!` has conflicting model selection logic
- No unified execution function exists
- Implementation that was attempted was lost due to file corruption

## Problem Statement

The current PSM scoring system has several architectural issues:

### Current Issues:
1. **Duplicate Model Selection Logic**: `score_precursor_isotope_traces_in_memory!` has its own PSM count-based model selection (1K vs 100K vs >100K) that conflicts with the model comparison framework
2. **Three Separate Functions**: 
   - `score_precursor_isotope_traces_in_memory!` (with embedded model selection)
   - `score_precursor_isotope_traces_in_memory_with_comparison!` (with model comparison)
   - `score_precursor_isotope_traces_out_of_memory!` (for large datasets)
3. **Inconsistent Model Selection**: The standard function shouldn't do model selection since it only runs when psms_count > 100K (no model comparison)
4. **Code Duplication**: Feature definitions and model configurations are duplicated across functions
5. **Poor Separation of Concerns**: Model selection mixed with model execution

### Desired Architecture:
- **Single Execution Function**: One `score_precursor_isotope_traces_in_memory` that takes a model configuration
- **Separate Model Selection**: Dedicated function that chooses the model before execution
- **Clean Interface**: Model type specified explicitly, not inferred from PSM count

## Proposed Architecture

### New Function Structure:
```
score_precursor_isotope_traces()  # Main entry point
├── select_psm_scoring_model()    # NEW: Model selection logic
├── score_precursor_isotope_traces_in_memory()   # NEW: Unified in-memory execution 
└── score_precursor_isotope_traces_out_of_memory() # Existing: Large dataset handling
```

### Model Selection Flow (CORRECTED):
1. **Count PSMs** → Determine processing approach (in-memory vs out-of-memory)
2. **Select Model** → Choose model configuration based on corrected three-case logic
3. **Execute Model** → Run the selected model configuration
4. **Return Results** → Consistent output regardless of model used

### Corrected Three-Case Logic:
1. **PSMs ≥ max_psms_in_memory** (typically >> 100K):
   - Use **out-of-memory** processing
   - Use default/advanced XGBoost (no model selection)

2. **PSMs < max_psms_in_memory AND ≥ 100K**:
   - Use **in-memory** processing
   - Use default/advanced XGBoost (no model selection needed)
   - This is the standard high-performance path

3. **PSMs < 100K**:
   - Use **in-memory** processing  
   - Run **automatic model comparison** (SimpleXGBoost vs ProbitRegression vs SuperSimplified)
   - Select the best performing model

## Detailed Implementation Plan

### Phase 1: Fix Probit Regression CV Folds ❌ NOT STARTED

#### 1.1 Fix `probit_regression_scoring_cv!` function (lines 256-400)
**Required Changes**:
- Remove `n_folds::Int64 = 3` parameter from function signature (line 261)
- Remove CV fold assignment loop (lines 263-270)
- Add validation that `:cv_fold` column exists
- Update to handle library's 0,1 CV fold values (not 1,2,3)
- Update all fold iteration logic to use detected CV folds

#### 1.2 Update model_comparison.jl ProbitRegression configuration
**Required Changes**:
- Remove `:n_folds => 3` from ProbitRegression hyperparameters
- Update `train_probit_model` to not pass n_folds parameter

### Phase 2: Create Core Functions ❌ NOT STARTED

#### 2.1 Add MAX_FOR_MODEL_SELECTION constant
**Location**: After imports in score_psms.jl (line ~21)
```julia
const MAX_FOR_MODEL_SELECTION = 100_000
```

#### 2.2 Create `select_psm_scoring_model` function
**Location**: After constants in score_psms.jl
**Implementation**:
```julia
function select_psm_scoring_model(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    validation_split_ratio::Float64,
    qvalue_threshold::Float64,
    output_dir::String
)
    psms_count = size(best_psms, 1)
    
    if psms_count >= MAX_FOR_MODEL_SELECTION
        # Return default/advanced XGBoost config
        return create_default_advanced_xgboost_config()
    else
        # Run model comparison and return best model
        return run_model_comparison_and_select_best(
            best_psms, file_paths, precursors, match_between_runs,
            max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore, validation_split_ratio,
            qvalue_threshold, output_dir
        )
    end
end
```

#### 2.3 Create unified execution function

```julia
"""
    score_precursor_isotope_traces_in_memory(psms::DataFrame,
                                           file_paths::Vector{String},
                                           precursors::LibraryPrecursors,
                                           model_config::ModelConfig,
                                           match_between_runs::Bool,
                                           max_q_value_xgboost_rescore::Float32,
                                           max_q_value_xgboost_mbr_rescore::Float32,
                                           min_PEP_neg_threshold_xgboost_rescore::Float32) -> Models

Unified in-memory PSM scoring function that executes the specified model.

# Arguments
- `model_config`: ModelConfig specifying which model to use and its hyperparameters

# Supported Models:
- SimpleXGBoost: Full feature set, standard hyperparameters
- ProbitRegression: Linear probit model with CV folds from library
- SuperSimplified: Minimal 5-feature XGBoost model
"""
function score_precursor_isotope_traces_in_memory(...)
    # Route to appropriate model implementation
    if model_config.model_type == :xgboost
        return train_xgboost_model_in_memory(...)
    elseif model_config.model_type == :probit  
        return train_probit_model_in_memory(...)
    else
        error("Unsupported model type: $(model_config.model_type)")
    end
end
```

#### 2.4 Create helper functions

**`create_default_advanced_xgboost_config()`**:
```julia
function create_default_advanced_xgboost_config()
    return ModelConfig(
        "AdvancedXGBoost",
        :xgboost,
        REDUCED_FEATURE_SET,  # From model_comparison.jl
        Dict(
            :colsample_bytree => 0.5,
            :colsample_bynode => 0.5,
            :min_child_weight => 5,
            :gamma => 1.0,
            :subsample => 0.25,
            :max_depth => 10,
            :eta => 0.05,
            :iter_scheme => [100, 200, 200]
        )
    )
end
```

**`train_xgboost_model_in_memory()`**:
```julia
function train_xgboost_model_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features and hyperparams from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    if match_between_runs
        append!(features, [:MBR_rv_coefficient, :MBR_best_irt_diff, :MBR_num_runs,
                          :MBR_max_pair_prob, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio])
    end
    
    hp = model_config.hyperparams
    return sort_of_percolator_in_memory!(
        best_psms, file_paths, features, match_between_runs;
        max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
        min_PEP_neg_threshold_xgboost_rescore,
        colsample_bytree = hp[:colsample_bytree],
        colsample_bynode = hp[:colsample_bynode],
        min_child_weight = hp[:min_child_weight],
        gamma = hp[:gamma],
        subsample = hp[:subsample],
        max_depth = hp[:max_depth],
        eta = hp[:eta],
        iter_scheme = hp[:iter_scheme]
    )
end
```

**`train_probit_model_in_memory()`**:
```julia
function train_probit_model_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    
    probit_regression_scoring_cv!(
        best_psms,
        file_paths,
        features,
        match_between_runs
        # Note: n_folds parameter removed after Phase 1 fix
    )
    
    # Write results back to files
    dropVectorColumns!(best_psms)
    for (ms_file_idx, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    
    return nothing  # Probit doesn't return models
end
```

### Phase 3: Refactor Main Entry Point ❌ NOT STARTED

#### 3.1 ❌ TODO: Update `score_precursor_isotope_traces()` function
**Current Status**: Still has old logic with `enable_model_comparison` parameter and conditional logic

**Current Issues**:
- Still has `enable_model_comparison` parameter (should be removed)
- Uses old conditional logic instead of corrected three-case logic
- Calls `score_precursor_isotope_traces_in_memory_with_comparison!` (should be removed)
- Calls old `score_precursor_isotope_traces_in_memory!` function

**Required Refactoring**:
```julia
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    max_psms_in_memory::Int64;
    validation_split_ratio::Float64 = 0.2,
    qvalue_threshold::Float64 = 0.01,
    output_dir::String = "."
)
    # Step 1: Count PSMs and determine processing approach
    psms_count = get_psms_count(file_paths)
    
    if psms_count >= max_psms_in_memory
        # Case 1: Out-of-memory processing with default XGBoost
        best_psms = sample_psms_for_xgboost(second_pass_folder, psms_count, max_psms_in_memory)
        return score_precursor_isotope_traces_out_of_memory!(...)
    else
        # In-memory processing - load PSMs first
        best_psms = load_psms_for_xgboost(second_pass_folder)
        
        if psms_count >= MAX_FOR_MODEL_SELECTION  # 100K
            # Case 2: In-memory with default/advanced XGBoost (no comparison)
            model_config = create_default_advanced_xgboost_config()
        else
            # Case 3: In-memory with automatic model comparison (<100K)
            model_config = select_psm_scoring_model(
                best_psms, file_paths, precursors, match_between_runs,
                max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
                min_PEP_neg_threshold_xgboost_rescore, validation_split_ratio,
                qvalue_threshold, output_dir
            )
        end
        
        @user_info "Selected PSM scoring model: $(model_config.name)"
        
        # Execute selected model using unified function
        return score_precursor_isotope_traces_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
        )
    end
end
```

**Key Changes Needed**:
1. Remove `enable_model_comparison` and related parameters
2. Implement corrected three-case logic
3. Call new unified execution function
4. Add `create_default_advanced_xgboost_config()` helper function

### Phase 4: Clean Up Legacy Functions ❌ NOT STARTED

#### 4.1 Remove Legacy Functions
**Functions to Remove**:
1. **`score_precursor_isotope_traces_in_memory!`** (lines 402-667 in score_psms.jl)
   - Has embedded model selection logic that conflicts with new architecture
   - Replace with new unified function

2. **Eventually: `score_precursor_isotope_traces_in_memory_with_comparison!`** (model_comparison.jl)
   - Keep temporarily for reference during implementation
   - Remove after confirming new implementation works

#### 4.2 Update Call Sites
- Main entry point will call new unified function
- No other call sites found in codebase

## Implementation Order

### Implementation Steps (commit after each):
1. **Make initial commit** - Save current state before changes
2. **Phase 1**: Fix probit regression CV folds
   - Update probit_regression_scoring_cv! 
   - Update model_comparison.jl
   - Commit: "fix: Update probit regression to use library-assigned CV folds"
3. **Phase 2**: Create core functions
   - Add MAX_FOR_MODEL_SELECTION constant
   - Add all new functions
   - Commit: "feat: Add model selection and unified execution functions"
4. **Phase 3**: Refactor main entry point
   - Update score_precursor_isotope_traces
   - Commit: "refactor: Update main entry point with three-case logic"
5. **Phase 4**: Clean up
   - Remove legacy function
   - Commit: "cleanup: Remove legacy score_precursor_isotope_traces_in_memory!"
6. **Test**: Verify implementation works

## Critical Dependencies

- **Model comparison logic** in `model_comparison.jl` is already implemented ✅
- **Probit regression CV fix** is already complete ✅
- **ModelConfig structures** are defined in `model_comparison.jl` ✅
- **Feature set constants** are defined in `model_comparison.jl` ✅

---

## Summary

This refactoring will create a clean PSM scoring architecture with:
- **Automatic model comparison** for datasets <100K PSMs (no enable/disable option)
- **Consistent protein-based CV folds** from library (values 0,1) for all models
- **Clean separation** between model selection and execution
- **Three-case logic** that correctly handles memory vs model selection thresholds:
  - ≥ max_psms_in_memory: Out-of-memory processing
  - < max_psms_in_memory AND ≥ 100K: In-memory with default XGBoost
  - < 100K: In-memory with automatic model comparison

**Current Status**: Plan ready for implementation. No code has been written yet.

## Benefits of This Refactoring

### 1. **Clean Separation of Concerns**
- Model selection logic isolated and testable
- Model execution is pure function of configuration
- No hidden model selection based on PSM count

### 2. **Eliminates Code Duplication**
- Single feature set definitions
- Unified model execution logic
- Shared hyperparameter configurations

### 3. **Improved Testability**
- Model selection can be tested independently
- Model execution can be tested with specific configurations
- Clear interfaces between components

### 4. **Better Maintainability**
- Adding new models only requires updating model factory
- Feature sets managed centrally
- Clear flow: select → configure → execute

### 5. **Consistent Behavior**
- No hidden model selection surprises
- Predictable behavior across dataset sizes
- Explicit model choice logging

## Implementation Order

1. **Phase 1**: Create model selection function
2. **Phase 2**: Create unified execution function  
3. **Phase 3**: Update main entry point
4. **Phase 4**: Implement model-specific functions
5. **Phase 5**: Centralize model configurations
6. **Phase 6**: Remove legacy functions and update call sites

## Backward Compatibility

- Main entry point `score_precursor_isotope_traces()` keeps same signature
- Behavior changes only in edge cases (removes hidden model selection)
- All model comparison functionality preserved and enhanced

## Testing Strategy

- Unit tests for model selection logic
- Unit tests for each model implementation
- Integration tests for full pipeline
- Regression tests to ensure same results as before