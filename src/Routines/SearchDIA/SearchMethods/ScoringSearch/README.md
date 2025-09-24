# ScoringSearch Module

ScoringSearch is the 7th stage in the Pioneer DIA search pipeline, responsible for machine learning-based PSM rescoring, FDR control, and protein group analysis. This module features adaptive model selection and supports both in-memory and out-of-memory processing for datasets of varying sizes.

## Overview

ScoringSearch performs three main functions:
1. **PSM Scoring**: Machine learning models (LightGBM or Probit Regression) rescore PSMs
2. **FDR Control**: Calculate q-values and filter PSMs based on false discovery rate thresholds
3. **Protein Inference**: Group peptides into minimal protein sets and calculate protein-level scores

## Data Flow Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────────┐
│   SecondPass    │───▶│   ScoringSearch  │───▶│ IntegrateChromato-  │
│   PSM Files     │    │                  │    │ gramSearch          │
└─────────────────┘    └──────────────────┘    └─────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    ScoringSearch Pipeline                          │
│                                                                     │
│  ┌─────────────────┐  ┌──────────────────┐  ┌──────────────────┐   │
│  │  1. PSM Scoring │─▶│ 2. FDR Control   │─▶│ 3. Protein       │   │
│  │     & Model     │  │    & Filtering   │  │    Inference     │   │
│  │    Selection    │  │                  │  │                  │   │
│  └─────────────────┘  └──────────────────┘  └──────────────────┘   │
└─────────────────────────────────────────────────────────────────────┘
```

## Model Selection Strategy

ScoringSearch uses different strategies based on dataset size:

```
PSM Count Decision Tree:

                    Count PSMs
                        │
           ┌────────────┼────────────┐
           │            │            │
    < 1,000 PSMs   1K - 100K PSMs   > 100K PSMs
           │            │            │
           ▼            ▼            ▼
   SimpleLightGBM   Model Comparison   AdvancedLightGBM
     (Direct)      (4 models tested)    (Direct)
                        │
                        ▼
            ┌─────────────────────────────┐
            │  Test 4 Models:             │
            │  • SimpleLightGBM            │
            │  • AdvancedLightGBM         │
            │  • ProbitRegression        │
            │  • SuperSimplified         │
            │                            │
            │  Select best performer     │
            │  (most targets at q ≤ thresh) │
            └─────────────────────────────┘
```

### Memory Strategy

```
                    PSM Count vs Memory Limit
                              │
                 ┌────────────┼────────────┐
                 │            │            │
        < max_psms_in_memory  │     ≥ max_psms_in_memory
              (100K)          │           (100K)
                 │            │            │
                 ▼            │            ▼
          In-Memory Processing│    Out-of-Memory Processing
          (All PSMs loaded)   │    (Sample subset for training)
                              │
                    Uses same ModelConfig
                    and hyperparameters
```

## Model Configurations

### Available Models (model_config.jl)

| Model | Type | Features | Use Case | Hyperparameters |
|-------|------|----------|----------|-----------------|
| **SimpleLightGBM** | LightGBM | REDUCED_FEATURE_SET (40+ features) | Default for small datasets | Conservative (depth=4, eta=0.1) |
| **AdvancedLightGBM** | LightGBM | ADVANCED_FEATURE_SET (50+ features) | Default for large datasets | Aggressive (depth=10, eta=0.05) |
| **ProbitRegression** | Linear | REDUCED_FEATURE_SET | Fast alternative | Linear model (max_iter=30) |
| **SuperSimplified** | LightGBM | MINIMAL_FEATURE_SET (5 features) | Minimal overfitting | Conservative (depth=4, eta=0.1) |

### Feature Sets

- **ADVANCED_FEATURE_SET**: 50+ features including all spectral, RT, MS1, and quality metrics
- **REDUCED_FEATURE_SET**: 40+ core features for balanced performance  
- **MINIMAL_FEATURE_SET**: 5 essential features (spectral contrast, residuals, error norms, intensity explained)
- **MBR Features**: Automatically appended when match_between_runs=true

## Function Call Flow

### Main Entry Point

```julia
score_precursor_isotope_traces(
    second_pass_folder,
    file_paths,
    precursors,
    match_between_runs,
    max_q_value_lightgbm_rescore,
    max_q_value_lightgbm_mbr_rescore,
    min_PEP_neg_threshold_lightgbm_rescore,
    max_psms_in_memory,
    q_value_threshold
)
```

### Processing Flow Diagram

```
score_precursor_isotope_traces()
├── get_psms_count(file_paths)
├── if psms_count >= max_psms_in_memory:
│   ├── sample_psms_for_lightgbm()
│   ├── create_default_advanced_lightgbm_config()
│   └── score_precursor_isotope_traces_out_of_memory!()
│       └── sort_of_percolator_out_of_memory!()
│           └── [For each CV fold]
│               └── train_booster() → LightGBMModel
└── else (in-memory processing):
    ├── load_psms_for_lightgbm()
    ├── if psms_count >= 100K:
    │   └── create_default_advanced_lightgbm_config()
    └── else (< 100K):
        └── select_psm_scoring_model()
            ├── create_model_configurations()
            └── [For each model config]
                ├── score_precursor_isotope_traces_in_memory()
                └── count_passing_targets()
    └── score_precursor_isotope_traces_in_memory()
        ├── if model_config.model_type == :lightgbm:
        │   └── train_lightgbm_model_in_memory()
        │       └── sort_of_percolator_in_memory!()
        │           └── [For each CV fold × iteration]
        │               └── train_booster() → LightGBMModel
        └── elif model_config.model_type == :probit:
            └── train_probit_model_in_memory()
                └── probit_regression_scoring_cv!()
                    └── [For each CV fold]
                        └── Pioneer.ProbitRegression()
```

## LightGBM Training Pipeline (percolatorSortOf.jl)

### sort_of_percolator_in_memory!() Flow

```
sort_of_percolator_in_memory!()
├── Sort PSMs by [pair_id, isotopes_captured]
├── Initialize probability arrays and CV fold mappings
├── For each CV fold:
│   ├── Split train/test data
│   ├── For each iteration in iter_scheme [100, 200, 200]:
│   │   ├── get_training_data_for_iteration!()
│   │   │   ├── First iteration: Use all training PSMs
│   │   │   └── Later iterations: 
│   │   │       ├── Convert worst targets to decoys (PEP filtering)
│   │   │       └── Filter to high-confidence PSMs (q-value filtering)
│   │   ├── train_booster()
│   │   │   └── LightGBMModel(hyperparams...)
│   │   ├── Predict on train/test sets
│   │   ├── Calculate q-values
│   │   └── update_mbr_features!() [if match_between_runs]
│   └── Store fold predictions
├── Handle MBR transfer candidates [if match_between_runs]
└── Return trained models (Dict{UInt8, LightGBMModelVector})
```

### Key Training Features

1. **Cross-Validation**: Uses CV folds from spectral library (typically 2 folds with values 0,1)
2. **Iterative Training**: 3-stage scheme [100, 200, 200] rounds with progressive filtering
3. **Negative Mining**: Converts low-scoring targets to negatives in later iterations
4. **MBR Integration**: Optional match-between-runs features added in final iteration
5. **Dynamic Training Data**: Training set filtered by q-value thresholds between iterations

## Probit Regression Pipeline

### probit_regression_scoring_cv!() Flow

```
probit_regression_scoring_cv!()
├── Validate CV fold column exists
├── Initialize probability estimates array
├── Filter features to available columns
├── For each CV fold (0, 1):
│   ├── Split train/test data
│   ├── Check sufficient training data (≥100 targets, ≥100 decoys)
│   ├── Pioneer.ProbitRegression()
│   │   └── Linear probit model with logistic link
│   ├── Pioneer.ModelPredictProbs!()
│   │   └── Predict probabilities on test set
│   └── Store fold predictions
├── Assign final probabilities to DataFrame
├── Create MBR columns for compatibility [if match_between_runs]
├── Validate and clamp probability values
└── Clean up temporary columns
```

### Probit vs LightGBM Differences

| Aspect | LightGBM | Probit Regression |
|--------|---------|------------------|
| **Model Type** | Gradient boosted trees | Linear logistic model |
| **Training** | Multi-stage iterative | Single-pass |
| **Features** | Non-linear interactions | Linear combinations |
| **Speed** | Slower (minutes) | Faster (seconds) |
| **Overfitting** | Can overfit small data | More robust |
| **Interpretability** | Feature importance | Linear coefficients |
| **MBR Handling** | Separate MBR iterations | Creates compatibility columns |
| **Hyperparameters** | Many (depth, eta, etc.) | Few (max_iter) |

## Integration with ScoringSearch Pipeline

### PSM Processing Output

Both LightGBM and Probit models produce standardized outputs:

```julia
# Required columns for downstream processing
psms[!, :prob]                    # Primary PSM probabilities
psms[!, :MBR_prob]               # MBR probabilities (if applicable)
psms[!, :MBR_is_best_decoy]      # MBR transfer metadata
psms[!, :MBR_transfer_candidate] # MBR filtering flags
```

### 23-Step ScoringSearch Pipeline Integration

The trained models integrate into the broader ScoringSearch pipeline:

1. **Steps 1-3**: Model training and PSM scoring (this README's focus)
2. **Steps 4-10**: PSM processing, merging, and FDR control
3. **Steps 11-23**: Protein inference and scoring

## Configuration Parameters

### Key Parameters

```json
{
    "optimization": {
        "machine_learning": {
            "max_psms_in_memory": 100000,
            "enable_model_comparison": true,
            "max_q_value_lightgbm_rescore": 0.01,
            "max_q_value_lightgbm_mbr_rescore": 0.20,
            "min_PEP_neg_threshold_lightgbm_rescore": 0.90
        }
    },
    "global_settings": {
        "scoring": {
            "q_value_threshold": 0.01
        },
        "match_between_runs": true
    }
}
```

### Parameter Effects

- **max_psms_in_memory**: Determines in-memory vs out-of-memory processing
- **enable_model_comparison**: Controls automatic model selection for medium datasets
- **q_value_threshold**: Used for model performance evaluation during comparison
- **match_between_runs**: Enables MBR features and processing

## Performance Characteristics

### Model Training Times (Approximate)

| Dataset Size | Model | Training Time | Memory Usage |
|-------------|--------|---------------|--------------|
| 1K PSMs | ProbitRegression | ~5 seconds | Low |
| 1K PSMs | SimpleLightGBM | ~30 seconds | Low |
| 50K PSMs | AdvancedLightGBM | ~5 minutes | Medium |
| 200K PSMs | AdvancedLightGBM (OOM) | ~10 minutes | Constant |

### Scaling Behavior

- **In-Memory**: Linear scaling with PSM count until memory limit
- **Out-of-Memory**: Constant memory usage, longer I/O times
- **Model Comparison**: 4x training time overhead for comparison phase

## Files and Structure

```
ScoringSearch/
├── README.md                      # This documentation
├── ScoringSearch.jl              # Main 23-step pipeline orchestrator
├── score_psms.jl                 # PSM scoring entry point and model selection
├── model_config.jl               # Model configurations and feature sets
├── utils.jl                      # Legacy protein group analysis functions
├── utils_protein_ml.jl           # ML-enhanced protein scoring
├── protein_inference_pipeline.jl # Modern protein inference pipeline
└── scoring_interface.jl          # Type-safe file reference operations
```

### Related Files

```
src/utils/ML/
└── percolatorSortOf.jl          # Core LightGBM training implementation
```

## Debugging and Diagnostics

### Key Log Messages

- **Model Selection**: Shows comparison results and selected model
- **Training Progress**: Progress bars during model training (suppressed during comparison)
- **CV Fold Statistics**: Target/decoy counts per fold
- **Feature Importance**: Optional detailed feature rankings

### Common Issues

1. **Model Comparison Failures**: Check feature availability in PSM data
2. **Memory Issues**: Reduce max_psms_in_memory or use out-of-memory processing
3. **Poor Performance**: Inspect CV fold balance and feature distributions
4. **Probit Convergence**: Check for sufficient training data per fold

## Future Enhancements

Potential areas for improvement:

1. **Additional Models**: Support for other ML algorithms (Random Forest, Neural Networks)
2. **Hyperparameter Tuning**: Automatic hyperparameter optimization
3. **Feature Selection**: Automatic feature importance-based selection
4. **Ensemble Methods**: Combining multiple model predictions
5. **Online Learning**: Incremental model updates for very large datasets