# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the PrecursorScoringSearch module.

## PrecursorScoringSearch Overview

PrecursorScoringSearch is the precursor-level rescoring stage in the Pioneer DIA search pipeline. It performs machine learning-based PSM rescoring and precursor-level FDR control. The module features adaptive model selection for datasets with 1,000-100,000 PSMs, automatically choosing between LightGBM, probit regression, and simplified models based on validation performance.

## Module Structure

```
PrecursorScoringSearch/
├── PrecursorScoringSearch.jl      # Main search method implementation
├── score_psms.jl                 # PSM scoring entry point with automatic in-memory/OOM selection and model comparison
├── model_config.jl               # Model configurations for comparison (SimpleLightGBM, AdvancedLightGBM, ProbitRegression, SuperSimplified)
├── utils.jl                      # Shared q-value / spline helpers
└── scoring_interface.jl          # Type-safe file reference operations
```

## Key Components

### Main Workflow (PrecursorScoringSearch.jl)

PrecursorScoringSearch implements the precursor rescoring and precursor-level filtering workflow:

**Phase 1: Model Training & PSM Scoring (Steps 1-3)**
1. **Model Training**: Adaptive model selection with comparison framework (when enabled)
2. **Probability Computation**: Calculate precursor probabilities with MBR filtering
3. **Best Trace Selection**: Identify optimal isotope traces per precursor

**Phase 2: PSM Processing & FDR Control (Steps 4-10)**
4. **Quantification Processing**: Apply best trace indicators and column filtering
5. **Merge by Global Probability**: Combine PSMs sorted by global_prob for global q-values
6. **Global Precursor Q-values**: Calculate experiment-wide FDR for unique precursors
7. **Re-merge by Precursor Probability**: Resort and merge by prec_prob for run-specific q-values
8. **Experiment-wide Precursor Q-values**: Calculate FDR across all precursor instances
9. **PSM Filtering**: Apply q-value thresholds to retain high-confidence PSMs
10. **Q-value Recalculation**: Re-compute experiment-wide q-values after filtering

Downstream protein handling now runs in later stages:
- `ProteinInferenceSearch` annotates integrated passing precursor tables with inferred protein groups
- `ProteinScoringSearch` builds protein-group tables, fits the protein probit model, and writes protein scores back

### PSM Scoring (score_psms.jl)

**Unified Entry Point**: `score_precursor_isotope_traces` automatically chooses processing strategy:
- **In-Memory**: Datasets within the derived in-memory row budget
- **Out-of-Memory**: Large datasets use sampling and streaming approaches

**Model Selection Strategy**:
- **Adaptive Selection**: For 1,000-100,000 PSMs with model comparison enabled
- **Default LightGBM**: For other dataset sizes or when comparison disabled
- **Cross-validation**: Maintains CV fold consistency for all downstream analysis

**Features**:
- Spectral matching scores (spectral contrast, residuals, GOF)
- Retention time errors and predictions
- Mass accuracy metrics
- MS1 isotope pattern features
- Match-between-runs features (when enabled)

### Model Comparison Framework (score_psms.jl and model_config.jl)

**Adaptive Model Selection**: For datasets with 1,000-100,000 PSMs, PrecursorScoringSearch automatically selects the best-performing model through comparison testing with suppressed output for clean results.

**Available Models** (defined in model_config.jl):
1. **SimpleLightGBM** (Default for small datasets)
   - Reduced feature set (40+ features)
   - Conservative hyperparameters for robust performance
   - Suitable for datasets with limited PSMs

2. **AdvancedLightGBM** (Default for >100K PSMs)
   - Full ADVANCED_FEATURE_SET (50 features including all available metrics)
   - Aggressive hyperparameters for maximum performance
   - Deeper trees (max_depth=10) and lower learning rate (eta=0.05)

3. **ProbitRegression** 
   - Linear probit model using reduced feature set
   - Fast training with interpretable coefficients
   - Good baseline, especially for smaller datasets

4. **SuperSimplified**
   - Minimal feature set (5 core features: spectral contrast, residuals, error norms, intensity explained)
   - LightGBM with conservative parameters
   - Fastest training, minimal overfitting risk

**Selection Process** (in select_psm_scoring_model):
1. **Model Training**: Each model trained on full dataset with output suppression
2. **Performance Evaluation**: Count targets passing user-defined q-value threshold
3. **Model Selection**: Model with highest target count selected
4. **Final Training**: Selected model retrained with progress bars visible

**Feature Sets** (defined in model_config.jl):
- **ADVANCED_FEATURE_SET**: 50 features including all spectral, RT, MS1, and quality metrics
- **REDUCED_FEATURE_SET**: 40+ core features for balanced performance
- **MINIMAL_FEATURE_SET**: 5 essential spectral matching features
- **Cross-run transfer features**: Currently disabled in production runs

**Clean Output Design**:
- Progress bars suppressed during comparison using `show_progress=false` parameter
- stdout redirected to devnull to suppress LightGBM output
- Only essential results shown: model name and target count at q-value threshold
- Final training shows normal progress bars for user feedback

**Integration**:
- Automatic for datasets with 1,000-100,000 PSMs
- Uses user-defined q_value_threshold from parameters (not hardcoded 0.01)
- Falls back to AdvancedLightGBM for >100K PSMs
- Falls back to SimpleLightGBM if all models fail during comparison

**Performance Evaluation**:
- **Primary Metric**: Number of target PSMs passing q-value threshold
- **Clean Reporting**: Shows results in format "ModelName: X IDs at q ≤ threshold"
- **Selection Display**: Clear indication of selected model with checkmark

### Downstream Protein Stages

- `ProteinInferenceSearch` handles type-safe protein inference and quantification eligibility flags
- `ProteinScoringSearch` handles protein-group feature extraction, semi-supervised probit fitting, and protein-level q-values
- `scoring_interface.jl` provides file reference operations used by precursor rescoring only

## Data Structures

### Protein Group Key
```julia
(protein_name::String, target::Bool, entrap_id::UInt8)
```

### Protein Group Features
- `pg_score`: Protein group score (log-sum or probit)
- `n_peptides`: Number of unique peptides
- `peptide_coverage`: Fraction of possible peptides observed
- `n_possible_peptides`: Total peptides in library
- `global_pg_score`: Log-odds combination of scores across all files

## Memory Management

- In-memory vs out-of-memory precursor scoring is selected from the derived row budget
- File reference abstractions avoid unnecessary table copies
- ML training dominates memory usage in this stage

## Configuration Parameters

### Core Scoring Parameters
```json
"optimization": {
    "machine_learning": {
        "max_in_memory_table_mb": 2000,
        "min_trace_prob": 0.01,
        "spline_points": 100,
        "q_value_interpolation_points_per_bin": 100,
        "min_PEP_neg_threshold_itr": 0.90
    }
},
"proteinScoring": {
    "min_peptides": 2,
    "write_qc_plots": false,
    "log_feature_importance": false
},
"global": {
    "scoring": {
        "q_value_threshold": 0.01  // Used for model comparison target counting
    }
}
```

### Model Comparison Behavior
- **Automatic Selection**: Enabled by default for datasets with 1,000-100,000 PSMs
- **Dataset Size Thresholds**:
  - < 1,000 PSMs: Uses SimpleLightGBM directly (no comparison)
  - 1,000-100,000 PSMs: Automatic model comparison and selection
  - > 100,000 PSMs (in-memory): Uses AdvancedLightGBM directly
  - Above the derived in-memory row budget from `optimization.machine_learning.max_in_memory_table_mb`: Out-of-memory processing with default LightGBM
- **Q-value Threshold**: Uses the user-defined `q_value_threshold` from `global.scoring`
- **Clean Output**: Progress bars and verbose output suppressed during comparison

### Downstream Protein Handling
Protein handling now runs in separate post-integration steps:
```json
"proteinScoring": {
    "min_peptides": 2,
    "write_qc_plots": true,
    "log_feature_importance": true
}
```

- `ProteinInferenceSearch` annotates integrated passing precursors with inferred protein groups
- `ProteinScoringSearch/protein_inference_pipeline.jl` builds protein-group tables from those annotated precursors
- The model fitting, QC plotting, and protein-level q-value logic live under `ProteinScoringSearch/`.

## Common Issues and Solutions

### Model Comparison Issues

#### Issue: Model comparison not triggering
**Causes**: 
- PSM count outside 1,000-100,000 range
- `enable_model_comparison: false` in parameters  
- Using out-of-memory processing (>100K PSMs)
**Solution**: Verify PSM count and parameter settings. Model comparison only works with in-memory processing.

#### Issue: Model training failures during comparison
**Symptoms**: "No models trained successfully - falling back to existing logic"
**Solution**: Check feature availability in PSM data. Missing features cause training failures.

#### Issue: Poor model validation performance
**Symptoms**: Very low number of targets passing q-value threshold
**Solution**: Inspect validation split quality and feature distributions. Ensure balanced target/decoy ratios.

### Pipeline Issues

#### Issue: Q-value calculation errors
**Symptoms**: NaN or infinite q-values
**Solution**: Verify target/decoy balance and score distributions. Check FDR scale factors.

## Testing

### Integration Testing
```julia
SearchDIA("./data/ecoli_test/ecoli_test_params.json")
```

### Model Comparison Testing
```julia
# Test with model comparison enabled
params = GetSearchParams("lib.arrow", "data_files.json", "output/")
params.optimization.machine_learning.enable_model_comparison = true
SearchDIA(params)
```

### Unit Testing Individual Components
```julia
# Test model comparison framework
include("test/UnitTests/PrecursorScoringSearch/test_model_comparison.jl")
```

### Key Outputs to Verify

**Precursor/PSM Files**:
- `trace_prob`, `prec_prob`, `global_prob` - Model predictions and aggregated probabilities
- `qval`, `global_qval`, `pep` - FDR estimates at different levels
- `passing_psms` output tables are the input to later integration and protein stages

**Model Comparison Outputs** (when enabled):
- `model_comparison_report.csv` - Detailed performance metrics
- Log messages showing model selection process
- Training time and feature count reporting

## Performance Considerations

- File I/O dominates for large precursor rescoring runs
- In-memory training scales with precursor-table size and feature count
- Out-of-memory processing trades memory for additional sampling and streaming overhead
- Thread-safe operations in `score_psms.jl`

## Recent Changes

### Trait-Based ML Scoring System (February 2025)

The PSM scoring system was refactored to use a trait-based architecture for better composability and extensibility. The core implementation lives in `src/utils/ML/` with PrecursorScoringSearch using it via `sort_of_percolator!()`.

**New Architecture**:
```julia
# 6 independent traits combined into ScoringConfig
ScoringConfig{M,P,T,F,I,B}
├── model::PSMScoringModel         # LightGBMScorer, ProbitScorer
├── pairing::PairingStrategy       # RandomPairing, NoPairing
├── training_data::TrainingDataStrategy  # QValueNegativeMining, AllDataSelection
├── feature_selection::FeatureSelectionStrategy  # IterativeFeatureSelection
├── iteration_scheme::IterationScheme  # FixedIterationScheme
└── mbr_update::MBRUpdateStrategy  # PairBasedMBR, NoMBR
```

**Key Files**:
- `src/utils/ML/scoring_traits.jl` - Trait type definitions
- `src/utils/ML/scoring_config.jl` - ScoringConfig combining traits
- `src/utils/ML/percolator_generic.jl` - Main `percolator_scoring!()` entry point
- `src/utils/ML/psm_container.jl` - AbstractPSMContainer data abstraction
- `src/utils/ML/pairing.jl` - Pairing strategy implementations

**Integration with PrecursorScoringSearch**:
- `sort_of_percolator!()` in `percolatorSortOf.jl` builds a `ScoringConfig` from parameters
- Delegates to `percolator_scoring!(psms, config)` for actual scoring
- Model comparison in `score_psms.jl` uses `build_scoring_config()` from `model_config.jl`

**See**: `src/utils/ML/CLAUDE.md` for detailed trait-based architecture documentation.

### Automatic Model Selection Feature (January 2025)
- **Four-Model Comparison**: Added automatic comparison of SimpleLightGBM, AdvancedLightGBM, ProbitRegression, and SuperSimplified models
- **Advanced Feature Set**: Added ADVANCED_FEATURE_SET with 50 features for maximum performance model
- **User-Defined Q-value**: Model selection now uses q_value_threshold from parameters instead of hardcoded 0.01
- **Clean Output Design**: 
  - Implemented stdout redirection to suppress LightGBM progress bars
  - Added `show_progress` parameter to control ProgressBar display
  - Removed `colsample_bynode` parameter (not supported by LightGBM)
  - Fixed duplicate parameter issues in train_booster calls
- **Separation of Concerns**: Separated PSM scoring from file I/O operations
  - Removed file writing from `sort_of_percolator_in_memory!`
  - Added dedicated `write_scored_psms_to_files!` function
  - Eliminated need for file backup/restore during model comparison
- **Model Configuration**: Created `model_config.jl` with ModelConfig struct and feature set definitions
- **Deprecated Code**: Marked `model_comparison.jl` as deprecated (functionality moved to score_psms.jl)

### Pipeline Refinement
- **Comprehensive Documentation**: Detailed breakdown of precursor rescoring and FDR-control flow
- **Phase-Based Organization**: Organized steps into logical phases for scoring and precursor filtering
- **Memory-Efficient Processing**: Enhanced file reference system with pipeline operations
- **Q-value Recalculation**: Added step 10 for post-filtering q-value updates

### Enhanced Pipeline Operations
- **TransformPipeline**: Composable file transformation operations with automatic sort state tracking
- **Stream Merging**: Generic N-key heap-based merging supporting arbitrary sort keys
- **File Reference System**: Type-safe PSMFileReference and ProteinGroupFileReference abstractions
- **Column Management**: Automatic addition/removal of temporary columns during processing

### Protein Inference Modernization
- **Type-Safe Inference**: Protein inference now lives in `ProteinInferenceSearch`
- **Composable Pipeline**: Protein-group table building now lives under `ProteinScoringSearch/protein_inference_pipeline.jl`
- **Stage Separation**: Protein inference and protein scoring now run after chromatogram integration
