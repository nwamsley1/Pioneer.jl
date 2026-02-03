# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the ScoringSearch module.

## ScoringSearch Overview

ScoringSearch is the 7th stage in the Pioneer DIA search pipeline. It performs machine learning-based PSM rescoring, FDR control, and protein group analysis with optional ML-enhanced protein scoring. The module features adaptive model selection for datasets with 1,000-100,000 PSMs, automatically choosing between LightGBM, probit regression, and simplified models based on validation performance.

## Module Structure

```
ScoringSearch/
├── ScoringSearch.jl              # Main search method implementation with 23-step pipeline
├── score_psms.jl                 # PSM scoring entry point with automatic in-memory/OOM selection and model comparison
├── model_config.jl               # Model configurations for comparison (SimpleLightGBM, AdvancedLightGBM, ProbitRegression, SuperSimplified)
├── model_comparison.jl           # DEPRECATED - Model comparison framework (functionality moved to score_psms.jl)
├── utils.jl                      # Protein group analysis and helper functions (legacy)
├── utils_protein_ml.jl           # ML-enhanced protein scoring (when enabled)
├── protein_inference_pipeline.jl # Modern composable protein inference pipeline
└── scoring_interface.jl          # Type-safe file reference operations
```

## Key Components

### Main Workflow (ScoringSearch.jl)

ScoringSearch implements a comprehensive 23-step pipeline:

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

**Phase 3: Protein Inference & Scoring (Steps 11-23)**
11. **Protein Peptide Counting**: Count possible peptides per protein for feature calculation
12. **Protein Inference**: Group peptides into minimal protein sets using parsimony
13. **CV Fold Mapping**: Build protein-to-CV-fold mapping for consistent scoring
14. **Probit Regression**: Refine protein scores with additional features
15. **Global Protein Scores**: Calculate max scores across all runs per protein
16-22. **Protein Q-value Processing**: Sort, merge, and calculate protein-level FDR
23. **PSM Score Updates**: Backpropagate final protein scores to PSMs

### PSM Scoring (score_psms.jl)

**Unified Entry Point**: `score_precursor_isotope_traces` automatically chooses processing strategy:
- **In-Memory**: Datasets ≤ `max_psms_in_memory` (typically 100,000 PSMs)
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

**Adaptive Model Selection**: For datasets with 1,000-100,000 PSMs, ScoringSearch automatically selects the best-performing model through comparison testing with suppressed output for clean results.

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
- **MBR Features**: Automatically appended when match_between_runs=true

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

### Protein Group Analysis

**Modern Pipeline (protein_inference_pipeline.jl)**:
- **`perform_protein_inference_pipeline`**: Composable pipeline approach
- **`apply_inference_to_dataframe`**: Wrapper around type-safe `infer_proteins()`
- **`group_psms_by_protein`**: Aggregates PSMs into protein groups
- **`add_inferred_protein_column`**: Updates PSMs with protein assignments
- **`add_quantification_flag`**: Marks peptides for quantification use

**Legacy Functions (utils.jl)**:
- **`get_protein_groups`**: Main entry point for protein analysis (legacy)
- **`perform_protein_inference`**: File-by-file protein inference (legacy)
- **`getProteinGroupsDict`**: Creates protein groups from PSMs (legacy)
- **`writeProteinGroups`**: Outputs protein groups with features
- **`update_psms_with_probit_scores`**: Updates PSMs with refined scores
- **`merge_sorted_protein_groups`**: Memory-efficient merging

**Type-Safe Operations (scoring_interface.jl)**:
- File reference-based operations for safer data handling
- Abstracts file operations from core algorithms

### Protein Scoring Flow

1. **Initial Scoring**: Log-sum of peptide probabilities
2. **Probit Regression**: Uses features like peptide coverage, n_peptides
3. **Global Scoring**: Max score across all runs
4. **PSM Update**: Backpropagates refined scores to PSMs

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

- File-by-file protein inference reduces memory usage
- Heap-based merging for large datasets
- Bidirectional PSM-PG file mappings maintained

## Configuration Parameters

### Core Scoring Parameters
```json
"optimization": {
    "machine_learning": {
        "max_psms_in_memory": 100000,
        "min_trace_prob": 0.01,
        "spline_points": 100,
        "interpolation_points": 100
    }
},
"protein_inference": {
    "min_peptides": 2
},
"global_settings": {
    "scoring": {
        "q_value_threshold": 0.01  // Used for model comparison target counting
    },
    "match_between_runs": true
}
```

### Model Comparison Behavior
- **Automatic Selection**: Enabled by default for datasets with 1,000-100,000 PSMs
- **Dataset Size Thresholds**:
  - < 1,000 PSMs: Uses SimpleLightGBM directly (no comparison)
  - 1,000-100,000 PSMs: Automatic model comparison and selection
  - > 100,000 PSMs (in-memory): Uses AdvancedLightGBM directly
  - > max_psms_in_memory: Out-of-memory processing with default LightGBM
- **Q-value Threshold**: Uses the user-defined `q_value_threshold` from global_settings.scoring
- **Clean Output**: Progress bars and verbose output suppressed during comparison

### ML Protein Scoring (Optional)
When enabled via parameters:
```json
"machine_learning": {
    "use_ml_protein_scoring": true,
    "n_top_precursors": 5
}
```

Uses LightGBM with top N precursor scores as features, implemented in `utils_protein_ml.jl`.

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

### Legacy Issues

#### Issue: BoundsError in perform_protein_inference
**Solution**: Function must return 3 values: `(pg_count, psm_to_pg_path, pg_to_psm_path)`

#### Issue: Missing probit scores in PSMs
**Solution**: Use `update_psms_with_probit_scores` after probit regression

#### Issue: Memory overflow with many files
**Solution**: Out-of-memory probit regression automatically triggered

### Pipeline Issues

#### Issue: Empty protein groups after inference
**Symptoms**: "No protein groups created during protein inference"
**Solution**: Check min_peptides parameter and PSM quality. Lower threshold or investigate upstream filtering.

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
include("test/UnitTests/ScoringSearch/test_model_comparison.jl")

# Test protein inference pipeline
include("test/UnitTests/ScoringSearch/test_protein_inference.jl")
```

### Key Outputs to Verify

**PSM Files**:
- `trace_prob`, `prec_prob`, `global_prob` - Model predictions and aggregated probabilities
- `qval`, `global_qval`, `pep` - FDR estimates at different levels
- `pg_score`, `global_pg_score` - Protein group scores propagated back
- `inferred_protein_group` - Protein inference assignments

**Protein Group Files**:
- `pg_score`, `global_pg_score` - Initial and refined protein scores
- `pg_qval`, `global_pg_qval`, `pg_pep` - Protein-level FDR estimates
- `n_peptides`, `peptide_coverage` - Protein inference features
- `n_possible_peptides` - Library-based peptide counts

**Model Comparison Outputs** (when enabled):
- `model_comparison_report.csv` - Detailed performance metrics
- Log messages showing model selection process
- Training time and feature count reporting

## Performance Considerations

- Protein inference is O(n_peptides * n_proteins)
- Probit regression scales with n_protein_groups
- File I/O dominates for large experiments
- Thread-safe operations in score_psms.jl

## Recent Changes (2025-02)

### OOM Percolator Pipeline Refinement (February 2025)
- **Per-Fold Q-Value Computation**: Fixed `build_mbr_tracker_for_fold` to compute q-values on CV fold subsets (matching in-memory behavior)
- **MBR_num_runs Consistency**: Both modes now exclude current run from count if it has passing PSMs
- **Transfer Candidate Logic**: Verified equivalent implementation between in-memory and OOM modes
- **Comprehensive Documentation**: Added `/docs/plans/in-memory-vs-oom-percolator-documentation.md`
- **Result Alignment**: OOM and in-memory modes now produce ~0.3% difference (within expected variance)

### Key OOM Functions
| Function | Location | Purpose |
|----------|----------|---------|
| `apply_models_to_files!` | `percolatorSortOf.jl:857` | Apply trained models file-by-file |
| `build_mbr_tracker_for_fold` | `percolatorSortOf.jl:960` | Build per-fold MBR tracker |
| `apply_mbr_from_tracker!` | `percolatorSortOf.jl:1099` | Apply MBR features from tracker |

## Recent Changes (2025-01)

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

### 23-Step Pipeline Refinement
- **Comprehensive Documentation**: Detailed step-by-step breakdown of entire ScoringSearch pipeline
- **Phase-Based Organization**: Organized steps into logical phases (Model Training, PSM Processing, Protein Inference)
- **Memory-Efficient Processing**: Enhanced file reference system with pipeline operations
- **Q-value Recalculation**: Added step 10 for post-filtering q-value updates

### Enhanced Pipeline Operations
- **TransformPipeline**: Composable file transformation operations with automatic sort state tracking
- **Stream Merging**: Generic N-key heap-based merging supporting arbitrary sort keys
- **File Reference System**: Type-safe PSMFileReference and ProteinGroupFileReference abstractions
- **Column Management**: Automatic addition/removal of temporary columns during processing

### Protein Inference Modernization
- **Type-Safe Inference**: Migrated to `infer_proteins()` using `ProteinKey` and `PeptideKey` types
- **Composable Pipeline**: Added `protein_inference_pipeline.jl` with modular operations
- **CV Fold Mapping**: Protein-to-CV-fold mapping for consistent cross-validation
- **Legacy Preservation**: Maintained backward compatibility with existing `utils.jl` functions