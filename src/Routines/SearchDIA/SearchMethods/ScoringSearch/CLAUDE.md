# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the ScoringSearch module.

## ScoringSearch Overview

ScoringSearch is the 7th stage in the Pioneer DIA search pipeline. It performs machine learning-based PSM rescoring, FDR control, and protein group analysis with optional ML-enhanced protein scoring. The module features adaptive model selection for datasets with 1,000-100,000 PSMs, automatically choosing between XGBoost, probit regression, and simplified models based on validation performance.

## Module Structure

```
ScoringSearch/
├── ScoringSearch.jl              # Main search method implementation with 23-step pipeline
├── score_psms.jl                 # PSM scoring entry point with automatic in-memory/OOM selection
├── model_comparison.jl           # Model comparison framework for adaptive model selection
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
- **Default XGBoost**: For other dataset sizes or when comparison disabled
- **Cross-validation**: Maintains CV fold consistency for all downstream analysis

**Features**:
- Spectral matching scores (spectral contrast, residuals, GOF)
- Retention time errors and predictions
- Mass accuracy metrics
- MS1 isotope pattern features
- Match-between-runs features (when enabled)

### Model Comparison Framework (model_comparison.jl)

**Adaptive Model Selection**: For datasets with 1,000-100,000 PSMs, ScoringSearch can automatically select the best-performing model through validation testing.

**Available Models**:
1. **SimpleXGBoost** (Current Default)
   - Full reduced feature set (40+ features)
   - Standard XGBoost hyperparameters
   - 3-fold cross-validation training

2. **ProbitRegression** 
   - Linear probit model with same feature set
   - Fast training, interpretable coefficients
   - Good baseline performance

3. **SuperSimplified**
   - Minimal feature set (5 core features)
   - Same XGBoost architecture as SimpleXGBoost
   - Fastest training, reduced overfitting risk

**Selection Process**:
1. **Data Split**: 80% training, 20% validation with target/decoy stratification
2. **Model Training**: All three models trained on training set
3. **Validation Testing**: Models evaluated on held-out validation data
4. **Performance Metrics**: Primary metric is number of targets passing q-value ≤ 0.01
5. **Model Selection**: Best model retrained on full dataset (100% of PSMs)

**Feature Sets**:
- **Reduced Feature Set**: Core peptide properties, RT features, spectral features, quality metrics, MS1 features
- **Minimal Feature Set**: 5 essential features (spectral contrast, residuals, error norms, intensity explained)
- **MBR Features**: Automatically added when match_between_runs=true

**Integration**:
- Enabled via `enable_model_comparison: true` in parameters
- Only active for in-memory processing (1K-100K PSMs)
- Falls back to default XGBoost for other dataset sizes
- Generates detailed performance reports in CSV format

**Performance Evaluation**:
- **Primary**: Number of target PSMs passing q-value threshold
- **Secondary**: AUC, accuracy, sensitivity, specificity
- **Efficiency**: Training time and feature count tracking

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
        "interpolation_points": 100,
        
        // Model Comparison (NEW)
        "enable_model_comparison": true,
        "validation_split_ratio": 0.2,
        "qvalue_threshold": 0.01,
        "min_psms_for_comparison": 1000,
        "max_psms_for_comparison": 100000
    }
},
"protein_inference": {
    "min_peptides": 2
},
"global_settings": {
    "scoring": {
        "q_value_threshold": 0.01
    },
    "match_between_runs": true
}
```

### Model Comparison Configuration
- **enable_model_comparison**: Enable adaptive model selection (default: true)
- **validation_split_ratio**: Fraction for validation split (default: 0.2)
- **qvalue_threshold**: Q-value threshold for target counting (default: 0.01)
- **min_psms_for_comparison**: Minimum PSMs to enable comparison (default: 1000)
- **max_psms_for_comparison**: Maximum PSMs to enable comparison (default: 100000)

### ML Protein Scoring (Optional)
When enabled via parameters:
```json
"machine_learning": {
    "use_ml_protein_scoring": true,
    "n_top_precursors": 5
}
```

Uses EvoTrees/XGBoost with top N precursor scores as features, implemented in `utils_protein_ml.jl`.

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
- `prob`, `prec_prob`, `global_prob` - Model predictions and aggregated probabilities
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

## Recent Changes (2025-01)

### Model Comparison Framework Implementation
- **Adaptive Model Selection**: Added comparison of 3 models (SimpleXGBoost, ProbitRegression, SuperSimplified) for 1K-100K PSM datasets
- **Validation-Based Selection**: 80/20 train/validation split with performance-based model selection
- **Feature Set Optimization**: Defined reduced (40+ features) and minimal (5 features) feature sets
- **Automatic Integration**: Seamless integration with existing scoring pipeline via `score_precursor_isotope_traces`
- **Performance Reporting**: Detailed CSV reports with training times, AUC, and target yield metrics

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