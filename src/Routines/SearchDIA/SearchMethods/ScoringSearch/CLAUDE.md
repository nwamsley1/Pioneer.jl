# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the ScoringSearch module.

## ScoringSearch Overview

ScoringSearch is the 7th stage in the Pioneer DIA search pipeline. It performs XGBoost-based PSM rescoring, FDR control, and protein group analysis with optional ML-enhanced protein scoring.

## Module Structure

```
ScoringSearch/
├── ScoringSearch.jl       # Main search method implementation
├── score_psms.jl          # XGBoost model training and PSM scoring
├── utils.jl               # Protein group analysis and helper functions
└── utils_protein_ml.jl    # ML-enhanced protein scoring (when enabled)
```

## Key Components

### Main Workflow (ScoringSearch.jl)

The scoring search performs these steps:

1. **XGBoost Model Training**: Trains CV models on second pass PSMs
2. **PSM Rescoring**: Applies models to score all PSMs
3. **Best Trace Selection**: Identifies optimal isotope traces
4. **Q-value Calculation**: Computes FDR at PSM and protein levels
5. **Protein Inference**: Groups peptides into proteins
6. **Probit Regression**: Refines protein scores with additional features
7. **Global Scoring**: Calculates experiment-wide protein scores

### PSM Scoring (score_psms.jl)

- Implements cross-validation XGBoost training
- Features include spectral scores, RT errors, mass errors
- Handles both in-memory and out-of-memory training
- Maintains CV fold consistency for downstream analysis

### Protein Group Analysis (utils.jl)

Key functions:

- **`get_protein_groups`**: Main entry point for protein analysis
- **`perform_protein_inference`**: File-by-file protein inference
- **`getProteinGroupsDict`**: Creates protein groups from PSMs
- **`writeProteinGroups`**: Outputs protein groups with features
- **`update_psms_with_probit_scores`**: Updates PSMs with refined scores
- **`merge_sorted_protein_groups`**: Memory-efficient merging

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
- `global_pg_score`: Max score across all files

## Memory Management

- File-by-file protein inference reduces memory usage
- Heap-based merging for large datasets
- Bidirectional PSM-PG file mappings maintained

## ML Protein Scoring (Optional)

When enabled via parameters:
```json
"machine_learning": {
    "use_ml_protein_scoring": true,
    "n_top_precursors": 5
}
```

Uses XGBoost with top N precursor scores as features, implemented in `utils_protein_ml.jl`.

## Common Issues and Solutions

### Issue: BoundsError in perform_protein_inference
**Solution**: Function must return 3 values: `(pg_count, psm_to_pg_path, pg_to_psm_path)`

### Issue: Missing probit scores in PSMs
**Solution**: Use `update_psms_with_probit_scores` after probit regression

### Issue: Memory overflow with many files
**Solution**: Out-of-memory probit regression automatically triggered

## Testing

Integration test via:
```julia
SearchDIA("./data/ecoli_test/ecoli_test_params.json")
```

Key outputs to verify:
- PSM files have `pg_score`, `global_pg_score`, `inferred_protein_group`
- Protein group files have features and scores
- Q-value calculations produce reasonable FDR estimates

## Performance Considerations

- Protein inference is O(n_peptides * n_proteins)
- Probit regression scales with n_protein_groups
- File I/O dominates for large experiments
- Thread-safe operations in score_psms.jl

## Recent Changes

- Moved protein inference functions to module scope for accessibility
- Added `update_psms_with_probit_scores` for proper score propagation
- Enhanced bidirectional file mapping between PSMs and protein groups