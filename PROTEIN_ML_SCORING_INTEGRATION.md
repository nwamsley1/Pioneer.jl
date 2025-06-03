# Protein Group ML Scoring Integration Guide

## Overview

This implementation adds machine learning-based protein group scoring using XGBoost random forests. The model uses the top N precursor scores as features along with additional protein-level statistics.

## Key Features

1. **Cross-validation consistency**: Uses the same CV fold assignments as precursor scoring (from `LibraryIon.jl`)
2. **Random Forest model**: Uses XGBoost's random forest implementation for robust scoring
3. **Top-N precursor features**: Uses scores from the best N precursors per protein group
4. **Additional features**: Number of peptides, number of precursors, mean/median/std of precursor scores

## Integration Steps

### 1. Add the new modules to the project

The implementation consists of two new files:
- `src/utils/ML/proteinGroupScoring.jl` - Core ML scoring functionality
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_protein_ml.jl` - Integration wrapper

### 2. Modify ScoringSearch to use ML scoring

In `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`, make these changes:

```julia
# Add to imports at the top
include("utils_protein_ml.jl")

# In the ScoringSearchParameters struct, add:
use_ml_protein_scoring::Bool
n_top_precursors::Int
num_protein_trees::Int

# In the ScoringSearchParameters constructor, add:
Bool(ml_params.use_ml_protein_scoring),  # Default: true
Int64(ml_params.n_top_precursors),       # Default: 5
Int64(ml_params.num_protein_trees),      # Default: 100

# In summarize_results!, replace the protein scoring section (around line 273):
# Old:
protein_inference_dict = get_protein_groups(
    getPassingPsms(getMSData(search_context)),
    getPassingProteins(getMSData(search_context)),
    passing_proteins_folder,
    temp_folder,
    getPrecursors(getSpecLib(search_context)),
    min_peptides = params.min_peptides
)

# New:
protein_inference_dict = get_protein_groups_with_ml(
    getPassingPsms(getMSData(search_context)),
    getPassingProteins(getMSData(search_context)),
    passing_proteins_folder,
    temp_folder,
    getPrecursors(getSpecLib(search_context)),
    min_peptides = params.min_peptides,
    use_ml_scoring = params.use_ml_protein_scoring,
    n_top_precursors = params.n_top_precursors,
    num_parallel_tree = params.num_protein_trees
)
```

### 3. Add parameters to configuration

In your parameter JSON files, add to the machine learning section:

```json
{
  "optimization": {
    "machine_learning": {
      "use_ml_protein_scoring": true,
      "n_top_precursors": 5,
      "num_protein_trees": 100,
      // ... existing parameters
    }
  }
}
```

## How It Works

1. **Feature Extraction**: For each protein group, extracts:
   - Top N precursor scores (sorted by score)
   - Number of unique peptides
   - Number of total precursors
   - Mean, median, and standard deviation of precursor scores

2. **CV Fold Assignment**: Protein groups inherit CV folds from their constituent peptides
   - Uses the most common fold among peptides
   - Consistent with the library precursor CV fold assignment

3. **Model Training**: 
   - Trains separate random forest models for each CV fold
   - Uses out-of-fold data for unbiased scoring
   - Single round of boosting with multiple trees (typical for RF)

4. **Scoring**: 
   - Applies the appropriate model based on protein's CV fold
   - Falls back to traditional scoring if insufficient data

## Benefits

- **Better calibration**: ML model can learn complex relationships between precursor quality and protein confidence
- **Robustness**: Random forests are less prone to overfitting than gradient boosting
- **Consistency**: Uses same CV scheme as PSM scoring for proper error estimation
- **Flexibility**: Can easily add more features or adjust model parameters

## Debugging

The implementation exports protein features to `temp_folder/protein_features.arrow` for analysis. This includes:
- All features used for training
- Original and ML scores
- CV fold assignments

## Future Enhancements

1. Add more protein-level features:
   - Sequence coverage
   - Number of unique peptides vs shared peptides
   - Hydrophobicity/physicochemical properties

2. Consider iterative scoring:
   - Use initial protein scores as features
   - Multiple rounds of training

3. Experiment with different models:
   - Gradient boosting instead of random forest
   - Different hyperparameters