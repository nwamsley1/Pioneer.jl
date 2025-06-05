# Protein Group Scoring Enhancement Plan

## Current State

Pioneer.jl currently scores protein groups using a simple log-sum approach:
- **Formula**: `pg_score = -sum(log(1 - peptide_prob))` for all peptides
- **Location**: `getProteinGroupsDict()` in `utils.jl:640-730`
- **Limitation**: Only considers individual peptide probabilities, ignoring other valuable features

## Proposed Enhancement: Multi-Feature Protein Scoring

### Objective
Develop a more sophisticated protein group scoring system that leverages multiple features from child peptides to improve protein-level FDR control and quantification accuracy.

## Implementation Plan

### Phase 1: Feature Engineering

#### 1.1 Peptide-Level Features to Extract
- **Count Features**:
  - `n_peptides`: Total number of unique peptides
  - `n_psms`: Total number of PSMs supporting the protein
  - `n_shared_peptides`: Number of peptides shared with other proteins
  - `n_unique_peptides`: Number of protein-specific peptides

- **Score Distribution Features**:
  - `top_peptide_score`: Highest peptide probability
  - `median_peptide_score`: Median of all peptide probabilities
  - `mean_peptide_score`: Mean peptide probability
  - `std_peptide_score`: Standard deviation of peptide scores
  - `min_peptide_score`: Lowest peptide probability

- **Coverage Features**:
  - `sequence_coverage`: Percentage of protein sequence covered
  - `spectral_count`: Total MS2 spectra matched
  - `avg_spectral_count_per_peptide`: Spectral count normalized by peptide count

- **Quality Features**:
  - `avg_mass_error`: Average mass error of peptides
  - `avg_rt_deviation`: Average RT deviation from predicted
  - `charge_state_diversity`: Number of different charge states observed

### Phase 2: Scoring Models

#### 2.1 Traditional Enhancement (Non-ML)
Implement weighted log-sum scoring:
```julia
pg_score = α * log_sum_score + 
           β * log(n_unique_peptides + 1) + 
           γ * log(n_psms + 1) +
           δ * (sequence_coverage / 100)
```

#### 2.2 Machine Learning Integration
Leverage existing XGBoost infrastructure:

1. **Training Data Collection**:
   - Use protein groups from FDR-controlled ScoringSearch
   - Extract features for target/decoy proteins
   - Create balanced training sets per CV fold

2. **Model Architecture**:
   - Reuse `percolatorSortOf.jl` patterns
   - Random forest with 100-500 trees
   - Cross-validation consistent with PSM scoring

3. **Feature Vector**:
   ```julia
   features = [
       n_peptides,
       n_unique_peptides,
       n_shared_peptides,
       top_peptide_score,
       median_peptide_score,
       mean_peptide_score,
       std_peptide_score,
       sequence_coverage,
       spectral_count,
       avg_mass_error,
       avg_rt_deviation,
       charge_state_diversity,
       log_sum_score  # Include original score as feature
   ]
   ```

### Phase 3: Implementation Details

#### 3.1 Code Structure
```julia
# New file: src/Routines/SearchDIA/SearchMethods/ScoringSearch/proteinGroupML.jl

struct ProteinGroupFeatures
    n_peptides::Int32
    n_unique_peptides::Int32
    n_shared_peptides::Int32
    top_peptide_score::Float32
    median_peptide_score::Float32
    mean_peptide_score::Float32
    std_peptide_score::Float32
    sequence_coverage::Float32
    spectral_count::Int32
    avg_mass_error::Float32
    avg_rt_deviation::Float32
    charge_state_diversity::Int32
    log_sum_score::Float32
end

function extractProteinGroupFeatures(
    protein_groups::Dictionary,
    psm_data::NamedTuple,
    precursors::LibraryPrecursors,
    protein_sequences::Dict{String, String}
) -> Dict{NamedTuple, ProteinGroupFeatures}
```

#### 3.2 Integration Points

1. **Modify `getProteinGroupsDict()`**:
   - Collect additional statistics during aggregation
   - Store peptide-level details for feature extraction

2. **Add to `ScoringSearch.jl`**:
   - New parameter: `use_ml_protein_scoring::Bool`
   - Call ML scoring after traditional scoring
   - Maintain backward compatibility

3. **Update `writeProteinGroups()`**:
   - Add `ml_pg_score` column when ML is enabled
   - Include feature columns for interpretability

#### 3.3 Memory-Efficient Processing
Following existing patterns:
- Process protein groups in chunks
- Use Arrow file streaming
- Leverage thread-local storage for feature calculation

### Phase 4: Validation & Testing

#### 4.1 Unit Tests
```julia
# test/UnitTests/proteinGroupML.jl
- Test feature extraction accuracy
- Validate score improvements
- Check memory usage
```

#### 4.2 Integration Tests
- Compare traditional vs ML scoring on E. coli dataset
- Validate protein FDR control
- Benchmark performance impact

#### 4.3 Quality Metrics
- Protein-level FDR curves
- Number of proteins at 1% FDR
- Correlation with known protein abundances
- Reproducibility across technical replicates

### Phase 5: Configuration

Add to search parameters JSON:
```json
{
  "protein_scoring": {
    "method": "ml",  // "traditional", "weighted", "ml"
    "min_peptides": 2,
    "ml_config": {
      "n_trees": 100,
      "max_depth": 6,
      "learning_rate": 0.1,
      "feature_subset": ["all"]  // or specific feature names
    },
    "weighted_config": {
      "alpha": 1.0,   // log-sum weight
      "beta": 0.5,    // unique peptide weight
      "gamma": 0.3,   // PSM count weight
      "delta": 0.2    // coverage weight
    }
  }
}
```

## Benefits

1. **Improved Sensitivity**: Better distinguish true proteins from false positives
2. **Robustness**: Less susceptible to single high-scoring contaminant peptides
3. **Interpretability**: Feature importance analysis reveals decision factors
4. **Flexibility**: Easy to add domain-specific features

## Risks & Mitigation

1. **Performance Impact**: 
   - Mitigation: Parallelize feature extraction, cache results

2. **Overfitting**: 
   - Mitigation: Cross-validation, feature selection, regularization

3. **Backward Compatibility**: 
   - Mitigation: Feature flag, default to current method

## Future Extensions

1. **Deep Learning**: Graph neural networks for protein complexes
2. **Bayesian Models**: Uncertainty quantification in scores
3. **Transfer Learning**: Pre-trained models across organisms
4. **Interactive Features**: Peptide-peptide interaction features

## References

- Cox, J., & Mann, M. (2008). MaxQuant enables high peptide identification rates
- Nesvizhskii, A. I., & Aebersold, R. (2005). Interpretation of shotgun proteomic data
- The, M., et al. (2016). Fast and accurate protein false discovery rates on large-scale proteomics data sets