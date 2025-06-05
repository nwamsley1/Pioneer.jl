# Probit Regression Implementation for Protein Group Scoring

## Overview
This document describes the implementation of probit regression for protein group scoring in Pioneer.jl, replacing the previous Linear Discriminant Analysis (LDA) approach.

## Implementation Details

### Location
The implementation is in `/src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl` within the `get_protein_groups` function (lines 1035-1224).

### Features Used
The probit regression model uses four features:
1. **pg_score** - The protein group score (sum of log(1-p) for constituent peptides)
2. **peptide_coverage** - Number of matched peptides / number of possible peptides  
3. **n_possible_peptides** - Total number of possible peptides for the protein
4. **log_n_possible_peptides** - Log-transformed version of n_possible_peptides

### Key Implementation Steps

1. **Feature Extraction and Standardization** (lines 1013-1033)
   ```julia
   feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_n_possible_peptides]
   X = Matrix{Float64}(all_protein_groups[:, feature_names])
   y = all_protein_groups.target
   
   # Standardize features
   X_mean = mean(X, dims=1)
   X_std = std(X, dims=1)
   X_std[X_std .== 0] .= 1.0  # Avoid division by zero
   X_standardized = (X .- X_mean) ./ X_std
   ```

2. **Probit Model Fitting** (lines 1038-1056)
   ```julia
   # Add intercept column
   X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
   X_df = DataFrame(X_with_intercept, [:intercept; feature_names])
   
   # Initialize coefficients
   β = zeros(Float64, size(X_with_intercept, 2))
   
   # Create data chunks for parallel processing
   n_chunks = max(1, Threads.nthreads())
   chunk_size = max(1, ceil(Int, length(y) / n_chunks))
   data_chunks = Iterators.partition(1:length(y), chunk_size)
   
   # Fit probit model
   β_fitted = Pioneer.ProbitRegression(β, X_df, y, data_chunks, max_iter=30)
   
   # Get predicted probabilities
   prob_scores = zeros(Float64, length(y))
   Pioneer.ModelPredictProbs!(prob_scores, X_df, β_fitted, data_chunks)
   ```

3. **Performance Evaluation** (lines 1058-1096)
   - Classification accuracy
   - AUC (Area Under Curve)
   - Feature coefficients
   - Score separation statistics

4. **Q-Value Calculation** (lines 1098-1166)
   - Sorts protein groups by probit probability scores
   - Calculates FDR and converts to q-values with monotonicity enforcement
   - Reports targets passing at 1%, 5%, and 10% FDR thresholds

5. **Comparison with pg_score-only baseline** (lines 1168-1223)
   - Calculates q-values using only pg_score for comparison
   - Reports improvement metrics at different FDR thresholds

### Handling Protein Groups with Multiple Proteins

The implementation properly handles protein groups containing multiple proteins separated by semicolons (lines 795-819):
```julia
# Handle protein groups with multiple proteins separated by semicolons
n_possible_peptides = zeros(Int64, length(keys_array))
for (i, k) in enumerate(keys_array)
    # Split the protein group name by semicolons
    protein_names_in_group = split(k[:protein_name], ';')
    
    # Union of all peptide sets from proteins in the group
    all_possible_peptides = Set{String}()
    for individual_protein in protein_names_in_group
        # Create key for each individual protein
        individual_key = (protein_name = String(individual_protein), 
                        target = k[:target], 
                        entrap_id = k[:entrap_id])
        # Get the set of peptides for this protein and union with existing
        if haskey(protein_to_possible_peptides, individual_key)
            union!(all_possible_peptides, protein_to_possible_peptides[individual_key])
        end
    end
    
    # Count unique peptides across all proteins in the group
    n_possible_peptides[i] = max(length(all_possible_peptides), 1)
end
```

### Advantages over LDA
1. **Binary Classification Focus**: Probit regression is specifically designed for binary outcomes (target/decoy)
2. **Probability Interpretation**: Direct probability outputs are more interpretable than LDA discriminant scores
3. **Robustness**: Less sensitive to multicollinearity among features
4. **Parallel Processing**: Utilizes Julia's threading capabilities for efficient computation

### Output and Logging
The implementation provides comprehensive logging including:
- Feature statistics and correlations
- Model coefficients
- Classification performance metrics
- Q-value analysis at multiple FDR thresholds
- Comparison with baseline approach
- Scatter plots of key features (saved to desktop)

## Usage
The probit regression is automatically used when `get_protein_groups` is called during the ScoringSearch phase of the SearchDIA pipeline. No additional configuration is required beyond the standard Pioneer.jl parameters.

## Dependencies
- Uses the existing `probitRegression.jl` implementation in `src/utils/ML/`
- Requires SpecialFunctions.jl for error function calculations
- Utilizes DataFrames.jl for data manipulation
- Uses Plots.jl for visualization