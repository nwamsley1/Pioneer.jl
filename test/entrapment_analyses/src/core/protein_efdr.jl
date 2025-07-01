"""
Protein-specific EFDR calculation functions.

These functions are designed to work directly with protein data without requiring
precursor indices or library DataFrames.
"""

"""
    add_protein_efdr_columns!(protein_results::DataFrame;
                             method_types::Vector=[CombinedEFDR, PairedEFDR],
                             score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:pg_score, :qval)],
                             r::Float64=1.0)

Add empirical FDR columns to protein results based on entrapment analysis.

This function is specifically designed for protein-level data where:
- entrap_id is directly available in the DataFrame
- No precursor_idx mapping is needed
- Original target scores have already been added

# Arguments
- `protein_results::DataFrame`: Protein results with columns:
  - `entrap_id`: Entrapment group ID (0 = original)
  - Score columns specified in `score_qval_pairs`
  - Original target score columns (e.g., `pg_score_original_target`)
- `method_types`: EFDR calculation methods to use
- `score_qval_pairs`: Vector of (score_col, qval_col) tuples
- `r`: Ratio parameter for EFDR calculation

# Notes
- Adds columns named as: {score_col}_{method_name}_efdr
- Example: pg_score_combined_efdr, pg_score_paired_efdr
"""
function add_protein_efdr_columns!(protein_results::DataFrame;
                                  method_types::Vector=[CombinedEFDR, PairedEFDR],
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:pg_score, :qval)],
                                  r::Float64=1.0)
    
    # Check for required columns
    if !hasproperty(protein_results, :entrap_id)
        error("DataFrame must have :entrap_id column")
    end
    
    # Extract entrapment labels directly
    entrap_labels = protein_results.entrap_id
    
    # Process each score/qval pair
    for (score_col, qval_col) in score_qval_pairs
        # Check columns exist
        if !hasproperty(protein_results, score_col)
            @warn "Column $score_col not found in DataFrame, skipping..."
            continue
        end
        if !hasproperty(protein_results, qval_col)
            @warn "Column $qval_col not found in DataFrame, skipping..."
            continue
        end
        
        # Get original target score column name
        original_target_col = Symbol(String(score_col) * "_original_target")
        if !hasproperty(protein_results, original_target_col)
            @warn "Column $original_target_col not found. Make sure to run add_original_target_protein_scores! first."
            continue
        end
        
        # Extract vectors
        scores = Float64.(protein_results[!, score_col])
        original_target_scores = Float64.(protein_results[!, original_target_col])
        qvals = Float64.(protein_results[!, qval_col])
        
        # Apply each method type
        for method_type in method_types
            # Determine method name for column
            method_name = if method_type == CombinedEFDR
                "combined"
            elseif method_type == PairedEFDR
                "paired"
            else
                lowercase(string(method_type))
            end
            
            # Create EFDR column name
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            
            # Create method instance and calculate EFDR
            method = method_type(scores, original_target_scores, entrap_labels, qvals; r=r)
            efdr_values = calculate_efdr(method)
            
            # Add to dataframe
            protein_results[!, efdr_col] = efdr_values
        end
    end
    
    return nothing
end

"""
    compare_protein_efdr_methods(protein_results::DataFrame, qval_col::Symbol, score_col::Symbol)

Compare different EFDR methods for protein data.

Similar to `compare_efdr_methods` but works directly with protein data without
requiring precursor indices or library DataFrames.

# Arguments
- `protein_results::DataFrame`: Protein results with entrap_id and EFDR columns
- `qval_col::Symbol`: Q-value column to compare against
- `score_col::Symbol`: Score column used for EFDR calculation

# Returns
- DataFrame with comparison results at different thresholds
"""
function compare_protein_efdr_methods(protein_results::DataFrame, qval_col::Symbol, score_col::Symbol)
    # Get method names
    combined_efdr_col = Symbol(String(score_col) * "_combined_efdr")
    paired_efdr_col = Symbol(String(score_col) * "_paired_efdr")
    
    # Check required columns exist
    required_cols = [qval_col, combined_efdr_col, paired_efdr_col, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    
    # Extract entrapment labels
    entrap_labels = protein_results.entrap_id
    
    # Define thresholds to test
    thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    
    # Initialize results
    results = []
    
    for threshold in thresholds
        # Q-value based selection
        qval_mask = protein_results[!, qval_col] .<= threshold
        qval_n = sum(qval_mask)
        qval_actual_fdr = qval_n > 0 ? sum(qval_mask .& (entrap_labels .!= 0)) / qval_n : 0.0
        
        # Combined EFDR based selection
        combined_mask = protein_results[!, combined_efdr_col] .<= threshold
        combined_n = sum(combined_mask)
        combined_actual_fdr = combined_n > 0 ? sum(combined_mask .& (entrap_labels .!= 0)) / combined_n : 0.0
        
        # Paired EFDR based selection
        paired_mask = protein_results[!, paired_efdr_col] .<= threshold
        paired_n = sum(paired_mask)
        paired_actual_fdr = paired_n > 0 ? sum(paired_mask .& (entrap_labels .!= 0)) / paired_n : 0.0
        
        push!(results, (
            threshold = threshold,
            qval_n = qval_n,
            qval_actual_fdr = qval_actual_fdr,
            combined_n = combined_n,
            combined_efdr = threshold,  # By definition
            combined_actual_fdr = combined_actual_fdr,
            paired_n = paired_n,
            paired_efdr = threshold,  # By definition
            paired_actual_fdr = paired_actual_fdr
        ))
    end
    
    return DataFrame(results)
end

"""
    calculate_protein_efdr_calibration_error(protein_results::DataFrame, qval_col::Symbol, efdr_col::Symbol)

Calculate calibration error for protein EFDR estimates.

# Arguments
- `protein_results::DataFrame`: Protein results with entrap_id and EFDR columns
- `qval_col::Symbol`: Q-value column
- `efdr_col::Symbol`: EFDR column to evaluate

# Returns
- Tuple of (calibration_data::DataFrame, mean_error::Float64)
"""
function calculate_protein_efdr_calibration_error(protein_results::DataFrame, qval_col::Symbol, efdr_col::Symbol)
    # Check required columns
    required_cols = [qval_col, efdr_col, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    
    # Extract entrapment labels
    entrap_labels = protein_results.entrap_id
    
    # Define thresholds
    thresholds = 0.001:0.001:0.1
    
    # Calculate calibration at each threshold
    calibration_data = []
    
    for threshold in thresholds
        mask = protein_results[!, efdr_col] .<= threshold
        n_selected = sum(mask)
        
        if n_selected > 0
            actual_fdr = sum(mask .& (entrap_labels .!= 0)) / n_selected
            error = abs(actual_fdr - threshold)
            
            push!(calibration_data, (
                threshold = threshold,
                estimated_fdr = threshold,
                actual_fdr = actual_fdr,
                n_selected = n_selected,
                error = error
            ))
        end
    end
    
    cal_df = DataFrame(calibration_data)
    mean_error = isempty(cal_df) ? 0.0 : mean(cal_df.error)
    
    return cal_df, mean_error
end