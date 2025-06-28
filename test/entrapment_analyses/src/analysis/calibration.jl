"""
    calculate_efdr_calibration_error(df::DataFrame, qval_col::Symbol, efdr_col::Symbol,
                                   library_precursors::DataFrame;
                                   n_bins=20)

Calculate calibration error between estimated and actual FDR.

Returns:
- calibration_data: DataFrame with bin-wise comparison of estimated vs actual FDR
- mean_calibration_error: Mean absolute error across all bins
"""
function calculate_efdr_calibration_error(df::DataFrame, qval_col::Symbol, efdr_col::Symbol,
                                        library_precursors::DataFrame;
                                        n_bins::Int=20)
    # Sort by EFDR
    sorted_indices = sortperm(df[!, efdr_col])
    sorted_df = df[sorted_indices, :]
    
    # Get entrapment labels
    entrap_labels = [library_precursors.entrapment_group_id[pid] 
                     for pid in sorted_df.precursor_idx]
    
    # Bin the data
    bin_size = div(nrow(sorted_df), n_bins)
    calibration_data = DataFrame()
    
    for i in 1:n_bins
        start_idx = (i-1) * bin_size + 1
        end_idx = min(i * bin_size, nrow(sorted_df))
        
        if start_idx <= end_idx
            bin_data = sorted_df[start_idx:end_idx, :]
            bin_labels = entrap_labels[start_idx:end_idx]
            
            # Average estimated EFDR in bin
            avg_estimated = mean(bin_data[!, efdr_col])
            
            # Actual FDR up to this point
            all_labels_up_to_here = entrap_labels[1:end_idx]
            actual_fdr = sum(all_labels_up_to_here .> 0) / length(all_labels_up_to_here)
            
            push!(calibration_data, (
                bin = i,
                estimated_efdr = avg_estimated,
                actual_fdr = actual_fdr,
                error = abs(avg_estimated - actual_fdr),
                n_items = end_idx - start_idx + 1
            ))
        end
    end
    
    mean_calibration_error = mean(calibration_data.error)
    
    return calibration_data, mean_calibration_error
end