using DataFrames
using Statistics

"""
    analyze_efdr_at_threshold(df::DataFrame, qval_col::Symbol, efdr_col::Symbol, threshold::Float64)

Analyze EFDR performance at a specific q-value threshold.

Returns a NamedTuple with:
- threshold: The q-value threshold used
- n_passing: Number of identifications passing threshold
- empirical_fdr: The empirical FDR at this threshold
- n_targets: Number of targets passing
- n_entrapments: Number of entrapments passing
"""
function analyze_efdr_at_threshold(df::DataFrame, qval_col::Symbol, efdr_col::Symbol, 
                                 threshold::Float64, library_precursors::DataFrame)
    # Sort by q-value
    sorted_indices = sortperm(df[!, qval_col])
    sorted_df = df[sorted_indices, :]
    
    # Find identifications passing threshold
    passing_mask = sorted_df[!, qval_col] .<= threshold
    passing_df = sorted_df[passing_mask, :]
    
    if nrow(passing_df) == 0
        return (threshold=threshold, n_passing=0, empirical_fdr=0.0, 
                n_targets=0, n_entrapments=0)
    end
    
    # Get entrapment labels
    entrap_labels = [library_precursors.entrapment_group_id[pid] 
                     for pid in passing_df.precursor_idx]
    
    n_targets = sum(entrap_labels .== 0)
    n_entrapments = sum(entrap_labels .> 0)
    
    # Get EFDR at the threshold (last passing value)
    empirical_fdr = passing_df[end, efdr_col]
    
    return (
        threshold = threshold,
        n_passing = nrow(passing_df),
        empirical_fdr = empirical_fdr,
        n_targets = n_targets,
        n_entrapments = n_entrapments
    )
end

"""
    compare_efdr_methods(df::DataFrame, qval_col::Symbol, score_col::Symbol, 
                        library_precursors::DataFrame; 
                        thresholds=[0.001, 0.01, 0.05, 0.1])

Compare different EFDR methods at multiple q-value thresholds.
"""
function compare_efdr_methods(df::DataFrame, qval_col::Symbol, score_col::Symbol,
                            library_precursors::DataFrame;
                            thresholds::Vector{Float64}=[0.001, 0.01, 0.05, 0.1])
    
    combined_col = Symbol(String(score_col) * "_combined_efdr")
    paired_col = Symbol(String(score_col) * "_paired_efdr")
    
    results = DataFrame()
    
    for threshold in thresholds
        # Standard q-value analysis
        qval_sorted = sortperm(df[!, qval_col])
        first_failing = findfirst(df[qval_sorted, qval_col] .> threshold)
        if isnothing(first_failing)
            # All pass the threshold
            qval_passing = df[qval_sorted, :]
        else
            qval_passing = df[qval_sorted[1:first_failing-1], :]
        end
        qval_entrap_labels = [library_precursors.entrapment_group_id[pid] 
                             for pid in qval_passing.precursor_idx]
        
        # Combined EFDR analysis
        combined_result = analyze_efdr_at_threshold(df, qval_col, combined_col, 
                                                   threshold, library_precursors)
        
        # Paired EFDR analysis
        paired_result = analyze_efdr_at_threshold(df, qval_col, paired_col,
                                                 threshold, library_precursors)
        
        push!(results, (
            threshold = threshold,
            qval_n = nrow(qval_passing),
            qval_actual_fdr = sum(qval_entrap_labels .> 0) / max(1, length(qval_entrap_labels)),
            combined_n = combined_result.n_passing,
            combined_efdr = combined_result.empirical_fdr,
            paired_n = paired_result.n_passing,
            paired_efdr = paired_result.empirical_fdr
        ))
    end
    
    return results
end

"""
    print_efdr_comparison_table(comparison_df::DataFrame)

Pretty-print the EFDR comparison results.
"""
function print_efdr_comparison_table(comparison_df::DataFrame)
    println("\nEFDR Method Comparison")
    println("="^80)
    println("Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR")
    println("-"^80)
    
    for row in eachrow(comparison_df)
        @printf("%9.3f | %9d | %10.4f | %12d | %13.4f | %10d | %11.4f\n",
                row.threshold, row.qval_n, row.qval_actual_fdr,
                row.combined_n, row.combined_efdr,
                row.paired_n, row.paired_efdr)
    end
    println("="^80)
end

"""
    calculate_efdr_calibration_error(df::DataFrame, qval_col::Symbol, efdr_col::Symbol,
                                   library_precursors::DataFrame;
                                   n_bins=20)

Calculate calibration error between estimated and actual FDR.
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

export analyze_efdr_at_threshold, compare_efdr_methods, print_efdr_comparison_table, 
       calculate_efdr_calibration_error