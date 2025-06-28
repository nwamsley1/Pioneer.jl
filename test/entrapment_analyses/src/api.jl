using Arrow
using DataFrames
using Printf
using Markdown

"""
    run_efdr_analysis(prec_results_path::String, library_precursors_path::String;
                      output_dir::String="efdr_out",
                      method_types=[CombinedEFDR, PairedEFDR],
                      score_qval_pairs=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                      r_lib::Float64=1.0,
                      plot_formats=[:png, :pdf],
                      verbose::Bool=true)

Run empirical FDR analysis on proteomics data with entrapment sequences.

# Arguments
- `prec_results_path`: Path to Arrow file containing precursor results
- `library_precursors_path`: Path to Arrow file containing library precursors with entrapment info
- `output_dir`: Directory to save all output files (default: "efdr_out")
- `method_types`: EFDR calculation methods to use
- `score_qval_pairs`: Vector of (score_col, qval_col) tuples to analyze
- `r_lib`: Ratio of library entrapments to real entrapments
- `plot_formats`: File formats for saving plots
- `verbose`: Whether to print progress messages

# Returns
- `results`: NamedTuple containing:
  - `filtered_data`: DataFrame with EFDR columns added
  - `comparison_results`: Dict of comparison DataFrames for each score/qval pair
  - `calibration_results`: Dict of calibration results for each EFDR column
  - `output_files`: Vector of paths to generated output files
"""
function run_efdr_analysis(prec_results_path::String, library_precursors_path::String;
                          output_dir::String="efdr_out",
                          method_types::Vector=[CombinedEFDR, PairedEFDR],
                          score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                          r_lib::Float64=1.0,
                          plot_formats::Vector{Symbol}=[:png, :pdf],
                          verbose::Bool=true)
    
    # Create output directory
    mkpath(output_dir)
    output_files = String[]
    
    # Load data
    verbose && println("Loading data...")
    prec_results = DataFrame(Arrow.Table(prec_results_path))
    library_precursors = DataFrame(Arrow.Table(library_precursors_path))
    
    verbose && println("Loaded $(nrow(prec_results)) precursor results")
    verbose && println("Loaded $(nrow(library_precursors)) library precursors")
    
    # Check for target column and filter
    original_rows = nrow(prec_results)
    if hasproperty(prec_results, :target)
        non_target_rows = sum(.!prec_results.target)
        if non_target_rows > 0
            verbose && @warn "Filtering out $non_target_rows non-target (decoy) rows before EFDR calculation"
            filter!(x -> x.target, prec_results)
        end
    else
        verbose && @warn "No 'target' column found. Assuming all rows are targets."
    end
    
    # Add mod_key column if not present
    if !hasproperty(library_precursors, :mod_key)
        verbose && println("Adding mod_key column to library...")
        library_precursors[!, :mod_key] = map(x -> getModKey(x), library_precursors.structural_mods)
    end
    
    # Assign entrapment pairs
    verbose && println("Assigning entrapment pairs...")
    assign_entrapment_pairs!(library_precursors)
    
    # Add entrap_pair_ids to results
    add_entrap_pair_ids!(prec_results, library_precursors)
    
    # Add original target scores
    verbose && println("Adding original target scores...")
    score_cols = [pair[1] for pair in score_qval_pairs]
    add_original_target_scores!(prec_results, library_precursors, score_cols)
    
    # Calculate EFDR columns
    verbose && println("Calculating EFDR columns...")
    add_efdr_columns!(prec_results, library_precursors;
                     method_types=method_types,
                     score_qval_pairs=score_qval_pairs,
                     r=r_lib)
    
    # Generate comparison results
    comparison_results = Dict{Tuple{Symbol,Symbol}, DataFrame}()
    for (score_col, qval_col) in score_qval_pairs
        verbose && println("Comparing EFDR methods for $score_col/$qval_col...")
        comparison_df = compare_efdr_methods(prec_results, qval_col, score_col, 
                                           library_precursors)
        comparison_results[(score_col, qval_col)] = comparison_df
    end
    
    # Calculate calibration errors
    calibration_results = Dict{Symbol, Tuple{DataFrame, Float64}}()
    for (score_col, qval_col) in score_qval_pairs
        for method_type in method_types
            method_name = method_type == CombinedEFDR ? "combined" : "paired"
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            
            verbose && println("Calculating calibration error for $efdr_col...")
            cal_data, cal_error = calculate_efdr_calibration_error(
                prec_results, qval_col, efdr_col, library_precursors
            )
            calibration_results[efdr_col] = (cal_data, cal_error)
        end
    end
    
    # Save plots
    verbose && println("Generating plots...")
    save_efdr_plots(prec_results, output_dir;
                   score_qval_pairs=score_qval_pairs,
                   method_types=method_types,
                   formats=plot_formats)
    
    # Add plot files to output list
    for (score_col, _) in score_qval_pairs
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_$(score_col).$(format)"))
        end
    end
    if length(score_qval_pairs) > 1
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_all.$(format)"))
        end
    end
    
    # Generate markdown report
    verbose && println("Generating markdown report...")
    report_path = joinpath(output_dir, "efdr_analysis_report.md")
    generate_markdown_report(prec_results, library_precursors, comparison_results, 
                           calibration_results, score_qval_pairs, method_types,
                           original_rows, output_dir, report_path)
    push!(output_files, report_path)
    
    # Save processed data
    verbose && println("Saving processed data...")
    data_path = joinpath(output_dir, "prec_results_with_efdr.arrow")
    Arrow.write(data_path, prec_results)
    push!(output_files, data_path)
    
    verbose && println("\nAnalysis complete! Output saved to: $output_dir")
    
    return (
        filtered_data = prec_results,
        comparison_results = comparison_results,
        calibration_results = calibration_results,
        output_files = output_files
    )
end

"""
    generate_markdown_report(prec_results, library_precursors, comparison_results,
                           calibration_results, score_qval_pairs, method_types,
                           original_rows, output_dir, output_path)

Generate a comprehensive markdown report of the EFDR analysis.
"""
function generate_markdown_report(prec_results::DataFrame, library_precursors::DataFrame,
                                comparison_results::Dict, calibration_results::Dict,
                                score_qval_pairs::Vector{Tuple{Symbol,Symbol}},
                                method_types::Vector, original_rows::Int,
                                output_dir::String, output_path::String)
    
    open(output_path, "w") do io
        # Header
        println(io, "# Empirical FDR Analysis Report")
        println(io, "\nGenerated: $(now())")
        println(io)
        
        # Data summary
        println(io, "## Data Summary")
        println(io)
        println(io, "- Original precursor results: $(original_rows)")
        filtered_rows = original_rows - nrow(prec_results)
        if filtered_rows > 0
            println(io, "- Filtered non-target rows: $(filtered_rows)")
        end
        println(io, "- Analyzed precursor results: $(nrow(prec_results))")
        println(io, "- Library precursors: $(nrow(library_precursors))")
        println(io)
        
        # Library composition
        println(io, "### Library Composition")
        entrap_counts = combine(groupby(library_precursors, :entrapment_group_id), nrow => :count)
        sort!(entrap_counts, :entrapment_group_id)
        println(io)
        println(io, "| Entrapment Group | Count | Percentage |")
        println(io, "|-----------------|-------|------------|")
        total_lib = nrow(library_precursors)
        for row in eachrow(entrap_counts)
            percentage = row.count / total_lib * 100
            group_name = row.entrapment_group_id == 0 ? "Target (0)" : "Entrapment ($(row.entrapment_group_id))"
            @printf(io, "| %s | %d | %.2f%% |\n", group_name, row.count, percentage)
        end
        println(io)
        
        # Results composition
        println(io, "### Results Composition")
        result_entrap_labels = [library_precursors.entrapment_group_id[pid] 
                               for pid in prec_results.precursor_idx]
        result_counts = combine(groupby(DataFrame(entrapment_group_id=result_entrap_labels), 
                                      :entrapment_group_id), nrow => :count)
        sort!(result_counts, :entrapment_group_id)
        println(io)
        println(io, "| Entrapment Group | Count | Percentage |")
        println(io, "|-----------------|-------|------------|")
        total_results = nrow(prec_results)
        for row in eachrow(result_counts)
            percentage = row.count / total_results * 100
            group_name = row.entrapment_group_id == 0 ? "Target (0)" : "Entrapment ($(row.entrapment_group_id))"
            @printf(io, "| %s | %d | %.2f%% |\n", group_name, row.count, percentage)
        end
        println(io)
        
        # EFDR comparison results
        println(io, "## EFDR Method Comparison")
        println(io)
        
        for (score_col, qval_col) in score_qval_pairs
            println(io, "### $(String(score_col)) / $(String(qval_col))")
            println(io)
            
            comparison_df = comparison_results[(score_col, qval_col)]
            
            println(io, "| Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR |")
            println(io, "|-----------|-----------|------------|--------------|---------------|------------|-------------|")
            
            for row in eachrow(comparison_df)
                @printf(io, "| %.3f | %d | %.4f | %d | %.4f | %d | %.4f |\n",
                       row.threshold, row.qval_n, row.qval_actual_fdr,
                       row.combined_n, row.combined_efdr,
                       row.paired_n, row.paired_efdr)
            end
            println(io)
        end
        
        # Calibration results
        println(io, "## Calibration Analysis")
        println(io)
        println(io, "Mean absolute calibration errors:")
        println(io)
        println(io, "| Method | Mean Calibration Error |")
        println(io, "|--------|----------------------|")
        
        for (efdr_col, (_, cal_error)) in calibration_results
            @printf(io, "| %s | %.4f |\n", String(efdr_col), cal_error)
        end
        println(io)
        
        # Plots section
        println(io, "## Plots")
        println(io)
        println(io, "The following plots have been generated:")
        println(io)
        
        for (score_col, _) in score_qval_pairs
            println(io, "- ![EFDR Comparison for $(score_col)](efdr_comparison_$(score_col).png)")
        end
        
        if length(score_qval_pairs) > 1
            println(io, "- ![All EFDR Comparisons](efdr_comparison_all.png)")
        end
        
        println(io)
        println(io, "## Analysis Parameters")
        println(io)
        println(io, "- EFDR Methods: ", join([m == CombinedEFDR ? "Combined" : "Paired" for m in method_types], ", "))
        println(io, "- Score/Q-value pairs analyzed: ", join(["$s/$q" for (s,q) in score_qval_pairs], ", "))
        println(io, "- Output directory: `$output_dir`")
    end
end

# Export the main API function
export run_efdr_analysis