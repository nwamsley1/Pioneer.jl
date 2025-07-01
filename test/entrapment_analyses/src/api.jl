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
    
    # Separate global and per-file analyses
    verbose && println("Separating global and per-file analyses...")
    global_scores = Symbol[]
    perfile_scores = Symbol[]
    
    # Identify which scores are global vs per-file based on column name
    for (score_col, qval_col) in score_qval_pairs
        if occursin("global", String(score_col))
            push!(global_scores, score_col)
        else
            push!(perfile_scores, score_col)
        end
    end
    
    # Create global results dataframe if we have global scores
    global_results_df = nothing
    if !isempty(global_scores)
        verbose && println("Creating global results dataframe...")
        # Use the first global score for selection (typically global_prob)
        global_results_df = create_global_results_df(prec_results; score_col=global_scores[1])
        # Add entrap_pair_ids to global results
        add_entrap_pair_ids!(global_results_df, library_precursors)
        verbose && println("Global dataframe has $(nrow(global_results_df)) unique precursors")
    end
    
    # Process per-file scores on original dataframe
    if !isempty(perfile_scores)
        verbose && println("Processing per-file scores...")
        add_original_target_scores!(prec_results, library_precursors, perfile_scores)
        
        # Calculate EFDR for per-file scores
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        if !isempty(perfile_pairs)
            add_efdr_columns!(prec_results, library_precursors;
                             method_types=method_types,
                             score_qval_pairs=perfile_pairs,
                             r=r_lib)
        end
    end
    
    # Process global scores on global dataframe
    if !isnothing(global_results_df) && !isempty(global_scores)
        verbose && println("Processing global scores...")
        add_original_target_scores!(global_results_df, library_precursors, global_scores)
        
        # Calculate EFDR for global scores
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        if !isempty(global_pairs)
            add_efdr_columns!(global_results_df, library_precursors;
                             method_types=method_types,
                             score_qval_pairs=global_pairs,
                             r=r_lib)
        end
    end
    
    # For downstream analysis, use the appropriate dataframe
    # If we have both, we'll need to handle them separately
    analysis_df = if !isnothing(global_results_df) && !isempty(perfile_scores)
        # We have both - need to handle this case
        # For now, prioritize the original dataframe and note this limitation
        verbose && @warn "Both global and per-file scores detected. Comparison results will be generated separately."
        prec_results
    elseif !isnothing(global_results_df)
        global_results_df
    else
        prec_results
    end
    
    # Generate comparison results
    comparison_results = Dict{Tuple{Symbol,Symbol}, DataFrame}()
    
    # Process per-file comparisons
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            verbose && println("Comparing EFDR methods for $score_col/$qval_col (per-file)...")
            comparison_df = compare_efdr_methods(prec_results, qval_col, score_col, 
                                               library_precursors)
            comparison_results[(score_col, qval_col)] = comparison_df
        end
    end
    
    # Process global comparisons
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                verbose && println("Comparing EFDR methods for $score_col/$qval_col (global)...")
                comparison_df = compare_efdr_methods(global_results_df, qval_col, score_col, 
                                                   library_precursors)
                comparison_results[(score_col, qval_col)] = comparison_df
            end
        end
    end
    
    # Calculate calibration errors
    calibration_results = Dict{Symbol, Tuple{DataFrame, Float64}}()
    
    # Calibration for per-file scores
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            for method_type in method_types
                method_name = method_type == CombinedEFDR ? "combined" : "paired"
                efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                
                verbose && println("Calculating calibration error for $efdr_col (per-file)...")
                cal_data, cal_error = calculate_efdr_calibration_error(
                    prec_results, qval_col, efdr_col, library_precursors
                )
                calibration_results[efdr_col] = (cal_data, cal_error)
            end
        end
    end
    
    # Calibration for global scores
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                for method_type in method_types
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                    
                    verbose && println("Calculating calibration error for $efdr_col (global)...")
                    cal_data, cal_error = calculate_efdr_calibration_error(
                        global_results_df, qval_col, efdr_col, library_precursors
                    )
                    calibration_results[efdr_col] = (cal_data, cal_error)
                end
            end
        end
    end
    
    # Save plots
    verbose && println("Generating plots...")
    
    # Generate plots for per-file scores
    if !isempty(perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        save_efdr_plots(prec_results, output_dir;
                       score_qval_pairs=perfile_pairs,
                       method_types=method_types,
                       formats=plot_formats)
    end
    
    # Generate plots for global scores
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        save_efdr_plots(global_results_df, output_dir;
                       score_qval_pairs=global_pairs,
                       method_types=method_types,
                       formats=plot_formats)
    end
    
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
    # Use analysis_df for report generation (could be either dataframe)
    report_df = if !isnothing(global_results_df) && isempty(perfile_scores)
        global_results_df
    else
        prec_results
    end
    generate_markdown_report(report_df, library_precursors, comparison_results, 
                           calibration_results, score_qval_pairs, method_types,
                           original_rows, output_dir, report_path)
    push!(output_files, report_path)
    
    # Save processed data
    verbose && println("Saving processed data...")
    
    # Save per-file results if we processed them
    saved_dfs = Dict{String, DataFrame}()
    if !isempty(perfile_scores)
        data_path = joinpath(output_dir, "prec_results_with_efdr.arrow")
        Arrow.write(data_path, prec_results)
        push!(output_files, data_path)
        saved_dfs["per_file"] = prec_results
    end
    
    # Save global results if we processed them
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_data_path = joinpath(output_dir, "global_results_with_efdr.arrow")
        Arrow.write(global_data_path, global_results_df)
        push!(output_files, global_data_path)
        saved_dfs["global"] = global_results_df
    end
    
    verbose && println("\nAnalysis complete! Output saved to: $output_dir")
    
    # Return both dataframes if we have them
    return_data = if length(saved_dfs) == 2
        # Merge both dataframes or return separately
        saved_dfs
    elseif haskey(saved_dfs, "global")
        saved_dfs["global"]
    else
        saved_dfs["per_file"]
    end
    
    return (
        filtered_data = return_data,
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

"""
    run_protein_efdr_analysis(protein_results_path::String;
                              output_dir::String="efdr_out",
                              method_types=[CombinedEFDR, PairedEFDR],
                              score_qval_pairs=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
                              r_lib::Float64=1.0,
                              plot_formats=[:png, :pdf],
                              verbose::Bool=true)

Run empirical FDR analysis on protein-level data with entrapment sequences.

# Arguments
- `protein_results_path`: Path to Arrow file containing protein results
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

# Notes
- This function is specifically for protein-level data
- Unlike precursor analysis, no library file is needed (entrap_id is in results)
- Automatically detects global vs per-file scores based on column names
"""
function run_protein_efdr_analysis(protein_results_path::String;
                                  output_dir::String="efdr_out",
                                  method_types::Vector=[CombinedEFDR, PairedEFDR],
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
                                  r_lib::Float64=1.0,
                                  plot_formats::Vector{Symbol}=[:png, :pdf],
                                  verbose::Bool=true)
    
    # Create output directory
    mkpath(output_dir)
    output_files = String[]
    
    # Load data
    verbose && println("Loading protein data...")
    protein_results = DataFrame(Tables.columntable(Arrow.Table(protein_results_path)))
    is_sorted = issorted(protein_results, [:pg_score, :entrap_id], rev = [true, false])
    @info is_sorted ? "Protein results are sorted by pg_score and entrap_id" :
        "Protein results are NOT sorted by pg_score and entrap_id, sorting now"
    sort!(protein_results, [:pg_score, :entrap_id], rev = [true, false])
    
    verbose && println("Loaded $(nrow(protein_results)) protein results")
    
    # Verify this is protein data by checking for expected columns
    protein_cols = [:protein, :entrap_id, :pg_score]
    if !all(col -> hasproperty(protein_results, col), protein_cols)
        error("Input file does not appear to be protein data. Expected columns: $protein_cols")
    end
    
    # Check for target column and filter
    original_rows = nrow(protein_results)
    if hasproperty(protein_results, :target)
        non_target_rows = sum(.!protein_results.target)
        if non_target_rows > 0
            verbose && @warn "Filtering out $non_target_rows non-target (decoy) rows before EFDR calculation"
            filter!(x -> x.target, protein_results)
        end
    else
        throw("No 'target' column found.")
    end
    
    # Count proteins with entrapment versions
    proteins_with_entrapment = 0
    for protein_group in groupby(protein_results, [:species,:protein,:entrap_id,:file_name])
        unique_entrap_ids = unique(protein_group.entrap_id)
        if 0 in unique_entrap_ids && any(x -> x > 0, unique_entrap_ids)
            proteins_with_entrapment += 1
        end
    end
    verbose && println("Found $proteins_with_entrapment unique proteins with entrapment versions")
    
    # Separate global and per-file analyses
    verbose && println("Separating global and per-file analyses...")
    global_scores = Symbol[]
    perfile_scores = Symbol[]
    
    # Identify which scores are global vs per-file based on column name
    for (score_col, _) in score_qval_pairs
        if occursin("global", String(score_col))
            push!(global_scores, score_col)
        else
            push!(perfile_scores, score_col)
        end
    end
    
    # Create global results dataframe if we have global scores
    global_results_df = nothing
    if !isempty(global_scores)
        verbose && println("Creating global protein results dataframe...")
        # Use the first global score for selection (typically global_pg_score)
        global_results_df = create_global_protein_results_df(protein_results; score_col=global_scores[1])
        verbose && println("Global dataframe has $(nrow(global_results_df)) unique proteins")
    end
    
    # Process per-file scores on original dataframe
    if !isempty(perfile_scores)
        verbose && println("Processing per-file protein scores...")
        add_original_target_protein_scores!(protein_results, perfile_scores)
        
        # Calculate EFDR for per-file scores
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        if !isempty(perfile_pairs)
            add_protein_efdr_columns!(protein_results;
                                     method_types=method_types,
                                     score_qval_pairs=perfile_pairs,
                                     r=r_lib)
        end
    end
    
    # Process global scores on global dataframe
    if !isnothing(global_results_df) && !isempty(global_scores)
        verbose && println("Processing global protein scores...")
        add_original_target_protein_scores!(global_results_df, global_scores)
        
        # Calculate EFDR for global scores
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        if !isempty(global_pairs)
            add_protein_efdr_columns!(global_results_df;
                                     method_types=method_types,
                                     score_qval_pairs=global_pairs,
                                     r=r_lib)
        end
    end
    
    # Generate comparison results
    comparison_results = Dict{Tuple{Symbol,Symbol}, DataFrame}()
    
    # Process per-file comparisons
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            verbose && println("Comparing EFDR methods for $score_col/$qval_col (per-file)...")
            comparison_df = compare_protein_efdr_methods(protein_results, qval_col, score_col)
            comparison_results[(score_col, qval_col)] = comparison_df
        end
    end
    
    # Process global comparisons
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                verbose && println("Comparing EFDR methods for $score_col/$qval_col (global)...")
                comparison_df = compare_protein_efdr_methods(global_results_df, qval_col, score_col)
                comparison_results[(score_col, qval_col)] = comparison_df
            end
        end
    end
    
    # Calculate calibration errors
    calibration_results = Dict{Symbol, Tuple{DataFrame, Float64}}()
    
    # Calibration for per-file scores
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            for method_type in method_types
                method_name = method_type == CombinedEFDR ? "combined" : "paired"
                efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                
                if hasproperty(protein_results, efdr_col)
                    verbose && println("Calculating calibration error for $efdr_col (per-file)...")
                    cal_data, cal_error = calculate_protein_efdr_calibration_error(
                        protein_results, qval_col, efdr_col
                    )
                    calibration_results[efdr_col] = (cal_data, cal_error)
                end
            end
        end
    end
    
    # Calibration for global scores
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                for method_type in method_types
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                    
                    if hasproperty(global_results_df, efdr_col)
                        verbose && println("Calculating calibration error for $efdr_col (global)...")
                        cal_data, cal_error = calculate_protein_efdr_calibration_error(
                            global_results_df, qval_col, efdr_col
                        )
                        calibration_results[efdr_col] = (cal_data, cal_error)
                    end
                end
            end
        end
    end
    
    # Save plots
    verbose && println("Generating plots...")
    
    # Generate plots for per-file scores
    if !isempty(perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        save_efdr_plots(protein_results, output_dir;
                       score_qval_pairs=perfile_pairs,
                       method_types=method_types,
                       formats=plot_formats)
    end
    
    # Generate plots for global scores
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        save_efdr_plots(global_results_df, output_dir;
                       score_qval_pairs=global_pairs,
                       method_types=method_types,
                       formats=plot_formats)
    end
    
    # Add plot files to output list
    for (score_col, _) in score_qval_pairs
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_$(score_col).$(format)"))
        end
    end
    
    # Generate markdown report
    verbose && println("Generating markdown report...")
    report_path = joinpath(output_dir, "protein_efdr_analysis_report.md")
    
    # For report, we'll use a simplified version
    open(report_path, "w") do io
        println(io, "# Protein-Level Empirical FDR Analysis Report")
        println(io, "\nGenerated: $(now())")
        println(io)
        println(io, "## Data Summary")
        println(io, "- Original protein results: $(original_rows)")
        println(io, "- Analyzed protein results: $(nrow(protein_results))")
        println(io, "- Unique proteins with entrapment versions: $proteins_with_entrapment")
        println(io)
        
        # Entrapment composition
        if hasproperty(protein_results, :entrap_id)
            println(io, "### Entrapment Composition")
            entrap_counts = combine(groupby(protein_results, :entrap_id), nrow => :count)
            sort!(entrap_counts, :entrap_id)
            println(io)
            println(io, "| Entrap ID | Count | Percentage |")
            println(io, "|-----------|-------|------------|")
            total = nrow(protein_results)
            for row in eachrow(entrap_counts)
                percentage = row.count / total * 100
                group_name = row.entrap_id == 0 ? "Target (0)" : "Entrapment ($(row.entrap_id))"
                @printf(io, "| %s | %d | %.2f%% |\n", group_name, row.count, percentage)
            end
        end
        
        println(io)
        println(io, "## EFDR Analysis Results")
        println(io)
        println(io, "See generated plots for EFDR comparison results.")
        println(io)
        println(io, "## Analysis Parameters")
        println(io, "- EFDR Methods: ", join([m == CombinedEFDR ? "Combined" : "Paired" for m in method_types], ", "))
        println(io, "- Score/Q-value pairs analyzed: ", join(["$s/$q" for (s,q) in score_qval_pairs], ", "))
        println(io, "- Output directory: `$output_dir`")
    end
    push!(output_files, report_path)
    
    # Save processed data
    verbose && println("Saving processed data...")
    
    # Save per-file results if we processed them
    saved_dfs = Dict{String, DataFrame}()
    if !isempty(perfile_scores)
        data_path = joinpath(output_dir, "protein_results_with_efdr.arrow")
        Arrow.write(data_path, protein_results)
        push!(output_files, data_path)
        saved_dfs["per_file"] = protein_results
    end
    
    # Save global results if we processed them
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_data_path = joinpath(output_dir, "global_protein_results_with_efdr.arrow")
        Arrow.write(global_data_path, global_results_df)
        push!(output_files, global_data_path)
        saved_dfs["global"] = global_results_df
    end
    
    verbose && println("\nProtein analysis complete! Output saved to: $output_dir")
    
    # Return results
    return_data = if length(saved_dfs) == 2
        saved_dfs
    elseif haskey(saved_dfs, "global")
        saved_dfs["global"]
    else
        saved_dfs["per_file"]
    end
    
    return (
        filtered_data = return_data,
        comparison_results = comparison_results,
        calibration_results = calibration_results,
        output_files = output_files
    )
end

# Export the main API functions
export run_efdr_analysis, run_protein_efdr_analysis