# Basic usage example for EntrapmentAnalysis module

# Add parent directory to load path
push!(LOAD_PATH, dirname(dirname(@__DIR__)))

using EntrapmentAnalysis

# Example 1: Basic EFDR analysis with default parameters
println("Running basic EFDR analysis...")

# Provide paths to your data files
prec_results_path = "path/to/your/precursor_results.arrow"
library_path = "path/to/your/library_precursors.arrow"

# Run analysis with default parameters
results = run_efdr_analysis(prec_results_path, library_path)

println("Analysis complete! Results saved to: efdr_out/")

# Example 2: Custom parameters
println("\nRunning EFDR analysis with custom parameters...")

results_custom = run_efdr_analysis(
    prec_results_path, 
    library_path;
    output_dir = "custom_efdr_output",
    method_types = [CombinedEFDR, PairedEFDR],
    score_qval_pairs = [(:global_prob, :global_qval), (:prec_prob, :qval)],
    r_lib = 1.0,
    plot_formats = [:png, :pdf, :svg],
    verbose = true
)

# Access the results
filtered_data = results_custom.filtered_data
comparison_results = results_custom.comparison_results
calibration_results = results_custom.calibration_results

# Example 3: Working with the data directly (handling global vs per-file scores)
println("\nExample of direct data manipulation with separate global analysis...")

using DataFrames, Arrow

# Load your data
prec_results = DataFrame(Arrow.Table(prec_results_path))
library_precursors = DataFrame(Arrow.Table(library_path))

# Add mod_key column if needed
if !hasproperty(library_precursors, :mod_key)
    library_precursors[!, :mod_key] = map(x -> getModKey(x), library_precursors.structural_mods)
end

# Assign entrapment pairs
assign_entrapment_pairs!(library_precursors)

# Add entrap_pair_ids to results
add_entrap_pair_ids!(prec_results, library_precursors)

# Create separate dataframe for global score analysis
println("Creating global results dataframe...")
global_results_df = create_global_results_df(prec_results; score_col=:global_prob)
add_entrap_pair_ids!(global_results_df, library_precursors)
println("Global dataframe has $(nrow(global_results_df)) unique precursors")

# Process per-file scores on original dataframe
add_original_target_scores!(prec_results, library_precursors, [:prec_prob])
add_efdr_columns!(prec_results, library_precursors;
                 method_types = [CombinedEFDR, PairedEFDR],
                 score_qval_pairs = [(:prec_prob, :qval)],
                 r = 1.0)

# Process global scores on global dataframe
add_original_target_scores!(global_results_df, library_precursors, [:global_prob])
add_efdr_columns!(global_results_df, library_precursors;
                 method_types = [CombinedEFDR, PairedEFDR],
                 score_qval_pairs = [(:global_prob, :global_qval)],
                 r = 1.0)

println("EFDR columns added to both dataframes!")

# Example 4: Plotting only
println("\nExample of generating plots from existing data...")

# If you already have a dataframe with EFDR columns
save_efdr_plots(prec_results, "my_plots";
               score_qval_pairs = [(:global_prob, :global_qval)],
               method_types = [CombinedEFDR, PairedEFDR],
               formats = [:png])

# Example 5: Saving results from both analyses
println("\nExample of saving both per-file and global analysis results...")

# Save per-file results
Arrow.write("per_file_results_with_efdr.arrow", prec_results)

# Save global results
Arrow.write("global_results_with_efdr.arrow", global_results_df)

println("Both dataframes saved!")

println("\nAll examples completed!")