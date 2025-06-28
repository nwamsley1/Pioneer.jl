# Example: EFDR analysis on yeast data with entrapment sequences
# This example shows how to use the EntrapmentAnalysis module with real data paths

# Add parent directory to load path
push!(LOAD_PATH, dirname(dirname(@__DIR__)))

using EntrapmentAnalysis

# Define input file paths (update these to your actual file locations)
prec_results_path = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/precursors_long.arrow"
library_path = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_entrap.poin/precursors_table.arrow"

# Check if files exist
if !isfile(prec_results_path)
    @warn "Precursor results file not found at: $prec_results_path"
    println("Please update the path to your precursor results file")
end

if !isfile(library_path)
    @warn "Library file not found at: $library_path" 
    println("Please update the path to your library file")
end

# Run the complete EFDR analysis
println("Starting EFDR analysis on yeast data...")

results = run_efdr_analysis(
    prec_results_path,
    library_path;
    output_dir = "yeast_efdr_analysis",
    method_types = [CombinedEFDR, PairedEFDR],
    score_qval_pairs = [(:global_prob, :global_qval), (:prec_prob, :qval)],
    r_lib = 1.0,  # Will be calculated automatically from the data
    plot_formats = [:png, :pdf],
    verbose = true
)

println("\nAnalysis complete!")
println("Results saved to: yeast_efdr_analysis/")
println("\nGenerated files:")
for file in results.output_files
    println("  - $file")
end

# Access specific results
println("\nAccessing analysis results...")

# Get comparison results for global scores
global_comparison = results.comparison_results[(:global_prob, :global_qval)]
println("\nGlobal score EFDR comparison at 1% FDR threshold:")
row_1pct = filter(row -> row.threshold == 0.01, global_comparison)[1, :]
println("  Q-value IDs: $(row_1pct.qval_n)")
println("  Actual FDR: $(round(row_1pct.qval_actual_fdr, digits=4))")
println("  Combined EFDR IDs: $(row_1pct.combined_n)")
println("  Combined EFDR: $(round(row_1pct.combined_efdr, digits=4))")
println("  Paired EFDR IDs: $(row_1pct.paired_n)")
println("  Paired EFDR: $(round(row_1pct.paired_efdr, digits=4))")

# Get calibration results
println("\nCalibration errors:")
for (efdr_col, (_, cal_error)) in results.calibration_results
    println("  $efdr_col: $(round(cal_error, digits=4))")
end

println("\nCheck the markdown report for detailed analysis: yeast_efdr_analysis/efdr_analysis_report.md")