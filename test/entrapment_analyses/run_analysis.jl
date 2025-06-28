#!/usr/bin/env julia

# Simple script to run EFDR analysis using the EntrapmentAnalysis module
# Usage: julia run_analysis.jl <prec_results.arrow> <library.arrow> [output_dir]

# Add module to load path
push!(LOAD_PATH, @__DIR__)

using EntrapmentAnalysis

function main(args)
    if length(args) < 2
        println("Usage: julia run_analysis.jl <prec_results.arrow> <library.arrow> [output_dir]")
        println()
        println("Arguments:")
        println("  prec_results.arrow - Path to precursor results Arrow file")
        println("  library.arrow      - Path to library precursors Arrow file")
        println("  output_dir         - Output directory (optional, default: efdr_out)")
        return 1
    end
    
    prec_results_path = args[1]
    library_path = args[2]
    output_dir = length(args) >= 3 ? args[3] : "efdr_out"
    
    # Check if files exist
    if !isfile(prec_results_path)
        error("Precursor results file not found: $prec_results_path")
    end
    
    if !isfile(library_path)
        error("Library file not found: $library_path")
    end
    
    println("Running EFDR analysis...")
    println("  Precursor results: $prec_results_path")
    println("  Library: $library_path")
    println("  Output directory: $output_dir")
    println()
    
    # Run analysis
    results = run_efdr_analysis(
        prec_results_path,
        library_path;
        output_dir = output_dir,
        verbose = true
    )
    
    println("\nAnalysis complete!")
    println("Results saved to: $output_dir/")
    
    return 0
end

# Run if called as script
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end