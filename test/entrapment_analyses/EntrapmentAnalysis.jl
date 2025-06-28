module EntrapmentAnalysis

using DataFrames
using Arrow
using Tables
using Printf
using Dates
using CSV
using Plots
using Statistics

# Include all source files
include("src/core/efdr_methods.jl")
include("src/core/entrapment_pairing.jl") 
include("src/core/scoring.jl")
include("src/analysis/efdr_analysis.jl")
include("src/analysis/calibration.jl")
include("src/plotting/efdr_plots.jl")
include("src/api.jl")

# Export main API
export run_efdr_analysis

# Export types
export EFDRMethod, CombinedEFDR, PairedEFDR

# Export core functions
export calculate_efdr, add_efdr_columns!
export assign_entrapment_pairs!, add_entrap_pair_ids!
export add_original_target_scores!, get_complement_score
export getModKey

# Export analysis functions  
export compare_efdr_methods, calculate_efdr_calibration_error
export analyze_efdr_at_threshold, print_efdr_comparison_table

# Export plotting functions
export plot_efdr_comparison, plot_efdr_vs_qval
export plot_multiple_efdr_comparisons, save_efdr_plots

# Re-export commonly used functions from dependencies
export DataFrame, @sprintf

end # module