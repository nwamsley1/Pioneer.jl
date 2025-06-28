# Load required packages and functions
using DataFrames
using Arrow
using Tables
using Test
using Printf
using Dates
using CSV
include("entrapment_helper_funcs.jl")
include("efdr_funcs.jl")

prec_results = DataFrame(
    Tables.columntable(
        Arrow.Table(
    "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/precursors_long.arrow"
)))
library_precursors = DataFrame(Arrow.Table(
    "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_entrap.poin/precursors_table.arrow"
))

# Create output directory
output_dir = joinpath(dirname(@__FILE__), "efdr_out")
mkpath(output_dir)

# Open summary file for writing
summary_file = open(joinpath(output_dir, "efdr_analysis_summary.md"), "w")

# Helper function to write to both console and file
function write_both(io::IO, str::String)
    println(str)
    println(io, str)
end

# Write header
write_both(summary_file, "# EFDR Analysis Summary")
write_both(summary_file, "")
write_both(summary_file, "Date: $(Dates.now())")
write_both(summary_file, "")

# Check the results
write_both(summary_file, "## Initial Data Summary")
write_both(summary_file, "")
write_both(summary_file, "- Total pairs created: $(length(unique(library_precursors.pair_id)))")
write_both(summary_file, "- Rows without pairs: $(sum(ismissing.(library_precursors.pair_id)))")

# Add mod_key column to library_precursors dataframe
library_precursors[!, :mod_key] = map(x -> getModKey(x), library_precursors.structural_mods)
assign_entrapment_pairs!(library_precursors)
const r_lib = sum(library_precursors[!,:entrapment_group_id].!=0)/sum(library_precursors[!,:entrapment_group_id].==0)

# Check the results after entrapment pairing
write_both(summary_file, "")
write_both(summary_file, "## After Entrapment Pairing")
write_both(summary_file, "")
write_both(summary_file, "- Total entrapment pairs created: $(length(unique(library_precursors.entrap_pair_id)))")
write_both(summary_file, "- Rows without entrapment pairs: $(sum(ismissing.(library_precursors.entrap_pair_id)))")
write_both(summary_file, "- Entrapment ratio (r): $(round(r_lib, digits=4))")

filter!(x->x.target, prec_results)
add_efdr_columns!(prec_results, library_precursors;
                    method_types = [CombinedEFDR, PairedEFDR],
                    score_qval_pairs = [(:global_prob, :global_qval), (:prec_prob, :qval)],
                    r = r_lib)

# Include plotting and analysis functions
include("efdr_plotting.jl")
include("efdr_analysis.jl")
using Plots

# 1. Create comparison plot
write_both(summary_file, "")
write_both(summary_file, "## EFDR Analysis Results")
write_both(summary_file, "")
write_both(summary_file, "### Global Score Analysis")
write_both(summary_file, "")
write_both(summary_file, "Creating EFDR comparison plots...")

p = plot_efdr_comparison(prec_results, :global_prob, :global_qval;
                        title="EFDR Comparison: Yeast Entrapment Analysis",
                        xlims=(0, 0.015),
                        ylims=(0, 0.015))
display(p)

# Save the plot
savefig(p, joinpath(output_dir, "efdr_comparison_global_prob.png"))
savefig(p, joinpath(output_dir, "efdr_comparison_global_prob.pdf"))
write_both(summary_file, "")
write_both(summary_file, "Plots saved:")
write_both(summary_file, "- ![Global EFDR Comparison](efdr_comparison_global_prob.png)")

# 2. Analyze EFDR performance at standard thresholds
write_both(summary_file, "")
write_both(summary_file, "#### Global EFDR Performance at Standard Thresholds")
write_both(summary_file, "")
comparison = compare_efdr_methods(prec_results, :global_qval, :global_prob, 
                                library_precursors;
                                thresholds=[0.001, 0.01, 0.05])

# Format comparison table for markdown
write_both(summary_file, "| Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR |")
write_both(summary_file, "|-----------|-----------|------------|--------------|---------------|------------|-------------|")
for row in eachrow(comparison)
    write_both(summary_file, @sprintf("| %.3f | %d | %.4f | %d | %.4f | %d | %.4f |",
            row.threshold, row.qval_n, row.qval_actual_fdr,
            row.combined_n, row.combined_efdr,
            row.paired_n, row.paired_efdr))
end

# 3. Calculate calibration error
write_both(summary_file, "")
write_both(summary_file, "#### Global Calibration Error")
write_both(summary_file, "")
cal_combined, error_combined = calculate_efdr_calibration_error(
    prec_results, :global_qval, :global_prob_combined_efdr, library_precursors)
cal_paired, error_paired = calculate_efdr_calibration_error(
    prec_results, :global_qval, :global_prob_paired_efdr, library_precursors)

write_both(summary_file, "Mean calibration error:")
write_both(summary_file, "- Combined EFDR: $(round(error_combined, digits=4))")
write_both(summary_file, "- Paired EFDR: $(round(error_paired, digits=4))")

# 4. Print summary statistics
write_both(summary_file, "")
write_both(summary_file, "### Summary Statistics")
write_both(summary_file, "")
write_both(summary_file, "- Total precursors analyzed: $(nrow(prec_results))")
write_both(summary_file, "- Total entrapment pairs: $(length(unique(skipmissing(library_precursors.entrap_pair_id))))")

# Show EFDR at specific q-value thresholds
write_both(summary_file, "")
write_both(summary_file, "#### Global EFDR at Specific Q-value Thresholds")
write_both(summary_file, "")
for threshold in [0.001, 0.01, 0.05]
    idx = findfirst(prec_results[sortperm(prec_results.global_qval), :global_qval] .>= threshold)
    if !isnothing(idx)
        sorted_df = prec_results[sortperm(prec_results.global_qval), :]
        combined_efdr = sorted_df[idx, :global_prob_combined_efdr]
        paired_efdr = sorted_df[idx, :global_prob_paired_efdr]
        write_both(summary_file, "At q-value threshold $threshold:")
        write_both(summary_file, "- Combined EFDR: $(round(combined_efdr, digits=4))")
        write_both(summary_file, "- Paired EFDR: $(round(paired_efdr, digits=4))")
        write_both(summary_file, "")
    end
end

# 5. Create comparison plot for local scores
write_both(summary_file, "### Local Score Analysis")
write_both(summary_file, "")
write_both(summary_file, "Creating EFDR comparison plots for local scores...")
p_local = plot_efdr_comparison(prec_results, :prec_prob, :qval;
                              title="Local Score EFDR Comparison: Yeast Entrapment Analysis",
                              xlims=(0, 0.015),
                              ylims=(0, 0.015))
display(p_local)

# Save the local plot
savefig(p_local, joinpath(output_dir, "efdr_comparison_prec_prob.png"))
savefig(p_local, joinpath(output_dir, "efdr_comparison_prec_prob.pdf"))
write_both(summary_file, "")
write_both(summary_file, "Plots saved:")
write_both(summary_file, "- ![Local EFDR Comparison](efdr_comparison_prec_prob.png)")

# 6. Analyze local EFDR performance at standard thresholds
write_both(summary_file, "")
write_both(summary_file, "#### Local EFDR Performance at Standard Thresholds")
write_both(summary_file, "")
comparison_local = compare_efdr_methods(prec_results, :qval, :prec_prob, 
                                       library_precursors;
                                       thresholds=[0.001, 0.01, 0.05])

# Format comparison table for markdown
write_both(summary_file, "| Threshold | Q-val IDs | Actual FDR | Combined IDs | Combined EFDR | Paired IDs | Paired EFDR |")
write_both(summary_file, "|-----------|-----------|------------|--------------|---------------|------------|-------------|")
for row in eachrow(comparison_local)
    write_both(summary_file, @sprintf("| %.3f | %d | %.4f | %d | %.4f | %d | %.4f |",
            row.threshold, row.qval_n, row.qval_actual_fdr,
            row.combined_n, row.combined_efdr,
            row.paired_n, row.paired_efdr))
end

# 7. Calculate calibration error for local scores
write_both(summary_file, "")
write_both(summary_file, "#### Local Calibration Error")
write_both(summary_file, "")
cal_combined_local, error_combined_local = calculate_efdr_calibration_error(
    prec_results, :qval, :prec_prob_combined_efdr, library_precursors)
cal_paired_local, error_paired_local = calculate_efdr_calibration_error(
    prec_results, :qval, :prec_prob_paired_efdr, library_precursors)

write_both(summary_file, "Mean local calibration error:")
write_both(summary_file, "- Combined EFDR: $(round(error_combined_local, digits=4))")
write_both(summary_file, "- Paired EFDR: $(round(error_paired_local, digits=4))")

# Show local EFDR at specific q-value thresholds
write_both(summary_file, "")
write_both(summary_file, "#### Local EFDR at Specific Q-value Thresholds")
write_both(summary_file, "")
for threshold in [0.001, 0.01, 0.05]
    idx = findfirst(prec_results[sortperm(prec_results.qval), :qval] .>= threshold)
    if !isnothing(idx)
        sorted_df = prec_results[sortperm(prec_results.qval), :]
        combined_efdr = sorted_df[idx, :prec_prob_combined_efdr]
        paired_efdr = sorted_df[idx, :prec_prob_paired_efdr]
        write_both(summary_file, "At q-value threshold $threshold:")
        write_both(summary_file, "- Combined EFDR: $(round(combined_efdr, digits=4))")
        write_both(summary_file, "- Paired EFDR: $(round(paired_efdr, digits=4))")
        write_both(summary_file, "")
    end
end

# 8. Create combined plot showing all four EFDR curves
write_both(summary_file, "### Combined Analysis")
write_both(summary_file, "")
write_both(summary_file, "Creating combined plot with all EFDR curves...")
p_all = plot_efdr_vs_qval(prec_results, :global_qval,
                          [:global_prob_combined_efdr, :global_prob_paired_efdr],
                          title="Global vs Local EFDR Comparison",
                          labels=["Global Combined", "Global Paired"],
                          colors=[:blue, :red],
                          xlims=(0, 0.015),
                          ylims=(0, 0.015))

# Add local EFDR curves to the same plot
plot!(p_all, prec_results[sortperm(prec_results.qval), :qval],
      prec_results[sortperm(prec_results.qval), :prec_prob_combined_efdr],
      label="Local Combined", color=:darkblue, linestyle=:dash, linewidth=2)
plot!(p_all, prec_results[sortperm(prec_results.qval), :qval],
      prec_results[sortperm(prec_results.qval), :prec_prob_paired_efdr],
      label="Local Paired", color=:darkred, linestyle=:dash, linewidth=2)

display(p_all)
savefig(p_all, joinpath(output_dir, "efdr_comparison_all.png"))
savefig(p_all, joinpath(output_dir, "efdr_comparison_all.pdf"))

write_both(summary_file, "")
write_both(summary_file, "Plots saved:")
write_both(summary_file, "- ![Combined EFDR Comparison](efdr_comparison_all.png)")

write_both(summary_file, "")
write_both(summary_file, "---")
write_both(summary_file, "")
write_both(summary_file, "All analysis files saved to: `$output_dir`")

# Save calibration data to separate CSV files
CSV.write(joinpath(output_dir, "calibration_global_combined.csv"), cal_combined)
CSV.write(joinpath(output_dir, "calibration_global_paired.csv"), cal_paired)
CSV.write(joinpath(output_dir, "calibration_local_combined.csv"), cal_combined_local)
CSV.write(joinpath(output_dir, "calibration_local_paired.csv"), cal_paired_local)

# Save comparison tables to CSV
CSV.write(joinpath(output_dir, "efdr_comparison_global.csv"), comparison)
CSV.write(joinpath(output_dir, "efdr_comparison_local.csv"), comparison_local)

close(summary_file)

println("\n=== Analysis Complete ===")
println("All results saved to: $output_dir")
println("- Main summary: efdr_analysis_summary.md")
println("- Plots: *.png and *.pdf files")
println("- Data tables: *.csv files")