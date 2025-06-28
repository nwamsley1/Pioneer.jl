# Example/test script for EFDR plotting
using Pkg
Pkg.activate(".")
using DataFrames, Test, Random
using Plots

# Include necessary files
include("entrapment_helper_funcs.jl")
include("efdr_funcs.jl") 
include("efdr_plotting.jl")

# Create more realistic test data
function create_realistic_test_data(n=1000)
    Random.seed!(42)
    
    # Create library with mix of targets and entrapments
    n_pairs = div(n, 3)
    entrapment_groups = Int[]
    pair_ids = UInt32[]
    
    for i in 1:n_pairs
        # Target
        push!(entrapment_groups, 0)
        push!(pair_ids, UInt32(i))
        # Entrapment from group 1
        push!(entrapment_groups, 1)
        push!(pair_ids, UInt32(i))
        # Sometimes add entrapment from group 2
        if rand() > 0.5 && length(entrapment_groups) < n
            push!(entrapment_groups, 2)
            push!(pair_ids, UInt32(i))
        end
    end
    
    # Pad to desired length if needed
    while length(entrapment_groups) < n
        push!(entrapment_groups, 0)
        push!(pair_ids, UInt32(0))  # Unpaired
    end
    
    library_precursors = DataFrame(
        entrapment_group_id = entrapment_groups[1:n],
        entrap_pair_id = pair_ids[1:n]
    )
    
    # Create scoring results with realistic distributions
    # Targets generally have better scores than entrapments
    scores = Float64[]
    qvals = Float64[]
    
    for i in 1:n
        if entrapment_groups[i] == 0
            # Targets: better scores
            score = 0.5 + 0.5 * rand()  # 0.5 to 1.0
            qval = rand() * 0.1  # 0 to 0.1
        else
            # Entrapments: worse scores on average
            score = 0.2 + 0.6 * rand()  # 0.2 to 0.8
            qval = 0.01 + rand() * 0.15  # 0.01 to 0.16
        end
        push!(scores, score)
        push!(qvals, qval)
    end
    
    prec_results = DataFrame(
        precursor_idx = 1:n,
        global_prob = scores,
        global_qval = qvals,
        prec_prob = scores .* (0.9 .+ 0.1 .* rand(n)),  # Slightly different scores
        qval = qvals .* (0.8 .+ 0.4 .* rand(n))  # Slightly different q-values
    )
    
    return prec_results, library_precursors
end

# Generate test data
println("Generating test data...")
prec_results, library_precursors = create_realistic_test_data(2000)

# Calculate EFDR values
println("Calculating EFDR values...")
add_efdr_columns!(prec_results, library_precursors)

# Create plots
println("Creating plots...")

# 1. Basic comparison plot for global scores
p1 = plot_efdr_comparison(prec_results, :global_prob, :global_qval;
                         title="Global Score EFDR Comparison",
                         xlims=(0, 0.15),
                         ylims=(0, 0.3))
display(p1)

# 2. Comparison with custom styling
p2 = plot_efdr_comparison(prec_results, :prec_prob, :qval;
                         title="Precursor Score EFDR Comparison",
                         colors=[:darkblue, :darkred],
                         xlims=(0, 0.2),
                         ylims=(0, 0.4),
                         legend=:topleft)
display(p2)

# 3. Combined plot for both score types
p3 = plot_multiple_efdr_comparisons(prec_results,
                                   [(:global_prob, :global_qval), (:prec_prob, :qval)])
display(p3)

# 4. Manual plot with all EFDR curves on one axis
p4 = plot_efdr_vs_qval(prec_results, :global_qval,
                      [:global_prob_combined_efdr, :global_prob_paired_efdr],
                      title="All EFDR Methods vs Global Q-value",
                      labels=["Combined EFDR", "Paired EFDR"],
                      colors=[:blue, :red],
                      xlims=(0, 0.1),
                      ylims=(0, 0.25))
display(p4)

# Save plots to files (uncomment to use)
# save_efdr_plots(prec_results, "efdr_plots/")

println("\nPlotting examples completed!")

# Print some statistics
println("\nEFDR Statistics:")
for col in [:global_prob_combined_efdr, :global_prob_paired_efdr]
    println("$col:")
    println("  Mean: $(round(mean(prec_results[!, col]), digits=4))")
    println("  Max:  $(round(maximum(prec_results[!, col]), digits=4))")
    println("  At FDR 0.01: $(round(prec_results[sortperm(prec_results.global_qval), col][findfirst(prec_results[sortperm(prec_results.global_qval), :global_qval] .>= 0.01)], digits=4))")
end