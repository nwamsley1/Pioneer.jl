using Plots
using DataFrames
using Statistics

"""
    plot_efdr_vs_qval(df::DataFrame, qval_col::Symbol, efdr_cols::Vector{Symbol};
                      title="Empirical FDR vs Q-value",
                      xlabel="Q-value",
                      ylabel="Empirical FDR",
                      labels=nothing,
                      colors=nothing,
                      xlims=(0, 0.1),
                      ylims=(0, 0.2),
                      legend=:bottomright,
                      diagonal=true)

Plot empirical FDR values against original q-values.

# Arguments
- `df`: DataFrame containing the data
- `qval_col`: Symbol of the q-value column to plot on x-axis
- `efdr_cols`: Vector of EFDR column symbols to plot
- `title`: Plot title
- `xlabel`: X-axis label
- `ylabel`: Y-axis label
- `labels`: Custom labels for each EFDR series (defaults to column names)
- `colors`: Custom colors for each series
- `xlims`: X-axis limits
- `ylims`: Y-axis limits
- `legend`: Legend position
- `diagonal`: Whether to add diagonal reference line
"""
function plot_efdr_vs_qval(df::DataFrame, qval_col::Symbol, efdr_cols::Vector{Symbol};
                          title="Empirical FDR vs Q-value",
                          xlabel="Q-value", 
                          ylabel="Empirical FDR",
                          labels=nothing,
                          colors=nothing,
                          xlims=(0, 0.015),
                          ylims=(0, 0.015),
                          legend=:bottomright,
                          diagonal=true)
    
    # Sort by q-value
    sorted_indices = sortperm(df[!, qval_col])
    sorted_df = df[sorted_indices, :]
    
    # Create plot
    p = plot(title=title, xlabel=xlabel, ylabel=ylabel, 
             xlims=xlims, ylims=ylims, legend=legend,
             size=(600, 500), dpi=300)
    
    # Add diagonal reference line if requested
    if diagonal
        max_val = min(xlims[2], ylims[2])
        plot!(p, [0, max_val], [0, max_val], 
              label="y=x", linestyle=:dash, color=:gray, alpha=0.5)
    end
    
    # Default colors if not provided
    if isnothing(colors)
        colors = [:blue, :red, :green, :orange, :purple]
    end
    
    # Plot each EFDR column
    for (i, efdr_col) in enumerate(efdr_cols)
        label = isnothing(labels) ? String(efdr_col) : labels[i]
        color = colors[mod1(i, length(colors))]
        
        plot!(p, sorted_df[!, qval_col], sorted_df[!, efdr_col],
              label=label, color=color, linewidth=2, alpha=0.8)
    end
    
    return p
end

"""
    plot_efdr_comparison(df::DataFrame, score_col::Symbol, qval_col::Symbol;
                        methods=[CombinedEFDR(), PairedEFDR()],
                        kwargs...)

Plot comparison of different EFDR methods for a given score/qval pair.
Automatically extracts the appropriate column names based on the score column and methods.

# Arguments
- `df`: DataFrame with EFDR columns already calculated
- `score_col`: The score column (e.g., :global_prob)
- `qval_col`: The q-value column (e.g., :global_qval)
- `methods`: Vector of EFDR methods to compare
- `kwargs`: Additional keyword arguments passed to plot_efdr_vs_qval
"""
function plot_efdr_comparison(df::DataFrame, score_col::Symbol, qval_col::Symbol;
                             method_types::Vector=[CombinedEFDR, PairedEFDR],
                             kwargs...)
    
    # Build EFDR column names based on score_col and method types
    efdr_cols = Symbol[]
    labels = String[]
    
    for method_type in method_types
        method_name = if method_type == CombinedEFDR
            "combined"
        elseif method_type == PairedEFDR
            "paired"
        else
            error("Unknown EFDR method type: $method_type")
        end
        
        efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
        
        # Check if column exists
        if !hasproperty(df, efdr_col)
            error("Column $efdr_col not found. Make sure to run add_efdr_columns! first.")
        end
        
        push!(efdr_cols, efdr_col)
        push!(labels, titlecase(method_name) * " EFDR")
    end
    
    # Default title if not provided
    default_kwargs = Dict(
        :title => "EFDR Comparison for $(String(score_col))",
        :labels => labels
    )
    
    # Merge with user kwargs (user kwargs take precedence)
    merged_kwargs = merge(default_kwargs, kwargs)
    
    return plot_efdr_vs_qval(df, qval_col, efdr_cols; merged_kwargs...)
end

"""
    plot_multiple_efdr_comparisons(df::DataFrame, score_qval_pairs::Vector{Tuple{Symbol,Symbol}};
                                  method_types=[CombinedEFDR, PairedEFDR],
                                  layout=nothing,
                                  kwargs...)

Create a grid of plots comparing EFDR methods for multiple score/qval pairs.

# Arguments
- `df`: DataFrame with EFDR columns already calculated
- `score_qval_pairs`: Vector of (score_col, qval_col) tuples
- `method_types`: EFDR method types to compare
- `layout`: Plot layout (defaults to automatic grid)
- `kwargs`: Additional keyword arguments passed to each subplot
"""
function plot_multiple_efdr_comparisons(df::DataFrame, 
                                      score_qval_pairs::Vector{Tuple{Symbol,Symbol}};
                                      method_types::Vector=[CombinedEFDR, PairedEFDR],
                                      layout=nothing,
                                      kwargs...)
    
    n_plots = length(score_qval_pairs)
    
    # Determine layout if not provided
    if isnothing(layout)
        n_cols = ceil(Int, sqrt(n_plots))
        n_rows = ceil(Int, n_plots / n_cols)
        layout = (n_rows, n_cols)
    end
    
    # Create subplots
    plots = []
    for (score_col, qval_col) in score_qval_pairs
        p = plot_efdr_comparison(df, score_col, qval_col; 
                               method_types=method_types, 
                               legend=(length(plots) == 0 ? :bottomright : false),
                               kwargs...)
        push!(plots, p)
    end
    
    # Combine plots
    return plot(plots..., layout=layout, size=(800 * layout[2], 600 * layout[1]))
end

"""
    save_efdr_plots(df::DataFrame, output_dir::String;
                   score_qval_pairs=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                   method_types=[CombinedEFDR, PairedEFDR],
                   formats=[:png, :pdf])

Save EFDR comparison plots to files.

# Arguments
- `df`: DataFrame with EFDR columns
- `output_dir`: Directory to save plots
- `score_qval_pairs`: Score/qval pairs to plot
- `method_types`: EFDR method types to compare
- `formats`: File formats to save (e.g., [:png, :pdf])
"""
function save_efdr_plots(df::DataFrame, output_dir::String;
                        score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                        method_types::Vector=[CombinedEFDR, PairedEFDR],
                        formats::Vector{Symbol}=[:png, :pdf])
    
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Save individual plots
    for (score_col, qval_col) in score_qval_pairs
        p = plot_efdr_comparison(df, score_col, qval_col; method_types=method_types)
        
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_$(score_col).$(format)")
            savefig(p, filename)
            println("Saved: $filename")
        end
    end
    
    # Save combined plot if multiple pairs
    if length(score_qval_pairs) > 1
        p_combined = plot_multiple_efdr_comparisons(df, score_qval_pairs; method_types=method_types)
        
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_all.$(format)")
            savefig(p_combined, filename)
            println("Saved: $filename")
        end
    end
end

# Export functions
export plot_efdr_vs_qval, plot_efdr_comparison, plot_multiple_efdr_comparisons, save_efdr_plots