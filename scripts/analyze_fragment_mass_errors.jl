#!/usr/bin/env julia

"""
    analyze_fragment_mass_errors.jl

Script to analyze fragment mass error distributions from Pioneer.jl FirstPassSearch output.

This script investigates how mass error distributions vary with:
- Fragment intensity (absolute and relative to max fragment)
- Retention time
- Fragment type (b, y, p ions)
- Fragment charge state
- Precursor characteristics

Usage:
    julia analyze_fragment_mass_errors.jl <path_to_fragment_errors.arrow>

The script generates:
- Mass error distribution plots
- Intensity vs mass error scatter plots
- Retention time vs mass error plots
- Statistical summaries and diagnostic information
"""

using Arrow
using DataFrames
using Plots
using StatsPlots
using Statistics
using Printf
using Dates

# Ion type mapping for better visualization
const ION_TYPE_NAMES = Dict(
    0 => "unknown",
    1 => "b",
    2 => "y",
    3 => "p"
)

function load_fragment_errors(file_path::String)
    """Load fragment mass error data from Arrow file."""
    if !isfile(file_path)
        error("File not found: $file_path")
    end

    df = DataFrame(Arrow.Table(file_path))
    println("Loaded $(nrow(df)) fragment mass errors from $file_path")

    # Add derived columns for analysis
    df.ion_type_name = [get(ION_TYPE_NAMES, ion_type, "unknown") for ion_type in df.ion_type]
    df.relative_intensity = df.fragment_intensity ./ df.max_fragment_intensity
    df.log10_intensity = log10.(max.(df.fragment_intensity, 1.0))
    df.abs_ppm_error = abs.(df.ppm_error)

    return df
end

function basic_statistics(df::DataFrame)
    """Print basic statistics about the fragment mass errors."""
    println("\n" * "="^60)
    println("FRAGMENT MASS ERROR STATISTICS")
    println("="^60)

    n_total = nrow(df)
    n_files = length(unique(df.ms_file_idx))
    n_precursors = length(unique(df.precursor_idx))
    n_scans = length(unique(df.scan_idx))

    println("Dataset Overview:")
    println("  Total fragments: $(n_total)")
    println("  MS files: $(n_files)")
    println("  Unique precursors: $(n_precursors)")
    println("  Unique scans: $(n_scans)")

    println("\nMass Error Distribution:")
    println("  Mean PPM error: $(@sprintf("%.3f", mean(df.ppm_error)))")
    println("  Median PPM error: $(@sprintf("%.3f", median(df.ppm_error)))")
    println("  Std dev PPM error: $(@sprintf("%.3f", std(df.ppm_error)))")
    println("  MAD PPM error: $(@sprintf("%.3f", mad(df.ppm_error, normalize=true)))")
    println("  95th percentile |error|: $(@sprintf("%.3f", quantile(df.abs_ppm_error, 0.95)))")

    println("\nRetention Time Range:")
    println("  Min RT: $(@sprintf("%.2f", minimum(df.retention_time))) min")
    println("  Max RT: $(@sprintf("%.2f", maximum(df.retention_time))) min")

    println("\nIntensity Range:")
    println("  Min intensity: $(@sprintf("%.0f", minimum(df.fragment_intensity)))")
    println("  Max intensity: $(@sprintf("%.0f", maximum(df.fragment_intensity)))")
    println("  Median intensity: $(@sprintf("%.0f", median(df.fragment_intensity)))")

    println("\nIon Type Distribution:")
    for (ion_type, count) in sort(collect(countmap(df.ion_type_name)))
        pct = 100 * count / n_total
        println("  $ion_type ions: $count ($(@sprintf("%.1f", pct))%)")
    end

    println("\nCharge State Distribution:")
    for (charge, count) in sort(collect(countmap(df.fragment_charge)))
        pct = 100 * count / n_total
        println("  Charge $charge: $count ($(@sprintf("%.1f", pct))%)")
    end
end

function plot_mass_error_distributions(df::DataFrame, output_dir::String=".")
    """Create comprehensive mass error distribution plots."""

    println("\nGenerating mass error distribution plots...")

    # Overall mass error histogram
    p1 = histogram(df.ppm_error,
                   bins=100,
                   title="Fragment Mass Error Distribution",
                   xlabel="PPM Error",
                   ylabel="Frequency",
                   alpha=0.7,
                   label="All fragments")
    vline!([mean(df.ppm_error)], color=:red, linewidth=2, label="Mean")
    vline!([median(df.ppm_error)], color=:orange, linewidth=2, label="Median")

    # Mass error by ion type
    p2 = @df df violin(:ion_type_name, :ppm_error,
                       title="Mass Error by Ion Type",
                       xlabel="Ion Type",
                       ylabel="PPM Error",
                       alpha=0.7)

    # Mass error vs retention time
    p3 = scatter(df.retention_time, df.ppm_error,
                 alpha=0.3,
                 markersize=2,
                 title="Mass Error vs Retention Time",
                 xlabel="Retention Time (min)",
                 ylabel="PPM Error",
                 label="")

    # Mass error vs intensity (log scale)
    p4 = scatter(df.log10_intensity, df.ppm_error,
                 alpha=0.3,
                 markersize=2,
                 title="Mass Error vs Fragment Intensity",
                 xlabel="Log10(Intensity)",
                 ylabel="PPM Error",
                 label="")

    # Combine plots
    combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))

    # Save plot
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_file = joinpath(output_dir, "fragment_mass_error_analysis_$(timestamp).png")
    savefig(combined_plot, output_file)
    println("Saved mass error distribution plots to: $output_file")

    return combined_plot
end

function analyze_intensity_dependence(df::DataFrame, output_dir::String=".")
    """Analyze mass error dependence on fragment intensity."""

    println("\nAnalyzing intensity dependence...")

    # Create intensity bins for analysis
    n_bins = 10
    intensity_quantiles = quantile(df.log10_intensity, range(0, 1, length=n_bins+1))
    df.intensity_bin = cut(df.log10_intensity, intensity_quantiles, labels=false)

    # Calculate statistics by intensity bin
    intensity_stats = combine(groupby(df, :intensity_bin)) do sdf
        DataFrame(
            mean_intensity = mean(sdf.log10_intensity),
            mean_error = mean(sdf.ppm_error),
            median_error = median(sdf.ppm_error),
            mad_error = mad(sdf.ppm_error, normalize=true),
            std_error = std(sdf.ppm_error),
            n_fragments = nrow(sdf)
        )
    end

    println("\nMass Error Statistics by Intensity Bin:")
    println("Bin | Mean Log10(Int) | Mean Error | Median Error | MAD Error | N Fragments")
    println("-"^75)
    for row in eachrow(intensity_stats)
        println("$(row.intensity_bin) | $(@sprintf("%12.2f", row.mean_intensity)) | " *
                "$(@sprintf("%10.3f", row.mean_error)) | " *
                "$(@sprintf("%12.3f", row.median_error)) | " *
                "$(@sprintf("%9.3f", row.mad_error)) | " *
                "$(@sprintf("%11d", row.n_fragments))")
    end

    # Plot error vs intensity with trend
    p1 = scatter(df.log10_intensity, df.abs_ppm_error,
                 alpha=0.2,
                 markersize=1,
                 title="Absolute Mass Error vs Fragment Intensity",
                 xlabel="Log10(Intensity)",
                 ylabel="Absolute PPM Error",
                 label="")

    # Add binned medians
    scatter!(intensity_stats.mean_intensity,
             abs.(intensity_stats.median_error),
             color=:red,
             markersize=4,
             label="Binned Medians")

    # Relative intensity analysis
    p2 = scatter(df.relative_intensity, df.abs_ppm_error,
                 alpha=0.2,
                 markersize=1,
                 title="Absolute Mass Error vs Relative Intensity",
                 xlabel="Relative Intensity (fraction of max)",
                 ylabel="Absolute PPM Error",
                 label="")

    # Combine plots
    intensity_plot = plot(p1, p2, layout=(1,2), size=(1200, 400))

    # Save plot
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_file = joinpath(output_dir, "fragment_intensity_analysis_$(timestamp).png")
    savefig(intensity_plot, output_file)
    println("Saved intensity analysis plots to: $output_file")

    return intensity_stats, intensity_plot
end

function analyze_rt_dependence(df::DataFrame, output_dir::String=".")
    """Analyze mass error dependence on retention time."""

    println("\nAnalyzing retention time dependence...")

    # Create RT bins for analysis
    n_bins = 20
    rt_quantiles = quantile(df.retention_time, range(0, 1, length=n_bins+1))
    df.rt_bin = cut(df.retention_time, rt_quantiles, labels=false)

    # Calculate statistics by RT bin
    rt_stats = combine(groupby(df, :rt_bin)) do sdf
        DataFrame(
            mean_rt = mean(sdf.retention_time),
            mean_error = mean(sdf.ppm_error),
            median_error = median(sdf.ppm_error),
            mad_error = mad(sdf.ppm_error, normalize=true),
            std_error = std(sdf.ppm_error),
            n_fragments = nrow(sdf)
        )
    end

    println("\nMass Error Statistics by Retention Time:")
    println("Bin | Mean RT (min) | Mean Error | Median Error | MAD Error | N Fragments")
    println("-"^75)
    for row in eachrow(rt_stats)
        println("$(row.rt_bin) | $(@sprintf("%13.2f", row.mean_rt)) | " *
                "$(@sprintf("%10.3f", row.mean_error)) | " *
                "$(@sprintf("%12.3f", row.median_error)) | " *
                "$(@sprintf("%9.3f", row.mad_error)) | " *
                "$(@sprintf("%11d", row.n_fragments))")
    end

    # Plot error vs RT with trend
    p1 = scatter(df.retention_time, df.ppm_error,
                 alpha=0.2,
                 markersize=1,
                 title="Mass Error vs Retention Time",
                 xlabel="Retention Time (min)",
                 ylabel="PPM Error",
                 label="")

    # Add binned medians
    scatter!(rt_stats.mean_rt,
             rt_stats.median_error,
             color=:red,
             markersize=4,
             label="Binned Medians")

    # Plot absolute error vs RT
    p2 = scatter(df.retention_time, df.abs_ppm_error,
                 alpha=0.2,
                 markersize=1,
                 title="Absolute Mass Error vs Retention Time",
                 xlabel="Retention Time (min)",
                 ylabel="Absolute PPM Error",
                 label="")

    # Add binned MADs
    scatter!(rt_stats.mean_rt,
             rt_stats.mad_error,
             color=:red,
             markersize=4,
             label="Binned MADs")

    # Combine plots
    rt_plot = plot(p1, p2, layout=(1,2), size=(1200, 400))

    # Save plot
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_file = joinpath(output_dir, "fragment_rt_analysis_$(timestamp).png")
    savefig(rt_plot, output_file)
    println("Saved RT analysis plots to: $output_file")

    return rt_stats, rt_plot
end

function main()
    """Main analysis function."""

    if length(ARGS) != 1
        println("Usage: julia analyze_fragment_mass_errors.jl <path_to_fragment_errors.arrow>")
        exit(1)
    end

    file_path = ARGS[1]
    output_dir = dirname(file_path)

    println("Fragment Mass Error Analysis")
    println("Input file: $file_path")
    println("Output directory: $output_dir")

    # Load data
    df = load_fragment_errors(file_path)

    # Basic statistics
    basic_statistics(df)

    # Generate plots and analyses
    plot_mass_error_distributions(df, output_dir)
    intensity_stats, _ = analyze_intensity_dependence(df, output_dir)
    rt_stats, _ = analyze_rt_dependence(df, output_dir)

    # Summary insights
    println("\n" * "="^60)
    println("ANALYSIS INSIGHTS")
    println("="^60)

    println("\n1. Overall Mass Accuracy:")
    println("   Mean systematic error: $(@sprintf("%.3f", mean(df.ppm_error))) ppm")
    println("   Random error (MAD): $(@sprintf("%.3f", mad(df.ppm_error, normalize=true))) ppm")

    intensity_range = maximum(intensity_stats.mad_error) - minimum(intensity_stats.mad_error)
    rt_range = maximum(rt_stats.mad_error) - minimum(rt_stats.mad_error)

    println("\n2. Intensity Dependence:")
    println("   MAD error range across intensity bins: $(@sprintf("%.3f", intensity_range)) ppm")
    if intensity_range > 1.0
        println("   ⚠️  Significant intensity-dependent mass error detected!")
    else
        println("   ✓  Mass error appears relatively independent of intensity")
    end

    println("\n3. Retention Time Dependence:")
    println("   MAD error range across RT bins: $(@sprintf("%.3f", rt_range)) ppm")
    if rt_range > 1.0
        println("   ⚠️  Significant RT-dependent mass error detected!")
    else
        println("   ✓  Mass error appears relatively independent of retention time")
    end

    # Ion type comparison
    ion_stats = combine(groupby(df, :ion_type_name)) do sdf
        DataFrame(
            median_error = median(sdf.ppm_error),
            mad_error = mad(sdf.ppm_error, normalize=true)
        )
    end

    println("\n4. Ion Type Comparison:")
    for row in eachrow(ion_stats)
        println("   $(row.ion_type_name) ions: median = $(@sprintf("%.3f", row.median_error)) ppm, " *
                "MAD = $(@sprintf("%.3f", row.mad_error)) ppm")
    end

    println("\n✓ Analysis complete! Check output directory for plots.")
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end