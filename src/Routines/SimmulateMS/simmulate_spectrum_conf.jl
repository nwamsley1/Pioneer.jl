using Distributions, Plots, Statistics

"""
    simulate_spectrum_with_confidence(
        spectral_lib, 
        prec_idx; 
        n_ions=1000, 
        n_simulations=100, 
        confidence_level=0.95
    )

Simulate multiple spectra for a precursor and calculate empirical confidence bands.

Parameters:
- spectral_lib: FragmentIndexLibrary containing the fragment data
- prec_idx: Precursor index/ID to simulate
- n_ions: Total number of ions in each simulation (default: 1000)
- n_simulations: Number of simulations to run (default: 100)
- confidence_level: Confidence level for bands (default: 0.95)

Returns:
- A NamedTuple with simulation results and statistics
"""
function simulate_spectrum_with_confidence(
    spectral_lib::FragmentIndexLibrary,
    prec_idx::Integer;
    n_ions::Int = 1000,
    n_simulations::Int = 100,
    confidence_level::Float64 = 0.95
)
    # Extract the fragment lookup table from the spectral library
    lookup = getFragmentLookupTable(spectral_lib)
    
    # Get the fragments for this precursor
    frag_range = getPrecFragRange(lookup, prec_idx)
    fragments = getFragments(lookup)[frag_range]
    
    # Skip if no fragments
    if isempty(fragments)
        @warn "No fragments found for precursor $prec_idx"
        return nothing
    end
    
    # Extract true intensities and m/z values
    intensities = Float32[]
    for frag in fragments
        if typeof(frag) <: DetailedFrag
            push!(intensities, Float32(frag.intensity))
        else
            spline_data = getSplineData(lookup)
            push!(intensities, getIntensity(frag, spline_data))
        end
    end
    
    # Normalize intensities to create probability distribution
    total_intensity = sum(intensities)
    if total_intensity == 0
        #@warn "Total intensity is zero for precursor $prec_idx"
        return nothing
    end
    
    probs = intensities ./ total_intensity
    mz_values = [getMz(frag) for frag in fragments]
    
    # Matrix to store simulation results: n_simulations × n_fragments
    n_fragments = length(fragments)
    simulation_results = zeros(Int, n_simulations, n_fragments)
    
    # Run simulations
    for i in 1:n_simulations
        # Create multinomial distribution and sample
        dist = Multinomial(n_ions, probs)
        counts = rand(dist)
        simulation_results[i, :] = counts
    end
    
    # Calculate statistics
    alpha = 1 - confidence_level
    mean_counts = mean(simulation_results, dims=1)[:]
    lower_bounds = [quantile(simulation_results[:, i], alpha/2) for i in 1:n_fragments]
    upper_bounds = [quantile(simulation_results[:, i], 1-alpha/2) for i in 1:n_fragments]
    
    # Return results
    return (
        mz_values = mz_values,
        true_intensities = intensities,
        normalized_intensities = probs,
        mean_counts = mean_counts,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        all_simulations = simulation_results
    )
end

"""
    plot_spectrum_with_confidence(sim_results; title="Spectrum with Confidence Bands")

Plot a simulated spectrum with empirical confidence bands.

Parameters:
- sim_results: Output from simulate_spectrum_with_confidence
- title: Title for the plot (default: "Spectrum with Confidence Bands")
- sort_by_mz: Whether to sort by m/z values (default: true)

Returns:
- A Plots.jl plot object
"""
function plot_spectrum_with_confidence(
    sim_results; 
    title="Spectrum with Confidence Bands",
    sort_by_mz=true
)
    if isnothing(sim_results)
        @warn "No simulation results to plot"
        return plot()
    end
    
    # Extract data
    mz_values = sim_results.mz_values
    mean_counts = sim_results.mean_counts
    lower_bounds = sim_results.lower_bounds
    upper_bounds = sim_results.upper_bounds
    
    # Sort by m/z if requested
    if sort_by_mz
        idx = sortperm(mz_values)
        mz_values = mz_values[idx]
        mean_counts = mean_counts[idx]
        lower_bounds = lower_bounds[idx]
        upper_bounds = upper_bounds[idx]
    end
    
    # Create plot
    p = plot(
        xlabel="m/z", 
        ylabel="Ion Count",
        title=title,
        legend=false
    )
    
    # Add confidence bands as ribbons
    for i in 1:length(mz_values)
        plot!(
            p,
            [mz_values[i], mz_values[i]],
            [lower_bounds[i], upper_bounds[i]],
            color=:blue,
            alpha=0.3,
            linewidth=2
        )
    end
    
    # Add mean points
    scatter!(
        p,
        mz_values,
        mean_counts,
        markersize=4,
        markerstrokewidth=0,
        color=:blue
    )
    
    # Add stem lines for means
    for i in 1:length(mz_values)
        plot!(
            p,
            [mz_values[i], mz_values[i]],
            [0, mean_counts[i]],
            color=:black,
            alpha=0.2,
            linewidth=1
        )
    end
    
    return p
end

"""
    plot_empirical_distribution(sim_results, fragment_idx; bins=20)
    
Plot the empirical distribution of counts for a specific fragment across simulations.

Parameters:
- sim_results: Output from simulate_spectrum_with_confidence
- fragment_idx: Index of the fragment to analyze
- bins: Number of histogram bins (default: 20)

Returns:
- A Plots.jl plot object
"""
function plot_empirical_distribution(
    sim_results, 
    fragment_idx::Int; 
    bins=20
)
    if isnothing(sim_results)
        @warn "No simulation results to plot"
        return plot()
    end
    
    # Extract data for the specific fragment
    mz = sim_results.mz_values[fragment_idx]
    counts = sim_results.all_simulations[:, fragment_idx]
    
    # Calculate the expected distribution (Binomial)
    n_ions = sum(sim_results.all_simulations[1, :])  # Total ions per simulation
    prob = sim_results.normalized_intensities[fragment_idx]  # Probability
    
    # Create histogram of empirical distribution
    plot_obj = histogram(
        counts, 
        bins=bins,
        normalize=:pdf,
        alpha=0.6,
        label="Empirical",
        title="Distribution for Fragment at m/z = $(round(mz, digits=2))",
        xlabel="Ion Count",
        ylabel="Probability Density"
    )
    
    # Overlay theoretical binomial distribution
    x_range = range(minimum(counts), maximum(counts), length=100)
    binom_dist = Binomial(n_ions, prob)
    binom_pdf = [pdf(binom_dist, round(Int, x)) for x in x_range]
    
    plot!(
        plot_obj,
        x_range,
        binom_pdf,
        line=:solid,
        linewidth=2,
        color=:red,
        label="Theoretical Binomial"
    )
    
    # Add vertical line for the mean
    vline!(
        plot_obj,
        [mean(counts)],
        color=:black,
        linestyle=:dash,
        label="Mean"
    )
    
    return plot_obj
end

using KernelDensity

"""
    plot_spectrum_with_distributions(sim_results; title="Spectrum with Intensity Distributions")

Create a spectrum plot where each peak is represented by a vertical band showing
the distribution of intensities across simulations.

Parameters:
- sim_results: Output from simulate_spectrum_with_confidence
- title: Plot title
- band_width: Width of the vertical bands (default: 2.0)
- density_points: Number of points to use for density estimation (default: 50)
- sort_by_mz: Whether to sort by m/z values (default: true)

Returns:
- A Plots.jl plot object
"""
function plot_spectrum_with_distributions(
    sim_results; 
    title="Spectrum with Intensity Distributions",
    band_width=2.0,
    density_points=50,
    sort_by_mz=true
)
    if isnothing(sim_results)
        @warn "No simulation results to plot"
        return plot()
    end
    
    # Extract data
    mz_values = sim_results.mz_values
    all_simulations = sim_results.all_simulations
    
    # Sort by m/z if requested
    if sort_by_mz
        idx = sortperm(mz_values)
        mz_values = mz_values[idx]
        all_simulations = all_simulations[:, idx]
    end
    
    # Create base plot
    p = plot(
        xlabel="m/z", 
        ylabel="Ion Count",
        title=title,
        legend=false
    )
    
    # For each fragment, create a vertical density band
    for i in 1:length(mz_values)
        mz = mz_values[i]
        counts = all_simulations[:, i]
        
        # Skip fragments with very low counts
        if maximum(counts) < 5
            continue
        end
        
        # Estimate the density of counts
        try
            # Use kernel density estimation
            k = kde(counts)
            
            # Normalize density to max height
            density_y = range(0, maximum(counts), length=density_points)
            density_x = zeros(length(density_y))
            
            # Create scaled density values at each height
            for j in 1:length(density_y)
                # Scale density to create width at each height
                density_at_y = pdf(k, density_y[j])
                max_density = maximum(pdf(k, counts))
                # Scale width proportionally to density
                density_x[j] = density_at_y / max_density * band_width
            end
            
            # Create left and right boundaries for the density shape
            left_x = mz .- density_x
            right_x = mz .+ density_x
            
            # Fill the density shape
            for j in 1:length(density_y)-1
                # Plot a polygon segment for each vertical slice
                plot!(
                    p,
                    [left_x[j], left_x[j+1], right_x[j+1], right_x[j], left_x[j]],
                    [density_y[j], density_y[j+1], density_y[j+1], density_y[j], density_y[j]],
                    seriestype=:shape,
                    color=:blue,
                    alpha=1.0,
                    linewidth=0
                )
            end
            
            # Add stem line for the mean
            mean_count = mean(counts)
            plot!(
                p,
                [mz, mz],
                [0, mean_count],
                color=:black,
                alpha=0.4,
                linewidth=1
            )
            
            # Add mean point
            scatter!(
                p,
                [mz],
                [mean_count],
                markersize=4,
                markerstrokewidth=0,
                color=:blue
            )
        catch e
            # If density estimation fails, fall back to simple confidence interval
            println("Density estimation failed for m/z $mz: $e")
            mean_count = mean(counts)
            lower = quantile(counts, 0.025)
            upper = quantile(counts, 0.975)
            
            # Plot confidence interval line
            plot!(
                p,
                [mz, mz],
                [lower, upper],
                color=:blue,
                alpha=0.4,
                linewidth=2
            )
            
            # Add mean point
            scatter!(
                p,
                [mz],
                [mean_count],
                markersize=4,
                markerstrokewidth=0,
                color=:blue
            )
        end
    end
    
    return p
end