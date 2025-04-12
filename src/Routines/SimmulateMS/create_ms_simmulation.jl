"""
    create_ms_simulation(
        spec_lib::FragmentIndexLibrary;
        precursor_proportion::Float64 = 0.1,
        abundance_mean::Float64 = 0.0,
        abundance_std::Float64 = 1.0,
        rt_std::Float64 = 0.5,
        seed::Union{Int, Nothing} = nothing
    )

Create a simulation by randomly selecting precursors and assigning properties.

Parameters:
- spec_lib: Spectral library containing precursor information
- precursor_proportion: Proportion of precursors to select (default: 0.1)
- abundance_mean: Mean of log-normal distribution for abundances (default: 0.0)
- abundance_std: Standard deviation of log-normal distribution (default: 1.0)
- rt_std: Standard deviation for retention time variation (default: 0.5)
- seed: Random seed for reproducibility

Returns:
- MSSimulation object containing the selected precursors
"""
function create_ms_simulation(
    spec_lib::FragmentIndexLibrary;
    precursor_proportion::Float64 = 0.1,
    abundance_mean::Float64 = 0.0,
    abundance_std::Float64 = 1.0,
    rt_std::Float64 = 0.5,
    seed::Union{Int, Nothing} = nothing
)
    # Set random seed if provided
    if !isnothing(seed)
        Random.seed!(seed)
    end
    
    # Get precursor information from the library
    precursors = getPrecursors(spec_lib)
    n_total_precursors = length(precursors)
    n_select = round(Int, n_total_precursors * precursor_proportion)
    
    # Randomly select precursors
    selected_indices = randperm(n_total_precursors)[1:n_select]
    
    # Get necessary data from the library
    mzs = getMz(precursors)
    charges = getCharge(precursors)
    irts = getIrt(precursors)
    sequences = getSequence(precursors)
    
    # Create log-normal distribution for abundances
    abundance_dist = LogNormal(abundance_mean, abundance_std)
    
    # Create normal distribution for retention time noise
    rt_noise_dist = Normal(0.0, rt_std)
    
    # Create simulated precursors
    simulated_precursors = Vector{SimulatedPrecursor}(undef, n_select)
    
    for (i, idx) in enumerate(selected_indices)
        # Generate abundance from log-normal distribution
        abundance = rand(abundance_dist)
        
        # Generate retention time with Gaussian noise
        true_rt = irts[idx]
        rt_noise = rand(rt_noise_dist)
        simulated_rt = true_rt + rt_noise
        
        # Create the simulated precursor
        simulated_precursors[i] = SimulatedPrecursor(
            UInt32(idx),
            mzs[idx],
            charges[idx],
            abundance,
            simulated_rt,
            true_rt,
            sequences[idx]
        )
    end
    
    # Sort by retention time for easier access later
    sort!(simulated_precursors, by = p -> p.retention_time)
    
    # Create and return the simulation object
    return MSSimulation(
        simulated_precursors,
        spec_lib,
        Dict(
            "precursor_proportion" => precursor_proportion,
            "abundance_mean" => abundance_mean,
            "abundance_std" => abundance_std,
            "rt_std" => rt_std,
            "n_selected_precursors" => n_select
        )
    )
end

"""
    simulate_spectrum_from_precursors(
        simulation::MSSimulation,
        rt::Float64;
        rt_window::Float64 = 0.5,
        n_ions_per_precursor::Int = 1000
    )

Simulate a mass spectrum at a specific retention time.

Parameters:
- simulation: MSSimulation object
- rt: Retention time at which to simulate
- rt_window: Width of retention time window to consider
- n_ions_per_precursor: Number of ions to simulate per precursor

Returns:
- A NamedTuple with precursor information and simulated fragment spectra
"""
function simulate_spectrum_from_precursors(
    simulation::MSSimulation,
    rt::Float64;
    rt_window::Float64 = 0.5,
    n_ions_per_precursor::Int = 1000
)
    # Find precursors within the retention time window
    rt_min = rt - rt_window/2
    rt_max = rt + rt_window/2
    
    active_precursors = filter(
        p -> rt_min <= p.retention_time <= rt_max,
        simulation.precursors
    )
    
    if isempty(active_precursors)
        @warn "No precursors found in RT window $(rt_min)-$(rt_max)"
        return (precursors = [], spectra = [])
    end
    
    # Generate spectra for each active precursor
    precursor_spectra = []
    
    for precursor in active_precursors
        # Simulate the fragment spectrum using our existing function
        spectrum = simulate_spectrum(
            simulation.spectral_library,
            Int(precursor.precursor_id),
            n_ions = n_ions_per_precursor
        )
        
        # Scale intensities by the precursor abundance
        scaled_spectrum = (
            mz_values = spectrum.mz_values,
            intensities = spectrum.intensities * precursor.abundance,
            counts = spectrum.counts
        )
        
        push!(precursor_spectra, scaled_spectrum)
    end
    
    return (
        precursors = active_precursors,
        spectra = precursor_spectra
    )
end

"""
    plot_simulated_precursors(simulation::MSSimulation; title="Simulated Precursors")

Plot the simulated precursors with retention time vs m/z, colored by abundance.

Parameters:
- simulation: MSSimulation object
- title: Title for the plot

Returns:
- A Plots.jl plot object
"""
function plot_simulated_precursors(simulation::MSSimulation; title="Simulated Precursors")
    precursors = simulation.precursors
    
    # Extract data
    rts = [p.retention_time for p in precursors]
    mzs = [p.mz for p in precursors]
    log_abundances = log10.([p.abundance for p in precursors])
    
    # Create scatter plot
    p = scatter(
        rts, 
        mzs,
        marker_z = log_abundances,
        xlabel = "Retention Time",
        ylabel = "m/z",
        title = title,
        colorbar_title = "Log10(Abundance)",
        markersize = 4,
        alpha = 0.7,
        legend = false
    )
    
    return p
end

"""
    plot_rt_profile(simulation::MSSimulation; bin_width=1.0, title="Retention Time Profile")

Plot the distribution of precursors across retention time.

Parameters:
- simulation: MSSimulation object
- bin_width: Width of retention time bins
- title: Title for the plot

Returns:
- A Plots.jl plot object
"""
function plot_rt_profile(
    simulation::MSSimulation; 
    bin_width=1.0, 
    title="Retention Time Profile"
)
    precursors = simulation.precursors
    rts = [p.retention_time for p in precursors]
    
    p = histogram(
        rts,
        bins = round(Int, (maximum(rts) - minimum(rts))/bin_width),
        xlabel = "Retention Time",
        ylabel = "Count",
        title = title,
        alpha = 0.7,
        legend = false
    )
    
    return p
end