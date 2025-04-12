using Distributions, Random

"""
    SimulatedPrecursor

Structure to hold information about a simulated precursor.

Fields:
- precursor_id: Identifier from the spectral library
- mz: Mass-to-charge ratio
- charge: Charge state
- abundance: Simulated abundance
- retention_time: Simulated retention time
- true_retention_time: Original retention time from library
- sequence: Peptide sequence
"""
struct SimulatedPrecursor
    precursor_id::UInt32
    mz::Float32
    charge::UInt8
    abundance::Float64
    retention_time::Float64
    true_retention_time::Float32
    sequence::String
end

struct SimulatedSpectrum{T<:AbstractFloat}
    mz_array::Vector{T}
    intensity_array::Vector{T}
    precursor_idxs::Vector{UInt32}
    rt::T
    centerMz::T
    isolationWidthMz::T
    ms_order::UInt8
end

"""
    MSSimulation

Container for all simulation data.

Fields:
- precursors: Vector of simulated precursors
- spectral_library: Reference to the original library
- parameters: Dictionary of simulation parameters
"""
struct MSSimulation
    precursors::Vector{SimulatedPrecursor}
    spectral_library::FragmentIndexLibrary
    parameters::Dict{String, Any}

    # Internal constructor that ensures precursors are sorted by retention time
    function MSSimulation(
        precursors::Vector{SimulatedPrecursor},
        spectral_library::FragmentIndexLibrary,
        parameters::Dict{String, <:Any}
    )
        # Sort precursors by retention time
        sorted_precursors = sort(precursors, by = p -> p.retention_time)
        
        # Return new instance with sorted precursors
        new(sorted_precursors, spectral_library, parameters)
    end

        
end

# First, properly import the function
import Base: issorted

# Then, extend it with your new method
issorted(x::MSSimulation) = issorted(x.precursors, by=p->p.retention_time)

"""
    plot_spectrum(
        spectrum::SimulatedSpectrum{T};
        plot_type=:stem,
        intensity_scale=:linear,
        top_n=nothing,
        min_intensity_pct=0.0,
        xlim=nothing,
        title=nothing,
        color=:blue,
        alpha=0.7,
        label_peaks=false,
        label_n_peaks=5
    ) where {T<:AbstractFloat}

Plot a mass spectrum from a SimulatedSpectrum object.

Parameters:
- spectrum: The SimulatedSpectrum to plot
- plot_type: :stem (default), :bar, or :scatter
- intensity_scale: :linear (default) or :log
- top_n: Only plot the top N most intense peaks (default: all)
- min_intensity_pct: Only plot peaks above this percentage of max intensity
- xlim: Tuple with (min, max) m/z range to display
- title: Custom title (default: auto-generated based on spectrum)
- color: Color for the peaks
- alpha: Transparency for peaks
- label_peaks: Whether to label the most intense peaks with their m/z values
- label_n_peaks: Number of peaks to label if label_peaks is true

Returns:
- A Plots.jl plot object
"""
function plot_spectrum(
    spectrum::SimulatedSpectrum{T};
    plot_type=:stem,
    intensity_scale=:linear,
    top_n=nothing,
    min_intensity_pct=0.0,
    xlim=nothing,
    title=nothing,
    color=:blue,
    alpha=0.7,
    label_peaks=false,
    label_n_peaks=5
) where {T<:AbstractFloat}
    
    # Check if spectrum is empty
    if isempty(spectrum.mz_array)
        return plot(
            xlabel="m/z",
            ylabel="Intensity",
            title="Empty Spectrum",
            legend=false
        )
    end
    
    # Filter peaks based on intensity threshold
    max_intensity = maximum(spectrum.intensity_array)
    min_intensity = max_intensity * (min_intensity_pct / 100.0)
    
    # Create indices mask for filtering
    mask = spectrum.intensity_array .>= min_intensity
    
    # Apply filtering
    mz_values = spectrum.mz_array[mask]
    intensities = spectrum.intensity_array[mask]
    
    # Further filter to top N peaks if specified
    if !isnothing(top_n) && top_n < length(mz_values)
        idx = sortperm(intensities, rev=true)[1:top_n]
        mz_values = mz_values[idx]
        intensities = intensities[idx]
        # Re-sort by m/z for display
        sort_idx = sortperm(mz_values)
        mz_values = mz_values[sort_idx]
        intensities = intensities[sort_idx]
    end
    
    # Apply log scale if requested
    if intensity_scale == :log
        intensities = log10.(intensities .+ 1.0)  # Add 1 to avoid log(0)
        y_label = "log₁₀(Intensity + 1)"
    else
        y_label = "Intensity"
    end
    
    # Generate automatic title if none provided
    if isnothing(title)
        ms_level = spectrum.ms_order
        if ms_level == 0x01
            title = "MS1 Spectrum at RT=$(spectrum.rt)"
        else
            title = "MS2 Spectrum at RT=$(spectrum.rt) for m/z=$(spectrum.centerMz)"
        end
    end
    
    # Create the plot with the requested plot type
    p = plot(
        xlabel="m/z",
        ylabel=y_label,
        title=title,
        legend=false
    )
    
    if plot_type == :stem
        # Stem plot (lines from baseline to points)
        for i in 1:length(mz_values)
            plot!(
                p,
                [mz_values[i], mz_values[i]],
                [0, intensities[i]],
                color=color,
                alpha=alpha,
                linewidth=1
            )
        end
        # Add points at the top
        scatter!(
            p,
            mz_values,
            intensities,
            markersize=3,
            markerstrokewidth=0,
            color=color
        )
    elseif plot_type == :bar
        # Bar plot
        bar!(
            p,
            mz_values,
            intensities,
            bar_width=(maximum(mz_values) - minimum(mz_values))/length(mz_values)/2,
            color=color,
            alpha=alpha
        )
    else  # Default to scatter
        # Scatter plot
        scatter!(
            p,
            mz_values,
            intensities,
            markersize=4,
            markerstrokewidth=0,
            color=color,
            alpha=alpha
        )
    end
    
    # Label the top N peaks if requested
    if label_peaks && !isempty(mz_values)
        # Get indices of top peaks to label
        top_indices = sortperm(intensities, rev=true)[1:min(label_n_peaks, length(intensities))]
        
        for idx in top_indices
            annotate!(
                p,
                mz_values[idx],
                intensities[idx] * 1.05,  # Slightly above the peak
                text("$(round(mz_values[idx], digits=1))", 8, :center)
            )
        end
    end
    
    # Set x-axis limits if provided
    if !isnothing(xlim)
        xlims!(p, xlim)
    end
    
    return p
end

"""
    ScanDefinition{T<:AbstractFloat}

Defines a single scan within an MS cycle.

Fields:
- ms_order: The MS level (MS1() or MS2())
- center_mz: Center m/z for isolation (missing for MS1)
- isolation_width: Width of isolation window (missing for MS1)
- duration_ms: Duration of the scan in milliseconds
"""
struct ScanDefinition{T<:AbstractFloat}
    ms_order::MsOrder
    center_mz::Union{T, Missing}
    isolation_width::Union{T, Missing}
    duration_ms::T
end

# Constructor for MS1 scans
function ScanDefinition(::MS1, duration_ms::T) where {T<:AbstractFloat}
    ScanDefinition{T}(MS1(), missing, missing, duration_ms)
end

# Constructor for MS2 scans
function ScanDefinition(::MS2, center_mz::T, isolation_width::T, duration_ms::T) where {T<:AbstractFloat}
    ScanDefinition{T}(MS2(), center_mz, isolation_width, duration_ms)
end

"""
    SimulatedMSRun{T<:AbstractFloat}

Container for a simulated MS run.

Fields:
- rt_points: Vector of retention time points
- spectra: Vector of simulated spectra at each RT point
- cycle_definition: Definition of the MS cycle used
- parameters: Additional parameters used for the simulation
"""
struct SimulatedMSRun{T<:AbstractFloat}
    rt_points::Vector{T}
    spectra::Vector{SimulatedSpectrum{T}}
    cycle_definition::Vector{ScanDefinition{T}}
    parameters::Dict{String, Any}
end