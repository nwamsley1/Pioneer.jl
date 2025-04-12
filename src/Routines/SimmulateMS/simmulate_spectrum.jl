using Distributions, Plots
abstract type MsOrder end
struct MS1 <: MsOrder end
struct MS2 <: MsOrder end
# Option 1: Using [] syntax
import Base: getindex
getindex(::MS1) = UInt8(1)
getindex(::MS2) = UInt8(2)

"""
    addIsotopes!(
        intensity_array::Vector{T},
        mz_array::Vector{T},
        frag::DetailedFrag{T},
        iso_splines::IsotopeSplineModel,
        n_isotopes::Int,
        ms_order::MsOrder
    ) where {T<:AbstractFloat}

Add isotope peaks for a fragment to mz and intensity arrays.

Parameters:
- intensity_array: Vector to store intensities
- mz_array: Vector to store m/z values
- frag: Fragment ion
- iso_splines: Isotope spline model for calculating relative abundances
- n_isotopes: Number of isotope peaks to generate
- ms_order: MS1 or MS2 level
"""
function addIsotopes!(
    intensity_array::Vector{T},
    mz_array::Vector{T},
    base_intensity::Float16,
    mono_mz::U,
    charge::UInt8,
    nsulfur::UInt8,
    iso_splines::IsotopeSplineModel,
    n_isotopes::Int,
) where {T,U<:AbstractFloat}
    
    for iso_idx in range(1, n_isotopes)
        # Calculate isotope m/z
        isotope_mz = mono_mz + iso_idx * (NEUTRON/charge)
        
        # Calculate isotope intensity using splines
        isotope_intensity = iso_splines(
            min(nsulfur, 5),  # Cap sulfur count at 5 for model limits
            iso_idx,
            mono_mz * charge  # Mass (not m/z)
        ) * base_intensity
        
        # Add to arrays
        push!(mz_array, isotope_mz)
        push!(intensity_array, isotope_intensity)
    end
end

"""
    simulate_spectrum(spectral_lib, prec_idx; n_ions=1000)

Simulate a spectrum for a precursor by sampling fragment ions according to their 
relative intensities using a multinomial distribution.

Parameters:
- spectral_lib: FragmentIndexLibrary containing the fragment data
- prec_idx: Precursor index/ID to simulate
- n_ions: Total number of ions to simulate (default: 1000)

Returns:
- A NamedTuple with mz_values, intensities, and the raw count data
"""
function simulate_spectrum(
    spectral_lib::FragmentIndexLibrary,
    iso_splines::IsotopeSplineModel,
    prec_idx::Integer,
    n_isotopes::Int,
    n_ions::Int,
    ::MS2
)
    # Extract the fragment lookup table from the spectral library
    lookup = getFragmentLookupTable(spectral_lib)
    
    # Get the fragments for this precursor
    frag_range = getPrecFragRange(lookup, prec_idx)
    fragments = getFragments(lookup)[frag_range]
    
    # Skip if no fragments
    if isempty(fragments)
        @warn "No fragments found for precursor $prec_idx"
        return (mz_values=Float32[], intensities=Float32[], counts=Int[])
    end

    # Extract intensities
    intensity_array = Float32[]
    mz_array = Float32[]
    for frag in fragments
        intensity = getIntensity(frag)
        if isnan(intensity)
            intensity = zero(Float16)
        end
        addIsotopes!(
            intensity_array,
            mz_array,
            intensity,
            getMz(frag),
            getFragCharge(frag),
            getSulfurCount(frag),
            iso_splines,
            n_isotopes
        )
    end

    # Normalize intensities to create probability distribution
    total_intensity = sum( intensity_array )
    if total_intensity == 0
        #@warn "Total intensity is zero for precursor $prec_idx"
        return (mz_values=Float32[],  intensities=Float32[], counts=Int[])
    end
    
    fragment_probs = intensity_array ./ total_intensity
    # Create multinomial distribution and sample
    dist = Multinomial(n_ions, fragment_probs)
    fragment_counts = rand(dist)
    
    # Return simulated spectrum
    return (
        mz_values = mz_array[fragment_counts .> 0],
        intensities = fragment_counts[fragment_counts .> 0],
        counts = counts
    )
end


"""
    simulate_spectrum(spectral_lib, prec_idx; n_ions=1000)

Simulate a spectrum for a precursor by sampling fragment ions according to their 
relative intensities using a multinomial distribution.

Parameters:
- spectral_lib: FragmentIndexLibrary containing the fragment data
- prec_idx: Precursor index/ID to simulate
- n_ions: Total number of ions to simulate (default: 1000)

Returns:
- A NamedTuple with mz_values, intensities, and the raw count data
"""
function simulate_spectrum(
    spectral_lib::FragmentIndexLibrary,
    iso_splines::IsotopeSplineModel,
    prec_idx::Integer,
    n_isotopes::Int,
    n_ions::Int,
    ::MS1
)
    precursors = getPrecursors(spectral_lib)
    # Extract intensities
    intensity_array = Float32[]
    mz_array = Float32[]
    addIsotopes!(
        intensity_array,
        mz_array,
        one(Float16),
        getMz(precursors)[prec_idx],
        getCharge(precursors)[prec_idx],
        getSulfurCount(precursors)[prec_idx],
        iso_splines,
        n_isotopes
    )

    # Normalize intensities to create probability distribution
    total_intensity = sum( intensity_array )
    if total_intensity == 0
        #@warn "Total intensity is zero for precursor $prec_idx"
        return (mz_values=Float32[],  intensities=Float32[], counts=Int[])
    end
    
    ion_probs = intensity_array ./ total_intensity
    
    # Create multinomial distribution and sample
    dist = Multinomial(n_ions,  ion_probs)
    ion_counts = rand(dist)
    
    # Return simulated spectrum
    return (
        mz_values = mz_array[ion_counts .> 0],
        intensities = ion_counts[ion_counts .> 0],
        counts = counts
    )
end

#=
function simulate_spectrum(
    spectral_lib::FragmentIndexLibrary,
    prec_idx::Integer,
    ms_order::MS1,
    n_ions::Int = 1000
)
    # Extract the fragment lookup table from the spectral library
    lookup = getFragmentLookupTable(spectral_lib)
    
    # Get the fragments for this precursor
    frag_range = getPrecFragRange(lookup, prec_idx)
    fragments = getFragments(lookup)[frag_range]
    
    # Skip if no fragments
    if isempty(fragments)
        @warn "No fragments found for precursor $prec_idx"
        return (mz_values=Float32[], intensities=Float32[], counts=Int[])
    end
    
    # Extract intensities
    intensities = Float32[]
    for frag in fragments
        
        # For DetailedFrag, intensity is stored directly
        if typeof(frag) <: DetailedFrag
            push!(intensities, Float32(frag.intensity))
        else
            # For SplineDetailedFrag, we need to calculate using splines
            spline_data = getSplineData(lookup)
            push!(intensities, getIntensity(frag, spline_data))
        end
    end
    
    # Normalize intensities to create probability distribution
    total_intensity = sum(intensities)
    if total_intensity == 0
        @warn "Total intensity is zero for precursor $prec_idx"
        return (mz_values=Float32[], intensities=Float32[], counts=Int[])
    end
    
    probs = intensities ./ total_intensity
    
    # Create multinomial distribution and sample
    dist = Multinomial(n_ions, probs)
    counts = rand(dist)
    
    # Extract m/z values
    mz_values = [getMz(frag) for frag in fragments]
    
    # Return simulated spectrum
    return (
        mz_values = mz_values,
        intensities = intensities,
        counts = counts
    )
end
=#

"""
    plot_spectrum(spectrum; title="Simulated Spectrum")

Plot a simulated spectrum with proper formatting.

Parameters:
- spectrum: The output from simulate_spectrum
- title: Title for the plot (default: "Simulated Spectrum")

Returns:
- A Plots.jl plot object
"""
function plot_spectrum(spectrum; title="Simulated Spectrum")
    p = scatter(
        spectrum.mz_values, 
        spectrum.counts,
        xlabel="m/z", 
        ylabel="Ion Count",
        title=title,
        legend=false,
        markersize=4,
        markerstrokewidth=0,
        alpha=0.7
    )
    
    # Add stem lines for better visibility
    for i in 1:length(spectrum.mz_values)
        plot!(
            p,
            [spectrum.mz_values[i], spectrum.mz_values[i]],
            [0, spectrum.counts[i]],
            color=:black,
            alpha=0.3,
            linewidth=1
        )
    end
    
    return p
end

"""
    merge_peaks(mz_array::Vector{T}, intensity_array::Vector{T}, ppm_tolerance::T) where {T<:AbstractFloat}

Merge peaks that are within a specified ppm tolerance.

When peaks are merged:
- Intensities are summed
- The m/z value becomes a weighted average based on intensity

Parameters:
- mz_array: Array of m/z values
- intensity_array: Array of intensities
- ppm_tolerance: Parts-per-million tolerance for merging

Returns:
- Tuple of (merged_mz_array, merged_intensity_array)
"""
function merge_peaks(mz_array::Vector{T}, intensity_array::Vector{T}, ppm_tolerance::T) where {T<:AbstractFloat}
    if isempty(mz_array)
        return mz_array, intensity_array
    end
    
    # Sort by m/z for efficiency
    idx = sortperm(mz_array)
    sorted_mz = mz_array[idx]
    sorted_intensities = intensity_array[idx]
    
    # Arrays for merged peaks
    merged_mz = T[]
    merged_intensities = T[]
    
    # Start with the first peak
    current_mz = sorted_mz[1]
    current_intensity = sorted_intensities[1]
    
    for i in 2:length(sorted_mz)
        # Check if current peak is within tolerance of the accumulated peak
        ppm_diff = abs(sorted_mz[i] - current_mz) / current_mz * T(1_000_000)
        
        if ppm_diff <= ppm_tolerance
            # Merge peaks using weighted average for m/z
            total_intensity = current_intensity + sorted_intensities[i]
            current_mz = (current_mz * current_intensity + 
                         sorted_mz[i] * sorted_intensities[i]) / total_intensity
            current_intensity = total_intensity
        else
            # Record the accumulated peak and start a new one
            push!(merged_mz, current_mz)
            push!(merged_intensities, current_intensity)
            current_mz = sorted_mz[i]
            current_intensity = sorted_intensities[i]
        end
    end
    
    # Don't forget the last accumulated peak
    push!(merged_mz, current_mz)
    push!(merged_intensities, current_intensity)
    
    return merged_mz, merged_intensities
end

"""
    simulate_precursor_spectra(
        spectral_lib::FragmentIndexLibrary,
        precursors::Vector{SimulatedPrecursor},
        ms_order::UInt8;
        n_ions_per_precursor::Int = 1000
    ) where {T<:AbstractFloat}

Simulate spectra for a list of precursors.

Parameters:
- spectral_lib: Fragment library containing the spectral data
- precursors: Vector of precursors to simulate
- ms_order: MS level (1 for MS1, 2 for MS2, etc.)
- center_mz: Center of m/z isolation window (for MS2+)
- isolation_width: Width of m/z isolation window (for MS2+)
- min_intensity: Minimum intensity threshold to include a peak
- n_ions_per_precursor: Number of ions to simulate per precursor

Returns:
- Tuple of (mz_array, intensity_array) for all simulated peaks
"""
function simulate_precursor_spectra(
    spectral_lib::FragmentIndexLibrary,
    iso_splines::IsotopeSplineModel,
    precursors::Vector{SimulatedPrecursor},
    ms_order::MsOrder;
    n_ions_per_precursor::Int = 1000,
    ppm_tolerance::T = T(5.0)
) where {T<:AbstractFloat}
    # Arrays to hold all fragment data
    all_mz_values = T[]
    all_intensities = T[]
        
    # Simulate fragments for each matching precursor
    for precursor in precursors
        # Get precursor ID as integer
        precursor_id = Int(precursor.precursor_id)
        
        # Simulate the fragment spectrum
        simulated_spectrum = simulate_spectrum(
            spectral_lib,
            iso_splines,
            precursor_id,
            5,#n_isotopes,
            n_ions_per_precursor,
            ms_order
        )
        # Skip if empty
        if isempty(simulated_spectrum.mz_values)
            continue
        end
        
        # Scale intensity by precursor abundance
        scaled_intensities = T.(simulated_spectrum.intensities .* 
                                precursor.abundance)
        # Add fragments with intensity above threshold
        for i in 1:length(simulated_spectrum.mz_values)
            push!(all_mz_values, T(simulated_spectrum.mz_values[i]))
            push!(all_intensities, scaled_intensities[i])
            # Merge peaks within ppm tolerance
            all_mz_values, all_intensities = merge_peaks(all_mz_values, all_intensities, ppm_tolerance)
        end
    end
    
    return all_mz_values, all_intensities
end



"""
    filter_precursors(
        precursors::Vector{SimulatedPrecursor},
        rt::T,
        rt_tolerance::T,
        center_mz::T,
        isolation_width::T,
        ms_order::MsOrder
    ) where {T<:AbstractFloat}

Filter precursors based on retention time and optionally m/z windows,
using the MsOrder trait for type-specific behavior.

Returns:
- Vector of filtered precursors
"""
function filter_precursors(
    precursors::Vector{SimulatedPrecursor},
    rt::T,
    rt_tolerance::T,
    center_mz::T,
    isolation_width::T,
    ms_order::MsOrder
) where {T<:AbstractFloat}
    # Define the retention time window
    rt_min = rt - rt_tolerance/2
    rt_max = rt + rt_tolerance/2
    
    # First filter by retention time (common for all MS orders)
    rt_filtered = filter(
        p -> rt_min <= p.retention_time <= rt_max,
        precursors
    )
    
    # Return early if no precursors match RT window
    isempty(rt_filtered) && return rt_filtered
    
    # Use the MsOrder trait to determine if additional filtering is needed
    return _filter_by_ms_order(rt_filtered, center_mz, isolation_width, ms_order)
end

# MS1 doesn't need m/z filtering
function _filter_by_ms_order(
    precursors::Vector{SimulatedPrecursor},
    center_mz::T,
    isolation_width::T,
    ::MS1
) where {T<:AbstractFloat}
    return precursors
end

# MS2 needs additional m/z filtering
function _filter_by_ms_order(
    precursors::Vector{SimulatedPrecursor},
    center_mz::T,
    isolation_width::T,
    ::MS2
) where {T<:AbstractFloat}
    mz_min = center_mz - isolation_width/2
    mz_max = center_mz + isolation_width/2
    
    return filter(
        p -> mz_min <= p.mz <= mz_max,
        precursors
    )
end

"""
    simulate_combined_spectrum(
        ms_sim::MSSimulation,
        rt::T,
        rt_tolerance::T,
        center_mz::T,
        isolation_width::T;
        ms_order::UInt8 = 0x01,
        min_intensity::T = T(0.0),
        n_ions_per_precursor::Int = 1000,
        ppm_tolerance::T = T(5.0)
    ) where {T<:AbstractFloat}

Simulate a mass spectrum by combining fragments from precursors within specified windows.
Peaks within the specified ppm tolerance are merged.

Parameters:
- ms_sim: MSSimulation object containing precursors
- rt: Center of retention time window
- rt_tolerance: Width of retention time window (full width)
- center_mz: Center of m/z isolation window
- isolation_width: Width of m/z isolation window (full width)
- ms_order: MS level (1 for MS1, 2 for MS2, etc.)
- min_intensity: Minimum intensity threshold to include a peak
- n_ions_per_precursor: Number of ions to simulate per precursor
- ppm_tolerance: Parts-per-million tolerance for merging peaks

Returns:
- SimulatedSpectrum containing combined peaks from all matching precursors
"""
function simulate_combined_spectrum(
    ms_sim::MSSimulation,
    iso_splines::IsotopeSplineModel,
    rt::T,
    rt_tolerance::T,
    center_mz::T,
    isolation_width::T,
    ms_order::MsOrder;
    n_ions_per_precursor::Int = 1000,
    ppm_tolerance::T = T(5.0)
) where {T<:AbstractFloat}
    
    # Filter precursors using helper function
    matching_precursors = filter_precursors(
        ms_sim.precursors,
        rt,
        rt_tolerance,
        center_mz,
        isolation_width,
        ms_order
    )


    # No matching precursors, return empty spectrum
    if isempty(matching_precursors)
        return SimulatedSpectrum{T}(
            Vector{T}(),
            Vector{T}(),
            Vector{UInt32}(),
            rt,
            center_mz,
            isolation_width,
            ms_order[]  # Get UInt8 value from trait
        )
    end

    # Simulate spectra for matching precursors
    all_mz_values, all_intensities = simulate_precursor_spectra(
        ms_sim.spectral_library,
        iso_splines,
        matching_precursors,
        ms_order;
        n_ions_per_precursor = n_ions_per_precursor,
        ppm_tolerance = ppm_tolerance
    )
    
    # If no peaks were generated, return an empty spectrum
    if isempty(all_mz_values)
        return SimulatedSpectrum{T}(
            Vector{T}(),
            Vector{T}(),
            map(x->x.precursor_id, matching_precursors),
            rt,
            center_mz,
            isolation_width,
            ms_order[]
        )
    end

    # Create and return the simulated spectrum
    return SimulatedSpectrum{T}(
        all_mz_values,
        all_intensities,
        map(x->x.precursor_id, matching_precursors),
        rt,
        center_mz,
        isolation_width,
        ms_order[]
    )
end