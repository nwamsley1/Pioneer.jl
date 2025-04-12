#Need to load a .pion formmated spectral library and use this to simmulate data. 
#using Distributions

using Distributions 
# Example: 5 possible fragments with their probabilities
fragment_probs = [0.2, 0.15, 0.3, 0.25, 0.1]

# Create a multinomial distribution with 1000 trials
dist = Multinomial(1000, fragment_probs)

# Sample from the distribution
result = rand(dist)

SPEC_LIB_DIR = "/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/Astral_tag6_SageLibrary_aligned.pion"

# Add this before calling loadSpectralLibrary
using Pioneer 
import Pioneer: DetailedFrag, StandardFragmentLookup
# Add this constructor method
function StandardFragmentLookup(
    frags::Vector{Pioneer.DetailedFrag{T}}, 
    prec_frag_ranges::Vector{UInt64}
) where {T<:AbstractFloat}
    # Either map the fragments to the correct type
    converted_frags = map(f -> DetailedFrag{T}(
        f.prec_id,
        f.mz,
        f.intensity,
        f.ion_type,
        f.is_y,
        f.is_b,
        f.is_p,
        f.is_isotope,
        f.frag_charge,
        f.ion_position,
        f.prec_charge,
        f.rank,
        f.sulfur_count
    ), frags)
    
    # Then create a StandardFragmentLookup with the converted fragments]
    StandardFragmentLookup{T}(converted_frags, prec_frag_ranges)
end
SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, 30.0);


include("Routines/SimmulateMS/simmulate_spectrum.jl")
# Select a precursor index to simulate
prec_idx = 1  # Replace with your desired precursor index

# Simulate spectrum for a precursor
spectrum = simulate_spectrum(SPEC_LIB, 1, n_ions=1000)

include("Routines/SimmulateMS/structs.jl")
include("Routines/SimmulateMS/create_ms_simmulation.jl")
include("Routines/SimmulateMS/simmulate_spectrum.jl")
include("Routines/SimmulateMS/build_run.jl")
include("Routines/SimmulateMS/plot_run.jl")
# Create a simulation with 10% of the precursors
ms_sim = create_ms_simulation(
    SPEC_LIB,
    precursor_proportion = 1.0,
    abundance_mean = 10.0,
    abundance_std = 1.0,
    rt_std = 0.5,
    seed = 42
);

# Define the cycle based on the MS table data pattern
cycle_def = create_cycle_definition([
    # First MS1 scan
    (MS1(), 115.0),  
    
    # First set of MS2 scans
    (MS2(), 519.25f0, 38.5f0, 65.0f0),
    (MS2(), 552.25f0, 28.5f0, 65.0f0),
    (MS2(), 578.25f0, 24.5f0, 65.0f0),
    (MS2(), 600.75f0, 21.5f0, 65.0f0),
    (MS2(), 622.25f0, 22.5f0, 65.0f0),
    (MS2(), 641.25f0, 16.5f0, 65.0f0),
    (MS2(), 656.25f0, 14.5f0, 65.0f0),
    (MS2(), 668.75f0, 11.5f0, 65.0f0),
    (MS2(), 679.75f0, 11.5f0, 65.0f0),
    (MS2(), 692.75f0, 15.5f0, 65.0f0),
    
    # Second MS1 scan
    (MS1(), 115.0f0),
    
    # Second set of MS2 scans
    (MS2(), 708.25f0, 16.5f0, 65.0f0),
    (MS2(), 722.25f0, 12.5f0, 65.0f0),
    (MS2(), 735.75f0, 15.5f0, 65.0f0),
    (MS2(), 748.75f0, 11.5f0, 65.0f0),
    (MS2(), 762.75f0, 17.5f0, 65.0f0),
    (MS2(), 780.25f0, 18.5f0, 65.0f0),
    (MS2(), 802.75f0, 27.5f0, 65.0f0),
    (MS2(), 831.75f0, 31.5f0, 65.0f0),
    (MS2(), 870.25f0, 46.5f0, 65.0f0),
    (MS2(), 946.5f0, 107.0f0, 65.0f0),
    
    # Third MS1 scan to complete the cycle
    (MS1(), 115.0f0)
])

_iso_splines_= parseIsoXML(joinpath(@__DIR__,"../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml"));
   
# Simulate the entire run from RT 10 to 30 minutes
ms_run = simulate_ms_run(
    ms_sim,
    _iso_splines_,
    cycle_def,
    Float32(20.0),   # Start RT
    Float32(21.0),   # End RT
    n_cycles = 100,  # 100 cycles over the RT range
    rt_tolerance = Float32(0.1), #
    n_ions_per_precursor = 1000,
    ppm_tolerance = Float32(5.0) #merge peaks 
)

# Create an overview plot
overview_plot = plot_ms_run_overview(ms_run)
display(overview_plot)

# Basic plot of an MS2 spectrum
p1 = plot_spectrum(ms_run.spectra[300])
display(p1)

         # Simulate an MS2 spectrum with peak merging
ms2_spectrum = simulate_combined_spectrum(
    ms_sim,
    _iso_splines_,                  # Averagine splines 
    Float32(20.0),                # RT center
    Float32(0.1),                 # RT window
    Float32(720.4),               # Precursor m/z center
    Float32(10.0),                 # Isolation width
    MS2(),              # MS2 level
    ppm_tolerance = Float32(3.0)            # Merge peaks within 5 ppm
)

# Basic plot of an MS2 spectrum
p1 = plot_spectrum(ms2_spectrum)
display(p1)

# Log scale with top 20 peaks labeled
p2 = plot_spectrum(
    ms2_spectrum,
    intensity_scale=:log,
    label_peaks=true,
    label_n_peaks=20,
    color=:purple
)
display(p2)

# Visualize the MS2 spectrum
p = scatter(
    ms2_spectrum.mz_array,
    ms2_spectrum.intensity_array,
    xlabel = "m/z",
    ylabel = "Intensity",
    title = "Simulated MS2 Spectrum at RT=$(ms2_spectrum.rt) for m/z=$(ms2_spectrum.centerMz)",
    legend = false
)
display(p)
#histogram(log.([x.abundance for x in ms_sim.precursors]))
# Plot the selected precursors
precursor_plot = plot_simulated_precursors(ms_sim)
display(precursor_plot)

# Plot the retention time profile
rt_plot = plot_rt_profile(ms_sim, bin_width=0.5)
display(rt_plot)

# Simulate a spectrum at a specific retention time
rt = 20.0  # Choose a retention time with precursors
result = simulate_spectrum_from_precursors(
    ms_sim, 
    rt,
    rt_window = 1.0,
    n_ions_per_precursor = 1000
)

# Display information about precursors in this window
for (i, precursor) in enumerate(result.precursors)
    println("Precursor $i: m/z = $(precursor.mz), RT = $(precursor.retention_time), Abundance = $(precursor.abundance)")
    
    # Plot the spectrum for this precursor
    spectrum = result.spectra[i]
    p = scatter(
        spectrum.mz_values,
        spectrum.intensities,
        title = "Spectrum for $(precursor.sequence)",
        xlabel = "m/z",
        ylabel = "Intensity",
        legend = false
    )
    display(p)
end


using Plots
using Distributions

# Parameters for the Gaussian distribution
μ = 500    # mean (in milliseconds)
σ = 100    # standard deviation (in milliseconds)
dist = Normal(μ, σ)

# Range of x values (time in milliseconds)
x_range = 100:1:900  # from 100ms to 900ms

# First plot: Gaussian distribution curve
p1 = plot(x_range, x -> 1000*pdf(dist, x), 
    title = "Peptide Elution Example", 
    xlabel = "Time (ms)", 
    ylabel = "ions/ms",
    legend = false,
    lw = 2,
    color = :blue)



    # Define points x1 and x2 for shading
x1 = 300  # Example value
x2 = 350  # Example value

# Second plot: Shaded area under the curve from x1 to x2
p2 = plot(x_range, x -> 1000*pdf(dist, x), 
    title = "Peptide Elution Example", 
    xlabel = "Time (ms)", 
    ylabel = "ions/ms",
    legend = false,
    lw = 2,
    color = :blue)

# Add shaded area from x1 to x2
x_fill = x1:1:x2
plot!(p2, x_fill, x -> 1000*pdf(dist, x), 
    fill = (0, 0.5, :green),
    alpha = 0.3)

# Add vertical lines at x1 and x2
vline!(p2, [x1, x2], color = :red, linestyle = :dash, label = "")

# Display both plots side by side
plot(p1, p2, layout = (1, 2), size = (900, 400))

# To save the plots to a file
savefig("gaussian_plots.png")
