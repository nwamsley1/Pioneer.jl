SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, 30.0);


include("Routines/SimmulateMS/simmulate_spectrum.jl")
# Select a precursor index to simulate
prec_idx = 1  # Replace with your desired precursor index

# Simulate spectrum for a precursor
spectrum = simulate_spectrum(SPEC_LIB, 1, n_ions=1000)

# Plot the spectrum
p = plot_spectrum(spectrum, title="Simulated Spectrum for Precursor 1")
display(p)



include("Routines/SimmulateMS/simmulate_spectrum_conf.jl")
# Run multiple simulations for a precursor
sim_results = simulate_spectrum_with_confidence(
    SPEC_LIB,
    1,  # Precursor index
    n_ions=20,
    n_simulations=10000,
    confidence_level=0.95
)


# Create the distribution spectrum plot
p_dist = plot_spectrum_with_distributions(
    sim_results,
    title="MS/MS Spectrum with Intensity Distributions",
    band_width=1.5  # Adjust this to control the width of bands
)
display(p_dist)

# Plot the spectrum with confidence bands
p1 = plot_spectrum_with_confidence(
    sim_results,
    title="Spectrum with 95% Confidence Bands for Precursor 1"
)
display(p1)

# Analyze the distribution of a specific fragment (e.g., the most intense one)
most_intense_idx = argmax(sim_results.mean_counts)
p2 = plot_empirical_distribution(sim_results, most_intense_idx)
display(p2)