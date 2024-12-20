
"""
Base configuration parameters for library building
"""
struct LibraryBuildParams
    fasta_params::Dict{String, Any}
    nce_params::Dict{String, Any}
    library_params::Dict{String, Any}
    fasta_paths::Vector{String}
    fasta_names::Vector{String}
    out_dir::String
    lib_name::String
    predict_fragments::Bool
end

"""
Fragment boundary model for m/z range prediction
"""
struct FragBoundModel
    low_mass_poly::ImmutablePolynomial  # Predicts min m/z
    high_mass_poly::ImmutablePolynomial # Predicts max m/z
end

# Accessor functions
function get_bounds(model::FragBoundModel, prec_mz::AbstractFloat)
    min_mz = model.low_mass_poly(prec_mz)
    max_mz = model.high_mass_poly(prec_mz)
    return (min_mz, max_mz)
end