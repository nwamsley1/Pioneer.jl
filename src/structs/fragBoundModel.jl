struct FragBoundModel
    low_mass::ImmutablePolynomial{Float32}
    high_mass::ImmutablePolynomial{Float32}
end

function (fbm::FragBoundModel)(x::AbstractFloat)
    return fbm.low_mass(x), fbm.high_mass(x)
end

# Accessor functions
function get_bounds(model::FragBoundModel, prec_mz::AbstractFloat)
    min_mz = model.low_mass_poly(prec_mz)
    max_mz = model.high_mass_poly(prec_mz)
    return (min_mz, max_mz)
end
