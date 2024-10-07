struct FragBoundModel
    low_mass::ImmutablePolynomial{Float32}
    high_mass::ImmutablePolynomial{Float32}
end

function (fbm::FragBoundModel)(x::AbstractFloat)
    return fbm.low_mass(x), fbm.high_mass(x)
end