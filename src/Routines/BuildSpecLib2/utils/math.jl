
# src/utils/math.jl
"""
Calculate parts per million difference between two masses.
"""
function get_ppm(mass1::T, mass2::T)::T where T <: AbstractFloat
    return abs(mass1 - mass2) / ((mass1 + mass2) / 2.0) * 1e6
end

"""
Get PPM tolerance window for a given mass.
"""
function get_ppm_window(mass::T, tolerance_ppm::T)::Tuple{T,T} where T <: AbstractFloat
    window = mass * tolerance_ppm / 1e6
    return (mass - window, mass + window)
end

"""
Calculate theoretical isotope distribution.
"""
function calculate_isotope_distribution(
    n_carbons::Int,
    n_sulfurs::Int,
    charge::Int;
    max_isotopes::Int=3
)::Vector{Float32}
    
    # Constants for natural isotope abundances
    C13_ABUNDANCE = 0.0107
    S34_ABUNDANCE = 0.0442
    
    distribution = zeros(Float32, max_isotopes)
    distribution[1] = 1.0  # monoisotopic peak
    
    # Calculate probability of each isotope peak
    for i in 2:max_isotopes
        # Carbon contribution
        c_prob = binomial(n_carbons, i-1) * (C13_ABUNDANCE^(i-1)) * ((1-C13_ABUNDANCE)^(n_carbons-(i-1)))
        
        # Sulfur contribution (if present)
        s_prob = 0.0
        if n_sulfurs > 0
            s_prob = binomial(n_sulfurs, i-1) * (S34_ABUNDANCE^(i-1)) * ((1-S34_ABUNDANCE)^(n_sulfurs-(i-1)))
        end
        
        distribution[i] = Float32((c_prob + s_prob) / charge)
    end
    
    return distribution
end