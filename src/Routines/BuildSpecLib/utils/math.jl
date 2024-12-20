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
