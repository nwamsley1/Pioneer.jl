"""
Type alias for m/z to eV interpolation functions.
Uses `GriddedInterpolation` with linear interpolation and line extrapolation.
"""
const InterpolationTypeAlias = Interpolations.Extrapolation{
    Float32,
    1,
    Interpolations.GriddedInterpolation{
        Float32,
        1,
        Vector{Float32},
        Gridded{Linear{Throw{OnGrid}}},
        Tuple{Vector{Float32}},
    },
    Gridded{Linear{Throw{OnGrid}}},
    Line{Nothing},
}

const CHARGE_ADJUSTMENT_FACTORS = Float64[1, 0.9, 0.85, 0.8, 0.75]
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)
