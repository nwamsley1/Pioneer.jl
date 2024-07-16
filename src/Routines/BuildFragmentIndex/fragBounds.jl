
#MS_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/astral_test/20230324_OLEP08_200ng_30min_E10H50Y40_180K_2Th3p5ms_01.arrow")
struct FragBoundModel
    low_mass::ImmutablePolynomial{Float32}
    high_mass::ImmutablePolynomial{Float32}
end

function (fbm::FragBoundModel)(x::AbstractFloat)
    return fbm.low_mass(x), fbm.high_mass(x)
end

"""
    getFragBounds(
        center_mass::AbstractVector{Float32},
        isolation_width::AbstractVector{Float32},

        low_frag_mass::AbstractVector{Float32},
        high_frag_mass::AbstractVector{Float32}
    ) where {T<:AbstractFloat}

Given the isolation centers and widths of MS2 scans from a raw file and the corresponding
scan frange for these MS2 spectra, learns a linear model to predict the maximum and minimum 
size (m/z) fragment peak given the precursor m/z. 

Outputs: FragBoundModel
- Essentially a tuple of polynomials. One to predict the high fragment m/z given the precursor m/z
and the other to predict the low fragment m/z
"""
function getFragBounds(
    center_mass::AbstractVector{Union{Missing, Float32}},
    isolation_width::AbstractVector{Union{Missing, Float32}},
    ms_order::AbstractVector{Union{Missing, Int32}},
    low_frag_mass::AbstractVector{Union{Missing, Float32}},
    high_frag_mass::AbstractVector{Union{Missing, Float32}}
    )

    #Filter out non-ms2 scans
    ms2_indicator = ms_order .== 2
    center_mass = center_mass[ms2_indicator]
    isolation_width = isolation_width[ms2_indicator]
    low_frag_mass = low_frag_mass[ms2_indicator]
    high_frag_mass = high_frag_mass[ms2_indicator]

    w_center_low = hcat(ones(length(center_mass)), center_mass .-  isolation_width./2)
    w_center_high = hcat(ones(length(center_mass)),center_mass .+ isolation_width./2)
    return FragBoundModel(
        ImmutablePolynomial(w_center_high\low_frag_mass),
        ImmutablePolynomial(w_center_low\high_frag_mass)
    ), minimum(w_center_low[:,2]), maximum(w_center_high[:,2])
end
