"""
    get_fragment_bounds(
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
function get_fragment_bounds(
    center_mass::AbstractVector{<:Union{Missing, Float32}},
    isolation_width::AbstractVector{<:Union{Missing, Float32}},
    ms_order::AbstractVector{UInt8},
    low_frag_mass::AbstractVector{Float32},
    high_frag_mass::AbstractVector{Float32}
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

function get_fragment_bounds(
    auto_detect_frag_bounds::Bool,
    frag_bounds_detection_raw_file_path::String,
    default_frag_bounds::Tuple{Float32, Float32},
    default_precursor_bounds::Tuple{Float32, Float32})::@NamedTuple{frag_bounds::FragBoundModel, prec_mz_min::Float32, prec_mz_max::Float32}
    
    if auto_detect_frag_bounds
        if isfile(frag_bounds_detection_raw_file_path)
            try
            MS_TABLE = Arrow.Table(frag_bounds_detection_raw_file_path)
            frag_bounds, prec_mz_min, prec_mz_max =  get_fragment_bounds(
                MS_TABLE[:centerMz],
                MS_TABLE[:isolationWidthMz],
                MS_TABLE[:msOrder],
                MS_TABLE[:lowMz],
                MS_TABLE[:highMz]
            )
            prec_mz_min -= 1.0f0
            prec_mz_max += 1.0f0
            return (frag_bounds = frag_bounds, prec_mz_min = prec_mz_min, prec_mz_max = prec_mz_max)
            catch# e 
                #throw(e)
                @warn "failed to estimate fragbounds from the example raw file.
                Using default values. Frag bounds: $default_frag_bounds, precursor bounds: $default_precursor_bounds"
            end
        else
            @warn "Could not find example file $frag_bounds_detection_raw_file_path to determine frag bounds. 
            Using default values. Frag bounds: $default_frag_bounds, precursor bounds: $default_precursor_bounds"
        end
    end
    frag_bounds = FragBoundModel(
        ImmutablePolynomial(first(default_frag_bounds)),
        ImmutablePolynomial(last(default_frag_bounds)) 
    )
    return (frag_bounds = frag_bounds, 
            prec_mz_min = first(default_precursor_bounds), 
            prec_mz_max = last(default_precursor_bounds))
end