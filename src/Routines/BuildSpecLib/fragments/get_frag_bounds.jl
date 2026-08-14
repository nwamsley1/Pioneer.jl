# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

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

"""
Coefficients for a manual fragment-bound rule, as `(slope, intercept)` pairs
applied to the precursor m/z.

Expressed against the **low edge** of the isolation window, matching the
regression in the auto path above, which fits `high_mass` on `centre - width/2`.
A rule written against the window centre is wider by `slope × width/2` — about
14 m/z on a 14 Th window — so the distinction is not cosmetic.

`thermo_auto_documented` is Thermo's own rule, `last_mass = z × centre + 10`
with DIA's z = 2. `thermo_auto` is what two Thermo instruments were measured
doing (Astral and Exploris 480, both slope 2.04 to four decimals); its intercept
is method-specific — 24.1 on one, 50.0 on the other — so it is offered rather
than assumed. Neither is a large win: against a flat ceiling the measured rule
recovers 0.36% of intensity, and exactly nothing at charge 2, because 2 × centre
is very nearly the largest fragment a 2+ precursor at that centre can produce.

SCIEX and Bruker hold one fixed MS2 range for every window, which is what
`constant` already does — there is deliberately no preset for them.
"""
const FRAG_BOUND_PRESETS = Dict{String, @NamedTuple{
    low_slope::Float32, low_intercept::Float32,
    high_slope::Float32, high_intercept::Float32,
}}(
    # Both bounds flat: frag_mz_min / frag_mz_max are used as-is.
    "constant" => (low_slope = 0.0f0, low_intercept = 0.0f0,
                   high_slope = 0.0f0, high_intercept = 0.0f0),
    "thermo_auto" => (low_slope = 0.0f0, low_intercept = 0.0f0,
                      high_slope = 2.04f0, high_intercept = 24.1f0),
    "thermo_auto_documented" => (low_slope = 0.0f0, low_intercept = 0.0f0,
                                 high_slope = 2.0f0, high_intercept = 10.0f0),
)

"""
    frag_bound_polynomials(spec, default_frag_bounds) -> (low, high)

The pair of polynomials a `frag_bounds` spec describes.

`nothing`, or a spec resolving to flat bounds, gives degree-0 polynomials at
`default_frag_bounds` — byte-identical to the behaviour before this key existed,
which is what keeps every params file written to date building the same library.
"""
function frag_bound_polynomials(
    spec::Union{Nothing, @NamedTuple{
        low_slope::Float32, low_intercept::Float32,
        high_slope::Float32, high_intercept::Float32,
    }},
    default_frag_bounds::Tuple{Float32, Float32},
)
    lo_const, hi_const = first(default_frag_bounds), last(default_frag_bounds)
    spec === nothing && return (ImmutablePolynomial(lo_const), ImmutablePolynomial(hi_const))

    # A zero slope keeps the configured constant rather than falling back to the
    # rule's intercept, so "constant" and an absent key agree.
    low = spec.low_slope == 0.0f0 ?
        ImmutablePolynomial(lo_const) :
        ImmutablePolynomial(Float32[spec.low_intercept, spec.low_slope])
    high = spec.high_slope == 0.0f0 ?
        ImmutablePolynomial(hi_const) :
        ImmutablePolynomial(Float32[spec.high_intercept, spec.high_slope])
    return (low, high)
end

"""
    parse_frag_bounds_spec(value) -> spec or nothing

Read the optional `library_params.frag_bounds` key.

Accepts nothing (flat bounds — the behaviour before this key existed), a preset
name, or an object giving `low`/`high` slope and intercept explicitly. Anything
else throws, listing what is valid: a misspelled preset must not quietly build a
library with the wrong fragment ceiling.
"""
function parse_frag_bounds_spec(value)
    value === nothing && return nothing

    if value isa AbstractString
        key = lowercase(strip(String(value)))
        haskey(FRAG_BOUND_PRESETS, key) || throw(InvalidParametersError(
            "Unknown frag_bounds preset \"$value\". Valid names: " *
            join(sort(collect(keys(FRAG_BOUND_PRESETS))), ", ") *
            ". Omit the key for constant bounds, or give explicit " *
            "{\"low\": {...}, \"high\": {...}} coefficients.",
            Dict{String, Any}("frag_bounds" => value)))
        return FRAG_BOUND_PRESETS[key]
    end

    if value isa AbstractDict
        # A preset and hand-written coefficients together is ambiguous; say so
        # rather than silently honouring one of them.
        if haskey(value, "preset")
            throw(InvalidParametersError(
                "frag_bounds takes either a preset name or explicit low/high " *
                "coefficients, not both. Give the name as a bare string, or " *
                "drop \"preset\" and keep the coefficients.",
                Dict{String, Any}("frag_bounds" => value)))
        end
        coeff(side) = begin
            side_val = get(value, side, nothing)
            side_val isa AbstractDict || throw(InvalidParametersError(
                "frag_bounds.$side must be an object with numeric \"slope\" and " *
                "\"intercept\".", Dict{String, Any}("frag_bounds" => value)))
            slope = get(side_val, "slope", 0.0)
            intercept = get(side_val, "intercept", 0.0)
            (slope isa Real && intercept isa Real) || throw(InvalidParametersError(
                "frag_bounds.$side slope and intercept must be numbers.",
                Dict{String, Any}("frag_bounds" => value)))
            (Float32(slope), Float32(intercept))
        end
        low_slope, low_intercept = coeff("low")
        high_slope, high_intercept = coeff("high")
        return (low_slope = low_slope, low_intercept = low_intercept,
                high_slope = high_slope, high_intercept = high_intercept)
    end

    throw(InvalidParametersError(
        "frag_bounds must be a preset name or an object with low/high " *
        "coefficients; got $(typeof(value)).",
        Dict{String, Any}("frag_bounds" => value)))
end

function get_fragment_bounds(
    auto_detect_frag_bounds::Bool,
    frag_bounds_detection_raw_file_path::String,
    default_frag_bounds::Tuple{Float32, Float32},
    default_precursor_bounds::Tuple{Float32, Float32},
    frag_bounds_spec::Union{Nothing, @NamedTuple{
        low_slope::Float32, low_intercept::Float32,
        high_slope::Float32, high_intercept::Float32,
    }} = nothing)::@NamedTuple{frag_bounds::FragBoundModel, prec_mz_min::Float32, prec_mz_max::Float32}

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
                @user_warn "failed to estimate fragbounds from the example raw file.
                Using default values. Frag bounds: $default_frag_bounds, precursor bounds: $default_precursor_bounds"
            end
        else
            @user_warn "Could not find example file $frag_bounds_detection_raw_file_path to determine frag bounds. 
            Using default values. Frag bounds: $default_frag_bounds, precursor bounds: $default_precursor_bounds"
        end
    end
    low_poly, high_poly = frag_bound_polynomials(frag_bounds_spec, default_frag_bounds)
    # frag_mz_min / frag_mz_max stop being "the bounds" and become the
    # instrument's absolute limits once a bound has a slope. Clamping to them is
    # what keeps a sloped ceiling from running past anything the instrument can
    # record at the top of a wide precursor range.
    frag_bounds = FragBoundModel(
        low_poly,
        high_poly,
        first(default_frag_bounds),
        last(default_frag_bounds),
    )
    return (frag_bounds = frag_bounds,
            prec_mz_min = first(default_precursor_bounds),
            prec_mz_max = last(default_precursor_bounds))
end