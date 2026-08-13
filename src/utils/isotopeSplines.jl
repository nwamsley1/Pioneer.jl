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
    DecodeCoefficients(encoded::String) -> Vector{Float64}

Base64-decode a string of packed Float64 spline coefficients (from Goldfarb XML).
"""
function DecodeCoefficients(encoded::String)
    return reinterpret(Float64, Base64.base64decode(encoded))
end

"""
    CubicSpline{N, T}

Piecewise cubic spline interpolant with N packed coefficients (groups of 4).
Callable: `spline(mass)` evaluates via Horner's method with `@fastmath`.
"""
struct CubicSpline{N, T<:AbstractFloat}
    coeffs::SVector{N, T}
    first::T
    last::T
    bin_width::T
    inv_bin_width::T
end

function (s::CubicSpline)(t::U) where {U<:AbstractFloat}
    @inbounds @fastmath begin
        if t < s.first
            return s.coeffs[1]
        end

        idx = floor(Int32, (t - s.first)*s.inv_bin_width)
        u = t - (s.first + s.bin_width*idx)

        coeff = idx*4 + 1
        a0 = s.coeffs[coeff]
        a1 = s.coeffs[coeff + 1]
        a2 = s.coeffs[coeff + 2]
        a3 = s.coeffs[coeff + 3]

        x = muladd(a3, u, a2)
        x = muladd(x, u, a1)
        x = muladd(x, u, a0)
    end
    return x
end

"""
    QuadTransmission{T}

Quartic quadrupole transmission model with overhang and steepness parameters.
Callable: `qtf(window_center, window_half_width, x)` returns transmission probability.
"""
struct QuadTransmission{T<:AbstractFloat}
    overhang::T
    b::T
end

function (qtf::QuadTransmission)(window_center::T, window_half_width::T, x::T) where {T<:AbstractFloat}
    return one(T)/(one(T) + abs((x - (window_center))/(window_half_width + qtf.overhang))^(T(2)*qtf.b))
end

"""
    IsotopeSplineModel{T}

Container for isotope probability splines indexed by [sulfur_count][isotope_index].
Loaded from Goldfarb et al. 2018 XML via [`parseIsoXML`](@ref).
Callable: `model(sulfur_count, isotope_idx, mass)` returns isotope probability.
"""
# 40 is hardcoded rather than carried as a type parameter. parseIsoXML already builds these as
# `CubicSpline{40, Float32}` with `SVector{40, Float32}` coefficients, so the parameter only ever
# held the constant 40 while forcing every signature that mentions the model to spell it out. A
# different coefficient count in the XML still fails loudly, at the SVector{40} conversion.
struct IsotopeSplineModel{T<:Real}
    splines::Vector{Vector{CubicSpline{40, T}}}
end

function (p::IsotopeSplineModel)(S, I, x)
    return p.splines[S::Int64 + 1][I::Int64 + 1](x::Float32)
end

"""
    parseIsoXML(iso_xml_path::String) -> IsotopeSplineModel{Float32}

Load isotope probability splines from the Goldfarb et al. 2018 XML file.
Returns an `IsotopeSplineModel` indexed by [sulfur_count][isotope_index].
"""
function parseIsoXML(iso_xml_path::String)
    xdoc = parse_file(iso_xml_path)

    # Determine dimensions: max sulfur count and max isotope index
    max_S, max_iso = 0, 0
    for model in LightXML.root(xdoc)["model"]
        if haskey(attributes_dict(model), "S")
            max_S = max(max_S, parse(Int64, attributes_dict(model)["S"]) + 1)
            max_iso = max(max_iso, parse(Int64, attributes_dict(model)["isotope"]) + 1)
        end
    end

    # Pre-allocate with zero splines
    splines = Vector{Vector{CubicSpline{40, Float32}}}()
    zero_coeffs = SVector{40, Float32}(zeros(Float32, 40))
    for i in 1:max_S
        push!(splines, [CubicSpline(zero_coeffs, 0.0f0, 0.0f0, 0.0f0, 0.0f0) for _ in 1:max_iso])
    end

    # Fill splines from XML
    for model in LightXML.root(xdoc)["model"]
        if haskey(attributes_dict(model), "S")
            S = parse(Int64, attributes_dict(model)["S"])
            iso = parse(Int64, attributes_dict(model)["isotope"])
            knots = collect(Float32.(DecodeCoefficients(content(model["knots"][1]))))[1:end - 1]
            coefficients = Float32.(DecodeCoefficients(content(model["coefficients"][1])))
            bin_width = (last(knots) - first(knots)) / (length(knots) - 1)
            splines[S+1][iso+1] = CubicSpline(
                SVector{length(coefficients), Float32}(coefficients),
                Float32(first(knots)),
                Float32(last(knots)),
                Float32(bin_width),
                Float32(1.0f0 / bin_width)
            )
        end
    end

    return IsotopeSplineModel(splines)
end

"""
    isotope{T, I}

Lightweight representation of an isotopic species (fragment or precursor) for
the Goldfarb abundance calculation.  Lowercase name to avoid collision with
`Isotope <: LibraryIon` in `src/utils/isotopes.jl`.

Fields: `mass` (neutral), `sulfurs`, `iso` (isotope state).
"""
struct isotope{T<:AbstractFloat,I<:Int}
    mass::T
    sulfurs::I
    iso::I
end

import Base: -
function -(a::isotope{T, I}, b::isotope{T, I}) where {T<:Real,I<:Integer}
    return isotope(
        a.mass - b.mass,
        a.sulfurs - b.sulfurs,
        a.iso - b.iso
    )
end

#############################################################################
# getFragAbundance! — core Goldfarb algorithm
#
# Four method overloads:
#   1. (isotopes, iso_splines, frag, prec, pset)           — isolation-set mode
#   2. (frag_isotopes, prec_isotopes, iso_splines, frag, prec) — transmission mode
#   3. (isotopes, iso_splines, prec_mz, ..., frag::LibraryFragmentIon, pset) — wrapper for (1)
#   4. (frag_isotopes, precursor_transmission, iso_splines, ..., frag::LibraryFragmentIon) — wrapper for (2)
#############################################################################

"""
    getFragAbundance!(isotopes, iso_splines, frag, prec, pset)

Compute relative fragment isotope abundances in place using the method of
Goldfarb et al. 2018 (ACS Omega 3(9):11383-11391).

`isotopes[1]` = M+0, `isotopes[2]` = M+1, etc.  Result is **not** normalised
to sum to one -- the caller is responsible for normalisation.

# Arguments
- `isotopes::Vector{T}` -- output; filled with unnormalised fragment isotope abundances.
- `iso_splines::IsotopeSplineModel` -- splines returning P(isotope | sulfurs, mass).
- `frag::isotope{T,I}` -- fragment neutral mass and sulfur count.
- `prec::isotope{T,I}` -- precursor neutral mass and sulfur count.
- `pset::Tuple{I,I}` -- (first, last) isolated precursor isotope indices.
"""
function getFragAbundance!(isotopes::Vector{T},
                            iso_splines::IsotopeSplineModel,
                            frag::isotope{T, I},
                            prec::isotope{T, I},
                            pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    min_p, max_p = first(pset), last(pset)
    frag_sulfurs = min(frag.sulfurs, 5)
    comp_sulfurs = min(prec.sulfurs - frag.sulfurs, 5)
    frag_mass = Float32(frag.mass)
    comp_mass = Float32(prec.mass - frag.mass)
    @inbounds @fastmath for f in range(0, min(length(isotopes)-1, max_p))
        # Eq. 5, Goldfarb et al. 2018 pg. 11389
        complement_prob = 0.0
        f_i = iso_splines(frag_sulfurs, f, frag_mass)

        for p in range(max(f, min_p), max_p)
            complement_prob += iso_splines(comp_sulfurs, p - f, comp_mass)
        end

        isotopes[f+1] = f_i*complement_prob
    end
    return nothing
end

"""
    getFragAbundance!(frag_isotopes, prec_isotopes, iso_splines, frag, prec)

Transmission-mode variant: `prec_isotopes[i]` gives the transmission probability
for each precursor isotope (from quadrupole model) instead of a discrete set.
"""
function getFragAbundance!(frag_isotopes::Vector{T},
                            prec_isotopes::Vector{T},
                            iso_splines::IsotopeSplineModel,
                            frag::isotope{T, I},
                            prec::isotope{T, I}
                        ) where {T<:Real,I<:Integer}
    frag_sulfurs = min(frag.sulfurs, 5)
    comp_sulfurs = min(prec.sulfurs - frag.sulfurs, 5)
    frag_mass = Float32(frag.mass)
    comp_mass = Float32(prec.mass - frag.mass)
    @inbounds @fastmath for f in range(0, length(prec_isotopes)-1)
        complement_prob = 0.0
        f_i = iso_splines(frag_sulfurs, f, frag_mass)
        for p in range(max(f, 0), length(prec_isotopes) - 1)
            complement_prob += iso_splines(comp_sulfurs, p - f, comp_mass) * prec_isotopes[p + 1]
        end
        frag_isotopes[f+1] = f_i*complement_prob
    end
    return nothing
end

# Convenience wrapper: LibraryFragmentIon + isolation set -> isotope structs
function getFragAbundance!(isotopes::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32},
                            pset::Tuple{I, I}) where {I<:Integer}
    getFragAbundance!(
        isotopes,
        iso_splines,
        isotope(getMz(frag)*getFragCharge(frag), Int64(getSulfurCount(frag)), 0),
        isotope(prec_mz*prec_charge, Int64(prec_sulfur_count), 0),
        pset
    )
end

# Convenience wrapper: LibraryFragmentIon + transmission vector -> isotope structs
function getFragAbundance!(frag_isotopes::Vector{Float32},
                            precursor_transmission::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32})
    getFragAbundance!(
        frag_isotopes,
        precursor_transmission,
        iso_splines,
        isotope((getMz(frag) - Float32(PROTON))*getFragCharge(frag), Int64(getSulfurCount(frag)), 0),
        isotope((prec_mz - Float32(PROTON))*prec_charge, Int64(prec_sulfur_count), 0)
    )
end

#############################################################################
# getFragIsotopes! — normalised fragment isotope intensities
#
# Three method overloads:
#   1. (isotopes, iso_splines, prec_mz, ..., frag::LibraryFragmentIon, pset)
#   2. (isotopes, iso_splines, prec_mz, ..., frag::SplineDetailedFrag, knots, nce, pset)
#   3. (frag_isotopes, precursor_transmission, iso_splines, ..., frag::LibraryFragmentIon)
#############################################################################

"""
    getFragIsotopes!(isotopes, iso_splines, prec_mz, prec_charge, prec_sulfur_count, frag, prec_isotope_set)

Fill `isotopes` with fragment isotope intensities scaled by `frag.intensity`,
normalised so they sum to `frag.intensity`.
"""
function getFragIsotopes!(isotopes::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32},
                            prec_isotope_set::Tuple{Int64, Int64})
    fill!(isotopes, zero(eltype(isotopes)))
    total_fragment_intensity = frag.intensity

    getFragAbundance!(isotopes, iso_splines, prec_mz, prec_charge,
                      prec_sulfur_count, frag, prec_isotope_set)

    iso_sum = sum(isotopes)
    @inbounds @fastmath for i in reverse(range(1, length(isotopes)))
        isotopes[i] = total_fragment_intensity*isotopes[i]/iso_sum
    end
end

"""
    getFragIsotopes!(isotopes, iso_splines, prec_mz, prec_charge, prec_sulfur_count, frag::SplineDetailedFrag, knots, nce, prec_isotope_set)

Variant for `SplineDetailedFrag`: intensity is evaluated from the spline model
at the given NCE rather than taken from a stored field.
"""
function getFragIsotopes!(isotopes::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::SplineDetailedFrag{N, Float32},
                            knots::NTuple{M, Float32},
                            nce::Float32,
                            prec_isotope_set::Tuple{Int64, Int64}) where {M, N}
    fill!(isotopes, zero(eltype(isotopes)))
    total_fragment_intensity = getIntensity(frag, knots, 3, nce)

    getFragAbundance!(isotopes, iso_splines, prec_mz, prec_charge,
                      prec_sulfur_count, frag, prec_isotope_set)

    iso_sum = sum(isotopes)
    @inbounds @fastmath for i in reverse(range(1, length(isotopes)))
        isotopes[i] = total_fragment_intensity*isotopes[i]/iso_sum
    end
end

"""
    getFragIsotopes!(frag_isotopes, precursor_transmission, iso_splines, prec_mz, prec_charge, prec_sulfur_count, frag)

Transmission-mode variant: uses measured precursor isotope transmission
probabilities instead of a discrete isolation set.
"""
function getFragIsotopes!(frag_isotopes::Vector{Float32},
                          precursor_transmission::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32})
    fill!(frag_isotopes, zero(eltype(frag_isotopes)))
    total_fragment_intensity = frag.intensity

    getFragAbundance!(frag_isotopes, precursor_transmission, iso_splines,
                      prec_mz, prec_charge, prec_sulfur_count, frag)

    @inbounds @fastmath for i in reverse(range(1, length(frag_isotopes)))
        frag_isotopes[i] = total_fragment_intensity*frag_isotopes[i]
    end
end

#############################################################################
# getPrecursorIsotopeSet — which precursor isotopes fall in the isolation window
#############################################################################

# Maximum isotope index considered (M+0 .. M+MAX_PRECURSOR_ISOTOPE)
const MAX_PRECURSOR_ISOTOPE = 5

"""
    getPrecursorIsotopeSet(prec_mz, prec_charge, min_prec_mz, max_prec_mz)

Given quadrupole isolation bounds, return `(first_iso, last_iso)` indicating
which precursor isotopes (M+0..M+$MAX_PRECURSOR_ISOTOPE) fall within the window.
Returns `(-1, -1)` if no isotope is captured.
"""
function getPrecursorIsotopeSet(prec_mz::Float32,
                                prec_charge::UInt8,
                                min_prec_mz::Float32,
                                max_prec_mz::Float32)
    first_iso, last_iso = -1, -1
    @fastmath for iso_count in range(0, MAX_PRECURSOR_ISOTOPE)
        iso_mz = iso_count*C13_C12_MASS_DIFF/prec_charge + prec_mz
        if (iso_mz > min_prec_mz) & (iso_mz < max_prec_mz)
            if first_iso < 0
                first_iso = iso_count
            end
            last_iso = iso_count
        end
    end
    return (first_iso, last_iso)
end

"""
    getPrecursorIsotopeSet(prec_mz, prec_charge, qtf::QuadTransmissionFunction)

Convenience overload extracting m/z bounds from a `QuadTransmissionFunction`.
"""
function getPrecursorIsotopeSet(prec_mz::Float32,
                                prec_charge::UInt8,
                                qtf::QuadTransmissionFunction)
    return getPrecursorIsotopeSet(prec_mz, prec_charge,
                                  getPrecMinBound(qtf), getPrecMaxBound(qtf))
end

#############################################################################
# getPrecursorIsotopeTransmission! / getPrecursorFractionTransmitted!
#############################################################################

"""
    getPrecursorIsotopeTransmission!(prec_isotope_transmission, prec_mono_mz, prec_charge, qtf)

Fill `prec_isotope_transmission[i]` with the quadrupole transmission probability
for the i-th precursor isotope.
"""
function getPrecursorIsotopeTransmission!(
                                            prec_isotope_transmission::Vector{Float32},
                                            prec_mono_mz::Float32,
                                            prec_charge::UInt8,
                                            qtf::QuadTransmissionFunction)
    fill!(prec_isotope_transmission, zero(Float32))
    prec_iso_mz = prec_mono_mz
    @inbounds @fastmath for i in range(1, length(prec_isotope_transmission))
        prec_isotope_transmission[i] = qtf(prec_iso_mz)
        prec_iso_mz += Float32(C13_C12_MASS_DIFF/prec_charge)
    end
end

"""
    getPrecursorFractionTransmitted!(iso_splines, precursor_isotopes, qtf, prec_mono_mz, prec_charge, sulfur_count)

Compute the total fraction of precursor signal transmitted through the quadrupole,
allocating a temporary transmission buffer internally.
"""
function getPrecursorFractionTransmitted!(
    iso_splines::IsotopeSplineModel{Float32},
    precursor_isotopes::Tuple{I, I},
    qtf::QuadTransmissionFunction,
    prec_mono_mz::Float32,
    prec_charge::UInt8,
    sulfur_count::UInt8,
) where {I<:Real}

    last_iso = floor(Int, last(precursor_isotopes))
    last_iso < 1 && return 0.0f0
    precursor_transmission = zeros(Float32, last_iso)
    getPrecursorIsotopeTransmission!(precursor_transmission, prec_mono_mz, prec_charge, qtf)
    return _precursor_fraction_transmitted(
        precursor_transmission, iso_splines, precursor_isotopes,
        prec_mono_mz, prec_charge, sulfur_count,
    )
end

"""
    getPrecursorFractionTransmitted!(precursor_transmission, iso_splines, precursor_isotopes, qtf, prec_mono_mz, prec_charge, sulfur_count)

Pre-allocated variant: reuses `precursor_transmission` buffer to avoid allocation.
"""
function getPrecursorFractionTransmitted!(
    precursor_transmission::AbstractVector{Float32},
    iso_splines::IsotopeSplineModel{Float32},
    precursor_isotopes::Tuple{I, I},
    qtf::QuadTransmissionFunction,
    prec_mono_mz::Float32,
    prec_charge::UInt8,
    sulfur_count::UInt8,
) where {I<:Real}

    getPrecursorIsotopeTransmission!(precursor_transmission, prec_mono_mz, prec_charge, qtf)
    return _precursor_fraction_transmitted(
        precursor_transmission, iso_splines, precursor_isotopes,
        prec_mono_mz, prec_charge, sulfur_count,
    )
end

@inline function _precursor_fraction_transmitted(
    precursor_transmission::AbstractVector{Float32},
    iso_splines::IsotopeSplineModel{Float32},
    precursor_isotopes::Tuple{I, I},
    prec_mono_mz::Float32,
    prec_charge::UInt8,
    sulfur_count::UInt8,
) where {I<:Real}
    probability = 0.0f0
    first_iso = floor(Int, first(precursor_isotopes))
    last_iso = floor(Int, last(precursor_isotopes))
    last_iso < 1 && return probability
    first_idx = max(first_iso, 1)
    last_idx = min(last_iso, length(precursor_transmission))
    last_idx < first_idx && return probability

    precursor_mass = (prec_mono_mz * prec_charge) - prec_charge
    sulfur_idx = min(Int64(sulfur_count), 5)
    @inbounds @fastmath for iso in first_idx:last_idx
        probability += iso_splines(sulfur_idx, iso - 1, precursor_mass) * precursor_transmission[iso]
    end
    return probability
end
