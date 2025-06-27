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



function DecodeCoefficients(encoded::String)
    return reinterpret(Float64, Base64.base64decode(encoded))
end

struct CubicSpline{N, T<:AbstractFloat} 
    coeffs::SVector{N, T}
    first::T
    last::T
    bin_width::T
end

function (s::CubicSpline)(t::U) where {U<:AbstractFloat}
    @inbounds @fastmath begin
        if t < s.first
            u = 0.0f0
            idx = 0
        else
            idx = floor(Int32, 
                        (t - s.first)/s.bin_width
                        )
            u = (t - (s.first + s.bin_width*(idx)))
        end
        x = zero(U)
        coeff = (idx)*4 + 1
        c = one(U)
        x += s.coeffs[coeff]*c
        c *= u
        coeff += 1
        x += s.coeffs[coeff]*c
        c *= u
        coeff += 1
        x += s.coeffs[coeff]*c
        c *= u
        coeff += 1
        x += s.coeffs[coeff]*c
    end
    return x
end

struct QuadTransmission{T<:AbstractFloat}
    overhang::T
    b::T
end

function (qtf::QuadTransmission)(window_center::T, window_half_width::T, x::T) where {T<:AbstractFloat}
    return one(T)/(one(T) + abs((x - (window_center))/(window_half_width + qtf.overhang))^(T(2)*qtf.b))
end


struct IsotopeSplineModel{N, T<:Real}
    splines::Vector{Vector{CubicSpline{N, T}}}
end

function (p::IsotopeSplineModel)(S, I, x)
    return p.splines[S::Int64 + 1][I::Int64 + 1](x::Float32)
end


function buildPolynomials(coefficients::Vector{T}, order::I) where {T<:Real, I<:Integer}
    return SVector{length(coefficients)}(coefficients)#polynomials
end

function parseIsoXML(iso_xml_path::String)
    #From LightXML.jl
    xdoc = parse_file(iso_xml_path)

    max_S, max_iso = 0, 0
    for model in LightXML.root(xdoc)["model"]
        #Use only sulfur-specific models
        if (haskey(attributes_dict(model),"S"))
            if parse(Int64, attributes_dict(model)["S"])+1 > max_S
                max_S = parse(Int64, attributes_dict(model)["S"])+1
            end
            if parse(Int64, attributes_dict(model)["isotope"])+1 > max_iso
                max_iso = parse(Int64, attributes_dict(model)["isotope"])+1
            end
        end
    end

    #Pre-allocate splines 
    #splines = Vector{Vector{PolynomialSpline{Float32}}}()
    splines = Vector{Vector{CubicSpline{40, Float32}}}()
    for i in range(1, max_S)
        push!(splines, [])
        for j in range(1, max_iso)
            push!(
                splines[i], 
                CubicSpline(
                    @SVector[x for x in zeros(Float32, 40)],
                    0.0f0,
                    0.0f0,
                    0.0f0)
            )
        end
    end

    #Fill Splines 
    for model in LightXML.root(xdoc)["model"]
        if (haskey(attributes_dict(model),"S"))
            S = parse(Int64, attributes_dict(model)["S"])
            iso =  parse(Int64, attributes_dict(model)["isotope"]) 
            knots = collect(Float32.(DecodeCoefficients(content(model["knots"][1]))))[1:end - 1]
            A = hcat(ones(length(knots)), knots)
            x = A\collect(range(0, length(knots) - 1))
            spline = Float32.(DecodeCoefficients(content(model["coefficients"][1])))
            splines[S+1][iso+1] = CubicSpline(
                SVector{length(spline), Float32}(spline),
                Float32(first(knots)),
                Float32(last(knots)),
                Float32((last(knots)-first(knots))/(length(knots) - 1))
            )
        end
    end

    return IsotopeSplineModel(splines)

end

struct isotope{T<:AbstractFloat,I<:Int}
    mass::T
    sulfurs::I
    iso::I
end

import Base.-
function -(a::isotope{T, I}, b::isotope{T, I}) where {T<:Real,I<:Integer}
    return isotope(
        a.mass - b.mass,
        a.sulfurs - b.sulfurs,
        a.iso - b.iso
    )
end

"""
    getFragAbundance!(isotopes::Vector{Float64}, iso_splines::IsotopeSplineModel{Float64}, frag::isotope{T, I}, prec::isotope{T, I}, pset::Tuple{I, I}) where {T<:Real,I<:Integer}

Get the relative intensities of fragment isotopes starting at M+0. Fills `isotopes` in place. isotopes[1] is M+0, isotopes[2] is M+1, etc. 
Based on Goldfarb et al. 2018 Approximating Isotope Distributions of Biomolecule Fragments 
CS Omega 2018, 3, 9, 11383-11391
Publication Date:September 19, 2018
https://doi.org/10.1021/acsomega.8b01649

### Input

- `isotopes::Vector{Float64}`: -- Vector to hold relative abundances of fragment isotopes. 
- `iso_splines::IsotopeSplineModel{Float64}` -- Splines from Goldfarb et. al. that return isotope probabilities given the number of sulfurs and average mass 
- `frag::isotope{T, I}` -- The fragment isotope
- `prec::isotope{T, I}` -- The precursor isotope
- `pset::Tuple{I, I}` -- The first and last precursor isotope that was isolated. (1, 3) would indicate the M+1 through M+3 isotopes were isolated and fragmented.

### Output

Fills `isotopes` in place with the relative abundances of the fragment isotopes. Does not normalize to sum to one!

### Notes

- See methods from Goldfarb et al. 2018

### Algorithm 

### Examples 

"""
function getFragAbundance!(isotopes::Vector{T}, 
                            iso_splines::IsotopeSplineModel, 
                            frag::isotope{T, I}, 
                            prec::isotope{T, I}, 
                            pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope
    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    for f in range(0, min(length(isotopes)-1, max_p)) #Fragment cannot take an isoopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018

        #Splines don't go above five sulfurs
        f_i = coalesce(iso_splines(min(frag.sulfurs, 5), f, Float32(frag.mass)), 0.0) #Probability of fragment isotope in state 'f' assuming full precursor distribution 

        for p in range(max(f, min_p), max_p) #Probabilities of complement fragments 
            #Splines don't go above five sulfurs 
            complement_prob += coalesce(iso_splines(
                                                            min(prec.sulfurs - frag.sulfurs, 5), 
                                                            p - f, 
                                                            Float32(prec.mass - frag.mass)), 
                                        0.0)
        end

        isotopes[f+1] = f_i*complement_prob
    end
    return nothing
    #return isotopes#isotopes./sum(isotopes)
end

function getFragAbundance!(frag_isotopes::Vector{T}, #Probability of each fragment isotopes (initialized as zeros)
                            prec_isotopes::Vector{T}, #Probability of each precursor isotope (determined by quadrupole transmission function)
                            iso_splines::IsotopeSplineModel,
                            frag::isotope{T, I}, 
                            prec::isotope{T, I}
                        ) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    for f in range(0, length(prec_isotopes)-1) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018
        #Splines don't go above five sulfurs
        f_i = iso_splines(min(frag.sulfurs, 5), f, Float32(frag.mass)) #Probability of fragment isotope in state 'f' assuming full precursor distribution 
        for p in range(max(f, 0), length(prec_isotopes) - 1) #Probabilities of complement fragments 
            #Splines don't go above five sulfurs 
            complement_prob += iso_splines(
                                            min(prec.sulfurs - frag.sulfurs, 5), 
                                            p - f, 
                                            Float32(prec.mass - frag.mass)
                                        )*prec_isotopes[p + 1]
                                        
        end
        frag_isotopes[f+1] = f_i*complement_prob
    end
    return nothing
end


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
        isotope(frag.mz*frag.frag_charge, Int64(frag.sulfur_count), 0),
        isotope(prec_mz*prec_charge, Int64(prec_sulfur_count), 0),
        pset
        )
end


function getFragIsotopes!(isotopes::Vector{Float32}, 
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32}, 
                            prec_isotope_set::Tuple{Int64, Int64})
    #Reset relative abundances of isotopes to zero 
    fill!(isotopes, zero(eltype(isotopes)))

    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity = frag.intensity

    getFragAbundance!(isotopes, 
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag, 
                    prec_isotope_set)

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    #total_fragment_intensity /= sum(isotopes)
    iso_sum = sum(isotopes)
    #isotopes /= isotopes[1]
    #total_fragment_intensity = total_fragment_intensity*sum(isotopes)
    for i in reverse(range(1, length(isotopes)))
        isotopes[i] = total_fragment_intensity*isotopes[i]/iso_sum
    end
end

function getFragIsotopes!(isotopes::Vector{Float32}, 
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::SplineDetailedFrag{N, Float32}, 
                            knots::NTuple{M, Float32},
                            nce::Float32,
                            prec_isotope_set::Tuple{Int64, Int64}) where {M, N}
    #Reset relative abundances of isotopes to zero 
    fill!(isotopes, zero(eltype(isotopes)))

    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity = getIntensity(frag, knots, 3, nce)

    getFragAbundance!(isotopes, 
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag, 
                    prec_isotope_set)

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    #total_fragment_intensity /= sum(isotopes)
    iso_sum = sum(isotopes)
    #isotopes /= isotopes[1]
    #total_fragment_intensity = total_fragment_intensity*sum(isotopes)
    for i in reverse(range(1, length(isotopes)))
        isotopes[i] = total_fragment_intensity*isotopes[i]/iso_sum
    end
end

function getFragAbundance!(frag_isotopes::Vector{Float32}, 
                            precursor_transmition::Vector{Float32},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32})
    getFragAbundance!(
        frag_isotopes,
        precursor_transmition,
        iso_splines,
        isotope((frag.mz - Float32(PROTON))*frag.frag_charge, Int64(frag.sulfur_count), 0),
        isotope((prec_mz - Float32(PROTON))*prec_charge, Int64(prec_sulfur_count), 0)
        )
end

function getFragIsotopes!(frag_isotopes::Vector{Float32}, 
                          precursor_transmition::Vector{Float32},
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32})
    #Reset relative abundances of isotopes to zero 
    fill!(frag_isotopes, zero(eltype(frag_isotopes)))
    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity = frag.intensity

    getFragAbundance!(
                    frag_isotopes, 
                    precursor_transmition,
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag
                    )

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    for i in reverse(range(1, length(frag_isotopes)))
        frag_isotopes[i] = total_fragment_intensity*frag_isotopes[i]
    end
end
#=
function getFragIsotopes!(
                            ::PartialPrec,
                            frag_isotopes::Vector{Float32}, 
                            precursor_transmition::Vector{Float32},
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::SplineDetailedFrag{N, Float32}, 
                            knots::NTuple{M, Float32},
                            nce::Float32) where {M, N}
    #Reset relative abundances of isotopes to zero 
    fill!(frag_isotopes, zero(eltype(frag_isotopes)))
    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity =  getIntensity(frag, knots, 3, nce)

    getFragAbundance!(
                    frag_isotopes, 
                    precursor_transmition,
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag
                    )

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    for i in reverse(range(1, length(frag_isotopes)))
        frag_isotopes[i] = total_fragment_intensity*frag_isotopes[i]
    end
end
=#
"""
    getPrecursorIsotopeSet(prec_mz::T, prec_charge::U, window::Tuple{T, T})where {T<:Real,U<:Unsigned}

Given the quadrupole isolation window and the precursor mass and charge, calculates which precursor isotopes were isolated

### Input

- `prec_mz::T`: -- Precursor mass-to-charge ratio
- `prec_charge::U` -- Precursor charge state 
- ` window::Tuple{T, T}` -- The lower and upper m/z bounds of the quadrupole isolation window


### Output

A Tuple of two integers. (1, 3) would indicate the M+1 through M+3 isotopes were isolated and fragmented.

### Notes

- See methods from Goldfarb et al. 2018

### Algorithm 

### Examples 

"""
function getPrecursorIsotopeSet(prec_mz::Float32, 
                                prec_charge::UInt8, 
                                min_prec_mz::Float32, 
                                max_prec_mz::Float32;
                                max_iso::Int64 = 5)
    first_iso, last_iso = -1, -1
    for iso_count in range(0, max_iso) #Arbitrary cutoff after 5 
        iso_mz = iso_count*NEUTRON/prec_charge + prec_mz
        if (iso_mz > min_prec_mz) & (iso_mz < max_prec_mz) 
            if first_iso < 0
                first_iso = iso_count
            end
            last_iso = iso_count
        end
    end
    return (first_iso, last_iso)
end


function getPrecursorIsotopeSet(prec_mz::Float32, 
                                prec_charge::UInt8, 
                                qtf::QuadTransmissionFunction;
                                max_iso::Int64 = 5)
    first_iso, last_iso = -1, 0
    min_prec_mz = getPrecMinBound(qtf)
    max_prec_mz = getPrecMaxBound(qtf)
    for iso_count in range(0, max_iso) #Arbitrary cutoff after 5 
        iso_mz = iso_count*NEUTRON/prec_charge + prec_mz
        if (iso_mz > min_prec_mz) & (iso_mz < max_prec_mz) 
            if first_iso < 0
                first_iso = iso_count
            end
            last_iso = iso_count
        end
    end
    return (first_iso, last_iso)
end

function getPrecursorIsotopeTransmission!(
                                            prec_isotope_transmission::Vector{Float32},
                                            prec_mono_mz::Float32,
                                            prec_charge::UInt8,
                                            qtf::QuadTransmissionFunction)
    fill!(prec_isotope_transmission, zero(Float32))
    prec_iso_mz = prec_mono_mz
    for i in range(1, length(prec_isotope_transmission))
        prec_isotope_transmission[i] = qtf(prec_iso_mz)
        prec_iso_mz += Float32(NEUTRON/prec_charge)
    end
end


function getPrecursorFractionTransmitted!(
    iso_splines::IsotopeSplineModel{40, Float32},
    precursor_isotopes::Tuple{I, I},
    qtf::QuadTransmissionFunction,
    prec_mono_mz::Float32,
    prec_charge::UInt8,
    sulfur_count::UInt8,
    ) where {I<:Real}

    precursor_transmission = zeros(Float32, last(precursor_isotopes))
    getPrecursorIsotopeTransmission!(precursor_transmission, prec_mono_mz, prec_charge, qtf)
    probability = 0.0f0
    for iso in range(first(precursor_isotopes), last(precursor_isotopes))
        probability += iso_splines(min(Int64(sulfur_count), 5), Int64(iso-1), (prec_mono_mz*prec_charge) - prec_charge) * precursor_transmission[iso]
    end

    return probability
end


#=
function correctPrecursorAbundance(
    abundance::Float32,
    isotope_splines::IsotopeSplineModel{40, Float32},
    precursor_isotopes::Tuple{I, I},
    precursor_mass::Float32,
    sulfur_count::UInt8,
    ) where {I<:Real}
    #println("precursor_isotopes $precursor_isotopes")
    probability = 0.0f0
    for i in range(first(precursor_isotopes), last(precursor_isotopes))
        #println(isotope_splines(min(Int64(sulfur_count), 5), Int64(i), precursor_mass))
        probability += isotope_splines(min(Int64(sulfur_count), 5), Int64(i), precursor_mass)
    end
    return abundance/probability
end

function correctPrecursorAbundances!(
    abundances::AbstractVector{Float32},
    isotope_splines::IsotopeSplineModel{40, Float32},
    precursor_isotopes::AbstractVector{Tuple{I,I}},
    precursor_idxs::AbstractVector{UInt32},
    precursor_mzs::AbstractVector{Float32},
    precursor_charges::AbstractVector{UInt8},
    sulfur_counts::AbstractVector{UInt8}) where {I<:Real}
    for i in range(1, length(abundances))
        prec_isotopes = precursor_isotopes[i]
        prec_idx = precursor_idxs[i]
        sulfur_count = sulfur_counts[prec_idx]
        prec_mz = precursor_mzs[prec_idx]
        prec_charge = precursor_charges[prec_idx]
        abundances[i] = correctPrecursorAbundance(
            abundances[i],
            isotope_splines,
            prec_isotopes,
            prec_mz*prec_charge,
            sulfur_count
        )
    end
end

function correctPrecursorAbundance(
    abundance::Float32,
    precursor_transmission::Vector{Float32},
    isotope_splines::IsotopeSplineModel{40, Float32},
    precursor_isotopes::Tuple{I, I},
    precursor_mass::Float32,
    sulfur_count::UInt8,
    ) where {I<:Real}
    #println("precursor_isotopes $precursor_isotopes")
    probability = 0.0f0
    
    #p = sum(precursor_transmission)
    #abundance = abundance/p
    for i in range(0, length(precursor_transmission) - 1)
        #println(isotope_splines(min(Int64(sulfur_count), 5), Int64(i), precursor_mass))
        probability += isotope_splines(min(Int64(sulfur_count), 5), Int64(i), precursor_mass)*precursor_transmission[i + 1]
    end
    println("probability $probability sum(precursor_transmission) ", sum(precursor_transmission)/5)
    return abundance/probability
end

function correctPrecursorAbundances!(
    abundances::AbstractVector{Float32},
    precursor_transmission::Vector{Float32},
    qtf::QuadTransmission,
    isotope_splines::IsotopeSplineModel{40, Float32},
    window_mzs::AbstractVector{Float32},
    precursor_isotopes::AbstractVector{Tuple{I,I}},
    precursor_idxs::AbstractVector{UInt32},
    precursor_mzs::AbstractVector{Float32},
    precursor_charges::AbstractVector{UInt8},
    sulfur_counts::AbstractVector{UInt8}) where {I<:Real}

    for i in range(1, length(abundances))
        prec_isotopes = precursor_isotopes[i]
        prec_idx = precursor_idxs[i]
        sulfur_count = sulfur_counts[prec_idx]
        prec_mz = precursor_mzs[prec_idx]
        prec_charge = precursor_charges[prec_idx]
        window_mz = window_mzs[i]

        getPrecursorIsotopeTransmission!(
            precursor_transmission,
            prec_mz,
            prec_charge,
            window_mz,
            qtf
        )

        abundances[i] = correctPrecursorAbundance(
            abundances[i],
            precursor_transmission,
            isotope_splines,
            prec_isotopes,
            prec_mz*prec_charge,
            sulfur_count
        )

        println("prec_isotopes $prec_isotopes precursor_transmission $precursor_transmission prec_mz $prec_mz window_mz $window_mz")
    end
end
=#