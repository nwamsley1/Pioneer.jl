using LightXML, Base64, Polynomials

function DecodeCoefficients(encoded::String)
    return reinterpret(Float64, Base64.base64decode(encoded))
end
struct CubicPolynomial{T<:Real}
    x0::T
    x1::T
    x2::T
    x3::T
end

function (a::CubicPolynomial{T})(x::T) where {T<:Real}
    return @fastmath a.x0 + a.x1*x + a.x2*x^2 + a.x3*x^3
end
struct PolynomialSpline{T<:Real}
    polynomials::Vector{CubicPolynomial{T}}
    slope::T
    intercept::T
    knots::Vector{T}
end



function (p::PolynomialSpline)(x)
    #idx = searchsortedfirst(p.knots, x)
    
    #idx = max(ceil(Int, first(p.knots) + x*last(p.knots)), 1)
    idx = ceil(Int, p.intercept + x*p.slope)
    #if (idx == 1) | (idx > length(p.knots))
    #    return missing
    #end
    if (idx == 0)
        return max(p.polynomials[1](x - p.knots[1]), zero(Float32))
    elseif (idx > length(p.knots))
        return p.polynomials[end](x - p.knots[end])
    else
        return p.polynomials[idx](x - p.knots[idx])
    end
end

struct IsotopeSplineModel{T<:Real}
    splines::Vector{Vector{PolynomialSpline{T}}}
end

function (p::IsotopeSplineModel)(S, I, x)
    return p.splines[S::Int64 + 1][I::Int64 + 1](x::Float32)
end


function buildPolynomials(coefficients::Vector{T}, order::I) where {T<:Real, I<:Integer}
    order += 1
    n_knots = length(coefficients)÷order
    polynomials = Vector{CubicPolynomial{T}}()
    for n in range(1, n_knots)
        start = ((n - 1)*order + 1)
        stop = (n)*order
        push!(polynomials, CubicPolynomial(coefficients[start], coefficients[start+1], coefficients[start+2], coefficients[start+3]))
    end
    return polynomials
end

function parseIsoXML(iso_xml_path::String)
    #From LightXML.jl
    xdoc = parse_file(iso_xml_path)

    max_S, max_iso = 0, 0
    for model in root(xdoc)["model"]
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
    splines = Vector{Vector{PolynomialSpline{Float32}}}()
    for i in range(1, max_S)
        push!(splines, [])
        for j in range(1, max_iso)
            push!(
                splines[i], 
                PolynomialSpline(
                                buildPolynomials(Float32[0, 0, 0], 3),
                                zero(Float32),
                                zero(Float32),
                                Float32[0]
                                            )
            )
        end
    end

    #Fill Splines 
    for model in root(xdoc)["model"]
        if (haskey(attributes_dict(model),"S"))
            S = parse(Int64, attributes_dict(model)["S"])
            iso =  parse(Int64, attributes_dict(model)["isotope"]) 

            polynomials =                 buildPolynomials(
                collect(Float32.(DecodeCoefficients(content(model["coefficients"][1])))),
                parse(Int64, attributes_dict(model)["order"]) - 1
                )

            knots = collect(Float32.(DecodeCoefficients(content(model["knots"][1]))))[1:end - 1]
            A = hcat(ones(length(knots)), knots)
            x = A\collect(range(0, length(knots) - 1))
            splines[S+1][iso+1] = PolynomialSpline(
                polynomials,
                Float32(last(x)),
                Float32(first(x)),
                knots
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
    getFragAbundance(iso_splines::IsotopeSplineModel{Float64}, frag::isotope{T, I}, prec::isotope{T, I}, pset::Tuple{I, I}) where {T<:Real,I<:Integer}

Get the relative intensities of fragment isotopes starting at M+0. Returns `isotopes` where isotopes[1] is M+0, isotopes[2] is M+1, etc. 
Based on Goldfarb et al. 2018 Approximating Isotope Distributions of Biomolecule Fragments 
CS Omega 2018, 3, 9, 11383-11391
Publication Date:September 19, 2018
https://doi.org/10.1021/acsomega.8b01649

### Input

- `iso_splines::IsotopeSplineModel{Float64}` -- Splines from Goldfarb et. al. that return isotope probabilities given the number of sulfurs and average mass 
- `frag::isotope{T, I}` -- The fragment isotope
- `prec::isotope{T, I}` -- The precursor isotope
- `pset::Tuple{I, I}` -- The first and last precursor isotope that was isolated. (1, 3) would indicate the M+1 through M+3 isotopes were isolated and fragmented.

### Output

Returns `isotopes` where isotopes[1] is M+0, isotopes[2] is M+1, etc. Does not normalize to sum to one

### Notes

- See methods from Goldfarb et al. 2018

### Algorithm 

### Examples 

"""
function getFragAbundance(iso_splines::IsotopeSplineModel{T}, 
                            frag::isotope{T, I}, 
                            prec::isotope{T, I}, 
                            pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope

    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    isotopes = zeros(Float64, max_p + 1)
    for f in range(0, max_p) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0f0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018

        f_i = coalesce(iso_splines(min(frag.sulfurs, 5), f, Float32(frag.mass)), 0.0f0) #Probability of fragment isotope in state 'f' assuming full precursor distribution 

        for p in range(max(f, min_p), max_p) #Probabilities of complement fragments 
            complement_prob += coalesce(iso_splines(min(prec.sulfurs - frag.sulfurs, 5), p - f, Float32(prec.mass - frag.mass)), 0.0f0)
        end
        isotopes[f+1] = f_i*complement_prob
    end

    return isotopes#./sum(isotopes)
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
                            iso_splines::IsotopeSplineModel{T}, 
                            frag::isotope{T, I}, 
                            prec::isotope{T, I}, 
                            pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope
    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    for f in range(0, min(length(isotopes)-1, max_p)) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
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

    #return isotopes#isotopes./sum(isotopes)
end

function getFragAbundance!(isotopes::Vector{Float32}, 
                            iso_splines::IsotopeSplineModel{Float32},
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
                            iso_splines::IsotopeSplineModel{Float32}, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::LibraryFragmentIon{Float32}, 
                            prec_isotope_set::Tuple{Int64, Int64})
    fill!(isotopes, zero(eltype(isotopes)))

    monoisotopic_intensity = frag.intensity
    getFragAbundance!(isotopes, 
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag, 
                    prec_isotope_set)

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    for i in reverse(range(1, length(isotopes)))
        isotopes[i] = max(monoisotopic_intensity*isotopes[i]/first(isotopes), zero(Float32))
    end
end



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
                                max_prec_mz::Float32)
    first_iso, last_iso = -1, -1
    for iso_count in range(0, 5) #Arbitrary cutoff after 5 
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