using LightXML, Base64, Polynomials

function DecodeCoefficients(encoded::String)
    return reinterpret(Float64, Base64.base64decode(encoded))
end

struct PolynomialSpline{T<:Real}
    polynomials::Vector{Polynomial{T, :x}}
    knots::Vector{T}
end

function (p::PolynomialSpline)(x)
    idx = searchsortedfirst(p.knots, x)
    if (idx == 1) | (idx > length(p.knots))
        return missing
    end
    return p.polynomials[idx-1](x - p.knots[idx - 1])
end

struct IsotopeSplineModel{T<:Real}
    splines::Vector{Vector{PolynomialSpline{T}}}
end

function (p::IsotopeSplineModel)(S, I, x)
    return p.splines[S::Int64 + 1][I::Int64 + 1](x::Float64)
end


function buildPolynomials(coefficients::Vector{T}, order::I) where {T<:Real, I<:Integer}
    order += 1
    n_knots = length(coefficients)÷order
    polynomials = Vector{Polynomial{T, :x}}()
    for n in range(1, n_knots)
        start = ((n - 1)*order + 1)
        stop = (n)*order
        push!(polynomials, Polynomial(coefficients[start:stop], :x))
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
    splines = Vector{Vector{PolynomialSpline{Float64}}}()
    for i in range(1, max_S)
        push!(splines, [])
        for j in range(1, max_iso)
            push!(
                splines[i], 
                PolynomialSpline(
                                buildPolynomials(Float64[0, 0, 0], 3),
                                Float64[0]
                                            )
            )
        end
    end

    #Fill Splines 
    for model in root(xdoc)["model"]
        if (haskey(attributes_dict(model),"S"))
            S = parse(Int64, attributes_dict(model)["S"])
            iso =  parse(Int64, attributes_dict(model)["isotope"]) 
            splines[S+1][iso+1] = PolynomialSpline(
                buildPolynomials(
                collect(DecodeCoefficients(content(model["coefficients"][1]))),
                parse(Int64, attributes_dict(model)["order"]) - 1
                ),
                collect(DecodeCoefficients(content(model["knots"][1])))
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

function getFragAbundance(iso_splines::IsotopeSplineModel{Float64}, frag::isotope{T, I}, prec::isotope{T, I}, pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope

    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    isotopes = zeros(Float64, max_p + 1)

    for f in range(0, max_p) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018

        f_i = iso_splines(frag.sulfurs, f, frag.mass) #Probability of fragment isotope in state 'f' assuming full precursor distribution 

        for p in range(max(f, min_p), max_p) #Probabilities of complement fragments 
            complement_prob += iso_splines(prec.sulfurs - frag.sulfurs, p - f, prec.mass - frag.mass)
        end

        isotopes[f+1] = f_i*complement_prob
    end

    return isotopes./sum(isotopes)
end

#given the precursor window and precursor mass and charge, calculate which precursor fragments were isolated. 
function getPrecursorIsotopeSet(prec_mz::T, prec_charge::U, window::Tuple{T, T})where {T<:Real,U<:Unsigned}
    first_iso, last_iso = -1, -1
    for iso_count in range(0, 5) #Arbitrary cutoff after 5 
        iso_mz = iso_count*NEUTRON/prec_charge + prec_mz
        if (iso_mz > first(window)) & (iso_mz < last(window)) 
            if first_iso < 0
                first_iso = iso_count
            end
            last_iso = iso_count
        end
    end
    return (first_iso, last_iso)
end
#=
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")

#When the full precursor distribution is isolated the fragment isotope distribution should look like its natural distribution 
getFragAbundance(iso_splines, isotope(2500.0, 0, 0), isotope(3000.0, 0, 0), [0, 1])
getFragAbundance(iso_splines, isotope(2500.0, 0, 0), isotope(3000.0, 0, 0), [0, 1, 2, 3, 4, 5, 6, 7])
[iso_splines(0,x,2500.0) for x in [0, 1, 2, 3, 4, 5, 6, 7]]

=#
