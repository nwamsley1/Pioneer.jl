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
function getFragAbundance(iso_splines::IsotopeSplineModel{Float64}, frag::isotope{T, I}, prec::isotope{T, I}, pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope

    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    isotopes = zeros(Float64, max_p + 1)
    for f in range(0, max_p) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018

        f_i = coalesce(iso_splines(min(frag.sulfurs, 5), f, Float64(frag.mass)), 0.0) #Probability of fragment isotope in state 'f' assuming full precursor distribution 

        for p in range(max(f, min_p), max_p) #Probabilities of complement fragments 
            complement_prob += coalesce(iso_splines(min(prec.sulfurs - frag.sulfurs, 5), p - f, Float64(prec.mass - frag.mass)), 0.0)
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
function getFragAbundance!(isotopes::Vector{Float64}, iso_splines::IsotopeSplineModel{Float64}, frag::isotope{T, I}, prec::isotope{T, I}, pset::Tuple{I, I}) where {T<:Real,I<:Integer}
    #Approximating Isotope Distributions of Biomolecule Fragments, Goldfarb et al. 2018 
    min_p, max_p = first(pset), last(pset) #Smallest and largest precursor isotope

    #placeholder for fragment isotope distributions
    #zero to isotopic state of largest precursor 
    for f in range(0, min(length(isotopes)-1, max_p)) #Fragment cannot take an isotopic state grater than that of the largest isolated precursor isotope
        complement_prob = 0.0 #Denominator in 5) from pg. 11389, Goldfarb et al. 2018

        #Splines don't go above five sulfurs
        f_i = coalesce(iso_splines(min(frag.sulfurs, 5), f, Float64(frag.mass)), 0.0) #Probability of fragment isotope in state 'f' assuming full precursor distribution 

        for p in range(max(f, min_p), max_p) #Probabilities of complement fragments 
            #Splines don't go above five sulfurs 
            complement_prob += coalesce(iso_splines(min(prec.sulfurs - frag.sulfurs, 5), p - f, Float64(prec.mass - frag.mass)), 0.0)
        end
        isotopes[f+1] = f_i*complement_prob
    end

    return isotopes#isotopes./sum(isotopes)
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


isotopes = [0.0, 0.0]
getFragAbundance!(
    isotopes,
    iso_splines,
    isotope(Float64(frag.frag_mz*frag.frag_charge), Int64(frag.sulfur_count), 0),
    isotope(Float64(prec.mz*prec.charge), Int64(prec.sulfur_count), 0),
    prec_isotope_set)

p = plot([0, 0], [0, 0])
isotopes = [0.0, 0.0]
prec_isotope_set = (1, 5)
prec = precursors[1000000]
for frag in f_det[1000000]

    getFragAbundance!(
        isotopes,
        iso_splines,
        isotope(Float64(frag.frag_mz*frag.frag_charge), Int64(frag.sulfur_count), 0),
        isotope(Float64(prec.mz*prec.charge), Int64(prec.sulfur_count), 0),
        prec_isotope_set)

    intensity = frag.intensity
    isotopes = isotopes./sum(isotopes)
    frag_mz = frag.frag_mz
    iso_mz = frag.frag_mz + NEUTRON/frag.frag_charge
    plot!([frag_mz, frag_mz], 
    [0.0, intensity*isotopes[1]], 
    show = true, color = :black, legend = false, alpha = 0.5)

    plot!([iso_mz, iso_mz], 
    [0.0, intensity*isotopes[2]], 
    show = true, color = :red, legend = false, alpha = 0.5)

end
plot()
#=
prec = precursors[1]
for frag in f_det[1]
    isotopes = getFragAbundance(iso_splines, 
                        isotope(Float64(frag.frag_mz*frag.frag_charge), Int64(frag.sulfur_count), 0), 
                        isotope(Float64(prec.mz*prec.charge), Int64(prec.sulfur_count), 0), 
                        (1, 3))

    println(isotopes[1:2]./sum(isotopes[1:2]))
end
=#
#=
struct LibraryFragmentS{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    frag_charge::UInt8
    is_y_ion::Bool
    ion_position::UInt8
    ion_index::UInt8
    intensity::Float32
    prec_charge::UInt8
    prec_id::UInt32
    rank::UInt8
    sulfur_count::UInt8
end

#getIntensity(f::LibraryFragment) = f.intensity
#isyIon(f::LibraryFragment) = f.is_y_ion
#getIonIndex(f::LibraryFragment) = f.ion_index
#getIonPosition(f::LibraryFragment) = f.ion_position
#getFragCharge(f::LibraryFragment) = f.frag_charge
#getRank(f::LibraryFragment) = f.rank
#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8))

f_det_sulfur = Vector{Vector{LibraryFragmentS{Float32}}}()
for prec in ProgressBar(f_det)
    push!(f_det_sulfur, Vector{LibraryFragmentS{Float32}}())
    precursor = precursors[first(prec).prec_id]
    sulfurs = count(r"[CM]", precursor.sequence)
    sulfurs = 0
    for frag in prec
        if frag.is_y_ion
            sulfurs = count(r"[CM]", reverse(precursor.sequence)[1:frag.ion_position])
        else
            sulfurs = count(r"[CM]", precursor.sequence[1:frag.ion_position])
        end
        push!(last(f_det_sulfur), 
            LibraryFragmentS(
                frag.frag_mz,
                frag.frag_charge,
                frag.is_y_ion,
                frag.ion_position,
                frag.ion_index,
                frag.intensity,
                frag.prec_charge,
                frag.prec_id,
                frag.rank,
                UInt8(sulfurs)
            ))
    end
end

struct LibraryPrecursorS{T<:AbstractFloat}
    iRT::T
    mz::T
    total_intensity::T
    base_peak_intensity::T
    isDecoy::Bool
    charge::UInt8
    pep_id::UInt32
    prot_ids::Vector{UInt32}
    accession_numbers::String
    sequence::String
    missed_cleavages::UInt8
    variable_mods::UInt8
    length::UInt8
    sulfur_count::UInt8
end

precursors_sulfur = Vector{LibraryPrecursorS{Float32}}()
for precursor in ProgressBar(precursors)
    sulfurs = count(r"[CM]", precursor.sequence)
    push!(precursors_sulfur, 
            LibraryPrecursorS(
                precursor.iRT,
                precursor.mz,
                precursor.total_intensity,
                precursor.base_peak_intensity,
                precursor.isDecoy,
                precursor.charge,
                precursor.pep_id,
                precursor.prot_ids,
                precursor.accession_numbers,
                precursor.sequence,
                precursor.missed_cleavages,
                precursor.variable_mods,
                precursor.length,
                UInt8(sulfurs),
            )
    )
end

f_det_sulfur_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_f_det_sulfur.jld2"
jldsave("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_f_det_sulfur.jld2";
f_det_sulfur)

jldsave("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_precursors_sulfur.jld2";
precursors_sulfur)

f_det_sulfur = load("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_f_det_sulfur.jld2");
f_det_sulfur = f_det_sulfur["f_det_sulfur"];

precursors_sulfur = load("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_precursors_sulfur.jld2");
precursors_sulfur =precursors_sulfur["precursors_sulfur"];


f_det = Vector{Vector{LibraryFragment{Float32}}}()
for prec in ProgressBar(f_det_sulfur)
    push!(f_det, Vector{LibraryFragmentS{Float32}}())
    for frag in prec
        push!(last(f_det), 
            LibraryFragment(
                frag.frag_mz,
                frag.frag_charge,
                frag.is_y_ion,
                false,
                frag.ion_position,
                frag.ion_index,
                frag.intensity,
                frag.prec_charge,
                frag.prec_id,
                frag.rank,
                frag.sulfur_count
            ))
    end
end

jldsave("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123_sulfur/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_f_det.jld2";
f_det)

precursors = Vector{LibraryPrecursor{Float32}}()
for precursor in ProgressBar(precursors_sulfur)
    #sulfurs = count(r"[CM]", precursor.sequence)
    push!(precursors, 
            LibraryPrecursor(
                precursor.iRT,
                precursor.mz,
                precursor.total_intensity,
                precursor.base_peak_intensity,
                precursor.isDecoy,
                precursor.charge,
                precursor.pep_id,
                precursor.prot_ids,
                precursor.accession_numbers,
                precursor.sequence,
                precursor.missed_cleavages,
                precursor.variable_mods,
                precursor.length,
                precursor.sulfur_count
            )
    )
end

jldsave("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123_sulfur/HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_precursors.jld2";
precursors)
=#

#=
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml");

#When the full precursor distribution is isolated the fragment isotope distribution should look like its natural distribution 
getFragAbundance(iso_splines, isotope(2500.0, 0, 0), isotope(3000.0, 0, 0), [0, 1])
getFragAbundance(iso_splines, isotope(2500.0, 0, 0), isotope(3000.0, 0, 0), [0, 1, 2, 3, 4, 5, 6, 7])
[iso_splines(0,x,2500.0) for x in [0, 1, 2, 3, 4, 5, 6, 7]]




=#
