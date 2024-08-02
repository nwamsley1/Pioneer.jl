#C(12):12.000000    98.8922
#C(13):13.0033548378    1.1078
#O(16):15.99491462219   99.7628
#O(17):16.99913150  0.0372
#O(18):17.9991604   0.20004
#S(32):31.97207069  95.018
#S(33):32.97145850  0.750
#S(34):33.96786683    4.215
#S(36):35.96708088    0.017
#N(14):14.0030740052    99.6337
#N(15):15.0001088984    0.3663

struct QRoots
    S::Vector{Vector{ComplexF64}}
    C::Vector{Vector{ComplexF64}}
    H::Vector{Vector{ComplexF64}}
    O::Vector{Vector{ComplexF64}}
    N::Vector{Vector{ComplexF64}}
    npeaks::UInt32
    function QRoots(npeaks::Integer)
        #=
        
        =#

        ###Get values from 
        #Sroots = roots(Polynomial([0.9499, 0.0075, 0.0425, 0.0, 0.0001]))
        Sroots = [ -0.09928332843567012 - 4.863784343504499im
                    -0.09928332843567012 + 4.863784343504499im
                    0.09928332843567106 - 20.03405391373626im
                    0.09928332843567106 + 20.03405391373626im]
        Sroots = [Complex.(Sroots.^n) for n in range(1, npeaks)]
        #Croots = roots(Polynomial([0.9893, 0.0107]))
        Croots = [-92.45794392523365]
        Croots = [Complex.(Croots.^n) for n in range(1, npeaks)]
        #Hroots = roots(Polynomial([0.999885, 0.000115]))
        Hroots = [-8694.652173913044]
        Hroots = [Complex.(Hroots.^n) for n in range(1, npeaks)]
        #Oroots = roots(Polynomial([0.99757, 0.00038, 0.00205]))
        Oroots = [ -0.09268292682926828 - 22.05925932732548im
                -0.09268292682926828 + 22.05925932732548im]
        Oroots = [Complex.(Oroots.^n) for n in range(1, npeaks)]
        #Nroots = roots(Polynomial([0.99636, 0.00364]))
        Nroots = [-273.72527472527474]
        Nroots = [Complex.(Nroots.^n) for n in range(1, npeaks)]
        new(Sroots, Croots, Hroots, Oroots, Nroots, UInt32(npeaks))
    end
end

mutable struct Composition
    C::Int32
    H::Int32
    N::Int32
    O::Int32
    S::Int32
end

import Base.+

+(y::Composition, x::Composition) = Composition(x.C + y.C, x.H + y.H, x.N + y.N, x.O + y.O, x.S + y.S)

Composition() = Composition(zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32))

struct Isotope{T<:AbstractFloat} <: Ion{T}
    mass::T
    intensity::T
    prec_idx::UInt32
end
#test 
getMZ(i::Isotope{T}) where {T<:AbstractFloat} = i.mass
getIntensity(i::Isotope{T}) where {T<:AbstractFloat} = i.intensity
getPrecID(i::Isotope{T}) where {T<:AbstractFloat} = i.prec_idx
Isotope{Float32}() = Isotope(zero(Float32), zero(Float32), zero(UInt32))

function getMonoMass(comp::Composition,charge::I) where {I<:Integer}
    return (comp.C*12.000000 + 
           comp.H*1.007825 +
           comp.N*14.0030740052 +
           comp.O*15.99491462219 + 
           comp.S*31.97207069 + PROTON*charge)/charge
end



function getIsotopes(comp::Composition, roots::QRoots, npeaks::Int, charge::I, prec_id::J; precision::DataType = Float32) where {I,J<:Integer}
    npeaks = min(npeaks, roots.npeaks)
    Ψ = zeros(ComplexF64, npeaks - 1)
    for n in range(1, npeaks - 1)
        Ψ[n] = sum(comp.C./(roots.C[n])) + sum(comp.H./(roots.H[n])) + sum(comp.N./(roots.N[n])) + sum(comp.O./(roots.O[n])) + sum(comp.S./(roots.S[n]))
    end
    q = zeros(Float32, npeaks)
    #C, 
    q[1] = prod([0.9893, 0.999885, 0.99757, 0.99636, 0.9499].^[comp.C, comp.H, comp.N, comp.O, comp.S])

    for i in range(2, npeaks)
        qΨ = 0.0
        for k in range(1, i-1)
            qΨ += q[i - k]Ψ[k]
        end
        q[i] = -(1/(i - 1))*qΨ
    end

    mass = Float32(getMonoMass(comp, charge))
    isotopes = Vector{Isotope{precision}}(undef, npeaks)
    for i in eachindex(isotopes)
        isotopes[i] = Isotope(mass, q[i], UInt32(prec_id))
        mass += Float32((NEUTRON/charge))
    end
    return isotopes
end


const aa_to_composition::Dict{Char, Composition} = Dict{Char, Composition}(
        'A' => Composition(3, 5, 1, 1, 0),
        'R' => Composition(6, 12, 4, 1, 0),
        'N' => Composition(4, 6, 2, 2, 0),
        'D' => Composition(4, 5, 1, 3, 0),
        'C' => Composition(3, 5, 1, 1, 1),
        'E' => Composition(5, 7, 1, 3, 0),
        'Q' => Composition(5, 8, 2, 2, 0),
        'G' => Composition(2, 3, 1, 1, 0),
        'H' => Composition(6, 7, 3, 1, 0),
        'I' => Composition(6, 11, 1, 1, 0),
        'L' => Composition(6, 11, 1, 1, 0),
        'K' => Composition(6, 12, 2, 1, 0),
        'M' => Composition(5, 9, 1, 1, 1),
        'F' => Composition(9, 9, 1, 1, 0),
        'P' => Composition(5, 7, 1, 1, 0),
        'S' => Composition(3, 5, 1, 2, 0),
        'T' => Composition(4, 7, 1, 2, 0),
        'W' => Composition(11, 10, 2, 1, 0),
        'Y' => Composition(9, 9, 1, 2, 0),
        'V' => Composition(5, 9, 1, 1, 0),
        'U' => Composition(3, 5, 1, 1, 0),
        'O' => Composition(12, 19, 3, 2, 0)
        )
const mod_to_composition::Dict{String, Composition} = Dict{String, Composition}(
    "ox" => Composition(0, -1, 0, 1, 0)
)
function getElementalComposition(seq::String)
    comp = Composition()
    in_mod = false
    #for AA in eachindex(seq)
    AA = 1
    while AA <= length(seq)
        if seq[AA] == '('
            mod = ""
            i = AA + 1
            while seq[i] != ')'
                mod *= seq[i]
                i += 1
            end
            comp += mod_to_composition[mod]
            AA = i + 1
            continue
        else
            comp += aa_to_composition[seq[AA]]
        end 
        AA += 1
    end
    #Comp + water/H2O
    return comp +  Composition(zero(UInt32), UInt32(2), zero(UInt32), UInt32(1), zero(UInt32))
end

function getIsotopes(seqs::Vector{String}, ids::Vector{I}, charges::Vector{J}, roots::QRoots, npeaks::Int; precision::DataType = Float32) where {I,J<:Integer}
    isotopes = UnorderedDictionary{UInt32, Vector{Isotope{precision}}}()
    for i in eachindex(seqs)
        insert!(isotopes, ids[i], getIsotopes(getElementalComposition(String(seqs[i])), roots, npeaks, charges[i], ids[i], precision = precision))
    end
    isotopes
end
