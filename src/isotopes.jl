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

struct Isotope{T<:AbstractFloat} <: IonType
    mass::T
    intensity::T
    prec_idx::UInt32
end

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

sequence = "M(ox)QVDQEEPHTEEQQQQPQTPAENK"
for match in findall(r"[A-Z](\(.*?\))+", sequence)
    println(match)
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
#=
#=
        Sroots = roots(Polynomial([0.9499, 0.0075, 0.0425, 0.0, 0.0001]))
        Croots = roots(Polynomial([0.9893, 0.0107]))
        Hroots = roots(Polynomial([0.999885, 0.000115]))
        Oroots = roots(Polynomial([0.99757, 0.00038, 0.00205]))
        Nroots = roots(Polynomial([0.99636, 0.00364]))
Ψ1 = Complex(sum(115 ./Croots)) + sum(4 ./Sroots) + sum(47 ./ Oroots) + Complex(sum(35 ./ Nroots)) + Complex(sum(2 ./Hroots))
Ψ2 = Complex(sum(115 ./(Croots.^2))) + sum(4 ./(Sroots.^2)) + sum(47 ./(Oroots.^2)) + Complex(sum(35 ./ (Nroots.^2))) + Complex(sum(2 ./Hroots.^2))
Ψ3 = Complex(sum(115 ./Croots.^3)) + sum(4 ./Sroots.^3) + sum(47 ./ Oroots.^3) + Complex(sum(35 ./ Nroots.^3)) + Complex(sum(2 ./Hroots.^3))
Ψ4 = Complex(sum(115 ./Croots.^4)) + sum(4 ./Sroots.^4) + sum(47 ./ Oroots.^4) + Complex(sum(35 ./ Nroots.^4)) + Complex(sum(2 ./Hroots.^4))
Ψ5 = Complex(sum(115 ./Croots.^5)) + sum(4 ./Sroots.^5) + sum(47 ./ Oroots.^5) + Complex(sum(35 ./ Nroots.^5)) + Complex(sum(2 ./Hroots.^5))
q = zeros(6)
q0 = prod([0.9499,  0.9893, 0.99757, 0.99636, 0.999885].^[4, 115, 47, 35, 0])
q1 = -(1/1)*q0*Ψ1
q2 = -(1/2)*(q1*Ψ1 + q0*Ψ2)
q3 = -(1/3)*(q2*Ψ3 + q1*Ψ2 + q0*Ψ1)
q4 = -(1/4)*(q3*Ψ4 + q2*Ψ3 + q1*Ψ2 + q0*Ψ1)
q5 = -(1/5)*(q4*Ψ5 + q3*Ψ4 + q2*Ψ3 + q1*Ψ2 + q0*Ψ1)

q = zeros(6) 
q[1] = prod([0.9499,  0.9893, 0.99757, 0.99636].^[4, 115, 47, 35])
q[2] = -(1/1)*q[1]*Ψ1
q[3] = -(1/2)*(q[2]*Ψ1 + q[1]*Ψ2)
q[4] = -(1/3)*(q[3]*Ψ1 + q[2]*Ψ2 + q[1]*Ψ3)
q[5] = -(1/4)*(q[4]*Ψ1 + q[3]*Ψ2 + q[2]*Ψ3 + q[1]*Ψ4)
q[6] = -(1/5)*(q[5]*Ψ1 + q[4]*Ψ2 + q[3]*Ψ3 + q[2]*Ψ4 + q[1]*Ψ5)
const AA_to_composition::Dict{Char, Composition} = Dict{Char, Composition}(
        'A' => Composition(2, 3, 1, 1, 0),
        'R' => Composition(6, 12, 4, 1, 0),
        'N' => Composition(4, 6, 2, 2, 0),
        'D' => Composition(4, 5, 1, 3, 0),
        'C' => Composition(3, 5, 1, 1, 1),
        'E' => Composition(5, 7, 1, 3, 0),
        'Q' => Composition(5, 8, 2, 2, 0),
        'G' => Composition(2, 3, 1, 0, 0),
        'H' => Composition(6, 7, 3, 1, 0),
        'I' => Composition(6, 11, 1, 1, 0),
        'L' => Composition(6, 11, 1, 1, 0),
        'K' => Composition(6, 12, 2, 1, 0),
        'M' => Composition(5, 9, 1, 1, 1),
        'F' => Composition(6, 9, 1, 1, 0),
        'P' => Composition(5, 7, 1, 1, 0),
        'S' => Composition(3, 5, 1, 2, 0),
        'T' => Composition(4, 7, 1, 2, 0),
        'W' => Composition(11, 10, 2, 1, 0),
        'Y' => Composition(6, 9, 1, 2, 0),
        'V' => Composition(5, 9, 1, 1, 0),
        'U' => Composition(3, 5, 1, 1, 0),
        'O' => Composition(12, 19, 3, 2, 0),
        )
import Base.+
+(y::Composition, x::Composition) = Composition(x.C + y.C, x.H + y.H, x.N + y.N, x.O + y.O, x.S + y.S)
mutable struct Composition
    C::UInt32
    H::UInt32
    N::UInt32
    O::UInt32
    S::UInt32
end
Composition() = Composition(zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32), zero(UInt32))
function getAAComposition(seq::String)
    comp = Composition()
    for AA in eachindex(seq)
        comp += AA_to_composition[seq[AA]]
    end
    return comp +  Composition(zero(UInt32), UInt32(2), zero(UInt32), UInt32(1), zero(UInt32))
end
C34H51N7O14
Int64(A.C)
Int64(A.H)
Int64(A.N)
Int64(A.O)
Int64(A.S)
isotopic_distribution("C115H179N35O47S1", 0.90)

[0.263097, 0.20978, 0.163548]./[0.225466, 0.179775,  0.140155]

DataFrame(a[2:end,:], a[1,:])
struct IsotopeProbsGenerator
    C::Binomal{Float64}
    S::Binomal{Float64}
end

struct Isotope{T<:AbstractFloat}
    mass::T
    prob::T
end

struct IsotopeProbs
    C::Tuple{Isotope{T},Isotope{T},Isotope{T}}
    S::Tuple{Isotope{T},Isotope{T},Isotope{T}}
end

function getIsotopes
    (
        C[1]*S[1],
        C[2]*S[1] + C[1]*S[2],
        C[3]*S[1] + C[1]*S[3] + C[2]*S[2],

    )
end


function getIsotopes
    (
        C[1]*S[1],
        C[2]*S[1],
        C[1]*S[2],
        C[3]*S[1],
        C[1]*S[3],
        C[2]*S[2],
        C[4]

    )
end

#=(
C[1]*N[1]*O[1]*S[1],

C[2]*N[1]*O[1]*S[1] + C[1]*N[2]*O[1]*S[1] + C[1]*N[2]*O[1]*S[1] + C[1]*N[2]*O[1]*S[1],

C[2]*(N[2]*O[1]*S[1] + N[1]*O[2]*S[1] + N[1]*O[1]*S[2]) +
N[2]*(C[1]*O[2]*S[1] + C[1]*O[1]*S[2]) + 
C[3]*N[1]*O[1]*S[1] + C[1]*N[3]*O[1]*S[1] + C[1]*N[1]*O[3]*S[1] + C[1]*N[1]*O[1]*S[3],

C[2]*N[2]*(O[1]*S[2] + O[2]*S[2]) + 
C[2]*O[2]*(S[1]*N[2] + S[2]*N[1]) + 
C[2]*S[2]*(O[1]*N[2] + O[2]*N[1]) + 

C[3]*(N[2]*O[1]*S[1] + N[1]*O[2]*S[1] + N[1]*O[1]*S[2]) + 
S[3]*(N[2]*O[1]*S[1] + N[1]*O[2]*S[1] + N[1]*O[1]*S[2]) + 
N[3]*(N[2]*O[1]*S[1] + N[1]*O[2]*S[1] + N[1]*O[1]*S[2]) + 
O[3]*(N[2]*O[1]*S[1] + N[1]*O[2]*S[1] + N[1]*O[1]*S[2])

)=#

"""
Isotopic distributions. 
"""

#include("isotopes-data.jl")

# User Interface.
# ---------------


export isotopic_distribution, masses, simulate, formula


"""
    masses(input::String)
Calculates the average, monoistopic and nominal masses for the chemical formula given as an input. The result is returned in a dictionary with the following entries: "Monoiotopic", Average" and "Nominal".
# Examples
```julia-repl
julia> masses("C254 H377 N65 O75 S6")
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
Dict{String,Float64} with 3 entries:
  "Monoisotopic" => 5729.6
  "Average"      => 5733.55
  "Nominal"      => 5727.0
```
"""
function masses(input::String)
    f = formula(input)
    masses(f)
end

"""
    masses(f::Dict{String,Int})
Calculates the average, monoistopic and nominal masses for the chemical formula dictionary, such as prodcued by MassJ.formula. The result is returned in a dictionary with the following entries: "Monoiotopic", Average" and "Nominal".
# Examples
```julia-repl
julia> masses("C254 H377 N65 O75 S6")
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
Dict{String,Float64} with 3 entries:
  "Monoisotopic" => 5729.6
  "Average"      => 5733.55
  "Nominal"      => 5727.0
```
"""
function masses(f::Dict{String,Int})
    m = Dict{String, Float64}()

    # exact mass
    avg = 0.0
    mono = 0.0
    nomi = 0.0
    for (key,val) in f
        avg += AverageMasses[key] * val
        mono += first(Elements[key]).m * val
        nomi += Int( map(round,first(Elements[key]).m) ) * val
    end
    m["Average"] = avg
    m["Monoisotopic"] = mono
    m["Nominal"] = nomi
    m
end


"""
    simulate(I::Array{Union{Float64, Int, String}}, ∆mz::Real; model::Symbol=:gauss, Npoints::Int=1000)
From an isotopic distribution and a peak width returns a mass spectrum (MSScan). The number of points of the resulting mass spectrum is passed as an optional argument. Peak shape are :gauss (default), :lorentz, :voight.
# Examples
```julia-repl
julia>  a = simulate(I, 0.4)
MassJ.MSscan(1, 0.0, 30898.192348114364, [5727.102517458742 ..., "", "", 0.0)
```
"""
function simulate(I::Array{Union{Float64, Int, String}}, ∆mz::Real; model::Symbol=:gauss , Npoints::Int=1000)
    r =  maximum(I[2:end,1]) - minimum(I[2:end,1])
    start =  minimum(I[2:end,1]) - r/2.
    stop  =  maximum(I[2:end,1]) + r/2.
    r = stop - start
    
    p  = zeros(Npoints)
    mz = zeros(Npoints)
    peak = zeros(Npoints)  
    
    mz = [start + (x * r / length(mz)) for x=1:length(mz)]
    
    for (i,val) in enumerate(I[2:end, 1])
        ind = num2pnt(mz, val)
        x0 = mz[ind]
        y0 = I[2:end,2][i]
        param = [∆mz, x0, y0, 0.0]
        if model == :gauss
            peak = [gauss(x,param) for x in mz]
        elseif model == :lorentz
            peak = [lorentz(x,param) for x in mz]
        elseif model == :voigt
            peak = [voigt(x,param) for x in mz]
        else
            ErrorException("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
        end            
        p = p .+ peak
    end
    p /= maximum(p)
    p *= 100.0
    basePeakMz = mz[num2pnt(p, 100.0)]
    MSscan(1, 0.0, sum(p), mz, p, 0, basePeakMz, 100.0, 0.0, "", "", 0.0)
end

"""
    isotopic_distribution(input::String, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{MassJ.Isotope,1}} = MassJ.Elements)
Calculates the isotopic distribution of input formula for which the overall probabilities equals p_target using the isospec algorithm. The charge state is entered as an optional argument. The peaks detection threshold tau is by default set to 10%.
# Examples
```julia-repl
julia> isotopic_distribution("C254 H377 N65 O75 S6", 0.5)
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
9×15 Array{Union{Float64, Int, String},2}:
     "Masses"   "Probability"     "12C"   "13C"   "32S"   "34S"   "33S"   "36S"    "14N"   "15N"    "16O"   "18O"   "17O"     "1H"   "2H"
 5731.61       0.112302        252       2       6       0       0       0       65       0       75       0       0       377      0    
 5732.61       0.102878        251       3       6       0       0       0       65       0       75       0       0       377      0    
 5730.6        0.0814037       253       1       6       0       0       0       65       0       75       0       0       377      0    
 5733.61       0.0704028       250       4       6       0       0       0       65       0       75       0       0       377      0    
 5734.62       0.0383896       249       5       6       0       0       0       65       0       75       0       0       377      0    
 5733.6        0.0301637       252       2       5       1       0       0       65       0       75       0       0       377      0    
 5729.6        0.0293871       254       0       6       0       0       0       65       0       75       0       0       377      0    
 5734.61       0.0276323       251       3       5       1       0       0       65       0       75       0       0       377      0    
``julia> isotopic_distribution("C254 H377 N65 O75 S6", 0.5, charge = +7)
Dict("S" => 6,"C" => 254,"N" => 65,"H" => 377,"O" => 75)
9×15 Array{Union{Float64, Int, String},2}:
    "Masses"   "Probability"     "12C"   "13C"   "32S"   "34S"   "33S"   "36S"    "14N"   "15N"    "16O"   "18O"   "17O"     "1H"   "2H"
 818.801      0.112302        252       2       6       0       0       0       65       0       75       0       0       377      0    
 818.944      0.102878        251       3       6       0       0       0       65       0       75       0       0       377      0    
 818.658      0.0814037       253       1       6       0       0       0       65       0       75       0       0       377      0    
 819.088      0.0704028       250       4       6       0       0       0       65       0       75       0       0       377      0    
 819.231      0.0383896       249       5       6       0       0       0       65       0       75       0       0       377      0    
 819.086      0.0301637       252       2       5       1       0       0       65       0       75       0       0       377      0    
 818.514      0.0293871       254       0       6       0       0       0       65       0       75       0       0       377      0    
 819.23       0.0276323       251       3       5       1       0       0       65       0       75       0       0       377      0    
```
"""
function isotopic_distribution(input::String, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{Isotope,1}} = Elements)
    f = formula(input)
    println(f)
    isotopic_distribution(formula(input), p_target, charge = charge, tau = tau)
end

function isotopic_distribution(form::Dict{String,Int}, p_target::Real; charge::Int = +1, tau::Real = 0.1, Elements::Dict{String,Array{Isotope,1}} = Elements)
    c_form = collect(form)
    sort_formula(x) = 
        begin 
            if length(Elements[x]) > 1
                return Elements[x][2].f * form[x]
            else
                return 0.0
            end
        end
    sort!(c_form, by = x -> sort_formula(x[1]), rev=true)    
    I = isospec(form, c_form, p_target, tau, Elements)    
    cI = collect(I)
    
    N = 0
    for (key, val) in c_form
        N += length(Elements[key])
    end
    N += 2                              # columns
    M = length(I)                       # rows
    result = Matrix{Union{String, Float64, Int}}(undef, M+1,N)
    labels = Vector{String}(undef, N)
    labels[1] = "Masses"
    labels[2] = "Probability"
    for i=1:M
        temp_dist = Dict{String,Array{Int,1}}()
        j = 1
        for k = 1:length(c_form)
            elem = c_form[k][1]
            temp_vect = [ cI[i][1][j:j+length(Elements[elem])-1]... ]
            result[i,j+2:j+2+length(Elements[elem])-1] = temp_vect
            j = j + length(Elements[elem])
            temp_dist[elem] = temp_vect
        end
        masses =  isotopologue_mass(c_form, temp_dist, Elements)
        result[i,1] = masses / abs(charge)
        result[i,2] = cI[i][2]
    end
    j = 3
    for i=1:length(c_form)
        elem = c_form[i][1]
        isot = length(Elements[elem])
        for m=1:length(Elements[elem])
            labels[j] = string(Elements[elem][m].A) * elem
            j+=1
        end
    end
    result[M+1,:] = labels
    result[end:-1:1,:]
end


"""
    formula(formula::String)
Private function that reads the input chemical formula and sorts the atoms. It returns a dictionary in which the different entries represent the atoms and the values are the number of times the atoms have been found in the formula.
"""
function formula(formula::String)
    if !occursin(r"[A-Z-(]", formula[1:1])                                        # is Formula starting with capital ?
        error("Invalid chemical formula")
    end    
    form_dict = Dict{String,Int}() 
    i = 1
    j = 1
    while i <= length(formula)
        if occursin(r"[A-Z]", formula[i:i])                                     # is formula[i] is a capital
            temp = formula[i:i]
            prev = (i:i)
            if i < length(formula) 
                if occursin(r"[a-z]", formula[i+1:i+1])                         # is the capital is followed by a lower-case ? 
                    temp = formula[i:i+1]                 
                    prev = (i:i+1)
                end
            end
            if haskey(Elements, temp)                                   # valid element --> looking for numbers 
                if last(prev) == length(formula)                                # last element
                    if !haskey(form_dict, formula[prev])                           # already found ?
                        form_dict[formula[prev]] = 1                               # no --> indice = 1
                    else
                        form_dict[formula[prev]] = form_dict[formula[prev]] + 1 # yes --> add 1 to indice
                    end 
                elseif last(prev) < length(formula)                             # not last element
                    if occursin(r"[0-9]", formula[last(prev):last(prev)+1])        # indice after element ?               
                        next = findnext(r"\d+", formula, last(prev)+1)
                        if !haskey(form_dict, formula[prev])                          # already found ?
                            form_dict[formula[prev]] = parse(Int,formula[next])  
                        else 
                            form_dict[formula[prev]] = form_dict[formula[prev]] + parse(Int,formula[next])
                        end
                    else                                                       # no indice after element
                        if !haskey(form_dict, formula[prev])                          # already found ?
                            form_dict[formula[prev]] = 1                                # no --> indice = 1
                        else
                            form_dict[formula[prev]] = form_dict[formula[prev]] + 1 
                        end    
                    end    
                end
            else
                error("Invalid chemical element: $temp")
            end
        elseif occursin(r"[(]", formula[i:i])                                  # formula contains a ( --> isotopes or group
            j = i+1
            mult = 1
            l = j
            while !occursin(r"[)]", formula[l+1:l+1])         # find closing ) and check for an indice --> mult factor
                l+=1
            end
            if l+2 > length(formula)                       # if ) is the last char of the formula, mult = 1
                mult = 1
            else
                if occursin(r"[1-9]", formula[l+2:l+2])   #  indice found after )
                    k = l+2
                    if k < length(formula)
                        while occursin(r"[0-9]", formula[k:k])
                            k += 1
                        end
                        mult = parse(Int,formula[l+2:k-1])
                    else
                        mult = parse(Int,formula[l+2:k])
                    end
                else
                    mult = 1
                end
            end
            #println("mult :$mult")
            while !occursin(r"[)]", formula[j:j])       # repeat until the end
                l = j
                if occursin(r"[0-9]", formula[j:j])                      # number found after ( --> isotope
                    while occursin(r"[0-9]", formula[l:l])
                        l += 1
                    end                    
                    prev = formula[j:l]
                    if haskey(Elements, prev)               # element found
                        j=l
                        if occursin(r"[)]", formula[j+1:j+1])        # no indice after element --> indice = 1
                            next = "1"
                        elseif occursin(r"[1-9]", formula[j+1:j+1])
                            l = j+1
                            while occursin(r"[0-9]", formula[l:l])
                                l+=1
                            end
                            next = formula[j+1:l-1]
                            j = l-1
                        elseif occursin(r"[A-Z]", formula[j+1:j+1])    # found following element --> storing previous one
                            next = "1"
                        end
                        if !haskey(form_dict, prev)
                            form_dict[prev] = parse(Int,next) * mult                                # no --> indice = 1
                        else
                            form_dict[prev] = form_dict[prev] +  (parse(Int,next) * mult) 
                        end                                
                    else
                        error("Invalid chemical element: $prev")
                    end
                elseif occursin(r"[A-Z]", formula[j:j])               # letter found after ( --> chemical group
                    l = j               
                    prev = formula[j:j]
                    next = ""
                    if occursin(r"[a-z]", formula[j+1:j+1])                 # is the capital is followed by a lower-case ? 
                        prev = prev * formula[j+1:j+1]
                        j += 1
                        println(prev)
                        if occursin(r"[)]", formula[j+1:j+1])               # capital letter followed by ) indice -> 1
                            next = "1"
                        elseif occursin(r"[1-9]", formula[j+1:j+1])            # finding indices
                            l = j+1
                            while occursin(r"[0-9]", formula[l:l])
                                l += 1
                            end
                            next = formula[j+1:l-1]
                            j = l-1
                        elseif  occursin(r"[A-Z]", formula[j+1:j+1])        # is followed by another element
                            next = "1"
                        end
                    elseif occursin(r"[1-9]", formula[j+1:j+1])            # finding indices
                        l = j+1
                        while occursin(r"[0-9]", formula[l:l])
                            l += 1
                        end
                        next = formula[j+1:l-1]
                        j = l-1
                    elseif occursin(r"[A-Z]", formula[j+1:j+1])            # found another element next
                        next = "1"
                    elseif occursin(r"[)]", formula[j+1:j+1])               # last element before )
                        next = "1"
                    end
                    if haskey(Elements, prev)                             # element found

                        if !haskey(form_dict, prev)                           # already found ?
                            form_dict[prev] = parse(Int,next) * mult         # no --> indice
                        else
                            form_dict[prev] = form_dict[prev] +  (parse(Int,next) * mult) # yes --> add 1 to indice
                        end
                    else
                        error("Invalid chemical element: $prev")
                    end
                end
                j+=1
            end
            i = j
        end
        i += 1
    end

    form_dict
end



function isotopologue_probability(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{Isotope,1}})
    p = 1.0                                                                                        # isotopologue probability
    Natoms = sum(collect(values(formula)))
    if Natoms < 20
        return low_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements)
    else

        return high_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements)
    end
    
end

function low_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{Isotope,1}})
    p = 1.0
    for (key, val) in formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        
        for i=1:length(distri[key])
            p = p * E[i].f^distri[key][i]
        end
        # getting coefficients
        coefs = factorial(BigInt(val)) 
        for el in distri[key]
            coefs /= factorial(BigInt(el))
        end
        p *= coefs
    end
    p = Float64(p)                            # returns mass and probability for the input isotopologue
end

function high_masses(formula::Dict{String,Int}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{Isotope,1}})
    # log (factorial(val) ) - Somme log(factorial(distri[key])) + somme(distri[key] * log(E[i].f) )
    log_p = 0.0
    for (key,val) in formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        
        for i=1:length(distri[key])
            log_p = log_p + distri[key][i]*E[i].logf
        end
        if val <= 256
            log_p = log_p + log_factorial[val]
        else
            log_p  = log_p + stirling(val)
        end
        for i in eachindex(distri[key])
            el = distri[key][i]
            if el != 0
                if el <= 256
                    log_p = log_p - log_factorial[el]
                else
                    log_p = log_p - stirling(el)
                end
            end
            
        end
    end
    exp(log_p)
end

function stirling(n::Int)
    x = n + 1
    return (x - 0.5)*log(x) - x + 0.5*log(2*pi) + 1.0/(12.0*x)
end


function isotopologue_mass(c_formula::Array{Pair{String,Int},1}, distri::Dict{String,Array{Int,1}}, Elements::Dict{String,Array{Isotope,1}})
    m = 0.0                                      # isotopologue mass
    for (key, val) in c_formula
        E = Elements[key]
        if sum(distri[key]) != val
            error("Distribution does not match chemical formula")
        end
        for i=1:length(distri[key])
            m = m + E[i].m * distri[key][i]
        end
    end
    m                                      # returns the mass of the input isotopologue
end


function most_probable_isotopologue(formula::Dict{String,Int}, Elements::Dict{String,Array{Isotope,1}})
    distri=Dict{String,Array{Int,1}}()        
    # enter hill climbing
    for (key, val) in formula
        E = Elements[key]
        if length(E) > 1
            # # start with the mean of the multimodal distribution
            P0 = Int.(map(round, val .* [el.f for el in E]))
            counts = 1
            while sum(P0) != val               # checking for mismatch between distri & formula
                for i=length(P0):-1:1
                    if P0[i] > 0               # found a non negative element                 
                        if sum(P0) > val       
                            P0[i] -= 1
                        elseif sum(P0) < val   
                            P0[i] +=1
                        end
                    end
                end
                if counts < 10
                    break
                end
                counts += 1
            end
            # Objective function isotopologue_probability(Dict(key => val), x::AbstractArray)
            f(x) = isotopologue_probability(Dict(key => val), Dict(key => x), Elements)
            
            # call hill_climbing(P, f) --> returns P
            P = hill_climbing(P0, f)
            distri[key] = P
        else
            # element with single isotopes
            distri[key] = Int.(val .* [el.f for el in E])
        end
    end
    distri
end


function hill_climbing(P::AbstractArray, f::Function)
    max_count = 5000
    counts = 0
    Pabs = P
    Pnext = P
    # find all neighbours
    n = fill(0, length(P))
    if length(P) > 1
        n[1] = -1
        n[2] = +1
        neighbours = unique(permutations(n) )
    end
    imax = 0
    while counts < max_count
        counts += 1
        for i in eachindex(neighbours)   # find neighbours that increase the function f
            temp = P .+ neighbours[i] 
            if !in(-1, temp)                 # Pnext has no negative values
                Pnext = temp
                if f(Pnext) >= f(P)     # Pnext is better than Plocal
                    imax = i                   
                end
            end
        end
        if imax != 0
            if !in(-1, P .+ neighbours[imax]  )
                P = P .+ neighbours[imax]       # move to the best neighbour
            end
        end
        
        if f(P) > f(Pabs)         
            Pabs = P
        else
            break
        end
            
    end
    if counts >= max_count
        alert("Couldn't find the most probable subisotopologue.")
    end
    return Pabs
end


function subgenerator!(formula::Dict{String,Int}, D_max::Dict{String,Array{Int,1}},  V::PriorityQueue, τ::Real, Elements::Dict{String,Array{Isotope,1}})
    key = first(formula)[1]
    val = first(formula)[2]
    E = Elements[key]
    P_max = D_max[key]
    fP_max = isotopologue_probability(formula, D_max, Elements)
    
    f(x) = isotopologue_probability(Dict(key => val), Dict(key => x), Elements)
    
    if isempty(V)
        V[P_max] = f(P_max)
    end
    PQ = deepcopy(V)
    n = fill(0, length(P_max))
    if length(P_max) > 1                                       # ensure there at least 2 isotopes for this element
        n[1] = -1
        n[2] = +1
        neighbours = unique( permutations(n) )
        while !isempty(PQ)
            beta = dequeue!(PQ)
            for i in eachindex(neighbours)
                temp = beta .+ neighbours[i]
                if !in(-1, temp)                       
                    Pnext = temp
                    if !haskey(V, Pnext)                      # if n not in V
                        fPnext = f(Pnext)
                        if fPnext >= (τ * fP_max)          # if prob above threshold
                            V[Pnext] = fPnext                 # then keep configuration
                            PQ[Pnext] = fPnext
                        end
                    end
                end
            end
        end
    end
    V
end

function ordered_isotopologue(c_formula::Array{Pair{String,Int},1}, α::Dict{String,Array{Int,1}}, V::Array{DataStructures.PriorityQueue,1}, I::PriorityQueue, tau::Real, Elements::Dict{String,Array{Isotope,1}} ) 
    S = Array{PriorityQueue}(undef,0)

    i = 1
    for (key, val) in c_formula
        form = Dict(key => c_formula[i][2])
        temp = Dict(key => α[key])
        subgenerator!(form, temp, V[i], tau, Elements)
        push!(S, V[i])
        i += 1
    end
    I = S[1]
    for i=2:length(S)
        cI = collect(I)
        I = PriorityQueue()
        cSi = collect(S[i])
        for el in cI
            for j=1:length(S[i])
                conf = (el[1]..., cSi[j][1]...)
                prob = el[2] * cSi[j][2]
                I[conf] = prob
            end
        end
    end
    V, I
end

function trim!( I::PriorityQueue, p_target::Real, summ::Real)
    next = peek(I)
    somme = 0.0
    while true
        if !isempty(I)
            somme = sum( collect( values(I) )[2:end] )
        end
        if somme <= p_target
            break
        end
        if isempty(I)
            break
        else
            dequeue!(I)
        end       
    end
    I
end


function isospec(formula::Dict{String,Int}, c_formula::Array{Pair{String,Int},1}, p_target::Real, tau::Real, Elements::Dict{String,Array{Isotope,1}} )
    α = most_probable_isotopologue( formula, Elements )
    I = PriorityQueue()
    N_elem = length(formula)
    V = Vector{PriorityQueue}(undef,0)
    for i =1:N_elem
        push!(V, PriorityQueue())
    end
    count = 1
    summ = 0.0
    while true
        s = @sprintf "%i \t tau: %e" count tau                             # used for debugging
        #print(s)                                                          # used for debugging
        
        V,I = ordered_isotopologue(c_formula, α, V, I, tau, Elements)

        summ = 0.0
        for (key,val) in I
            summ += val
        end
        u = @sprintf "\tprob.: %f \t diff.: %e" summ (p_target-summ)        # used for debugging
        #println(u)                                                         # used for debugging

        a = (summ - 1.0)/tau
        tau = (p_target-1.0) / a
        
        count += 1
        if count >= 5
            break
        end
        if summ >= p_target
            break
        end
    end
    trim!(I, p_target, summ)
    I
end


=#
=#