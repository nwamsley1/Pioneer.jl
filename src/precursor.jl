##########
#Global Constants
##########
const H2O = Float32(18.010565)
const PROTON = Float32(1.0072764)
const NEUTRON = Float32(1.00335)
const AA_to_mass =
Dict{Char, Float32}(
        'A' => 71.03711,
        'R' => 156.10111,
        'N' => 114.04293,
        'D' => 115.02694,
        'C' => 103.00919,
        'E' => 129.04259,
        'Q' => 128.05858,
        'G' => 57.02146,
        'H' => 137.05891,
        'I' => 113.08406,
        'L' => 113.08406,
        'K' => 128.09496,
        'M' => 131.04049,
        'F' => 147.06841,
        'P' => 97.05276,
        'S' => 87.03203,
        'T' => 101.04768,
        'W' => 186.07931,
        'Y' => 163.06333,
        'V' => 99.06841,
        'U' => 150.95363,
        'O' => 237.14773
)

const default_mods = 
Dict{String, Float32}(
    "Carb" => 57.021464
)
##########
#Representation of Amino Acid
##########
struct AA 
    aa::Char
    mass::Float32
    #Constructor for amino acid. Restric inputs
    #to valid amino acid symbols and assign correct mass
    function AA(aa::Char)
        m = try
            AA_to_mass[aa]
        catch
            throw(ErrorException("The character $aa cannot be interpreted as an amino acid!"))
        end
        return new(aa, m)
    end
end

#Getter methods
getMass(aa::AA) = aa.mass
getAA(aa::AA) = aa.aa

export AA
##########
#Modification. Simple representation of mass modifiction
##########
struct Mod
    mass::Float32
end

function Mod(mod::String, mods_dict::Dict{String, Float32})
    """
    Given a string, first parse by the regular expression  
        Example: "K[+8.014199]" or "C[Carb]"
        
    In the first case, "K[+8.014199]", "K[+8.014199]" is the modification name
        and 8.014199 is the modification mass. 

    In the second case, "C[Carb]" is the modification name
        and "Carb" is a key to the dictionary "mods_dict".
        mods_dict["Carb"] returns the modification mass. 

    If the Mod string can't be parsed, returns an error. 
    """
    m = match(r"^[A-Z]\[(.*)\]$", mod)

    try
        if m == nothing
            Mod(0.0)
        elseif startswith(m[1], "+")
            Mod(parse(Float32, m[1][2:end])) #"K[+8.014199]"
        else 
            Mod(mods_dict[m[1]]) #getAA("C[Carb]")
        end
    catch
        throw(ErrorException("$mod could not be parsed as given"))
    end 
end

#Optionally parse mods without a mods_dict
Mod(mod::String) = Mod(0.0)
Mod(name::Char) = Mod(0.0)
Mod() = Mod(0.0)
#Getter Functions
getMass(mod::Mod) = mod.mass    
export Mod

##########
#Residue. Implementation of amino acid with custom mass modifications
##########
struct Residue
    mass::Float32
end

#Residue(AA('A'))
function Residue(aa::AA)
    Residue(getMass(aa))
end

#Residue('A')
Residue(aa::Char) = Residue(AA(aa))

Residue(aa::AA, mod::Mod) = Residue(getMass(mod)+getMass(aa))

Residue(residue::String, mods_dict::Dict{String, Float32}) = Residue(AA(residue[1]), Mod(residue, mods_dict))

Residue(residue::String) = Residue(residue, Dict{String, Float32}())

Residue(residue::String, mod_mass::Float32) = Residue(getMass(AA(residue[1])) + mod_mass)

Residue(residue::Char, mod_mass::Float32) = Residue(getMass(AA(residue)) + mod_mass)

function getResidues(sequence::String, mods_dict::Dict{String, Float32})
    map(residue -> Residue(sequence[residue], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", sequence))
end

#Getter methods
getMass(residue::Residue) = residue.mass

export Residue
export getMass
export getMod
export getAA
export default_mods

##########
#Frag, Transition, Precursor
##########
#Abstract Frag type and concrete daugher types Transition and Precursor
abstract type Ion end

struct Transition <: Ion
    frag_mz::Float32
    prec_mz::Float32
    pep_id::Int32
    ion_type::Char
    charge::Int32
    iosotope::Int32
    ind::Int32
end

struct Precursor <: Ion
    mz::Float32
    charge::Int32
    iosotope::Int32
    pep_id::Int32
end

function PrecursorMZ(residues::Array{Residue, 1}, charge::Int32, isotope::Int32)
    (sum(map(residue->getMass(residue), residues))+(PROTON*charge + H2O + isotope*NEUTRON))/charge
end

#Constructor for precursor
function Precursor(residues::Array{Residue, 1}, charge::Int32, isotope::Int32, pep_id::Int32)
    Precursor(
              PrecursorMZ(residues, charge, isotope),
              charge, 
              isotope,
              pep_id
            )
end



#Getter Functions
getFragMZ(ion::Ion) = ion.frag_mz
getPrecMZ(transition::Transition) = transition.prec_mz
getMZ(precursor::Precursor) = precursor.mz
getPepID(ion::Ion) = ion.pep_id
getIonType(transition::Transition) = transition.ion_type
getIonType(precursor::Precursor) = 'p'
getCharge(ion::Ion) = ion.charge
getIsotope(ion::Ion) = ion.isotope

function getFragIons(residues::Array{Residue, 1}, prec_mz::Float32, pep_id::Int32, modifier::Float32, ion_type::Char, start::Int32, charge::Int32, isotope::Int32)
    function __getFragIons__(residues::Array{Residue, 1}, modifier::Float32, charge::Int32, isotope::Int32)
        enumerate((
                    cumsum(
                            map(residue->getMass(residue), residues) .+ (modifier + isotope*NEUTRON)
                          )
                    )/charge)
    end
    map(frag -> FragInd(frag[2], prec_mz, pep_id, ion_type, charge, isotope, Int32(frag[1])+start), __getFragIons__(residues[start:end], modifier, charge, isotope))
end

#Shorthands for getting precursors, b ions, and y ions
#Vector of precursors for given charge states and isotopes
getPrecursors(residues::Array{Residue, 1}, 
              charges::Vector{Int32}, 
              isotopes::Vector{Int32}, 
              pep_id::Int32) = Precursor.(residues, charges, isotopes, pep_id)

getBIons(
        residues::Array{Residue, 1}, 
        prec_mz::Float32, 
        prec_id::Int32, 
        start::Int32,
        charge::Int32, 
        isotope::Int32
        ) = getFragIons(
                        residues,
                        prec_mz,
                        prec_id,
                        PROTON*charge,
                        'b',
                        start,
                        charge,
                        isotope)
getYIons(
        residues::Array{Residue, 1}, 
        prec_mz::Float32, 
        prec_id::Int32, 
        start::Int32,
        charge::Int32, 
        isotope::Int32
        ) = getFragIons(
                        residues,
                        prec_mz,
                        prec_id,
                        PROTON*charge + H2O,
                        'y',
                        start,
                        charge,
                        isotope)


function getFragIons(residues::Array{Residue, 1}, prec_mz::Float32, prec_id::Int32, charge::Int32, isotope::Int32, y_start::Int32, b_start::Int32)
    vcat(
        getBIons(residues, prec_mz, prec_id, b_start, charge, isotope),
        getYIons(residues, prec_mz, prec_id, y_start, charge, isotope)
        )
end

getFragIons(residues::Array{Residue, 1},
            pep_id::Int32, charge::Int32, isotope::Int32, y_start::Int32, b_start::Int32
            ) = getFragIons(residues,
                           PrecursorMZ(residues, charge, isotope),
                           pep_id, charge, isotope, y_start, b_start)


#Get all b and y ions at each charge state and isotope specified. 
getFragIons(residues::Array{Residue, 1}, pep_id::Int32,
            charges::Array{Int32}, 
            isotopes::Array{Int32}, 
            y_start::Int32, b_start::Int32) = getFragIons.(residues, pep_id,
                                                           charges, isotopes, y_start, b_start)

#Same as above for when no isotopes specified
getFragIons(residues::Array{Residue, 1}, 
            pep_id::Int32, 
            charges::Array{Int32}, 
            y_start::Int32, b_start::Int32) = getFragIons.(residues, pep_id, charges, [Int32(0)], y_start, b_start)

#Can provide a list of named tuples to specify exactly which fragments to get


export Frag
export getMZ
export getCharge
export getType
export getIso

##########
#Peptide
##########
abstract type AbstractPeptide end

mutable struct TargetPeptide <: AbstractPeptide
    sequence::String,
    unmodified_sequence::String
    transitions::Array{Transition, 1}
    precursors::Array{Precursors, 1}
    mods::Array{String, 1}
    isotope_label::String
    function TargetPeptide(sequence::String, charge::Int32, isotope_label::String, mods_dict::Dict{String, Float32})
        new(
            sequence, #sequence
            replace(sequence, r"(\[.*?\])"=>"") #unmodified_sequence
            Array{Transition, 1}(), #transitions
            Array{Precursor, 1}(), #precursors
            map(mod -> sequence[mod], findall(r"[A-Z]\[.*?\]", sequence)) #Mods
            isotope_label #isotope_label
        )
    end
end

getSequence(pep::AbsractPeptide) = pep.sequence

function frag!(residues::Array{Residue, 1}, prec_id::Int32, charges::Array{Int32, 1}, isotopes::Array{Int32, 1}, y_start::Int32, b_start::Int32)
    peptide.transitions = getFragIons(residues, 
                                      prec_id, charges, isotopes, y_start, b_start
                                     )
    peptide.precursors = getPrecursors(residues,
                                       charges, isotopes, pep_id
                                       )
end

frag!(peptide::AbstractPeptide, mods_dict::Dict{String, Float32}, 
      prec_id::Int32, charges::Array{Int32, 1}, 
      isotopes::Array{Int32, 1}, y_start::Int32, b_start::Int32) = frag!(getResidues(getSequence(peptide), mods_dict), 
                                                                         prec_id, charges, 
                                                                         isotopes, y_start, b_start)

#getMZ(Frag(map(mod -> Residue(sequence[mod], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", sequence)), 'p', charge)),
# mutable struct Peptide <: AbstractPeptide
#     sequence::String,
#     unmodified_sequence::String
#     transitions::Array{Transition, 1}
#     precursors::Array{Precursors, 1}
#     mods::Array{String, 1}
#     function Peptide(sequence::String, charge::Int32, mods_dict::Dict{String, Float32})
#         new(
#             sequence,
#             Array{Frag, 1}(),
#             charge,
#             getMZ(Frag(map(mod -> Residue(sequence[mod], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", sequence)), 'p', charge)),
#             map(mod -> sequence[mod], findall(r"[A-Z]\[.*?\]", sequence))
#         )
#     end
# end

function getFrag(residues::Vector{Residue}, frag::NamedTuple)
    #get combinations of b, y, and p ions and charges that don't violate the filters
    #Loop through them. 
    if frag.ion_type == 'b'
        Frag(residues[1:frag.ind], frag.ion_type, frag.charge)
    elseif frag.ion_type == 'y'
        Frag(reverse(reverse(residues)[1:frag.ind]), frag.ion_type, frag.charge)
    elseif frag.ion_type == 'p'
        Frag(residues, frag.ion_type, frag.charge, frag.isotope)
    else
        throw(ErrorException(string("Ion type ", frag.ion_type," not recognized")))
    end
end




function frag!(peptide::Peptide, frags::Vector{NamedTuple}, mods_dict::Dict{String, Float32})
    peptide.fragments = map(frag -> getFrag(getResidues(peptide, mods_dict), frag), frags)
end

#function frag!(precursor::Peptide, frag_filter::String)
    #frag_filter "b(1-2;1-N);y(1-2;1-N);p(2-4;0-2)"
    #

#end
export getFrag
export Peptide
export frag!
export getResidues
export PROTON
export H2O
export NEUTRON
export getFragIons
export FragInd
print(Peptide("C[Carb]TIDEK[+8.014199]", Int32(2), default_mods))
print("hello")


# function Transition(residues::Vector{Residue}, frag_ions::NamedTuple)
#     if frag_ions.type=='b'
#         Transition()
#     elseif frag_ions.type=='y'
#     end
# end
# # struct Frag
# #     charge::Int32
# #     type::Char
# #     mz::Float32
# #     isotope::Int32
# #     function Frag(residues::Array{Residue, 1}, type::Char, charge::Int32, isotope::Int32)
# #         if type=='b'
# #             new(charge, type, (sum(residue->getMass(residue), residues) + PROTON*charge + isotope*NEUTRON)/charge, isotope)
# #         elseif type∈('y','p')
# #             new(charge, type, (sum(residue->getMass(residue), residues) + PROTON*charge + H2O + isotope*NEUTRON)/charge, isotope)
# #         #Could add functionality for a/x/c/z ions here
# #         end
# #     end
# end