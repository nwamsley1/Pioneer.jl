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
    name::String
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
            Mod(
                mod,
                0.0
            )
        elseif startswith(m[1], "+")
            Mod(
                mod,                    #"K[+8.014199]"
                parse(Float32, m[1][2:end]) #8.014199
                )
        else 
            Mod(
                mod,                #getAA("C[Carb]")
                mods_dict[m[1]]         #57.021464
                )
        end
    catch
        throw(ErrorException("$mod could not be parsed as given"))
    end 
end

#Optionally parse mods without a mods_dict
Mod(mod::String) = Mod(mod, Dict{String, Float32}())
Mod(name::Char, mass::Float32) = Mod(string(name), mass)
#Empty modification
Mod() = Mod("", 0.0)

#Getter Functions
getMass(mod::Mod) = mod.mass    
getName(mod::Mod) = mod.name
export Mod

##########
#Residue. Implementation of amino acid with custom mass modifications
##########
struct Residue
    aa::AA
    mod::Mod
    mass::Float32
end

#Residue(AA('A'))
function Residue(aa::AA)
    Residue(aa, Mod(), getMass(aa))
end

#Residue('A')
Residue(aa::Char) = Residue(AA(aa))

Residue(aa::AA, mod::Mod) = Residue(aa::AA, mod, getMass(mod)+getMass(aa))

function Residue(residue::String, mods_dict::Dict{String, Float32}) 
    if length(residue)>1
        Residue(AA(residue[1]), Mod(residue, mods_dict))
    else
        Residue(AA(residue[1]), Mod())
    end
end
function Residue(residue::String)
    if length(residue)>1
        Residue(AA(residue[1]), Mod(residue))
    else
        Residue(AA(residue[1]), Mod())
    end
end


Residue(residue::String, mod_mass::Float32) = Residue(AA(residue[1]), Mod(residue, mod_mass))

Residue(residue::Char, mod_mass::Float32) = Residue(AA(residue), Mod(residue, mod_mass))
#    """
#    Residue('K', )
#    
#    """
#    Residue(
#            AA(residue), 
#            Mod(join([residue,"[+", string(mod_mass),"]"]), mod_mass),
#            mod_mass
#            )
#end
#Getter methods
getMass(residue::Residue) = residue.mass
getMod(residue::Residue) = residue.mod
getAA(residue::Residue) = residue.aa

export Residue
export getMass
export getMod
export getAA
export default_mods

##########
#Frag
##########
struct Frag
    charge::Int32
    type::Char
    mz::Float32
    isotope::Int32
    function Frag(residues::Array{Residue, 1}, type::Char, charge::Int32, isotope::Int32)
        if type=='b'
            new(charge, type, (sum(residue->getMass(residue), residues) + PROTON*charge + isotope*NEUTRON)/charge, isotope)
        elseif type∈('y','p')
            new(charge, type, (sum(residue->getMass(residue), residues) + PROTON*charge + H2O + isotope*NEUTRON)/charge, isotope)
        #Could add functionality for a/x/c/z ions here
        end
    end
end

Frag(residues::Array{Residue, 1}, type::Char, charge::Int32) = Frag(residues, type, charge, Int32(0))

getCharge(frag::Frag) = frag.charge
getMZ(frag::Frag) = frag.mz
getType(frag::Frag) = frag.type
getIso(frag::Frag) = frag.isotope

export Frag
export getMZ
export getCharge
export getType
export getIso

##########
#Precursor
##########
mutable struct Peptide
    sequence::String
    fragments::Array{Frag, 1}
    charge::Int32
    mz::Float32
    mods::Array{String, 1}
    function Peptide(sequence::String, charge::Int32, mods_dict::Dict{String, Float32})
        new(
            sequence,
            Array{Frag, 1}(),
            charge,
            getMZ(Frag(map(mod -> Residue(sequence[mod], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", sequence)), 'p', charge)),
            map(mod -> sequence[mod], findall(r"[A-Z]\[.*?\]", sequence))
        )
    end

end

getSequence(peptide::Peptide) = peptide.sequence

function getResidues(peptide::Peptide, mods_dict::Dict{String, Float32})
    map(mod -> Residue(getSequence(peptide)[mod], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", getSequence(peptide)))
end

function getResidues(peptide::String, mods_dict::Dict{String, Float32})
    map(mod -> Residue(peptide[mod], mods_dict), findall(r"[A-Z]\[.*?\]|[A-Z]", peptide))
end

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

struct FragInd
    frag_mz::Float32
    prec_mz::Float32
    prec_id::Int32
    ion_type::Char
    charge::Int32
    iosotope::Int32
    ind::Int32
end

function getFragIons(residues::Vector{Residue}, prec_mz::Float32, prec_id::Int32, modifier::Float32, ion_type::Char, charge::Int32, isotope::Int32)
    function __getFragIons__(residues::Vector{Residue}, modifier::Float32, charge::Int32, isotope::Int32)
        (cumsum(map(residue->getMass(residue), residues) .+ (modifier + isotope*NEUTRON)))/charge
    end

    map(frag -> FragInd(frag[2], prec_mz, prec_id, ion_type, charge, isotope, Int32(frag[1])), enumerate(__getFragIons__(residues, modifier, charge, isotope)))
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
