##########
#Global Constants
##########
const H2O::Float32 = Float32(18.010565)
const PROTON::Float32 = Float32(1.0072764)
const NEUTRON::Float32 = Float32(1.00335)
const ion_mods = Dict{Char, Float32}(
        'y' => H2O,
        'b' => Float32(0),
        'p' => H2O,
    )
export PROTON
export H2O
export NEUTRON

const default_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464)
)

export default_mods
##########
#Representation of Amino Acid
##########
struct AA 
    aa::Char
    mass::Float32
    #Constructor for amino acid. Restric inputs
    #to valid amino acid symbols and assign correct mass
    function AA(aa::Char)
        AA_to_mass::Dict{Char, Float32} = Dict{Char, Float32}(
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
    mass::Float32 #Mass of modification
end

"""
    Mod(mod::String, mods_dict::Dict{String, Float32})
Given a string, first parse by the regular expression  
    Examples: "K", "K[+8.014199]" or "C[Carb]"
    
In the second case, "K[+8.014199]", "K[+8.014199]" is the modification name
    and 8.014199 is the modification mass. 

In the second case, "C[Carb]" is the modification name
    and "Carb" is a key to the dictionary "mods_dict".
    mods_dict["Carb"] returns the modification mass. 

If the Mod string can't be parsed, returns an error. 
"""
function Mod(mod::String, mods_dict::Dict{String, Float32})

    try
        #Single, unmodified amino acid
        if length(mod) == 1
            Mod(0.0)
        #Modified amino acid, so interpret modification
        else
            m = match(r"^[A-Z]\[(.*)\]$", mod)
            if startswith(m[1], "+")
                Mod(parse(Float32, m[1][2:end])) #"K[+8.014199]"
            else 
                Mod(mods_dict[m[1]]) #getAA("C[Carb]")
            end
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

function getResidues(sequence::String, mods_dict::Dict{String, Float32} = default_mods)
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
    prec_id::Int32 #Integer referencing the precursor
    ion_type::Char #'y' or 'b'. Could use others. 
    ind::UInt8 #for b3+1 ion, the "ind" is 3
    charge::UInt8
    iosotope::UInt8#diference in number of neutrons between the monoisotopic
end

struct Precursor <: Ion
    residues::Array{Residue, 1}
    mz::Float32
    charge::UInt32
    iosotope::UInt32 #diference in number of neutrons between the monoisotopic
    prec_id::Int32 
    pep_id::Int32
end

function Precursor()
    Precursor(Array{Residue, 1}(), Float32(0), UInt8(0), UInt8(0), Int32(0), Int32(0))
end
"""
    getIonMZ(residues::Array{Residue, 1}, charge::Int32; modifier::Float32 = PROTON*charge + H2O, isotope::Int32 = Int32(0))

Get the mz ratio of an ion given a list of the amino acid residues, the charge, 
the isotope state, and modifier. 

The modifier should depend on the kind of ion (y, b, precursor etc.). 
The default value `PROTON*charge + H2O` is correct for precursors and y ions. For b ions, `PROTON*charge` 
should be used instead. 

The isotope state is the difference in the number of isotopes from the monoisotopic
state. For a monoisotopic precurosr, `isotope` should equal the default value of zero. 
"""
function getIonMZ(residues::Array{Residue, 1}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))
    (sum(map(residue->getMass(residue), residues))+ #Sum residue masses
    (modifier + isotope*NEUTRON) #Add Proton and H2O
    )/charge #Divide mass by charge
end
export getIonMZ
#Constructor for precursor
"""
    Precursor(residues::Array{Residue, 1}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))

Constructor for the `Precursor` struct. Given a list of amino acid residues, a charge, and
an isotope state, makes a precursor object with the correct mz. 
(link to PrecursorMZ)
"""
function Precursor(residues::Array{Residue, 1}, charge::UInt8, isotope::UInt8 = UInt8(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))
    Precursor(
              residues,
              getIonMZ(residues, charge, isotope = isotope),
              charge, 
              isotope,
              prec_id,
              pep_id
            )
end

"""
    Precursor(sequence::String, mods_dict::Dict{String, Float32}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0)) 

Alternate constructor for the `Precursor` struct. Can accept a string representation of a peptide and a `mods_dict` 
and convert to residues `Array{Residue, 1}`. 
(link to getResidues())
"""
Precursor(sequence::String, mods_dict::Dict{String, Float32}, 
            charge::UInt8, isotope::UInt8 = UInt8(0), 
            prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0)
        ) = Precursor(getResidues(sequence, mods_dict), charge, isotope, prec_id, pep_id)

#Getter Functions
getFragMZ(ion::Ion) = ion.frag_mz
getPrecMZ(transition::Transition) = transition.prec_mz
getPrecCharge(transition::Transition) = transition.prec_id
getMZ(precursor::Precursor) = precursor.mz
getPepID(ion::Ion) = ion.pep_id
getIonType(transition::Transition) = transition.ion_type
getIonType(precursor::Precursor) = 'p'
getCharge(ion::Ion) = ion.charge
getIsotope(ion::Ion) = ion.isotope
getResidues(precursor::Precursor) = precursor.residues
getPrecID(ion::Ion) = ion.prec_id

function getIonModifier(ion_type::Char, charge::UInt8)::Float32
    ion_mods = Dict{Char, Float32}(
        'y' => H2O,
        'b' => Float32(0),
        'p' => H2O,
    )
    return ion_mods[ion_type] + PROTON*charge
end

"""
    Transition(residues::Array{Residue, 1}, prec_mz::Float32, ion_type::Char, charge::Int32, ind::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0))

Constructor for the `Transition` struct. Given a list of amino acids, the precursor mz, ion type (b, y, etc.), charge,
fragment position, isotope, and precursor id, returns a Transition with the correct fragment mz.  
(link to other usages)
"""
function Transition(residues::Array{Residue, 1}, prec_mz::Float32, ion_type::Char, charge::UInt8; ind::UInt8 = UInt8(length(residues)), isotope::UInt8 = UInt8(0), prec_id::Int32 = Int32(0))

    function getFragIonMZ(residues::Array{Residue, 1}, modifier::Float32, charge::UInt8, isotope::UInt8)              
        (sum(map(residue->getMass(residue), residues)) .+ (modifier + isotope*NEUTRON))/charge
    end
    print("hello")
    if ion_type == 'b'
        Transition(
                    getIonMZ(residues[1:ind], charge, modifier = getIonModifier(ion_type, charge), isotope = isotope), #frag_mz
                    prec_mz, prec_id,ion_type,ind,charge,isotope
                )       
    elseif ion_type == 'y'
        #Default modifier gives the correct fragment mz for a 'y' ion.
        #Only need to reverse the sequence to get y ions. 
        Transition(
                    getIonMZ(reverse(residues)[1:ind], charge, isotope = isotope),
                    prec_mz, prec_id,ion_type,ind,charge,isotope
        )     
    else
        throw(ErrorException(string("Ion type ", ion_type," not recognized")))
    end
end

"""
    Transition(sequence::String, mods_dict::Dict{String, Float32}, precursor::Precursor, ion_type::Char, charge::Int32, isotope::Int32 = Int32(0),  prec_id::Int32 = Int32(0))

Constructor for the `Transition` struct. Accepts a sequence string and mods_dict as objects. Can
get the residues array from these using getResidues. 
(link to getResidues)
"""
function Transition(sequence::String, prec_mz::Float32, ion_type::Char, charge::UInt8; mods_dict::Dict{String, Float32} = default_mods, isotope::UInt8 = UInt8(0), prec_id::Int32 = Int32(0))
    print("hellow")
    residues = getResidues(sequence, mods_dict)
    Transition(residues, 
                                                    prec_mz, 
                                                    ion_type, 
                                                    charge, 
                                                    ind = UInt8(length(residues)),
                                                    isotope = isotope, 
                                                    prec_id = prec_id)
end
#Change order of arguments in transition to reflect y2+1 type, ind, charge
function Transition(sequence::String, prec_mz::Float32, ion_type::Char, charge::UInt8, ind::UInt8; mods_dict::Dict{String, Float32} = default_mods, isotope::UInt8 = UInt8(0), prec_id::Int32 = Int32(0))
    Transition(getResidues(sequence, mods_dict), 
                                                    prec_mz, 
                                                    ion_type, 
                                                    charge, 
                                                    ind = ind,
                                                    isotope = isotope, 
                                                    prec_id = prec_id)
end
#unction getIonMZ(residues::Array{Residue, 1}, charge::Int32; modifier::Float32 = PROTON*charge + H2O, isotope::Int32 = Int32(0))
#    (sum(map(residue->getMass(residue), residues))+ #Sum residue masses
#    (modifier + isotope*NEUTRON) #Add Proton and H2O
#    )/charge #Divide mass by charge
#end
"""
    getFragIons(residues::Array{Residue, 1}, prec_mz::Float32, prec_id::Int32, modifier::Float32, ion_type::Char, start::Int32, charge::Int32, isotope::Int32)

Returns an `Array{Transition, 1}`. For the specific ion type and charge state (yn+2 for example) gets all fragment ions from start to the end. 
"""
function getAllIonMZ(residues::Array{Residue, 1}, charge::UInt8; start::Int = 3, modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))::Array{Float32, 1}
    (@view(cumsum((map(residue->getMass(residue), residues)))[start:end]).+ #Sum residue masses
    (modifier + isotope*NEUTRON) #Add Proton and H2O
    )./charge #Divide mass by charge
end
function getFragIons(residues::Array{Residue, 1}, prec_mz::Float32, prec_id::Int32, ion_type::Char, start::Int, charge::UInt8, isotope::UInt8)::Array{Transition, 1}
    #use cumsum
    #function getSlices(start::Int, residues::Array{Residue, 1})
    #    enumerate([residues[slice] for slice in [range(1, start+i) for i in range(0, length(residues)-start)]])
    #end

    #map(frag -> Transition(getIonMZ(frag[2], charge, modifier = getIonModifier(ion_type, charge), isotope = isotope), #frag_mz
    #                       prec_mz, prec_id, ion_type, charge,
    #                       Int32(frag[1])+start - 1, #ind. 3 for y3+n ion.
    #                    isotope), 
    #                    getSlices(start, residues)
    #    )
    #function getSlices(start::Int, residues::Array{Residue, 1})
    #    enumerate([residues[slice] for slice in [range(1, start+i) for i in range(0, length(residues)-start)]])
    #end
    
    map(frag -> Transition(frag[2]::Float32, #frag_mz
                            prec_mz::Float32, prec_id::Int32, ion_type::Char, charge::UInt8,
                            (frag[1]+start - 1)::Int, #ind. 3 for y3+n ion.
                            #ststart -1,
                        isotope::UInt8), 
                        enumerate(getAllIonMZ(residues, charge, start = start, modifier = ion_mods[ion_type] + PROTON*charge, isotope = isotope))
                        #getAllIonMZ(residues, charge, start = start, modifier = ion_mods[ion_type] + PROTON*charge, isotope = isotope)
        )
    #ion_mzs = Array{Float32, 1}(undef, length(residues)-start)
    #ion_mzs = getAllIonMZ(residues, charge, start = start, modifier = ion_mods[ion_type] + PROTON*charge, isotope = isotope);
    #ret = Array{Transition, 1}(undef, length(ion_mzs));
    #for frag in enumerate(ion_mzs)
    #    ret[frag[1]] = Transition(frag[2]::Float32, #frag_mz
    #                    prec_mz::Float32, prec_id::Int32, ion_type::Char, charge::Int32,
    ##                    (frag[1]+start - 1)::Int, #ind. 3 for y3+n ion.
    #                    #ststart -1,
    #                    isotope::Int32);
    #end
    #return ret
end


#Shorthands for getting precursors, b ions, and y ions
#Vector of precursors for given charge states and isotopes
#Really need to test this. 
getPrecursors(residues::Array{Residue, 1}, 
              charges::Vector{UInt8}, 
              isotopes::Vector{UInt8}, 
              pep_id::Int32) = Precursor.(residues, charges, isotopes, pep_id)
export getPrecursors

#Get all b and y ions 
function getFragIons(residues::Array{Residue, 1}, prec_mz::Float32, prec_id::Int32; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)
    vcat(
        getFragIons(residues, prec_mz,prec_id, 'b', b_start, charge, isotope),#,
        getFragIons(reverse(residues), prec_mz,prec_id, 'y', y_start, charge, isotope)
        #getFragIons(residues, prec_mz,prec_id, 'y', y_start, charge, isotope)
        )
end

getFragIons(precursor::Precursor; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3
            ) = getFragIons(getResidues(precursor),
                            getMZ(precursor),
                            getPrecID(precursor),
                            charge = charge, isotope = isotope, y_start = y_start, b_start = b_start
            )

function getFragIons(precursor::Precursor, charges::Array{UInt8, 1}; isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)
    collect(Iterators.flatten(map(charge -> getFragIons(precursor, charge = charge, isotope = isotope, y_start = y_start, b_start = b_start), charges)))           
end
#Can provide a list of named tuples to specify exactly which fragments to get
export Ion
export Transition
export Precursor
export getFragMZ
export getPrecMZ
export getMZ
export getPepID
export getIonType
export getCharge
export getIsotope
export getFragIons
export getBIons
export getYIons
export getPrecursors
export getFragIons
export getFragIonMZ

##########
#Peptide
##########
abstract type AbstractPrecursor end

mutable struct TargetPrecursor <: AbstractPrecursor
    sequence::String
    unmodified_sequence::String
    residues::Array{Residue, 1}
    transitions::Vector{Transition}
    precursor::Precursor
    #Some sort of dataframe with rows for scans and columns for transitions
    mods::Array{String, 1}
    isotope_label::String
    function TargetPeptide(sequence::String, mods_dict::Dict{String, Float32}; charge::Int32 = Int32(2), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))
        residues = getResidues(sequence, mods_dict)
        new(
            sequence, #sequence
            replace(sequence, r"(\[.*?\])"=>""), #unmodified_sequence
            residues, #residues
            Array{Transition, 1}(), #transitions
            Precursor(residues, charge, prec_id, pep_id), #mz for each isotope
            map(mod -> sequence[mod], findall(r"[A-Z]\[.*?\]", sequence)), #Mods
            isotope_label #isotope_label.
        )
    end
end

getResidues(precursor::AbstractPrecursor) = precursor.residues
getSequence(pep::AbstractPrecursor) = pep.sequence
getMZ(precursor::AbstractPrecursor) = getMZ(precursor.precursor)
getPrecID(precursor::AbstractPrecursor) = getPrecID(precursor.precursor)
getCharge(precursor::AbstractPrecursor) = getCharge(precursor.precursor)
getPepID(precursor::AbstractPrecursor) = getPepID(precursor.precursor)

function frag!(precursor::AbstractPrecursor, charges::Array{Int32, 1}; y_start::Int = 3, b_start::Int = 3, isotopes::Array{Int32, 1} = Array{Int32, 1}([0]))

    peptide.transitions = getFragIons.(getResidues(precursor), 
                                      getID(precursor), 
                                      getMZ(precursor), 
                                      charges, isotopes, y_start, b_start
                                     )

end

#Transition(residues::Array{Residue, 1}, prec_mz::Float32, ion_type::Char, charge::Int32, ind::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0))
function addFrag!(precursor::AbstractPrecursor, ion_type::Char, charge::Int32, isotope::Int32, ind::Int32)
    push!(precursor.transitions, Transition(getResidues(precursor), getMZ(precursor), ion_type, charge, ind, isotope, prec_id))
end

export Transition
export Precursor
export AbstractPeptide
export TargetPeptide

#print(Peptide("C[Carb]TIDEK[+8.014199]", Int32(2), default_mods))
#print("hello")
