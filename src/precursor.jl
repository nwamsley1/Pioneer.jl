##########
#Global Constants
##########
const H2O::Float32 = Float32(18.010565)
const PROTON::Float32 = Float32(1.0072764)
const NEUTRON::Float32 = Float32(1.00335)
const default_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464)
)
export PROTON, H2O, NEUTRON, default_mods
##########
#Representation of Amino Acid
##########
"""
AA

Type that represents an amino acid

### Fields

- aa::Char -- single character AA code 
- mass::Float32 -- mass of an amino acid

### Examples

- `AA(aa::Char)` -- Inner constructor

"""
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
"""
MzFeature

Type that represents a mass modification

### Fields

- mass:Float32 -- mass of the modification

### Examples

- `Mod(mass::Float32)` -- default constructor
- `Mod(mod::String, mods_dict::Dict{String, Float32} = Dict{String, Float32}()) -- Constructor 
that parses a string representation of an amino acid to get the mass of any modification of that amino acid`
- `Mod(name::Char) = Mod(0.0)` -- constructor for unmodified AA
- `Mod() = Mod(0.0)` -- constructor for no input

"""
struct Mod
    mass::Float32
end

"""
    Mod(mod::String, mods_dict::Dict{String, Float32})

Parses a string that is assumed to represent an amino acid and calculates its mass

### Input

- `mod::String`: -- String representation of a potentially modified amino acid. Constructor enforces formatting. 
- `mods_dict` -- Dictionary for named modifications that could appear in `mod`

### Output

Object of type "Mod" with field mass::Float32 that represents the mass of the modification

### Notes

Acceptable `mod` arguments match a regular expression defined in the method
There can be three types of modifications. 
- 1) There is no modification as in "K"
- 2) The mass is explicitly stated as in "K[+8.014199]"
- 3) The mass is named with a valid key for `mods_dict` as in "C[Carb]".
In this case mods_dict["Carb"] should return the appropriate mass

"""
function Mod(mod::String, mods_dict::Dict{String, Float32} = Dict{String, Float32}())

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

#Mod(string::String) = Mod(string, Dict{String, Float32}())
#Optionally parse mods without a mods_dict
Mod(name::Char) = Mod(0.0)
Mod() = Mod(0.0)
#Getter Functions
getMass(mod::Mod) = mod.mass    
export Mod

##########
#Residue. Implementation of amino acid with custom mass modifications
##########
"""
Residue

Type that represents a (potentially modified) amino acid within a peptide

### Fields

- mass:Float32 -- mass of the amino acid

### Examples

- `Residue(aa::AA)` -- default constructor
- `Residue(aa::Char) = Residue(AA(aa))`
- `Residue(aa::AA, mod::Mod) = Residue(getMass(mod)+getMass(aa))`
- `Residue(residue::String, mods_dict::Dict{String, Float32}) = Residue(AA(residue[1]), Mod(residue, mods_dict))`
- `Residue(residue::String) = Residue(residue, Dict{String, Float32}())`
- `Residue(residue::String, mod_mass::Float32) = Residue(getMass(AA(residue[1])) + mod_mass)`
- `Residue(residue::Char, mod_mass::Float32) = Residue(getMass(AA(residue)) + mod_mass)`

### Getter methods

- getMass(residue::Residue) = residue.mass

### See Also

- `getResidues(sequence::String, 
               mods_dict::Dict{String, Float32} = default_mods)` -- Gets a vector of residues given a string
   representation of an amino acid sequence
"""
struct Residue
    mass::Float32
end

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
#Ion, Transition, Precursor
##########
"""
MzFeature

Type that represents an m/z ratio with a ppm tolerance

### Fields

- mono::Float32 -- MZ with upper and lower bounds given a ppm tolerance
- low::Float32 -- Identifier of the precursor ion (parent ion of the fragment/transition)
- high::Float32 -- Type of transition. For example 'b' for b ion or 'y' for y ion

### Examples

- `MzFeature(mono::Float32; ppm::Float32 = Float32(20))` -- default internal constructor
- `MzFeature() = MzFeature(Float32(0.0))` -- constructor for a default/placeholder MzFeature

### GetterMethods

- getMZ(mz_feature::MzFeature) = mz_feature.mono
- getLow(mz_feature::MzFeature) = mz_feature.low
- getHigh(mz_feature::MzFeature) = mz_feature.high
"""
struct MzFeature
    mono::Float32
    low::Float32
    high::Float32
    function MzFeature(mono::Float32; ppm::Float32 = Float32(20))
        new(mono, Float32(mono*(1 - ppm/1e6)), Float32(mono*(1 + ppm/1e6)))
    end
end

MzFeature() = MzFeature(Float32(0.0))

getMZ(mz_feature::MzFeature) = mz_feature.mono
getLow(mz_feature::MzFeature) = mz_feature.low
getHigh(mz_feature::MzFeature) = mz_feature.high

"""
    Ion

Abstract type that represents an ion 

Types that inherit from `Ion` should implement the following. 

- getMZFeature(ion::Ion) = ion.mz
- getMZ(ion::Ion) = getMZ(getMZFeature(ion))
- getLow(ion::Ion) = getLow(getMZFeature(ion))
- getHigh(ion::Ion) = getHigh(getMZFeature(ion))
- getPrecID(ion::Ion) = ion.prec_id
- getCharge(ion::Ion) = ion.charge
- getIsotope(ion::Ion) = ion.isotope

"""
abstract type Ion end

getMZFeature(ion::Ion) = ion.mz
getMZ(ion::Ion) = getMZ(getMZFeature(ion))
getLow(ion::Ion) = getLow(getMZFeature(ion))
getHigh(ion::Ion) = getHigh(getMZFeature(ion))
getPrecID(ion::Ion) = ion.prec_id
getCharge(ion::Ion) = ion.charge
getIsotope(ion::Ion) = ion.isotope

"""
Transition <: Ion

Type that represents transition (fragment ion of a peptide)

### Fields

- mz::MzFeature -- MZ with upper and lower bounds given a ppm tolerance
- prec_id::UInt32 -- Identifier of the precursor ion (parent ion of the fragment/transition)
- ion_type::Char  -- Type of transition. For example 'b' for b ion or 'y' for y ion
- ind::UInt8 -- Position of fragment ion with reference to parent ion. (A b5+2 ion should have an ind equal to 5)
- charge::UInt8 -- Charge of the fragment ion
- isotope::UInt8 -- Difference in number of neutrons from the monoisotopic fragment 

### Examples

- `Transition(frag_mz::Float32, prec_id::UInt32, ion_type::Char, ind::UInt8, 
              charge::UInt8, 
              isotope::UInt8; 
              ppm = Float32(20))` -- default internal constructor
- `Transition(residues::Vector{Residue}; ion_type::Char = 'y', charge::UInt8 = UInt8(1), 
              ind::UInt8 = UInt8(length(residues)), 
              isotope::UInt8 = UInt8(0), 
              prec_id::UInt32 = UInt32(0))` -- constructor that calculates the appropriate mz

### GetterMethods

- getIonType(transition::Transition) = transition.ion_type
"""
struct Transition <: Ion
    mz::MzFeature
    prec_id::UInt32 #Integer referencing the precursor
    ion_type::Char #'y' or 'b'. Could use others. 
    ind::UInt8 #for b3+1 ion, the "ind" is 3
    charge::UInt8
    iosotope::UInt8 #diference in number of neutrons between the monoisotopic
    function Transition(frag_mz::Float32, prec_id::UInt32, ion_type::Char, ind::UInt8, charge::UInt8, isotope::UInt8; ppm = Float32(20))
        new(MzFeature(frag_mz, ppm = ppm), 
                       prec_id, ion_type, ind, charge, isotope
                       )
                       
    end
    #Internal constructor that makes lower 
    #and upper mz bounds given an optionall ppm tolerance. 
end

"""
    Transition(residues::Vector{Residue}; ion_type::Char = 'y', charge::UInt8 = UInt8(1), ind::UInt8 = UInt8(length(residues)), isotope::UInt8 = UInt8(0), prec_id::UInt32 = UInt32(0))
    
    Constructor for the `Transition` struct. Given a list of amino acid residues, ion type (b or y), charge state, and isotopic state makes a transition with the correct mz. 
"""
function Transition(residues::Vector{Residue}; ion_type::Char = 'y', charge::UInt8 = UInt8(1), ind::UInt8 = UInt8(length(residues)), isotope::UInt8 = UInt8(0), prec_id::UInt32 = UInt32(0))
    #function Transition(frag_mz::Float32, prec_id::UInt32, ion_type::Char, ind::UInt8, charge::UInt8, isotope::UInt8; ppm = Float32(20))
    #function getFragIonMZ(residues::Vector{Residue}, modifier::Float32, charge::UInt8, isotope::UInt8)              
    #    (sum(map(residue->getMass(residue), residues)) .+ (modifier + isotope*NEUTRON))/charge
    #end
    if ion_type == 'b'
        Transition(
                    getIonMZ(residues[1:ind], charge, modifier =  getBIonModifier(charge), isotope = isotope), #frag_mz
                    prec_id,ion_type,ind,charge,isotope
                )       
    elseif ion_type == 'y'
        #Default modifier gives the correct fragment mz for a 'y' ion.
        #Only need to reverse the sequence to get y ions. 
        Transition(
                    getIonMZ(reverse(residues)[1:ind], charge, modifier =  getYIonModifier(charge), isotope = isotope),
                    prec_id,ion_type,ind,charge,isotope
        )     
    else
        throw(ErrorException(string("Ion type ", ion_type," not recognized")))
    end
end

getIonType(transition::Transition) = transition.ion_type

"""
Precursor <: Ion

Type that represents a precursor (A peptide parent ion)

### Fields

- residues::Vector{Residue} -- List of amino acid residues of the precursor in their appropriate order
- mz::MzFeature -- MZ with upper and lower bounds given a ppm tolerance
- prec_id::UInt32 -- Identifier of the precursor ion (parent ion of the fragment/transition)
- charge::UInt8 -- Charge of the fragment ion
- isotope::UInt8 -- Difference in number of neutrons from the monoisotopic fragment 
- pep_id::UInt32 -- Identifier of the peptide from which precursor is derived

### Examples

- `Precursor(residues::Vector{Residue}, mz::Float32, charge::UInt8, 
             isotope::UInt8, 
             pep_id::UInt32,
             prec_id::UInt32; 
             ppm = Float32(20))` -- default internal constructor
- `Precursor()` -- constructor for null/empty precursor
- `Precursor(residues::Vector{Residue}, charge::UInt8, 
             isotope::UInt8 = UInt8(0), 
             pep_id::UInt32 = UInt32(0),
             prec_id::UInt32 = UInt32(0)` -- Constructor that calculates mz without having to supply it
- `Precursor(sequence::String; mods_dict::Dict{String, Float32} = Dict{String, Float32}(), charge::UInt8 = UInt8(2), 
             isotope::UInt8 = UInt8(0), 
             pep_id::UInt32 = UInt32(0),
             prec_id::UInt32 = UInt32(0))` -- Constructor that accepts a string representation of a peptide

### GetterMethods

- getResidues(precursor::Precursor) = precursor.residues
"""
struct Precursor <: Ion
    residues::Vector{Residue}
    mz::MzFeature
    charge::UInt8
    iosotope::UInt8
    pep_id::UInt32
    prec_id::UInt32
    """
        Precursor(residues::Array{Residue, 1}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))

    Inner constructor for the `Precursor` struct. Given a list of amino acid residues, an mz, charge, and
    an isotope state, makes a precursor object with the correct mz.
    """
    function Precursor(residues::Vector{Residue}, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32, prec_id::UInt32; ppm = Float32(20))
        new(residues, MzFeature(mz, ppm = ppm), 
                       charge, isotope, pep_id, prec_id)
            
    end
    #Internal constructor that makes lower 
    #and upper mz bounds given an optionall ppm tolerance. 
end

"""
    Precursor(residues::Array{Residue, 1}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))

Constructor for the `Precursor` struct. Given a list of amino acid residues, a charge, and
an isotope state, makes a precursor object with the correct mz. 
(link to Precursor)
"""
function Precursor(residues::Vector{Residue}, charge::UInt8, isotope::UInt8 = UInt8(0), pep_id::UInt32 = UInt32(0), prec_id::UInt32 = UInt32(0))
    Precursor(
              residues,
              getIonMZ(residues, charge, isotope = isotope),
              charge, 
              isotope,
              pep_id,
              prec_id
            )
end

"""
    Precursor(sequence::String, mods_dict::Dict{String, Float32}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0)) 

Alternate constructor for the `Precursor` struct. Can accept a string representation of a peptide and a `mods_dict` 
and convert to residues `Array{Residue, 1}`. 
(link to Precusor)
"""
Precursor(sequence::String; mods_dict::Dict{String, Float32} = Dict{String, Float32}(), 
            charge::UInt8 = UInt8(2), isotope::UInt8 = UInt8(0), 
            pep_id::UInt32 = UInt32(0), prec_id::UInt32 = UInt32(0)
        ) = Precursor(getResidues(sequence, mods_dict), charge, isotope, pep_id, prec_id)

"""
    Precursor()
    Constructor for an "empty" or "default" precursor
"""
Precursor() = Precursor(Vector{Residue}(), MzFeature(), UInt8(0), UInt8(0), UInt32(0), UInt32(0))

getResidues(precursor::Precursor) = precursor.residues

import Base.length
length(precursor::Precursor) = length(getResidues(precursor))

"""
    getBIonModifier(charge::UInt8)
Mass modification that needs to be added to b-ions
"""
function getBIonModifier(charge::UInt8)
    PROTON*charge
end

"""
    getYIonModifier(charge::UInt8)
Mass modification that needs to be added to y-ions
"""
function getYIonModifier(charge::UInt8)
    PROTON*charge + H2O
end

"""
    getIonMZ(residues::Vector{Residue}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))::Float32

Get the mz ratio of an ion

### Input

- `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion
- `charge::UInt8` -- Charge of the ion
- `modifier::Float32=PROTON*charge + H2O` -- Added to the mass of the ion
- `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion. 

### Output

A Float32 representing the mass-to-charge ratio (m/z) of an ion

### Notes

The `modifier` argument ought to depend on the kind of ion. For B ions PROTON*charge is appropriate,
but for 'y' or precursor ions, PROTON*charge + H2O would be appropriate.

### Algorithm 

Sum the amino acid residue masses, add `modifier` + isotope*NEUTRON and then divide the total by the charge. 

### Examples 

#Gets the b6+1 ion MZ
```julia-repl
julia> getIonMZ(getResidues("PEPTIDE")[1:6], UInt8(1), modifier = PROTON)
653.314f0
```
#Gets the y6+1 ion MZ
```julia-repl
julia> getIonMZ(reverse(getResidues("PEPTIDE"))[1:6], UInt8(1))
703.3144f0
```

"""
function getIonMZ(residues::Vector{Residue}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))
    #modifier should be "PROTON*charge + H2O" for 'p' and 'y' ions
    #and "PROTON*charge" for 'b' ions. 
    (sum(map(residue->getMass(residue), residues))+ #Sum residue masses
    (modifier + isotope*NEUTRON) #Add Proton and H2O
    )/charge #Divide mass by charge
end

"""
    getIonMZ(residues::Vector{Residue}, ion_type::Char, charge::UInt8; isotope::UInt8 = UInt8(0))

Alternate getIonMZ method that chooses the correct mass modifier for 'b', 'y', and 'p' ions respectively. 

### Input

    - `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion
    - `ion_type::Char` -- Type of fragment ion. Currently supports, 'b', 'y', and 'p'. 
    - `charge::UInt8` -- Charge of the ion
    - `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion. 

### Notes 
    See main method
        getIonMZ(residues::Vector{Residue}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))
    for more details
"""
function getIonMZ(residues::Vector{Residue}, ion_type::Char, charge::UInt8; isotope::UInt8 = UInt8(0))
    #modifier should be "PROTON*charge + H2O" for 'p' and 'y' ions
    #and "PROTON*charge" for 'b' ions. 
    if ion_type == 'b'
        (sum(map(residue->getMass(residue), residues))+ #Sum residue masses
        (getBIonModifier(charge::UInt8) + isotope*NEUTRON) #Add Proton and H2O
        )/charge #Divide mass by charge
    elseif ion_type ∈ ('y','p')
        (sum(map(residue->getMass(residue), residues))+ #Sum residue masses
        (getYIonModifier(charge::UInt8) + isotope*NEUTRON) #Add Proton and H2O
        )/charge #Divide mass by charge
    else
        throw(ErrorException(string("Ion type ", ion_type," not recognized")))
    end
end 

export getIonMZ
#Constructor for precursor


"""
    getPrecursors(residues::Vector{Residue}; charges::Vector{UInt8} = UInt8[1,2],isotopes::Vector{UInt8}=UInt8[0],pep_id::UInt32=UInt32(0)) 
    
    Alternate constructor for the `Precursor` struct that Can accept a string representation of a peptide
    ### Input

    ### Output 

    ### Notes


    (link to getResidues())
 """
function getPrecursors(residues::Vector{Residue}; charges::Vector{UInt8} = UInt8[1,2],isotopes::Vector{UInt8}=UInt8[0],pep_id::UInt32=UInt32(0))
    [Precursor(residues, charge, isotope, pep_id) for charge in charges for isotope in isotopes]
end

export getPrecursors



"""
    getIonSeries(residues::Vector{Residue}, charge::UInt8; start::Int = 3, modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))

Gets the m/z's for an ion series as a Vector{Float32}. 

### Input

- `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion
- `charge::UInt8` -- Charge of the fragment ions
- `start::Int=3`  -- Index of first ion the the series to compute.
- `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion. 

### Output

A Vector{Float32} wich each m/z in the ion series

### Notes

- The `modifier` argument ought to depend on the kind of ion. For a 'b' ion series PROTON*charge is appropriate,
but for a 'y' ion series, PROTON*charge + H2O would be appropriate. 

- Will not allow the index of the ion to be equal to or less than the charge. For example,
b2+2 ions could only be calculated in error and are therefore excluded even if `start` 
is set to 2. 

- If `start` exceeds length(residues)-1, then only the N-1 ion is calculated, that is,
the highest mass ion in the series. 

### Examples 

#Gets the y3+2 through y6+2 ions
```julia-repl
julia> getIonSeries(reverse(getResidues("PEPTIDE")), UInt8(2), start = 3)
4-element view(::Vector{Float32}, 3:6) with eltype Float32:
 188.58934
 239.11317
 287.63956
 352.16086
```
#Gets the b4+2 through b6+2 ions
```julia-repl
julia> getIonSeries(getResidues("PEPTIDE"), UInt8(2), start = 4, modifier = PROTON*UInt8(2))
3-element view(::Vector{Float32}, 4:6) with eltype Float32:
 213.10518
 269.6472
 327.16064
```

"""
function getIonSeries(residues::Vector{Residue}, charge::UInt8; start::Int = 3, modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))
    N = length(residues)
    cumulative_mass = zeros(Float32, N)
    cumulative_mass[1] = getMass(residues[1])
    for i in 2:N
        cumulative_mass[i] = cumulative_mass[i-1] + getMass(residues[i])
    end
    @views begin
        #min(max(charge+1, start), N-1):N-1 has three purposes
        #1) Exclude the last residue (N-1 instead of N) because that would just be the precursor ion
        #2) The minimum fragment length must be greater than the fragment charge
        #3) The fragment length can't exceed N-1 (the length of the peptide minus 1)
        @. ((cumulative_mass + (modifier + isotope*NEUTRON)) / charge)[min(max(charge+1, start), N-1):N-1]
    end 
end

"""
    getFragIons(residues::Vector{Residue}; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)

Uses calls to `getIonSeries` to concatenate both the b and y ion series together in a single Vector{Float32}

### Input

- `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion
- `charge::UInt8` -- Charge of the fragment ions
- `b_start::Int=3`  -- Index of first ion the the b-ion series to compute. 
- `y_start::Int=3`  -- Index of first ion the the y-ion series to compute. 
- `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion. 

### Output

A Vector{Float32} wich each m/z in the ion series

### Notes

- The `modifier` argument ought to depend on the kind of ion. For a 'b' ion series PROTON*charge is appropriate,
but for a 'y' ion series, PROTON*charge + H2O would be appropriate. 

- Will not allow the index of the ion to be equal to or less than the charge. For example,
b2+2 ions could only be calculated in error and are therefore excluded even if `start` 
is set to 2. 

- If `start` exceeds length(residues)-1, then only the N-1 ion is calculated, that is,
the highest mass ion in the series. 

### Examples 

#Gets the b3+1-b6+1 and y3+1-y6+1 ions
```julia-repl
julia> getFragIons(reverse(getResidues("PEPTIDE")), b_start = 3, y_start = 3)
8-element Vector{Float32}:
 358.16083
 459.2085
 556.2612
 685.30383
 342.16595
 443.21365
 556.29767
 671.3246
```

### See Also

Alternate convience methods
- `getFragIons(residues::Vector{Residue}; charge::UInt8 = UInt8(1),
              isotope::UInt8 = UInt8(0), 
              y_start::Int = 3, 
              b_start::Int = 3)` - Default method

- `getFragIons(precursor::Precursor; charge::UInt8 = UInt8(1), 
              isotope::UInt8 = UInt8(0), 
              y_start::Int = 3, 
              b_start::Int = 3)` - Can supply a `Precursor` rather than a  `Vector{Residues}` input

- `getFragIons(precursor::Precursor,charges::Vector{UInt8}, 
              isotopes::Vector{UInt8}; 
              y_start::Int = 3, 
              b_start::Int = 3)` - Gets b and y ion seriers for multiple charge and isotopic states
"""
function getFragIons(residues::Vector{Residue}; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)
    vcat(getIonSeries(residues, charge, start = b_start, modifier =  getBIonModifier(charge), isotope = isotope),
    getIonSeries(reverse(residues), charge, start = y_start, modifier =  getYIonModifier(charge), isotope = isotope))
end

getFragIons(precursor::Precursor; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3
            ) = getFragIons(getResidues(precursor),
                            charge = charge, isotope = isotope, y_start = y_start, b_start = b_start
            )

function getFragIons(precursor::Precursor,charges::Vector{UInt8}, isotopes::Vector{UInt8}; y_start::Int = 3, b_start::Int = 3)
    vec(hcat([getFragIons(getResidues(precursor), charge = charge, isotope = isotope, y_start = y_start, b_start = b_start) for charge in charges for isotope in isotopes]...))
end 

export getAllIonMz
function getTransitionSeries(residues::Vector{Residue}, charge::UInt8 = UInt8(1); ion_type::Char = 'y', prec_id::UInt32 = UInt32(0), isotope::UInt8 = UInt8(0), start::Int = 3, modifier::Float32 = H2O + PROTON*charge)
        map(transition -> Transition(transition[2]::Float32
                            , prec_id, ion_type, charge,
                            UInt8(transition[1]+start - 1), #ind. 3 for y3+n ion.
                            #ststart -1,
                        isotope), 
                        enumerate(
                                getIonSeries(residues, 
                                              charge, isotope = isotope,
                                              start = start, 
                                              modifier = modifier
                                              )
                                )
        )
end

function getTransitions(precursor::Precursor; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)
    vcat(getTransitionSeries(getResidues(precursor), charge, ion_type = 'b', prec_id = getPrecID(precursor), isotope = isotope, start = b_start, modifier = PROTON*charge),
    getTransitionSeries(reverse(getResidues(precursor)), charge, ion_type = 'y', prec_id = getPrecID(precursor), isotope = isotope, start = y_start, modifier = H2O + PROTON*charge))
end

function getTransitions(precursor::Precursor, charges::Vector{UInt8}, isotopes::Vector{UInt8}; y_start::Int = 3, b_start::Int = 3)
    vec(hcat([getTransitions(precursor, charge = charge, isotope = isotope, y_start = y_start, b_start = b_start) for charge in charges for isotope in isotopes]...))
end

export getTransitions


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
export getResidues
export MzFeature
export getHigh, getLow, MzFeature
##########
#Peptide
##########
#=
abstract type AbstractPrecursor end

mutable struct TargetPrecursor <: AbstractPrecursor
    sequence::String
    unmodified_sequence::String
    residues::Vector{Residue}
    transitions::Vector{Transition}
    precursor::Precursor
    #Some sort of dataframe with rows for scans and columns for transitions
    mods::Vector{String}
    isotope_label::String
    function TargetPeptide(sequence::String, mods_dict::Dict{String, Float32}; charge::Int32 = Int32(2), pep_id::Int32 = Int32(0))
        residues = getResidues(sequence, mods_dict)
        new(
            sequence, #sequence
            replace(sequence, r"(\[.*?\])"=>""), #unmodified_sequence
            residues, #residues
            Vector{Transition}(), #transitions
            Precursor(residues, charge, pep_id), #mz for each isotope
            map(mod -> sequence[mod], findall(r"[A-Z]\[.*?\]", sequence)), #Mods
            isotope_label #isotope_label.
        )
    end
end

getResidues(precursor::AbstractPrecursor) = precursor.residues
getSequence(pep::AbstractPrecursor) = pep.sequence
getMZ(precursor::AbstractPrecursor) = getMZ(precursor.precursor)
getCharge(precursor::AbstractPrecursor) = getCharge(precursor.precursor)
getPepID(precursor::AbstractPrecursor) = getPepID(precursor.precursor)

function frag!(precursor::AbstractPrecursor, charges::Vector{UInt8}; y_start::Int = 3, b_start::Int = 3, isotopes::Vector{Int32} = Vector{Int32}([0]))

    peptide.transitions = getFragIons.(getResidues(precursor), 
                                      getID(precursor), 
                                      getMZ(precursor), 
                                      charges, isotopes, y_start, b_start
                                     )

end
export Transition
export Precursor
export AbstractPeptide
export TargetPeptide
=#

