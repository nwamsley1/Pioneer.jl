##########
#Global Constants
##########
const H2O::Float64 = Float64(18.010565)
const PROTON::Float64 = Float64(1.0072764)
const NEUTRON::Float64 = Float64(1.00335)

const default_mods::Dict{String, Float64} = Dict{String, Float64}(
    "Carb" => Float64(57.021464))

const AA_to_mass::Dict{Char, Float64} = Dict{Char, Float64}(
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
    

function getMass(residue::Char)::Float64
    if haskey(AA_to_mass, residue)
        AA_to_mass[residue]
    else
        0.0
    end
end

function getMass(sequence::String)::Float64
    mass = zero(Float64)
    for aa in sequence
        mass += getMass(aa)
    end
    return mass
end


"""
    parseMods(mods_string::AbstractString)::Base.RegexMatchIterator

Given the meta data line (starting in "Comment") for a preursor, parse out the fields. 

### Input
-'mods_string::AbstractString' -- String with modifications for hte precursor "1(7,M,Oxidation)(13,K,AnExampleMod)"

### Output
-'::Base.RegexMatchIterator}' iterator of RegexMatch for each modification parsed out of the string

### Examples
julia> collect(
            parseMods("1(7,M,Oxidation)(13,K,AnExampleMod)")
        )
2-element Vector{RegexMatch}:
 RegexMatch("7,M,Oxidation")
 RegexMatch("13,K,AnExampleMod")
"""
function parseMods(mods_string::AbstractString)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, mods_string)
end

function parseMods(mods_string::Missing)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, "")
end



"""
    getModName(mod_string::AbstractString)::String

Given the a modification string "7,M,Oxidation", parse out the name

### Input
-'mod_string::AbstractString' -- Modification string

### Output
-'::String' modifictaion name

### Examples
julia> getModName("7,M,Oxidation")
"Oxidation"
"""
function getModName(mod_string::AbstractString)::String
    match(r"[^,]+(?=$)", mod_string).match
end

"""
    getModIndex(mod_string::AbstractString)::UInt8

Given the a modification string "7,M,Oxidation", parse out the index

### Input
-'mod_string::AbstractString' -- Modification string

### Output
-'::UInt8' modifictaion index

### Examples
julia> getModName("7,M,Oxidation")
0x07
"""
function getModIndex(mod_string::AbstractString)::UInt8
    parse(UInt8, match(r"^[0-9]+(?=,)", mod_string).match)
end

function getModMass(mods::String, mod_to_mass::Dict{String, Float64})
    mass = zero(Float64)
    for mod in parseMods(mods)
        mass += mod_to_mass[getModName(mod.match)]
    end
    return mass
end

function getModMass(mod::Missing)
    return zero(Float64)
end

function getMZ(sequence::String, mods::String, charge::UInt8, mod_to_mass::Dict{String, Float64})
    mass = getMass(sequence) + getModMass(mods, mod_to_mass)
    return (mass + PROTON*charge + H2O)/charge
end

function getMZ(sequence::String, mods::Missing, charge::UInt8, mod_to_mass::Dict{String, Float64})
    mass = getMass(sequence)
    return (mass + PROTON*charge + H2O)/charge
end

function getMZs(sequences::AbstractVector{String}, mods::AbstractVector{Union{Missing, String}}, charges::AbstractVector{UInt8}, mod_to_mass::Dict{String, Float64})
    mz = Vector{Float32}(undef, length(sequences))
    for i in range(1, length(sequences))
        mz[i] = getMZ(sequences[i], mods[i], charges[i], mod_to_mass)
    end
    return mz
end