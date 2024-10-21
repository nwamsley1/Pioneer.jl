
"""
parseInternalIon(frag_name::String)

Parse out the non-sequence information from the internal
ion annotation. See Examples 

### Input
-'frag_name::String'

### Output
-'::SubString{String}'

### Examples
julia> parseInternalIon("\"Int/PDTAENLPAmY-CO/7\"")
"-CO"
julia> parseInternalIon("\"Int/PDTAENLPAmY-CO^2+i/7\"")
"-CO^2+i"
"""
function parseInternalIon(
    frag_name::AbstractString
    ) #Parse an internal ion annotation to get rid of the sequence part "\"Int/PDTAENLPAmY-CO/7\"" => "-CO^2+i"
    #internal_ion_regex = r"[+\-^].*$"
    internal_ion_regex = r"[+\-^].*$"
    frag_name = split(frag_name, '/')[2] #"\"Int/PDTAENLPAmY-CO/7\"" => "PDTAENLPAmY-CO"
    frag_match = match(internal_ion_regex, frag_name) #"PDTAENLPAmY-CO" => RegexMatch("-H2O")
    if frag_match===nothing #No losses, isotopes, or >1 charge state 
        return ""
    else
        return frag_match.match*"\""
    end
end

"""
    getIonAnnotationSet(lib_path::String)

Extracts unique ion annotations from an .msp (MSPiano) library and returns a set 
of unique annotations. Special parsing required for internal ions. For example
"\"Int/PDTAENLPAmY-CO^2+i/7\"" poses a problem because there are many unique sequence
and internal annotation could have. We just want to extract the loss/isotope/charge information.
So we would parse "\"Int/PDTAENLPAmY-CO^2+i/7\"" => "-CO^2+i". 

### Input
-'lib_path::String' -- "path/to/lib.msp"

### Output
A Set of strings 
-'::Set{String}'


### Examples 
test_set = getIonAnnotationSet(test_lib)
julia> test_set
Set{String} with 4364 elements:
  "\"b9-H2O^2+3i\""
  "\"y17-CH3SOH-H2O^2\""
  "\"y16-CH3SOH^4\""
  "\"y17-H2O-NH3^4+i\""

"""
function getIonAnnotationSet(
            frag_names::AbstractVector{String}
            )
    unique_ion_types = Set{String}()
   
    for frag_name in frag_names
        if startswith(frag_name, "Int") #startswith("\"Int/PDTAENLPAmY-CO/7\"", "\"Int")
            push!(unique_ion_types, "Int/"*parseInternalIon(frag_name))
        else frag_name ∉ unique_ion_types #Regular and immonium ions can enter without pre-processing 
            push!(unique_ion_types, frag_name)
        end
    end
    return unique_ion_types
end

"""
    getIonAnnotationDict(ion_annotation_set::Set{String})

Given a set of strings representing ion annotations, map
each one to a unique 16-bit unsigned integer. Return a Julia
base Dict mapping the annotations to the integers. 

### Input
-'ion_annotation_set::String' -- ["\"b9-H2O^2+3i\"",...,"\"y10-CH3SOH-H2O+i\""]

### Output
A base julia Dict mappign ion annotations to unique 16-bit unsigned integers. 
-'::Dict{String, UInt32}'

### Examples 
test_set = getIonAnnotationSet(test_lib)
julia> getIonAnnotationDict(test_set)
Dict{String, UInt32} with 4364 entries:
  "\"b20-H2O^2+3i\""     => 0x000003c6
  "\"b13^3+i\""          => 0x00000217
  "\"y5-C2H5NOS\""       => 0x00000f85
  "\"b9-H2O^2+3i\""      => 0x00000639
  "\"b12-CH3SOH-H2O^3\"" => 0x000001a5
"""
function getIonAnnotationDict(
    ion_annotation_set::Set{String}
    )
    annotations = sort(collect(ion_annotation_set))
    annotations_to_id = Dict(
                            zip(
                                annotations, 
                                range(one(UInt32), UInt16(length(annotations)))
                            )
                        )

    return annotations_to_id
end

"""
    countSulfurLoss(molecular_formula::SubString{String},isloss::Bool)

Given a set of strings representing ion annotations, map
each one to a unique 16-bit unsigned integer. Return a Julia
base Dict mapping the annotations to the integers. 

### Input
-'molecular_formula::AbstractString' -- A molecualr formula like "-2CH3SOH", or "+CH3". Must have sign in front

### Output
8-bit integer representing the number of sulfurs added or subtracted from 
the fragment by the loss 
-'::Int8'

### Examples
julia> countSulfurLoss("+2CH3SOH")
2
julia> countSulfurLoss("-2CH3SOH")
-2
"""
function countSulfurLoss(
    molecular_formula::AbstractString, #"2CH3SOH"
    )
    #Formula multiplier
    multiplier = one(Int8) #Times the formula is repeated #"2CH3SOH" => 2 or "CH3SOH" => 1
    multiplier_string = match(r"[+-]([0-9]+)", molecular_formula) #Get times the formula is repeated if >1
    if multiplier_string!==nothing
        multiplier = parse(Int8, first(multiplier_string.captures))
    end
    sulfur_sign = startswith(molecular_formula, '-') ? -one(Int8) : one(Int8)
    #Sulfur Count 
    sulfur = match(r"S([0-9]+)*",molecular_formula) #Get number of sulfurs in the formula
    if sulfur===nothing #No sulfurs in the formula 
        return zero(Int8)
    elseif first(sulfur.captures) === nothing #Ony one sulfur in the formula 
        return multiplier*sulfur_sign
    else #Multiple sulfurs in the formula 
        return multiplier*sulfur_sign*parse(Int8, first(sulfur.captures))
    end
end

"""
    getNumeric(substring::SubString{String}, default::UInt8)

Finds the first numeric subset of a string. Parses to an 8-bit unsigned integer. 
If no numeric subset is round returns a default value. 

### Input
-'substring::AbstractString'
-'defualt::UInt8'

### Output
-'::UInt8'

### Examples
julia> getNumeric("x101x", zero(UInt8))
0x65
julia> getNumeric("x101x102", zero(UInt8))
0x65
julia> getNumeric("xxxxx", zero(UInt8))
0x00
julia> getNumeric("xxxxx", one(UInt8))
0x01
"""
function getNumeric(substring::AbstractString, default::R) where {R<:Real}
    numeric_capture = match(r"\d+", substring)
    if numeric_capture !== nothing
        return parse(R, numeric_capture.match)
    else
        return default
    end
end

"""
    getNumeric(immonium_table_path::String)

Given a tab delmited table, maps immonium ion annotations to 
the respective number of Sulfur's in the ion. See Examples below. 

### Input
-'immonium_table_path::String'

### Output
-'::Dict{String, String}'

### Examples
julia> getImmoniumToFormulaDict("data/immonium.txt")
Dict{String, Int8} with 64 entries:
  "IQA" => 0x0
  "IWF" => 0x0
  "ILC" => 0x0
  "IQB" => 0x0
"""
function getImmoniumToSulfurDict(
    immonium_table_path::String
    )
    immonium_to_sulfur_count = Dict{String, Int8}()
    open(immonium_table_path) do file
        for l in eachline(file)
            ion_name, formula = split(l, '\t')
            immonium_to_sulfur_count[ion_name] = countSulfurLoss(formula)
        end
    end

    return immonium_to_sulfur_count
end