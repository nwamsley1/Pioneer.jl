# src/fragments/fragment_annotation.jl

"""
    get_ion_annotation_set(frag_names::AbstractVector{String})::Set{String}

Extract unique ion annotations from fragment names, with special handling for internal ions.

For regular ions, adds the annotation directly to the set. For internal ions, processes
them to extract only the loss/isotope/charge information. For example:
"Int/PDTAENLPAmY-CO^2+i/7" => "Int/-CO^2+i"

Parameters:
- frag_names::AbstractVector{String}: Vector of fragment annotations/names

Returns:
- Set{String}: Set of unique ion annotations, with internal ions processed to their
  generic form.

Examples:
```julia
# Regular annotations remain as-is
annotations = ["b9-H2O^2+3i", "y17-H2O-NH3^4+i"]
get_ion_annotation_set(annotations)  # Keeps these exactly as input

# Internal ions get processed
annotations = ["Int/PDTAENLPAmY-CO/7", "Int/PEPTIDE-H2O^2+i/4"]
get_ion_annotation_set(annotations)
# Returns Set(["Int/-CO", "Int/-H2O^2+i"])  # Sequence removed
```

Notes:
- Internal ion processing preserves:
  - Neutral losses (e.g., -H2O, -NH3)
  - Charge states (e.g., ^2+)
  - Isotope states (e.g., +i)
- Regular and immonium ions are preserved exactly as input
"""
function get_ion_annotation_set(frag_names::AbstractVector{String})::Set{String}
    unique_ion_types = Set{String}()
   
    for frag_name in frag_names
        if startswith(frag_name, "\"Int")  # Internal ion
            # Extract modification info and create generic form
            push!(unique_ion_types, "Int/" * parse_internal_ion(frag_name))
        else  # Regular or immonium ion
            push!(unique_ion_types, frag_name)
        end
    end
    
    return unique_ion_types
end

"""
    parse_internal_ion(frag_name::AbstractString)::String

Extract modification information from internal ion annotation.
Processes internal ion strings to get only the loss/isotope/charge information.

Example:
"\"Int/PDTAENLPAmY-CO/7\"" => "-CO"
"\"Int/PDTAENLPAmY-CO^2+i/7\"" => "-CO^2+i"

Notes:
- Extracts neutral losses (e.g., -CO, -H2O)
- Extracts charge state if present (e.g., ^2+)
- Extracts isotope state if present (e.g., +i)
- Removes sequence and position information
"""
function parse_internal_ion(frag_name::AbstractString)::String
    # Match modifications/charge/isotope after the sequence
    internal_ion_regex = r"[+\-^].*$"
    
    # Split to get the middle part containing sequence and mods
    parts = split(frag_name, '/')
    sequence_and_mods = parts[2]  # "PDTAENLPAmY-CO^2+i"
    
    # Find modifications/charge/isotope
    mod_match = match(internal_ion_regex, sequence_and_mods)
    
    if isnothing(mod_match)
        return ""  # No modifications or special states
    else
        return mod_match.match * "\""  # Include closing quote
    end
end

"""
    create_ion_annotation_index(ion_annotation_set::Set{String})::Dict{String,UInt16}

Create a mapping from ion annotations to unique indices.
Useful for efficient storage and lookup of ion types.

Parameters:
- ion_annotation_set::Set{String}: Set of unique ion annotations

Returns:
- Dict{String,UInt16}: Mapping from annotations to indices

Example:
```julia
annotations = Set(["b9^2", "y8-H2O", "Int/-CO^2+i"])
index = create_ion_annotation_index(annotations)
# Returns Dict("b9^2" => 0x0001, "y8-H2O" => 0x0002, ...)
```
"""
function create_ion_annotation_index(ion_annotation_set::Set{String})::Dict{String,UInt16}
    # Sort for deterministic indexing
    sorted_annotations = sort(collect(ion_annotation_set))
    
    # Create mapping from annotation to index
    return Dict(
        annotation => UInt16(i)
        for (i, annotation) in enumerate(sorted_annotations)
    )
end

"""
    parse_fragment_annotation(annotation_pieces::Vector{RegexMatch},internal::Bool,immonium::Bool)

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
function parse_fragment_annotation(frag_regex_iterator::UniSpecFragRegexIterator,
                        immonium_to_sulfur_count::Dict{String, Int8},
                        internal::Bool,
                        immonium::Bool)
    annotation_pieces = getAnnotationPieces(frag_regex_iterator)
    #Initialize Data 
    base_type = '_'
    frag_index = one(UInt8)
    charge = one(UInt8)
    isotope = zero(UInt8)
    sulfur_diff = zero(Int8)
    is_gain_loss = false
    #= Example 
    julia>annotation_pices 
    4-element Vector{RegexMatch}:
    RegexMatch("b11")
    RegexMatch("-2H2O")
    RegexMatch("^3")
    RegexMatch("+2i")
    =#
    for annotation_piece in annotation_pieces
        if occursin('i', annotation_piece.match) #"+2i" isotope_stage
            match = annotation_piece.match
            if occursin('/', annotation_piece.match) #"/2i" isotope stage
                match = first(split(annotation_piece.match, '/'))
            end
            if startswith(match, '+') #"+2i" isotope stage
                #Isotope stage can't be higher than 9
                isotope = getNumeric(match, one(UInt8))
                continue
            elseif startswith(match, '-') #"-2i" isotope stage
                #Isotope stage can't be lower than -9
                isotope = -getNumeric(match, one(UInt8))
                continue
            end
        end
        
        if startswith(annotation_piece.match, '^') #"^2" charge state
            #Charge state can't be higher than 9
            charge = getNumeric(annotation_piece.match, one(UInt8))
        elseif startswith(annotation_piece.match, '-') #neutral loss "-2H2O"
            sulfur_diff += countSulfurLoss(
                annotation_piece.match
            )
            is_gain_loss = true
        elseif startswith(annotation_piece.match, '+') #neutral gain "+NH3"
            sulfur_diff += countSulfurLoss(
                annotation_piece.match
            ) 
            is_gain_loss = true
        elseif iszero(count(r"^[A-Z]+$", annotation_piece.match)) #not an immonium ion 
            base_type = first(string(annotation_piece.match))
            frag_index = getNumeric(annotation_piece.match, one(UInt8))
        else #is an immonium ion
            #Get number of sulfurs in the immonium ion?
            sulfur_diff += immonium_to_sulfur_count[string(annotation_piece.match)]
        end
    end

    return PioneerFragAnnotation(
        base_type,
        frag_index,
        charge,
        isotope,
        internal,
        immonium,
        is_gain_loss,
        sulfur_diff
    )
end

function parse_fragment_annotation(frag_regex_iterator::GenericFragRegexIterator,
                        immonium_to_sulfur_count::Dict{String, Int8},
                        internal::Bool,
                        immonium::Bool)
    annotation_pieces = getAnnotationPieces(frag_regex_iterator)
    #Initialize Data 
    base_type = '_'
    frag_index = one(UInt8)
    charge = one(UInt8)
    isotope = zero(UInt8)
    sulfur_diff = zero(Int8)
    is_gain_loss = false
    #= Example 
    julia>annotation_pices 
    4-element Vector{RegexMatch}:
    RegexMatch("b11")
    RegexMatch("-2H2O")
    RegexMatch("^3")
    RegexMatch("+2i")
    =#
    for annotation_piece in annotation_pieces
        if startswith(annotation_piece.match, '+') #"^2" charge state
            #Charge state can't be higher than 9
            charge = getNumeric(annotation_piece.match, one(UInt8))
        elseif iszero(count(r"^[A-Z]+$", annotation_piece.match)) #not an immonium ion 
            base_type = first(string(annotation_piece.match))
            frag_index = getNumeric(annotation_piece.match, one(UInt8))
        end
    end

    return PioneerFragAnnotation(
        base_type,
        frag_index,
        charge,
        isotope,
        internal,
        immonium,
        is_gain_loss,
        sulfur_diff
    )
end
#=
golden_regex = r"(?<=[^-+\^]|^)[^\"].*?(?=[-+\^]|$|\")"
collect(eachmatch(golden_regex, "\"ICCAM-NH3^2+i\""))
collect(eachmatch(golden_regex, "\"IRF\""))
collect(eachmatch(golden_regex,  "\"b11-2H2O+2i\""))
collect(eachmatch(golden_regex,  "\"b11-2H2O^3+2i\""))
collect(eachmatch(golden_regex,  "-CO^2+i"))
=#
function getAnnotationPieces(frag_annotation::UniSpecFragAnnotation)
    golden_regex = r"(?<=[^-+\^]|^)[^\"].*?(?=[-+\^]|$|\")"
    return UniSpecFragRegexIterator(eachmatch(golden_regex, getAnnotation(frag_annotation)))
end
function getAnnotationPieces(frag_annotation::GenericFragAnnotation)
    golden_regex = r"(?<=[^-+\^]|^)[^\"].*?(?=[-+\^]|$|\")"
    return GenericFragRegexIterator(eachmatch(golden_regex, getAnnotation(frag_annotation)))
end

"""
    parse_fragment_annotation(annotation::String, immonium_to_sulfur)

Extracts unique ion annotations from an .msp (MSPiano) library and returns a set 
of unique annotations. Special parsing required for internal ions. For example
"\"Int/PDTAENLPAmY-CO^2+i/7\"" poses a problem because there are many unique sequence
and internal annotation could have. We just want to extract the loss/isotope/charge information.
So we would parse "\"Int/PDTAENLPAmY-CO^2+i/7\"" => "\"Int/-CO^2+i\"". 

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
function parse_fragment_annotation(frag_annotation::UniSpecFragAnnotation;
                        immonium_to_sulfur_count::Dict{String, Int8}=Dict{String, Int8}())::PioneerFragAnnotation
    annotation = getAnnotation(frag_annotation)
    is_internal_fragment = false
    if occursin(r"Int", annotation)
        is_internal_fragment = true
        annotation = parseInternalIon(annotation)
    end
    return parse_fragment_annotation(
        getAnnotationPieces(frag_annotation),
        immonium_to_sulfur_count,
        is_internal_fragment,
        occursin(r"I[A-Z]{1,}", annotation),#Is immonium ion 
    )
end

function parse_fragment_annotation(frag_annotation::GenericFragAnnotation;
                        immonium_to_sulfur_count::Dict{String, Int8}=Dict{String, Int8}())::PioneerFragAnnotation
    return parse_fragment_annotation(
        getAnnotationPieces(frag_annotation),
        immonium_to_sulfur_count,
        false,
        false,#Is immonium ion 
    )
end

"""
getAnnotationToID(answer_dict::Dict{String, PioneerFragAnnotation})

Takes a dictionary mapping fragment annotations to `PioneerFragAnnotation` type. Returns 
a vector mapping integer values to the `PioneerFragAnnotation`s. Also returns a dictionary
mapping the fragment annotations (strings) to their integer keys 

### Input
-'answer_dict::Dict{String, PioneerFragAnnotation}' -- maps fragment annotations to `PioneerFragAnnotation`

### Output
A Set of strings 
-'::Vector{PioneerFragAnnotation}'
-'::Dict{String, UInt16}'

### Examples 
julia> test_dict
Dict{String, PioneerFragAnnotation} with 4 entries:
  "\"p-2CH3SOH-H2O\"" => PioneerFragAnnotation('p', 0x00, 0x01, 0x00, false, false, -2)
  "\"IEA\""           => PioneerFragAnnotation('_', 0x00, 0x01, 0x00, false, true, 0)
  "\"y1-NH3\""        => PioneerFragAnnotation('y', 0x01, 0x01, 0x00, false, false, 0)
  "\"Int/-CO^2+i\""   => PioneerFragAnnotation('_', 0x00, 0x02, 0x01, true, false, 0)

julia> id_to_annotation, annotation_to_id = getAnnotationToID(test_dict);

julia> id_to_annotation
4-element Vector{PioneerFragAnnotation}:
 PioneerFragAnnotation('p', 0x00, 0x01, 0x00, false, false, -2)
 PioneerFragAnnotation('_', 0x00, 0x01, 0x00, false, true, 0)
 PioneerFragAnnotation('y', 0x01, 0x01, 0x00, false, false, 0)
 PioneerFragAnnotation('_', 0x00, 0x02, 0x01, true, false, 0)

julia> annotation_to_id
Dict{String, UInt16} with 4 entries:
  "\"p-2CH3SOH-H2O\"" => 0x0001
  "\"IEA\""           => 0x0002
  "\"y1-NH3\""        => 0x0003
  "\"Int/-CO^2+i\""   => 0x0004
"""
function getAnnotationToID(answer_dict::Dict{String, PioneerFragAnnotation})
    id_to_annotation = Vector{PioneerFragAnnotation}(undef, length(answer_dict))
    annotation_to_id = Dict{String, UInt16}()
    i = one(UInt16)
    for (key, value) in pairs(answer_dict)
        id_to_annotation[i] = value
        annotation_to_id[key] = i
        i += one(UInt16)
    end
    return id_to_annotation, annotation_to_id
end


function get_altimeter_ion_dict(ion_table_path::String)
    ion_index_to_name = Dict{Int32, String}()
    open(ion_table_path) do file
        for (i, l) in enumerate(eachline(file))
            ion_name, _, _, _ = split(l, '\t')
            ion_index_to_name[Int32(i-1)] = ion_name
        end
    end

    return ion_index_to_name
end
