"""
    parseAnnotation(annotation_pieces::Vector{RegexMatch},internal::Bool,immonium::Bool)

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
function parseAnnotation(frag_regex_iterator::UniSpecFragRegexIterator,
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
        if endswith(annotation_piece.match,'i') #"+2i" isotope_stage
            #making an assumption here that there are no "-"
            isotope = getNumeric(annotation_piece.match, one(UInt8))
        elseif startswith(annotation_piece.match, '^') #"^2" charge state
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

function parseAnnotation(frag_regex_iterator::GenericFragRegexIterator,
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
    parseAnnotation(annotation::String, immonium_to_sulfur)

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
function parseAnnotation(frag_annotation::UniSpecFragAnnotation;
                        immonium_to_sulfur_count::Dict{String, Int8}=Dict{String, Int8}())::PioneerFragAnnotation
    annotation = getAnnotation(frag_annotation)
    is_internal_fragment = false
    if occursin(r"Int", annotation)
        is_internal_fragment = true
        annotation = parseInternalIon(annotation)
    end
    return parseAnnotation(
        getAnnotationPieces(frag_annotation),
        immonium_to_sulfur_count,
        is_internal_fragment,
        occursin(r"I[A-Z]{2,}", annotation),#Is immonium ion 
    )
end

function parseAnnotation(frag_annotation::GenericFragAnnotation;
                        immonium_to_sulfur_count::Dict{String, Int8}=Dict{String, Int8}())::PioneerFragAnnotation
    return parseAnnotation(
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