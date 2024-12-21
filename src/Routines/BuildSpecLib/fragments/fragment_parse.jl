# src/fragments/fragment_parse.jl

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
    countSulfurs(plain_sequence::AbstractString, 
                    mods_iterator::Base.RegexMatchIterator,
                    mods_to_sufur_diff::Dict{String, Int8})::Int8

Counts the number of Sulfurs in the sequence (assumex M/S are the only Sulfur containing characters). 
Adds or subtracts any sulfurs contributed by the modifications. 

### Input
-'plain_sequence::AbstractString' -- Plaine peptide sequence
-'mods_iterator::Base.RegexMatchIterator' -- Iteratres through sequence modifications that might contribute/subtract sulfurs
-'mods_to_sufur_diff::Dict{String, Int8}' -- maps modification name to the number of sulfurs gained or lost 

### Output
-'::UInt8' Number of Sulfurs in the precursor (if mods_to_sulfur_diff is incorrectly specified there 
might be an error because UInt8 doesn't allow for negative values)

### Examples
julia> countSulfurs(
           "LLDYRTIMHDENK",
           Base.RegexMatchIterator(r"(?<=\\().*?(?=\\))", "1(7,M,Oxidation)(13,K,AnExampleMod)", false),
           Dict("AnExampleMod" => Int8(2))
       )
3
"""
function count_sulfurs!(seq_idx_to_sulfur::Vector{UInt8},
                        plain_sequence::AbstractString, 
                        mods_iterator::Base.RegexMatchIterator,
                        mods_to_sulfur_diff::Dict{String, Int8})::Int8
    sulfur_count = zero(UInt8)
    for (idx, aa) in enumerate(plain_sequence)
        has_sulfur = (aa=='C')|(aa=='M')
        seq_idx_to_sulfur[idx] = UInt8(has_sulfur)
        sulfur_count += has_sulfur
    end

    for mod in mods_iterator
        mod_string = getModName(mod.match)
        if haskey(mods_to_sulfur_diff, mod_string)
            n_sulfur = mods_to_sulfur_diff[mod_string]
            seq_idx_to_sulfur[getModIndex(mod.match)] += n_sulfur
            sulfur_count += n_sulfur
        end
    end
    return sulfur_count
end


function fill_isotope_mods!(seq_idx_to_iso_mod_mass::Vector{Float32},
                        mods_iterator::Base.RegexMatchIterator,
                        iso_mod_to_mass::Dict{String, Float32})

    for mod in mods_iterator
        mod_string = getModName(mod.match)
        if haskey(iso_mod_to_mass, mod_string)
            seq_idx_to_iso_mod_mass[getModIndex(mod.match)] += iso_mod_to_mass[mod_string]
        end
    end
    return nothing
end


function get_fragment_indices(base_type::Char, index::UInt8, sequence_length::UInt8)
    #This ion types run from the N-terminus to the right
    if (base_type=='a')|(base_type=='b')|(base_type=='c')
        return one(UInt8), index
    elseif base_type=='p' #precursor ion covers the entire sequence 
        return one(UInt8), sequence_length
    #Other ion types run from the C-terminus to the left 
    else
        return sequence_length - index + one(UInt8), sequence_length
    end
end

function apply_isotope_mod(charge::UInt8, seq_idx_to_iso_mod::Vector{Float32}, start_idx::UInt8, stop_idx::UInt8)::Float32
    mod_mass = zero(Float32)
    for i in range(start_idx, stop_idx)
        mod_mass += seq_idx_to_iso_mod[i]
    end
    return mod_mass/charge
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
function get_immonium_sulfur_dict(
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


"""
    parse_koina_fragments(
        precursor_table::Arrow.Table,
        fragment_table::Arrow.Table,
        annotation_type::Type{<:FragAnnotation},  # Changed from instance to type
        ion_annotation_set::Set{String},
        frag_name_to_idx::Dict{String, UInt16},
        precursor_batch_size::Int64,
        immonium_data_path::String,
        out_dir::String,
        mods_to_sulfur_diff::Dict{String, Int8},
        iso_mod_to_mass::Dict{String, Float32},
        model_type::KoinaModelType
    )

Process Koina prediction outputs into Pioneer's fragment format.
"""
function parse_koina_fragments(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    annotation_type::FragAnnotation,  # Changed parameter type
    ion_annotation_set::Set{String},
    frag_name_to_idx::Dict{String, UInt16},
    precursor_batch_size::Int64,
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    model_type::KoinaModelType
)
    # Create output paths
    fragment_table_path = joinpath(out_dir, "fragments_table.arrow")
    fragment_indices_path = joinpath(out_dir, "prec_to_frag.arrow")

    # Clean existing files
    rm(fragment_table_path, force=true)
    rm(fragment_indices_path, force=true)

    # Load immonium ion data
    immonium_to_sulfur_count = get_immonium_sulfur_dict(immonium_data_path)

    # Process annotations
    ion_annotation_to_features_dict = Dict{String, PioneerFragAnnotation}()
    annotation_type = typeof(annotation_type)
    for ion in ion_annotation_set
        # Create the annotation instance here
        frag_annotation = annotation_type(ion)  # Now properly constructing instance
        ion_annotation_to_features_dict[ion] = parse_fragment_annotation(
            frag_annotation,
            immonium_to_sulfur_count=immonium_to_sulfur_count
        )
    end
    # Process fragments in batches
    precursor_idx, frag_idx = one(UInt32), one(UInt64)
    total_precursors = length(precursor_table[1])
    
    pbar = ProgressBar(total=ceil(Int64, total_precursors/precursor_batch_size))
    
    while precursor_idx <= total_precursors
        precursor_idx, frag_idx = append_pioneer_lib_batch(
            UInt32(precursor_idx),
            UInt64(frag_idx),
            fragment_table_path,
            fragment_indices_path,
            precursor_batch_size,
            precursor_table,
            fragment_table,
            ion_annotation_to_features_dict,
            frag_name_to_idx,
            mods_to_sulfur_diff,
            iso_mod_to_mass,
            model_type
        )
        update(pbar)
    end
    
    # Write final index
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = UInt64[frag_idx - 1],))
    )

    return ion_annotation_to_features_dict
end

"""
Process fragments in batches to manage memory usage.
"""
function process_fragments_batched(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    ion_annotation_to_features_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    batch_size::Int64,
    fragment_table_path::String,
    fragment_indices_path::String,
    model_type::KoinaModelType
)
    precursor_idx, frag_idx = one(UInt32), one(UInt64)
    total_precursors = length(precursor_table[1])
    
    pbar = ProgressBar(total=ceil(Int64, total_precursors/batch_size))
    
    while precursor_idx <= total_precursors
        precursor_idx, frag_idx = append_pioneer_lib_batch(
            UInt32(precursor_idx),
            UInt64(frag_idx),
            fragment_table_path,
            fragment_indices_path,
            batch_size,
            precursor_table,
            fragment_table,
            ion_annotation_to_features_dict,
            frag_name_to_idx,
            mods_to_sulfur_diff,
            iso_mod_to_mass,
            model_type
        )
        update(pbar)
    end
    
    # Write final index
    Arrow.append(
        fragment_indices_path,
        DataFrame((start_idx = UInt64[frag_idx - 1],))
    )
end

"""
Process a batch of fragments for standard intensity models.
"""
function append_pioneer_lib_batch(
    first_prec_idx::UInt32,
    first_frag_idx::UInt64,
    fragment_table_path::String,
    fragment_indices_path::String,
    n_precursors_in_batch::Int64,
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    ::Union{InstrumentSpecificModel, InstrumentAgnosticModel}
)
    # Count fragments in batch
    n = first_frag_idx
    pid = first_prec_idx
    precs_encountered = 1
    
    while n <= length(fragment_table[:precursor_idx])
        if pid != fragment_table[:precursor_idx][n]
            if precs_encountered == n_precursors_in_batch
                n = n - 1
                pid = fragment_table[:precursor_idx][n]
                break
            end
            pid = fragment_table[:precursor_idx][n]
            precs_encountered += 1
        end
        n += 1
    end
    
    n = min(n, length(fragment_table[:precursor_idx]))
    n_frags = n - first_frag_idx + 1
    n_precursors = precs_encountered
    
    # Allocate storage
    prec_frags = Vector{PioneerFrag}(undef, n_frags)
    pid_to_frag_idxs = Vector{UInt64}(undef, n_precursors)
    prec_sulfur_count = Vector{UInt8}(undef, n_precursors)
    
    # Process each fragment
    process_batch!(
        prec_frags,
        pid_to_frag_idxs,
        prec_sulfur_count,
        first_frag_idx,
        precursor_table,
        fragment_table,
        ion_annotation_to_data_dict,
        frag_name_to_idx,
        mods_to_sulfur_diff,
        iso_mod_to_mass
    )
    
    # Write results
    Arrow.append(fragment_table_path, DataFrame(prec_frags))
    Arrow.append(fragment_indices_path, DataFrame(start_idx = pid_to_frag_idxs))
    
    return first_prec_idx + n_precursors_in_batch, first_frag_idx + n_frags
end

"""
Process a batch of fragments for spline coefficient models.
"""
function append_pioneer_lib_batch(
    first_prec_idx::UInt32,
    first_frag_idx::UInt64,
    fragment_table_path::String,
    fragment_indices_path::String,
    n_precursors_in_batch::Int64,
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    ::SplineCoefficientModel
)
    # Similar to standard batch processing but handles coefficients
    n = first_frag_idx
    pid = first_prec_idx
    precs_encountered = 1
    
    while n <= length(fragment_table[:precursor_idx])
        if pid != fragment_table[:precursor_idx][n]
            if precs_encountered == n_precursors_in_batch
                n = n - 1
                pid = fragment_table[:precursor_idx][n]
                break
            end
            pid = fragment_table[:precursor_idx][n]
            precs_encountered += 1
        end
        n += 1
    end
    
    n = min(n, length(fragment_table[:precursor_idx]))
    n_frags = n - first_frag_idx + 1
    n_precursors = precs_encountered
    
    # Get coefficient tuple type
    coef_type = eltype(fragment_table[:coefficients])
    n_coef = length(coef_type.parameters)
    
    # Allocate storage for spline fragments
    prec_frags = Vector{PioneerSplineFrag{n_coef}}(undef, n_frags)
    pid_to_frag_idxs = Vector{UInt64}(undef, n_precursors)
    prec_sulfur_count = Vector{UInt8}(undef, n_precursors)
    
    # Process batch with coefficients
    process_spline_batch!(
        prec_frags,
        pid_to_frag_idxs,
        prec_sulfur_count,
        first_frag_idx,
        precursor_table,
        fragment_table,
        ion_annotation_to_data_dict,
        frag_name_to_idx,
        mods_to_sulfur_diff,
        iso_mod_to_mass
    )
    
    # Write results
    Arrow.append(fragment_table_path, DataFrame(prec_frags))
    Arrow.append(fragment_indices_path, DataFrame(start_idx = pid_to_frag_idxs))
    
    return first_prec_idx + n_precursors_in_batch, first_frag_idx + n_frags
end

# src/fragments/fragment_parse.jl

"""
Process a batch of fragments, creating Pioneer fragment entries.

Updates the provided vectors in place with fragment information.
"""
function process_batch!(
    prec_frags::Vector{PioneerFrag},
    pid_to_frag_idxs::Vector{UInt64},
    prec_sulfur_count::Vector{UInt8},
    first_frag_idx::UInt64,
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32}
)
    # Track working variables
    seq_idx_to_sulfur = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    last_pid = zero(UInt32)
    precursor_sulfur_count = zero(UInt8)
    precursor_length = zero(UInt8)
    batch_pid = 0
    
    for frag_idx in 1:length(prec_frags)
        # Get fragment info
        actual_frag_idx = first_frag_idx + frag_idx - 1
        frag_annotation = fragment_table[:annotation][actual_frag_idx]
        pid = fragment_table[:precursor_idx][actual_frag_idx]

        # New precursor encountered - update precursor info
        if pid != last_pid
            batch_pid += 1
            last_pid = pid
            
            # Reset and update sulfur tracking
            fill!(seq_idx_to_sulfur, zero(UInt8))
            fill!(seq_idx_to_iso_mod, zero(Float32))
            
            # Calculate sulfur count for new precursor
            prec_sulfur_count[batch_pid] = count_sulfurs!(
                seq_idx_to_sulfur,
                precursor_table[:sequence][pid],
                parseMods(precursor_table[:mods][pid]),
                mods_to_sulfur_diff
            )
            
            # Handle isotope modifications
            iso_mods_iterator = parseMods(precursor_table[:isotope_mods][pid])
            fill_isotope_mods!(
                seq_idx_to_iso_mod,
                iso_mods_iterator,
                iso_mod_to_mass
            )
            
            precursor_sulfur_count = prec_sulfur_count[batch_pid]
            precursor_length = UInt8(length(precursor_table[:sequence][pid]))
            pid_to_frag_idxs[batch_pid] = actual_frag_idx
        end

        # Get fragment data and parse annotation
        frag_data = nothing
        frag_name = nothing
        
        # Handle different fragment types
        if startswith(frag_annotation, "Int")
            frag_name = "Int/" * parse_internal_ion(frag_annotation)
        else
            frag_name = frag_annotation
        end
        #println("frag_annotation $frag_annotation")
        #println("frag_name $frag_name")
        frag_data = ion_annotation_to_data_dict[frag_annotation]

        # Get basic fragment info
        frag_mz = fragment_table[:mz][actual_frag_idx]
        frag_intensity = fragment_table[:intensities][actual_frag_idx]

        # Calculate sequence bounds
        start_idx, stop_idx = get_fragment_indices(
            frag_data.base_type,
            frag_data.frag_index,
            precursor_length
        )

        # Calculate sulfur count for fragment
        sulfur_count = zero(UInt8)

        if !frag_data.immonium
            for i in start_idx:stop_idx
                sulfur_count += seq_idx_to_sulfur[i]
            end
        end

        # Handle internal fragments
        if frag_data.internal
            # Parse position information
            start_idx = parse(UInt8, first(match(r"/([0-9]+)", frag_annotation).captures)) + one(UInt8)
            internal_ion_seq = match(r"(?<=/)[A-Za-z]+(?=[^A-Za-z])", frag_annotation).match
            internal_ion_length = length(internal_ion_seq)
            stop_idx = UInt8(start_idx + internal_ion_length - one(UInt8))
        end

        # Add any sulfur difference from modifications
        sulfur_count += frag_data.sulfur_diff

        # Adjust m/z for isotope modifications
        frag_mz += apply_isotope_mod(
            frag_data.charge,
            seq_idx_to_iso_mod,
            start_idx,
            stop_idx
        )

        # Create fragment entry
        prec_frags[frag_idx] = PioneerFrag(
            # Basic info
            frag_mz,
            Float16(frag_intensity),
            frag_name_to_idx[frag_annotation],
            
            # Ion type flags
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            
            # Other type flags
            ((frag_data.base_type == 'a') |
             (frag_data.base_type == 'x') |
             (frag_data.base_type == 'c') |
             (frag_data.base_type == 'z')),
            frag_data.neutral_diff,
            
            # Fragment details
            frag_data.frag_index,
            frag_data.charge,
            frag_data.isotope,
            frag_data.internal,
            frag_data.immonium,
            
            # Sequence and sulfur info
            (start_idx, stop_idx),
            min(sulfur_count, precursor_sulfur_count)
        )
    end
end

"""
Process a batch of fragments with spline coefficients.

Similar to process_batch! but handles spline coefficients for intensity modeling.
"""
function process_spline_batch!(
    prec_frags::Vector{PioneerSplineFrag{N}},
    pid_to_frag_idxs::Vector{UInt64},
    prec_sulfur_count::Vector{UInt8},
    first_frag_idx::UInt64,
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    ion_annotation_to_data_dict::Dict{String, PioneerFragAnnotation},
    frag_name_to_idx::Dict{String, UInt16},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32}
) where N
    # Working arrays for tracking modifications and sulfur
    seq_idx_to_sulfur = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    
    # Tracking variables
    last_pid = zero(UInt32)
    precursor_sulfur_count = zero(UInt8)
    precursor_length = zero(UInt8)
    batch_pid = 0

    for frag_idx in 1:length(prec_frags)
        actual_frag_idx = first_frag_idx + frag_idx - 1
        frag_annotation = fragment_table[:annotation][actual_frag_idx]
        pid = fragment_table[:precursor_idx][actual_frag_idx]

        # Handle new precursor
        if pid != last_pid
            batch_pid += 1
            frag_idx_start = actual_frag_idx
            last_pid = pid

            # Reset tracking arrays
            fill!(seq_idx_to_sulfur, zero(UInt8))
            fill!(seq_idx_to_iso_mod, zero(Float32))

            # Calculate sulfur count for new precursor
            prec_sulfur_count[batch_pid] = count_sulfurs!(
                seq_idx_to_sulfur,
                precursor_table[:sequence][pid],
                parseMods(precursor_table[:mods][pid]),
                mods_to_sulfur_diff
            )

            # Process isotope modifications
            iso_mods_iterator = parseMods(precursor_table[:isotope_mods][pid])
            fill_isotope_mods!(
                seq_idx_to_iso_mod,
                iso_mods_iterator,
                iso_mod_to_mass
            )

            # Update precursor tracking
            precursor_sulfur_count = prec_sulfur_count[batch_pid]
            precursor_length = UInt8(length(precursor_table[:sequence][pid]))
            pid_to_frag_idxs[batch_pid] = actual_frag_idx
        end

        # Get fragment data
        frag_name = if startswith(frag_annotation, "Int")
            "Int/" * parse_internal_ion(frag_annotation)
        else
            frag_annotation
        end
        frag_data = ion_annotation_to_data_dict[frag_name]

        # Get basic fragment information
        frag_mz = fragment_table[:mz][actual_frag_idx]
        frag_coef = fragment_table[:coefficients][actual_frag_idx]
        frag_intensity = fragment_table[:intensities][actual_frag_idx]

        # Calculate sequence bounds
        start_idx, stop_idx = get_fragment_indices(
            frag_data.base_type,
            frag_data.frag_index,
            precursor_length
        )

        # Calculate sulfur count
        sulfur_count = zero(UInt8)
        if !frag_data.immonium
            for i in start_idx:stop_idx
                sulfur_count += seq_idx_to_sulfur[i]
            end
        end

        # Handle internal fragments
        if frag_data.internal
            start_idx = parse(UInt8, first(match(r"/([0-9]+)", frag_annotation).captures)) + one(UInt8)
            internal_ion_seq = match(r"(?<=/)[A-Za-z]+(?=[^A-Za-z])", frag_annotation).match
            internal_ion_length = length(internal_ion_seq)
            stop_idx = UInt8(start_idx + internal_ion_length - one(UInt8))
        end

        # Add modification-based sulfur changes
        sulfur_count += frag_data.sulfur_diff

        # Adjust m/z for isotope modifications
        frag_mz += apply_isotope_mod(
            frag_data.charge,
            seq_idx_to_iso_mod,
            start_idx,
            stop_idx
        )

        # Create spline fragment entry
        prec_frags[frag_idx] = PioneerSplineFrag(
            # Basic info
            frag_mz,
            frag_coef,  # Spline coefficients
            Float16(frag_intensity),
            frag_name_to_idx[frag_name],
            
            # Ion type flags
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            ((frag_data.base_type == 'a') |
             (frag_data.base_type == 'x') |
             (frag_data.base_type == 'c') |
             (frag_data.base_type == 'z')),
            frag_data.neutral_diff,

            # Fragment details
            frag_data.frag_index,
            frag_data.charge,
            frag_data.isotope,
            frag_data.internal,
            frag_data.immonium,

            # Sequence and sulfur info
            (start_idx, stop_idx),
            min(sulfur_count, precursor_sulfur_count)
        )
    end
end