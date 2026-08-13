# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

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
julia> getImmoniumToFormulaDict("assets/immonium.txt")
Dict{String, Int8} with 64 entries:
  "IQ" => 0x0
  "IW" => 0x0
  "IL" => 0x0
  "IQ" => 0x0
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
    ion_annotation_to_data_dict::Dict{Int32, PioneerFragAnnotation},
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
        mods_to_sulfur_diff,
        iso_mod_to_mass
    )
    
    # Write results
    Arrow.append(fragment_table_path, DataFrame(prec_frags))
    Arrow.append(fragment_indices_path, DataFrame(start_idx = pid_to_frag_idxs))
    
    return first_prec_idx + n_precursors_in_batch, first_frag_idx + n_frags
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
    ion_annotation_to_data_dict::Dict{Int32, PioneerFragAnnotation},
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
        rank = Float16(255)
        if pid != last_pid
            rank -= one(Float16)
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
        #frag_name = if startswith(frag_annotation, "Int")
        #    "Int/" * parse_internal_ion(frag_annotation)
        #else
        #    frag_annotation
        #end
        frag_data = ion_annotation_to_data_dict[frag_annotation]

        # Get basic fragment information
        frag_mz = fragment_table[:mz][actual_frag_idx]
        frag_coef = fragment_table[:coefficients][actual_frag_idx]
        frag_intensity = rank #- one(Float32)#fragment_table[:intensities][actual_frag_idx]

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
        #if frag_data.internal
        #    start_idx = parse(UInt8, first(match(r"/([0-9]+)", frag_annotation).captures)) + one(UInt8)
        #    internal_ion_seq = match(r"(?<=/)[A-Za-z]+(?=[^A-Za-z])", frag_annotation).match
        #    internal_ion_length = length(internal_ion_seq)
        #    stop_idx = UInt8(start_idx + internal_ion_length - one(UInt8))
        #end

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
            UInt16(frag_annotation),
            
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


function parse_altimeter_fragments(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    annotation_type::FragAnnotation,  # Changed parameter type
    ion_dictionary::Dict{Int32, String},
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
    ion_annotation_to_features_dict = Dict{Int32, PioneerFragAnnotation}()
    annotation_type = typeof(annotation_type)
    for (ion_idx, ion_name) in ion_dictionary
        # Create the annotation instance here
        frag_annotation = annotation_type(ion_name)  # Now properly constructing instance
        ion_annotation_to_features_dict[ion_idx] = parse_fragment_annotation(
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
            mods_to_sulfur_diff,
            iso_mod_to_mass,
            model_type
        )
        update(pbar)
    end
    
    # Discard the per-batch incrementally-built `prec_to_frag.arrow` and
    # rebuild it from the on-disk raw fragments table. The per-batch writer
    # in `process_spline_batch!` indexes `pid_to_frag_idxs` by *appearance
    # order* (it bumps `batch_pid` only on observed pid transitions), so any
    # precursor whose fragments were ALL filtered out upstream — possible
    # since §9 moved the metadata + rank-cap filters into
    # `filter_fragments!(::SplineCoefficientModel)` — ends up with an
    # uninitialized slot (a stray 0) that `getDetailedFrags` later reads as
    # an arbitrary fragment range and walks past the end of `frag_mz`.
    #
    # Rebuilding from `raw_fragments.arrow` (which still carries
    # `precursor_idx` per fragment, and whose row order matches
    # `fragments_table.arrow` 1:1 since `process_spline_batch!` is a pure
    # 1:1 copy) guarantees one start index per logical precursor, with
    # empty CSR ranges (`start[i+1] == start[i]`) for precursors that
    # have zero surviving fragments.
    rebuild_prec_to_frag_index!(
        fragment_indices_path,
        joinpath(out_dir, "raw_fragments.arrow"),
        length(precursor_table[1]),
    )

    serialize_to_jls(joinpath(out_dir, "frag_name_to_idx.jls"), ion_dictionary)
    serialize_to_jls(joinpath(out_dir, "ion_annotations.jls"), ion_annotation_to_features_dict)
    return ion_annotation_to_features_dict
end

"""
    rebuild_prec_to_frag_index!(prec_to_frag_path, raw_fragments_path, n_precursors)

Recompute `prec_to_frag.arrow` from `raw_fragments.arrow[:precursor_idx]`,
producing one start index per logical precursor 1..n_precursors plus a
final boundary at `length(raw_fragments) + 1`. Precursors with zero
fragments get the same start index as the next non-empty precursor
(empty CSR range). Overwrites `prec_to_frag_path`.
"""
function rebuild_prec_to_frag_index!(prec_to_frag_path::String,
                                     raw_fragments_path::String,
                                     n_precursors::Integer)
    raw_pids = Arrow.Table(raw_fragments_path).precursor_idx
    n_frags = length(raw_pids)
    starts = fill(zero(UInt64), n_precursors + 1)
    @inbounds for i in 1:n_frags
        pid = Int(raw_pids[i])
        if 1 <= pid <= n_precursors && starts[pid] == 0
            starts[pid] = UInt64(i)
        end
    end
    starts[end] = UInt64(n_frags + 1)
    # Back-fill from the right so empty-range precursors inherit the next
    # valid start. Readers compute the stop index as
    # `starts[pid+1] - 1`, so an empty range emerges naturally as
    # `starts[pid] == starts[pid+1]`.
    next = starts[end]
    for i in n_precursors:-1:1
        if starts[i] == 0
            starts[i] = next
        else
            next = starts[i]
        end
    end
    rm(prec_to_frag_path; force = true)
    Arrow.write(prec_to_frag_path, DataFrame(start_idx = starts))
    return nothing
end

"""
    build_detailed_frags_from_raw(precursor_table, fragment_table, annotation_type,
        ion_dictionary, immonium_data_path, out_dir, mods_to_sulfur_diff,
        iso_mod_to_mass, ::SplineCoefficientModel) -> (detailed_frags, prec_to_frag)

Fused replacement for `parse_altimeter_fragments` + `getDetailedFrags`: decode the
raw Koina fragment table directly into the final `Vector{SplineCompactFrag}` (plus
the CSR precursor→fragment index) in a single pass, WITHOUT the intermediate
`fragments_table.arrow` re-encode/re-read.

It reproduces exactly the per-fragment work of `process_spline_batch!` (sulfur count
over the fragment's sequence span, capped at the precursor sulfur count; isotope-mod
m/z adjustment; ion-type flags) and the per-precursor positional rank of
`getDetailedFrags`. `fragment_table` is assumed grouped by ascending `precursor_idx`
(predict_fragments emits/filters in that order), so a single sequential pass yields
contiguous per-precursor ranges; the CSR is built first-occurrence-only +
right-back-filled, identical to `rebuild_prec_to_frag_index!`. `detailed_frags[fi]`
is 1:1 with raw row `fi` (no filtering happens here — it already happened in
`filter_fragments!`), so positions match the two-step path.

Still serializes `frag_name_to_idx.jls` + `ion_annotations.jls` (as the two-step path
did). Returns `(detailed_frags, prec_to_frag)` to be handed to `buildPionLib` in
memory; nothing fragment-related is written to disk here.
"""
function build_detailed_frags_from_raw(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    annotation_type::FragAnnotation,
    ion_dictionary::Dict{Int32, String},
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    ::SplineCoefficientModel
)
    # Annotation feature dict — identical to parse_altimeter_fragments.
    immonium_to_sulfur_count = get_immonium_sulfur_dict(immonium_data_path)
    ion_annotation_to_data_dict = Dict{Int32, PioneerFragAnnotation}()
    atype = typeof(annotation_type)
    for (ion_idx, ion_name) in ion_dictionary
        ion_annotation_to_data_dict[ion_idx] = parse_fragment_annotation(
            atype(ion_name); immonium_to_sulfur_count = immonium_to_sulfur_count)
    end

    coef_col = fragment_table[:coefficients]
    coef_type = eltype(coef_col)
    N = length(coef_type.parameters)
    T = eltype(coef_type)
    n_frags  = length(fragment_table[:mz])
    n_precursors = length(precursor_table[1])
    detailed = Vector{SplineCompactFrag{N, T}}(undef, n_frags)
    prec_to_frag = fill(zero(UInt64), n_precursors + 1)

    # Function barrier: the Arrow.Table columns are abstractly typed at this call
    # site (Arrow.Table getindex returns Any), so pass them to a kernel that
    # specializes on their concrete element types (ChainedVector{Int32}, ...) for
    # a type-stable hot loop — the pattern getDetailedFrags gets from its signature.
    _fill_detailed_from_raw!(
        detailed, prec_to_frag,
        fragment_table[:annotation], fragment_table[:precursor_idx],
        fragment_table[:mz], coef_col,
        precursor_table[:sequence], precursor_table[:mods],
        precursor_table[:isotope_mods], precursor_table[:precursor_charge],
        ion_annotation_to_data_dict, mods_to_sulfur_diff, iso_mod_to_mass,
        n_frags, n_precursors)

    serialize_to_jls(joinpath(out_dir, "frag_name_to_idx.jls"), ion_dictionary)
    serialize_to_jls(joinpath(out_dir, "ion_annotations.jls"), ion_annotation_to_data_dict)
    return detailed, prec_to_frag
end

# Type-stable kernel for build_detailed_frags_from_raw. Reproduces
# process_spline_batch!'s per-fragment decode (sulfur over the fragment's sequence
# span capped at the precursor sulfur; isotope-mod m/z shift; ion flags) and
# getDetailedFrags's per-precursor positional rank, emitting SplineCompactFrag at
# detailed[fi] (1:1 with raw row fi). Builds the CSR prec_to_frag first-occurrence-
# only + right-back-filled, identical to rebuild_prec_to_frag_index!.
function _fill_detailed_from_raw!(
    detailed::Vector{SplineCompactFrag{N, T}},
    prec_to_frag::Vector{UInt64},
    frag_ann, frag_pid, frag_mzs, coef_col,
    seqs, mods_col, iso_mod_col, charge_col,
    ion_dict::Dict{Int32, PioneerFragAnnotation},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    n_frags::Int, n_precursors::Int) where {N, T}

    seq_idx_to_sulfur  = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    last_pid = zero(UInt32)
    prec_sulfur_count = 0
    prec_length = zero(UInt8)
    prec_charge = zero(UInt8)
    rank = 0

    @inbounds for fi in 1:n_frags
        pid = frag_pid[fi]
        if pid != last_pid
            last_pid = pid
            rank = 0
            if prec_to_frag[pid] == 0
                prec_to_frag[pid] = UInt64(fi)
            end
            fill!(seq_idx_to_sulfur, zero(UInt8))
            fill!(seq_idx_to_iso_mod, zero(Float32))
            seq = seqs[pid]
            prec_sulfur_count = Int(count_sulfurs!(
                seq_idx_to_sulfur, seq, parseMods(mods_col[pid]), mods_to_sulfur_diff))
            fill_isotope_mods!(seq_idx_to_iso_mod, parseMods(iso_mod_col[pid]), iso_mod_to_mass)
            prec_length = UInt8(length(seq))
            prec_charge = UInt8(charge_col[pid])
        end
        rank += 1

        frag_data = ion_dict[frag_ann[fi]]
        start_idx, stop_idx = get_fragment_indices(
            frag_data.base_type, frag_data.frag_index, prec_length)

        sulfur_count = 0
        if !frag_data.immonium
            for i in start_idx:stop_idx
                sulfur_count += seq_idx_to_sulfur[i]
            end
        end
        sulfur_count += frag_data.sulfur_diff

        frag_mz = frag_mzs[fi] + apply_isotope_mod(
            frag_data.charge, seq_idx_to_iso_mod, start_idx, stop_idx)

        detailed[fi] = SplineCompactFrag(
            UInt32(pid),
            frag_mz,
            coef_col[fi],
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            frag_data.isotope > 0,
            frag_data.charge,
            frag_data.frag_index,
            prec_charge,
            UInt8(rank),
            UInt8(min(sulfur_count, prec_sulfur_count)),
        )
    end

    # CSR boundaries: final boundary + right-back-fill so empty precursors inherit
    # the next non-empty start (identical to rebuild_prec_to_frag_index!).
    prec_to_frag[end] = UInt64(n_frags + 1)
    next = prec_to_frag[end]
    for i in n_precursors:-1:1
        if prec_to_frag[i] == 0
            prec_to_frag[i] = next
        else
            next = prec_to_frag[i]
        end
    end
    return nothing
end

"""
    build_detailed_frags_from_raw(precursor_table, fragment_table, annotation_type,
        immonium_data_path, out_dir, mods_to_sulfur_diff, iso_mod_to_mass,
        ::InstrumentAgnosticModel) -> (detailed_frags, prec_to_frag)

Prosit analog of the spline fused decode: decode the topN-filtered
`raw_fragments.arrow` directly into `Vector{CompactFrag}` (constant scalar
intensity) + the CSR precursor→fragment index, in one pass. Differs from the
spline path only in that Prosit carries STRING annotations (parsed via the
`GenericFragAnnotation` parser, cached) and a scalar `:intensities` column that
already holds *total* abundance (the mono→total conversion happened in
`filter_fragments!(::InstrumentAgnosticModel)`). Per-fragment sulfur is recomputed
here with the same `count_sulfurs!`/`get_fragment_indices` used at predict time
(deterministic — keeps this decode symmetric with the spline path, no sulfur column
in `raw_fragments.arrow`). Serializes `ion_annotations.jls` for parity.
"""
function build_detailed_frags_from_raw(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    annotation_type::FragAnnotation,
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    ::InstrumentAgnosticModel
)
    n_frags = length(fragment_table[:mz])
    n_precursors = length(precursor_table[1])
    detailed = Vector{CompactFrag{Float32}}(undef, n_frags)
    prec_to_frag = fill(zero(UInt64), n_precursors + 1)
    annotation_cache = Dict{String, PioneerFragAnnotation}()

    _fill_compact_from_raw!(
        detailed, prec_to_frag,
        fragment_table[:annotation], fragment_table[:precursor_idx],
        fragment_table[:mz], fragment_table[:intensities],
        precursor_table[:sequence], precursor_table[:mods],
        precursor_table[:isotope_mods], precursor_table[:precursor_charge],
        annotation_cache, mods_to_sulfur_diff, iso_mod_to_mass,
        n_frags, n_precursors)

    serialize_to_jls(joinpath(out_dir, "ion_annotations.jls"), annotation_cache)
    return detailed, prec_to_frag
end

# Type-stable kernel for the Prosit (InstrumentAgnosticModel) fused decode.
# Mirrors _fill_detailed_from_raw! (same per-precursor sulfur/iso-mod setup, same
# per-fragment span sulfur, same CSR back-fill, same per-precursor positional rank)
# but reads a scalar `:intensities` (already total abundance) and emits CompactFrag
# with a STRING-annotation lookup instead of the Int32 ion-dict.
function _fill_compact_from_raw!(
    detailed::Vector{CompactFrag{T}},
    prec_to_frag::Vector{UInt64},
    frag_ann, frag_pid, frag_mzs, frag_ints,
    seqs, mods_col, iso_mod_col, charge_col,
    annotation_cache::Dict{String, PioneerFragAnnotation},
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    n_frags::Int, n_precursors::Int) where {T}

    seq_idx_to_sulfur  = zeros(UInt8, 255)
    seq_idx_to_iso_mod = zeros(Float32, 255)
    last_pid = zero(UInt32)
    prec_sulfur_count = 0
    prec_length = zero(UInt8)
    prec_charge = zero(UInt8)
    rank = 0

    @inbounds for fi in 1:n_frags
        pid = frag_pid[fi]
        if pid != last_pid
            last_pid = pid
            rank = 0
            if prec_to_frag[pid] == 0
                prec_to_frag[pid] = UInt64(fi)
            end
            fill!(seq_idx_to_sulfur, zero(UInt8))
            fill!(seq_idx_to_iso_mod, zero(Float32))
            seq = seqs[pid]
            prec_sulfur_count = Int(count_sulfurs!(
                seq_idx_to_sulfur, seq, parseMods(mods_col[pid]), mods_to_sulfur_diff))
            fill_isotope_mods!(seq_idx_to_iso_mod, parseMods(iso_mod_col[pid]), iso_mod_to_mass)
            prec_length = UInt8(length(seq))
            prec_charge = UInt8(charge_col[pid])
        end
        rank += 1

        frag_data = _agnostic_annotation!(annotation_cache, frag_ann[fi])
        start_idx, stop_idx = get_fragment_indices(
            frag_data.base_type, frag_data.frag_index, prec_length)

        sulfur_count = 0
        if !frag_data.immonium
            for i in start_idx:stop_idx
                sulfur_count += seq_idx_to_sulfur[i]
            end
        end
        sulfur_count += frag_data.sulfur_diff

        frag_mz = frag_mzs[fi] + apply_isotope_mod(
            frag_data.charge, seq_idx_to_iso_mod, start_idx, stop_idx)

        detailed[fi] = CompactFrag(
            UInt32(pid),
            frag_mz,
            Float16(frag_ints[fi]),
            frag_data.base_type == 'y',
            frag_data.base_type == 'b',
            frag_data.base_type == 'p',
            frag_data.isotope > 0,
            frag_data.charge,
            frag_data.frag_index,
            prec_charge,
            UInt8(rank),
            UInt8(min(sulfur_count, prec_sulfur_count)),
        )
    end

    prec_to_frag[end] = UInt64(n_frags + 1)
    next = prec_to_frag[end]
    for i in n_precursors:-1:1
        if prec_to_frag[i] == 0
            prec_to_frag[i] = next
        else
            next = prec_to_frag[i]
        end
    end
    return nothing
end