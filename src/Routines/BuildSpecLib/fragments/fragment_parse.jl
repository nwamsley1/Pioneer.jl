# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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