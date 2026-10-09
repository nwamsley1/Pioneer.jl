# Library precursor metadata access.
#
# LibraryPrecursors wraps an Arrow table of precursor information and provides
# typed getter methods for each column. Also handles CV fold assignment for
# cross-validation (protein-group-based to prevent data leakage).

abstract type LibraryPrecursors end

"""
    PrecursorTextIds

Integer stand-ins for the library's per-precursor text, so intermediate tables can carry
`precursor_idx` alone and compare peptides without strings. Each ID is the rank of its string
among the library's distinct values in `String` sort order, so sorting IDs sorts the strings.

- `sequence_id`: the peptide sequence.
- `mods_id`: the `(structural_mods, isotopic_mods)` pair, with `missing` equal to `""`.
- `species_id` / `species_names`: the canonical species string (`_canonical_species`).
- `accession_set_id`: the `accession_numbers` string. `accession_set_members[id]` lists its
  `;`-separated accessions, sorted and without duplicates, as IDs into `accession_names`
  (the distinct accessions, sorted).
"""
struct PrecursorTextIds
    sequence_id::Vector{UInt32}
    mods_id::Vector{UInt32}
    species_id::Vector{UInt32}
    species_names::Vector{String}
    accession_set_id::Vector{UInt32}
    accession_set_members::Vector{Vector{UInt32}}
    accession_names::Vector{String}
end

struct StandardLibraryPrecursors <: LibraryPrecursors
    data::Arrow.Table
    n::Int64
    accession_numbers_to_pid::Dictionary{String, UInt32}
    pid_to_cv_fold::Vector{UInt8}
    inferred_num_variable_modifications::Union{Nothing, Vector{UInt8}}
    text_ids::Base.RefValue{Union{Nothing, PrecursorTextIds}}
    function StandardLibraryPrecursors(
        precursor_table::Arrow.Table,
        inferred_num_variable_modifications::Union{Nothing, Vector{UInt8}} = nothing
    )
        try
            n = length(precursor_table[:sequence])
            if inferred_num_variable_modifications !== nothing
                length(inferred_num_variable_modifications) == n ||
                    throw(ArgumentError(
                        "inferred_num_variable_modifications must match the precursor table length"
                    ))
            end
            accession_numbers = precursor_table[:accession_numbers]
            accession_keys = String[_accession_key(accs) for accs in accession_numbers]

            unique_proteins = unique(accession_keys)
            accession_number_to_pgid = Dictionary(
                unique_proteins, range(one(UInt32), UInt32(length(unique_proteins)))
            )

            # Assign CV folds by protein group (all precursors of a protein get the same fold)
            pg_to_cv_fold = Dictionary{String, UInt8}()
            cv_folds = UInt8[0, 1]
            Random.seed!(1776)
            for pg in unique_proteins
                insert!(pg_to_cv_fold, pg, rand(cv_folds))
            end
            pid_to_cv_fold = Vector{UInt8}(undef, n)
            for pid in range(1, n)
                pid_to_cv_fold[pid] = pg_to_cv_fold[accession_keys[pid]]
            end
            if length(keys(accession_number_to_pgid)) <= 1
                # Only one (or zero) unique protein in the library — protein-grouped
                # CV folds collapse to a single fold, which defeats the purpose.
                # Fall back to per-precursor random fold assignment.
                n_unique = length(keys(accession_number_to_pgid))
                @user_warn "Library has $n_unique unique protein accession$(n_unique == 1 ? "" : "s"); cannot split CV folds by protein. Falling back to per-precursor random fold assignment."
                for pid in range(1, n)
                    pid_to_cv_fold[pid] = rand(cv_folds)
                end
            end
            new(
                precursor_table,
                n,
                accession_number_to_pgid,
                pid_to_cv_fold,
                inferred_num_variable_modifications,
                Ref{Union{Nothing, PrecursorTextIds}}(nothing)
            )
        catch e
            @user_warn "Failed to load precursor table"
            throw(e)
        end
    end
end

_accession_key(accs::AbstractString) = String(accs)
_accession_key(accs) = join(string.(collect(accs)), ";")

@inline _count_mox(seq::AbstractString) = UInt8(count("Unimod:35", seq))
@inline _count_mox(::Missing) = zero(UInt8)

function _count_variable_modifications(
    structural_mods::Union{Missing, AbstractString},
    variable_mod_names::AbstractSet{String}
)::UInt8
    (ismissing(structural_mods) || isempty(structural_mods) ||
     isempty(variable_mod_names)) && return zero(UInt8)

    n_variable = 0
    for mod_match in eachmatch(r"(?<=\().*?(?=\))", structural_mods)
        mod_name_match = match(r"[^,]+(?=$)", mod_match.match)
        mod_name_match === nothing && continue
        String(mod_name_match.match) in variable_mod_names || continue
        n_variable += 1
    end
    n_variable <= typemax(UInt8) ||
        throw(ArgumentError("A precursor cannot have more than 255 variable modifications"))
    return UInt8(n_variable)
end

_count_variable_modifications(::Any, ::Nothing) = nothing

function _configured_variable_mod_names(config)::Union{Nothing, Set{String}}
    config isa AbstractDict || return nothing
    variable_mods = get(config, "variable_mods", nothing)
    variable_mods isa AbstractDict || return nothing
    names = get(variable_mods, "name", nothing)
    names isa AbstractVector || return nothing
    variable_names = Set(String.(names))
    fixed_mods = get(config, "fixed_mods", nothing)
    fixed_names = if fixed_mods isa AbstractDict
        configured_fixed_names = get(fixed_mods, "name", nothing)
        configured_fixed_names isa AbstractVector ?
            Set(String.(configured_fixed_names)) : Set{String}()
    else
        Set{String}()
    end
    # An old table cannot distinguish fixed and variable instances when a
    # name appears in both lists. Treat that configuration as unavailable and
    # use the silent M-oxidation fallback instead of misclassifying fixed mods.
    isempty(intersect(variable_names, fixed_names)) || return nothing
    return variable_names
end

function SetPrecursors(
    precursor_table::Arrow.Table;
    variable_mod_names::Union{Nothing, AbstractSet{String}} = nothing
)
    inferred_num_variable_modifications = if !hasproperty(
        precursor_table,
        :num_variable_modifications
    ) && variable_mod_names !== nothing
        UInt8[
            _count_variable_modifications(mods, variable_mod_names)
            for mods in precursor_table[:structural_mods]
        ]
    else
        nothing
    end
    return StandardLibraryPrecursors(
        precursor_table,
        inferred_num_variable_modifications
    )
end

Base.length(lp::LibraryPrecursors) = lp.n

# ============================================================================
# Precursor column getters
# ============================================================================

getCvFold(lp::LibraryPrecursors, precursor_idx::I) where {I<:Integer} = lp.pid_to_cv_fold[precursor_idx]
getProteinGroupId(lp::LibraryPrecursors, accession_numbers::String)::UInt32 = lp.accession_numbers_to_pid[accession_numbers]
getProteomeIdentifiers(lp::LibraryPrecursors)::Arrow.List{S,Int32,Array{UInt8,1}} where {S<:AbstractString} = lp.data[:proteome_identifiers]
getAccessionNumbers(lp::LibraryPrecursors)::Arrow.List{String, Int32, Vector{UInt8}} = lp.data[:accession_numbers]
getSequence(lp::LibraryPrecursors)::Arrow.List{S,Int32,Array{UInt8,1}} where {S<:AbstractString} = lp.data[:sequence]
getStructuralMods(lp::LibraryPrecursors)::Arrow.List{Union{Missing, String}, Int32, Vector{UInt8}} = lp.data[:structural_mods]
getCharge(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:prec_charge]
getIsDecoy(lp::LibraryPrecursors)::Arrow.BoolVector{Bool} = lp.data[:is_decoy]
getEntrapmentGroupId(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:entrapment_group_id]
getMz(lp::LibraryPrecursors)::Arrow.Primitive{Float32, Vector{Float32}} = lp.data[:mz]
getLength(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:length]
getMissedCleavages(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:missed_cleavages]
getStartIdx(lp::LibraryPrecursors) = lp.data[:start_idx]

@inline _format_start_idx(start::Integer) = string(start)
@inline _format_start_idx(starts) = join(starts, ';')

"""Species as written to the outputs: the distinct proteome identifiers, sorted and `;`-joined."""
_canonical_species(proteome) = join(sort(unique(split(coalesce(proteome, ""), ';'))), ';')

# Replace each value by its rank among the distinct values in sort order.
function _rank_ids(values, ::Type{K}) where {K}
    first_seen = Dict{K, UInt32}()
    ids = Vector{UInt32}(undef, length(values))
    @inbounds for i in eachindex(values)
        ids[i] = get!(first_seen, convert(K, values[i]), UInt32(length(first_seen) + 1))
    end
    distinct = Vector{K}(undef, length(first_seen))
    for (value, id) in first_seen
        distinct[id] = value
    end
    order = sortperm(distinct)
    rank = Vector{UInt32}(undef, length(distinct))
    rank[order] = UInt32.(1:length(distinct))
    @inbounds for i in eachindex(ids)
        ids[i] = rank[ids[i]]
    end
    return ids, distinct[order]
end

"""
    getPrecursorTextIds(lp::LibraryPrecursors) -> PrecursorTextIds

Build (once, on first use) and return the integer IDs for the library's text columns.
"""
function getPrecursorTextIds(lp::StandardLibraryPrecursors)
    ids = lp.text_ids[]
    ids === nothing || return ids
    sequence_id, _ = _rank_ids(getSequence(lp), String)
    structural_mods, isotopic_mods = getStructuralMods(lp), getIsotopicMods(lp)
    mods = [(coalesce(structural_mods[i], ""), coalesce(isotopic_mods[i], "")) for i in 1:length(lp)]
    mods_id, _ = _rank_ids(mods, Tuple{String, String})
    canonical = Dict{Union{Missing, String}, String}()
    species = [get!(() -> _canonical_species(p), canonical, p) for p in getProteomeIdentifiers(lp)]
    species_id, species_names = _rank_ids(species, String)
    accession_set_id, accession_sets = _rank_ids(getAccessionNumbers(lp), String)
    accession_names = sort!(unique!(String[a for set in accession_sets for a in split(set, ';')]))
    accession_rank = Dict(name => UInt32(i) for (i, name) in enumerate(accession_names))
    accession_set_members = [sort!(unique!(UInt32[accession_rank[a] for a in split(set, ';')]))
                             for set in accession_sets]
    lp.text_ids[] = PrecursorTextIds(sequence_id, mods_id, species_id, species_names,
                                     accession_set_id, accession_set_members, accession_names)
    return lp.text_ids[]
end

"""
    PrecursorOutputText(precursors, file_names, protein_group_names = String[])

The lookups needed to put text back on precursor output rows: library text by
`:precursor_idx`, file names by `:ms_file_idx`, and protein-group names by the `pg_id` in
`:inferred_protein_group`.
"""
struct PrecursorOutputText{P<:LibraryPrecursors}
    precursors::P
    text_ids::PrecursorTextIds
    file_names::Vector{String}
    protein_group_names::Vector{String}
end

PrecursorOutputText(precursors::LibraryPrecursors, file_names::AbstractVector{<:AbstractString},
                    protein_group_names::AbstractVector{<:AbstractString} = String[]) =
    PrecursorOutputText(precursors, getPrecursorTextIds(precursors), String.(file_names),
                        String.(protein_group_names))

const PRECURSOR_TEXT_COLUMNS_BEFORE_NORMALIZED = (:accession_numbers, :species)
const PRECURSOR_TEXT_COLUMNS_AFTER_NORMALIZED =
    (:structural_mods, :isotopic_mods, :sequence, :peptide_start_positions, :file_name)
const PRECURSOR_TEXT_COLUMNS =
    (PRECURSOR_TEXT_COLUMNS_BEFORE_NORMALIZED..., PRECURSOR_TEXT_COLUMNS_AFTER_NORMALIZED...)

function getNumEnzymaticTermini(lp::LibraryPrecursors)
    hasproperty(lp.data, :num_enzymatic_termini) &&
        return lp.data[:num_enzymatic_termini]

    # Libraries built before enzymatic specificity was recorded were fully
    # specific, so both termini are enzymatic by construction.
    return fill(UInt8(2), length(lp))
end

"""
    getNumVariableModifications(lp::LibraryPrecursors)

Return exact variable-modification counts stored in the library or inferred
from its build configuration. Returns `nothing` for legacy libraries without
either source of metadata; callers silently retain the historical
M-oxidation-only behavior for those libraries.
"""
getNumVariableModifications(lp::LibraryPrecursors) =
    hasproperty(lp.data, :num_variable_modifications) ?
        lp.data[:num_variable_modifications] :
        lp.inferred_num_variable_modifications

@inline _num_variable_modifications_at(::Nothing, structural_mods, precursor_idx::Integer) =
    _count_mox(structural_mods[precursor_idx])
@inline _num_variable_modifications_at(values, structural_mods, precursor_idx::Integer) =
    UInt8(values[precursor_idx])
getIrt(lp::LibraryPrecursors)::Arrow.Primitive{Float32, Vector{Float32}} = lp.data[:irt]
# Ion-mobility columns exist only in libraries built with `library_params.im_model`;
# `nothing` otherwise.
getCCS(lp::LibraryPrecursors) = hasproperty(lp.data, :ccs) ? lp.data[:ccs] : nothing
getInvIonMobility(lp::LibraryPrecursors) =
    hasproperty(lp.data, :inv_ion_mobility) ? lp.data[:inv_ion_mobility] : nothing
getSulfurCount(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:sulfur_count]
"Isotopic mods per precursor; all missing for libraries built without isotope groups, which leave the column out."
getIsotopicMods(lp::LibraryPrecursors) = hasproperty(lp.data, :isotopic_mods) ?
    lp.data[:isotopic_mods]::Arrow.List{Union{Missing, String}, Int32, Vector{UInt8}} : AllMissingStrings(length(lp))

"A read-only column of `n` missing strings."
struct AllMissingStrings <: AbstractVector{Union{Missing, String}}
    n::Int
end
Base.size(v::AllMissingStrings) = (v.n,)
Base.getindex(v::AllMissingStrings, i::Int) = (@boundscheck checkbounds(v, i); missing)
getBasePepId(lp::LibraryPrecursors)::Arrow.Primitive{UInt32, Vector{UInt32}} = lp.data[:base_pep_id]

# ============================================================================
# Pair index helpers (target/decoy pairing)
# ============================================================================

getPairIdx(lp::LibraryPrecursors) = lp.data[:pair_id]

function extract_pair_idx(pair_idx_column::Arrow.Primitive{Union{Missing, UInt32}, Array{UInt32,1}}, idx)
    value = pair_idx_column[idx]
    return ismissing(value) ? zero(UInt32) : convert(UInt32, value)
end

function extract_pair_idx(pair_idx_column::Arrow.Primitive{UInt32, Vector{UInt32}}, idx)
    return pair_idx_column[idx]
end

function extract_pair_idx(pair_idx_column, idx)
    value = pair_idx_column[idx]
    if isa(value, UInt32)
        return value
    elseif ismissing(value)
        return zero(UInt32)
    else
        return convert(UInt32, value)
    end
end
