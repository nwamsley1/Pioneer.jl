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

"""
    _register_protein_ambiguities!(candidates_by_id, inference_result, last_id)

Register ambiguous peptide assignments in memory and return the peptide-to-ID lookup plus the
final assigned ID.
"""
function _register_protein_ambiguities!(
    candidates_by_id::Dict{UInt32, Vector{P}},
    inference_result::InferenceResult{P, Q},
    last_id::UInt32
) where {P, Q}
    peptide_to_id = Dictionary{Q, UInt32}()
    peptide_keys = sort!(collect(keys(inference_result.ambiguous_peptide_to_proteins)))

    for peptide_key in peptide_keys
        last_id == typemax(UInt32) && error("Protein ambiguity ID space exhausted")
        last_id += one(UInt32)
        insert!(peptide_to_id, peptide_key, last_id)

        candidates = sort!(unique(copy(
            inference_result.ambiguous_peptide_to_proteins[peptide_key]
        )))
        candidates_by_id[last_id] = candidates
    end

    return peptide_to_id, last_id
end

function _collect_final_protein_groups!(
    final_groups::Set{P},
    inference_result::InferenceResult{P}
) where {P}
    for protein in values(inference_result.peptide_to_protein)
        push!(final_groups, protein)
    end
    for candidates in values(inference_result.ambiguous_peptide_to_proteins)
        union!(final_groups, candidates)
    end
    return final_groups
end

"""
    count_protein_peptide_opportunities(
        accession_numbers,
        sequences,
        is_decoys,
        entrap_ids,
        final_groups
    )

Classify distinct library peptide sequences as unique or shared relative to the
retained final protein groups. Sharing among accessions within one final group
is unique-to-group; only mappings to multiple final groups are shared.
"""
function count_protein_peptide_opportunities(
    accession_numbers::AbstractVector{<:AbstractString},
    sequences::AbstractVector{<:AbstractString},
    is_decoys::AbstractVector{Bool},
    entrap_ids::AbstractVector{UInt8},
    final_groups::Set{ProteinKey{String}};
    common_precursor_mask::Union{Nothing, AbstractVector{Bool}} = nothing
)
    return _count_protein_peptide_opportunities(
        accession_numbers, sequences, is_decoys, entrap_ids, final_groups,
        accessions -> (String(strip(a)) for a in split(accessions, ';')),
        group -> (String(strip(a)) for a in split(group.name, ';'));
        common_precursor_mask = common_precursor_mask)
end

"""
Generic core: `row_accessions(accession_numbers[i])` and `group_accessions(group)` give the
accessions of a library row and of a final group, as comparable keys.
"""
function _count_protein_peptide_opportunities(
    accession_numbers::AbstractVector,
    sequences::AbstractVector,
    is_decoys::AbstractVector{Bool},
    entrap_ids::AbstractVector{UInt8},
    final_groups::Set{G},
    row_accessions,
    group_accessions;
    common_precursor_mask::Union{Nothing, AbstractVector{Bool}} = nothing
) where {G<:ProteinKey}
    n_rows = length(accession_numbers)
    length(sequences) == n_rows ||
        throw(ArgumentError("sequences must match accession_numbers length"))
    length(is_decoys) == n_rows ||
        throw(ArgumentError("is_decoys must match accession_numbers length"))
    length(entrap_ids) == n_rows ||
        throw(ArgumentError("entrap_ids must match accession_numbers length"))
    common_precursor_mask === nothing ||
        length(common_precursor_mask) == n_rows ||
        throw(ArgumentError(
            "common_precursor_mask must match accession_numbers length"
        ))

    A = eltype(accession_numbers) <: AbstractString ? String : eltype(accession_numbers)
    K = eltype(sequences) <: AbstractString ? String : eltype(sequences)
    # Accession sets and single accessions share a key type (String, or UInt32 IDs).
    accession_to_groups = Dict{Tuple{A, Bool, UInt8}, Vector{G}}()
    for group in final_groups
        for accession in group_accessions(group)
            key = (accession, group.is_target, group.entrap_id)
            push!(get!(accession_to_groups, key, G[]), group)
        end
    end
    for groups in values(accession_to_groups)
        sort!(unique!(groups))
    end

    group_cache =
        Dict{Tuple{A, Bool, UInt8}, Vector{G}}()
    peptide_to_groups =
        Dict{Tuple{K, Bool, UInt8}, Vector{G}}()
    common_peptide_to_groups =
        Dict{Tuple{K, Bool, UInt8}, Vector{G}}()

    @inbounds for i in eachindex(accession_numbers, sequences, is_decoys, entrap_ids)
        target = !is_decoys[i]
        entrap_id = entrap_ids[i]
        accessions = convert(A, accession_numbers[i])
        groups = get!(group_cache, (accessions, target, entrap_id)) do
            mapped_groups = G[]
            for accession in row_accessions(accessions)
                append!(
                    mapped_groups,
                    get(
                        accession_to_groups,
                        (accession, target, entrap_id),
                        G[]
                    )
                )
            end
            sort!(unique!(mapped_groups))
        end

        isempty(groups) && continue
        peptide_key = (convert(K, sequences[i]), target, entrap_id)
        if !haskey(peptide_to_groups, peptide_key)
            peptide_to_groups[peptide_key] = groups
        elseif peptide_to_groups[peptide_key] != groups
            peptide_to_groups[peptide_key] =
                sort!(unique!(vcat(peptide_to_groups[peptide_key], groups)))
        end
        is_common = common_precursor_mask === nothing || common_precursor_mask[i]
        if is_common
            if !haskey(common_peptide_to_groups, peptide_key)
                common_peptide_to_groups[peptide_key] = groups
            elseif common_peptide_to_groups[peptide_key] != groups
                common_peptide_to_groups[peptide_key] = sort!(unique!(vcat(
                    common_peptide_to_groups[peptide_key],
                    groups
                )))
            end
        end
    end

    unique_counts = Dict(group => 0 for group in final_groups)
    shared_counts = Dict(group => 0 for group in final_groups)
    common_unique_counts = Dict(group => 0 for group in final_groups)
    for (peptide_key, groups) in pairs(peptide_to_groups)
        counts = length(groups) == 1 ? unique_counts : shared_counts
        for group in groups
            counts[group] += 1
        end
        if length(groups) == 1
            common_groups = get(
                common_peptide_to_groups,
                peptide_key,
                G[]
            )
            groups[1] in common_groups && (common_unique_counts[groups[1]] += 1)
        end
    end

    opportunities = Dict(
        group => ProteinPeptideOpportunityCounts(
            unique_counts[group],
            shared_counts[group],
            common_unique_counts[group]
        )
        for group in final_groups
    )
    return opportunities
end

function count_protein_peptide_opportunities(
    precursors::LibraryPrecursors,
    final_groups::Set{PGKey},
    protein_group_members::Vector{Vector{UInt32}}
)
    missed_cleavages = getMissedCleavages(precursors)
    num_enzymatic_termini = getNumEnzymaticTermini(precursors)
    num_variable_modifications = getNumVariableModifications(precursors)
    structural_mods = getStructuralMods(precursors)
    common_precursor_mask = BitVector(undef, length(precursors))
    @inbounds for precursor_idx in eachindex(common_precursor_mask)
        common_precursor_mask[precursor_idx] = _is_common_peptide(
            missed_cleavages[precursor_idx],
            num_enzymatic_termini[precursor_idx],
            _num_variable_modifications_at(
                num_variable_modifications,
                structural_mods,
                precursor_idx
            )
        )
    end
    # Integer keys: library accession sets and sequences by ID, groups by their member accessions.
    text_ids = getPrecursorTextIds(precursors)
    return _count_protein_peptide_opportunities(
        text_ids.accession_set_id,
        text_ids.sequence_id,
        getIsDecoy(precursors),
        getEntrapmentGroupId(precursors),
        final_groups,
        set_id -> text_ids.accession_set_members[set_id],
        group -> protein_group_members[group.name];
        common_precursor_mask = common_precursor_mask
    )
end

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add the per-precursor library columns protein inference needs. Text (sequence, accessions,
species, modifications) is not added: it is looked up from the library by `:precursor_idx`.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"

    op = function(df)
        all_is_decoys = getIsDecoy(precursors)::AbstractVector{Bool}
        all_entrap_ids = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}
        all_base_pep_ids = getBasePepId(precursors)::AbstractVector{UInt32}
        all_structural_mods = getStructuralMods(precursors)::AbstractVector{Union{Missing, String}}
        all_num_variable_modifications = getNumVariableModifications(precursors)

        precursor_idx = df.precursor_idx::AbstractVector{UInt32}
        n_rows = length(precursor_idx)

        is_decoy_vec = Vector{Bool}(undef, n_rows)
        for i in 1:n_rows
            is_decoy_vec[i] = all_is_decoys[precursor_idx[i]]
        end
        df.is_decoy = is_decoy_vec

        entrap_vec = Vector{UInt8}(undef, n_rows)
        for i in 1:n_rows
            entrap_vec[i] = all_entrap_ids[precursor_idx[i]]
        end
        df.entrap_id = entrap_vec

        base_pep_ids = Vector{UInt32}(undef, n_rows)
        for i in 1:n_rows
            base_pep_ids[i] = all_base_pep_ids[precursor_idx[i]]
        end
        df.base_pep_id = base_pep_ids

        num_variable_modifications = Vector{UInt8}(undef, n_rows)
        for i in 1:n_rows
            num_variable_modifications[i] = _num_variable_modifications_at(
                all_num_variable_modifications,
                all_structural_mods,
                precursor_idx[i]
            )
        end
        df.num_variable_modifications = num_variable_modifications

        return df
    end

    return desc => op
end

"""
    IntegerInferenceInputs(precursors)

Library lookups for running protein inference on integers: each precursor's sequence ID, and its
accession set as a `ProteinGroupRegistry` group of accession IDs.
"""
struct IntegerInferenceInputs{D, E}
    sequence_id::Vector{UInt32}
    accession_set_id::Vector{UInt32}
    set_group::Vector{UInt32}              # accession set ID -> registry group ID
    is_decoy::D
    entrap_id::E
    registry::ProteinGroupRegistry
    accession_names::Vector{String}
end

function IntegerInferenceInputs(precursors::LibraryPrecursors)
    text_ids = getPrecursorTextIds(precursors)
    registry = ProteinGroupRegistry(length(text_ids.accession_names))
    set_group = UInt32[intern_group!(registry, members) for members in text_ids.accession_set_members]
    return IntegerInferenceInputs(text_ids.sequence_id, text_ids.accession_set_id, set_group,
        getIsDecoy(precursors), getEntrapmentGroupId(precursors), registry, text_ids.accession_names)
end

# (sequence ID, accession set ID, is_decoy, entrap_id) of a precursor: what inference reads of it.
@inline _inference_key(inputs::IntegerInferenceInputs, pid::Integer) = (
    inputs.sequence_id[pid], inputs.accession_set_id[pid], Bool(inputs.is_decoy[pid]), UInt8(inputs.entrap_id[pid]))

"""Infer protein groups for a set of `_inference_key` tuples. Group names are registry IDs."""
function _infer_integer_proteins(inputs::IntegerInferenceInputs, keys)
    proteins = ProteinKey{UInt32}[ProteinKey(inputs.set_group[acc], !dec, ent) for (_, acc, dec, ent) in keys]
    peptides = PeptideKey{UInt32}[PeptideKey(seq, !dec, ent) for (seq, _, dec, ent) in keys]
    isempty(proteins) &&
        return InferenceResult(Dictionary{PeptideKey{UInt32}, ProteinKey{UInt32}}())
    return infer_proteins(proteins, peptides; names = inputs.registry)
end

"""
    _protein_group_ids(results, inputs) -> (results, pg_names, pg_members)

Replace registry group IDs in inference results by `pg_id`s: the rank of the group's name
(its sorted `;`-joined accessions) among all final group names, keyed by name alone. Sorting
`pg_id`s therefore sorts names, so every order downstream matches a sort by name. Returns the
converted results, `pg_names[pg_id]` and `pg_members[pg_id]` (sorted accession IDs).
"""
function _protein_group_ids(results::AbstractVector{<:InferenceResult}, inputs::IntegerInferenceInputs)
    group_ids = Set{UInt32}()
    for result in results
        foreach(protein -> push!(group_ids, protein.name), values(result.peptide_to_protein))
        for candidates in values(result.ambiguous_peptide_to_proteins)
            foreach(protein -> push!(group_ids, protein.name), candidates)
        end
    end
    ids = collect(group_ids)
    names = [join(view(inputs.accession_names, inputs.registry.members[id]), ';') for id in ids]
    pg_names = sort!(unique(names))
    rank = Dict(name => UInt32(i) for (i, name) in enumerate(pg_names))
    to_pg = Dict(id => rank[name] for (id, name) in zip(ids, names))
    pg_members = Vector{Vector{UInt32}}(undef, length(pg_names))
    for (id, name) in zip(ids, names)
        pg_members[rank[name]] = inputs.registry.members[id]
    end
    remap(protein) = PGKey(to_pg[protein.name], protein.is_target, protein.entrap_id)
    converted = map(results) do result
        assigned = Dictionary{PeptideKey{UInt32}, PGKey}()
        for (peptide, protein) in pairs(result.peptide_to_protein)
            insert!(assigned, peptide, remap(protein))
        end
        ambiguous = Dictionary{PeptideKey{UInt32}, Vector{PGKey}}()
        for (peptide, candidates) in pairs(result.ambiguous_peptide_to_proteins)
            insert!(ambiguous, peptide, sort!(map(remap, candidates)))
        end
        InferenceResult(assigned, ambiguous)
    end
    return converted, pg_names, pg_members
end

"""
    add_inferred_protein_column(inference_result::InferenceResult, library_sequences)

Add inferred protein-group assignments to PSMs.
"""
function add_inferred_protein_column(inference_result::InferenceResult, library_sequences::AbstractVector)
    desc = "add_inferred_protein_column"

    op = function(df)
        sequences = view(library_sequences, df.precursor_idx::AbstractVector{UInt32})
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        inferred_proteins = Vector{Union{Missing, UInt32}}(undef, length(sequences))
        for i in eachindex(sequences, is_decoy, entrap_ids)
            pep_key = PeptideKey(sequences[i], !is_decoy[i], entrap_ids[i])
            if haskey(inference_result.peptide_to_protein, pep_key)
                inferred_proteins[i] = inference_result.peptide_to_protein[pep_key].name
            else
                inferred_proteins[i] = missing
            end
        end

        df.inferred_protein_group = inferred_proteins
        return df
    end

    return desc => op
end

"""
    add_quantification_flag(inference_result::InferenceResult, library_sequences)

Mark peptides assigned to inferred protein groups as usable for protein quant/scoring.
"""
function add_quantification_flag(inference_result::InferenceResult, library_sequences::AbstractVector)
    desc = "add_quantification_flag"

    op = function(df)
        sequences = view(library_sequences, df.precursor_idx::AbstractVector{UInt32})
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        use_for_quant = Vector{Bool}(undef, length(sequences))
        for i in eachindex(sequences, is_decoy, entrap_ids)
            pep_key = PeptideKey(sequences[i], !is_decoy[i], entrap_ids[i])
            use_for_quant[i] = haskey(inference_result.peptide_to_protein, pep_key)
        end

        df.use_for_protein_quant = use_for_quant
        return df
    end

    return desc => op
end

"""
    add_protein_ambiguity_id(peptide_to_id, library_sequences)

Annotate PSMs with the normalized ambiguous-peptide assignment ID. Zero denotes a peptide that
is not ambiguous between multiple retained protein groups.
"""
function add_protein_ambiguity_id(peptide_to_id::Dictionary{<:PeptideKey, UInt32}, library_sequences::AbstractVector)
    desc = "add_protein_ambiguity_id"

    op = function(df)
        sequences = view(library_sequences, df.precursor_idx::AbstractVector{UInt32})
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        ambiguity_ids = zeros(UInt32, length(sequences))
        for i in eachindex(sequences, is_decoy, entrap_ids)
            peptide_key = PeptideKey(sequences[i], !is_decoy[i], entrap_ids[i])
            ambiguity_ids[i] = get(peptide_to_id, peptide_key, zero(UInt32))
        end

        df.protein_ambiguity_id = ambiguity_ids
        return df
    end

    return desc => op
end

"""
    run_protein_inference!(search_context; passing_refs, global_inference=true)

Annotate passing precursor tables in place with inferred protein groups and
protein-quant eligibility flags.

`global_inference=true` (default) runs `infer_proteins` once over the union of
unique `(sequence, accession_numbers, is_decoy, entrap_id)` tuples from every
file, then applies the single result to every file. This produces a stable
peptide → group mapping across the experiment and pools decoys for the
downstream protein-level PEP fit.

`global_inference=false` runs inference per file (legacy behavior).

Inference runs on integers (`IntegerInferenceInputs`), and protein groups are written as
`pg_id`s (see `_protein_group_ids`). Returns the in-memory ambiguity mapping, theoretical
unique/shared peptide opportunity counts for the retained final protein groups, and the group
names by `pg_id`.
"""
function run_protein_inference!(
    search_context::SearchContext;
    passing_refs::Vector{PSMFileReference},
    global_inference::Bool = true,
)
    protein_ambiguity_candidates = Dict{UInt32, Vector{PGKey}}()
    final_protein_groups = Set{PGKey}()
    if isempty(passing_refs)
        return (
            protein_ambiguity_candidates = protein_ambiguity_candidates,
            protein_peptide_opportunities =
                Dict{PGKey, ProteinPeptideOpportunityCounts}(),
            protein_group_names = String[]
        )
    end

    precursors = getPrecursors(getSpecLib(search_context))
    inputs = IntegerInferenceInputs(precursors)

    if !global_inference
        annotation_pipeline = TransformPipeline() |>
            add_peptide_metadata(precursors)
        indexed_refs = collect(enumerate(passing_refs))
        @debug_l1 "Annotating passing PSM files with inferred protein groups and protein-quant flags (per-file)"

        # Infer every file first: pg_ids rank group names across all files.
        file_refs = PSMFileReference[]
        file_results = InferenceResult{ProteinKey{UInt32}, PeptideKey{UInt32}}[]
        for (_, psm_ref) in ProgressBar(indexed_refs)
            exists(psm_ref) || continue
            apply_pipeline!(psm_ref, annotation_pipeline)
            pidx = materialize_columns(psm_ref, [:precursor_idx])[!, :precursor_idx]::AbstractVector{UInt32}
            keys = unique(_inference_key(inputs, p) for p in pidx)
            push!(file_refs, psm_ref)
            push!(file_results, _infer_integer_proteins(inputs, keys))
        end
        pg_results, pg_names, pg_members = _protein_group_ids(file_results, inputs)

        last_ambiguity_id = zero(UInt32)
        for (psm_ref, inference_result) in zip(file_refs, pg_results)
            _collect_final_protein_groups!(final_protein_groups, inference_result)
            peptide_to_ambiguity_id, last_ambiguity_id = _register_protein_ambiguities!(
                protein_ambiguity_candidates,
                inference_result,
                last_ambiguity_id
            )

            update_pipeline = TransformPipeline() |>
                add_inferred_protein_column(inference_result, inputs.sequence_id) |>
                add_quantification_flag(inference_result, inputs.sequence_id) |>
                add_protein_ambiguity_id(peptide_to_ambiguity_id, inputs.sequence_id)
            apply_pipeline!(psm_ref, update_pipeline)
        end

        return (
            protein_ambiguity_candidates = protein_ambiguity_candidates,
            protein_peptide_opportunities = count_protein_peptide_opportunities(
                precursors,
                final_protein_groups,
                pg_members
            ),
            protein_group_names = pg_names
        )
    end

    @debug_l1 "Annotating passing PSM files with inferred protein groups and protein-quant flags (global)"

    # Pass 1: stream :precursor_idx from each passing PSM Arrow file and collect the distinct
    # (sequence ID, accession set ID, is_decoy, entrap_id) tuples from the library. No DataFrame
    # is materialized and no file is rewritten — peak per-file memory is one mmap'd UInt32 column.
    unique_set = Set{Tuple{UInt32, UInt32, Bool, UInt8}}()
    for psm_ref in ProgressBar(passing_refs)
        exists(psm_ref) || continue
        table = Arrow.Table(file_path(psm_ref))
        :precursor_idx in Tables.columnnames(table) || continue
        pidx = Tables.getcolumn(table, :precursor_idx)
        @inbounds for i in eachindex(pidx)
            push!(unique_set, _inference_key(inputs, pidx[i]))
        end
    end

    # Pass 2: one InferenceResult from the global tuple set, with groups as pg_ids.
    (inference_result,), pg_names, pg_members =
        _protein_group_ids([_infer_integer_proteins(inputs, unique_set)], inputs)
    _collect_final_protein_groups!(final_protein_groups, inference_result)
    peptide_to_ambiguity_id, _ = _register_protein_ambiguities!(
        protein_ambiguity_candidates,
        inference_result,
        zero(UInt32)
    )

    # Pass 3: compute peptide metadata + inferred protein group +
    # use_for_protein_quant + protein_ambiguity_id directly into row-aligned arrays and write them
    # as a single sidecar instead of rewriting the main file. The sidecar is
    # consolidated into the main file later (at MaxLFQ's sort).
    all_base_pep_ids = getBasePepId(precursors)::AbstractVector{UInt32}

    for psm_ref in ProgressBar(passing_refs)
        exists(psm_ref) || continue
        # Library text stays in the library (looked up by :precursor_idx downstream). Add:
        #   :entrap_id, :base_pep_id  (from library lookup)
        #   :inferred_protein_group (a pg_id), :use_for_protein_quant,
        #   :protein_ambiguity_id  (from inference)
        pidx = materialize_columns(psm_ref, [:precursor_idx])[!, :precursor_idx]::AbstractVector{UInt32}
        n = length(pidx)
        entrap_ids    = Vector{UInt8}(undef, n)
        base_pep_ids  = Vector{UInt32}(undef, n)
        inferred      = Vector{Union{Missing, UInt32}}(undef, n)
        use_for_quant = Vector{Bool}(undef, n)
        ambiguity_ids = zeros(UInt32, n)

        @inbounds for i in 1:n
            p = pidx[i]
            seq, _, dec, ent = _inference_key(inputs, p)
            entrap_ids[i]   = ent
            base_pep_ids[i] = all_base_pep_ids[p]
            pep_key = PeptideKey(seq, !dec, ent)
            if haskey(inference_result.peptide_to_protein, pep_key)
                inferred[i]       = inference_result.peptide_to_protein[pep_key].name
                use_for_quant[i]  = true
            else
                inferred[i]       = missing
                use_for_quant[i]  = false
                ambiguity_ids[i]  = get(
                    peptide_to_ambiguity_id,
                    pep_key,
                    zero(UInt32)
                )
            end
        end

        # :is_decoy is deliberately NOT emitted: it equals `decoy`, itself exactly `!target`, so it
        # was the third output column encoding one bit (0 mismatches across 235,194 rows). Nothing in
        # the search path read it; the only `is_decoy` readers are in BuildSpecLib, a different table.
        add_columns_via_sidecar!(psm_ref,
            :entrap_id              => entrap_ids,
            :base_pep_id            => base_pep_ids,
            :inferred_protein_group => inferred,
            :use_for_protein_quant  => use_for_quant,
            :protein_ambiguity_id   => ambiguity_ids;
            tag = "ProteinInference")
    end
    return (
        protein_ambiguity_candidates = protein_ambiguity_candidates,
        protein_peptide_opportunities = count_protein_peptide_opportunities(
            precursors,
            final_protein_groups,
            pg_members
        ),
        protein_group_names = pg_names
    )
end
