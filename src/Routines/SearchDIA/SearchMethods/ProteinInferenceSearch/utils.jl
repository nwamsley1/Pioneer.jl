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
    candidates_by_id::Dict{UInt32, Vector{ProteinKey}},
    inference_result::InferenceResult,
    last_id::UInt32
)
    peptide_to_id = Dictionary{PeptideKey, UInt32}()
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
    final_groups::Set{ProteinKey},
    inference_result::InferenceResult
)
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
    final_groups::Set{ProteinKey};
    common_precursor_mask::Union{Nothing, AbstractVector{Bool}} = nothing
)
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

    accession_to_groups =
        Dict{Tuple{String, Bool, UInt8}, Vector{ProteinKey}}()
    for group in final_groups
        for accession in split(group.name, ';')
            key = (strip(accession), group.is_target, group.entrap_id)
            push!(get!(accession_to_groups, key, ProteinKey[]), group)
        end
    end
    for groups in values(accession_to_groups)
        sort!(unique!(groups))
    end

    group_cache =
        Dict{Tuple{String, Bool, UInt8}, Vector{ProteinKey}}()
    peptide_to_groups =
        Dict{Tuple{String, Bool, UInt8}, Vector{ProteinKey}}()
    common_peptide_to_groups =
        Dict{Tuple{String, Bool, UInt8}, Vector{ProteinKey}}()

    @inbounds for i in eachindex(accession_numbers, sequences, is_decoys, entrap_ids)
        target = !is_decoys[i]
        entrap_id = entrap_ids[i]
        accessions = String(accession_numbers[i])
        groups = get!(group_cache, (accessions, target, entrap_id)) do
            mapped_groups = ProteinKey[]
            for accession in split(accessions, ';')
                append!(
                    mapped_groups,
                    get(
                        accession_to_groups,
                        (strip(accession), target, entrap_id),
                        ProteinKey[]
                    )
                )
            end
            sort!(unique!(mapped_groups))
        end

        isempty(groups) && continue
        peptide_key = (String(sequences[i]), target, entrap_id)
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
                ProteinKey[]
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
    final_groups::Set{ProteinKey}
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
    return count_protein_peptide_opportunities(
        getAccessionNumbers(precursors),
        getSequence(precursors),
        getIsDecoy(precursors),
        getEntrapmentGroupId(precursors),
        final_groups;
        common_precursor_mask = common_precursor_mask
    )
end

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add peptide sequence and protein information from the precursor library.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"

    op = function(df)
        all_sequences = getSequence(precursors)::AbstractVector{String}
        all_accessions = getAccessionNumbers(precursors)::AbstractVector{String}
        all_is_decoys = getIsDecoy(precursors)::AbstractVector{Bool}
        all_entrap_ids = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}
        all_species = getProteomeIdentifiers(precursors)::AbstractVector{<:AbstractString}
        all_base_pep_ids = getBasePepId(precursors)::AbstractVector{UInt32}
        all_structural_mods = getStructuralMods(precursors)::AbstractVector{Union{Missing, String}}
        all_isotopic_mods = getIsotopicMods(precursors)::AbstractVector{Union{Missing, String}}
        all_num_variable_modifications = getNumVariableModifications(precursors)

        precursor_idx = df.precursor_idx::AbstractVector{UInt32}
        n_rows = length(precursor_idx)

        sequences = Vector{String}(undef, n_rows)
        for i in 1:n_rows
            sequences[i] = all_sequences[precursor_idx[i]]
        end
        df.sequence = sequences

        accessions = Vector{String}(undef, n_rows)
        for i in 1:n_rows
            accessions[i] = all_accessions[precursor_idx[i]]
        end
        df.accession_numbers = accessions

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

        species = Vector{String}(undef, n_rows)
        for i in 1:n_rows
            species[i] = join(sort(unique(split(coalesce(all_species[precursor_idx[i]], ""), ';'))), ';')
        end
        df.species = species

        base_pep_ids = Vector{UInt32}(undef, n_rows)
        for i in 1:n_rows
            base_pep_ids[i] = all_base_pep_ids[precursor_idx[i]]
        end
        df.base_pep_id = base_pep_ids

        structural_mods = Vector{String}(undef, n_rows)
        for i in 1:n_rows
            structural_mods[i] = coalesce(all_structural_mods[precursor_idx[i]], "")
        end
        df.structural_mods = structural_mods

        isotopic_mods = Vector{String}(undef, n_rows)
        for i in 1:n_rows
            isotopic_mods[i] = coalesce(all_isotopic_mods[precursor_idx[i]], "")
        end
        df.isotopic_mods = isotopic_mods

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
    apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)

Apply the core protein inference algorithm to a prepared DataFrame.
"""
function apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)
    if nrow(df) == 0
        return InferenceResult(Dictionary{PeptideKey, ProteinKey}())
    end

    unique_pairs = unique(df, [:sequence, :accession_numbers, :is_decoy, :entrap_id])
    proteins_vec = Vector{ProteinKey}(undef, nrow(unique_pairs))
    peptides_vec = Vector{PeptideKey}(undef, nrow(unique_pairs))

    for (i, row) in enumerate(eachrow(unique_pairs))
        proteins_vec[i] = ProteinKey(row.accession_numbers, !row.is_decoy, row.entrap_id)
        peptides_vec[i] = PeptideKey(row.sequence, !row.is_decoy, row.entrap_id)
    end

    return infer_proteins(proteins_vec, peptides_vec)
end

"""
    add_inferred_protein_column(inference_result::InferenceResult)

Add inferred protein-group assignments to PSMs.
"""
function add_inferred_protein_column(inference_result::InferenceResult)
    desc = "add_inferred_protein_column"

    op = function(df)
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        inferred_proteins = Vector{Union{Missing, String}}(undef, length(sequences))
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
    add_quantification_flag(inference_result::InferenceResult)

Mark peptides assigned to inferred protein groups as usable for protein quant/scoring.
"""
function add_quantification_flag(inference_result::InferenceResult)
    desc = "add_quantification_flag"

    op = function(df)
        sequences = df.sequence::AbstractVector{String}
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
    add_protein_ambiguity_id(peptide_to_id)

Annotate PSMs with the normalized ambiguous-peptide assignment ID. Zero denotes a peptide that
is not ambiguous between multiple retained protein groups.
"""
function add_protein_ambiguity_id(peptide_to_id::Dictionary{PeptideKey, UInt32})
    desc = "add_protein_ambiguity_id"

    op = function(df)
        sequences = df.sequence::AbstractVector{String}
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

Returns the in-memory ambiguity mapping and theoretical unique/shared peptide
opportunity counts for the retained final protein groups.
"""
function run_protein_inference!(
    search_context::SearchContext;
    passing_refs::Vector{PSMFileReference},
    global_inference::Bool = true,
)
    protein_ambiguity_candidates = Dict{UInt32, Vector{ProteinKey}}()
    final_protein_groups = Set{ProteinKey}()
    if isempty(passing_refs)
        return (
            protein_ambiguity_candidates = protein_ambiguity_candidates,
            protein_peptide_opportunities =
                Dict{ProteinKey, ProteinPeptideOpportunityCounts}()
        )
    end

    precursors = getPrecursors(getSpecLib(search_context))
    annotation_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors)

    if !global_inference
        indexed_refs = collect(enumerate(passing_refs))
        @debug_l1 "Annotating passing PSM files with inferred protein groups and protein-quant flags (per-file)"

        last_ambiguity_id = zero(UInt32)

        for (_, psm_ref) in ProgressBar(indexed_refs)
            exists(psm_ref) || continue

            apply_pipeline!(psm_ref, annotation_pipeline)
            prepared_df = load_dataframe(psm_ref)
            inference_result = apply_inference_to_dataframe(prepared_df, precursors)
            _collect_final_protein_groups!(final_protein_groups, inference_result)
            peptide_to_ambiguity_id, last_ambiguity_id = _register_protein_ambiguities!(
                protein_ambiguity_candidates,
                inference_result,
                last_ambiguity_id
            )

            update_pipeline = TransformPipeline() |>
                add_inferred_protein_column(inference_result) |>
                add_quantification_flag(inference_result) |>
                add_protein_ambiguity_id(peptide_to_ambiguity_id)
            apply_pipeline!(psm_ref, update_pipeline)
        end

        return (
            protein_ambiguity_candidates = protein_ambiguity_candidates,
            protein_peptide_opportunities = count_protein_peptide_opportunities(
                precursors,
                final_protein_groups
            )
        )
    end

    @debug_l1 "Annotating passing PSM files with inferred protein groups and protein-quant flags (global)"

    # Pass 1: stream :precursor_idx from each passing PSM Arrow file and look
    # up the (sequence, accession_numbers, is_decoy, entrap_id) tuple directly
    # from the spec-lib precursor arrays. No DataFrame is materialized and no
    # file is rewritten — peak per-file memory is one mmap'd UInt32 column.
    all_seqs    = getSequence(precursors)::AbstractVector{String}
    all_accs    = getAccessionNumbers(precursors)::AbstractVector{String}
    all_decoys  = getIsDecoy(precursors)::AbstractVector{Bool}
    all_entraps = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}

    UniqueKey = Tuple{String, String, Bool, UInt8}
    unique_set = Set{UniqueKey}()
    for psm_ref in ProgressBar(passing_refs)
        exists(psm_ref) || continue
        table = Arrow.Table(file_path(psm_ref))
        :precursor_idx in Tables.columnnames(table) || continue
        pidx = Tables.getcolumn(table, :precursor_idx)
        @inbounds for i in eachindex(pidx)
            p = pidx[i]
            push!(unique_set, (
                String(all_seqs[p]),
                String(all_accs[p]),
                Bool(all_decoys[p]),
                UInt8(all_entraps[p]),
            ))
        end
    end

    # Pass 2: build a single InferenceResult from the global tuple set.
    n_unique = length(unique_set)
    proteins_vec = Vector{ProteinKey}(undef, n_unique)
    peptides_vec = Vector{PeptideKey}(undef, n_unique)
    let i = 0
        for t in unique_set
            i += 1
            sequence, accession_numbers, is_decoy_val, entrap_val = t
            proteins_vec[i] = ProteinKey(accession_numbers, !is_decoy_val, entrap_val)
            peptides_vec[i] = PeptideKey(sequence, !is_decoy_val, entrap_val)
        end
    end
    inference_result = if n_unique == 0
        InferenceResult(Dictionary{PeptideKey, ProteinKey}())
    else
        infer_proteins(proteins_vec, peptides_vec)
    end
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
        # By the time ProteinInference runs, IntegrateChromatograms has
        # already populated :sequence, :accession_numbers, :species,
        # :structural_mods, :isotopic_mods on the main file. We only need
        # to add the 5 columns it does NOT already provide:
        #   :entrap_id, :base_pep_id  (from library lookup)
        #   :inferred_protein_group, :use_for_protein_quant,
        #   :protein_ambiguity_id  (from inference)
        # The old code redundantly overwrote the IntegrateChromatograms cols.
        pidx = materialize_columns(psm_ref, [:precursor_idx])[!, :precursor_idx]::AbstractVector{UInt32}
        n = length(pidx)
        entrap_ids    = Vector{UInt8}(undef, n)
        base_pep_ids  = Vector{UInt32}(undef, n)
        inferred      = Vector{Union{Missing, String}}(undef, n)
        use_for_quant = Vector{Bool}(undef, n)
        ambiguity_ids = zeros(UInt32, n)

        @inbounds for i in 1:n
            p = pidx[i]
            seq = all_seqs[p]
            dec = all_decoys[p]
            ent = all_entraps[p]
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
            final_protein_groups
        )
    )
end
