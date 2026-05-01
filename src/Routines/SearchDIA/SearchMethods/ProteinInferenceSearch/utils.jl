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
        all_pair_ids = getPairIdx(precursors)

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

        pair_ids = Vector{UInt32}(undef, n_rows)
        for i in 1:n_rows
            pair_ids[i] = extract_pair_idx(all_pair_ids, precursor_idx[i])
        end
        df.pair_id = pair_ids

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
    run_protein_inference!(search_context; passing_refs, global_inference=true)

Annotate passing precursor tables in place with inferred protein groups and
protein-quant eligibility flags.

`global_inference=true` (default) runs `infer_proteins` once over the union of
unique `(sequence, accession_numbers, is_decoy, entrap_id)` tuples from every
file, then applies the single result to every file. This produces a stable
peptide → group mapping across the experiment and pools decoys for the
downstream protein-level PEP fit.

`global_inference=false` runs inference per file (legacy behavior).
"""
function run_protein_inference!(
    search_context::SearchContext;
    passing_refs::Vector{PSMFileReference},
    global_inference::Bool = true,
)
    isempty(passing_refs) && return nothing

    precursors = getPrecursors(getSpecLib(search_context))
    annotation_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors)

    if !global_inference
        indexed_refs = collect(enumerate(passing_refs))
        @user_info "Annotating passing PSM files with inferred protein groups and protein-quant flags (per-file)"

        for (_, psm_ref) in ProgressBar(indexed_refs)
            exists(psm_ref) || continue

            apply_pipeline!(psm_ref, annotation_pipeline)
            prepared_df = load_dataframe(psm_ref)
            inference_result = apply_inference_to_dataframe(prepared_df, precursors)

            update_pipeline = TransformPipeline() |>
                add_inferred_protein_column(inference_result) |>
                add_quantification_flag(inference_result)
            apply_pipeline!(psm_ref, update_pipeline)
        end

        return nothing
    end

    @user_info "Annotating passing PSM files with inferred protein groups and protein-quant flags (global)"

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

    # Pass 3: annotate each file with peptide metadata and apply the single
    # global inference result. Combined into one file rewrite per file.
    final_pipeline = annotation_pipeline |>
        add_inferred_protein_column(inference_result) |>
        add_quantification_flag(inference_result)
    for psm_ref in ProgressBar(passing_refs)
        exists(psm_ref) || continue
        apply_pipeline!(psm_ref, final_pipeline)
    end

    return nothing
end
