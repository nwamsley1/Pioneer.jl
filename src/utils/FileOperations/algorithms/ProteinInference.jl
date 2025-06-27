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
Protein inference algorithm wrappers.

Provides safe wrappers around protein inference algorithms
with proper file reference management and validation.
"""

using Arrow, DataFrames, Tables

#==========================================================
Algorithm Wrappers
==========================================================#

"""
    apply_protein_inference(psm_ref::PSMFileReference, output_path::String, 
                          precursors::Union{Nothing, DataFrame}, 
                          protein_inference_fn::Function;
                          min_peptides=2, batch_size=1_000_000)

Apply protein inference algorithm to PSMs, creating protein groups.
Uses the provided protein_inference_fn function (typically getProteinGroupsDict from ScoringSearch/utils.jl).
"""
function apply_protein_inference(psm_ref::PSMFileReference, 
                               output_path::String,
                               precursors,
                               protein_inference_fn::Function;
                               min_peptides::Int=2)
    validate_exists(psm_ref)
    
    # Validate required columns for protein inference
    required_cols = Set([:prob, :inferred_protein_group, :precursor_idx, :target, :entrapment_group_id])
    validate_required_columns(schema(psm_ref), required_cols)
    
    # For protein inference, we need all PSMs at once (not streaming)
    # This is because protein inference needs global view of peptide-protein relationships
    psms_df = DataFrame(Arrow.Table(file_path(psm_ref)))
    
    # Call provided protein inference function
    protein_groups_dict = protein_inference_fn(
        psms_df.prob,
        psms_df.inferred_protein_group,
        psms_df.precursor_idx,
        psms_df.target,
        psms_df.entrapment_group_id,
        precursors,
        min_peptides
    )
    
    # Convert to DataFrame format expected by ScoringSearch
    # Following the schema from write_protein_groups_arrow
    n_groups = length(protein_groups_dict)
    
    # Pre-allocate vectors for each column
    protein_names = Vector{String}(undef, n_groups)
    targets = Vector{Bool}(undef, n_groups)
    entrap_ids = Vector{UInt8}(undef, n_groups) 
    pg_scores = Vector{Float32}(undef, n_groups)
    n_peptides = Vector{Int64}(undef, n_groups)
    
    # Fill vectors from protein groups
    for (i, ((name, target, entrap_id), peptide_set)) in enumerate(protein_groups_dict)
        protein_names[i] = name
        targets[i] = target
        entrap_ids[i] = entrap_id
        # Basic score as negative log-sum of peptide probabilities
        pg_scores[i] = -sum(log(1.0f0 - psms_df.prob[pid]) for pid in peptide_set)
        n_peptides[i] = length(peptide_set)
    end
    
    # Create DataFrame
    protein_df = DataFrame(
        protein_name = protein_names,
        target = targets,
        entrapment_group_id = entrap_ids,
        pg_score = pg_scores,
        n_peptides = n_peptides
    )
    
    # Sort by score (descending)
    sort!(protein_df, :pg_score, rev=true)
    
    # Write to Arrow using writeArrow for Windows compatibility
    writeArrow(output_path, protein_df)
    
    return ProteinGroupFileReference(output_path)
end

"""
    update_psms_with_scores(psm_ref::PSMFileReference, protein_ref::ProteinGroupFileReference,
                          output_path::String; batch_size=100_000)

Update PSMs with protein group scores from protein reference.
Streaming operation that preserves memory efficiency.
"""
function update_psms_with_scores(psm_ref::PSMFileReference, 
                               protein_ref::ProteinGroupFileReference,
                               output_path::String;
                               batch_size::Int=100_000)
    validate_exists(psm_ref)
    validate_exists(protein_ref)
    
    # Load protein scores (typically smaller)
    protein_df = DataFrame(Arrow.Table(file_path(protein_ref)))
    
    # Create lookup dictionary
    protein_scores = Dict{Tuple{String,Bool,UInt8},NamedTuple}()
    for row in eachrow(protein_df)
        key = (row.protein_name, row.target, row.entrapment_group_id)
        protein_scores[key] = (
            pg_score = row.pg_score,
            global_pg_score = row.global_pg_score,
            pg_qval = row.pg_qval,
            global_pg_qval = row.global_pg_qval
        )
    end
    
    # Stream through PSMs, updating scores
    stream_transform(psm_ref, output_path, batch_size=batch_size) do batch
        # Add score columns
        batch.pg_score = Vector{Float32}(undef, nrow(batch))
        batch.global_pg_score = Vector{Float32}(undef, nrow(batch))
        batch.pg_qval = Vector{Float32}(undef, nrow(batch))
        batch.global_qval_pg = Vector{Float32}(undef, nrow(batch))
        
        # Update scores
        for i in 1:nrow(batch)
            key = (batch.inferred_protein_group[i], 
                   batch.target[i], 
                   batch.entrapment_group_id[i])
            
            if haskey(protein_scores, key)
                scores = protein_scores[key]
                batch.pg_score[i] = scores.pg_score
                batch.global_pg_score[i] = scores.global_pg_score
                batch.pg_qval[i] = scores.pg_qval
                batch.global_qval_pg[i] = scores.global_pg_qval
            else
                # This should not happen if protein inference was done correctly
                error("Missing protein scores for: $key")
            end
        end
        
        return batch
    end
end

