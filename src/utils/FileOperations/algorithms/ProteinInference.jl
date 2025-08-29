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
        batch.qlobal_pg_qval = Vector{Float32}(undef, nrow(batch))
        
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
                batch.qlobal_pg_qval[i] = scores.global_pg_qval
            else
                # This should not happen if protein inference was done correctly
                error("Missing protein scores for: $key")
            end
        end
        
        return batch
    end
end

