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
Data types for protein inference and scoring in the ScoringSearch module.

These types provide clear, type-safe structures for protein group analysis,
replacing complex NamedTuple-based dictionaries.
"""

using Dictionaries

"""
    ProteinKey

Unique identifier for a protein or protein group.

# Fields
- `name::String`: Protein name or semicolon-delimited list of indistinguishable proteins
- `is_target::Bool`: True if this is a target protein (not decoy)
- `entrap_id::UInt8`: Entrapment group identifier for FDR estimation
"""
struct ProteinKey
    name::String
    is_target::Bool
    entrap_id::UInt8
end

"""
    PeptideKey

Unique identifier for a peptide in the context of protein inference.

# Fields
- `sequence::String`: Peptide amino acid sequence
- `is_target::Bool`: True if from target proteins (not decoy)
- `entrap_id::UInt8`: Entrapment group identifier
"""
struct PeptideKey
    sequence::String
    is_target::Bool
    entrap_id::UInt8
end

"""
    ProteinFeatures

Statistical features of a protein group used for scoring and FDR control.

# Fields
- `n_peptides::UInt16`: Number of unique peptides observed
- `n_possible_peptides::UInt16`: Total possible peptides from in silico digestion
- `peptide_coverage::Float32`: Fraction of possible peptides observed
- `total_peptide_length::UInt16`: Sum of lengths of all observed peptides
- `log_n_possible_peptides::Float32`: Log of n_possible_peptides for regression
- `log_binom_coeff::Float32`: Log binomial coefficient for peptide sampling
"""
struct ProteinFeatures
    n_peptides::UInt16
    n_possible_peptides::UInt16
    peptide_coverage::Float32
    total_peptide_length::UInt16
    log_n_possible_peptides::Float32
    log_binom_coeff::Float32
end

"""
    ProteinGroup

Complete protein group with scoring and peptide information.

# Fields
- `key::ProteinKey`: Unique identifier
- `score::Float32`: Protein group score (higher is better)
- `peptides::Set{String}`: Set of unique peptide sequences
- `features::ProteinFeatures`: Statistical features for FDR control
"""
struct ProteinGroup
    key::ProteinKey
    score::Float32
    peptides::Set{String}
    features::ProteinFeatures
end

"""
    InferenceResult

Result of protein inference for a single file.

# Fields
- `peptide_to_protein::Dictionary{PeptideKey, ProteinKey}`: Maps peptides to inferred proteins
- `use_for_quant::Dictionary{PeptideKey, Bool}`: Whether peptide should be used for quantification
"""
struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
    use_for_quant::Dictionary{PeptideKey, Bool}
end

"""
    ProteinGroupBuilder

Mutable builder for constructing protein groups incrementally.

# Fields
- `key::ProteinKey`: Protein identifier
- `score::Float32`: Accumulated score (in log space during construction)
- `peptides::Set{String}`: Unique peptides
"""
mutable struct ProteinGroupBuilder
    key::ProteinKey
    score::Float32
    peptides::Set{String}
    
    function ProteinGroupBuilder(key::ProteinKey)
        new(key, 0.0f0, Set{String}())
    end
end

"""
    add_peptide!(builder::ProteinGroupBuilder, peptide::String, score::Float32)

Add a peptide to the protein group with its score.

# Arguments
- `builder`: ProteinGroupBuilder to update
- `peptide`: Peptide sequence
- `score`: Peptide probability score (0-1)
"""
function add_peptide!(builder::ProteinGroupBuilder, peptide::String, score::Float32)
    push!(builder.peptides, peptide)
    builder.score += log1p(-score)  # Accumulate in log space
end

"""
    finalize(builder::ProteinGroupBuilder, features::ProteinFeatures) -> ProteinGroup

Convert a builder to a final ProteinGroup with calculated features.

# Arguments
- `builder`: Completed ProteinGroupBuilder
- `features`: Calculated protein features

# Returns
- `ProteinGroup`: Finalized protein group with score converted from log space
"""
function finalize(builder::ProteinGroupBuilder, features::ProteinFeatures)::ProteinGroup
    return ProteinGroup(
        builder.key,
        -builder.score,  # Convert from log space
        builder.peptides,
        features
    )
end

"""
    FileMapping

Bidirectional mapping between PSM and protein group files.

# Fields
- `psm_to_pg::Dictionary{String, String}`: Maps PSM file paths to protein group paths
- `pg_to_psm::Dictionary{String, String}`: Maps protein group paths to PSM file paths
"""
mutable struct FileMapping
    psm_to_pg::Dictionary{String, String}
    pg_to_psm::Dictionary{String, String}
    
    FileMapping() = new(Dictionary{String, String}(), Dictionary{String, String}())
end

"""
    add_mapping!(mapping::FileMapping, psm_path::String, pg_path::String)

Add a bidirectional file mapping.
"""
function add_mapping!(mapping::FileMapping, psm_path::String, pg_path::String)
    insert!(mapping.psm_to_pg, psm_path, pg_path)
    insert!(mapping.pg_to_psm, pg_path, psm_path)
end

# Conversion utilities for backward compatibility

"""
    to_protein_key(nt::NamedTuple) -> ProteinKey

Convert old NamedTuple format to ProteinKey.
"""
function to_protein_key(nt::NamedTuple{(:protein_name, :target, :entrap_id), Tuple{String, Bool, UInt8}})
    return ProteinKey(nt.protein_name, nt.target, nt.entrap_id)
end

"""
    to_namedtuple(key::ProteinKey) -> NamedTuple

Convert ProteinKey to old NamedTuple format for compatibility.
"""
function to_namedtuple(key::ProteinKey)
    return (protein_name = key.name, target = key.is_target, entrap_id = key.entrap_id)
end

"""
    to_peptide_key(nt::NamedTuple) -> PeptideKey

Convert old NamedTuple format to PeptideKey.
"""
function to_peptide_key(nt::NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}})
    return PeptideKey(nt.peptide, !nt.decoy, nt.entrap_id)
end