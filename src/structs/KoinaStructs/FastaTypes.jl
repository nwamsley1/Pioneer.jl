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
Represents an entry in a FASTA file with associated metadata.
"""
struct FastaEntry
    accession::String
    description::String
    gene::String
    protein::String
    organism::String
    proteome::String
    sequence::String
    start_idx::UInt32
    structural_mods::Union{Missing,Vector{PeptideMod}}
    isotopic_mods::Union{Missing,Vector{PeptideMod}}
    charge::UInt8
    base_target_id::UInt32
    base_pep_id::UInt32
    entrapment_pair_id::UInt32
    is_decoy::Bool
end

# Accessor methods
get_id(entry::FastaEntry) = entry.accession
get_description(entry::FastaEntry) = entry.description
get_gene(entry::FastaEntry) = entry.gene
get_protein(entry::FastaEntry) = entry.protein
get_organism(entry::FastaEntry) = entry.organism
get_proteome(entry::FastaEntry) = entry.proteome
get_sequence(entry::FastaEntry) = entry.sequence
get_start_idx(entry::FastaEntry) = entry.start_idx
get_base_target_id(entry::FastaEntry) = entry.base_target_id
get_base_pep_id(entry::FastaEntry) = entry.base_pep_id
get_entrapment_pair_id(entry::FastaEntry) = entry.entrapment_pair_id
is_decoy(entry::FastaEntry) = entry.is_decoy
get_isotopic_mods(entry::FastaEntry) = entry.isotopic_mods
get_structural_mods(entry::FastaEntry) = entry.structural_mods
get_charge(entry::FastaEntry) = entry.charge

# ---------------------------------------------------------------------------
# Backward-compatible outer constructors
# ---------------------------------------------------------------------------

"""
    FastaEntry(acc, desc, gene, protein, organism, proteome, sequence,
               start_idx::UInt32, struct_mods, iso_mods,
               charge::UInt8, base_target_id::UInt32, base_pep_id::UInt32,
               entrapment_pair_id::Union{UInt8,UInt32}, is_decoy::Bool)

Backward-compatible convenience constructor (struct has no base_prec_id field).
Coerces `entrapment_pair_id` to UInt32.
"""
function FastaEntry(
    accession::AbstractString,
    description::AbstractString,
    gene::AbstractString,
    protein::AbstractString,
    organism::AbstractString,
    proteome::AbstractString,
    sequence::AbstractString,
    start_idx::UInt32,
    struct_mods::Union{Missing,Vector{PeptideMod}},
    iso_mods::Union{Missing,Vector{PeptideMod}},
    charge::UInt8,
    base_target_id::UInt32,
    base_pep_id::UInt32,
    entrapment_pair_id::Integer,
    is_decoy::Bool,
)
    return FastaEntry(
        String(accession), String(description), String(gene), String(protein), String(organism), String(proteome), String(sequence),
        start_idx, struct_mods, iso_mods, charge,
        base_target_id, base_pep_id, UInt32(entrapment_pair_id), is_decoy,
    )
end

"""
    FastaEntry(acc, desc, gene, protein, organism, proteome, sequence,
               start_idx::UInt32, struct_mods, iso_mods,
               charge::UInt8, base_target_id::UInt32, base_pep_id::UInt32,
               base_prec_id::UInt32,
               entrapment_pair_id::Union{UInt8,UInt32}, is_decoy::Bool)

Compatibility constructor that accepts a `base_prec_id` argument but ignores it.
"""
function FastaEntry(
    accession::AbstractString,
    description::AbstractString,
    gene::AbstractString,
    protein::AbstractString,
    organism::AbstractString,
    proteome::AbstractString,
    sequence::AbstractString,
    start_idx::UInt32,
    struct_mods::Union{Missing,Vector{PeptideMod}},
    iso_mods::Union{Missing,Vector{PeptideMod}},
    charge::UInt8,
    base_target_id::UInt32,
    base_pep_id::UInt32,
    base_prec_id::UInt32,
    entrapment_pair_id::Integer,
    is_decoy::Bool,
)
    return FastaEntry(
        String(accession), String(description), String(gene), String(protein), String(organism), String(proteome), String(sequence),
        start_idx, struct_mods, iso_mods, charge,
        base_target_id, base_pep_id, UInt32(entrapment_pair_id), is_decoy,
    )
end
