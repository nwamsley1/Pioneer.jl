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
    structral_mods::Union{Missing,Vector{PeptideMod}}
    isotopic_mods::Union{Missing,Vector{PeptideMod}}
    charge::UInt8
    base_seq_id::UInt32
    base_pep_id::UInt32
    base_prec_id::UInt32
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
get_base_seq_id(entry::FastaEntry) = entry.base_seq_id
get_base_pep_id(entry::FastaEntry) = entry.base_pep_id
get_base_prec_id(entry::FastaEntry) = entry.base_prec_id
get_entrapment_pair_id(entry::FastaEntry) = entry.entrapment_pair_id
is_decoy(entry::FastaEntry) = entry.is_decoy
get_isotopic_mods(entry::FastaEntry) = entry.isotopic_mods
get_structural_mods(entry::FastaEntry) = entry.structral_mods
get_charge(entry::FastaEntry) = entry.charge