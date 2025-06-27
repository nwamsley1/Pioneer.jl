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
    id::String
    description::String 
    proteome::String
    sequence::String
    structral_mods::Union{Missing,Vector{PeptideMod}}
    isotopic_mods::Union{Missing,Vector{PeptideMod}}
    charge::UInt8
    base_pep_id::UInt32
    base_prec_id::UInt32
    entrapment_group_id::UInt8
    is_decoy::Bool
end

# Accessor methods
get_id(entry::FastaEntry) = entry.id
get_description(entry::FastaEntry) = entry.description
get_proteome(entry::FastaEntry) = entry.proteome
get_sequence(entry::FastaEntry) = entry.sequence
get_base_pep_id(entry::FastaEntry) = entry.base_pep_id
get_base_prec_id(entry::FastaEntry) = entry.base_prec_id
get_entrapment_group_id(entry::FastaEntry) = entry.entrapment_group_id
is_decoy(entry::FastaEntry) = entry.is_decoy
get_isotopic_mods(entry::FastaEntry) = entry.isotopic_mods
get_structural_mods(entry::FastaEntry) = entry.structral_mods
get_charge(entry::FastaEntry) = entry.charge