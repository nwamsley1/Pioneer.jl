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
    base_pep_id::UInt32
    base_prec_id::UInt32
    entrapment_group_id::UInt8
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
get_base_pep_id(entry::FastaEntry) = entry.base_pep_id
get_base_prec_id(entry::FastaEntry) = entry.base_prec_id
get_entrapment_group_id(entry::FastaEntry) = entry.entrapment_group_id
is_decoy(entry::FastaEntry) = entry.is_decoy
get_isotopic_mods(entry::FastaEntry) = entry.isotopic_mods
get_structural_mods(entry::FastaEntry) = entry.structral_mods
get_charge(entry::FastaEntry) = entry.charge