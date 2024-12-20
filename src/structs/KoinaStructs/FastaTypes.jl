"""
Represents an entry in a FASTA file with associated metadata.
"""
struct FastaEntry
    id::String
    description::String 
    proteome::String
    sequence::String
    base_pep_id::UInt32
    entrapment_group_id::UInt8
    is_decoy::Bool
end

# Accessor methods
get_id(entry::FastaEntry) = entry.id
get_description(entry::FastaEntry) = entry.description
get_proteome(entry::FastaEntry) = entry.proteome
get_sequence(entry::FastaEntry) = entry.sequence
get_base_pep_id(entry::FastaEntry) = entry.base_pep_id
get_entrapment_group_id(entry::FastaEntry) = entry.entrapment_group_id
is_decoy(entry::FastaEntry) = entry.is_decoy
