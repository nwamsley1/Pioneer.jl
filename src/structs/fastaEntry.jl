struct FastaEntry
    uniprot_id::String
    description::String
    proteome_id::String
    sequence::String
    base_pep_id::UInt32
    entrapment_group_id::UInt8
    decoy::Bool
end

getID(fe::FastaEntry) = fe.uniprot_id
getDescription(fe::FastaEntry) = fe.description
getProteome(fe::FastaEntry) = fe.proteome_id
getSeq(fe::FastaEntry) = fe.sequence
getBasePepId(fe::FastaEntry) = fe.base_pep_id
getEntrapmentGroupId(fe::FastaEntry) = fe.entrapment_group_id
isDecoy(fe::FastaEntry) = fe.decoy