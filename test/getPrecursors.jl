
struct PeptideGroup
    prot_ids::Set{UInt32}
end
function addProt!(pg::PeptideGroup, prot_id::UInt32)
    push!(pg.prot_ids, prot_id)
end

struct Peptide
    sequence::String
    pep_group_id::UInt32
end

getSeq(p::Peptide) = p.sequence
getGroupID(p::Peptide) = p.pep_group_id

struct Protein
    name::String
    pep_group_ids::Set{UInt32}
end
Protein(name::String) = Protein(name, Set{UInt32}())
PeptideGroup() = PeptideGroup("", Set{UInt32}(), UInt32(0))

struct PrecursorTable
    id_prot::UnorderedDictionary{UInt32, Protein}
    prot_id::UnorderedDictionary{String, UInt32}

    id_pep_group::UnorderedDictionary{UInt32, PeptideGroup}
    pep_group_id::UnorderedDictionary{String, UInt32}

    pep_id::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptide group
end

getProteinID(p::PrecursorTable, protein::String) = p.prot_id[protein]
getPepGroupID(p::PrecursorTable, peptide::String) = p.pep_group_id[peptide]
getPepGroup(p::PrecursorTable, pep_id::UInt32) = p.id_pep_group[pep_id]

function addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)
    addProt!(getPepGroup(p, getPepGroupID(p, peptide)), getProteinID(p, protein))
end
precursor_table.id_pep_group[precursor_table.pep_group_id[peptide]]

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, Protein}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, PeptideGroup}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, Peptide}()
                                    )

function addProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)
    insert!(precursor_table.prot_id, protein, prot_id);
    insert!(precursor_table.id_prot, prot_id, Protein(protein));
end

function addPeptideGroup!(peptide_seq::String, pep_group_id::UInt32, protein::String, precursor_table::PrecursorTable)
        insert!(precursor_table.pep_group_id, peptide_seq, pep_group_id);
        insert!(precursor_table.id_pep_group, 
                    pep_group_id, 
                    PeptideGroup(Set(getProteinID(precursor_table, protein)))
                )
end


PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, Peptide}(),
                                    Set{String}(),
                                    Set{String}()
                                    )
function applyMods!(var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String, peptides::Vector{Peptide}, group_id::UInt32; n::Int = 3)
    applyVariableMods!(matchVarMods(var_mods, input_string),
                    input_string,
                    peptides,
                    group_id,
                    n
                    )
end