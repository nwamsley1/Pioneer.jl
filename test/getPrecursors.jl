
using Dictionaries 
fixed_mods = [(p=r"C", r="C[Carb]")]
var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
struct PeptideGroup
    prot_ids::Set{UInt32}
end

function addProt!(pg::PeptideGroup, prot_id::UInt32)
    push!(pg.prot_ids, prot_id)
end

getProtIDs(p::PeptideGroup) = p.prot_ids

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
addPepGroup!(p::Protein, pep_group_id::UInt32) = push!(p.pep_group_ids, pep_group_id)

PeptideGroup() = PeptideGroup("", Set{UInt32}(), UInt32(0))

struct PrecursorTable
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptide group
end

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, Protein}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, PeptideGroup}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, Peptide}()
                                    )

containsProt(p::PrecursorTable, protein::String) = isassigned(p.prot_to_id, protein)
containsPepGroup(p::PrecursorTable, peptide::String) = isassigned(p.pepGroup_to_id, peptide)
getProtID(p::PrecursorTable, protein::String) = p.prot_to_id[protein]
getPepGroupID(p::PrecursorTable, peptide::String) = p.pepGroup_to_id[peptide]
getPepGroup(p::PrecursorTable, pep_id::UInt32) = p.id_to_pepGroup[pep_id]
insertProtID!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.prot_to_id, protein, prot_id)
insertProt!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.id_to_prot, prot_id, Protein(protein))
insertPepGroupID!(p::PrecursorTable, peptide::String, pepGroup_id::UInt32) = insert!(p.pepGroup_to_id, peptide, pepGroup_id)
PeptideGroup(p::PrecursorTable, protein::String) = PeptideGroup(Set(getProtID(p, protein)))
insertPepGroup!(p::PrecursorTable, protein::String, pepGroup_id::UInt32) = insert!(p.id_to_pepGroup, pepGroup_id, PeptideGroup(p, protein))

addPepGroup!(p::PrecursorTable, prot_id::UInt32, pep_group_id::UInt32) = addPepGroup!(p.id_to_prot[prot_id], pep_group_id)

#Adds the protein_id to the pep group if not already included. 
function addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)
    addProt!(getPepGroup(p, getPepGroupID(p, peptide)), getProtID(p, protein))
end

#Adds the pep_group_id to the protein if not already included. 
function addPepGroupToProtein!(p::PrecursorTable, protein::String, peptide::String)
    addPepGroup!(p, getProtID(p, protein), getPepGroupID(p, peptide))
end

"""
    Adds a map from prot_id top protein
    and from protien to prot_id
"""
function newProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)
    insertProtID!(precursor_table, protein, prot_id);
    insertProt!(precursor_table, protein, prot_id);
end

"""
    Adds a map from pep_group_id to a PeptideGroup
    and from a PeptideGroup.name to a pep_group_id
"""
function newPeptideGroup!(peptide::String, pepGroup_id::UInt32, protein::String, precursor_table::PrecursorTable)
        insertPepGroupID!(precursor_table, peptide, pepGroup_id);
        insertPepGroup!(precursor_table, protein, pepGroup_id)
end


PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, Peptide}(),
                                    Set{String}(),
                                    Set{String}()
                                    )

function matchVarMods(patterns::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String)
    matches = Vector{Tuple{UnitRange{Int64}, String}}()
    for pattern in patterns
        regex_matches = findall(pattern[:p], input_string)
        for match in regex_matches
            push!(matches, (match, pattern[:r]))
        end
    end
    matches
end

function applyVariableMods!(matches::Vector{Tuple{UnitRange{Int64}, String}}, input_string::String, precursors::UnorderedDictionary{UInt32, Peptide}, group_id::UInt32, pep_id::UInt32, n::Int)
#out = Vector{String}(undef, length(all_combinations))
    # Apply the replacements to the input string for each combination
    for N in 1:min(n, length(matches))
        for (i, combination) in enumerate(combinations(matches, N))

            sort!(combination, by=match->match[1][end]);
            output_str = [""]#nput_string[1:combination[1][1][1]]
            index = 1
            for mod in combination
                push!(output_str, input_string[index:mod[1][1]])
                push!(output_str, mod[2])
                index = mod[1][1]+1
            end
            push!(output_str, input_string[index:end])
            #append!(out, [join(output_str)])
            pep_id += UInt32(1);
            insert!(precursors, pep_id, Peptide(input_string, group_id));
        end
    end
    pep_id += UInt32(1);
    insert!(precursors, pep_id, Peptide(input_string, group_id))
    pep_id
    #out
end

function applyMods!(var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String, peptides::UnorderedDictionary{UInt32, Peptide}, group_id::UInt32, pep_id::UInt32; n::Int = 3)
    applyVariableMods!(matchVarMods(var_mods, input_string),
                    input_string,
                    peptides,
                    group_id,
                    pep_id,
                    n
                    )
end

function fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    for mod in fixed_mods
        peptide = replace(peptide, mod[:p]=>mod[:r])
    end
    peptide
end

