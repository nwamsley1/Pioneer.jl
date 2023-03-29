include("./precursor.jl");
using Combinatorics, Dictionaries

export buildPrecursorTable!, getPrecursors!

struct PeptideGroup
    prot_ids::Set{UInt32}
end

PeptideGroup() = PeptideGroup(Set{UInt32}())
getProtIDs(p::PeptideGroup) = p.prot_ids
export getProtIds

function addProtID!(pg::PeptideGroup, prot_id::UInt32)
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

getPepGroupIDs(p::Protein) = p.pep_group_ids
getName(p::Protein) = p.name
Protein(name::String) = Protein(name, Set{UInt32}())
addPepGroup!(p::Protein, pep_group_id::UInt32) = push!(p.pep_group_ids, pep_group_id)

export getPepGroupIDs, getName
mutable struct PrecursorTable
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptide group
    precursors::Vector{Precursor} #Needs to be sortable by precursor mass, therfore, not an UnorderedDictioanry. 
end

getIDToProt(p::PrecursorTable) = p.id_to_prot
getProtToID(p::PrecursorTable) = p.prot_to_id
getIDToPepGroup(p::PrecursorTable) = p.id_to_pepGroup
getPepGroupToID(p::PrecursorTable) = p.pepGroup_to_id
getIDToPep(p::PrecursorTable) = p.id_to_pep
getPrecursors(p::PrecursorTable) = p.precursors

export getIDToProt, getIDToPep, getIDToPepGroup, getPepGroupToID, getIDToPep, getPrecursors

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, Protein}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, PeptideGroup}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, Peptide}(),
                                    Vector{Precursor}())

containsProt(p::PrecursorTable, protein::AbstractString) = isassigned(p.prot_to_id, protein)
containsProtID(p::PrecursorTable, prot_id::UInt32) = isassigned(p.id_to_prot, prot_id)
containsPepGroup(p::PrecursorTable, peptide::String) = isassigned(p.pepGroup_to_id, peptide)
containsPepGroupID(p::PrecursorTable, pepGroup_id::UInt32) = isassigned(p.id_to_pepGroup, pepGroup_id)
containsPepID(p::PrecursorTable, pep_id::UInt32) = isassigned(p.id_to_pep, pep_id)

getProtID(p::PrecursorTable, protein::String) = p.prot_to_id[protein]
getProt(p::PrecursorTable, prot_id::UInt32) = p.id_to_prot[prot_id]
getPepGroupID(p::PrecursorTable, peptide::String) = p.pepGroup_to_id[peptide]
getPepGroup(p::PrecursorTable, pepGroup_id::UInt32) = p.id_to_pepGroup[pepGroup_id]
getPep(p::PrecursorTable, pep_id::UInt32) = p.id_to_pep[pep_id]

insertProtID!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.prot_to_id, protein, prot_id)
insertProt!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.id_to_prot, prot_id, Protein(protein))
insertPepGroupID!(p::PrecursorTable, peptide::String, pepGroup_id::UInt32) = insert!(p.pepGroup_to_id, peptide, pepGroup_id)
insertPepGroup!(p::PrecursorTable, protein::String, pepGroup_id::UInt32) = insert!(p.id_to_pepGroup, pepGroup_id, PeptideGroup(Set(getProtID(p, protein))))

#PeptideGroup(p::PrecursorTable, protein::String) = PeptideGroup(Set(getProtID(p, protein)))


#addNewPepGroup!(p::PrecursorTable, prot_id::UInt32, pep_group_id::UInt32) = addPepGroup!(p.id_to_prot[prot_id], pep_group_id)

#Adds the protein_id to the pep group if not already included. 
function addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)
    addProtID!(getPepGroup(p, getPepGroupID(p, peptide)), getProtID(p, protein))
end

#Adds the pep_group_id to the protein if not already included. 
function addPepGroupToProtein!(p::PrecursorTable, protein::String, peptide::String)
    prot_id = getProtID(p, protein)
    addPepGroup!(getProt(p, prot_id), getPepGroupID(p, peptide))
end

"""
    Adds a map from prot_id top protein
    and from protien to prot_id
"""
function addNewProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)
    insertProtID!(precursor_table, protein, prot_id);
    insertProt!(precursor_table, protein, prot_id);
end

"""
    Adds a map from pep_group_id to a PeptideGroup
    and from a PeptideGroup.name to a pep_group_id
"""
function addNewPeptideGroup!(peptide::String, pepGroup_id::UInt32, protein::String, precursor_table::PrecursorTable)
        insertPepGroupID!(precursor_table, peptide, pepGroup_id);
        insertPepGroup!(precursor_table, protein, pepGroup_id)
end

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

function applyVariableMods!(matches::Vector{Tuple{UnitRange{Int64}, String}}, unmod_seq::String, peptides::UnorderedDictionary{UInt32, Peptide}, group_id::UInt32, pep_id::UInt32, n::Int)
    
    function applyMods(combination, unmod_seq::String)
        output_str = [""] #Build the modified sequence from scratch
        index = 1
        for mod in combination
            push!(output_str, unmod_seq[index:mod[1][1]]) #From the prior modification up until the location of the next mod
            push!(output_str, mod[2])                     #Add the modification
            index = mod[1][1]+1                           
        end
        push!(output_str, unmod_seq[index:end]) #From the last mod to the end of the sequence
        return output_str        
    end

    insertPeptide!(peptides::UnorderedDictionary{UInt32, Peptide}, 
                   seq::String, 
                   group_id::UInt32, 
                   pep_id::UInt32) = insert!(peptides, pep_id, Peptide(seq, group_id))

    #Peptide with 0 variable mods. 
    #pep_id += UInt32(1);
    insertPeptide!(peptides, unmod_seq, group_id, pep_id);
    pep_id += UInt32(1);
    # Apply the replacements to the input string for each combination
    for N in 1:min(n, length(matches)) #Apply 1:min(n, length(matches)) variable mods
        for combination in combinations(matches, N) #Each combination of "N" variable mods
            sort!(combination, by=match->match[1][end]); #Sort applied variable mods in order of appearance in "unmod_seq". 
            modified_seq = applyMods(combination, unmod_seq) #Get the modified sequence for the given combination of mods
            #pep_id += UInt32(1);
            insertPeptide!(peptides, join(modified_seq), group_id, pep_id);
            pep_id += UInt32(1);
        end
    end
    pep_id
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
export applyMods!

function fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    for mod in fixed_mods
        peptide = replace(peptide, mod[:p]=>mod[:r])
    end
    peptide
end

#For tomorrow. Make this function. Make it sort the precursors. Also make tests for this module.
function getPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})
    prec_id = UInt32(1)
    for (pep_id, peptide) in pairs(ptable.id_to_pep)
        for charge in charges, isotope in isotopes
            push!(ptable.precursors, Precursor(peptide.sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict))
            prec_id+= UInt32(1)
        end
    end
    sort!(ptable.precursors, by=p->getMZ(p))
end

function buildPrecursorTable!(ptable::PrecursorTable, 
                              fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                              var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                              n::Int, f_path::String)
    open(f_path) do f #"./data/NRF2_SIL.txt"
        pepGroup_id, pep_id, prot_id = UInt32(1), UInt32(1), UInt32(1)
        timetaken = @elapsed for (row, protein_peptide) in enumerate(eachline(f))
            protein, peptide = map(string, split(protein_peptide, "\t")); #Parse input "PROTEIN_NAME\tPEPTIDE_SEQUENCE"
            peptide = fixedMods(peptide, fixed_mods); #Apply fixed modifications
            if !containsProt(ptable, protein) #If the protien hasn't been encountered,
                #then add it to the hash table
                addNewProtein!(protein, prot_id, ptable);
                prot_id += UInt32(1);  
            end
            if !containsPepGroup(ptable, peptide) #If the peptide hasn't been encountered, then
                addNewPeptideGroup!(peptide, pepGroup_id, protein, ptable); #add it to the hash table,
                addPepGroupToProtein!(ptable, protein, peptide); #add a new peptide group to the protein,
                pep_id = applyMods!(var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                                    peptide,               #and increase the pep_id for each variable mod applied 
                                    ptable.id_to_pep,
                                    pepGroup_id,
                                    pep_id,
                                    n = n); #
                pepGroup_id += UInt32(1); 
            else #If this peptide has been encountered before, we don't need to apply the variable modes. Instead,
                addProteinToPepGroup!(ptable, protein, peptide); #Add the current protein to this peptide group
                addPepGroupToProtein!(ptable, protein, peptide); #Add the current peptide group to this protein
            end
        end
        println("Time to build precursor table ", timetaken);
    end
end
