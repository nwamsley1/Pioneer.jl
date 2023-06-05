
struct Protein
    name::String
    pepGroup_ids::Set{UInt32}
end

Protein(name::String) = Protein(name, Set{UInt32}())
getName(p::Protein) = p.name
getPepGroupIDs(p::Protein) = p.pepGroup_ids

struct PeptideGroup
    prot_ids::Set{UInt32}
    pep_ids::Set{UInt32}
    sequence::String
end
PeptideGroup(sequence::String) = PeptideGroup(Set{UInt32}(), Set{UInt32}(), sequence)

getProtIDs(pg::PeptideGroup) = pg.prot_ids
getPepIDs(pg::PeptideGroup) = pg.pep_ids
getSeq(pg::PeptideGroup) = pg.sequence
addProtID!(pg::PeptideGroup, prot_id::UInt32) = push!(pg.prot_ids, prot_id)
addPepID!(pg::PeptideGroup, pep_id::UInt32) = push!(pg.pep_ids, pep_id)

struct Peptide
    sequence::String
    pepGroup_id::UInt32
    prec_ids::Set{UInt32}
end

getSeq(p::Peptide) = p.sequence
getPepGroupID(p::Peptide) = p.pepGroup_id
getPrecIDs(p::Peptide) = p.prec_ids
addPrecID!(p::Peptide, prec_id::UInt32) = push!(p.prec_ids, prec_id)

struct SimplePrecursor
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end

mutable struct PrecursorTable
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    id_to_pep::UnorderedDictionary{UInt32, Peptide}
    id_to_prec::UnorderedDictionary{UInt32, SimplePrecursor}
    prot_to_id::UnorderedDictionary{String, UInt32}
    peptides::Set{String}
end

PrecursorTable() = PrecursorTable(
    UnorderedDictionary{UInt32, Protein}(),
    UnorderedDictionary{UInt32, PeptideGroup}(),
    UnorderedDictionary{UInt32, Peptide}(),
    UnorderedDictionary{UInt32, SimplePrecursor}(),
    UnorderedDictionary{String, UInt32}(),
    Set{String}()
)

hasProtein(pt::PrecursorTable, prot::String) = haskey(pt.prot_to_id, prot)

function addProtein!(pt::PrecursorTable, prot::Protein, prot_id::UInt32) 
    insert!(pt.prot_to_id, getName(prot), prot_id)
    insert!(pt.id_to_prot, prot_id, prot)
end

getProtein(pt::PrecursorTable, prot_id::UInt32) = pt.id_to_prot[prot_id]
getProtID(pt::PrecursorTable, prot_name::String) = pt.prot_to_id[prot_name]
addPepGroupToProtein!(pt::PrecursorTable, prot_id::UInt32, pepGroup_id::UInt32) = push!(getPepGroupIDs(getProtein(pt, prot_id)), pepGroup_id)

getPepGroup(pt::PrecursorTable, pepGroup_id) = pt.id_to_pepGroup[pepGroup_id]
addPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, sequence::String) = insert!(pt.id_to_pepGroup, pepGroup_id, PeptideGroup(sequence))
addProtToPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, prot_id::UInt32) = addProtID!(getPepGroup(pt, pepGroup_id), prot_id)

getPep(pt::PrecursorTable, pep_id) = pt.id_to_pep[pep_id]
addPep!(pt::PrecursorTable, pep_id::UInt32, pepGroup_id::UInt32, sequence::String) = insert!(pt.id_to_pep, pep_id, Peptide(sequence, pepGroup_id, Set{UInt32}()))
addPrecToPep!(pt::PrecursorTable, pep_id::UInt32, prec_id::UInt32) = addPrecID!(getPep(pt, pep_id), prec_id)

getPrec(pt::PrecursorTable, prec_id::UInt32) = pt.id_to_prec[prec_id]
addPrec!(pt::PrecursorTable, prec_id::UInt32, pep_id::UInt32, charge::UInt8, isotope::UInt8) = insert!(pt.id_to_prec, prec_id, SimplePrecursor(charge, isotope, pep_id))

function getProtID(ptable::PrecursorTable, prot_name::String, max_prot_id::UInt32)
    if !hasProtein(ptable, prot_name)
        max_prot_id += UInt32(1)
        prot_id = max_prot_id 
        addProtein!(ptable, Protein(prot_name), prot_id)
    else
        prot_id = getProtID(ptable, prot_name)
    end
    return max_prot_id, prot_id
end

function fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    for mod in fixed_mods
        peptide = replace(peptide, mod[:p]=>mod[:r])
    end
    peptide
end

include("applyMods.jl")

function buildPrecursorTable!(ptable::PrecursorTable, peptides_fasta::Vector{FastaEntry},  
                                fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                                var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                                n::Int)
    max_prot_id, pepGroup_id, max_pep_id, prec_id = UInt32(1), UInt32(1), UInt32(1), UInt32(1)
    
    #Most recently encoutnered peptide sequence
    #It is assumed that `peptides_fasta` is sorted by the `sequence` field
    #so duplicated peptides are adjacent in the list
    previous_peptide = ""
    for peptide in peptides_fasta
        #Current peptide differs from previous
        if getSeq(peptide) != previous_peptide

            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(getSeq(peptide), fixed_mods)

            pepGroup_id += UInt32(1)
            addPepGroup!(ptable, pepGroup_id, peptide_sequence_fixed_mods)

            max_prot_id, prot_id = getProtID(ptable, 
                                             getID(peptide), 
                                             max_prot_id)
            addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
            
            addProtToPepGroup!(ptable, pepGroup_id, prot_id)

            previous_peptide = getSeq(peptide)

            ######
            #Apply Variable modifications
            #######
            max_pep_id = applyMods!(ptable.id_to_pep,
                        var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                        getSeq(peptide),
                        pepGroup_id,
                        max_pep_id,
                        n = n); 
        else #Duplicated peptide, so add a protein to the peptide group

            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(getSeq(peptide), fixed_mods)

            max_prot_id, prot_id = getProtID(ptable, peptide_sequence_fixed_mods, max_prot_id)
            addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
            addProtToPepGroup!(ptable, pepGroup_id, prot_id)
        end
    end
end