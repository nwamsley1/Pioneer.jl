"""
    Protein

Type that represents a protein. 

### Fields

- name::String -- Name of the protein
- pepGroup_ids::UInt32-- Identifiers for the `PeptideGroup`s linked to this `Protein`

### Examples

- `Protein(name::String) = Protein(name, Set{UInt32}())` -- constructor for an placeholder 

### GetterMethods

- getName(p::Protein) = p.name
- getPepGroupIDs(p::Protein) = p.pepGroup_ids

### Methods

- addPepGroupID!(p::Protein, pepGroup_id::UInt32)

### Notes

Protein > PeptideGroup > Peptide > Precursor
"""
struct Protein
    name::String
    pepGroup_ids::Set{UInt32}
end

Protein(name::String) = Protein(name, Set{UInt32}())
getName(p::Protein) = p.name
getPepGroupIDs(p::Protein) = p.pepGroup_ids
addPepGroupID!(p::Protein, pepGroup_id::UInt32) = push!(getPepGroupIDs(p), pepGroup_id)

"""
    PeptideGroup

Type that represents an unmodified peptide. Minimal unit specifiable by single-charater amino-acid  
symbols alone. Cannot describe modifications. 

### Fields

- prot_ids::Set{UInt32} -- Set of identifiers for `Protein`s that contain this `Peptide` group
- pep_ids::Set{UInt32} -- Set of identifiers for `Peptide`s contained by this `PeptideGroup`
- sequence::String -- Sequence of the `PeptideGroup`

### Examples

- `PeptideGroup() = PeptideGroup(Set{UInt32}(), "")` -- constructor for an placeholder 

### GetterMethods

- getSeq(p::PeptideGroup) = p.sequence
- getProtIDs(p::PeptideGroup) = p.prot_ids
- getPepIDs(p::PeptideGroup) = p.pep_ids

### Methods

- getProtIDs(pg::PeptideGroup) = pg.prot_ids
- getPepIDs(pg::PeptideGroup) = pg.pep_ids
- getSeq(pg::PeptideGroup) = pg.sequence
- addProtID!(pg::PeptideGroup, prot_id::UInt32) = push!(pg.prot_ids, prot_id)
- addPepID!(pg::PeptideGroup, pep_id::UInt32) = push!(pg.pep_ids, pep_id)

### Notes

Protein > PeptideGroup > Peptide > Precursor

"""
struct PeptideGroup
    prot_ids::Set{UInt32}
    pep_ids::Set{UInt32}
    sequence::String
    decoy::Bool
end
PeptideGroup(sequence::String, decoy::Bool) = PeptideGroup(Set{UInt32}(), Set{UInt32}(), sequence, decoy)

getProtIDs(pg::PeptideGroup) = pg.prot_ids
getPepIDs(pg::PeptideGroup) = pg.pep_ids
getSeq(pg::PeptideGroup) = pg.sequence
isDecoy(pg::PeptideGroup) = pg.decoy
addProtID!(pg::PeptideGroup, prot_id::UInt32) = push!(pg.prot_ids, prot_id)
addPepID!(pg::PeptideGroup, pep_id::UInt32) = push!(pg.pep_ids, pep_id)

"""
    Peptide

Type that represents a peptide. Minimal describing a peptide of a particular mono-isotopic mass. A `Precursor`
is a `Peptide` specifying a specific mass and charge.

### Fields

- sequence::String -- Sequence of the peptide
- pepGroup_id::UInt32 -- Identifier of the `PeptideGroup` to which this `Peptide` belongs. 
- prec_ids::Set{UInt32} -- Set of identifiers for `Precursor`s associated with this `Peptide`


### GetterMethods

- getSeq(p::Peptide) = p.sequence
- getPepGroupID(p::Peptide) = pepGroup_id
- getPrecIDs(p::Peptide) = p.prec_ids

### Methods

- addPrecID!(p::Peptide, prec_id::UInt32)

### Notes

Protein > PeptideGroup > Peptide > Precursor

"""
struct Peptide
    sequence::String
    pepGroup_id::UInt32
    prec_ids::Set{UInt32}
    decoy::Bool
end

getSeq(p::Peptide) = p.sequence
getPepGroupID(p::Peptide) = p.pepGroup_id
getPrecIDs(p::Peptide) = p.prec_ids
isDecoy(p::Peptide) = p.decoy
addPrecID!(p::Peptide, prec_id::UInt32) = push!(p.prec_ids, prec_id)

struct SimplePrecursor
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end

abstract type PeptideDatabase end

abstract type PrecursorDatabase <: PeptideDatabase end

"""
    PrecursorDatabase

Data Structure that represents relations between precursors, peptides, peptide groups, and proteins

### Suggested Fields

An implementation of `PrecursorTable` should have these fields 

- `id_to_prot::UnorderedDictionary{UInt32, Protein}` -- Maps from a protien identifier to a Protein
- `id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}` -- Maps a PeptideGroup identifier to a PeptideGroup
- `id_to_pep::UnorderedDictionary{UInt32, Peptide}` -- Maps a peptide identifier to a peptide
- `id_to_prec::Dictionary{UInt32, Precursor}` -- Precursor has fields `pep_id` and `prec_id`. `pep_id`s are keys for `id_to_pep`
- `prot_to_id::UnorderedDictionary{String, UInt32}` -- 
- `peptides::Set{String}` -- 

### GetterMethods

- getProtein(pt::PrecursorTable, prot_id::UInt32)
- getProtID(pt::PrecursorTable, prot_name::String)
- getPepGroup(pt::PrecursorTable, pepGroup_id)
- getPrec(pt::PrecursorTable, prec_id::UInt32
- getProtID!(ptable::PrecursorTable, prot_name::String, max_prot_id::UInt32)

### Methods

- addPepGroupToProtein!(pt::PrecursorTable, prot_id::UInt32, pepGroup_id::UInt32)
- addPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, sequence::String)
- addProtToPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, prot_id::UInt32)
- addPep!(pt::PrecursorTable, pep_id::UInt32, pepGroup_id::UInt32, sequence::String)
- addPrecToPep!(pt::PrecursorTable, pep_id::UInt32, prec_id::UInt32)
- addPrec!(pt::PrecursorTable, prec_id::UInt32, pep_id::UInt32, charge::UInt8, isotope::UInt8)
### Notes

`PrecursorDatabase` keeps an account of relations between Proteins, PeptideGroups, Peptides, and
Precursors. Each Protein, PeptideGroup, Peptide, and Precursor has a unique UInt32 indentifier. 
Each also has a string name. See `Protein`, `PeptideGroup`, `Peptide` and `Precursor`. 
Proteins map to one or more PeptideGroups. PeptideGroups map to one or more Peptides, 
and each Peptide maps to one or more precursors.
"""
mutable struct PrecursorTable <: PrecursorDatabase
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

"""
    addProtein!(pt::PrecursorTable, prot::Protein, prot_id::UInt32) 

Adds a new protein to `pt`. Updates the pt.id_to_prot dictionary and the 
pt.prot_to_id dictionary. 
"""
function addProtein!(pt::PrecursorTable, prot::Protein, prot_id::UInt32) 
    insert!(pt.prot_to_id, getName(prot), prot_id)
    insert!(pt.id_to_prot, prot_id, prot)
end

getProtein(pt::PrecursorTable, prot_id::UInt32) = pt.id_to_prot[prot_id]
getProtID(pt::PrecursorTable, prot_name::String) = pt.prot_to_id[prot_name]
addPepGroupToProtein!(pt::PrecursorTable, prot_id::UInt32, pepGroup_id::UInt32) = push!(getPepGroupIDs(getProtein(pt, prot_id)), pepGroup_id)

getPepGroup(pt::PrecursorTable, pepGroup_id) = pt.id_to_pepGroup[pepGroup_id]
addPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, sequence::String, decoy::Bool) = insert!(pt.id_to_pepGroup, pepGroup_id, PeptideGroup(sequence, decoy))
addProtToPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, prot_id::UInt32) = addProtID!(getPepGroup(pt, pepGroup_id), prot_id)

getPep(pt::PrecursorTable, pep_id) = pt.id_to_pep[pep_id]
addPep!(pt::PrecursorTable, pep_id::UInt32, pepGroup_id::UInt32, sequence::String) = insert!(pt.id_to_pep, pep_id, Peptide(sequence, pepGroup_id, Set{UInt32}()))
addPrecToPep!(pt::PrecursorTable, pep_id::UInt32, prec_id::UInt32) = addPrecID!(getPep(pt, pep_id), prec_id)

getPrec(pt::PrecursorTable, prec_id::UInt32) = pt.id_to_prec[prec_id]
addPrec!(pt::PrecursorTable, prec_id::UInt32, pep_id::UInt32, charge::UInt8, isotope::UInt8) = insert!(pt.id_to_prec, prec_id, SimplePrecursor(charge, isotope, pep_id))

function getProtFromPepID(pt::PrecursorTable, pep_id::Int)
    getProtIDs(pt.id_to_pepGroup[getPepGroupID(pt.id_to_pep[pep_id])])
end
"""
    getProtID!(ptable::PrecursorTable, prot_name::String, max_prot_id::UInt32)

Checks if a protein named `prot_name` is in the ptable. 
- If yes, returns `max_prot_id` unchanged and the protein ID of the protein named `prot_name`
- If no, adds a new protein to the ptable with `prot_name` with identifier `max_prot_id` += 1. Returns
max_prot_id + 1, max_prot_id + 1. 
"""
function getProtID!(ptable::PrecursorTable, prot_name::String, max_prot_id::UInt32)
    if !hasProtein(ptable, prot_name) #No protein with `prot_name` found in `ptable`
        max_prot_id += UInt32(1) #Update max_prot_id since adding a new protein
        prot_id = max_prot_id #Current prot_id is the max_prot_id
        addProtein!(ptable, Protein(prot_name), prot_id) #Add a new protein
    else
        prot_id = getProtID(ptable, prot_name) #Retrieve the ID of the protein with name `prot_name`
    end
    return max_prot_id, prot_id
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
        sequence = getSeq(peptide)
        #Current peptide differs from previous
        if sequence != previous_peptide
            decoy = isDecoy(peptide)
            if decoy
                sequence = shufflefast(sequence)
            end
            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(sequence, fixed_mods)

            pepGroup_id += UInt32(1)
            addPepGroup!(ptable, pepGroup_id, peptide_sequence_fixed_mods, decoy)

            max_prot_id, prot_id = getProtID!(ptable, 
                                            getID(peptide), 
                                            max_prot_id)
            if !decoy
                addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
                addProtToPepGroup!(ptable, pepGroup_id, prot_id)
                previous_peptide = sequence
            end

            ######
            #Apply Variable modifications
            #######
            max_pep_id = applyMods!(ptable.id_to_pep,
                        var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                        peptide_sequence_fixed_mods,
                        pepGroup_id,
                        max_pep_id,
                        decoy,
                        n = n); 
        else #Duplicated peptide, so add a protein to the peptide group

            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(getSeq(peptide), fixed_mods)

            max_prot_id, prot_id = getProtID!(ptable, peptide_sequence_fixed_mods, max_prot_id)
            addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
            addProtToPepGroup!(ptable, pepGroup_id, prot_id)
        end
    end
end
