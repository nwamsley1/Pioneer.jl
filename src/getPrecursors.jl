
include("./precursor.jl");
using Combinatorics, Dictionaries

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

- addProtID!(pg::PeptideGroup, prot_id::UInt32)

### Notes

Protein > PeptideGroup > Peptide > Precursor

"""
struct PeptideGroup
    prot_ids::Set{UInt32}
    pep_ids::Set{UInt32}
    sequence::String
end

PeptideGroup() = PeptideGroup(Set{UInt32}(), Set{UInt32}(), "")
getSeq(pg::PeptideGroup) = pg.sequence
getProtIDs(pg::PeptideGroup) = pg.prot_ids
getPepIDs(pg::PeptideGroup) = pg.pep_ids

addProtID!(pg::PeptideGroup, prot_id::UInt32) push!(pg.prot_ids, prot_id)

export getSeq, getProtIDs, getPepIDs, addProtID!

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
end

getSeq(p::Peptide) = p.sequence
getPepGroupID(p::Peptide) = pepGroup_id
getPrecIDs(p::Peptide) = p.prec_ids

addPrecID!(p::Peptide, prec_id::UInt32) = push!(p.prec_ids, prec_id)

export getSeq, getPepGroupID, getPrecIDs, addPrecID!

"""
    Protein

Type that represents a protein. 

### Fields

- sequence::String -- Sequence of the peptide
- pep_group_id::UInt32-- Identifier of thes `PeptideGroup`s contained by this `Protein`

### Examples

- `Protein(name::String) = Protein(name, Set{UInt32}())` -- constructor for an placeholder 

### GetterMethods

- getName(p::Protein) = p.name
- getPepGroupIDs(p::Protein) = p.pep_group_ids

### Methods

- addPepGroupID!(p::Protein, pep_group_id::UInt32)

### Notes

Protein > PeptideGroup > Peptide > Precursor
"""
struct Protein
    name::String
    pep_group_ids::Set{UInt32}
end

Protein(name::String) = Protein(name, Set{UInt32}())

getName(p::Protein) = p.name
getPepGroupIDs(p::Protein) = p.pep_group_ids

addPepGroupID!(p::Protein, pep_group_id::UInt32) = push!(p.pep_group_ids, pep_group_id)

export getPepGroupIDs, getName, addPepGroupID!

abstract type PeptideDatabase end

"""
    PrecursorTable

Data Structure that represents relations between precursors, peptides, peptide groups, and proteins

### Suggested Fields

An implementation of `PrecursorTable` should have these fields 

- id_to_prot::UnorderedDictionary{UInt32, Protein}-- Maps from a protien identifier to a Protein
- prot_to_id::UnorderedDictionary{String, UInt32} -- Maps from a protein name to a protein identifier
- id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup} -- Maps a PeptideGroup identifier to a PeptideGroup
- pepGroup_to_id::UnorderedDictionary{String, UInt32} -- Maps a PeptideGroup name/sequence to an identifier
- id_to_pep::UnorderedDictionary{UInt32, Peptide} -- Maps a peptide identifier to a peptide
- id_to_prec::Dictionary{UInt32, Precursor} -- Precursor has fields `pep_id` and `prec_id`. `pep_id`s are keys for `id_to_pep`
- sorted_prec_ids::Vector{UInt32} -- 
- prec_id_to_transitions::Dictionary{UInt32, Vector{Transition}} --

### GetterMethods

- getIDToProt(p::PrecursorDatabase) = p.id_to_prot
- getProtToID(p::PrecursorDatabase) = p.prot_to_id
- getIDToPepGroup(p::PrecursorDatabase) = p.id_to_pepGroup
- getPepGroupToID(p::PrecursorDatabase) = p.pepGroup_to_id
- getIDToPep(p::PrecursorDatabase) = p.id_to_pep
- getIDToPrec(p::PrecursorDatabase) = p.id_to_prec
- getPrecIDToTransitions(p::PrecursorDatabase) = p.prec_id_to_transitions
- getPrecursorIDs(p::PrecursorDatabase) = p.sorted_prec_ids
- getProtID(p::PrecursorDatabase, protein::String)
- getProtID(p::PrecursorDatabase, protein::Protein)
- getPepGroup(p::PrecursorDatabase, pepGroup_id::UInt32)
- getPepGroupID(p::PrecursorDatabase, peptide::String)
- getPepGroupID(p::PrecursorDatabase, pepGroup::PeptideGroup)
- getPep(p::PrecursorDatabase, pep_id::UInt32)
- getPrecursor(p::PrecursorDatabase, prec_id::UInt32)
- getTransitions(p::PrecursorDatabase, prec_id::UInt32)
- getProtNamesFromPepSeq(p::PrecursorDatabase, peptide::String)
- getProtNamesFromPepSeq(p::PrecursorDatabase, pepGroup::PeptideGroup)
- getPepGroupsFromProt(p::PrecursorDatabase, protein::String)
- getPepSeqsFromProt(p::PrecursorDatabase, protein::String)
- getPepGroupsFromProt(p::PrecursorDatabase, prot_id::UInt32)
- getPepSeqsFromProt(p::PrecursorDatabase, prot_id::UInt32)
- getPepIDFromPrecID(p::PrecursorDatabase, prec_id::UInt32)

### Methods

- insertProtID!(p::PrecursorDatabase, protein::String, prot_id::UInt32)
- insertProt!(p::PrecursorDatabase, protein::String, prot_id::UInt32)
- insertPepGroupID!(p::PrecursorDatabase, peptide::String, pepGroup_id::UInt32)
- insertPepGroup!(p::PrecursorDatabase, protein::String, peptide::String, pepGroup_id::UInt32)
- setSortedPrecursorKeys!(p::PrecursorDatabase)
- precursorRangeQuery(p::PrecursorDatabase, window_center::Float32, left_precursor_tolerance::Float32, right_precursor_tolerance::Float32)
- addProteinToPepGroup!(pd::PrecursorDatabase, protein::String, peptide::String)
- addPepGroupToProtein!(pd::PrecursorDatabase, protein::String, peptide::String)
- addNewProtein!(pd::PrecursorDatabase, protein::String, prot_id::UInt32)
- addNewPeptideGroup!(pd::PrecursorDatabase, peptide::String, pepGroup_id::UInt32, protein::String)
- addPrecursors!(ptable::PrecursorDatabase, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})
- addTransitions!(ptable::PrecursorDatabase, transition_charges::Vector{UInt8}, transition_isotopes::Vector{UInt8},b_start::Int64,y_start::Int64, fragment_match_ppm::Float32)
### Notes

`PrecursorDatabase` keeps an account of relations between Proteins, PeptideGroups, Peptides, and
Precursors. Each Protein, PeptideGroup, Peptide, and Precursor has a unique UInt32 indentifier. 
Each also has a string name. See `Protein`, `PeptideGroup`, `Peptide` and `Precursor`. 
Proteins map to one or more PeptideGroups. PeptideGroups map to one or more Peptides, 
and each Peptide maps to one or more precursors.
"""
abstract type PrecursorDatabase <: PeptideDatabase end

"""
    PrecursorTable <: PrecursorDatabase

Inherits from `PrecursorDatabase`. Minimal implementation of methods required for `PrecursorDatabase`

### Fields

- id_to_prot::UnorderedDictionary{UInt32, Protein}-- Maps from a protien identifier to a `Protein``
- prot_to_id::UnorderedDictionary{String, UInt32} -- Maps from a protein name to a protein identifier
- id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup} -- Maps a peptide group identifier to a `PeptideGroup``
- pepGroup_to_id::UnorderedDictionary{String, UInt32} -- Maps a `PeptideGroup`` name/sequence to an identifier
- id_to_pep::UnorderedDictionary{UInt32, Peptide} -- Maps a peptide identifier to a `Peptide`
- id_to_prec::Dictionary{UInt32, Precursor} -- Maps a precursor identifier to a `Precursor`
- sorted_prec_ids::Vector{UInt32} -- Precursor ID's for precursors in `id_to_prec` sorted by m/z 
- prec_id_to_transitions::Dictionary{UInt32, Vector{Transition}} --  Maps a precursor identifier to a Vector{Transition} (fragment ions)

### Examples

- `PrecursorTable() = PrecursorTable(
    UnorderedDictionary{UInt32, Protein}(),
    UnorderedDictionary{String, UInt32}(),
    UnorderedDictionary{UInt32, PeptideGroup}(),
    UnorderedDictionary{String, UInt32}(),
    UnorderedDictionary{UInt32, Peptide}(),
    Vector{Precursor}())` -- constructor for a placeholder 

"""
mutable struct PrecursorTable <: PrecursorDatabase
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide}
    id_to_prec::Dictionary{UInt32, Precursor}
    sorted_prec_ids::Vector{UInt32}
    prec_id_to_transitions::Dictionary{UInt32, Vector{Transition}}
end

getIDToProt(p::PrecursorDatabase) = p.id_to_prot
getProtToID(p::PrecursorDatabase) = p.prot_to_id
getIDToPepGroup(p::PrecursorDatabase) = p.id_to_pepGroup
getPepGroupToID(p::PrecursorDatabase) = p.pepGroup_to_id
getIDToPep(p::PrecursorDatabase) = p.id_to_pep
getIDToPrec(p::PrecursorDatabase) = p.id_to_prec
getPrecIDToTransitions(p::PrecursorDatabase) = p.prec_id_to_transitions
getPrecursorIDs(p::PrecursorDatabase) = p.sorted_prec_ids




getPrecIDToTransitions(p::PrecursorDatabase) = p.prec_id_to_transitions

export getIDToProt, getIDToPep, getIDToPepGroup, getPepGroupToID, getIDToPep, getIDToPrec, getPrecursors, getPrecIDToTransitions

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, Protein}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, PeptideGroup}(),
                                    UnorderedDictionary{String, UInt32}(),
                                    UnorderedDictionary{UInt32, Peptide}(),
                                    #UnorderedDictionary{UInt32, SimplePrecursor}(),
                                    Dictionary{UInt32, Precursor}(),
                                    Vector{UInt32}(),
                                    Dictionary{UInt32, Vector{Transition}}())

containsProt(p::PrecursorDatabase, protein::AbstractString) = isassigned(getProtToID(p), protein)
containsProtID(p::PrecursorDatabase, prot_id::UInt32) = isassigned(getIDToProt(p), prot_id)
containsPepGroup(p::PrecursorDatabase, peptide::String) = isassigned(getPepGroupToIDp(p), peptide)
containsPepGroupID(p::PrecursorDatabase, pepGroup_id::UInt32) = isassigned(getIDToPepGroup(p), pepGroup_id)
containsPepID(p::PrecursorDatabase, pep_id::UInt32) = isassigned(getIDToPep(p), pep_id)
containsPep(p::PrecursorDatabase, sequence::String) = isassigned(getPepSeqToPepID(p), sequence)
containsPrecID(p::PrecursorDatabase, prec_id::UInt32) = isassigned(getPrecursors(p), prec_id)
hasTransitions(p::PrecursorDatabase, prec_id::UInt32) = isassigned(getPrecIDToTransitions(p), prec_id)


getProt(p::PrecursorDatabase, prot_id::UInt32) = getIDToProt(p)[prot_id]

getProtID(p::PrecursorDatabase, protein::String) = getProtToID(p)[protein]
getProtID(p::PrecursorDatabase, protein::Protein) = getProtToID(getName(p))[protein]

getPepGroup(p::PrecursorDatabase, pepGroup_id::UInt32) = getIDToPepGroup(p)[pepGroup_id]
getPepGroupID(p::PrecursorDatabase, peptide::String) = getPepGroupToID(p)[string(replace(peptide, r"\[[^\]]*\]"=>""))]
getPepGroupID(p::PrecursorDatabase, pepGroup::PeptideGroup) = getPepGroupToID(p)[getSeq(pepGroup)]

getPep(p::PrecursorDatabase, pep_id::UInt32) = getIDToPep(p)[pep_id]
getPrecursor(p::PrecursorDatabase, prec_id::UInt32) = getIDToPrec(p)[prec_id]

getTransitions(p::PrecursorDatabase, prec_id::UInt32) = getPrecIDToTrasnitions(p)[prec_id]


getProtNamesFromPepSeq(p::PrecursorDatabase, peptide::String) = Set([getName(getProt(p, prot_id)) for prot_id in getProtIDs(getPepGroup(p, getPepGroupID(p, peptide)))])
getProtNamesFromPepSeq(p::PrecursorDatabase, pepGroup::PeptideGroup) = Set([getName(getProt(p, prot_id)) for prot_id in getProtIDs(getSeq(pepGroup))])

getPepGroupsFromProt(p::PrecursorDatabase, protein::String) = [getPepGroup(p, ID) for ID in getPepGroupIDs(getProt(p, getProtID(p, protein)))]
getPepSeqsFromProt(p::PrecursorDatabase, protein::String) = [getSeq(pep_group) for pep_group in getPepGroupsFromProt(p, protein)]
getPepGroupsFromProt(p::PrecursorDatabase, prot_id::UInt32) = [getPepGroup(p, ID) for ID in getPepGroupIDs(getProt(p, prot_id))]
getPepSeqsFromProt(p::PrecursorDatabase, prot_id::UInt32) = [getSeq(pep_group) for pep_group in getPepGroupsFromProt(p, prot_id)]
getPepIDFromPrecID(p::PrecursorDatabase, prec_id::UInt32) = getPepID(getPrecursor(p, prec_id))

export getProtID, getProt, getPepGroupID, getPepGroup, getPep

insertProtID!(p::PrecursorDatabase, protein::String, prot_id::UInt32) = insert!(p.prot_to_id, protein, prot_id)
insertProt!(p::PrecursorDatabase, protein::String, prot_id::UInt32) = insert!(p.id_to_prot, prot_id, Protein(protein))
insertPepGroupID!(p::PrecursorDatabase, peptide::String, pepGroup_id::UInt32) = insert!(p.pepGroup_to_id, peptide, pepGroup_id)
insertPepGroup!(p::PrecursorDatabase, protein::String, peptide::String, pepGroup_id::UInt32) = insert!(p.id_to_pepGroup, pepGroup_id, PeptideGroup(Set(getProtID(p, protein)), Set(UInt32[]), peptide))
function insertPep!(p::PrecursorDatabase, sequence::String, pep_id::UInt32, pepGroup_id::UInt32)
    insert!(getPepSeqToPepID(p), sequence, pep_id)
    insert!(getIDToPep(ptable), pep_id, Peptide(sequence, pepGroup_id, Set(UInt32[])))
end

"""
    setSortedPrecursorKeys!(p::PrecursorDatabase)

Sorts precursors in the `PrecursorDatabase` by their m/z ratio. First sorts the Dictionary{UInt32, Precursor}
`p.id_to_prec` that maps precursor id's to `Precursor`s. Sets `p.sorted_prec_ids`
as a vector of the ID's for these precursors (sorted). 
"""
function setSortedPrecursorKeys!(p::PrecursorDatabase)
    function setSortedPrecursorIDs(p::PrecursorDatabase, keys::Vector{UInt32}) 
        p.sorted_prec_ids = keys
    end
    sort!(getIDToPrec(p), by = prec->getMZ(prec));
    setSortedPrecursorKeys(p, collect(keys(getIDToPrec(p))))
end

"""
    precursorRangeQuery(p::PrecursorDatabase, window_center::Float32, left_precursor_tolerance::Float32, right_precursor_tolerance::Float32)

Finds precursor IDs mapping to `Precursor`s with m/z ratios in the tolerance specified by `window_center`, `left_precursor_tolerance`, and `right_precursor_tolerance`
Assumes the p.id_to_prec dictionary is already sorted. see `setSortedPrecursorKeys!`
"""
function precursorRangeQuery(p::PrecursorDatabase, window_center::Float32, left_precursor_tolerance::Float32, right_precursor_tolerance::Float32)
    l_bnd, u_bnd = window_center - left_precursor_tolerance, window_center + right_precursor_tolerance
    start = searchsortedfirst(getPrecursorIDs(p), l_bnd,lt=(t,x)->getMZ(getPrecursor(p, t))<x)
    stop = searchsortedlast(getPrecursorIDs(p), u_bnd,lt=(x,t)->getMZ(getPrecursor(p, t))>x)
    return @view(getPrecursorIDsp)[start:stop]
end

"""
addProteinToPepGroup!(pd::PrecursorDatabase, protein::String, peptide::String)

For the PeptideGroup to which "peptide" belongs, add the protein to its `prot_ids`
"""
function addProteinToPepGroup!(pd::PrecursorDatabase, protein::String, peptide::String)
    addProtID!(getPepGroup(pd, getPepGroupID(pd, peptide)), getProtID(pd, protein))
end

"""
addPepGroupToProtein!(pd::PrecursorDatabase, protein::String, peptide::String)

For the protein to which the `peptide` belongs, add the correct peptide group id to its `pep_group_ids`
"""
function addPepGroupToProtein!(pd::PrecursorDatabase, protein::String, peptide::String)
    addPepGroup!(getProt(pd, getProtID(pd, protein)), getPepGroupID(pd, peptide))
end

"""
addNewProtein!(pd::PrecursorDatabase, protein::String, prot_id::UInt32)

Create a new `Protein` and protein id and add them to the lookup tables
"""
function addNewProtein!(pd::PrecursorDatabase, protein::String, prot_id::UInt32)
    insertProtID!(pd, protein, prot_id);
    insertProt!(pd, protein, prot_id);
end

"""
addNewPeptideGroup!(pd::PrecursorDatabase, peptide::String, pepGroup_id::UInt32, protein::String)

Create a new `PeptideGroup` and peptide group id and add them to the lookup tables
"""
function addNewPeptideGroup!(pd::PrecursorDatabase, peptide::String, pepGroup_id::UInt32, protein::String)
        insertPepGroupID!(pd, peptide, pepGroup_id);
        insertPepGroup!(pd, protein, peptide, pepGroup_id)
end

"""
    addPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})

Modifies a `PrecursorDatabase` to fill the `precursors` field. Gets a `Precursor` for each unique `Peptide` in `ptable`. Sorts the 
precursors in order of increasing mz.   

### Input

- `ptable::PrecursorTable` -- See `PrecursorTable`
- `charges::Vector{UInt8}` -- Get a precursor with these charges for each Peptide in `ptable`
- `isotopes::Vector{UInt8}` -- Get a precursor with these isotopes for each Peptide in `ptable`
- `mods_dict::Dict{String, Float32}` -- Dict mapping modifications to their masses

### Output
- Fills the `id_to_prec` and `sorted_prec_ids` fields in a `PrecursorDatabase` 

### Examples 
test_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464),
    "Harg" => Float32(10),
    "Hlys" => Float32(8),
)
fixed_mods = [(p=r"C", r="C[Carb]")]
var_mods = [(p=r"(K\$)", r="[Hlys]"), (p=r"(R\$)", r="[Harg]")]

testPtable = PrecursorTable()
buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "../data/peptide_lists/PROT_PEPTIDE_TEST1.txt")
getPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
"""
function addPrecursors!(ptable::PrecursorDatabase, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})
    prec_id = UInt32(1)
    for (pep_id, peptide) in pairs(getIDToPep(ptable))
        for charge in charges, isotope in isotopes
            insert!(getIDToPrec(ptable), 
                    prec_id,
                    Precursor(peptide.sequence, 
                                charge = charge, 
                                isotope = isotope, 
                                pep_id = pep_id, 
                                prec_id = prec_id, 
                                mods_dict = mods_dict)
                    )
            #insert!(ptable.simple_precursors, prec_id, SimplePrecursor(getMZ(ptable.precursors[end]), charge, isotope, pep_id));
            prec_id+= UInt32(1)
        end
    end
    setSortedPrecursorKeys!(ptable)
end


"""
    addTransitions!(ptable::PrecursorDatabase, transition_charges::Vector{UInt8}, transition_isotopes::Vector{UInt8},b_start::Int64,y_start::Int64, fragment_match_ppm::Float32)

Modifies a `PrecursorDatabase` to add `Transition`s to the 'prec_id_to_transitions' field corresponding to all `Precursor`s  in the `id_to_prec` field. Gets all transitions 
for the specified charge and isotopic states for all precursors in the PrecursorDatabase. 

### Input

- `ptable::PrecursorTable` -- See `PrecursorTable`
- `transition_charges::Vector{UInt8}` -- Get `Transition`s with these charges for each `Precursor` in `ptable`
- `transition_isotopes::Vector{UInt8}` -- Get `Transition`s with these isotopes for each Precursor` in `ptable`
- `b_start::Int64` -- Get transitions bn+x and higher. At `3` would exclude b2 and b1 ions
- `y_start::Int64` -- Same as `b_start` but for y ions
- `fragment_match_ppm::Float32` -- Fragment match tolerance for the precursors

### Output
- Fills the `prec_id_to_transitions` field in a `PrecursorDatabase` 

### Examples 

"""
function addTransitions!(ptable::PrecursorDatabase, transition_charges::Vector{UInt8}, transition_isotopes::Vector{UInt8},b_start::Int64,y_start::Int64, fragment_match_ppm::Float32)
    for (prec_id, precursor) in pairs(getIDToPrec(ptable))
        for charge in transition_charges, isotope in transition_isotopes
            if !hasTransitions(ptable, prec_id)
                insert!(getPrecIDToTransitions(ptable),
                        prec_id,
                        getTransitions(precursor, 
                                        charge = charge, 
                                        isotope = isotope, 
                                        b_start = b_start, 
                                        y_start = y_start,
                                        ppm = fragment_match_ppm)
                        )
            else
                append!(getTransitions(ptable, prec_id), 
                        getTransitions(precursor, 
                                        charge = charge, 
                                        isotope = isotope, 
                                        b_start = b_start, 
                                        y_start = y_start,
                                        ppm = fragment_match_ppm)
                        )
            end
        end
    end
end
