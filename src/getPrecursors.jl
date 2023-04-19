
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
    matchVarMods(patterns::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String)::Vector{Tuple{UnitRange{Int64}, String}}

Finds each match to patterns in the input_string. Returns a Vector{Tuple{UnitRange{Int64}, String}}. Each pattern in `patterns` includes a 
regular expresion `:p` and a replacement `:r`, which should substitute the matched pattern. For the output, the first element of each tuple is a UnitRange
corresponding to the indices of `input_string` that matched a pattern. The second element is a string that should be substituted at those indices.  

### Input

- `patterns::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: -- Patterns to match and their substitutions (could replace r"K\$" with "K[C-term K mod]")
- `input_string::String` -- A string to search for matches to `patterns`

### Output
- `::Vector{Tuple{UnitRange{Int64}, String}}`

Each UnitRange{Int64} corresponds to a range in `input_string` that matched a pattern. Each String is intended to replace those indices
in the `input_string`. 

### Algorithm 
- Searches the input string for each pattern seperately. Appends each match to the output as a Tuple{UnitRange{Int64}, String}. 

### Notes

- Different combinations of each modification need to be applied. That is why this function returns a Vector of modifications
rather than applying each modification to the `input_string` and returning that as output. 

### Examples 
var_mods = [(p=r"(E)", r="[modE]"), (p=r"(C)", r="[modC]")]
test_pep = "PEPEPCPEPEP"
@test Set(matchVarMods(var_mods, test_pep)) == Set([(2:2, "[modE]"),
                                                (4:4, "[modE]"),
                                                (6:6, "[modC]"),
                                                (8:8, "[modE]"),
                                                (10:10, "[modE]")])
"""
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

"""
    applyVariableMods!(peptides::UnorderedDictionary{UInt32, Peptide}, matches::Vector{Tuple{UnitRange{Int64}, String}}, unmod_seq::String, group_id::UInt32, pep_id::UInt32, n::Int)

Applies all combinations of n or fewer variable modifications to an input peptide sequence, `unmod_seq`. 

### Input

- `peptides::UnorderedDictionary{UInt32, Peptide}` -- Maps unique peptide identifiers to `Peptide`s
- `matches::Vector{Tuple{UnitRange{Int64}, String}}` -- Tuples of substitutions `String` to be substituted into `unmod_seq` at indices `UnitRange{Int64}`
- `unmod_seq::String` -- A string representing a peptide sequence
- `group_id::UInt32` -- Unique identifier for a `PeptideGroup`
- `pep_id::UInt32` -- Unique identifier for a `Peptide`
- `n::Int` -- Apply all combinations of n or fewer modifications (specified in `matches`). 

### Output
- Modifies `peptides`. Adds a new key value pair UInt32=>`Peptide` for each unique peptide modification. Also
adds the unmodified peptide. 

### Algorithm 
- Relies on `Combinatorics` package to find all combinations of n or fewer `matches` that can be applied to the peptide sequence `unmod_seq`. 

### Notes

- If there are only M matches but `n`>length(`matches`), will apply at most M matches
- If multiple mods can be applied to the same index, they will be appended. For example: "DRAGRACER[R mod 1][R mod 2]". 

### Examples 
- See `applyMods!(peptides::UnorderedDictionary{UInt32, Peptide}, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String, group_id::UInt32, pep_id::UInt32; n::Int = 3)`
"""
function applyVariableMods!(peptides::UnorderedDictionary{UInt32, Peptide}, matches::Vector{Tuple{UnitRange{Int64}, String}}, unmod_seq::String, group_id::UInt32, pep_id::UInt32, n::Int)
    
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
                   pep_id::UInt32) = insert!(peptides, pep_id, Peptide(seq, group_id, Set(UInt32[])))

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

"""
    applyMods!(peptides::UnorderedDictionary{UInt32, Peptide}, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, unmod_seq::String, group_id::UInt32, pep_id::UInt32; n::Int = 3)

Wrapper for `applyVariableMods!`. Applies all combinations of n or fewer variable modifications, `var_mods`, to an input peptide sequence, `unmod_seq`. 

### Input

- `peptides::UnorderedDictionary{UInt32, Peptide}` -- Maps unique peptide identifiers to `Peptide`s
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Pairs of match patterns and substitutions for variable modifications
- `unmod_seq::String` -- A string representing a peptide sequence
- `group_id::UInt32` -- Unique identifier for a `PeptideGroup`
- `pep_id::UInt32` -- Unique identifier for a `Peptide`
- `n::Int` -- Apply all combinations of n or fewer modifications (specified in `matches`). 

### Output
- Modifies `peptides`. Adds a new key value pair UInt32=>`Peptide` for each unique peptide modification. Also
adds the unmodified peptide. 

### Algorithm 
- Relies on `Combinatorics` package to find all combinations of n or fewer `matches` that can be applied to the peptide sequence `unmod_seq`. 

### Notes

- If there are only M matches but `n`>length(`matches`), will apply at most M matches
- If multiple mods can be applied to the same index, they will be appended. For example: "DRAGRACER[R mod 1][R mod 2]". 

### Examples 

test_peps_dict = UnorderedDictionary{UInt32, Peptide}()
group_id = UInt32(1)
pep_id = UInt32(1)
var_mods = [(p=r"(E)", r="[modE]"), (p=r"(C)", r="[modC]")]
test_pep = "PEPEPCPEPEP"
applyMods!(test_peps_dict, var_mods, test_pep, group_id, pep_id, n=1);
@test Set(map(pep->getSeq(pep), test_peps_dict)) == Set([
                                                            "PEPEPCPEPEP",
                                                            "PEPEPC[modC]PEPEP",
                                                            "PE[modE]PEPCPEPEP",
                                                            "PEPE[modE]PCPEPEP",
                                                            "PEPEPCPE[modE]PEP",
                                                            "PEPEPCPEPE[modE]P"
                                                            ])
"""
function applyMods!(peptides::UnorderedDictionary{UInt32, Peptide}, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, unmod_seq::String, group_id::UInt32, pep_id::UInt32; n::Int = 3)
    applyVariableMods!(peptides,
                    matchVarMods(var_mods, unmod_seq),
                    unmod_seq,
                    group_id,
                    pep_id,
                    n
                    )
end
export applyMods!

"""
    fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})

Applies fixed moficitations to a peptide sequence and returns the modified sequence.  

### Input

- `peptide::String` -- A peptide sequence. 
- `Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Pairs of match match patterns and substitutions for fixed modifications

### Output
- A peptide sequence with fixed modifications applied. 

### Algorithm 

for mod in fixed_mods
    peptide = replace(peptide, mod[:p]=>mod[:r])
end

### Examples 

fixed_mods = [(p=r"C", r="C[Carb]")]
@test fixedMods("CARACHTER", fixed_mods) == "C[Carb]ARAC[Carb]HTER"
"""
function fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    for mod in fixed_mods
        peptide = replace(peptide, mod[:p]=>mod[:r])
    end
    peptide
end

export fixedMods



#=
"""
    getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, n::Int, f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},mass_mods::Dict{String, Float32})

Build a `PrecursorTable` given a file_path to a tab delimited table of protein_name peptide_sequence pairs. Applies fixed and variable modifications. Gets precursors
for each combination of `charges` and `isotopes` supplied. 

### Input

- `fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Specification of fixed modifications to apply
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}` -- Specification of variable modifications to apply
- `n::Int` -- Will apply all combinations of `n` or fewer variable modifications to each sequence
- `f_path::String` -- File path fo the protein-peptide list. 
- `charges::Vector{UInt8}` -- For each peptide, gets precursors with these charge states 
- `isotopes::Vector{UInt8}` -- For each peptide, gets precursors with these isotopic states
- ` mass_mods::Dict{String, Float32}` -- Specifies mass for each modification in `fixed_mods` and `var_mods` 

### Output
- Returns a `PrecursorTable` struct. See documentation for `PrecursorTable`. Mapps identifiers between precursors,
pepties, peptide groups, and proteins. 

### Notes

### Examples 

"""
function getPrecursors(fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        n::Int, 
        f_path::String, charges::Vector{UInt8}, isotopes::Vector{UInt8},
        mass_mods::Dict{String, Float32})
    ptable = PrecursorTable()
    buildPrecursorTable!(ptable, fixed_mods, var_mods, n, f_path)
    addPrecursors!(ptable, charges, isotopes, mass_mods)
    return ptable
end
=#


