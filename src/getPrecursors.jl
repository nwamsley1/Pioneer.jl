include("./precursor.jl");
using Combinatorics, Dictionaries

export buildPrecursorTable!, getPrecursors!

struct SimplePrecursor
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end

"""
    PeptideGroup

Type that represents an unmodified peptide 

### Fields

- prot_ids::Set{UInt32} -- Set of identifiers for proteins that contain the peptide group. 
- sequence::String -- Sequence of the peptide group

### Examples

- `PeptideGroup() = PeptideGroup(Set{UInt32}(), "")` -- constructor for an placeholder 

### GetterMethods

- getSeq(p::PeptideGroup) = p.sequence
- getProtIDs(p::PeptideGroup) = p.prot_ids

### Methods

- addProtID!(pg::PeptideGroup, prot_id::UInt32)
"""
struct PeptideGroup
    prot_ids::Set{UInt32}
    sequence::String
end

PeptideGroup() = PeptideGroup(Set{UInt32}(), "")
getSeq(p::PeptideGroup) = p.sequence
getProtIDs(p::PeptideGroup) = p.prot_ids
export getProtIds

function addProtID!(pg::PeptideGroup, prot_id::UInt32)
    push!(pg.prot_ids, prot_id)
end

"""
    Peptide

Type that represents a peptide

### Fields

- sequence::String -- Sequence of the peptide
- pep_group_id::UInt32-- Identifier of the `PeptideGroup` to which the peptide belongs. 

### GetterMethods

- getSeq(p::Peptide) = p.sequence
- getGroupID(p::Peptide) = p.pep_group_id

### Methods

- addProtID!(pg::PeptideGroup, prot_id::UInt32)
"""
struct Peptide
    sequence::String
    pep_group_id::UInt32
end

getSeq(p::Peptide) = p.sequence
getGroupID(p::Peptide) = p.pep_group_id

"""
    Protein

Type that represents a protein

### Fields

- sequence::String -- Sequence of the peptide
- pep_group_id::UInt32-- Identifier of the `PeptideGroup` to which the peptide belongs. 

### Examples

- `Protein(name::String) = Protein(name, Set{UInt32}())` -- constructor for an placeholder 

### GetterMethods

- getPepGroupIDs(p::Protein) = p.pep_group_ids
- getName(p::Protein) = p.name

### Methods

- addPepGroup!(p::Protein, pep_group_id::UInt32)
"""
struct Protein
    name::String
    pep_group_ids::Set{UInt32}
end

Protein(name::String) = Protein(name, Set{UInt32}())

getPepGroupIDs(p::Protein) = p.pep_group_ids
getName(p::Protein) = p.name

addPepGroup!(p::Protein, pep_group_id::UInt32) = push!(p.pep_group_ids, pep_group_id)

export getPepGroupIDs, getName

"""
    PrecursorTable

Data Structure that represents relations between precursors, peptides, peptide groups, and proteins

### Fields

- id_to_prot::UnorderedDictionary{UInt32, Protein} -- Maps from a protien identifier to a Protein
- prot_to_id::UnorderedDictionary{String, UInt32} -- Maps from a protein name to a protein identifier
- id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup} -- Maps a PeptideGroup identifier to a PeptideGroup
- pepGroup_to_id::UnorderedDictionary{String, UInt32} -- Maps a PeptideGroup name/sequence to an identifier
- id_to_pep::UnorderedDictionary{UInt32, Peptide} -- Maps a peptide identifier to a peptide
- precursors::Vector{Precursor} -- Precursor has fields `pep_id` and `prec_id`. `pep_id`s are keys for `id_to_pep`

### Examples

- `PrecursorTable() = PrecursorTable(
    UnorderedDictionary{UInt32, Protein}(),
    UnorderedDictionary{String, UInt32}(),
    UnorderedDictionary{UInt32, PeptideGroup}(),
    UnorderedDictionary{String, UInt32}(),
    UnorderedDictionary{UInt32, Peptide}(),
    Vector{Precursor}())` -- constructor for a placeholder 

### GetterMethods

- getIDToProt(p::PrecursorTable)
- getProtToID(p::PrecursorTable) 
- getIDToPepGroup(p::PrecursorTable) 
- getPepGroupToID(p::PrecursorTable)
- getIDToPep(p::PrecursorTable)
- getPrecursors(p::PrecursorTable)
- getProtID(p::PrecursorTable, protein::String)
- getProt(p::PrecursorTable, prot_id::UInt32)
- getPepGroupID(p::PrecursorTable, peptide::String)
- getPepGroup(p::PrecursorTable, pepGroup_id::UInt32)
- getPep(p::PrecursorTable, pep_id::UInt32)
- getProtNamesFromPepSeq(p::PrecursorTable, peptide::String)
- getPepGroupsFromProt(p::PrecursorTable, protein::String)
- getPepSeqsFromProt(p::PrecursorTable, protein::String)
- getPepGroupsFromProt(p::PrecursorTable, prot_id::UInt32)
- getPepSeqsFromProt(p::PrecursorTable, prot_id::UInt32)

### Methods

- insertProtID!(p::PrecursorTable, protein::String, prot_id::UInt32)
- insertProt!(p::PrecursorTable, protein::String, prot_id::UInt32)
- insertPepGroupID!(p::PrecursorTable, peptide::String, pepGroup_id::UInt32)
- insertPepGroup!(p::PrecursorTable, protein::String, peptide::String, pepGroup_id::UInt32)

- addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)
- addPepGroupToProtein!(p::PrecursorTable, protein::String, peptide::String)
- addNewProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)
- addNewPeptideGroup!(peptide::String, pepGroup_id::UInt32, protein::String, precursor_table::PrecursorTable)

### Notes

`PrecursorTable` keeps an account of relations between Proteins, PeptideGroups, Peptides, and
Precursors. Each Protein, PeptideGroup, Peptide, and Precursor has a unique UInt32 indentifier. 
Each also has a string name. See `Protein`, `PeptideGroup`, `Peptide` and `Precursor`. 
Proteins map to one or more PeptideGroups. PeptideGroups map to one or more Peptides, 
and each Peptide maps to one or more precursors.
"""
mutable struct PrecursorTable
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptide group
    simple_precursors::UnorderedDictionary{UInt32, SimplePrecursor}
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
                                    UnorderedDictionary{UInt32, SimplePrecursor}(),
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
getSimplePrecFromID(p::PrecursorTable, prec_id::UInt32) = p.simple_precursors[prec_id]
getProtNamesFromPepSeq(p::PrecursorTable, peptide::String) = Set([getName(getProt(p, prot_id)) for prot_id in getProtIDs(getPepGroup(p, getPepGroupID(p, peptide)))])
getPepGroupsFromProt(p::PrecursorTable, protein::String) = [getPepGroup(p, ID) for ID in getPepGroupIDs(getProt(p, getProtID(p, protein)))]
getPepSeqsFromProt(p::PrecursorTable, protein::String) = [getSeq(pep_group) for pep_group in getPepGroupsFromProt(p, protein)]
getPepGroupsFromProt(p::PrecursorTable, prot_id::UInt32) = [getPepGroup(p, ID) for ID in getPepGroupIDs(getProt(p, prot_id))]
getPepSeqsFromProt(p::PrecursorTable, prot_id::UInt32) = [getSeq(pep_group) for pep_group in getPepGroupsFromProt(p, prot_id)]
export getProtID, getProt, getPepGroupID, getPepGroup, getPep

insertProtID!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.prot_to_id, protein, prot_id)
insertProt!(p::PrecursorTable, protein::String, prot_id::UInt32) = insert!(p.id_to_prot, prot_id, Protein(protein))
insertPepGroupID!(p::PrecursorTable, peptide::String, pepGroup_id::UInt32) = insert!(p.pepGroup_to_id, peptide, pepGroup_id)
insertPepGroup!(p::PrecursorTable, protein::String, peptide::String, pepGroup_id::UInt32) = insert!(p.id_to_pepGroup, pepGroup_id, PeptideGroup(Set(getProtID(p, protein)), peptide))

"""
    addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)

For the PeptideGroup to which "peptide" belongs, add the protein to its `prot_ids`
"""
function addProteinToPepGroup!(p::PrecursorTable, protein::String, peptide::String)
    addProtID!(getPepGroup(p, getPepGroupID(p, peptide)), getProtID(p, protein))
end

"""
    addPepGroupToProtein!(p::PrecursorTable, protein::String, peptide::String)

For the protein to which the `peptide` belongs, add the correct peptide group id to its `pep_group_ids`
"""
function addPepGroupToProtein!(p::PrecursorTable, protein::String, peptide::String)
    addPepGroup!(getProt(p, getProtID(p, protein)), getPepGroupID(p, peptide))
end

"""
    addNewProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)

Create a new `Protein` and protein id and add them to the lookup tables
"""
function addNewProtein!(protein::String, prot_id::UInt32, precursor_table::PrecursorTable)
    insertProtID!(precursor_table, protein, prot_id);
    insertProt!(precursor_table, protein, prot_id);
end

"""
    addNewPeptideGroup!(peptide::String, pepGroup_id::UInt32, protein::String, precursor_table::PrecursorTable)

Create a new `PeptideGroup` and peptide group id and add them to the lookup tables
"""
function addNewPeptideGroup!(peptide::String, pepGroup_id::UInt32, protein::String, precursor_table::PrecursorTable)
        insertPepGroupID!(precursor_table, peptide, pepGroup_id);
        insertPepGroup!(precursor_table, protein, peptide, pepGroup_id)
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

"""
    getPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})

Modifies a `PrecursorTable` to fill the `precursors` field. Gets a `Precursor` for each unique `Peptide` in `ptable`. Sorts the 
precursors in order of increasing mz.   

### Input

- `ptable::PrecursorTable` -- See `PrecursorTable`
- `charges::Vector{UInt8}` -- Get a precursor with these charges for each Peptide in `ptable`
- `isotopes::Vector{UInt8}` -- Get a precursor with these isotopes for each Peptide in `ptable`
- `mods_dict::Dict{String, Float32}` -- Dict mapping modifications to their masses

### Output
- Fills the `precursors` field in a `ptable` 

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
function getPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})
    prec_id = UInt32(1)
    for (pep_id, peptide) in pairs(ptable.id_to_pep)
        for charge in charges, isotope in isotopes
            push!(ptable.precursors, Precursor(peptide.sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict))
            insert!(ptable.simple_precursors, prec_id, SimplePrecursor(charge, isotope, pep_id));
            prec_id+= UInt32(1)
        end
    end
    sort!(ptable.precursors, by=p->getMZ(p))
end

getPepIDFromPrecID(p::PrecursorTable, prec_id::UInt32) = getSimplePrecFromID(p, prec_id).pep_id
"""
    buildPrecursorTable!(ptable::PrecursorTable, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, n::Int, f_path::String)

Builds a `PrecursorTable` given a path to a tab delimited text file where each line has a protein name in the first field and a peptide sequence in the second. 

### Input

- `ptable::PrecursorTable` -- See `PrecursorTable`
- `fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Specifies fixed modifications to apply to each peptide
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}` -- Specifies variable modifications to apply to each peptide
- `n::Int` -- Apply all combinations of `n` or fewer variable modifications to each peptide. 
- `f_path::String` -- Path to a tab delimited text file where each line has a protein name in the first field and a peptide sequence in the second. 

### Output
- Fills the `ptable` given the protein_name => peptide_sequence pairs. 

### Examples 
    testPtable = PrecursorTable()
    fixed_mods = [(p=r"C", r="C[Carb]")]
    var_mods = [(p=r"(K\$)", r="[Hlys]"), (p=r"(R\$)", r="[Harg]")]
    buildPrecursorTable!(testPtable, fixed_mods, var_mods, 2, "../data/NRF2_SIL.txt")
    getPrecursors!(testPtable, UInt8[2, 3, 4], UInt8[0], test_mods)
    @test length(getIDToPepGroup(testPtable)) == 260
    @test length(getPrecursors(testPtable)) == 260*2*3 #2 because each protein ends in a variably modifiable K or R. 3 because 3 charge states. 

    @test getProtNamesFromPepSeq(testPtable, "LAIEAGFR") == Set(["AKR1C1","AKR1C2","AKR1C3"]) 
    @test Set(getPepSeqsFromProt(testPtable, "CD8A")) == Set(["AAEGLDTQR", "TWNLGETVELK", "LGDTFVLTLSDFR"])
    @test Set(getPepSeqsFromProt(testPtable, "ZNF746")) == Set(["LLSLEGR", "GQPLPTPPAPPDPFK"])
"""
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
                pep_id = applyMods!(ptable.id_to_pep,
                                    var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                                    peptide,               #and increase the pep_id for each variable mod applied 
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
