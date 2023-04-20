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

Applies all combinations of n or fewer variable modifications to an input peptide sequence, `unmod_seq`. Two outputs:
- 1) Increments the supplied `pep_id` for each unique combination of mods and returns the final `pep_id`. 
- 2) Adds each unique modified peptide to the `peptides` dict. Does this inplace. 

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
    # #Apply all combinations of n or fewer variable mods
    for N in 1:min(n, length(matches)) #min(n, length(matches)) because there could be fewer than `n` matches
        for combination in combinations(matches, N) #Each combination of "N" variable mods
            #Maybe the sorting could somehow ber done all at once outside the loop. Could be more efficient?
            sort!(combination, by=match->match[1][end]); #Sort applied variable mods in order of appearance in "unmod_seq". 
            modified_seq = applyMods(combination, unmod_seq) #Get the modified sequence for the given combination of mods
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
