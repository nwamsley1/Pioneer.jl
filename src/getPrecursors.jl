using Dictionaries

struct Peptide
    sequence::String
    prot_ids::Set{UInt32}
    pep_group_id::UInt32
end

Peptide() = Peptide("", Set{UInt32}(), UInt32(0))

struct PrecursorTable
    prot_id_string::UnorderedDictionary{UInt32, String}
    pep_group_id_string::UnorderedDictionary{UInt32, String}
    pep_id::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptides
    #pred_id_prec::Dictionary{UInt32, Precursor} #Vector of precursors
    peptide_groups::Set{String} #Set of peptide groups (unmodified peptide sequences)
    proteins::Set{String} #Set of protein names
end

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, String}(),
                                    UnorderedDictionary{UInt32, Peptide}(),
                                    Set{String}(),
                                    Set{String}()
                                    )

function getPrecursors(file::String)

#readdlm("./data/NRF2_SIL.txt", '\t', String)
    readdlm(file, '\t', String)
end

test_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464),
    "Harg" => Float32(10),
    "Hlys" => Float32(8)
)

precursor_table = PrecursorTable()

totaltime, totallines = open("./data/NRF2_SIL.txt") do f
    linecounter = 0
    pep_group_id = UInt32(1)
    pep_id = UInt32(1)
    prot_id = UInt32(1)
    timetaken = @elapsed for l in enumerate(eachline(f))
        linecounter += UInt32(1)
        #println(linecounter)
        protein, peptide_seq = [string(x) for x in split(l[2], "\t")];
        #Apply fixed modifications
        peptide_seq = fixedMods(peptide_seq, fixed_mods)
        if !isassigned(precursor_table.prot_id, protein)
            prot_id += UInt32(1);
            addProtein!(protein, prot_id, precursor_table);
        end
        if !isassigned(precursor_table.pep_group_id, peptide_seq)
            pep_group_id += UInt32(1);
            addPeptideGroup!(peptide_seq, pep_group_id, prot_id, precursor_table)
            #Need to apply the variable mods and add peptides. 
            applyMods!(var_mods,
                       peptide_seq,
                       p.pep_id,
                       pep_group_id,
                       2)
        else
            #Checks if the current protein is assigned to the current peptide group
            #Adds the protien ID if not. 
            addProteinToPepGroup!(precursor_table, protein, peptide)
            break
        end
    end
    (timetaken, linecounter)
end

function fixedMods(peptide::String, fixed_mods::Vector{Pair{Regex, String}})
    replace("PEPTIKDEK", fixed_mods[1])
end
function variableMods(peptide::String, fixed_mods::Vector{Pair{Regex, String}})

end

totaltime, totallines = open("./data/NRF2_SIL.txt") do f
    linecounter = 0
    timetaken = @elapsed for l in enumerate(eachline(f))
        linecounter += 1
        println(linecounter)
        a, b = split(l[2], "\t")
    end
    (timetaken, linecounter)
end

using IterTools

# Define the regular expressions and replacement strings
using Combinatorics

# Define the regular expressions and replacement strings
# Define the regular expressions and replacement strings
var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
fixed_mods = [(p=r"C", r="C[Carb]")]
#replacements = ["[+11]","[+10]", "[+8]", "[+2]"]
# Define the input string
input_string = "PEPTICKDECK"
x = String[]

# Get all matches for each regular expression
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


@btime getMatches(regexes, replacements, input_string)

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
            insert!(precursors, Peptide(input_string, group_id), pep_id);
        end
    end
    pep_id += UInt32(1);
    insert!(precursors, Peptide(input_string, group_id), pep_id)
    pep_id
    #out
end
#How can we prevent double modification in some cases?
#If we want lysine formyltaion and acetylation as variable mods in the 
#same experiment, obviously those mods don't coeixt? Prabably just need to add a step that checks.  
function applyFixedMods(patterns::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string)
    for pattern in patterns
        input_string = replace(input_string, pattern[:p]=>pattern[:r])
    end
    input_string
end

function applyMods!(var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, input_string::String, peptides::UnorderedDictionary{UInt32, Peptide}, group_id::UInt32, pep_id::UInt32; n::Int = 3)
    applyVariableMods!(matchVarMods(var_mods, input_string),
                    input_string,
                    peptides,
                    group_id,
                    pep_id,
                    n
                    )
end
#What about groups of mutually exclusive variable mods?
