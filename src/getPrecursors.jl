using Dictionaries

struct Peptide
    proteins::Set{String}
    precursors::Vector{UInt32}
end

Peptide() = Peptide(String[], UInt32[])

struct PrecursorTable
    peptides::UnorderedDictionary{String, Peptide}
    PepID_to_String::UnorderedDictionary{UInt32, String}
    precursors::Vector{Precursor}
end

PrecursorTable() = PrecursorTable(
                                    UnorderedDictionary{String, Peptide}(),
                                    UnorderedDictionary{UInt32, String}(),
                                    Vector{Precursor}())

function getPrecursors(file::String)

#readdlm("./data/NRF2_SIL.txt", '\t', String)
    readdlm(file, '\t', String)
end

new_table = PrecursorTable()

totaltime, totallines = open("./data/NRF2_SIL.txt") do f
    linecounter = 0
    pep_id = UInt32(1)
    prec_id = UInt32(1)
    timetaken = @elapsed for l in enumerate(eachline(f))
        linecounter += UInt32(1)
        #println(linecounter)
        protein, peptide = [string(x) for x in split(l[2], "\t")];
        #println(protein)
        #println(peptide)
        if !isassigned(new_table.peptides, peptide)
            #println("test")
            #Function add peptide
            insert!(new_table.PepID_to_String, pep_id, peptide)
            #println("test2")
            insert!(new_table.peptides, 
                    peptide, #Peptide Name
                    Peptide(Set([protein]), UInt32[]))
            #println("test3")
            pep_id += UInt32(1)
        end
        if protein ∉ new_table.peptides[peptide].proteins
            push!(new_table.peptides[peptide].proteins, protein)
        end
       #getPrecursors!(new_table, peptide, pep_id, prec_id)
        for charge in UInt8[1]
            for isotope in UInt8[0]
                push!(new_table.peptides[peptide].precursors, prec_id)
                #insert!(new_table.precursors,
                #        prec_id,
                push!(new_table.precursors,
                        Precursor(peptide, charge = charge, isotope = isotope, 
                        pep_id = pep_id, prec_id = prec_id))
                prec_id += UInt32(1)
            end
        end

        #low, high = split(l[2], "\t")
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
regexes = [r"E", r"(K$)"]
replacements = ["[+10]", "[+8]"]

# Define the input string
input_string = "PEPTIKDEK"

# Get all matches for each regular expression
function getMatches(regexes::Vector{Regex}, replacements::Vector{String})
    matches = []
    for (regex, replacement) in zip(regexes, replacements)
        regex_matches = findall(regex, input_string)
        for match in regex_matches
            push!(matches, (match, replacement))
        end
    end
    matches
end

function getMods(all_combinations, input_string::String)
    out = Vector{String}(undef, length(all_combinations))
    # Apply the replacements to the input string for each combination
    for (i, combination) in enumerate(all_combinations)
        output_str = [""]#nput_string[1:combination[1][1][1]]
        index = 1
        for mod in combination
            push!(output_str, input_string[index:mod[1][1]]*mod[2])
            index = mod[1][1]
        end
        push!(output_str, input_string[index:end])
        #append!(out, [join(output_str)])
        out[i] = join(output_str)
    end
    out
end


function GetMods(regexes::Vector{Regex}, replacements::Vector{String}, input_string::String)
    getMods(collect(combinations(getMatches(regexes, replacements))),
    input_string)
end
