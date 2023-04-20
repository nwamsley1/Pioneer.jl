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
function buildPrecursorTable!(ptable::PrecursorDatabase,
                              fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                              var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                              n::Int, f_path::String)
    open(f_path) do f #"./data/NRF2_SIL.txt"

        max_pepGroup_id, max_pep_id, max_prot_id = UInt32(1), UInt32(1), UInt32(1)
        timetaken = @elapsed for (row, protein_peptide) in enumerate(eachline(f))

            line = map(string, split(protein_peptide, "\t"))
            protein, peptide = line[1:2]; #Parse input "PROTEIN_NAME\tPEPTIDE_SEQUENCE"

            if !containsProt(ptable, protein) #If the protien hasn't been encountered, then add it to the hash table
                addProtein!(ptable, protein, max_prot_id);
                max_prot_id += UInt32(1);  
            end
            if !containsPepGroup(ptable, peptide) #If the peptide group hasn't been encountered
                initial_pep_id = max_pep_id
                max_pep_id = applyMods!(ptable.id_to_pep,
                                    var_mods,              #and lastly, apply variable mods and ad them to the peptide hash table
                                    fixedMods(peptide, fixed_mods),               #and increase the pep_id for each variable mod applied 
                                    max_pepGroup_id,
                                    max_pep_id,
                                    n = n); 
                pep_ids = Set([UInt32(id) for id in initial_pep_id:(max_pep_id - 1)])

                addPepGroup!(ptable, peptide, max_pepGroup_id, pep_ids, protein); #Add peptide group to the hash table
                addPepGroupToProtein!(ptable, #Add the peptide group to the protein 
                                        getProtIDFromName(ptable, protein), 
                                        max_pepGroup_id);

                max_pepGroup_id += UInt32(1); 
            else #If this peptide has been encountered before, we don't need to apply the variable modes. Instead,
                prot_id = getProtIDFromName(ptable, protein)
                pepGroup_id = getPepGroupIDFromSequence(ptable, peptide)
                addProteinToPepGroup!(ptable, 
                                        prot_id,
                                        pepGroup_id); #Add the current protein to this peptide group
                addPepGroupToProtein!(ptable, 
                                        prot_id,
                                        pepGroup_id); #Add the peptide group to the current protein 
            end
        end
        println("Time to build precursor table ", timetaken);
    end
end
