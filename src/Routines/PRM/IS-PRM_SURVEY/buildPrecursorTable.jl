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
            line = map(string, split(protein_peptide, "\t"))
            protein, peptide = line[1:2]; #Parse input "PROTEIN_NAME\tPEPTIDE_SEQUENCE"
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
                if length(line)>2
                    transition_names = line[3:end]
                    insert!(getPepIDToTransitions(ptable),
                            pep_id,
                            transition_names)
                end
                pepGroup_id += UInt32(1); 
            else #If this peptide has been encountered before, we don't need to apply the variable modes. Instead,
                addProteinToPepGroup!(ptable, protein, peptide); #Add the current protein to this peptide group
                addPepGroupToProtein!(ptable, protein, peptide); #Add the current peptide group to this protein
            end
        end
        println("Time to build precursor table ", timetaken);
    end
end

"""
    addPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})

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
function addPrecursors!(ptable::PrecursorTable, charges::Vector{UInt8}, isotopes::Vector{UInt8}, mods_dict::Dict{String, Float32})
    prec_id = UInt32(1)
    for (pep_id, peptide) in pairs(ptable.id_to_pep)
        for charge in charges, isotope in isotopes
            push!(ptable.precursors, Precursor(peptide.sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict))
            insert!(ptable.simple_precursors, prec_id, SimplePrecursor(getMZ(ptable.precursors[end]), charge, isotope, pep_id));
            prec_id+= UInt32(1)
        end
    end
    sort!(ptable.precursors, by=p->getMZ(p))
end
