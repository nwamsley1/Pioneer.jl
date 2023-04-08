struct LightHeavyPair
    prot_ids::Set{UInt32}
    sequence::String
    light_pep_id::UInt32
    heavy_pep_id::UInt32
    light_scan_idxs::Vector{Int64}
    heavy_scan_idxs::Vector{Int64}
    scan_idxs::Set{Int64}
end

getProtIDs(p::LightHeavyPair) = p.prot_ids
getSeq(p::LightHeavyPair) = p.sequence
getLightID(p::LightHeavyPair) = p.light_pep_id
getLightScanIDs(p::LightHeavyPair) = p.light_scan_idxs
getHeavyID(p::LightHeavyPair) = p.heavy_pep_id
getHeavyScanIDs(p::LightHeavyPair) = p.heavy_scan_idxs
getScanIDs(p::LightHeavyPair) = p.scan_idxs

#addProtID(p::LightHeavyPair, prot_id::UInt32) = push!(getProtIDs(p), prot_id)

mutable struct ISPRMPrecursorTable
    ptable::PrecursorTable 
    id_to_light_heavy_pair::UnorderedDictionary{UInt32, LightHeavyPair}
    pep_id_to_light_heavy_pair_id::UnorderedDictionary{UInt32, UInt32}
    pep_id_to_transitions::UnorderedDictionary{UInt32, Vector{Transition}}
    precursors::Dictionary{UInt32, Precursor}
end

getPTable(p::ISPRMPrecursorTable) = p.ptable
getPrecursors(p::ISPRMPrecursorTable) = p.precursors
getIDToLightHeavyPair(p::ISPRMPrecursorTable) = p.id_to_light_heavy_pair
getPepIDToLightHeavyPairID(p::ISPRMPrecursorTable) = p.pep_id_to_light_heavy_pair_id
getPepIDToTransitions(p::ISPRMPrecursorTable) = p.pep_id_to_transitions
getLightHeavyPairFromPepID(p::ISPRMPrecursorTable, pep_id::UInt32) = getIDToLightHeavyPair(p)[getPepIDToLightHeavyPairID(p)[pep_id]]



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
function buildPrecursorTable!(ptable::ISPRMPrecursorTable, mods_dict::Dict{String, Float32}, f_path::String)

    function readheader(f::IOStream)
        Dict(str => ind for (ind, str) in enumerate(map(string, split(readline(f), "\t"))))
    end
    function parseTransitoinName(transition_name::String)
        return (ion_type = transition_name[1], ind = split(transition_name,"+")[1][2:end], charge = transition_name[end])
    end
    open(f_path) do f #"./data/NRF2_SIL.txt"

        header = readheader(f)

        pepGroup_id, pep_id, prot_id, light_heavy_pair_id = UInt32(1), UInt32(1), UInt32(1), UInt32(1)

        timetaken = @elapsed for (row, protein_peptide) in enumerate(eachline(f))

            line = map(string, split(protein_peptide, "\t"))

            protein = line[header["protein_name"]]
            #Add a check to make sure it is a heavy peptide to begin with. Warn and skip if not. 
            heavy_peptide = line[header["sequence"]]
            light_peptide = replace(heavy_peptide, r"\[H[a-z]{3}\]" =>"")
            unmodified_peptide = replace(heavy_peptide, r"\[[^\]]*\]")
            
            #if !isassigned(getPepGroupToID(ptable), peptide)
            #    insert!(getPepGroupToID(ptable), peptide, )
            #end
            #push!(ptable.precursors, Precursor(peptide.sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict))
            #insert!(ptable.simple_precursors, prec_id, SimplePrecursor(getMZ(ptable.precursors[end]), charge, isotope, pep_id));
            if !containsProt(ptable, protein) #If the protien hasn't been encountered,
                #then add it to the hash table
                addNewProtein!(protein, prot_id, ptable);
                prot_id += UInt32(1);  
            end
            #replace("PEPTIDEK[Hlys]", testregex =>"")
            if !containsPepGroup(ptable, light_peptide) #If the peptide hasn't been encountered, then
                addNewPeptideGroup!(unmodified_peptide, pepGroup_id, protein, ptable); #add it to the hash table,
                addPepGroupToProtein!(ptable, protein, unmodified_peptide); #add a new peptide group to the protein,
                
                insert!(getIDToPep(ptable), pep_id, Peptide(light_peptide, pepGroup_id))
                insert!(getIDToPep(ptable), pep_id + 1, Peptide(light_peptide, pepGroup_id))

                #Modify incase multiple protin ids. 
                insert!(getIDToLightHeavyPair(ptable),
                        light_heavy_pair_id,
                        LightHeavyPair(Set(prot_id), light_peptide, pep_id, pep_id + 1, Int64[], Int64[],  Set([Int64[]]))
                    )

                insert!(getIDtoPep(getPTable(ptable)),
                        pep_id,
                        Peptide(light_peptide, pepGroup_id))

                insert!(getIDtoPep(getPTable(ptable)),
                        pep_id + 1,
                        Peptide(heavy_peptide, pepGroup_id))

                insert!(getPepIDToLightHeavyPairID(ptable),
                        pep_id, 
                        light_heavy_pair_id)

                insert!(getPepIDToLightHeavyPairID(ptable),
                        pep_id + 1, 
                        light_heavy_pair_id)                        

                #Need to be able to deal with lines that are the same excepting
                #Charge and isotope. 
                insert!(getPrecursors(ptable),
                        prec_id, 
                        Precursor(light_peptide, 
                                    charge = UInt8(line[header["charge"]]),
                                    isotope = UInt8(line[header["isotope"]]),
                                    pep_id = pep_id, 
                                    prec_id = prec_id,
                                    mods_dict)
                            )
                insert!(getPrecursors(ptable),
                        prec_id + 1, 
                        Precursor(heavy_peptide, 
                                    charge = UInt8(line[header["charge"]]),
                                    isotope = UInt8(line[header["isotope"]]),
                                    pep_id = pep_id + 1, 
                                    prec_id = prec_id + 1,
                                    mods_dict)
                )

                for transition_name in line[header["transition_names"]:end]
                    transition = parseTransitionName(transition_name)
                    #Light peptide transitions
                    insert!(getPepIDToTransitions(ptable),
                            pep_id,
                            Transition(getResidues(getPrecursors(ptable)[prec_id]),
                                        charge = transition[:charge],
                                        isotope = UInt8(0),
                                        prec_id = prec_id)
                                )
                    #Heavy peptide transitions
                    insert!(getPepIDToTransitions(ptable),
                            pep_id + 1,
                            Transition(getResidues(getPrecursors(ptable)[prec_id + 1]),
                                        charge = transition[:charge],
                                        isotope = UInt8(0),
                                        prec_id = prec_id + !)
                                )

                end
                struct LightHeavyPair
                    prot_ids::Set{UInt32}
                    sequence::String
                    light_pep_id::UInt32
                    heavy_pep_id::UInt32
                    light_scan_idxs::Vector{Int64}
                    heavy_scan_idxs::Vector{Int64}
                    scan_idxs::Set{Int64}
                end
                mutable struct ISPRMPrecursorTable
                    ptable::PrecursorTable 
                    id_to_light_heavy_pair::UnorderedDictionary{UInt32, LightHeavyPair}
                    pep_id_to_light_heavy_pair_id::UnorderedDictionary{UInt32, UInt32}
                end
                if length(line)>2
                    transition_names = line[3:end]
                    insert!(getPepIDToTransitions(ptable),
                            pep_id,
                            transition_names)
                end
                pepGroup_id += UInt32(1); 
            else #If this peptide has been encountered before, just update the protein and peptide groups
                addProteinToPepGroup!(ptable, protein, unmodified_peptide); #Add the current protein to this peptide group
                addPepGroupToProtein!(ptable, protein, unmodified_peptide); #Add the current peptide group to this protein
            end
        end
        println("Time to build precursor table ", timetaken);
    end
end
