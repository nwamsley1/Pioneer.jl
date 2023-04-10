struct LightHeavyPeptidePair
    pepGroup_id::UInt32
    light_sequence::String
    heavy_sequence::String
    lh_prec_id_pairs::Set{Tuple{UInt32, UInt32}}
    light_scan_idxs::Vector{Int64}
    heavy_scan_idxs::Vector{Int64}
    integration_boundary_scan_idxs::Set{Int64}
end
LightHeavyPeptidePair(pepGroup_id::UInt32, light_sequence::String, heavy_sequence::String) = LightHeavyPeptidePair(pepGroup_id, light_sequence, heavy_sequence,
                                                                                                                   Set{Tuple{UInt32, UInt32}}(),Int64[], Int64[], Set(Int64[]))

getPepGroupID(p::LightHeavyPeptidePair) = p.pepGroup_id
getLightSeq(p::LightHeavyPeptidePair) = p.light_sequence
getHeavySeq(p::LightHeavyPeptidePair) = p.heavy_sequence
getPrecIDPairs(p::LightHeavyPeptidePair) = p.lh_prec_id_pairs
getLightScanIDs(p::LightHeavyPeptidePair) = p.light_scan_idxs
getHeavyScanIDs(p::LightHeavyPeptidePair) = p.light_scan_idxs
getIntegrationBoundaryScanIDs(p::LightHeavyPeptidePair) = p.integration_boundary_scan_idxs

#addProtID(p::LightHeavyPair, prot_id::UInt32) = push!(getProtIDs(p), prot_id)

mutable struct ISPRMPrecursorTable
    ptable::PrecursorTable 
    pep_id_to_light_heavy_pair::UnorderedDictionary{UInt32, LightHeavyPeptidePair}
    pep_sequence_to_pep_id::UnorderedDictionary{String, UInt32}
    prec_id_to_transitions::UnorderedDictionary{UInt32, Vector{Transition}}
    light_precursor_set::Set{SimplePrecursor}
    prec_id_to_precursors::Dictionary{UInt32, Precursor}
end

ISPRMPrecursorTable() = ISPRMPrecursorTable(PrecursorTable(), UnorderedDictionary{UInt32, LightHeavyPeptidePair}(), UnorderedDictionary{String, UInt32}(), 
UnorderedDictionary{UInt32, Vector{Transition}}(), Set{SimplePrecursor}(), Dictionary{UInt32, Precursor}())

getPTable(p::ISPRMPrecursorTable) = p.ptable
getPrecIDToLightHeavyPair(p::ISPRMPrecursorTable) = p.pep_id_to_light_heavy_pair
getPrecIDToTransitions(p::ISPRMPrecursorTable) = p.prec_id_to_transitions
getPepSequenceToPepID(p::ISPRMPrecursorTable) = p.pep_sequence_to_pep_id
getLightHeavyPairFromPepID(p::ISPRMPrecursorTable, pep_id::UInt32) = getPrecIDToLightHeavyPair(p)[pep_id]
getPrecursors(p::ISPRMPrecursorTable) = p.prec_id_to_precursors

isinPrecursorSet(p::ISPRMPrecursorTable, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32) = SimplePrecursor(mz, charge, isotope, pep_id) ∈ p.light_precursor_set
isinPrecursorSet(p::ISPRMPrecursorTable, prec::Precursor) = isinPrecursorSet(p, getMZ(prec), getCharge(prec), getIsotope(prec), getPepID(prec))
addToPrecursorSet!(p::ISPRMPrecursorTable, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32) = push!(p.light_precursor_set, SimplePrecursor(mz, charge, isotope, pep_id))
addToPrecursorSet!(p::ISPRMPrecursorTable, prec::Precursor) = push!(p.light_precursor_set, SimplePrecursor(getMZ(prec), getCharge(prec), getIsotope(prec), getPepID(prec)))
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

    function readHeader(f::IOStream)
        Dict(str => ind for (ind, str) in enumerate(map(string, split(readline(f), ","))))
    end
    
    function parseLine(row::String, header::Dict{String, Int})
        println("header ", header)
        row = map(string, split(row, ","))
        return row[header["protein_name"]], row[header["sequence"]], parse(UInt8, row[header["precursor_charge"]]), parse(UInt8, row[header["precursor_isotope"]]), row[header["transition_names"]:end]
    end

    function parseTransitionName(transition_name::String)
        return (ion_type = transition_name[1], ind = parse(UInt8, split(transition_name,"+")[1][2:end]), charge = parse(UInt8, transition_name[end]))
    end

    function getProtID!(ptable::ISPRMPrecursorTable, protein_name::String, max_prot_id::UInt32)
        if !containsProt(getPTable(ptable), protein_name) #If the protien hasn't been encountered,
            #then add it to the hash table\
            max_prot_id += UInt32(1);  
            prot_id = max_prot_id
            addNewProtein!(protein_name, prot_id, getPTable(ptable));
        else
            prot_id = getProtToID(getPTable(ptable))[protein_name]
        end
        return prot_id, max_prot_id
    end

    function getPepGroupID!(ptable::ISPRMPrecursorTable, protein_name::String, unmodified_sequence::String, max_pepGroup_id::UInt32)
        if !containsPepGroup(getPTable(ptable), unmodified_sequence) #If the peptide hasn't been encountered
            max_pepGroup_id += UInt32(1)
            pepGroup_id = max_pepGroup_id
            addNewPeptideGroup!(unmodified_sequence, pepGroup_id, protein_name, getPTable(ptable)); #add it to the hash table,
            addPepGroupToProtein!(getPTable(ptable), protein_name, unmodified_sequence); #add a new peptide group to the protein,
        else #If this peptide has been encountered before, just update the protein and peptide groups
            addProteinToPepGroup!(getPTable(ptable), protein_name, unmodified_sequence); #Add the current protein to this peptide group
            addPepGroupToProtein!(getPTable(ptable), protein_name, unmodified_sequence); #Add the current peptide group to this protein
            pepGroup_id = getPepGroupToID(getPTable(ptable))[unmodified_sequence]
        end
        return pepGroup_id, max_pepGroup_id
    end 

    function getPepID!(ptable::ISPRMPrecursorTable, light_sequence::String, heavy_sequence::String, max_pep_id::UInt32, pepGroup_id::UInt32)
        if !isassigned(getPepSequenceToPepID(ptable), light_sequence) #If the peptide hasn't been encountered
            max_pep_id += UInt32(1);  
            pep_id = max_pep_id 
            insert!(getPepSequenceToPepID(ptable), light_sequence, pep_id)
            insert!(getPrecIDToLightHeavyPair(ptable), 
                    pep_id,
                    LightHeavyPeptidePair(pepGroup_id, light_sequence, heavy_sequence))
        else
            pep_id = getPepSequenceToPepID(ptable)[light_sequence]
        end
        pep_id, max_pep_id
    end
    
    function addPrecursors!(ptable::ISPRMPrecursorTable, light_precursor::Precursor, heavy_precursor::Precursor, max_prec_id::UInt32, pep_id::UInt32)
        max_prec_id += UInt32(2)
        light_prec_id, heavy_prec_id = max_prec_id - UInt32(1), max_prec_id

        addToPrecursorSet!(ptable, light_precursor)

        insert!(getPrecursors(ptable), light_prec_id, light_precursor)
        insert!(getPrecursors(ptable), heavy_prec_id, heavy_precursor)

        push!(getPrecIDPairs(getLightHeavyPairFromPepID(ptable, pep_id)), (light_prec_id, heavy_prec_id))
        return light_prec_id, heavy_prec_id, max_prec_id
    end

    function addTransitions!(ptable::ISPRMPrecursorTable, transition_names::Vector{String}, prec_ids::Vector{UInt32})
        #for transition_name in line[header["transition_names"]:end]
        for prec_id in prec_ids 
            insert!(getPrecIDToTransitions(ptable),
                    prec_id,Vector{Transition}())
            for transition_name in transition_names
                transition = parseTransitionName(transition_name)
                push!(getPrecIDToTransitions(ptable)[prec_id],
                Transition(getResidues(getPrecursors(ptable)[prec_id]),
                            charge = transition[:charge],
                            isotope = UInt8(0),
                            prec_id = prec_id)
                    )
            end
        end
    end

    function main(ptable::ISPRMPrecursorTable, mods_dict::Dict{String, Float32}, f_path::String)
        open(f_path) do f #"./data/NRF2_SIL.txt"

            header = readHeader(f)

            max_prot_id, max_pepGroup_id, max_pep_id, max_prec_id = UInt32(0), UInt32(0), UInt32(0), UInt32(0)
            pepGroup_id, pep_id, prot_id, prec_id = UInt32(1), UInt32(1), UInt32(1), UInt32(1)

            timetaken = @elapsed for row in eachline(f)

                protein_name, heavy_sequence, charge, isotope, transition_names = parseLine(row, header)
                #Add a check to make sure it is a heavy peptide to begin with. Warn and skip if not. 
                light_sequence = replace(heavy_sequence, r"\[H[a-z]{3}\]" =>"")
                unmodified_sequence = replace(heavy_sequence, r"\[[^\]]*\]"=>"")

                light_precursor = Precursor(light_sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict)
                heavy_precursor = Precursor(heavy_sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict)

                prot_id, max_prot_id = getProtID!(ptable, protein_name, max_prot_id)
                pepGroup_id, max_pepGroup_id = getPepGroupID!(ptable, protein_name, unmodified_sequence, max_pepGroup_id)
                pep_id, max_pep_id = getPepID!(ptable, light_sequence, heavy_sequence, max_pep_id, pepGroup_id)

                if !isinPrecursorSet(ptable, light_precursor) #If the precursor hasn't been encountered 
                    light_prec_id, heavy_prec_id, max_prec_id = addPrecursors!(ptable, light_precursor, heavy_precursor, max_prec_id, pep_id)
                    addTransitions!(ptable, transition_names, [light_prec_id, heavy_prec_id])
                end
            end
            println("Time to build precursor table ", timetaken);
        end
    end

    main(ptable, mods_dict, f_path)
end

mods_dict = Dict("Carb" => Float32(57.021464),
                 "Harg" => Float32(10.008269),
                 "Hlys" => Float32(8.014199))