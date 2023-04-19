struct SimplePrecursor
    sequence::String
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end
mutable struct teststruct
    a::Tuple{Int64, Int64}
end
getSeq(sp::SimplePrecursor) = sp.sequence
getCharge(sp::SimplePrecursor) = sp.charge
getIsotope(sp::SimplePrecursor) = sp.isotope
getPepID(sp::SimplePrecursor) = sp.pep_id

struct ParModel
    par_model_coef::Matrix{Float64}
    dev_ratio::Float64
end
ParModel() = ParModel(zeros((2,2)),0.0)
mutable struct LightHeavyPrecursorPair
    pepGroup_id::UInt32
    light_sequence::String
    heavy_sequence::String
    light_prec_id::UInt32
    heavy_prec_id::UInt32
    #integration_bounds::Dictionary{UInt32, NamedTuple{(:lower_bound, :upper_bound), Tuple{Int64, Int64}}}
    par_model::Dictionary{UInt32, ParModel}
end

function setParModel(p::LightHeavyPrecursorPair; coef::Matrix{Float64} = zeros((2, 2)), dev_ratio::Float64 = 0.0, ms_file_idx::UInt32 = UInt32(0))
    if !isassigned(p.par_model, ms_file_idx)
        insert!(p.par_model, ms_file_idx, ParModel(coef, dev_ratio))
    else
        p.par_model[ms_file_idx] = ParModel(coef, dev_ratio)
    end
end

getPepGroupID(p::LightHeavyPrecursorPair) = p.pepGroup_id
getLightSeq(p::LightHeavyPrecursorPair) = p.light_sequence
getHeavySeq(p::LightHeavyPrecursorPair) = p.heavy_sequence
getLightPrecID(p::LightHeavyPrecursorPair) = p.light_prec_id
getHeavyPrecID(p::LightHeavyPrecursorPair) = p.heavy_prec_id
getLightScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs
getHeavyScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs
getIntegrationBounds(p::LightHeavyPrecursorPair) = p.integration_bounds

function setIntegrationBounds!(p::LightHeavyPrecursorPair, run_idx::Int64, bounds::NamedTuple{(:lower_bound, :upper_bound), Tuple{Int64, Int64}})
    if !isassigned(getIntegrationBounds(p), run_idx)
        insert!(getIntegrationBounds(p), run_idx, bounds)
    else
        getIntegrationBounds(p)[run_idx] = bounds
    end
end

inBounds(p::LightHeavyPrecursorPair, run_idx, i::Int64) = (i>=getIntegrationBounds(p)[run_idx][:lower_bound]) & (i<=getIntegrationBounds(p)[run_idx][:upper_bound]) 

#addProtID(p::LightHeavyPair, prot_id::UInt32) = push!(getProtIDs(p), prot_id)

mutable struct ISPRMPrecursorTable <: PrecursorDatabase
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide} #Map peptide IDs to peptide group
    prec_id_to_precursor::Dictionary{UInt32, Precursor} #Needs to be sortable by precursor mass, therfore, not an UnorderedDictioanry. 
    sorted_precursor_keys::Vector{UInt32}
    prec_id_to_transitions::UnorderedDictionary{UInt32, Vector{Transition}}
    lh_pair_id_to_light_heavy_pair::UnorderedDictionary{UInt32, LightHeavyPrecursorPair}
    pep_sequence_to_pep_id::UnorderedDictionary{String, UInt32}
    simple_precursor_set::Set{SimplePrecursor}
    prec_id_to_lh_pair_id::UnorderedDictionary{UInt32, UInt32}
end

ISPRMPrecursorTable() = ISPRMPrecursorTable(UnorderedDictionary{UInt32, Protein}(),
                                            UnorderedDictionary{String, UInt32}(),
                                            UnorderedDictionary{UInt32, PeptideGroup}(),
                                            UnorderedDictionary{String, UInt32}(),
                                            UnorderedDictionary{UInt32, Peptide}(),
                                            Dictionary{UInt32, Precursor}(),
                                            Vector{UInt32}(),
                                            UnorderedDictionary{UInt32, Vector{Transition}}(), 
                                            UnorderedDictionary{UInt32, LightHeavyPrecursorPair}(), 
                                            UnorderedDictionary{String, UInt32}(), 
                                            Set{SimplePrecursor}(),
                                            UnorderedDictionary{UInt32, UInt32}())

getIDToLightHeavyPair(p::ISPRMPrecursorTable) = p.lh_pair_id_to_light_heavy_pair
getPepSequenceToPepID(p::ISPRMPrecursorTable) = p.pep_sequence_to_pep_id
getLightHeavyPairFromPrecID(p::ISPRMPrecursorTable, prec_id::UInt32) = getPrecIDToLightHeavyPair(p)[prec_id]
getSimplePrecursorSet(p::ISPRMPrecursorTable) = p.simple_precursor_set

function addPrecIDToLHPairID(p::ISPRMPrecursorTable, prec_id::UInt32, lh_pair_id::UInt32)
    if !isassigned(p.prec_id_to_lh_pair_id, prec_id)
        insert!(p.prec_id_to_lh_pair_id, prec_id, lh_pair_id)
    end
end

isinPrecursorSet(p::ISPRMPrecursorTable, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32) = SimplePrecursor(mz, charge, isotope, pep_id) ∈ getSimplePrecursorSet(p)
isinPrecursorSet(p::ISPRMPrecursorTable, prec::Precursor) = isinPrecursorSet(p, getMZ(prec), getCharge(prec), getIsotope(prec), getPepID(prec))
#addToPrecursorSet!(p::ISPRMPrecursorTable, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32) = push!(p.simple_precursor_set, SimplePrecursor(mz, charge, isotope, pep_id))
#addToPrecursorSet!(p::ISPRMPrecursorTable, prec::Precursor) = push!(p.simple_precursor_set, SimplePrecursor(getMZ(prec), getCharge(prec), getIsotope(prec), getPepID(prec)))
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
        row = map(string, split(row, ","))
        return row[header["protein_name"]], row[header["sequence"]], parse(UInt8, row[header["precursor_charge"]]), parse(UInt8, row[header["precursor_isotope"]]), map(string, split(row[header["transition_names"]],";"))
    end 

    function parseTransitionName(transition_name::String)
        return (ion_type = transition_name[1], ind = parse(UInt8, split(transition_name,"+")[1][2:end]), charge = parse(UInt8, transition_name[end]))
    end

    function getProtID!(ptable::ISPRMPrecursorTable, protein_name::String, max_prot_id::UInt32)
        if !containsProt(ptable, protein_name) #If the protien hasn't been encountered,
            #then add it to the hash table\
            max_prot_id += UInt32(1);  
            prot_id = max_prot_id
            addNewProtein!(protein_name, prot_id, ptable);
        else
            prot_id = getProtToID(ptable)[protein_name]
        end
        return prot_id, max_prot_id
    end

    function getPepGroupID!(ptable::ISPRMPrecursorTable, protein_name::String, unmodified_sequence::String, max_pepGroup_id::UInt32)
        if !containsPepGroup(ptable, unmodified_sequence) #If the peptide hasn't been encountered
            max_pepGroup_id += UInt32(1)
            pepGroup_id = max_pepGroup_id
            addNewPeptideGroup!(unmodified_sequence, pepGroup_id, protein_name, ptable); #add it to the hash table,
            addPepGroupToProtein!(ptable, protein_name, unmodified_sequence); #add a new peptide group to the protein,
        else #If this peptide has been encountered before, just update the protein and peptide groups
            addProteinToPepGroup!(ptable, protein_name, unmodified_sequence); #Add the current protein to this peptide group
            addPepGroupToProtein!(ptable, protein_name, unmodified_sequence); #Add the current peptide group to this protein
            pepGroup_id = getPepGroupToID(ptable)[unmodified_sequence]
        end
        return pepGroup_id, max_pepGroup_id
    end 

    function getPepID!(ptable::ISPRMPrecursorTable, light_sequence::String, max_pep_id::UInt32, pepGroup_id::UInt32)
        if !isassigned(getPepSequenceToPepID(ptable), light_sequence) #If the peptide hasn't been encountered
            max_pep_id += UInt32(1);  
            pep_id = max_pep_id 
            insert!(getPepSequenceToPepID(ptable), light_sequence, pep_id)
            insert!(getIDToPep(ptable), 
                    pep_id,
                    Peptide(light_sequence, pepGroup_id, Set(UInt32[])))
            push!(getIDToPepGroup(ptable)[pepGroup_id].pep_ids, pep_id)
        else
            pep_id = getPepSequenceToPepID(ptable)[light_sequence]
        end
        pep_id, max_pep_id
    end
    
    function addPrecursors!(ptable::ISPRMPrecursorTable, light_precursor::Precursor, heavy_precursor::Precursor, pep_id::UInt32, max_lh_pair_id::UInt32, max_prec_id::UInt32)
        max_prec_id += UInt32(2)
        max_lh_pair_id += UInt32(1)
        light_prec_id, heavy_prec_id = max_prec_id - UInt32(1), max_prec_id

        #addToPrecursorSet!(ptable, light_precursor)

        insert!(getPrecursors(ptable), light_prec_id, light_precursor)
        insert!(getPrecursors(ptable), heavy_prec_id, heavy_precursor)
        insert!(getIDToLightHeavyPair(ptable),
                max_lh_pair_id,
                LightHeavyPrecursorPair(pep_id, 
                                        getSeq(getPep(ptable, getPepID(light_precursor))), 
                                        getSeq(getPep(ptable, getPepID(heavy_precursor))),
                                        light_prec_id, 
                                        heavy_prec_id, 
                                        #Dictionary{Int64, NamedTuple{(:lower_bound, :upper_bound), Tuple{Int64, Int64}}}(),
                                        Dictionary{UInt32, ParModel}())
                )
        push!(getPrecIDs(getIDToPep(ptable)[getPepID(light_precursor)]), max_lh_pair_id)
        addPrecIDToLHPairID(ptable, light_prec_id, max_lh_pair_id)
        addPrecIDToLHPairID(ptable, heavy_prec_id, max_lh_pair_id)
        return light_prec_id, heavy_prec_id, max_prec_id, max_lh_pair_id
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
                            ion_type = transition[:ion_type],
                            ind = transition[:ind],
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

            max_prot_id, max_pepGroup_id, max_pep_id, max_lh_pair_id, max_prec_id = UInt32(0), UInt32(0), UInt32(0), UInt32(0), UInt32(0)
            pepGroup_id, pep_id, prot_id, prec_id = UInt32(1), UInt32(1), UInt32(1), UInt32(1)

            timetaken = @elapsed for row in eachline(f)

                protein_name, heavy_sequence, charge, isotope, transition_names = parseLine(row, header)
                #Add a check to make sure it is a heavy peptide to begin with. Warn and skip if not. 
                light_sequence = replace(heavy_sequence, r"\[H[a-z]{3}\]" =>"")
                unmodified_sequence = replace(heavy_sequence, r"\[[^\]]*\]"=>"")

                prot_id, max_prot_id = getProtID!(ptable, protein_name, max_prot_id)
                pepGroup_id, max_pepGroup_id = getPepGroupID!(ptable, protein_name, unmodified_sequence, max_pepGroup_id)
                pep_id, max_pep_id = getPepID!(ptable, light_sequence, max_pep_id, pepGroup_id)
                precursor = SimplePrecursor(light_sequence, charge, isotope, pep_id)
                if precursor ∉ ptable.simple_precursor_set
                    push!(ptable.simple_precursor_set, precursor)
                    light_precursor = Precursor(light_sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict)
                    heavy_precursor = Precursor(heavy_sequence, charge = charge, isotope = isotope, pep_id = pep_id, prec_id = prec_id, mods_dict = mods_dict)
                    light_prec_id, heavy_prec_id, max_prec_id, max_lh_pair_id = addPrecursors!(ptable, light_precursor, heavy_precursor, pep_id, max_lh_pair_id, max_prec_id)
                    addTransitions!(ptable, transition_names, [light_prec_id, heavy_prec_id])
                end
            end
            println("Time to build precursor table ", timetaken);
        end
    end

    main(ptable, mods_dict, f_path);
    makeSortedPrecursorKeys!(ptable);
end
#=struct SimplePrecursor
    sequence::String
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end=#
mods_dict = Dict("Carb" => Float32(57.021464),
                 "Harg" => Float32(10.008269),
                 "Hlys" => Float32(8.014199))