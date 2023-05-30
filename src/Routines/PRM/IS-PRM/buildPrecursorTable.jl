"""
    SimplePrecursor

Minimal representation of a peptide with a specific mass, charge state, and isotopic state. 

### Fields

- sequence::String -- Single letter amino-acid representation of the precursor.
- charge::UInt8 -- Charge state
- isotope::UInt8 -- Isotope state. 0 is the monoisotope. 
- pep_id::UInt32 -- Unique identifier for the `Peptide` to which this `SimplePrecursor` corresponds

### Examples

- `PeptideGroup() = PeptideGroup(Set{UInt32}(), "")` -- constructor for an placeholder 

### GetterMethods

- getSeq(p::PeptideGroup) = p.sequence
- getCharge(sp::SimplePrecursor) = sp.charge
- getIsotope(sp::SimplePrecursor) = sp.isotope
- getPepID(sp::SimplePrecursor) = sp.pep_id

### Methods

None 

### Notes

- Immutable so can be used as a dictionary key or member of a set. 

"""
struct SimplePrecursor
    sequence::String
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end

getSeq(sp::SimplePrecursor) = sp.sequence
getCharge(sp::SimplePrecursor) = sp.charge
getIsotope(sp::SimplePrecursor) = sp.isotope
getPepID(sp::SimplePrecursor) = sp.pep_id

"""
    ParModel

Represents results of a par estimation model. 

### Fields

- par_model_coef::Matrix{Float64} -- Coefficients of the fit model
- goodness_of_fit::Float64 -- A measure of goodness of fit. Could be deviance ratio, r^2, whatever you want. 

### Examples

- `PeptideGroup() = ParModel() = ParModel(zeros((2,2)),0.0)` -- constructor for an placeholder 

### GetterMethods

- getCoef(pm::ParModel) = pm.par_model_coef
- getGoodnessOfFit(pm::ParModel) = pm.goodness_of_fit

### Methods

None 

### Notes

None

"""
struct ParModel{T} <: AbstractFloat
    par_model_coef::T
    goodness_of_fit::T
end

getCoef(pm::ParModel) = pm.par_model_coef
getGoodnessOfFit(pm::ParModel) = pm.goodness_of_fit
ParModel() = ParModel(0.0,0.0)

"""
    LightHeavyPrecursorPair

Represents a pair of precursors, one with a heavy isotopic label (internal standard) and one without.  

### Fields

- pepGroup_id::UInt32 -- Unique identifier for the `Peptide` to which this precursor pair belongs. 
- light_sequence::String -- Sequence of the unlabeled, "light" peptide 
- heavy_sequence::String -- Sequence of the isotope labeled, "heavy" peptide 
- light_prec_id::UInt32 -- Unique identifier of the unlabeled, "light" peptide 
- heavy_prec_id::UInt32 -- Unique identifier of the isotope labeled, "heavy" peptide 
- par_model::Dictionary{UInt32, ParModel} -- A table of `ParModel`. The raw file id could stand for the key and the value
                                            could be the model specific to that raw file.

### Examples

None

### GetterMethods

- getPepID(p::LightHeavyPrecursorPair) = p.pep_id
- getLightSeq(p::LightHeavyPrecursorPair) = p.light_sequence
- getHeavySeq(p::LightHeavyPrecursorPair) = p.heavy_sequence
- getLightPrecID(p::LightHeavyPrecursorPair) = p.light_prec_id
- getHeavyPrecID(p::LightHeavyPrecursorPair) = p.heavy_prec_id
- getLightScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs
- getHeavyScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs

### Methods

- setParModel(p::LightHeavyPrecursorPair; coef::Matrix{Float64} = zeros((2, 2)), dev_ratio::Float64 = 0.0, ms_file_idx::UInt32 = UInt32(0))

### Notes

None

"""
mutable struct LightHeavyPrecursorPair
    light_pep_id::UInt32 #change to heavy pep_id
    heavy_pep_id::UInt32
    light_sequence::String
    heavy_sequence::String
    light_prec_id::UInt32
    heavy_prec_id::UInt32
    par_model::Dictionary{UInt32, ParModel}
end

getPepID(p::LightHeavyPrecursorPair) = p.light_pep_id
getLightSeq(p::LightHeavyPrecursorPair) = p.light_sequence
getHeavySeq(p::LightHeavyPrecursorPair) = p.heavy_sequence
getLightPrecID(p::LightHeavyPrecursorPair) = p.light_prec_id
getHeavyPrecID(p::LightHeavyPrecursorPair) = p.heavy_prec_id
getLightScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs
getHeavyScanIDs(p::LightHeavyPrecursorPair) = p.light_scan_idxs

function setParModel(p::LightHeavyPrecursorPair; coef::T = 0.0, goodness_of_fit::T = 0.0, ms_file_idx::UInt32 = UInt32(0)) where T <: AbstractFloat
    if !isassigned(p.par_model, ms_file_idx)
        insert!(p.par_model, ms_file_idx, ParModel(coef,  goodness_of_fit))
    else
        p.par_model[ms_file_idx] = ParModel(coef,  goodness_of_fit)
    end
end

#getParModel(p::LightHeavyPrecursorPair) = 
"""
    PrecursorTable <: PrecursorDatabase

Inherits from `PrecursorDatabase`. Additional methods and fields for dealing with light/heavy precursor pairs. 

### Fields

All fields mandated for `PrecursorDatabase`. In addition the following: 

- lh_pair_id_to_light_heavy_pair::UnorderedDictionary{UInt32, LightHeavyPrecursorPair}
- pep_sequence_to_pep_id::UnorderedDictionary{String, UInt32}
- simple_precursor_set::Set{SimplePrecursor}
- prec_id_to_lh_pair_id::UnorderedDictionary{UInt32, UInt32}

### GetterMethods
In addition to mandatory field for `PrecursorDatabase`

- getIDToLightHeavyPair(p::ISPRMPrecursorTable) = p.lh_pair_id_to_light_heavy_pair
- getPepSeqToPepID(p::ISPRMPrecursorTable) = p.pep_sequence_to_pep_id
- getLightHeavyPairFromPrecID(p::ISPRMPrecursorTable, prec_id::UInt32) = getIDToLightHeavyPair(p)[prec_id]
- getSimplePrecursors(p::ISPRMPrecursorTable) = p.simple_precursor_set

### Methods

- addPrecIDToLHPairID(p::ISPRMPrecursorTable, prec_id::UInt32, lh_pair_id::UInt32)
- isinPrecursorSet(p::ISPRMPrecursorTable, mz::Float32, charge::UInt8, isotope::UInt8, pep_id::UInt32)
- isinPrecursorSet(p::ISPRMPrecursorTable, prec::Precursor)

### Examples

- `ISPRMPrecursorTable() = ISPRMPrecursorTable(UnorderedDictionary{UInt32, Protein}(),
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
                                            UnorderedDictionary{UInt32, UInt32}())` -- constructor for a placeholder 

"""
mutable struct ISPRMPrecursorTable{T<:AbstractFloat} <: PrecursorDatabase
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    prot_to_id::UnorderedDictionary{String, UInt32}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    pepGroup_to_id::UnorderedDictionary{String, UInt32}
    id_to_pep::UnorderedDictionary{UInt32, Peptide}
    id_to_prec::Dictionary{UInt32, Precursor{T}}
    sorted_prec_ids::Vector{UInt32}
    prec_id_to_transitions::Dictionary{UInt32, Vector{Transition}}
    #Specific to ISPRMPrecursorTable
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
                                            Dictionary{UInt32, Precursor{Float64}}(),
                                            Vector{UInt32}(),
                                            Dictionary{UInt32, Vector{Transition}}(), 
                                            UnorderedDictionary{UInt32, LightHeavyPrecursorPair}(), 
                                            UnorderedDictionary{String, UInt32}(), 
                                            Set{SimplePrecursor}(),
                                            UnorderedDictionary{UInt32, UInt32}())

getIDToLightHeavyPair(p::ISPRMPrecursorTable) = p.lh_pair_id_to_light_heavy_pair
getPepSeqToPepID(p::ISPRMPrecursorTable) = p.pep_sequence_to_pep_id
getPrecIDToLHPairID(p::ISPRMPrecursorTable) = p.prec_id_to_lh_pair_id
getLightHeavyPairFromPrecID(p::ISPRMPrecursorTable, prec_id::UInt32) = getIDToLightHeavyPair(p)[prec_id]
getSimplePrecursors(p::ISPRMPrecursorTable) = p.simple_precursor_set

function addPrecIDToLHPairID(p::ISPRMPrecursorTable, prec_id::UInt32, lh_pair_id::UInt32)
    if !isassigned(getPrecIDToLHPairID(p), prec_id)
        insert!(getPrecIDToLHPairID(p), prec_id, lh_pair_id)
    end
end

isinPrecursorSet(p::ISPRMPrecursorTable, mz::T, charge::UInt8, isotope::UInt8, pep_id::UInt32) where {T<:AbstractFloat} = SimplePrecursor(mz, charge, isotope, pep_id) ∈ getSimplePrecursorSet(p)
isinPrecursorSet(p::ISPRMPrecursorTable, prec::Precursor{T}) where {T<:AbstractFloat}= isinPrecursorSet(p, getMZ(prec), getCharge(prec), getIsotope(prec), getPepID(prec)) ∈ getSimplePrecursorSet(p)

"""
    buildPrecursorTable!(ptable::ISPRMPrecursorTable, mods_dict::Dict{String, Float32}, f_path::String)

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

"""
function buildPrecursorTable!(ptable::ISPRMPrecursorTable, mods_dict::Dict{String, T}, f_path::String) where {T<:AbstractFloat}

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

    function getProtID!(ptable::PrecursorDatabase, protein_name::String, max_prot_id::UInt32)
        if !containsProt(ptable, protein_name) #If the protien hasn't been encountered,
            #then add it to the hash table
            max_prot_id += UInt32(1); #Increase the counter
            prot_id = max_prot_id #Set the current protein group ID
            addProtein!(ptable, protein_name, prot_id);
        else #If the protein has been encountered, get its protien ID
            prot_id = getProtIDFromName(ptable, protein_name) #Set the current protein group ID
        end
        return prot_id, max_prot_id
    end

    function getPepGroupID!(ptable::PrecursorDatabase, protein_name::String, unmodified_sequence::String, max_pepGroup_id::UInt32)
        if !containsPepGroup(ptable, unmodified_sequence) #If the peptide group hasn't been encountered
            max_pepGroup_id += UInt32(1) #Increase the counter
            pepGroup_id = max_pepGroup_id #Set the current pepGroup_id
            addPepGroup!(ptable, unmodified_sequence, pepGroup_id, Set(UInt32[]), protein_name); #add it to the hash table,
            addPepGroupToProtein!(ptable, protein_name, unmodified_sequence); #add a new peptide group to the protein,
        else #If this peptide has been encountered before, just update the protein and peptide groups
            addProteinToPepGroup!(ptable, protein_name, unmodified_sequence); #Add the current protein to this peptide group
            addPepGroupToProtein!(ptable, protein_name, unmodified_sequence); #Add the current peptide group to this protein
            pepGroup_id = getPepGroupIDFromSequence(ptable, unmodified_sequence) #Set the current pepGroup_id
        end
        return pepGroup_id, max_pepGroup_id
    end 

    function getPepID!(ptable::PrecursorDatabase, light_sequence::String, heavy_sequence::String, max_pep_id::UInt32, pepGroup_id::UInt32)
        if !containsPep(ptable, light_sequence) #If the peptide hasn't been encountered
            light_pep_id = max_pep_id + UInt32(1)
            heavy_pep_id = max_pep_id + UInt32(2)
            insertPep!(ptable, light_sequence, light_pep_id, pepGroup_id)
            insertPep!(ptable, heavy_sequence, heavy_pep_id, pepGroup_id)
            addPepToPepGroup!(ptable, pepGroup_id, light_pep_id)
            addPepToPepGroup!(ptable, pepGroup_id, heavy_pep_id)
            max_pep_id += UInt32(2);  
            return light_pep_id, heavy_pep_id, max_pep_id
        else
            light_pep_id = getPepSeqToPepID(ptable)[light_sequence]
            heavy_pep_id = getPepSeqToPepID(ptable)[heavy_sequence]
            return light_pep_id, heavy_pep_id, max_pep_id
        end
    end
    
    function addPrecursors!(ptable::ISPRMPrecursorTable, light_precursor::Precursor, light_sequence::String, heavy_precursor::Precursor{T}, heavy_sequence::String, light_pep_id::UInt32, heavy_pep_id::UInt32, max_lh_pair_id::UInt32, max_prec_id::UInt32) where {T<:AbstractFloat}
        
        max_lh_pair_id += UInt32(1)
        light_prec_id = getPrecID(light_precursor)
        heavy_prec_id = getPrecID(heavy_precursor)
        if !containsPrecID(ptable, light_prec_id)
            insertPrecursor!(ptable, light_precursor)
        end
        if !containsPrecID(ptable, heavy_prec_id)
            insertPrecursor!(ptable, heavy_precursor)
        end
        insert!(getIDToLightHeavyPair(ptable),
                max_lh_pair_id,
                LightHeavyPrecursorPair(light_pep_id,
                                        heavy_pep_id,
                                        light_sequence, 
                                        heavy_sequence,
                                        light_prec_id, 
                                        heavy_prec_id, 
                                        #Dictionary{Int64, NamedTuple{(:lower_bound, :upper_bound), Tuple{Int64, Int64}}}(),
                                        Dictionary{UInt32, ParModel}())
                )
        #push!(getPrecIDs(getIDToPep(ptable)[getPepID(light_precursor)]), getPrecID(light_precursor))
        #push!(getPrecIDs(getIDToPep(ptable)[getPepID(light_precursor)]), max_lh_pair_id)
        addPrecIDToLHPairID(ptable, light_prec_id, max_lh_pair_id)
        addPrecIDToLHPairID(ptable, heavy_prec_id, max_lh_pair_id)
        push!(getPrecIDs(
                        getPep(ptable, getPepIDFromPrecID(ptable, light_prec_id))), 
                        light_prec_id)
        push!(getPrecIDs(
                            getPep(ptable, getPepIDFromPrecID(ptable, heavy_prec_id))), 
                            heavy_prec_id)
        return max_lh_pair_id
    end

    function addTransitions!(ptable::ISPRMPrecursorTable, transition_names::Vector{String}, prec_ids::Vector{UInt32})
        #for transition_name in line[header["transition_names"]:end]
        for prec_id in prec_ids 
            if hasTransitions(ptable, prec_id)
                continue
            end
            insert!(getPrecIDToTransitions(ptable),prec_id,Vector{Transition}())
            for transition_name in transition_names
                transition = parseTransitionName(transition_name)
                addTransition!(ptable, 
                                prec_id, 
                                Transition(getResidues(getPrecursor(ptable, prec_id)),
                                        ion_type = transition[:ion_type],
                                        ind = transition[:ind],
                                        charge = transition[:charge],
                                        isotope = UInt8(0),
                                        prec_id = prec_id)
                                )
            end
        end
    end

    function main(ptable::ISPRMPrecursorTable, mods_dict::Dict{String, T}, f_path::String) where {T<:AbstractFloat}
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
                light_pep_id, heavy_pep_id, max_pep_id = getPepID!(ptable, light_sequence, heavy_sequence, max_pep_id, pepGroup_id)
                precursor = SimplePrecursor(light_sequence, charge, isotope, light_pep_id)
                if precursor ∉ ptable.simple_precursor_set
                    push!(ptable.simple_precursor_set, precursor)

                    max_prec_id += UInt32(2)
                    light_prec_id, heavy_prec_id = max_prec_id - UInt32(1), max_prec_id
                    light_precursor = Precursor(light_sequence, charge = charge, isotope = isotope, pep_id = light_pep_id, prec_id = light_prec_id, mods_dict = mods_dict)
                    heavy_precursor = Precursor(heavy_sequence, charge = charge, isotope = isotope, pep_id = heavy_pep_id, prec_id = heavy_prec_id, mods_dict = mods_dict)

                    max_lh_pair_id = addPrecursors!(ptable, light_precursor, light_sequence, heavy_precursor, heavy_sequence, light_pep_id, heavy_pep_id, max_lh_pair_id, max_prec_id)
                    addTransitions!(ptable, transition_names, [light_prec_id, heavy_prec_id])
                end
            end

            setSortedPrecursorKeys!(ptable)
            println("Time to build precursor table ", timetaken);
        end
    end

    main(ptable, mods_dict, f_path);
end
#=struct SimplePrecursor
    sequence::String
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end=#
const mods_dict = Dict("Carb" => Float64(57.021464),
                 "Harg" => Float64(10.008269),
                 "Hlys" => Float64(8.014199),
                 "Hglu" => Float64(6))
                 