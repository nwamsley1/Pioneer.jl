function parseIsotopeMods!(
    iso_mods_json::Dict{String, Any},
    iso_mod_to_mass_float::Dict{String, <:AbstractFloat})

    iso_mods = Vector{@NamedTuple{p::Regex, r::String}}()
    iso_mod_to_mass_dict = Dict{String, String}()
    n_iso_mods = length(iso_mods_json["pattern"])
    for n in range(1, n_iso_mods)
        pattern = Regex(iso_mods_json["pattern"][n])
        mass = Float32(iso_mods_json["mass"][n])
        name = iso_mods_json["mod_names"][n]
        push!(iso_mods,
        (
            p=pattern, 
            r=name
        ))
        iso_mod_to_mass_dict[name] = string(mass)
    end
    for (key, value) in pairs(iso_mod_to_mass_dict)
        iso_mod_to_mass_float[key] = parse(Float32, value)
    end
    return iso_mods, iso_mod_to_mass_float, iso_mods_json["label_name"]
end 

function addIsotopeMods!(
    isotope_mods_column::Vector{Union{Missing, String}},
    accession_number_column::AbstractVector{String},
    precursor_mz_column::Vector{Float32},
    precursor_sequence_column::Vector{String},
    precursor_charge_column::Vector{UInt8},
    mod_group_name::String,
    isotope_mods::Vector{@NamedTuple{p::Regex, r::String}},
    iso_mod_to_mass::Dict{String, Float32})

    for i in ProgressBar(range(1, length(isotope_mods_column)))
        precursor_sequence = precursor_sequence_column[i]
        accession_number_column[i] *= "|"*mod_group_name
        iso_mods_string = ""
        #difference in precursors mass between unmodified and modified 
        mass_offset = 0.0f0
        for isotope_mod in isotope_mods
            #For each amino acid matching the isotope_mod 
            for match in eachmatch(isotope_mod[:p], precursor_sequence)
                #seq_idx_to_iso_mod[match.offset] += isotope_mod[:mass]
                #For heavy arginine at position 7 could be "(7,R)"
                iso_mods_string *= "("*string(match.offset)*","*isotope_mod[:r]*")"
                mass_offset += iso_mod_to_mass[isotope_mod[:r]]
            end
        end
        #Convert mass offset to m/z offset
        mz_offset = mass_offset/precursor_charge_column[i]
        precursor_mz_column[i] += mz_offset
        #Update isotope mods column
        if length(iso_mods_string) > 0
            isotope_mods_column[i] = iso_mods_string
        else
            isotope_mods_column[i] = missing
        end
    end
end

function addIsotopeModifiedPrecursors!(
    precursors_df::DataFrame,
    iso_mod_to_mass::Dict{String, Float32},
    isotope_mods_groups::Vector{Any})
    #One data frame for each isotope mod group + unmodified 
    precursors_dfs = [precursors_df]
    #For each isotope mod group get the modified precursors dataframe 
    for iso_mod_set in isotope_mods_groups
        #Parse the isotope mods for the isotope mod group 
        isotope_mods, iso_mod_to_mass, mod_group_name = parseIsotopeMods!(
            iso_mod_set,
            iso_mod_to_mass
            )
        #New dataframe for modified precursors 
        isotope_modified_precursors_df = copy(precursors_df)
        addIsotopeMods!(
            isotope_modified_precursors_df[!,:isotope_mods],
            isotope_modified_precursors_df[!,:accession_number],
            isotope_modified_precursors_df[!,:mz],
            isotope_modified_precursors_df[!,:sequence],
            isotope_modified_precursors_df[!,:precursor_charge],
            mod_group_name,
            isotope_mods,
            iso_mod_to_mass,
        )
        #Remove rows for precursors that did not have an isotope mod
        filter!(x->!ismissing(x.isotope_mods), isotope_modified_precursors_df);
        #Add new dataframe 
        push!(precursors_dfs,
        isotope_modified_precursors_df);
    end
    #Concatentae the precursor dataframes 
    return vcat(precursors_dfs...);
end