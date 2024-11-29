function countSulfurs( plain_sequence::AbstractString, 
                        mods::String,
                        mods_to_sulfur_diff::Dict{String, Int8})::Int8
    sulfur_count = zero(UInt8)
    for aa in plain_sequence
        sulfur_count += (aa=='C')|(aa=='M')
    end

    for mod in parseMods(mods)
        if haskey(mods_to_sulfur_diff, getModName(mod.match))
            n_sulfur = mods_to_sulfur_diff[mod_string]
            seq_idx_to_sulfur[getModIndex(mod.match)] += n_sulfur
            sulfur_count += n_sulfur
        end
    end
    return sulfur_count
end

function countSulfurs( plain_sequence::AbstractString, 
                        mods::Missing,
                        mods_to_sulfur_diff::Dict{String, Int8})::Int8
    sulfur_count = zero(UInt8)
    for aa in plain_sequence
        sulfur_count += (aa=='C')|(aa=='M')
    end
    return sulfur_count
end

function parseChronologerOutput(
    path_to_precursors::String,
    pion_lib_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    isotope_mods_groups::Vector{Any},
    rt_bin_tol::AbstractFloat)
    #Read chronologer output
    precursors_df = CSV.read(path_to_precursors,
                                DataFrame, 
                                types=Dict(
                                :upid => String,
                                :accession_number => String,
                                :sequence => String,
                                :mods => String,
                                :precursor_charge => UInt8,
                                :collision_energy => Float32,
                                :decoy => Bool,
                                :mz => Float32,
                                :iRT => Float32,
                                :entrapment_group_id => UInt8,
                                :base_pep_id => UInt32
                                ))
    #add/rename  columns
    #precursors_df[!,:base_pep_id] = UInt32.(collect(range(1, size(precursors_df, 1))))
    rename!(precursors_df, :iRT => :irt)
    rename!(precursors_df, :upid => :proteome_identifiers)
    precursors_df[!,:length] = zeros(UInt8, size(precursors_df, 1))
    for i in range(1, size(precursors_df,1))
        precursors_df[i,:length] = UInt8(length(precursors_df[i,:sequence]))
    end
    precursors_df[!,:sulfur_count] = zeros(UInt8, size(precursors_df, 1))
    for i in range(1, size(precursors_df,1))
        precursors_df[i,:sulfur_count] = countSulfurs(precursors_df[i,:sequence],
                                                    precursors_df[i,:mods],
                                                    mods_to_sulfur_diff)
    end
    precursors_df[!,:isotope_mods] = Vector{Union{Missing, String}}(undef, size(precursors_df, 1));
    for i in range(1, size(precursors_df, 1))
        precursors_df[i,:isotope_mods] = missing
    end
    #Duplicates and appends the precursors_df for each isotope_mods group. 
    #At this time should only be used for SILAC. 
    println("pre isotope: ", size(precursors_df))
    precursors_df = addIsotopeModifiedPrecursors!(
                        precursors_df,
                        iso_mod_to_mass,
                        isotope_mods_groups
                        )
    println("post isotope: ", size(precursors_df))
    ########
    #Sort 
    #Need same sort order as fragment index.
    #1) Sort entire data frame in ascending order of irt
    #2) Within irt bins of identical width, sort by mz
    sort!(precursors_df, :irt)
    start_idx, stop_idx = 1, 1
    start_irt, stop_irt = first(precursors_df[!,:irt]), first(precursors_df[!,:irt])

    #Get irt bins and sort by mz within these 
    for pid in range(1, size(precursors_df, 1))
        stop_idx = pid
        stop_irt = precursors_df[stop_idx,:irt]
        if ((stop_irt - start_irt) > rt_bin_tol) & (stop_idx > start_idx)
            stop_idx -= 1
            sort!(@view(precursors_df[start_idx:stop_idx,:]), :mz)
            start_idx = pid
            start_irt = precursors_df[pid,:irt]
        end
    end
    sort!(precursors_df[start_idx:end,:],:mz)
    dir, filename = splitdir(path_to_precursors)
    name, _ = splitext(filename)
    precursors_arrow_path = joinpath(pion_lib_dir, name*".arrow")
    Arrow.write(
        precursors_arrow_path,
        precursors_df
    )   
    return precursors_arrow_path
end

