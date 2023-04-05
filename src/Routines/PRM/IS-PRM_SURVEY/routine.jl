using JSON

##########
#Parse arguments
##########
    # Parse the first argument as an integer
    params = JSON.parse(read(ARGS[1], String))
    MS_DATA_DIR = ARGS[2]
    PRECURSOR_LIST_PATH = ARGS[3]

    MS_TABLE_PATHS = filter(file -> isfile(joinpath(dir_path, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))


    # Print the argument
    println("Parameters: $params")

    function parse_mods(fixed_mods)
    fixed_mods_parsed = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
    for mod in fixed_mods
        push!(fixed_mods_parsed, (p=Regex(mod[1]), r = mod[2]))
    end
    return fixed_mods_parsed
    end
    #Parse argments
    params = (
    right_precursor_tolerance = Float32(params["right_precursor_tolerance"]),
    left_precursor_tolerance = Float32(params["left_precursor_tolerance"]),
    precursor_rt_tolerance = Float32(params["precursor_rt_tolerance"]),
    b_ladder_start = UInt8(params["b_ladder_start"]),
    y_ladder_start = UInt8(params["y_ladder_start"]),
    precursor_charges = [UInt8(charge) for charge in params["precursor_charges"]],
    precursor_isotopes = [UInt8(isotope) for isotope in params["precursor_isotopes"]],
    transition_charges = [UInt8(isotope) for isotope in params["transition_charges"]],
    fragment_match_ppm = Float32(params["fragment_match_ppm"]),
    minimum_fragment_count = UInt8(params["minimum_fragment_count"]),
    precursort_rt_window = Float32(params["precursor_rt_window"]),
    max_variable_mods = Int(params["max_variable_mods"]),
    fixed_mods = parse_mods(params["fixed_mods"]),
    modification_masses = Dict{String, Float32}(k => Float32(v) for (k, v) in params["modification_masses"])
    )

    #Dict{String, Float32}(k => Float32(v) for (k, v) in params["modification_masses"])
    println("params again $params")
    println(params[:b_ladder_start])

##########
#Read Precursor Table
##########
    ptable = PrecursorTable()
    buildPrecursorTable!(ptable, 
                        params[:fixed_mods], 
                        params[:var_mods], 
                        params[:max_variable_mods], 
                        PRECURSOR_LIST_PATH)
    addPrecursors!(
                        ptable, 
                        params[:precursor_charges], 
                        params[:precursor_isotopes], 
                        params[:modification_masses]
                        )
##########
#Search Survey Runs
##########
MS_TABLES = Arrow.Table[]
combined_scored_psms = makePSMsDict(FastXTandem())
combined_fragment_matches = Vector{Vector{FragmentMatch}}()
    for (mx_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)

        push!(MS_TABLES, Arrow.Table(MS_TABLE_PATH))

        scored_psms, fragment_matches = SearchRAW(
                                                MS_TABLES[end], 
                                                getPrecursors(ptable), 
                                                selectTransitionsPRM, 
                                                params[:right_precursor_tolerance],
                                                params[:left_precursor_tolerance],
                                                params[:transition_charges],
                                                params[:transition_isotopes],
                                                params[:b_ladder_start],
                                                params[:y_ladder_start],
                                                params[:fragment_match_ppm],
                                                ms_file_idx
                                                )
        for key in keys(combined_scored_psms)
            append!(combined_scored_psms[key], scored_psms[key])
        end
        push!(combined_fragment_matches, fragment_matches)
    end

##########
#Get Best PSMs for Each Peptide
##########
    best_psms = FilterPSMs(combined_scored_psms, ptable, MS_TABLE, params[:minimum_fragment_count])
##########
#Get MS1 Peak Heights
##########
    #Peak heights are zero to begin with
    ms1_peak_heights = UnorderedDictionary(best_psms[!,:precursor_idx], zeros(Float32, len(best_psms[!,:precursor_idx])))
    for (ms_file_idx, MS_TABLE) in enumerate(MS_TABLES)
        getMS1PeakHeights!(
                            MS_TABLE[:retentionTime], 
                            MS_TABLE[:masses], 
                            MS_TABLE[:intensities], 
                            MS_TABLE[:msOrder], 
                            MS1_PEAK_HEIGHTS, 
                            best_psms[!,:retentionTime], 
                            best_psms[!,:precursor_idx], 
                            best_psms[!,:ms_file_idx],
                            getSimplePrecursors(ptable), 
                            Float32(0.25), 
                            params[:right_precursor_tolerance], 
                            params[:left_precursor_tolerance],
                            ms_file_idx)
    end

    #Add MS1 Heights to the best_psms DataFrame 
    transform!(best_psms, AsTable(:) => ByRow(psm -> ms1_peak_heights[psm[:precursor_idx]]) => :ms1_peak_height)
##########
#Get Chromatograms for the Best Precursors
##########
    precursor_chromatograms = initPrecursorChromatograms(best_psms) |> (best_psms -> fillPrecursorChromatograms(best_psms, 
                                                                                                                fragment_matches, 
                                                                                                                MS_TABLE, 
                                                                                                                params[:precursor_rt_tolerance]
                                                                                                                )
                                                                    )  
    #Names and charges for the "n" most intense fragment ions for each precursor
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestTransitions(getBestPSM(precursor_chromatograms[psm[:precursor_idx]]),
                                                                        params[:maximum_fragment_count])) => :best_transitions)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:precursor_idx]])[:name][psm[:best_transitions]]) => :transition_names)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)

    writeTransitionList(best_psms, joinpath(MS_DATA_DIR, "transition_list.csv"))
    writeIAPIMethod(best_psms, joinpath(MS_DATA_DIR, "iapi_method.csv"))
##########
#Get Chromatograms for the Best Precursors
##########

