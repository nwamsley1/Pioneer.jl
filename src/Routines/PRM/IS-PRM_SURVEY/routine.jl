#julia ./src/Routines/PRM/IS-PRM_SURVEY/routine.jl ./data/test.json ./data/parquet/ ./data/NRF2_SIL.txt
#julia ./src/Routines/PRM/IS-PRM_SURVEY/routine.jl ./data/test.json /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/survey_runs /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/KINOME_MAR23.txt
#julia ./src/Routines/PRM/IS-PRM_SURVEY/routine.jl ./data/test.json /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/survey_runs /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/KINOME_MAR23.txt
using JSON
using PrettyPrinting
using PDFmerger
##########
#Parse arguments
##########

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "params_json"
            help = "Path to a .json file with the parameters"
            required = true
        "data_dir"
            help = "Path to a folder with .arrow MS data tables"
            required = true
        "precursor_list"
            help = "Path to a tab delimited table of precursors"
            required = true
        "--make_plots", "-p"
            help = "Whether to make plots. Defaults to `true`"
            default = true
        "--print_params", "-s"
            help = "Whether to print the parameters from the json. Defaults to `false`"
            default = false 
    end

    return parse_args(s)
end

ARGS = parse_commandline()


params = JSON.parse(read(ARGS["params_json"], String))
MS_DATA_DIR = ARGS["data_dir"]
PRECURSOR_LIST_PATH = ARGS["precursor_list"]
#MS_DATA_DIR = "./data/parquet/"
#PRECURSOR_LIST_PATH = "./data/NRF2_SIL.txt"
#Get all files in the `MS_DATA_DIR` ending in ".arrow" and append their names to the `MS_DATA_DIR` path. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))]


# Print the argument
#println("Parameters: $params")
#println("MS_DATA_DIR: $MS_DATA_DIR")
println("Processing: "*string(length(MS_TABLE_PATHS))*" files")
#println("PRECURSOR_LIST_PATH: $PRECURSOR_LIST_PATH")

function parse_mods(fixed_mods)
    fixed_mods_parsed = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
    for mod in fixed_mods
        push!(fixed_mods_parsed, (p=Regex(mod[1]), r = mod[2]))
    end
    return fixed_mods_parsed
end

#Parse argments
#params = JSON.parse(read("./data/test.json", String))
params = (
right_precursor_tolerance = Float32(params["right_precursor_tolerance"]),
left_precursor_tolerance = Float32(params["left_precursor_tolerance"]),
precursor_rt_tolerance = Float64(params["precursor_rt_tolerance"]),
b_ladder_start = Int64(params["b_ladder_start"]),
y_ladder_start = Int64(params["y_ladder_start"]),
precursor_charges = [UInt8(charge) for charge in params["precursor_charges"]],
precursor_isotopes = [UInt8(isotope) for isotope in params["precursor_isotopes"]],
transition_charges = [UInt8(charge) for charge in params["transition_charges"]],
transition_isotopes = [UInt8(isotope) for isotope in params["transition_isotopes"]],
fragment_match_ppm = Float32(params["fragment_match_ppm"]),
minimum_fragment_count = UInt8(params["minimum_fragment_count"]),
fragments_to_select = UInt8(params["fragments_to_select"]),
precursort_rt_window = Float32(params["precursor_rt_window"]),
max_variable_mods = Int(params["max_variable_mods"]),
fixed_mods = parse_mods(params["fixed_mods"]),
variable_mods = parse_mods(params["variable_mods"]),
modification_masses = Dict{String, Float32}(k => Float32(v) for (k, v) in params["modification_masses"]),
ms_file_conditions = params["ms_file_conditions"]
)

#Dict{String, Float32}(k => Float32(v) for (k, v) in params["modification_masses"])
if ARGS["print_params"]
    for (k,v) in zip(keys(params), params)
        pprint(k)
        print(" => ")
        pprintln(v)
    end
end
##########
#Load Dependencies 
##########
using Arrow, DataFrames, Tables
using Plots
include("../../../precursor.jl")
include("../../../binaryRangeQuery.jl")
include("../../../matchpeaks.jl")
include("../../../getPrecursors.jl")
include("../../../PSM_TYPES/PSM.jl")
include("../../../PSM_TYPES/FastXTandem.jl")
#include("../../../searchSpectra.jl")
include("../../../Routines/PRM/getBestPSMs.jl")
include("../../../Routines/PRM/precursorChromatogram.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/plotPRM.jl")
include("../../../Routines/PRM/getMS1PeakHeights.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/buildPrecursorTable.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/initTransitions.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/selectTransitions.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/getBestTransitions.jl")
include("../../../SearchRAW.jl")
include("../../../Routines/PRM/IS-PRM_SURVEY/writeTables.jl")

#=
include("src/precursor.jl")
include("src/binaryRangeQuery.jl")
include("src/matchpeaks.jl")
include("src/getPrecursors.jl")
include("src/PSM_TYPES/PSM.jl")
include("src/PSM_TYPES/FastXTandem.jl")
#include("../../../searchSpectra.jl")
include("src/Routines/PRM/getBestPSMs.jl")
include("src/Routines/PRM/precursorChromatogram.jl")
include("src/Routines/PRM/plotPRM.jl")
include("src/Routines/PRM/getMS1PeakHeights.jl")
include("src/Routines/PRM/IS-PRM_SURVEY/initTransitions.jl")
include("src/Routines/PRM/IS-PRM_SURVEY/selectTransitions.jl")
include("src/Routines/PRM/IS-PRM_SURVEY/getBestTransitions.jl")
include("src/Routines/PRM/IS-PRM_SURVEY/buildPrecursorTable.jl")
include("src/SearchRAW.jl")
include("src/Routines/PRM/IS-PRM_SURVEY/writeTables.jl")
=#


##########
#Read Precursor Table
##########
@time begin 
    ptable = PrecursorTable()
    buildPrecursorTable!(ptable, 
                        params[:fixed_mods], 
                        params[:variable_mods], 
                        params[:max_variable_mods], 
                        PRECURSOR_LIST_PATH)
    addPrecursors!(
                        ptable, 
                        params[:precursor_charges], 
                        params[:precursor_isotopes], 
                        params[:modification_masses]
                        )
    addTransitions!(
                        ptable, 
                        params[:transition_charges],
                        params[:transition_isotopes],
                        params[:b_ladder_start],
                        params[:y_ladder_start],
                        params[:fragment_match_ppm]
                    )
    #println(getPrecIDToTransitions(ptable))
##########
#Search Survey Runs
##########
MS_TABLES = Dict{UInt32, Arrow.Table}()
combined_scored_psms = makePSMsDict(FastXTandem())
combined_fragment_matches = Dict{UInt32, Vector{FragmentMatch}}()
MS_RT = Dict{UInt32, Vector{Float32}}()
    for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)

        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        MS_RT[UInt32(ms_file_idx)] = MS_TABLE[:retentionTime]
        scored_psms, fragment_matches = SearchRAW(
                                                MS_TABLE, 
                                                ptable, 
                                                selectTransitions, 
                                                params[:right_precursor_tolerance],
                                                params[:left_precursor_tolerance],
                                                params[:transition_charges],
                                                params[:transition_isotopes],
                                                params[:b_ladder_start],
                                                params[:y_ladder_start],
                                                params[:fragment_match_ppm],
                                                UInt32(ms_file_idx)
                                                )
        for key in keys(combined_scored_psms)
            append!(combined_scored_psms[key], scored_psms[key])
        end
        combined_fragment_matches[UInt32(ms_file_idx)] = fragment_matches
    end

##########
#Get Best PSMs for Each Peptide
##########
    best_psms = getBestPSMs(combined_scored_psms, ptable, MS_RT, params[:minimum_fragment_count])
##########
#Get MS1 Peak Heights
##########
    #First key is ms_file_idx (identifier of the ms file), second key is pep_idx (peptide id)
    ms1_peak_heights = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, Float32}}()
    #Peak heights are zero to begin with
    precursor_idxs = unique(best_psms[!,:precursor_idx])
    for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        insert!(ms1_peak_heights, 
                UInt32(ms_file_idx), 
                UnorderedDictionary(precursor_idxs, zeros(Float32, length(precursor_idxs)))
                )

        getMS1PeakHeights!( ptable,
                            MS_TABLE[:retentionTime], 
                            MS_TABLE[:masses], 
                            MS_TABLE[:intensities], 
                            MS_TABLE[:msOrder], 
                            ms1_peak_heights[ms_file_idx], 
                            best_psms[!,:retention_time], 
                            best_psms[!,:precursor_idx], 
                            best_psms[!,:ms_file_idx],
                            Float32(0.25), 
                            params[:right_precursor_tolerance], 
                            params[:left_precursor_tolerance],
                            UInt32(ms_file_idx))
    end

    #Add MS1 Heights to the best_psms DataFrame 
    transform!(best_psms, AsTable(:) => ByRow(psm -> ms1_peak_heights[psm[:ms_file_idx]][psm[:precursor_idx]]) => :ms1_peak_height)
##########
#Get Chromatograms for the best precursors in each file. 
##########
    precursor_chromatograms = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, PrecursorChromatogram}}()
    for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        insert!(precursor_chromatograms, UInt32(ms_file_idx), initPrecursorChromatograms(best_psms, UInt32(ms_file_idx)) |> (best_psms -> fillPrecursorChromatograms!(best_psms, 
                                                                                                                    combined_fragment_matches[UInt32(ms_file_idx)], 
                                                                                                                    MS_TABLE, 
                                                                                                                    params[:precursor_rt_tolerance],
                                                                                                                    UInt32(ms_file_idx))
                                                                                                        )
                ) 
    end 

    #Names and charges for the "n" most intense fragment ions for each precursor
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestTransitions(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]]),
                                                                        params[:fragments_to_select])) => :best_transitions)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:name][psm[:best_transitions]]) => :transition_names)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)

##########
#Apply conditions to MS Files. 
##########
MS_FILE_ID_TO_NAME = Dict(
                            zip(
                            [UInt32(i) for i in 1:length(MS_TABLE_PATHS)], 
                            [splitpath(filepath)[end] for filepath in MS_TABLE_PATHS]
                            )
                        )

MS_FILE_ID_TO_CONDITION = Dict(
                                    zip(
                                    [key for key in keys(MS_FILE_ID_TO_NAME)], 
                                    ["NONE" for key in keys(MS_FILE_ID_TO_NAME) ]
                                    )
                                )

for (file_id, file_name) in MS_FILE_ID_TO_NAME 
    for (condition, value) in params[:ms_file_conditions]
        if occursin(condition, file_name)
            MS_FILE_ID_TO_CONDITION[file_id] = condition
        end
    end
end

transform!(best_psms, AsTable(:) => ByRow(psm -> MS_FILE_ID_TO_CONDITION[psm[:ms_file_idx]]) => :condition)
##########
#Write Method Files
##########
    #Get best_psm for each peptide across all ms_file
    best_psms = combine(sdf -> sdf[argmax(sdf.hyperscore), :], groupby(best_psms, :pep_idx)) 

    writeTransitionList(best_psms, joinpath(MS_DATA_DIR, "transition_list.csv"))
    writeIAPIMethod(best_psms, joinpath(MS_DATA_DIR, "iapi_method.csv"))

    println(" Scored "*string(size(best_psms)[1])*" precursors")
#end
end
##########
#Make Plots
##########
#=
if ARGS["make_plots"]
    for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        sample_name = split(MS_FILE_ID_TO_NAME[ms_file_idx], ".")[1]
        plotAllBestSpectra(precursor_chromatograms[ms_file_idx], 
                            ptable, 
                            MS_TABLE,
                            joinpath("./figures/best_spectra/", sample_name),
                            join(sample_name*"_best_spectra"*".pdf"))
        plotAllFragmentIonChromatograms(precursor_chromatograms[ms_file_idx], 
                                        ptable,
                                        joinpath("./figures/precursor_chromatograms/", sample_name),
                                        join(sample_name*"_precursor_chromatograms"*".pdf"))
    end
end=#
