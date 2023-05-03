using DataFrames
include("../src/precursor.jl")
include("../src/matchpeaks.jl")
include("../src/PSM_TYPES/FastXTandem.jl")
include("../src/getPrecursors.jl")
include("../src/binaryRangeQuery.jl")
include("precursor.jl")

include("../src/applyMods.jl")
include("applyMods.jl")
include("../src/Routines/PRM/IS-PRM-SURVEY/buildPrecursorTable.jl")
#include("./Routines/PRM/IS-PRM-SURVEY/buildPrecursorTable.jl")
include("Routines/PRM/IS-PRM-SURVEY/buildPrecursorTable.jl")

include("../src/Routines/PRM/IS-PRM/buildPrecursorTable.jl")
#include("./Routines/PRM/IS-PRM/buildPrecursorTable.jl")
include("Routines/PRM/IS-PRM/buildPrecursorTable.jl")

include("matchpeaks.jl")
include("binaryRangeQuery.jl")
include("PSM_TYPES/FastXTandem.jl")

include("../src/Routines/PRM/getMS1PeakHeights.jl")
#include("./Routines/PRM/getMS1PeakHeights.jl")
include("Routines/PRM/getMS1PeakHeights.jl")

include("../src/Routines/PRM/IS-PRM-SURVEY/getBestTransitions.jl")
#include("./Routines/PRM/IS-PRM-SURVEY/getBestTransitions.jl")
include("Routines/PRM/IS-PRM-SURVEY/getBestTransitions.jl")

include("../src/Routines/PRM/IS-PRM/getIntegrationBounds.jl")
#include("./Routines/PRM/IS-PRM/getIntegrationBounds.jl")
include("Routines/PRM/IS-PRM/getIntegrationBounds.jl")

include("../src/LFQ.jl")
include("LFQ.jl")

include("../src/Routines/PRM/precursorChromatogram.jl")
include("../src/Routines/PRM/IS-PRM/parEstimation.jl")

include("../src/SearchRAW.jl")
include("../src/Routines/PRM/IS-PRM/selectTransitions.jl")

include("../src/Routines/PRM/IS-PRM/getBestPSMs.jl")
##########
#Test Routine
##########
using Arrow, JSON, Tables, DataFrames, Plots
    params = JSON.parse(read("../data/test.json", String))
    function parse_mods(fixed_mods)
        fixed_mods_parsed = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
        for mod in fixed_mods
            push!(fixed_mods_parsed, (p=Regex(mod[1]), r = mod[2]))
        end
        return fixed_mods_parsed
    end
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
    MS_TABLE_PATHS = [ "../data/parquet/GAPDH_VGVNGFGR.arrow",
                       "../data/parquet/GSTM1_RPWFAGNK.arrow",
                       "../data/parquet/GSTM4_VAVWGNK.arrow"]
    TRANSITION_LIST_PATH = "../data/parquet/transition_list.csv"
    ptable = ISPRMPrecursorTable()
    buildPrecursorTable!(ptable, mods_dict, TRANSITION_LIST_PATH)
    combined_scored_psms = makePSMsDict(FastXTandem())
    combined_fragment_matches = Dict{UInt32, Vector{FragmentMatch}}()
    MS_RT = Dict{UInt32, Vector{Float32}}()
    #println("START")
    lk = ReentrantLock()
    #MS_TABLE_PATHS = ["/Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out/MA4953_SQ_Kinome_H358_Rep1_EWZ_Rerun.arrow","/Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out/MA5032_SQ_Kinome_HCC827_Rep1_EWZ.arrow"]
    for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        println("MS_TABLE_PATH ", MS_TABLE_PATH)
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
        lock(lk) do 
            for key in keys(combined_scored_psms)
                append!(combined_scored_psms[key], scored_psms[key])
            end
            combined_fragment_matches[UInt32(ms_file_idx)] = fragment_matches
        end
    end             
    best_psms = getBestPSMs(combined_scored_psms, ptable, MS_RT, UInt8(1))
    precursor_chromatograms = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, PrecursorChromatogram}}()
    for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        lock(lk) do 
        insert!(precursor_chromatograms, UInt32(ms_file_idx), initPrecursorChromatograms(best_psms, UInt32(ms_file_idx)) |> (best_psms -> fillPrecursorChromatograms!(best_psms, 
                                                                                                                    combined_fragment_matches[UInt32(ms_file_idx)], 
                                                                                                                    MS_TABLE, 
                                                                                                                    params[:precursor_rt_tolerance],
                                                                                                                    UInt32(ms_file_idx))
                                                                                                        )
                ) 
        end
    end 
    using Statistics, RobustModels
    for i in collect(eachindex(MS_TABLE_PATHS))
        println(i)
        getPARs(ptable, precursor_chromatograms[i], minimum_scans = 3, ms_file_idx = UInt32(i))
    end
    transform!(best_psms, AsTable(:) => ByRow(psm -> getPAR(ptable, psm[:precursor_idx], psm[:ms_file_idx])) => [:par, :goodness_of_fit, :isotope])
    quant = select(best_psms[(best_psms.isotope .== "light"),:], [:ms_file_idx, :sequence, :protein_names, :par, :goodness_of_fit])
    filter!(row -> (coalesce(row.par, 0.0) > 1e-8) && (coalesce(row.goodness_of_fit, 0.0) < 0.10), quant)
    ##########
    include("./Routines/PRM/IS-PRM/parEstimation.jl")

    include("../src/Routines/PRM/IS-PRM-SURVEY/getBestPSMs.jl")
    #include("./Routines/PRM/IS-PRM-SURVEY/getBestPSMs.jl")
    include("Routines/PRM/IS-PRM-SURVEY/getBestPSMs.jl")
#include("./PSM_TYPES/FastXTandem.jl")
#include("getPrecursors.jl")