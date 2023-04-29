#julia ./src/Routines/PRM/IS-PRM/routine.jl ./data/test.json /Users/n.t.wamsley/Desktop/test_folder ./data/NRF2_SIL.txt\
#julia ./src/Routines/PRM/IS-PRM/routine.jl ./data/IS-PRM_TEST.json /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/survey_runs/transition_list.csv
#julia --threads 24 ./src/Routines/PRM/IS-PRM/routine.jl ./data/IS-PRM_TEST.json /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out /Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/survey_runs/transition_list.csv
using JSON
using PrettyPrinting
using PDFmerger
using CSV
using Arrow, DataFrames, Tables, Plots
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
        "transition_list"
            help = "Path to a tab delimited table of transitions"
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
TRANSITION_LIST_PATH = ARGS["transition_list"]
#MS_DATA_DIR = "./data/parquet/"
#PRECURSOR_LIST_PATH = "./data/NRF2_SIL.txt"
#Get all files in the `MS_DATA_DIR` ending in ".arrow" and append their names to the `MS_DATA_DIR` path. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))]

#=
MS_TABLE_PATHS = [ "./data/parquet/GAPDH_VGVNGFGR.arrow",
"./data/parquet/GSTM1_RPWFAGNK.arrow",
"./data/parquet/GSTM4_VAVWGNK.arrow"]
=#
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
include("../../../precursor.jl")
include("../../../binaryRangeQuery.jl")
include("../../../matchpeaks.jl")
include("../../../getPrecursors.jl")
include("../../../PSM_TYPES/PSM.jl")
include("../../../PSM_TYPES/FastXTandem.jl")
include("../../../Routines/PRM/IS-PRM/getBestPSMs.jl")
include("../../../Routines/PRM/precursorChromatogram.jl")
include("../../../Routines/PRM/getMS1PeakHeights.jl")
include("../../../Routines/PRM/IS-PRM/initTransitions.jl")
include("../../../Routines/PRM/IS-PRM/getBestTransitions.jl")
include("../../../SearchRAW.jl")
include("../../../Routines/PRM/IS-PRM/buildPrecursorTable.jl")
include("../../../Routines/PRM/IS-PRM/selectTransitions.jl")
include("../../../Routines/PRM/IS-PRM/getBestPSMs.jl")
include("../../../Routines/PRM/IS-PRM/getIntegrationBounds.jl")
include("../../../Routines/PRM/IS-PRM/parEstimation.jl")
include("../../../LFQ.jl")
include("../../../Routines/PRM/plotPRM.jl")
##########
#Read Precursor Table
##########
#@time begin 
ptable = ISPRMPrecursorTable()
buildPrecursorTable!(ptable, mods_dict, TRANSITION_LIST_PATH)
#println(ptable.pepGroup_to_id)
scan_adresses = Dict{UInt32, Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}}()
MS_TABLES = Dict{UInt32, Arrow.Table}()
println("N threads ", Threads.nthreads())
for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    #println("MS_TABLE_PATH ", MS_TABLE_PATH)
    scan_adresses[UInt32(ms_file_idx)] = getScanAdresses(MS_TABLE[:msOrder])
    #setIntegrationBounds!(ptable, ms_file_idx, scan_adresses[UInt32(ms_file_idx)], MS_TABLES[UInt32(ms_file_idx)][:precursorMZ])
end

##########
#Search Survey Runs
##########

#MS_TABLES = Dict{UInt32, Arrow.Table}()
combined_scored_psms = makePSMsDict(FastXTandem())
combined_fragment_matches = Dict{UInt32, Vector{FragmentMatch}}()
MS_RT = Dict{UInt32, Vector{Float32}}()
println("START")
lk = ReentrantLock()
#MS_TABLE_PATHS = ["/Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out/MA4967_SQ_Kinome_KYSE70_Rep3_EWZ.arrow",
#"/Users/n.t.wamsley/RIS_Temp/EWZ_KINOME/parquet_out/MA5032_SQ_Kinome_HCC827_Rep1_EWZ.arrow"]
@time begin
Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
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
end
println("SEARCHED")
##########
#Get Best PSMs for Each Peptide
##########
    @time begin
        best_psms = getBestPSMs(combined_scored_psms, ptable, MS_RT, UInt8(1))
    end
##########
#Get MS1 Peak Heights
##########
    #First key is ms_file_idx (identifier of the ms file), second key is pep_idx (peptide id)
    ms1_peak_heights = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, Float32}}()
    #Peak heights are zero to begin with
    precursor_idxs = unique(best_psms[!,:precursor_idx])
    @time begin
    Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        lock(lk) do 
        insert!(ms1_peak_heights, 
                UInt32(ms_file_idx), 
                UnorderedDictionary(precursor_idxs, zeros(Float32, length(precursor_idxs)))
                )
        end
        getMS1PeakHeights!(ms1_peak_heights[ms_file_idx],
             ptable,
                            MS_TABLE[:retentionTime], 
                            MS_TABLE[:masses], 
                            MS_TABLE[:intensities], 
                            MS_TABLE[:msOrder], 
                             
                            best_psms[!,:retention_time], 
                            best_psms[!,:precursor_idx], 
                            best_psms[!,:ms_file_idx],
                            Float32(0.25), 
                            params[:right_precursor_tolerance], 
                            params[:left_precursor_tolerance],
                            UInt32(ms_file_idx))
    end
    end
    #Add MS1 Heights to the best_psms DataFrame 
    transform!(best_psms, AsTable(:) => ByRow(psm -> ms1_peak_heights[psm[:ms_file_idx]][psm[:precursor_idx]]) => :ms1_peak_height)
    
##########
#Get Chromatograms for the best precursors in each file. 
##########
    precursor_chromatograms = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, PrecursorChromatogram}}()
    Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
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

    #Names and charges for the "n" most intense fragment ions for each precursor
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestTransitions(getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]]),
                                                                        params[:fragments_to_select])) => :best_transitions)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:name][psm[:best_transitions]]) => :transition_names)
    transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)
    
##########
#Get PARs. 
##########
using GLMNet
using Statistics
println("start getPARs")
@time begin
    #=for i in collect(eachindex(MS_TABLE_PATHS))
        println(i)
        getPARs(ptable, scan_adresses[i], precursor_chromatograms[i], minimum_scans = 3, ms_file_idx = UInt32(i))
    end=#
    #Threads.@threads 
    for i in collect(eachindex(MS_TABLE_PATHS))
        println(i)
        getPARs(ptable, scan_adresses[i], precursor_chromatograms[i], minimum_scans = 3, ms_file_idx = UInt32(i))
    end
end
println("end getPARs")
#=for i in ptable.lh_pair_id_to_light_heavy_pair
    if i.light_sequence == "AANLLLGYK"
        println(i)
        println(precursor_chromatograms[UInt32(2)][0x000004dd])
        println(precursor_chromatograms[UInt32(2)][0x000004de])
    end
end=#
transform!(best_psms, AsTable(:) => ByRow(psm -> getPAR(ptable, psm[:precursor_idx], psm[:ms_file_idx])) => [:par, :dev_ratio, :isotope])

#Get the peak area ratios into the dataframe. 
#Then make the plots with the PARs.
#Then figure out how to summarize to protien level by LFQ and make a protien table

##########
#Write Protein Peptide Quant Table
##########
using StatsBase
out = Dict(
    :protein => String[],
    :peptides => String[],
    :log2_abundance => Float64[],
    :experiments => UInt32[],
)

quant = select(best_psms[(best_psms.isotope .== "light"),:], [:ms_file_idx, :sequence, :protein_names, :par, :isotope, :dev_ratio])
println("start LFQ")
out_folder = joinpath(MS_DATA_DIR, "tables/")
if !isdir(out_folder)
    mkpath(out_folder)
end
MS_FILE_ID_TO_NAME = Dict(
                            zip(
                            [UInt32(i) for i in 1:length(MS_TABLE_PATHS)], 
                            [splitpath(filepath)[end] for filepath in MS_TABLE_PATHS]
                            )
                        )
transform!(quant, AsTable(:) => ByRow(psm -> MS_FILE_ID_TO_NAME[psm[:ms_file_idx]]) => :file_name)
CSV.write(joinpath(out_folder, "peptide.csv"), quant)
# Filter the rows based on the conditions on featurea and featureb
filter!(row -> (coalesce(row.par, 0.0) > 1e-12) && (coalesce(row.dev_ratio, 0.0) < 0.10), quant)
@time begin
    for (protein, data) in pairs(groupby(quant, :protein_names))

        getProtAbundance(string(protein[:protein_names]), 
                            collect(data[!,:sequence]), 
                            collect(data[!,:ms_file_idx]), 
                            (collect(data[!,:par])),
                            out[:protein],
                            out[:peptides],
                            out[:experiments],
                            out[:log2_abundance]
                        )
        #println(protein[:parent])
    end
end

out = DataFrame(out)

transform!(out, AsTable(:) => ByRow(psm -> MS_FILE_ID_TO_NAME[psm[:experiments]]) => :file_name)
#out_folder = joinpath(MS_DATA_DIR, "tables/")
if !isdir(out_folder)
    mkpath(out_folder)
end
CSV.write(joinpath(out_folder, "protein.csv"), out)
##########
#Make Plots
##########
#end
for (i, path) in collect(enumerate(MS_TABLE_PATHS))
    sample_name = split(splitpath(path)[end], ".")[1]
    plotAllPairedFragmentIonChromatograms(ptable, precursor_chromatograms[i], UInt32(i), joinpath(MS_DATA_DIR,"figures/paired_spectra",sample_name),
    join(sample_name*"_precursor_chromatograms"*".pdf"))
    println(sample_name)
end
#end