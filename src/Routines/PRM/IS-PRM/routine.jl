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
include("../../../Routines/PRM/plotPRM.jl")
include("../../../Routines/PRM/getMS1PeakHeights.jl")
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
include("src/Routines/PRM/IS-PRM/initTransitions.jl")
include("src/Routines/PRM/IS-PRM/getBestTransitions.jl")
include("src/SearchRAW.jl")
=#
include("src/Routines/PRM/IS-PRM/buildPrecursorTable.jl")
include("src/Routines/PRM/IS-PRM/selectTransitions.jl")
include("src/Routines/PRM/IS-PRM/getBestPSMs.jl")
include("src/Routines/PRM/IS-PRM/getIntegrationBounds.jl")
include("src/Routines/PRM/IS-PRM/parEstimation.jl")
##########
#Read Precursor Table
##########
@time begin 
    ptable = ISPRMPrecursorTable()
    buildPrecursorTable!(ptable, mods_dict, "data/parquet/transition_list.csv")
scan_adresses = Dict{UInt32, Vector{NamedTuple{(:scan_index, :ms1, :msn), Tuple{Int64, Int64, Int64}}}}()
MS_TABLES = Dict{UInt32, Arrow.Table}()
MS_TABLE_PATHS = [ "./data/parquet/GAPDH_VGVNGFGR.arrow",
"./data/parquet/GSTM1_RPWFAGNK.arrow",
"./data/parquet/GSTM4_VAVWGNK.arrow"]
for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)
    MS_TABLES[UInt32(ms_file_idx)] = Arrow.Table(MS_TABLE_PATH)
    scan_adresses[UInt32(ms_file_idx)] = getScanAdresses(MS_TABLES[UInt32(ms_file_idx)][:msOrder])
    #setIntegrationBounds!(ptable, ms_file_idx, scan_adresses[UInt32(ms_file_idx)], MS_TABLES[UInt32(ms_file_idx)][:precursorMZ])
end
##########
#Search Survey Runs
##########

MS_TABLES = Dict{UInt32, Arrow.Table}()
combined_scored_psms = makePSMsDict(FastXTandem())
combined_fragment_matches = Dict{UInt32, Vector{FragmentMatch}}()
    for (ms_file_idx, MS_TABLE_PATH) in enumerate(MS_TABLE_PATHS)

        MS_TABLES[UInt32(ms_file_idx)] = Arrow.Table(MS_TABLE_PATH)

        scored_psms, fragment_matches = SearchRAW(
                                                MS_TABLES[UInt32(ms_file_idx)], 
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
    best_psms = getBestPSMs(combined_scored_psms, ptable, MS_TABLES, UInt8(1))
 
##########
#Get MS1 Peak Heights
##########
    #First key is ms_file_idx (identifier of the ms file), second key is pep_idx (peptide id)
    ms1_peak_heights = UnorderedDictionary{UInt32, UnorderedDictionary{UInt32, Float32}}()
    #Peak heights are zero to begin with
    precursor_idxs = unique(best_psms[!,:precursor_idx])
    for (ms_file_idx, MS_TABLE) in MS_TABLES
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
    for (ms_file_idx, MS_TABLE) in MS_TABLES

        insert!(precursor_chromatograms, UInt32(ms_file_idx), initPrecursorChromatograms(best_psms, UInt32(ms_file_idx)) |> (best_psms -> fillPrecursorChromatograms!(best_psms, 
                                                                                                                    combined_fragment_matches[UInt32(ms_file_idx)], 
                                                                                                                    MS_TABLE, 
                                                                                                                    params[:precursor_rt_tolerance],
                                                                                                                    UInt32(ms_file_idx))
                                                                                                        )
                ) 
    end 
    
    #scan_adresses[1][precursor_chromatograms[1][183].last_scan_idx[2:end]]
    #precursor_chromatograms[1][184].last_scan_idx
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
for i in eachindex(MS_TABLES)
    getPARs(ptable, scan_adresses[i], precursor_chromatograms[i], minimum_scans = 3, ms_file_idx = UInt32(i))
end
#Get the peak area ratios into the dataframe. 
#Then make the plots with the PARs.
#Then figure out how to summarize to protien level by LFQ and make a protien table

##########
#Write Protein Peptide Quant Table
##########
#transform!(best_psms, AsTable(:) => ByRow(psm -> getBestPSM(precursor_chromatograms[psm[:ms_file_idx]][psm[:precursor_idx]])[:mz][psm[:best_transitions]]) => :transition_mzs)

function getPAR(ptable::ISPRMPrecursorTable, prec_id::UInt32, ms_file_idx::UInt32)
    lh_pair = ptable.lh_pair_id_to_light_heavy_pair[ptable.prec_id_to_lh_pair_id[prec_id]]
    isotope = "light"
    if lh_pair.heavy_prec_id == prec_id
        isotope = "heavy"
    end
    if length(lh_pair.par_model) == 0
        return (Missing, Missing, isotope)
    elseif !isassigned(lh_pair.par_model, ms_file_idx)
        return (missing, missing, isotope)
    else
        return (1/lh_pair.par_model[ms_file_idx].par_model_coef[1],
                lh_pair.par_model[ms_file_idx].dev_ratio,
                isotope
                )
    end
end
transform!(best_psms, AsTable(:) => ByRow(psm -> (getPAR(ptable, psm[:precursor_idx], psm[:ms_file_idx]))) => [:par, :dev_ratio, :isotope])

quant = select(best_psms[(best_psms.isotope .== "light"),:], [:ms_file_idx, :sequence, :protein_names, :par, :isotope, :dev_ratio])
getProtAbundance(prot[!, :peptide], prot[!, :file_idx], prot[!, :abundance])
combine(sdf -> getProtAbundance(sdf[!,:sequence], sdf[!,:ms_file_idx], sdf[!, :par]), groupby(quant, :protein_names))
combine(sdf -> typeof(sdf[!,:sequence]), groupby(quant, :protein_names))


quant = select(best_psms[(best_psms.isotope .== "light"),:], [:ms_file_idx, :sequence, :protein_names, :par, :isotope, :dev_ratio])
combine(groupby(quant, :protein_names), [:sequence, :ms_file_idx, :par] => (sequence, ms_file_idx, par) -> getProtAbundance(sequence, ms_file_idx, par))
a = groupby(quant, :protein_names)
for group in grouped
    getProtAbundance(group.sequence,group.ms_file_idx, group.par)
end
##########
#Make Plots
##########

if ARGS["make_plots"]
    for (ms_file_idx, MS_TABLE) in MS_TABLES
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
end

prot = DataFrame(Dict(
    :peptide => ["A","A","A","B","B","B","C","C","C"],
    :protein => split(repeat("A",9), ""),
    :file_idx => [1, 2, 3, 1, 2, 3, 1, 2, 3],
    :abundance => [10, 12, 14, 1, 1.2, 1.4, 100, 120, 140],
))

prot = DataFrame(Dict(
    :peptide => ["A","A","A","B","B","B","C","C","C"],
    :protein => split(repeat("A",9), ""),
    :file_idx => [1, 2, 3, 1, 2, 3, 1, 2, 3],
    :abundance => [10, 20, 40, 1, 2, 4, 100, 200, 400],
))

prot = DataFrame(Dict(
    :peptide => ["A","A","A","B","B","B","C","C","C"],
    :protein => split(repeat("A",9), ""),
    :file_idx => [1, 2, 3, 1, 2, 3, 1, 2, 3],
    :abundance => [10, 20, 40, 1, 2, 4, 100, 200, missing],
))

prot = DataFrame(Dict(
    :peptide => ["A"],
    :protein => split(repeat("A",1), ""),
    :file_idx => [1],
    :abundance => [10],
))


function getS(peptides::Vector{String}, experiments::Vector{Int}, abundance::Vector{Union{T, Missing}}) where T <: Real
    S = Matrix(zeros(length(unique(peptides)), length(unique(experiments))))
    peptides_dict = Dict(zip(unique(peptides), 1:length(unique(peptides))))
    experiments_dict = Dict(zip(unique(experiments), 1:length(unique(experiments))))
    for i in eachindex(peptides)
        S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = coalesce(abundance[i], 0.0)
    end
    return S
end

function getB(S::Matrix{Float64}, N::Int)
    B = zeros(N + 1)
    for i in 1:N
        for j in (i+1):N
                r_i_j = median(-log2.(S[:,i]) + log2.(S[:,j]))
                if (S[i, j] != 0.0)
                    B[i] = B[i] - r_i_j# M[i, j]
                    B[j] = B[j] + r_i_j
                end
        end 
    end
    B.*=2
    B[end] = sum(S[S.!=0])*N/length(S) #Normalizing factor
    B
end

function getA(N::Int)
    A = ones(N+1, N+1)
    for i in 1:(N)
        for j in 1:N
            if i == j
                A[i, j] = 2*(N - 1)
            else
                A[i, j] = -1*(N-1)
            end
        end
    end
    A[end, end] = 0
    A
end

function getProtAbundance(peptides::Vector{String}, experiments::Vector{Int}, abundance::Vector{Union{T, Missing}}) where T <: Real
    N = length(unique(experiments))
    S = getS(peptides, experiments, abundance)
    B = getB(S, N)
    A = getA(N)
    return (A\B)[1:(end - 1)]
end

function getProtAbundance(peptides::Vector{String}, experiments::Vector{Int}, abundance::Vector{T}) where T <: Real
    getProtAbundance(peptides, experiments, allowmissing(abundance))
end 

S = getS(prot[!, :peptide], prot[!, :file_idx], prot[!, :abundance])
N = 3
A = ones(N+1, N+1)
B = zeros(3 + 1)
N = 3

getProtAbundance(prot[!, :peptide], prot[!, :file_idx], prot[!, :abundance])

prot = DataFrame(Dict(
    :peptide => ["A","A","A","B","B","B","C","C"],
    :protein => split(repeat("A",8), ""),
    :file_idx => [1, 2, 3, 1, 2, 3, 1, 2],
    :abundance => [10, 12, 14, 1, 1.2, 1.4, 100, 120],
))

A*log.([10, 12, 14])
B

for prot in grouped
    getProtAbundance(prot[!, :sequence], prot[!, :ms_file_idx], prot[!, :par])
end
A, B = getProtAbundance(prot[!, :peptide], prot[!, :file_idx], prot[!, :abundance])
sum(A*log.([10, 12, 14]) .- B)
a = (-1)*ones(3, 3)
a[diagind(a)] .= 4
Diagonal(4*ones(3, 3))

function my_function(x::Vector{Int})
    # function implementation
end

df = DataFrame(group = [1, 1, 2, 2], values = [1, 2, 3, 4])
result = combine(groupby(df, :group), :values => (x -> my_function(x)) => :result)
