##########
#Import Libraries
##########
#Data Parsing/Printing
println("Importing Libraries...")

using ArgParse
using CSV, Arrow, Tables, DataFrames, JSON, JLD2, ProgressBars
using Plots, StatsPlots, PrettyPrinting, CategoricalArrays
#DataStructures 
using DataStructures, Dictionaries, Distributions, Combinatorics, StatsBase, LinearAlgebra, Random, LoopVectorization, SparseArrays
#Algorithms 
using Interpolations, XGBoost, SavitzkyGolay, NumericalIntegration, ExpectationMaximization, LsqFit, FastGaussQuadrature, GLM, StaticArrays
using Base.Order
using Base.Iterators: partition
using PDFmerger

##########
#Parse Arguments 
##########
#Example Usage 
#julia --threads 24 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/Projects/PROSIT/TEST_DATA/MOUSE_TEST /Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423 -s true 
#julia --threads 9 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/TEST_DATA/mzXML/ /Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HumanYeastEcoli/5ppm_15irt/ -s true
#julia --threads 9 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/TEST_DATA/mzXML/ /Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_start2 -s true -e TEST_EXP
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "params_json"
            help = "Path to a .json file with the parameters"
            required = true
        "data_dir"
            help = "Path to a folder with .arrow MS data tables"
            required = true
        "spec_lib_dir"
            help = "Path to a tab delimited table of transitions"
            required = true
        "--experiment_name", "-e"
            help = "Name of subdirectory for output"
            default = "EXPERIMENT"
        "--print_params", "-s"
            help = "Whether to print the parameters from the json. Defaults to `false`"
            default = false 
    end

    return parse_args(s)
end

println("Parsing Arguments...")
ARGS = parse_commandline();

params = JSON.parse(read(ARGS["params_json"], String));
#=
3-protome PC test March 4th 2024 for hupo
params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));

SPEC_LIB_DIR = "/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/HEIL_2023/"
#SPEC_LIB_DIR =  "C:\\Users\\n.t.wamsley\\PROJECTS\\HUPO_2023\\HUMAN_YEAST_ECOLI\\PIONEER\\LIB"
#MS_DATA_DIR = "C:\\Users\\n.t.wamsley\\PROJECTS\\HUPO_2023\\HUMAN_YEAST_ECOLI\\PIONEER\\RAW"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "HUPO_THREE_PROTEOME_Mar4_2024"
=#

println("ARGS ", ARGS)
MS_DATA_DIR = ARGS["data_dir"];
EXPERIMENT_NAME = ARGS["experiment_name"];
SPEC_LIB_DIR = ARGS["spec_lib_dir"];

#Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATH_TO_ID = Dictionary(MS_TABLE_PATHS, UInt32.(collect(range(1,length(MS_TABLE_PATHS)))))
MS_TABLE_ID_TO_PATH = Dictionary(UInt32.(collect(range(1,length(MS_TABLE_PATHS)))), MS_TABLE_PATHS)

MS_DATA_DIR = joinpath(MS_DATA_DIR, EXPERIMENT_NAME);
out_folder = joinpath(MS_DATA_DIR, "Search")
if !isdir(out_folder)
    mkpath(out_folder)
end
qc_plot_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS")
if !isdir(qc_plot_folder)
    mkpath(qc_plot_folder)
end
const rt_alignment_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS","rt_alignment")
if !isdir(rt_alignment_folder)
    mkpath(rt_alignment_folder)
end
const mass_err_estimation_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS","mass_error_estimation")
if !isdir(mass_err_estimation_folder)
    mkpath(mass_err_estimation_folder)
end
const results_folder = joinpath(MS_DATA_DIR, "Search", "RESULTS")
if !isdir(results_folder)
    mkpath(results_folder)
end
const params_folder = joinpath(MS_DATA_DIR, "Search", "PARAMS")
if !isdir(params_folder )
    mkpath(params_folder)
end

params_ = (
    expected_matches = Int64(params["expected_matches"]),
    isotope_err_bounds = Tuple([Int64(bound) for bound in params["isotope_err_bounds"]]),
    choose_most_intense = Bool(params["choose_most_intense"]),
    quadrupole_isolation_width = Float64(params["quadrupole_isolation_width"]),
    irt_err_sigma = params["irt_err_sigma"],

    presearch_params = Dict{String, Any}(k => v for (k, v) in params["presearch_params"]),
    first_search_params = Dict{String, Any}(k => v for (k, v) in params["first_search_params"]),
    quant_search_params = Dict{String, Any}(k => v for (k, v) in params["quant_search_params"]),
    frag_tol_params = Dict{String, Any}(k => v for (k, v) in params["frag_tol_params"]),
    irt_mapping_params = Dict{String, Any}(k => v for (k, v) in params["irt_mapping_params"]),
    integration_params = Dict{String, Any}(k => v for (k, v) in params["integration_params"]),
    deconvolution_params = Dict{String, Any}(k => v for (k, v) in params["deconvolution_params"]),
    summarize_first_search_params = Dict{String, Any}(k => v for (k, v) in params["summarize_first_search_params"]),
    qc_plot_params = Dict{String, Any}(k => v for (k, v) in params["qc_plot_params"])

    );

##########
#Load Dependencies 
##########
#Fragment Library Parsing

[include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in [
                                                                    "ArrayDict.jl",
                                                                    "Counter.jl",
                                                                    "Ion.jl",
                                                                    "LibraryIon.jl",
                                                                    "MatchIon.jl",
                                                                    "LibraryFragmentIndex.jl",
                                                                    "SparseArray.jl"]];

#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in [
                                                                    "ExponentialGaussianHybrid.jl",
                                                                    "isotopes.jl",
                                                                    "globalConstants.jl",
                                                                    "isotopeSplines.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "SpectralDeconvolution.jl",
                                                                    "percolatorSortOf.jl",
                                                                    "kdeRTAlignment.jl",
                                                                    "probitRegression.jl",
                                                                    "partitionThreadTasks.jl"]];

[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchpeaksLib.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl",
                                                                                    "integrateChroms.jl"]];
                                             

library_fragment_lookup_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"lib_frag_lookup_031424", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index_top5_7ppm_2hi_031424", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors_031424", file) != nothing, readdir(SPEC_LIB_DIR))][1]

println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
spec_load_time = @timed begin
    @time const f_index = load(f_index_path)["f_index"];
    prosit_lib["f_index"] = f_index;#["f_index"]
    @time const library_fragment_lookup_table = load(library_fragment_lookup_path)["lib_frag_lookup"]
    prosit_lib["f_det"] = library_fragment_lookup_table; #["f_det"];
    @time const precursors = load(precursors_path)["precursors"]
    prosit_lib["precursors"] = precursors;#["precursors"];
end

###########
#Load Pre-Allocated Data Structures. One of each for each thread. 
###########
@time begin
N = Threads.nthreads()
ionMatches = [[FragmentMatch{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
ionMisses = [[FragmentMatch{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
all_fmatches = [[FragmentMatch{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
IDtoCOL = [ArrayDict(UInt32, UInt16, length(precursors)) for _ in range(1, N)];
ionTemplates = [[DetailedFrag{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml");
scored_PSMs = [Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
unscored_PSMs = [[SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
spectral_scores = [Vector{SpectralScoresSimple{Float16}}(undef, 5000) for _ in range(1, N)];
precursor_weights = [zeros(Float32, length(precursors)) for _ in range(1, N)];
precs = [Counter(UInt32, UInt8,length(precursors)) for _ in range(1, N)];
complex_scored_PSMs = [Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
complex_unscored_PSMs = [[ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
complex_spectral_scores = [Vector{SpectralScoresComplex{Float16}}(undef, 5000) for _ in range(1, N)];
end;


###########
#File Names Parsing 
###########
file_names = first.(split.(basename.(MS_TABLE_PATHS), '.'))
split_file_names = split.(file_names, "_")



#If file names have inconsistnat number of delimiters, give up on parsing and use the entire file name
unique_split_file_name_lengths = unique(length.(split_file_names))
if unique_split_file_name_lengths == 1
    N = first(unique_split_file_name_lengths)
    M = length(file_names)
    split_strings = Array{String}(undef, (M, N))
    for i in 1:N
        for j in 1:M
            split_strings[j, i] = split_file_names[j][i]
        end
    end
    cols_to_keep = zeros(Bool, N)
    for (col_idx, col) in enumerate(eachcol(split_strings))
        if length(unique(col))>1
            cols_to_keep[col_idx] = true
        end
    end
    parsed_file_names = Vector{String}(undef, M)
    split_strings = split_strings[:, cols_to_keep]
    for i in 1:M
        parsed_file_names[i] = join(split_strings[i,:], "_")
    end
else
    parsed_file_names = file_names
end

const file_id_to_parsed_name = Dict(zip(1:M, parsed_file_names))
#Parsed file names 
const parsed_fnames = sort(collect(values(file_id_to_parsed_name)))
