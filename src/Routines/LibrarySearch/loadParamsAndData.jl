##########
#Import Libraries
##########
#Data Parsing/Printing
println("Importing Libraries...")

using ArgParse
using CSV, Arrow, Tables, DataFrames, JSON, JLD2, ProgressBars
using Plots, StatsPlots, PrettyPrinting, CategoricalArrays
#DataStructures 
using DataStructures, Dictionaries, Distributions, Combinatorics, StatsBase, LinearAlgebra, Random, LoopVectorization, SparseArrays, StaticArrays
#Algorithms 
using Interpolations, BSplineKit, Interpolations, XGBoost, SavitzkyGolay, NumericalIntegration, ExpectationMaximization, LsqFit, FastGaussQuadrature, GLM, LinearSolve
using Base.Order
using Base.Iterators: partition
using PDFmerger, Measures

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
library_fragment_lookup_table.prec_frag_ranges[end] = 0x29004baf:(0x29004be8 - 1)
3-protome PC test March 4th 2024 for hupo
julia --threads 15 ./src/Routines/LibrarySearch/MAIN.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/exploris_test /Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/UP000005640_9606_Apr20_24/pioneer_lib -s true 

params = JSON.parse(read("../OAT_103ISO/LibrarySearch.json", String));

SPEC_LIB_DIR ="../spec_lib/pioneer_lib"

MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/astral_test"


MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "EXPERIMENT"
=#

MS_DATA_DIR = ARGS["data_dir"];
EXPERIMENT_NAME = ARGS["experiment_name"];
SPEC_LIB_DIR = ARGS["spec_lib_dir"];

#Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATH_TO_ID = Dictionary(MS_TABLE_PATHS, UInt32.(collect(range(1,length(MS_TABLE_PATHS)))))
MS_TABLE_ID_TO_PATH = Dictionary(UInt32.(collect(range(1,length(MS_TABLE_PATHS)))), MS_TABLE_PATHS)

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
    qc_plot_params = Dict{String, Any}(k => v for (k, v) in params["qc_plot_params"]),
    normalization_params = Dict{String, Any}(k => v for (k, v) in params["normalization_params"]),
    benchmark_params = Dict{String, Any}(k => v for (k, v) in params["benchmark_params"])
    );


MS_DATA_DIR = joinpath(params_[:benchmark_params]["results_folder"], EXPERIMENT_NAME);
out_folder = joinpath(MS_DATA_DIR, "Search")
RESULTS_DIR = MS_DATA_DIR
if !isdir(out_folder)
    mkpath(out_folder)
end
qc_plot_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS")
if !isdir(qc_plot_folder)
    mkpath(qc_plot_folder)
else
    [rm(joinpath(qc_plot_folder, x)) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")]
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
#Write params to folder 
open(joinpath(params_folder, "config.json"),"w") do f 
    JSON.print(f, params)
end



##########
#Load Dependencies 
##########
#Fragment Library Parsing

[include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in [
                                                                    "ChromObject.jl",
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
                                                                    "uniformBasisCubicSpline.jl",
                                                                    "isotopeSplines.jl",
                                                                    "Normalization.jl",
                                                                    "SavitskyGolay.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "SpectralDeconvolution.jl",
                                                                    "percolatorSortOf.jl",
                                                                    "plotRTAlignment.jl",
                                                                    "probitRegression.jl",
                                                                    "partitionThreadTasks.jl"]];

[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl",
                                                                                    "integrateChroms.jl"]];
                                             

#library_fragment_lookup_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"lib_frag_lookup_031424", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index_top5_7ppm_1hi_rtmajor_031924", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors_031424", file) != nothing, readdir(SPEC_LIB_DIR))][1]

#println("Loading spectral libraries into main memory...")
#prosit_lib = Dict{String, Any}()
#@time detailed_frags = load("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_detailed_frags_032124.jld2")["detailed_frags"]
#@time prec_frag_ranges = load("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_prec_frag_ranges_032124.jld2")["prec_frag_ranges"]
#const library_fragment_lookup_table = LibraryFragmentLookup(detailed_frags, prec_frag_ranges)
#prosit_lib["f_det"] = library_fragment_lookup_table

#@time begin
f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
detailed_frags = load(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))["detailed_fragments"]
prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["precursor_to_fragment_indices"]
const library_fragment_lookup_table = LibraryFragmentLookup(detailed_frags, prec_frag_ranges)
last_range = library_fragment_lookup_table.prec_frag_ranges[end] #0x29004baf:(0x29004be8 - 1)
last_range = range(first(last_range), last(last_range) - 1)
library_fragment_lookup_table.prec_frag_ranges[end] = last_range
prosit_lib["f_det"] = library_fragment_lookup_table

#@time begin
#f_index_fragments = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_fragments_032124.arrow")
#f_index_rt_bins = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_rtbins_032124.arrow")
#f_index_frag_bins = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_fragbins_032124.arrow")


#precursor_bins = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_precursor_bins_031924.arrow");
#frag_bin_mzs = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_frag_bin_mzs_031924.arrow");

#detailed_frags = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_detailed_frags_031924.arrow")
#prec_frag_ranges = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_prec_frag_ranges_031924.arrow")
#Arrow.write("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/precursors_031924.arrow", DataFrame(precursors))
#Arrow.write("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/f_index_top5_7ppm_1hi_rtmajor_fragbins_031924.arrow", f_index_frag_bins)
#precursors = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/precursors_032124.arrow")#DataFrame(precursors)
precursors = Arrow.Table(joinpath(SPEC_LIB_DIR, "precursor_table.arrow"))#DataFrame(precursors)
f_index = FragmentIndex(
    f_index_frag_bins[:FragIndexBin],
    f_index_rt_bins[:FragIndexBin],
    f_index_fragments[:IndexFragment],
);
presearch_f_index = FragmentIndex(
    presearch_f_index_frag_bins[:FragIndexBin],
    presearch_f_index_rt_bins[:FragIndexBin],
    presearch_f_index_fragments[:IndexFragment],
);
prosit_lib["f_index"] = f_index;
prosit_lib["presearch_f_index"] = presearch_f_index;
prosit_lib["precursors"] = precursors;
###########
#Load Pre-Allocated Data Structures. One of each for each thread. 
###########
N = Threads.nthreads()
M = 250000
n_precursors = length(precursors[:mz])
ionMatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
ionMisses = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
all_fmatches = [[FragmentMatch{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors ) for _ in range(1, N)];
ionTemplates = [[DetailedFrag{Float32}() for _ in range(1, M)] for _ in range(1, N)];
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml");
scored_PSMs = [Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
unscored_PSMs = [[SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
spectral_scores = [Vector{SpectralScoresSimple{Float16}}(undef, 5000) for _ in range(1, N)];
precursor_weights = [zeros(Float32, n_precursors ) for _ in range(1, N)];
precs = [Counter(UInt32, UInt8,n_precursors ) for _ in range(1, N)];
chromatograms = [Vector{ChromObject}(undef, 5000) for _ in range(1, N)];
complex_scored_PSMs = [Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
complex_unscored_PSMs = [[ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
complex_spectral_scores = [Vector{SpectralScoresComplex{Float16}}(undef, 5000) for _ in range(1, N)];

###########
#File Names Parsing 
###########
file_names = first.(split.(basename.(MS_TABLE_PATHS), '.'))
split_file_names = split.(file_names, "_")



#If file names have inconsistnat number of delimiters, give up on parsing and use the entire file name
unique_split_file_name_lengths = unique(length.(split_file_names))
parsed_file_names = copy(split_file_names)
M = length(parsed_file_names)

if length(unique_split_file_name_lengths) == 1
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

const file_id_to_parsed_name = Dict(zip(1:M, [string(x) for x in parsed_file_names]))
#Parsed file names 
const parsed_fnames = sort(collect(values(file_id_to_parsed_name)))


