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
params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_012724_sulfur"
SPEC_LIB_DIR = "/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/HEIL_2023/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/PXD028735/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
#MS_TABLE_PATHS = MS_TABLE_PATHS[1:4]
EXPERIMENT_NAME = "TEST_y4b3_nOf5"
=#
#=
params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));
SPEC_LIB_DIR = ""pwd()*"\\..\\data\\nOf3_y4b3_012724_sulfur"#"/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_012724_sulfur"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = pwd()*"\\..\\data\\RAW"#"/Users/n.t.wamsley/TEST_DATA/PXD028735/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "TEST_y4b3_nOf5"
=#
#=
3-protome PC test March 4th 2024 for hupo
params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));
SPEC_LIB_DIR =  "C:\\Users\\n.t.wamsley\\PROJECTS\\HUPO_2023\\HUMAN_YEAST_ECOLI\\PIONEER\\LIB"
MS_DATA_DIR = "C:\\Users\\n.t.wamsley\\PROJECTS\\HUPO_2023\\HUMAN_YEAST_ECOLI\\PIONEER\\RAW"
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
out_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS")
if !isdir(out_folder)
    mkpath(out_folder)
end
out_folder = joinpath(MS_DATA_DIR, "Search", "RESULTS")
if !isdir(out_folder)
    mkpath(out_folder)
end
out_folder = joinpath(MS_DATA_DIR, "Search", "PARAMS")
if !isdir(out_folder)
    mkpath(out_folder)
end

nnls_params = Dict{String, Any}(k => v for (k, v) in params["nnls_params"]);

params_ = (
    expected_matches = Int64(params["expected_matches"]),
    frag_ppm_err = Float64(params["frag_ppm_err"]),
    frag_tol_quantile = Float32(params["frag_tol_quantile"]),
    frag_tol_presearch = Float64(params["frag_tol_presearch"]),
    intensity_filter_fraction = Float32(params["intensity_filter_fraction"]),
    isotope_err_bounds = Tuple([Int64(bound) for bound in params["isotope_err_bounds"]]),
    LsqFit_tol = Float64(params["LsqFit_tol"]),
    Lsq_max_iter = Int64(params["Lsq_max_iter"]),
    nnls_max_iter = Int64(nnls_params["nnls_max_iter"]),
    max_peaks = typeof(params["max_peaks"]) == Bool ? params["max_peaks"] : Int64(params["max_peaks"]),
    max_peak_width = Float64(params["max_peak_width"]),

    min_index_search_score = UInt8(params["min_index_search_score"]),
    min_index_search_score_presearch = UInt8(params["min_index_search_score_presearch"]),
    min_frag_count = Int64(params["min_frag_count"]),
    min_frag_count_presearch = Int64(params["min_frag_count_presearch"]),
    min_frag_count_integration = Int64(params["min_frag_count_integration"]),
    min_frag_count_index_search = Int64(params["min_frag_count_index_search"]),

    min_log2_matched_ratio_presearch = Float32(params["min_log2_matched_ratio_presearch"]),
  
    min_log2_matched_ratio_main_search = Float32(params["min_log2_matched_ratio_main_search"]),
    min_log2_matched_ratio_integration = Float32(params["min_log2_matched_ratio_integration"]),

    min_spectral_contrast = Float32(params["min_spectral_contrast"]),
    min_spectral_contrast_presearch = Float32(params["min_spectral_contrast_presearch"]),
    min_spectral_contrast_integration = Float32(params["min_spectral_contrast_integration"]),

    min_topn_of_m = Tuple([Int64(x) for x in params["min_topn_of_m"]]),
    min_topnOf3_integration = Int64(params["min_topnOf3_integration"]),
    min_weight = Float32(params["min_weight"]),

    most_intense = Bool(params["most_intense"]),
    n_quadrature_nodes = Int64(params["n_quadrature_nodes"]),
    n_frag_isotopes = Int64(params["n_frag_isotopes"]),
    nnls_tol = Float32(nnls_params["nnls_tol"]),
    quadrupole_isolation_width = Float64(params["quadrupole_isolation_width"]),
    regularize = Bool(nnls_params["regularize"]),
    rt_bounds = Tuple([Float64(rt) for rt in params["rt_bounds"]]),
    rt_tol = Float64(params["rt_tol"]),
    sample_rate = Float64(params["sample_rate"]),
    tail_distance = Float32(params["tail_distance"]),
    topN_presearch = Int64(params["topN_presearch"]),
    topN = Int64(params["topN"]),
    topN_index_search = Int64(params["topN_index_search"]),
    unmatched_penalty_factor = Float64(params["unmatched_penalty_factor"]),
    huber_δ = Float32(nnls_params["huber_delta"]),
    λ = Float32(nnls_params["lambda"]),
    γ = Float32(nnls_params["gamma"]),
);
#Search parameters
first_search_params = Dict(
    :collect_frag_errs => true,
    :expected_matches => params_[:expected_matches],
    :frag_ppm_err => params_[:frag_ppm_err],
    :frag_tol_presearch => params_[:frag_tol_presearch],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count_presearch],
    :min_frag_count_index_search => params_[:min_frag_count_index_search],
    :min_log2_matched_ratio => params_[:min_log2_matched_ratio_presearch],
    :min_index_search_score => params_[:min_index_search_score_presearch],
    :min_spectral_contrast => params_[:min_spectral_contrast_presearch],
    :min_topn_of_m => params_[:min_topn_of_m],
    :most_intense => params_[:most_intense],
    :nmf_tol => params_[:nnls_tol],
    :quadrupole_isolation_width => params_[:quadrupole_isolation_width],
    :regularize => params_[:regularize],
    :rt_bounds => params_[:rt_bounds],
    :rt_tol => Float64(1000),#params_[:rt_tol],
    :sample_rate => params_[:sample_rate],
    :topN => params_[:topN_presearch],
    :topN_index_search => params_[:topN_index_search],
    :λ => params_[:λ],
    :γ => params_[:γ],
    :huber_δ => params_[:huber_δ]
);
main_search_params = Dict(
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :isotope_err_bounds => params_[:isotope_err_bounds],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count],
    :min_frag_count_index_search => params_[:min_frag_count_index_search],
    :min_log2_matched_ratio_main_search => params_[:min_log2_matched_ratio_main_search],
    :min_index_search_score => params_[:min_index_search_score],
    :min_spectral_contrast => params_[:min_spectral_contrast],
    :min_topn_of_m => params_[:min_topn_of_m],
    :most_intense => params_[:most_intense],
    :nmf_tol => params_[:nnls_tol],
    :n_frag_isotopes => params_[:n_frag_isotopes],
    :quadrupole_isolation_width => params_[:quadrupole_isolation_width],
    :regularize => params_[:regularize],
    :rt_bounds => params_[:rt_bounds],
    :rt_tol => params_[:rt_tol],
    :topN => params_[:topN],
    :topN_index_search => params_[:topN_index_search],
    :λ => params_[:λ],
    :γ => params_[:γ],
    :huber_δ => params_[:huber_δ]
);
ms2_integration_params = Dict(
    :unmatched_penalty_factor => params_[:unmatched_penalty_factor],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :isotope_err_bounds => params_[:isotope_err_bounds],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :max_peak_width => params_[:max_peak_width],
    :min_frag_count => params_[:min_frag_count_integration],
    :min_log2_matched_ratio => params_[:min_log2_matched_ratio_integration],
    :min_spectral_contrast => params_[:min_spectral_contrast_integration],
    :min_topn_of_m => params_[:min_topn_of_m],
    :min_weight => params_[:min_weight],
    :most_intense => params_[:most_intense],
    :nmf_tol => params_[:nnls_tol],
    :n_frag_isotopes => params_[:n_frag_isotopes],
    :quadrupole_isolation_width => params_[:quadrupole_isolation_width],
    :regularize => params_[:regularize],
    :rt_bounds => params_[:rt_bounds],
    :rt_tol => params_[:rt_tol],
    :topN => params_[:topN],
    :topN_index_search => params_[:topN_index_search],
    :λ => params_[:λ],
    :γ => params_[:γ],
    :huber_δ => params_[:huber_δ]
);

if parse(Bool, lowercase(strip(ARGS["print_params"])))
    for (k,v) in zip(keys(params_), params_)
        pprint(k)
        print(" => ")
        pprintln(v)
    end
end

##########
#Load Dependencies 
##########
#Fragment Library Parsing

#[include(joinpath(pwd(), "src", jl_file)) for jl_file in ["IonType.jl","parseFasta.jl","PrecursorDatabase.jl"]];
#[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in ["IonType.jl"
#                                                                    #        ,"parseFasta.jl"
#                                                                    #"precursorDatabase.jl"
#                                                                    ]];

[include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in [
                                                                    "ArrayDict.jl",
                                                                    "Counter.jl",
                                                                    "Ion.jl",
                                                                    "LibraryIon.jl",
                                                                    "MatchIon.jl",
                                                                    "LibraryFragmentIndex.jl",
                                                                    "SparseArray.jl"]];
#[include(joinpath(pwd(), "src", "Routines","ParseProsit", jl_file)) for jl_file in ["parsePrositLib.jl"]];

#Generic files in src directory
#[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in ["IonType.jl","isotopes.jl"]]
#ML/Math Routines                                                                                    
[include(joinpath(pwd(), "src","ML", jl_file)) for jl_file in [
                                                                                            "percolatorSortOf.jl",
                                                                                            "kdeRTAlignment.jl",
                                                                                            "probitRegression.jl"]];
                                 
                                                                                            #Files needed for PSM scoring


#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in [
                                                                    "ExponentialGaussianHybrid.jl",
                                                                    "isotopes.jl",
                                                                    "globalConstants.jl",
                                                                    "isotopeSplines.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "SpectralDeconvolution.jl",
                                                                    
                                                                    "partitionThreadTasks.jl"]];

[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                        #"buildFragmentIndex.jl",
                                                                                    "matchpeaksLib.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    #"spectralDistanceMetrics.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl",
                                                                                    "integrateChroms.jl"]];
                                             
##########
#Load Spectral Library
#=
library_fragment_lookup_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"lib_frag_lookup", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index_7ppm_2hi", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors", file) != nothing, readdir(SPEC_LIB_DIR))][1]
=#
#library_fragment_lookup_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"lookup_table\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"index\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1]

#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
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
spec_load_time = @timed begin
    @time const f_index = load(f_index_path)["library_fragment_index"];
    prosit_lib["f_index"] = f_index;#["f_index"]
    @time const library_fragment_lookup_table = load(library_fragment_lookup_path)["library_fragment_lookup_table"]
    prosit_lib["f_det"] = library_fragment_lookup_table; #["f_det"];
    @time const precursors = load(precursors_path)["library_precursors"]
    prosit_lib["precursors"] = precursors;#["precursors"];
end
###########
#Load RAW File
#=
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = ["/Users/n.t.wamsley/TEST_DATA/PXD028735/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"]
MS_TABLE_PATHS = ["/Users/n.t.wamsley/TEST_DATA/mzXML/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"]
=#

###########
#Pre-Search
#Need to Estimate the following from a random sample of high-confidence targets
#1) Fragment Mass error/correction
#2) Fragment Mass tolerance
#3) iRT to RT conversion spline
###########



println("Begin processing: "*string(length(MS_TABLE_PATHS))*" files")
println("Starting Pre Search...")
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