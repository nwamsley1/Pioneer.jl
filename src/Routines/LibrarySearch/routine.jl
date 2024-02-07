
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
SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_012723_sulfur"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/PXD028735/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
#MS_TABLE_PATHS = MS_TABLE_PATHS[1:4]
EXPERIMENT_NAME = "TEST_y4b3_nOf5"
=#
#=
params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));
SPEC_LIB_DIR = pwd()*"\\..\\data\\nOf3_y4b3_012724_sulfur"#"/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_012724_sulfur"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = pwd()*"\\..\\data\\RAW"#"/Users/n.t.wamsley/TEST_DATA/PXD028735/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "TEST_y4b3_nOf5"
=#
println("ARGS ", ARGS)
MS_DATA_DIR = ARGS["data_dir"];
EXPERIMENT_NAME = ARGS["experiment_name"];
SPEC_LIB_DIR = ARGS["spec_lib_dir"];

#Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
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

    min_frag_count = Int64(params["min_frag_count"]),
    min_frag_count_presearch = Int64(params["min_frag_count_presearch"]),
    min_frag_count_integration = Int64(params["min_frag_count_integration"]),
    min_frag_count_index_search = Int64(params["min_frag_count_index_search"]),

    min_matched_ratio_presearch = Float32(params["min_matched_ratio_presearch"]),
  
    min_matched_ratio_main_search = Float32(params["min_matched_ratio_main_search"]),
    min_matched_ratio_integration = Float32(params["min_matched_ratio_integration"]),

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
    :min_matched_ratio => params_[:min_matched_ratio_presearch],
    :min_index_search_score => params_[:min_index_search_score],
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
    :min_matched_ratio_main_search => params_[:min_matched_ratio_main_search],
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
    :min_matched_ratio => params_[:min_matched_ratio_integration],
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

[include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in ["Ion.jl",
                                                                    
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
                                                                                            ]];
                                 
                                                                                            #Files needed for PSM scoring


#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in ["counter.jl",
                                                                    "exponentialGaussianHybrid.jl",
                                                                    "isotopes.jl",
                                                                    "globalConstants.jl",
                                                                    "isotopeSplines.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "spectralDeconvolution.jl"]];

[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch", jl_file)) for jl_file in [
                                                                                        #"buildFragmentIndex.jl",
                                                                                    "matchpeaksLib.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    #"spectralDistanceMetrics.jl",
                                                                                    "refinePSMs.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "logitRegression.jl"]];
                                             
##########
#Load Spectral Library
library_fragment_lookup_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"lookup_table\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"index\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
spec_load_time = @timed begin
    @time const f_index = load(f_index_path)["library_fragment_index"];
    prosit_lib["f_index"] = f_index#["f_index"]
    @time const library_fragment_lookup_table = load(library_fragment_lookup_path)["library_fragment_lookup_table"]
    prosit_lib["f_det"] = library_fragment_lookup_table #["f_det"];
    @time const precursors = load(precursors_path)["library_precursors"]
    prosit_lib["precursors"] = precursors#["precursors"];
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
N = 24
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
end
presearch_time = @timed begin
#init_frag_tol = 30.0 #Initial tolerance should probably be pre-determined for each different instrument and resolution. 
RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,MassErrorModel}()
lk = ReentrantLock()
#MS_TABLE_PATHS = MS_TABLE_PATHS[1:4]
#=
ms_file_idx = 1
MS_TABLE_PATH = MS_TABLE_PATHS[1]
=#
for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"
    @time RESULT =  firstSearch(
                                            MS_TABLE,
                                            prosit_lib["f_index"],
                                            prosit_lib["precursors"],
                                            prosit_lib["f_det"],
                                            x->x, #RT to iRT map'
                                            UInt32(ms_file_idx), #MS_FILE_IDX
                                            Laplace(zero(Float64), first_search_params[:frag_tol_presearch]),
                                            first_search_params,
                                            ionMatches,
                                            ionMisses,
                                            all_fmatches,
                                            IDtoCOL,
                                            ionTemplates,
                                            iso_splines,
                                            scored_PSMs,
                                            unscored_PSMs,
                                            spectral_scores,
                                            precursor_weights,
                                            precs,
                                            #scan_range = (100000, 110000)
                                            scan_range = (1, length(MS_TABLE[:masses]))
                                            );
    rtPSMs = vcat([first(result) for result in RESULT]...)
    all_matches = vcat([last(result) for result in RESULT]...)
    #@time begin

    @time refineFirstSearchPSMs!(rtPSMs, MS_TABLE, precursors)
    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end
    #Get Retention Times and Target/Decoy Status 
    ####################
    #Use best_psms to estimate 
    #1) RT to iRT curve and 
    #2) mass error (ppm) distribution 
    @time begin
    best_precursors = Set(rtPSMs[:,:precursor_idx]);
    best_matches = [match for match in all_matches if match.prec_id ∈ best_precursors];
    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];
    frag_ppm_intensities = [match.intensity for match in best_matches];
    end

    mass_err_model = ModelMassErrs(
        frag_ppm_intensities,
        frag_ppm_errs,
        40.0
    )
    
    #Model fragment errors with a mixture model of a uniform and laplace distribution 
    lock(lk) do 
        PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
        @time RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT_predicted], n = 50, bandwidth = 4.0);
        @time plotRTAlign(rtPSMs[:,:RT], rtPSMs[:,:iRT_predicted], RT_to_iRT_map, 
                    f_out = PLOT_PATH);
        RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map
        frag_err_dist_dict[ms_file_idx] = mass_err_model
    end
end

end
println("Finished presearch in ", presearch_time.time, " seconds")

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "rt_map_dict_020224.jld2"); RT_to_iRT_map_dict)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict_020224.jld2"); frag_err_dist_dict)
#MS_TABLE_PATHS = MS_TABLE_PATHS[1:2]
RT_to_iRT_map_dict = load("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\rt_map_dict_020224.jld2")["RT_to_iRT_map_dict"]
frag_err_dist_dict = load("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\frag_err_dist_dict_020224.jld2")["frag_err_dist_dict"]
###########
#Main PSM Search
###########
println("Begining Main Search...")

PSMs_Dict = Dictionary{String, DataFrame}()
RT_INDICES = Dictionary{String, retentionTimeIndex{Float32, Float32}}()

main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
#@profview for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        println("starting file $ms_file_idx")
        @time begin
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
        @time PSMs = vcat(mainLibrarySearch(
            MS_TABLE,
            prosit_lib["f_index"],
            prosit_lib["precursors"],
            #prosit_lib["f_det"],
            library_fragment_lookup_table,
            RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
            UInt32(ms_file_idx), #MS_FILE_IDX
            frag_err_dist_dict[ms_file_idx],
            16.1,
            main_search_params,
            ionMatches,
            ionMisses,
            all_fmatches,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            scored_PSMs,
            unscored_PSMs,
            spectral_scores,
            precursor_weights,
            precs,
            #scan_range = (100000, 200000),
            scan_range = (1, length(MS_TABLE[:masses]))
        #scan_range = (100000, 100010)
        )...);

        addMainSearchColumns!(PSMs, MS_TABLE, prosit_lib["precursors"]);
        getRTErrs!(PSMs);

        column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,
                        :iRT_error,:missed_cleavage,:Mox,:charge,:TIC,
                        :y_count,:err_norm,:spectrum_peak_count,:intercept]

        scoreMainSearchPSMs!(PSMs,
                                    column_names,
                                    n_train_rounds = 3,
                                    max_iter_per_round = 10,
                                    max_q_value = 0.01);
                                    
        getBestPSMs!(PSMs,
                        prosit_lib["precursors"],
                        max_q_value = 0.10);

        println("retained ", size(PSMs, 1), " psms")

        insert!(PSMs_Dict, 
            MS_TABLE_PATH, 
            PSMs[!,
                [:precursor_idx,:RT,:iRT_predicted,:prec_mz,:q_value,:prob]
                ]
        );
        end
end
println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_020224_M0.jld2"); PSMs_Dict)

rows_ = best_psms[!,:file_path] .== "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"
filter!(x -> x.q_value<=0.01, best_psms)
rows_ = best_psms[!,:file_path] .== "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"
precs_ = Set(best_psms[rows_,:precursor_idx])
Set(PSMs[!,:precursor_idx]) ∩ precs_
Set(PSMs_Dict.values[1][!,:precursor_idx]) ∩ precs_
x = 1
PSMs_Dict = load(joinpath("C:\\Users\\n.t.wamsley\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS", "PSMs_Dict_013024_M0.jld2"))["PSMs_Dict"]
@time begin
iRT_RT, RT_iRT = mapRTandiRT(PSMs_Dict)
precID_to_iRT = getPrecIDtoiRT(PSMs_Dict, RT_iRT)
RT_INDICES = makeRTIndices(PSMs_Dict,precID_to_iRT,iRT_RT)
end

BPSMS = Dict{Int64, DataFrame}()
PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];
scored_PSMs = [Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
unscored_PSMs = [[ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
spectral_scores = [Vector{SpectralScoresComplex{Float16}}(undef, 5000) for _ in range(1, N)];
end

quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    println("starting file $ms_file_idx")
    PSMS = vcat(integrateMS2_(MS_TABLE, 
                    prosit_lib["precursors"],
                    prosit_lib["f_det"],
                    RT_INDICES[MS_TABLE_PATH],
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    30.0f0,
                    ms2_integration_params,  
                    ionMatches,
                    ionMisses,
                    all_fmatches,
                    IDtoCOL,
                    ionTemplates,
                    iso_splines,
                    scored_PSMs,
                    unscored_PSMs,
                    spectral_scores,
                    precursor_weights,
                    scan_range = (1, length(MS_TABLE[:scanNumber])),
                    )...);

    filter!(x->x.weight>0.0, PSMS);
    addSecondSearchColumns!(PSMS, MS_TABLE, prosit_lib["precursors"]);
    getIsoRanks!(PSMS, MS_TABLE, 8.0036/2)
    addIntegrationFeatures!(PSMS)
    MS2_CHROMS = groupby(PSMS, [:precursor_idx,:iso_rank]);
    integratePrecursors(MS2_CHROMS, 
                                        n_quadrature_nodes = params_[:n_quadrature_nodes],
                                        intensity_filter_fraction = params_[:intensity_filter_fraction],
                                        α = 0.001f0,
                                        LsqFit_tol = params_[:LsqFit_tol],
                                        Lsq_max_iter = params_[:Lsq_max_iter],
                                        tail_distance = params_[:tail_distance]
                        )
    filter!(x -> x.best_scan, PSMS);
    PSMS[!,:ms_file_idx].=ms_file_idx
    BPSMS[ms_file_idx] = PSMS;
    GC.gc()
end

best_psms = vcat(values(BPSMS)...)

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_M0M1_q995to25_020224.jld2"); best_psms)
#best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_020124.jld2"))["best_psms"]
@time begin
getBestTrace!(best_psms)
addChromatogramFeatures!(PSMS, MS_TABLE, prosit_lib["precursors"],
                RT_iRT,
                precID_to_iRT);
#=
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_012524.jld2"); best_psms)
best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_012524.jld2"))["best_psms"];
=#
#best_psms = copy(MS2_CHROMS)
#getBestTrace!(best_psms)

#best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0M1_00_lambda0_topn2_total4_minimpute_siblingtest_5m_penalty50_011824.jld2"))["best_psms"]
features = [ 
    #:iRT_diff,
    :max_prob,
    :iso_rank,
    :max_city_fitted,
    :mean_city_fitted,
    #:max_scribe_fitted,
    :max_ions,
    :isotope_fraction,
    :assymetry,
    :fraction_censored,
    :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :Mox,
    :iRT_predicted,
    :iRT_error,
    :RT,
    #:RT_error,
    :iRT_diff,
    :longest_y,
    :y_count,
    :b_count,
    :base_width_min,
    :best_rank,
    :charge,
    :city_block,
    :city_block_fitted,
    :data_points,
    :entropy_score,
    :err_norm,
    :error,
    :hyperscore,
    #:intensity_explained,
    :ions_sum,
    #:log_sum_of_weights,
    :matched_ratio,
    :max_entropy,
    #:mean_log_probability,
    :max_spectral_contrast,
    :max_matched_ratio,
    :max_scribe_score,
    :missed_cleavage,
    #:ms1_ms2_diff,
    :peak_area,
    #:peak_area_ms1,
    :points_above_FWHM,
    :points_above_FWHM_01,
    :poisson,
    :prec_mz,
    :scribe,
    :scribe_corrected,
    :scribe_fitted,
    :sequence_length,
    :spectral_contrast,
    :spectral_contrast_corrected,
    #:spectrum_peaks,
    :topn,
    :total_ions,
    :weight,
    #:y_ladder,
    #:ρ,
    :log2_intensity_explained,
    :TIC,
    :adjusted_intensity_explained
    ]

best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
#best_psms = best_psms[best_psms[:,:file_path] .== "/Users/n.t.wamsley/TEST_DATA/mzXML/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow",:]
xgboost_time = @timed bst = rankPSMs2!(best_psms, 
                        features,
                        colsample_bytree = 1.0, 
                        min_child_weight = 5, 
                        gamma = 1, 
                        subsample = 0.5, 
                        n_folds = 2, 
                        num_round = 200, 
                        max_depth = 10, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        n_iters = 2);
best_psms = bst[3]
best_psms[!,:prob] = Float32.(best_psms[!,:prob])
getQvalues!(best_psms[!,:prob], best_psms[:,:decoy], best_psms[!,:q_value]);
value_counts(df, col) = combine(groupby(df, col), nrow)
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:file_path])
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers
) => :accession_numbers
);
end
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0M1_020224_q99min5max25.jld2"); best_psms)
end
histogram(best_psms[best_psms[!,:target], :prob], normalize = :probability, alpha = 0.5)
histogram!(best_psms[best_psms[!,:decoy], :prob], normalize = :probability, alpha = 0.5)

filter!(x->x.best_scan==true, best_psms );
#best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0M1_00_lambda0_topn2_total4_penalty50_clean5_012224.jld2"))["best_psms"]
#filter!(x->x.q_value<=0.01, best_psms );
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_020124_q975_5to40.jld2"); best_psms)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "ids_M0_010224.jld2"); IDs_PER_FILE)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "xgboost_M0_010224.jld2"); bst)



test = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_013024_ppm16.jld2"))