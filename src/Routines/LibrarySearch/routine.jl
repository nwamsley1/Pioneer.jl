
##########
#Import Libraries
##########
#Data Parsing/Printing
println("Importing Libraries...")

using ArgParse
using CSV, Arrow, Tables, DataFrames, JSON, JLD2, ProgressBars
using Plots, PrettyPrinting
#DataStructures 
using DataStructures, Dictionaries, Distributions, Combinatorics, StatsBase, LinearAlgebra, Random, LoopVectorization, SparseArrays
#Algorithms 
using Interpolations, XGBoost, SavitzkyGolay, NumericalIntegration, ExpectationMaximization, LsqFit, FastGaussQuadrature, GLM, StaticArrays
using Base.Order

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
SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123_sulfur"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/PXD028735/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATHS = MS_TABLE_PATHS[1:1]
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

    min_index_search_score = Float32(params["min_index_search_score"]),

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

    min_topnOf3 = Int64(params["min_topnOf3"]),
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
    :min_topnOf3 => params_[:min_topnOf3],
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
    :min_topnOf3 => params_[:min_topnOf3],
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
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :max_peak_width => params_[:max_peak_width],
    :min_frag_count => params_[:min_frag_count_integration],
    :min_matched_ratio => params_[:min_matched_ratio_integration],
    :min_spectral_contrast => params_[:min_spectral_contrast_integration],
    :min_topnOf3 => params_[:min_topnOf3_integration],
    :min_weight => params_[:min_weight],
    :most_intense => params_[:most_intense],
    :nmf_tol => params_[:nnls_tol],
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
integrate_ms1_params = Dict(
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
    :max_peak_width => params_[:max_peak_width],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count],
    :min_matched_ratio => params_[:min_matched_ratio],
    :min_spectral_contrast => params_[:min_spectral_contrast],
    :most_intense => params_[:most_intense],
    :nmf_tol => params_[:nnls_tol],
    :precursor_tolerance => params_[:prec_tolerance],
    :quadrupole_isolation_width => params_[:quadrupole_isolation_width],
    :regularize => params_[:regularize],
    :rt_bounds => params_[:rt_bounds],
    :rt_tol => params_[:rt_tol],
    :sample_rate => 1.0,
    :topN => params_[:topN],
    :λ => params_[:λ],
    :γ => params_[:γ]
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

[include(joinpath(pwd(), "src", jl_file)) for jl_file in ["IonType.jl","parseFasta.jl","PrecursorDatabase.jl"]];

[include(joinpath(pwd(), "src", "Routines","ParseProsit", jl_file)) for jl_file in ["buildPrositCSV.jl",
                                                                                    "parsePrositLib.jl"]];

#Generic files in src directory
[include(joinpath(pwd(), "src", jl_file)) for jl_file in ["precursor.jl","isotopes.jl"]];

#ML/Math Routines                                                                                    
[include(joinpath(pwd(), "src","ML", jl_file)) for jl_file in ["customsparse.jl","sparseNNLS.jl",
                                                                                            "percolatorSortOf.jl",
                                                                                            "kdeRTAlignment.jl",
                                                                                            "entropySimilarity.jl",
                                                                                            "EGH.jl"]];
                                                                                        
#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in ["counter.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "isotopeSplines.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch", jl_file)) for jl_file in [
                                                                                        "buildFragmentIndex.jl",
                                                                                    "matchpeaksLib.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "spectralDistanceMetrics.jl",
                                                                                    "refinePSMs.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                "integrateChroms.jl",
                                                                                "getCrossCorr.jl",
                                                                                    "queryFragmentIndex.jl",
                                                                                    "scratch_newchroms.jl"]];


                                                                                                                                
#Files needed for PSM scoring
[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","LibraryXTandem.jl","LibraryIntensity.jl"]]

[include(joinpath(pwd(), "src", "Routines", "PRM","IS-PRM",jl_file)) for jl_file in ["getScanPairs.jl"]]

##########
#Load Spectral Library
f_det_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_det\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
spec_load_time = @timed begin
    const f_index = load(f_index_path)["f_index"];
    prosit_lib["f_index"] = f_index#["f_index"]
    const f_det = load(f_det_path)["f_det"]
    prosit_lib["f_det"] = f_det#["f_det"];
    const precursors = load(precursors_path)["precursors"]
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
presearch_time = @timed begin
#init_frag_tol = 30.0 #Initial tolerance should probably be pre-determined for each different instrument and resolution. 
RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,Laplace{Float64}}()
lk = ReentrantLock()
#=
ms_file_idx = 1
MS_TABLE_PATH = MS_TABLE_PATHS[1]
=#
for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"

    sub_search_time = @timed RESULT =  firstSearch(
                                            MS_TABLE,
                                            prosit_lib["f_index"],
                                            prosit_lib["precursors"],
                                            prosit_lib["f_det"],
                                            x->x, #RT to iRT map'
                                            UInt32(ms_file_idx), #MS_FILE_IDX
                                            Laplace(zero(Float64), first_search_params[:frag_tol_presearch]),
                                            first_search_params,
                                            #scan_range = (100000, 110000)
                                            scan_range = (1, length(MS_TABLE[:masses]))
                                            );
    rtPSMs = vcat([first(result) for result in RESULT]...)
    all_matches = vcat([last(result) for result in RESULT]...)
    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end
    
    #Get Retention Times and Target/Decoy Status 
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_lib["precursors"][psm[:precursor_idx]]))) => :iRT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_lib["precursors"][psm[:precursor_idx]])) => :decoy);
    filter!(:entropy_score => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), rtPSMs);
    rtPSMs[!,:total_ions] = rtPSMs[!,:y_count] .+ rtPSMs[!,:b_count]
    #filter!(:entropy_sim => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), rtPSMs);
    rtPSMs[:,:target] =  rtPSMs[:,:decoy].==false


    FORM = FormulaTerm(
        (Term(:target),),
        (Term(:poisson),
         Term(:hyperscore),
         Term(:topn),
         Term(:scribe),
         Term(:spectral_contrast),
         Term(:total_ions))
    )

    model_fit = glm(FORM, rtPSMs, 
    Binomial(), 
    ProbitLink())
    Y′ = GLM.predict(model_fit, rtPSMs);
    println(size(rtPSMs))
    rtPSMs = rtPSMs[Y′.>0.75,:]
    println(size(rtPSMs))
    #Retain only targets
    rtPSMs = rtPSMs[rtPSMs[:,:decoy].==false,:];


    ####################
    #Use best_psms to estimate 
    #1) RT to iRT curve and 
    #2) mass error (ppm) distribution 
    best_precursors = Set(rtPSMs[:,:precursor_idx]);

    best_matches = [match for match in all_matches if match.prec_id ∈ best_precursors];
    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];


    #Model fragment errors with a mixture model of a uniform and laplace distribution 
    lock(lk) do 
        PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
        frag_err_dist = estimateErrorDistribution(frag_ppm_errs, Laplace{Float64}, 0.0, 3.0, first_search_params[:frag_tol_presearch],
                        f_out = PLOT_PATH );

        #Spline mapping RT to iRT
        RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50, bandwidth = 4.0);
        plotRTAlign(rtPSMs[:,:RT], rtPSMs[:,:iRT], RT_to_iRT_map, 
                    f_out = PLOT_PATH);

        RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map
        frag_err_dist_dict[ms_file_idx] = frag_err_dist
    end
end

end
println("Finished presearch in ", presearch_time.time, " seconds")
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "rt_map_dict_010324.jld2"); RT_to_iRT_map_dict)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict_010324.jld2"); frag_err_dist_dict)

RT_to_iRT_map_dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "rt_map_dict_010324.jld2"))["RT_to_iRT_map_dict"]
frag_err_dist_dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict_010324.jld2"))["frag_err_dist_dict"]
###########
#Main PSM Search
###########
println("Begining Main Search...")

PSMs_Dict = Dictionary{String, DataFrame}()
RT_INDICES = Dictionary{String, retentionTimeIndex{Float32, Float32}}()

Profile.Allocs.clear()
Profile.init()
Profile.Allocs.@profile [x for x in [x for x in range(1, 10000)]]
pprof(from_c = false)

results = Profile.Allocs.fetch()
last(sort(results.allocs, by=x->x.size))

retrieve = Profile.Allocs.retrieve()
main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
#@profview for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
        Profile.Allocs.@profile begin 
            @time PSMs = vcat(mainLibrarySearch(
            MS_TABLE,
            prosit_lib["f_index"],
            prosit_lib["precursors"],
            prosit_lib["f_det"],
            RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
            UInt32(ms_file_idx), #MS_FILE_IDX
            frag_err_dist_dict[ms_file_idx],
            16.1,
            main_search_params,
            scan_range = (100000, 200000),
            #scan_range = (1, length(MS_TABLE[:masses]))
        #scan_range = (100000, 100010)
        )...);
        end

        refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"], max_rt_error = 10.0);
        
        filter!(x->x.q_value<=0.30, PSMs);
        #bins = LinRange(0, 1, 100)
        #histogram(PSMs[(PSMs[!,:q_value].<=0.01),:entropy_score], normalize=:probability, bins = bins, alpha = 0.5)
        #histogram!(PSMs[(PSMs[!,:q_value].>0.1),:entropy_score], normalize=:probability, bins = bins, alpha = 0.5)
        #histogram(PSMs[(PSMs[!,:decoy].==false),:entropy_score], normalize=:probability, bins = bins, alpha = 0.5)
        #histogram!(PSMs[(PSMs[!,:decoy]),:entropy_score], normalize=:probability, bins = bins, alpha = 0.5)
        #############
        #Get best scan for every precursor
        sort!(PSMs,:RT); #Sorting before grouping is critical. 
        PSMs[!,:max_weight] .= zero(Float32)
        grouped_PSMs = groupby(PSMs, :precursor_idx);
        for i in range(1, length(grouped_PSMs))
            best_scan = argmax(grouped_PSMs[i].matched_ratio)
            grouped_PSMs[i].best_scan[best_scan] = true
        end
        filter!(x->x.best_scan, PSMs);

        transform!(PSMs, AsTable(:) => ByRow(psm -> 
        prosit_lib["precursors"][psm[:precursor_idx]].mz
        ) => :prec_mz
        );

        sort!(PSMs,:RT, rev = false);
        #Build RT index of precursors to integrate
        rt_index = buildRTIndex(PSMs);

        #RT_indices[MS_TABLE_PATH] = rt_index
        insert!(RT_INDICES, MS_TABLE_PATH, rt_index)#joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name);
        #jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name); PSMs);
        #Free up memory
        insert!(PSMs_Dict, MS_TABLE_PATH, PSMs);
end
end
println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_INDICES_isotopes_010224_isotopes_M0M1.jld2"); RT_INDICES)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_isotopes_010224_isotopes_M0M1.jld2"); PSMs_Dict)
#RT_INDICES = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_INDICES_isotopes_122023.jld2"));
#RT_INDICES = RT_INDICES["RT_INDICES"]

PSMs_Dict = nothing
GC.gc()
#PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_121423.jld2"));


BPSMS = Dict{Int64, DataFrame}()
PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];

quantitation_time = for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
    MS2_CHROMS = vcat(integrateMS2(MS_TABLE, 
                    prosit_lib["precursors"],
                    prosit_lib["f_det"],
                    RT_INDICES[MS_TABLE_PATH],
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    16.1,
                    ms2_integration_params, 
                    scan_range = (1, length(MS_TABLE[:scanNumber])),
                    #scan_range = (101357, 102357),
                    #scan_range = (101357, 110357)
                    )...);

    _refinePSMs!(MS2_CHROMS, MS_TABLE, prosit_lib["precursors"]);
    sort!(MS2_CHROMS,:RT); #Sorting before grouping is critical. 
    MS2_CHROMS_GROUPED = groupby(MS2_CHROMS, :precursor_idx);

    integratePrecursors(MS2_CHROMS_GROUPED, 
                                        n_quadrature_nodes = params_[:n_quadrature_nodes],
                                        intensity_filter_fraction = params_[:intensity_filter_fraction],
                                        α = 0.001f0,
                                        LsqFit_tol = params_[:LsqFit_tol],
                                        Lsq_max_iter = params_[:Lsq_max_iter],
                                        tail_distance = params_[:tail_distance])
    println("size(MS2_CHROMS) 1", size(MS2_CHROMS))
    filter!(x -> x.best_scan, MS2_CHROMS);
    println("size(MS2_CHROMS) 2", size(MS2_CHROMS))
    filter!(x->x.weight>0, MS2_CHROMS);
    println("size(MS2_CHROMS) 3", size(MS2_CHROMS))
    #for i in range(1, size(MS2_CHROMS)[1])
    #    if isinf(MS2_CHROMS[i,:mean_matched_ratio])
    #        MS2_CHROMS[i,:mean_matched_ratio] = Float16(6000)
    #    end
    #end

    transform!( MS2_CHROMS, AsTable(:) => ByRow(psm -> 
    prosit_lib["precursors"][psm[:precursor_idx]].mz
    ) => :prec_mz
    );
    MS2_CHROMS[:,:file_path].=MS_TABLE_PATH
    BPSMS[ms_file_idx] = MS2_CHROMS;
end

#=
for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH) 
        PSMs = BPSMS[ms_file_idx];
        isotopes = UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}()
        #lock(lk) do 
            isotopes = getIsotopes(PSMs[!,:sequence], 
                                    PSMs[!,:precursor_idx], 
                                    PSMs[!,:charge], QRoots(4), 4);
        #end
        prec_rt_table = sort(
                                collect(
                                        zip(
                                            PSMs[!,:tᵣ], 
                                            UInt32.(PSMs[!,:precursor_idx])
                                    )
                            ), by = x->first(x));
    
        #println("Integrating MS1 $ms_file_idx...")
        ms1_chroms = integrateMS1(MS_TABLE, 
                                        isotopes, 
                                        prec_rt_table, 
                                        UInt32(ms_file_idx), 
                                        Laplace(-1.0, 4.0),
                                        #frag_err_dist_dict[ms_file_idx], 
                                        integrate_ms1_params,
                                        scan_range = (1, length(MS_TABLE[:scanNumber])));

        #Integrate MS1 Chromatograms 
        gx, gw = gausslegendre(100);
        ms1_chrom_keys = keys(ms1_chroms);
        time_test = @timed transform!(PSMs, AsTable(:) => ByRow(psm -> integratePrecursor(ms1_chroms, 
                                            ms1_chrom_keys,
                                            gx,
                                            gw,
                                            UInt32(psm[:precursor_idx]), 
                                            [Float32(psm[:σ]),
                                            Float32(psm[:tᵣ]),
                                            Float32(psm[:τ]), 
                                            Float32(1e4)],
                                            isplot = false)) => [:peak_area_ms1,
                                            :GOF_ms1,
                                            :FWHM_ms1,
                                            :FWHM_01_ms1,
                                            :asymmetry_ms1,
                                            :points_above_FWHM_ms1,
                                            :points_above_FWHM_01_ms1,
                                            :σ_ms1,:tᵣ_ms1,:τ_ms1,:H_ms1]);

        pearson_corr!(PSMs, N = 500);
        PSMs[:,:ms1_ms2_diff] = abs.(PSMs[!,:tᵣ] .- PSMs[!,:tᵣ_ms1]);
        #########
        #Quality Filter
        PSMs[isnan.(coalesce.(PSMs[:,:GOF_ms1], 0.0)),:GOF_ms1].=missing;
        PSMs[:,"sequence_length"] = length.(replace.(PSMs[:,"sequence"], "M(ox)" => "M"));

        for i in range(1, size(PSMs)[1])
            if ismissing(PSMs[i,:ρ])
                continue
            end
            if isnan(PSMs[i,:ρ])
                PSMs[i,:ρ] = missing
            end
        end
end
=#

best_psms = vcat(values(BPSMS)...)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_010224_M0M1.jld2"); best_psms)
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_INDICES_110223.jld2"); RT_INDICES)
best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_010224_M0.jld2"))["best_psms"];
#RT_INDICES = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_INDICES_110123.jld2"));
#RT_INDICES = RT_INDICES["RT_INDICES"]
#MS2_CHROMS = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_110123.jld2"));
BPSMS = nothing
GC.gc()

filter!(x -> x.best_scan, best_psms);
filter!(x->x.weight>0, best_psms );

    #for i in range(1, size(MS2_CHROMS)[1])
    #    if isinf(MS2_CHROMS[i,:mean_matched_ratio])
    #        MS2_CHROMS[i,:mean_matched_ratio] = Float16(6000)
    #    end
    #end

    #transform!( MS2_CHROMS, AsTable(:) => ByRow(psm -> 
    #prosit_lib["precursors"][psm[:precursor_idx]].mz
    #) => :prec_mz
    #);

#=
for i in range(1, size(best_psms)[1])
    if isinf(best_psms[i,:mean_log_probability])
        println(i)
        best_psms[i,:mean_log_probability] = Float16(-6000)
    end
end
for name in names(best_psms)
    if eltype(best_psms[!,name]) == String
        continue
    else
        if any(isinf.(best_psms[!,name]))
            println(name)
            println(argmax(isinf.(best_psms[!,name])))
        end
    end
end
=#
for i in range(1, size(best_psms)[1])
    if ismissing(best_psms[i,:mean_log_probability])
        continue
    end
    if isinf(best_psms[i,:mean_log_probability])
        best_psms[i,:mean_log_probability] = Float16(6000)*sign(best_psms[i,:mean_log_probability] )
    end
end
for i in range(1, size(best_psms)[1])
    if isinf(best_psms[i,:GOF])
        println(i)
        best_psms[i,:GOF] = Float16(-6000)
    end
end
for i in range(1, size(best_psms)[1])
    if isinf(best_psms[i,:mean_matched_ratio])
        println(i)
        best_psms[i,:mean_matched_ratio] = Float16(6000)
    end
end
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers
) => :accession_numbers
);
for name in names(best_psms)
    if eltype(best_psms[:,name]) == String
        continue
    end
    #println("name $name")
    if any(skipmissing(isinf.(best_psms[:,name])))
        println("name $name is inf")
    end
end
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "all_psms_110223_mzXML.jld2"); best_psms)
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_unscored_nOf5_110123.jld2"); best_psms)
#best_psms = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_unscored_nOf5_110123.jld2"));
for i in range(1, size(best_psms)[1])
    if isinf(best_psms[i,:mean_log_probability])
        println(i)
        best_psms[i,:mean_log_probability] = Float16(-6000)
    end
end
#=
open(joinpath(MS_DATA_DIR, "Search", "RESULTS", "meta.txt"), "w") do f
    write(f, "Main Search Time $main_search_time \n")
end

open(joinpath(MS_DATA_DIR, "Search", "PARAMS","params.json"), "w") do f
    write(f, JSON.json(JSON.parse(read(ARGS["params_json"], String))))
end

filter!(x -> x.data_points > 2, best_psms)
=#
#best_psms = MS2_CHROMS
features = [ 
    :max_ions,
    :isotope_fraction,
    :assymetry,
    :fraction_censored,
    :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :Mox,
    :RT,
    :RT_error,
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
    :mean_log_probability,
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
best_psms[:,"sequence_length"] = length.(replace.(best_psms[:,"sequence"], "M(ox)" => "M"));
best_psms[!,:isotope_fraction] = best_psms[!,:isotope_count]./(best_psms[!,:b_count] .+ best_psms[!,:y_count])
#best_psms = best_psms[best_psms[:,:file_path] .== "/Users/n.t.wamsley/TEST_DATA/mzXML/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow",:]
xgboost_time = @timed bst = rankPSMs!(best_psms, 
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
                        train_fraction = 4.0/9.0,
                        n_iters = 2);

getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
value_counts(df, col) = combine(groupby(df, col), nrow)
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:file_path])
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0_010224.jld2"); best_psms)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "ids_M0_010224.jld2"); IDs_PER_FILE)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "xgboost_M0_010224.jld2"); bst)



IDs_M0M1 = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "ids_M0M1_010224.jld2"))["IDs_PER_FILE"]
bst_M0M1 = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "xgboost_M0M1_010224.jld2"))["bst"]

best_psms_passing = best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:decoy].==false),:]
pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)

pioneer_passing_fdr = Set("_".*best_psms[!,:stripped_sequence].*"_.".*string.(best_psms[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)


#@save "/Users/n.t.wamsley/TEST_DATA/best_psms_unfiltered_16ppm_huber10000_y4b3bestand1_cosineCorrected_102323.jld2" best_psms
       


#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_nOf5_110123.jld2"); best_psms)
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers
) => :accession_numbers
);
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_y4b3_nOf5_122123_isotopes.jld2"); best_psms)

using CSV, DataFrames, StatsBase, Plots, StatsPlots, Measures, JLD2, FASTX, CodecZlib, Loess, KernelDensity, Distributions, SavitzkyGolay,Interpolations

function parseFasta(fasta_path::String, parse_identifier::Function = x -> split(x,"|")[2])

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    #I/O for Fasta
    reader = getReader(fasta_path)

    #In memory representation of FastaFile
    #fasta = Vector{FastaEntry}()
    fasta = Vector{Tuple{String, String}}()
    @time begin
        for record in reader
                push!(fasta, 
                        (parse_identifier(FASTA.identifier(record)),
                         split(split(split(FASTA.description(record), ' ')[1], '|')[end], '_')[end],
                                #FASTA.sequence(record),
                                #false
                        )
                )
        end
    end

    return fasta
end

PIONEER = best_psms[(best_psms[!,:q_value].<=0.1).&(best_psms[!,:decoy].==false).&(best_psms[!,:data_points].>2),:];
PIONEER = PIONEER[(iszero.(PIONEER[:,:peak_area]).==false),:];
MEANS = combine(prec_group -> (mean_log2 = mean(log2.(prec_group[:,:peak_area])),
                               median_log2 = median(log2.(prec_group[:,:peak_area])),
                               CV = std(prec_group[:,:peak_area])/mean(prec_group[:,:peak_area])
                            ),
                groupby(PIONEER[:,[:precursor_idx,:peak_area]], 
                        [:precursor_idx])
    );
#Exclude precursors without enough data to calculate a global CV
MEANS = MEANS[isnan.(MEANS[:,:CV]).==false,:];

#Use only precursors with CV in the bottom 25%
MEANS = MEANS[MEANS[:,:CV].<quantile(MEANS[:,:CV], 0.25),:];

#Dictionary mapping precursor id's to named tuple with global log2_mean, log2_median, and CV
mean_dict = Dict(zip(MEANS[:,:precursor_idx], [NamedTuple(x) for x in eachrow(MEANS[:,[:mean_log2,:median_log2,:CV]])]));

PIONEER[:,:CV] = allowmissing(zeros(Float64, size(PIONEER)[1]));
#Apply precursor group CVs to the rows of the data
for i in range(1,size(PIONEER)[1])
    #if PIONEER[i,:file_path]=="/Users/n.t.wamsley/TEST_DATA/mzXML/LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.arrow"
    if haskey(mean_dict, PIONEER[i,:precursor_idx])
        PIONEER[i,:CV] = mean_dict[PIONEER[i,:precursor_idx]][:CV]
    else
        PIONEER[i,:CV] = missing
    end
end

PIONEER[:,:peak_area_norm] = PIONEER[:,:peak_area];
PIONEER[:,:peak_area_log2] = log2.(PIONEER[:,:peak_area]);

PIONEER_BEST = PIONEER[(PIONEER[:,:q_value].<=0.01).&( #1% FDR
                        PIONEER[:,:decoy].==false).&(
                        coalesce.(PIONEER[:,:CV], Inf).<median(skipmissing(PIONEER[:,:CV]))), #CV less than median
    :]

INTENSITIES = combine(file -> (median = median(file[:,:peak_area_log2]),mean = mean(file[:,:peak_area_log2])), groupby(PIONEER_BEST[:,[:file_path,:peak_area_log2]], :file_path));

#INTENSITIES = combine(file -> (median = median(file[:,:peak_area_log2]),mean = mean(file[:,:peak_area_log2])), groupby(PIONEER_BEST[:,[:file_path,:peak_area_log2]], :file_path));

INTENSITIES[:,:correction] = (INTENSITIES[:,:median] .- median(INTENSITIES[:,:median]));

NORM_CORRECTIONS = Dict(zip(INTENSITIES[:,:file_path],INTENSITIES[:,:correction]));

transform!(PIONEER, AsTable(:) => ByRow(precursor -> 2^(precursor[:peak_area_log2] .- NORM_CORRECTIONS[precursor[:file_path]])) => [:peak_area_norm]);
PIONEER[:,:peak_area_log2_norm] .= log2.(PIONEER[:,:peak_area_norm]);
#try
p = plot(legend=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n Original Response",
topmargin=5mm)

@df PIONEER[(PIONEER[:,:q_value].<=0.01).&(
             PIONEER[:,:decoy].==false),:] boxplot(p, (:file_path), :peak_area_log2, ylim = (10, 25), ylabel = "Log2 (Response)", xlabel = "Run ID", show = true)
#savefig(p, joinpath(MS_DATA_DIR, "Search", "RESULTS", "RAW_INTENSITIES.pdf"))

p = plot(legend=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n Normalized Response",
topmargin=5mm)
@df PIONEER[(PIONEER[:,:q_value].<=0.01).&(
            PIONEER[:,:decoy].==false),:] boxplot(p, (:file_path), :peak_area_log2_norm, ylim = (10, 25), ylabel = "Log2 (Response)", xlabel = "Run ID", show = true) 

#savefig(p, joinpath(MS_DATA_DIR, "Search", "RESULTS", "NORM_INTENSITIES.pdf"))            
##catch
#    println("couldn't printn normalization plots")
#end 
ACC_TO_SPEC = Dict(vcat([parseFasta("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/UP000000625_83333_Escherichia_coli.fasta.gz"),
            parseFasta("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/UP000002311_559292_Saccharomyces_cerevisiae.fasta.gz"),
            parseFasta("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/UP000005640_9606_human.fasta.gz")]...));

PIONEER[:,:species_set] = [Set([ACC_TO_SPEC[id] for id in ids]) for ids in split.(PIONEER[:,:accession_numbers],';')];
            PIONEER = PIONEER[length.(PIONEER[:,:species_set]).==1,:];
            PIONEER[:,:species] = first.(PIONEER[:,:species_set]);
value_counts(df, col) = combine(groupby(df, col), nrow)
IDs_PER_FILE = value_counts(PIONEER[(PIONEER[:,:q_value].<=0.01) .& (PIONEER[:,:decoy].==false),:], [:file_path])

CSV.write(joinpath(MS_DATA_DIR, "Search", "RESULTS", "id_per_file.csv"), IDs_PER_FILE)


function getCondition(row)
    file_name = split(split(row[:file_path], '/')[end], '.')[1]
    file_name_split = split(file_name, '_')
    condition, biological, technical = file_name_split[end - 3], file_name_split[end - 1], file_name_split[end]
    return condition, biological, technical
end

transform!(PIONEER, AsTable(:) => ByRow(precursor -> getCondition(precursor)) => [:condition, :biological, :technical]);

PIONEER_GROUPED = groupby(PIONEER[(PIONEER[:,:decoy].==false).&(PIONEER[:,:Mox].==0).&(PIONEER[:,:data_points].>2),:], [:species,:accession_numbers,:sequence,:charge,:condition]);
PIONEER_LOG2MEAN = combine(prec_group -> (log2_mean = log2(mean(prec_group[:,:peak_area_norm])), 
                                        CV = 100*std(prec_group[:,:peak_area_norm])/mean(prec_group[:,:peak_area_norm]),
                                        non_missing = length(prec_group[:,:peak_area_norm]),
                                        min_q_value = minimum(prec_group[:,:q_value])), 
    PIONEER_GROUPED );


#=
PIONEER_LOG2MEAN = combine(prec_group -> (log2_mean = log2(mean(trim(prec_group[:,:peak_area_norm], prop = 0.2))), 
                                        CV = 100*std(trim(prec_group[:,:peak_area_norm], prop = 0.2))/mean(trim(prec_group[:,:peak_area_norm], prop = 0.2)),
                                        non_missing = length(prec_group[:,:peak_area_norm]),
                                        min_q_value = minimum(prec_group[:,:q_value])), 
    PIONEER_GROUPED );
=#


filter!(:min_q_value => x -> x<=0.01, PIONEER_LOG2MEAN);

function getLog2Diff(prec_group)
    #println(size(prec_group))
    if size(prec_group)[1].!=2
        return (log2_diff = zeros(Float64, 1),
                 log2_mean = zeros(Float64, 1),
                 log2_a = zeros(Float64, 1),
                 log2_b = zeros(Float64, 1),
                 CV_a = zeros(Float64, 1),
                 CV_b = zeros(Float64, 1),
                 nonMissing_a = zeros(Int64, 1),
                 nonMissing_b = zeros(Int64, 1))
    end
    group_A = occursin.("A", prec_group[:,:condition])
    group_B = group_A.==false
    log2_diff = prec_group[:,:log2_mean][group_B] - prec_group[:,:log2_mean][group_A]
    out = (
     log2_diff = log2_diff, 
     log2_mean = (prec_group[:,:log2_mean][group_A] + prec_group[:,:log2_mean][group_B])/2,
     log2_a = prec_group[:,:log2_mean][group_A],
     log2_b = prec_group[:,:log2_mean][group_B],
     CV_a = prec_group[:,:CV][group_A],
     CV_b = prec_group[:,:CV][group_B],
     nonMissing_a = prec_group[:,:non_missing][group_A],
     nonMissing_b = prec_group[:,:non_missing][group_B]
        )
    return out
end

PIONEER_COMBINED = combine(prec_group -> getLog2Diff(prec_group), groupby(PIONEER_LOG2MEAN, [:species, :accession_numbers, :sequence, :charge]));
PIONEER_COMBINED = PIONEER_COMBINED[(PIONEER_COMBINED[:,:nonMissing_a].>3).&(PIONEER_COMBINED[:,:nonMissing_b].>3),:];
PIONEER_COMBINED = PIONEER_COMBINED[(PIONEER_COMBINED[:,:CV_a].<30).&(PIONEER_COMBINED[:,:CV_b].<30),:];

open(joinpath(MS_DATA_DIR, "Search", "RESULTS", "passing_benchmark.txt"), "w") do f
    #size = "size(PIONEER_COMBINED) "*string(size(PIONEER_COMBINED))
    write(f, "size(PIONEER_COMBINED) "*string(size(PIONEER_COMBINED)))
end

gdf = groupby(PIONEER_COMBINED , [:species])
nt = NamedTuple.(keys(gdf))
p = plot(legend=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n 6of9, CV<30%, 1% FDR \n 27498 Precursors",
topmargin=5mm)
i = 1
SPECIES_TO_LOG2FC = Dict("HUMAN" => 0.0,
                         "YEAST" => -1.0,
                         "ECOLI" => 2.0)
for (k,v) in pairs(gdf)
    density!(p, gdf[k][:,:log2_diff], label=nothing, color = i, bins = LinRange(-3, 3, 100), show = true, normalize = :pdf, alpha = 0.5)
    vline!([SPECIES_TO_LOG2FC[nt[i][:species]]], color = i, label = label="$(nt[i][:species])")
    i += 1
end

savefig(p, joinpath(MS_DATA_DIR, "Search", "RESULTS", "HIST.pdf"))

gdf = groupby(PIONEER_COMBINED , [:species])
nt = NamedTuple.(keys(gdf))
p = plot(legend=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n 6of9, CV<30%, 1% FDR \n 27498 Precursors",
topmargin=5mm)
i = 1

for (k,v) in pairs(gdf)
    violin!(p, gdf[k][:,:log2_diff], label="$(nt[i])", ylim = (-3, 4), span = [10.0])
    i += 1
end

savefig(p, joinpath(MS_DATA_DIR, "Search", "RESULTS", "VIOLIN.pdf"))

gdf = groupby(PIONEER_COMBINED , [:species])
nt = NamedTuple.(keys(gdf))
p = plot(legend=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n 6of9, CV<30%, 1% FDR \n 27498 Precursors",
topmargin=5mm)
i = 1
for (k,v) in pairs(gdf)
    plot!(p, gdf[k][:,:log2_mean], gdf[k][:,:log2_diff], color = i, label=nothing, xlim = (10, 30), ylim = (-3, 4), alpha = 0.1, seriestype=:scatter)
    hline!([SPECIES_TO_LOG2FC[nt[i][:species]]], color = i, label = label="$(nt[i][:species])")
  
    i += 1
end

savefig(p, joinpath(MS_DATA_DIR, "Search", "RESULTS", "SCATTER.pdf"))



####################
#Get Differential Expression/Benchmark

#=
features = [ :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :Mox,
    :RT,
    :RT_error,
    :b_ladder,
    :base_width_min,
    :best_rank,
    :charge,
    :city_block,
    :data_points,
    :entropy_sim,
    :err_norm_log2,
    :error,
    :hyperscore,
    :intensity_explained,
    :ions_sum,
    :log_sum_of_weights,
    :matched_ratio,
    :mean_log_entropy,
    :mean_log_probability,
    :mean_log_spectral_contrast,
    :mean_matched_ratio,
    :mean_scribe_score,
    :missed_cleavage,
    :ms1_ms2_diff,
    :peak_area,
    :peak_area_ms1,
    :points_above_FWHM,
    :points_above_FWHM_01,
    :poisson,
    :prec_mz,
    :scribe_score,
    :sequence_length,
    :spectral_contrast,
    :spectrum_peaks,
    :topn,
    :total_ions,
    :weight,
    :y_ladder,
    :ρ,

    :log2_base_peak_intensity,
    :TIC,
    :adjusted_intensity_explained
    ]


time_test = @timed filter!(x -> x.topn>1, best_psms);
time_test = @timed filter!(x -> x.best_rank<3, best_psms);
time_test = @timed filter!(x -> x.points_above_FWHM>0, best_psms);
time_test = @timed filter!(x -> x.data_points>2, best_psms);
#time_test = @timed filter!(x -> x.data_points>2, best_psms);
#time_test = @timed filter!(x -> x.points_above_FWHM>0, best_psms);

xgboost_time = @timed bst = rankPSMs!(best_psms, 
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

getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
best_psms_passing = best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false).&(ismissing.(best_psms[:,:peak_area]).==false),:]
=#