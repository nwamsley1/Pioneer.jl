
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
using Interpolations, XGBoost, SavitzkyGolay, NumericalIntegration, ExpectationMaximization, LsqFit, FastGaussQuadrature, GLM
##########
#Parse Arguments 
##########
#Example Usage 
#julia --threads 24 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/Projects/PROSIT/TEST_DATA/MOUSE_TEST /Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423 -s true 
#julia --threads 9 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/TEST_DATA/mzXML/ /Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HumanYeastEcoli/5ppm_15irt/ -s true
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
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATHS = MS_TABLE_PATHS[1:1]
=#
println("ARGS ", ARGS)
MS_DATA_DIR = ARGS["data_dir"];
SPEC_LIB_DIR = ARGS["spec_lib_dir"];

#Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
out_folder = joinpath(MS_DATA_DIR, "Search")
if !isdir(out_folder)
    mkpath(out_folder)
end
out_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS")
if !isdir(out_folder)
    mkpath(out_folder)
end

nnls_params = Dict{String, Any}(k => v for (k, v) in params["nnls_params"]);

params_ = (
    b_min_ind = Int64(params["b_min_ind"]),
    y_min_ind = Int64(params["y_min_ind"]),
    expected_matches = Int64(params["expected_matches"]),
    frag_ppm_err = Float64(params["frag_ppm_err"]),
    frag_tol_quantile = Float32(params["frag_tol_quantile"]),
    init_frag_tol = Float64(params["frag_tol_presearch"]),
    nnls_max_iter = Int64(nnls_params["nnls_max_iter"]),
    max_peaks = typeof(params["max_peaks"]) == Bool ? params["max_peaks"] : Int64(params["max_peaks"]),
    max_peak_width = Float64(params["max_peak_width"]),
    min_frag_count = Int64(params["min_frag_count"]),
    min_frag_count_presearch = Int64(params["min_frag_count_presearch"]),
    min_frag_count_index_search = Int64(params["min_frag_count_index_search"]),
    min_matched_ratio = Float32(params["min_matched_ratio"]),
    min_matched_ratio_presearch = Float32(params["min_matched_ratio_presearch"]),
    min_matched_ratio_index_search = Float32(params["min_matched_ratio_index_search"]),
    min_spectral_contrast = Float32(params["min_spectral_contrast"]),
    min_spectral_contrast_presearch = Float32(params["min_spectral_contrast_presearch"]),
    nnls_tol = Float32(nnls_params["nnls_tol"]),
    prec_tolerance = Float64(params["prec_tolerance"]),
    quadrupole_isolation_width = Float64(params["quadrupole_isolation_width"]),
    regularize = Bool(nnls_params["regularize"]),
    rt_bounds = Tuple([Float64(rt) for rt in params["rt_bounds"]]),
    rt_tol = Float64(params["rt_tol"]),
    sample_rate = Float64(params["sample_rate"]),
    topN_presearch = Int64(params["topN_presearch"]),
    topN = Int64(params["topN"]),
    topN_index_search = Int64(params["topN_index_search"]),
    huber_δ = Float32(nnls_params["huber_delta"]),
    λ = Float32(nnls_params["lambda"]),
    γ = Float32(nnls_params["gamma"]),
);
#Search parameters
first_search_params = Dict(
    :b_min_ind => params_[:b_min_ind],
    :y_min_ind => params_[:y_min_ind],
    :collect_frag_errs => true,
    :expected_matches => params_[:expected_matches],
    :frag_ppm_err => params_[:frag_ppm_err],
    :fragment_tolerance => params_[:init_frag_tol],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count_presearch],
    :min_frag_count_index_search => params_[:min_frag_count_index_search],
    :min_matched_ratio => params_[:min_matched_ratio_presearch],
    :min_matched_ratio_index_search => params_[:min_matched_ratio_index_search],
    :min_spectral_contrast => params_[:min_spectral_contrast_presearch],
    :nmf_tol => params_[:nnls_tol],
    :precursor_tolerance => params_[:prec_tolerance],
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
    :b_min_ind => params_[:b_min_ind],
    :y_min_ind => params_[:y_min_ind],
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count],
    :min_frag_count_index_search => params_[:min_frag_count_index_search],
    :min_matched_ratio => params_[:min_matched_ratio],
    :min_matched_ratio_index_search => params_[:min_matched_ratio_index_search],
    :min_spectral_contrast => params_[:min_spectral_contrast],
    :nmf_tol => params_[:nnls_tol],
    :precursor_tolerance => params_[:prec_tolerance],
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
#=
integrate_ms2_params = Dict(
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
    :max_peak_width => params_[:max_peak_width],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count],
    :min_matched_ratio => params_[:min_matched_ratio],
    :min_spectral_contrast => params_[:min_spectral_contrast],
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
=#
integrate_ms1_params = Dict(
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
    :max_peak_width => params_[:max_peak_width],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count],
    :min_matched_ratio => params_[:min_matched_ratio],
    :min_spectral_contrast => params_[:min_spectral_contrast],
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
[include(joinpath(pwd(), "src","ML", jl_file)) for jl_file in ["sparseNNLS.jl",
                                                                                            "percolatorSortOf.jl",
                                                                                            "kdeRTAlignment.jl",
                                                                                            "entropySimilarity.jl",
                                                                                            "EGH.jl"]];
                                                                                        
#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in ["counter.jl",
                                                                    "massErrorEstimation.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch", jl_file)) for jl_file in ["buildFragmentIndex.jl",
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
[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","LibraryXTandem.jl"]]

[include(joinpath(pwd(), "src", "Routines", "PRM","IS-PRM",jl_file)) for jl_file in ["getScanPairs.jl"]]

##########
#Load Spectral Library
#Need to find a way to speed this up.  
#SPEC_LIB_DIR = "/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HumanYeastEcoli/5ppm_15irt/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/"
prosit_lib_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
println("Loading spectral libraries into main memory...")

#SPEC_LIB_DIR = "/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HumanYeastEcoli/5ppm_15irt/"
spec_load_time = @timed begin
    #const frag_list = load_object(frags_list_path)
    #const precursors_list = load_object(precursors_list_path)
    #const frag_index = load_object(frag_index_path)
    prosit_lib = load(prosit_lib_path);
    #prosit_lib = load("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEXED_HumanYeastEcoli_NCE33_corrected_092823.jld2")
    
    #@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1
    #@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/prosit_mouse_33NCEcorrected_start1_5ppm_15irt.jld2" prosit_mouse_33NCEcorrected_start1_5ppm_15irt 
end
#frag_index = prosit_mouse_33NCEcorrected_start1_5ppm_15irt
#precursors_list = precursors_mouse_detailed_33NCEcorrected_start1
#frag_list = frags_mouse_detailed_33NCEcorrected_start1

###########
#Load RAW File
#=
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = ["/Users/n.t.wamsley/TEST_DATA/PXD028735/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"]
MS_TABLE_PATHS = ["/Users/n.t.wamsley/TEST_DATA/mzXML/LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"]
=#
#"/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA_duplicate.arrow",
#"/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA_duplicate_2.arrow",
#"/Users/n.t.wamsley/RIS_temp/MOUSE_DIA/ThermoRawFileToParquetConverter-main/parquet_out/MA5171_MOC1_DMSO_R01_PZ_DIA_duplicate_3.arrow"]

##########
#Set Search parameters
#=
first_search_params = Dict(
    :collect_frag_errs => true,
    :expected_matches => 1000000,
    :frag_ppm_err => 0.0,
    :fragment_tolerance => 30.0,
    :max_iter => 1000,
    :max_peaks => false,
    :min_frag_count => 7,
    :min_matched_ratio => Float32(0.6),
    :min_spectral_contrast => Float32(0.95),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.5,
    :regularize => false,
    :rt_bounds => (0.0, 200.0),
    :rt_tol => 200.0,
    :sample_rate => 0.01,
    :topN => 5,
    :λ => zero(Float32),
    :γ => zero(Float32)
)

main_search_params = Dict(
    :expected_matches => 1000000,
    :frag_tol_quantile => 0.975,
    :max_iter => 1000,
    :max_peaks => false,
    :min_frag_count => 4,
    :min_matched_ratio => Float32(0.45),
    :min_spectral_contrast => Float32(0.5),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.5,
    :regularize => false,
    :rt_bounds =>(-20.0, 200.0),
    :rt_tol => 20.0,
    :topN => 1000,
    :λ => zero(Float32),
    :γ => zero(Float32)
)

integrate_ms2_params = Dict(
    :expected_matches => 1000000,
    :frag_tol_quantile => 0.95,
    :max_iter => 1000,
    :max_peak_width => 2.0,
    :max_peaks => false,
    :min_frag_count => 4,
    :min_matched_ratio => Float32(0.45),
    :min_spectral_contrast => Float32(0.5),
    :nmf_tol => Float32(100),
    :precursor_tolerance => 5.0,
    :quadrupole_isolation_width => 8.0,
    :regularize => false,
    :rt_bounds => (0.0, 200.0),
    :rt_tol => 20.0,
    :sample_rate => 1.0,
    :topN => 1000,
    :λ => zero(Float32),
    :γ => zero(Float32)
)
integrate_ms1_params = Dict(
        :expected_matches => 1000000,
        :frag_tol_quantile => 0.95,
        :max_iter => 1000,
        :max_peak_width => 2.0,
        :max_peaks => false,
        :min_frag_count => 4,
        :min_matched_ratio => Float32(0.45),
        :min_spectral_contrast => Float32(0.5),
        :nmf_tol => Float32(100),
        :precursor_tolerance => 5.0,
        :quadrupole_isolation_width => 8.5,
        :regularize => false,
        :rt_tol => 20.0,
        :rt_bounds => (0.0, 200.0),
        :sample_rate => 1.0,
        :topN => 100,
        :λ => zero(Float32),
        :γ => zero(Float32)
)
MS_DATA_DIR = "/Users/n.t.wamsley/Projects/PROSIT/TEST_DATA/MOUSE_TEST"

MS_TABLE_PATHS = ["/Users/n.t.wamsley/RIS_temp/ThermoRawFileToParquetConverter-main/parquet_out//LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.arrow"]
MS_DATA_DIR = "/Users/n.t.wamsley/RIS_temp/ThermoRawFileToParquetConverter-main/parquet_out/"
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
init_frag_tol = 30.0 #Initial tolerance should probably be pre-determined for each different instrument and resolution. 
RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,Laplace{Float64}}()
lk = ReentrantLock()
Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"
    rtPSMs, all_matches =  firstSearch(
                                        MS_TABLE,
                                        prosit_lib["f_index"],
                                        prosit_lib["f_det"],
                                        x->x, #RT to iRT map'
                                        UInt32(ms_file_idx), #MS_FILE_IDX
                                        Laplace(zero(Float64), first_search_params[:fragment_tolerance]),
                                        first_search_params,
                                        scan_range = (1, length(MS_TABLE[:masses]))
                                        );

    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end
    
    #Get Retention Times and Target/Decoy Status 
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(prosit_lib["precursors"][psm[:precursor_idx]]))) => :iRT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> isDecoy(prosit_lib["precursors"][psm[:precursor_idx]])) => :decoy);
    filter!(:entropy_sim => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), rtPSMs);
    rtPSMs[:,:target] =  rtPSMs[:,:decoy].==false

    model_fit = glm(@formula(target ~ poisson + hyperscore + topn +
    scribe_score + topn + spectral_contrast + intensity_explained + total_ions), rtPSMs, 
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
        frag_err_dist = estimateErrorDistribution(frag_ppm_errs, Laplace{Float64}, 0.0, 3.0, 30.0,
                        f_out = PLOT_PATH );

        #Spline mapping RT to iRT
        RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50, bandwidth = 1.0);
        plotRTAlign(rtPSMs[:,:RT], rtPSMs[:,:iRT], RT_to_iRT_map, 
                    f_out = PLOT_PATH);

        RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map
        frag_err_dist_dict[ms_file_idx] = frag_err_dist
    end
end
end
println("Finished presearch in ", presearch_time.time, " seconds")

###########
#Main PSM Search
###########
println("Begining Main Search...")
BPSMS = Dict{Int64, DataFrame}()

main_search_time = @timed Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
#@profview for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)    
        sub_search_time = @timed PSMs = mainLibrarySearch(
                                                MS_TABLE,
                                                prosit_lib["f_index"],
                                                prosit_lib["f_det"],
                                                RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
                                                UInt32(ms_file_idx), #MS_FILE_IDX
                                                frag_err_dist_dict[ms_file_idx],
                                                main_search_params,
                                                #scan_range = (201389, 204389),
                                                
                                                scan_range = (1, length(MS_TABLE[:masses]))
                                            );

        @save "/Users/n.t.wamsley/TEST_DATA/PSMs_TEST_101623.jld2" PSMs;
        @load "/Users/n.t.wamsley/TEST_DATA/PSMs_TEST_101623.jld2" PSMs;
        #PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/PSMs_TEST_101423.csv"))
        #PSMs[isnan.(PSMs[:,:entropy_sim]),:entropy_sim] .= 0.0;
        #filter!(x -> x.topn>1, PSMs);
        #filter!(x -> x.best_rank<3, PSMs);
        filter!(x -> x.weight>10.0, PSMs);
        refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
        
        
        sort!(PSMs,:RT); #Sorting before grouping is critical. 
        test_chroms = groupby(PSMs, :precursor_idx);


        #Integrate MS2 and filter
        time_test = @timed integratePrecursors(test_chroms)
        time_test = @timed filter!(x -> x.best_scan, PSMs);
        time_test = @timed filter!(x -> x.FWHM.<100, PSMs);
        time_test = @timed filter!(x -> x.FWHM_01.<100, PSMs);
        time_test = @timed filter!(x -> x.data_points.>2, PSMs);

        #Get Precursor M/Z's 
        transform!(PSMs, AsTable(:) => ByRow(psm -> 
        prosit_lib["precursors"][psm[:precursor_idx]].mz
        ) => :prec_mz
        );

        #Correct RT prediction errors 
        for i in range(1, size(PSMs)[1])
            if ismissing(PSMs[i,:tᵣ])
                PSMs[i, :RT_error] = abs.(PSMs[i,:RT_pred] - PSMs[i,:RT])
                PSMs[i,:tᵣ] = PSMs[i,:RT]
            else
                PSMs[i,:RT_error] = abs.(PSMs[i,:tᵣ] .- PSMs[i,:RT_pred]);
            end
        end

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
    
        println("Integrating MS1 $ms_file_idx...")
        ms1_chroms = integrateMS1(MS_TABLE, 
                                        isotopes, 
                                        prec_rt_table, 
                                        UInt32(ms_file_idx), 
                                        frag_err_dist_dict[ms_file_idx], 
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
        #PSMs[coalesce.(PSMs[:,:GOF_ms1], -Inf).<=0.4,:peak_area_ms1].=missing
        #PSMs[coalesce.(PSMs[:,:GOF_ms1], -Inf).<=0.4,:ρ].=missing
        #best_psms[coalesce.(best_psms[:,:peak_area_ms1], 0.0).<=100.0,:peak_area_ms1].=missing
        #best_psms[coalesce.(best_psms[:,:peak_area_ms1], 0.0).<=100.0,:ρ].=missing
        PSMs[:,"sequence_length"] = length.(replace.(PSMs[:,"sequence"], "M(ox)" => "M"));

        sub_search_time = @timed lock(lk) do 
            println("ms_file_idx: $ms_file_idx has ", size(PSMs))
            BPSMS[ms_file_idx] = PSMs
        end
        println("Assigne best_psms_time for $ms_file_idx ", sub_search_time.time)
        #println("TEST length(prec_counts) ", length(prec_counts))
end
println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")

best_psms = vcat(values(BPSMS)...)

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

for i in range(1, size(best_psms)[1])
    if ismissing(best_psms[i,:ρ])
        continue
    end
    if isnan(best_psms[i,:ρ])
best_psms[i,:ρ] = missing
    end
end
#time_test = @timed filter!(x -> x.topn>1, best_psms);
#time_test = @timed filter!(x -> x.best_rank<3, best_psms);
       
time_test = @timed filter!(x -> x.data_points>2, best_psms);
time_test = @timed filter!(x -> x.points_above_FWHM>0, best_psms);

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

points_per_peak = value_counts(best_psms, :data_points)
bar(points_per_peak[!,:data_points], 
points_per_peak[!,:nrow]./sum(points_per_peak[!,:nrow]), alpha = 0.5)
points_per_peak = value_counts(best_psms_passing, :data_points)
bar!(points_per_peak[!,:data_points], 
points_per_peak[!,:nrow]./sum(points_per_peak[!,:nrow]), alpha = 0.5)


points_per_peak = value_counts(best_psms, :points_above_FWHM)
bar(points_per_peak[!,:points_above_FWHM], 
points_per_peak[!,:nrow]./sum(points_per_peak[!,:nrow]), alpha = 0.5)
points_per_peak = value_counts(best_psms_passing, :points_above_FWHM)
bar!(points_per_peak[!,:points_above_FWHM], 
points_per_peak[!,:nrow]./sum(points_per_peak[!,:nrow]), alpha = 0.5)

##########
#Regroup PSMs by file id 
#CSV.write("/Users/n.t.wamsley/TEST_DATA/PSMs_FIRSTTEST.csv", PSMs)
#PSMs = groupby(PSMs,:ms_file_idx);

########################
#Save PSMS 
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)

########################
#Quantification 
println("Begining Precursor Quantification ...")
integration_time = @timed Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)    
    
    best_psms = nothing
    #lock(lk) do 
    best_psms = BPSMS[ms_file_idx]     
    #end
    #Get Best PSM per precursor
    #best_psms = combine(sdf -> getBestPSM(sdf), groupby(PSMs[ms_file_idx][PSMs[ms_file_idx][:,:q_value].<=0.25,:], [:sequence,:charge]));
    println("SIZE best_psms $ms_file_idx", size(BPSMS[ms_file_idx]))
    transform!(best_psms, AsTable(:) => ByRow(psm -> 
                prosit_lib["precursors"][psm[:precursor_idx]].mz
                ) => :prec_mz
                );
    #=
    #Need to sort RTs 
    sort!(best_psms,:RT, rev = false);
    #Build RT index of precursors to integrate
    rt_index = buildRTIndex(best_psms);

    println("Integrating MS2 $ms_file_idx...")
    sub_search_time = @timed ms2_chroms = integrateMS2(MS_TABLE, 
                                    prosit_lib["f_det"],
                                    rt_index,
                                    UInt32(ms_file_idx), 
                                    frag_err_dist_dict[ms_file_idx],
                                    integrate_ms2_params, 
                                    scan_range = (0, length(MS_TABLE[:scanNumber]))
                                    #can_range = (101357, 110357)
                                    );
    println("ms2_chroms time for $ms_file_idx ", sub_search_time.time)
    println("ms2_chroms gc_stats for $ms_file_idx ", sub_search_time.gcstats)
    println("ms2_chroms gctime for $ms_file_idx ", sub_search_time.gctime)
    #Integrate MS2 Chromatograms 
    #Need to speed up this part. 
    sub_search_time = @timed transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms2_chroms, 
                                                UInt32(psm[:precursor_idx]), 
                                                (0.1f0, #Fraction peak height α
                                                0.15f0, #Distance from vertical line containing the peak maximum to the leading edge at α fraction of peak height
                                                0.15f0, #Distance from vertical line containing the peak maximum to the trailing edge at α fraction of peak height
                                                Float32(psm[:RT]), psm[:weight]), isplot = false)) => [:peak_area,:GOF,:FWHM,:FWHM_01,:asymmetry,:points_above_FWHM,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H]);
    println("integrate ms2 time for $ms_file_idx ", sub_search_time.time)
    println("integrate ms2 gctats for $ms_file_idx ", sub_search_time.gcstats)
    println("integrate ms2 gctime for $ms_file_idx ", sub_search_time.gctime)
    #Remove Peaks with 0 MS2 intensity or fewer than 6 points accross the peak. 
    #best_psms = best_psms[(ismissing.(best_psms[:,:peak_area]).==false).&(best_psms[:,:points_above_FWHM].>=1).&(best_psms[:,:points_above_FWHM_01].>=5),:];
    =#
    #best_psms = best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:];
    for i in range(1, size(best_psms)[1])
        if ismissing(best_psms[i,:tᵣ])
            best_psms[i, :RT_error] = abs.(best_psms[i,:RT_pred] - best_psms[i,:RT])
            best_psms[i,:tᵣ] = best_psms[i,:RT]
        else
            best_psms[i,:RT_error] = abs.(best_psms[i,:tᵣ] .- best_psms[i,:RT_pred]);
        end
    end

    #Get Predicted Isotope Distributions 
    #For some reason this requires a threadlock. Need to investiage further. 
    isotopes = UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}()
    lock(lk) do 
        isotopes = getIsotopes(best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:sequence], 
                                        best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:precursor_idx], 
                                        best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:charge], QRoots(4), 4);
    end

    prec_rt_table = sort(collect(zip(best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:tᵣ], 
    UInt32.(best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:precursor_idx]))), by = x->first(x));

    println("Integrating MS1 $ms_file_idx...")
    ms1_chroms = integrateMS1(MS_TABLE, 
                                    isotopes, 
                                    prec_rt_table, 
                                    UInt32(ms_file_idx), 
                                    frag_err_dist_dict[ms_file_idx], 
                                    integrate_ms1_params,
                                    scan_range = (1, length(MS_TABLE[:scanNumber])));

    #Integrate MS1 Chromatograms 
    gx, gw = gausslegendre(100)
    ms1_chrom_keys = keys(ms1_chroms)
    time_test = @timed transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms1_chroms, 
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
    
    #Get Correlation between MS1 and MS2 train_classes
    pearson_corr!(best_psms, N = 500)
    best_psms[:,:ms1_ms2_diff] = abs.(best_psms[!,:tᵣ] .- best_psms[!,:tᵣ_ms1]);
    #########
    #Quality Filter
    best_psms[isnan.(coalesce.(best_psms[:,:GOF_ms1], 0.0)),:GOF_ms1].=missing
    best_psms[coalesce.(best_psms[:,:GOF_ms1], -Inf).<=0.4,:peak_area_ms1].=missing
    best_psms[coalesce.(best_psms[:,:GOF_ms1], -Inf).<=0.4,:ρ].=missing
    #best_psms[coalesce.(best_psms[:,:peak_area_ms1], 0.0).<=100.0,:peak_area_ms1].=missing
    #best_psms[coalesce.(best_psms[:,:peak_area_ms1], 0.0).<=100.0,:ρ].=missing
    best_psms[:,"sequence_length"] = length.(replace.(best_psms[:,"sequence"], "M(ox)" => "M"));
    println("SIZE best_psms $ms_file_idx before filter ", size(best_psms))
    #Peak area is non-zero and non-missing
    #best_psms = best_psms[(best_psms[:,:peak_area].>0.0).&(ismissing.(best_psms[:,:peak_area]).==false),:]; 
    #FWHM must exceed minimum threshold
    #best_psms = best_psms[best_psms[:,:FWHM].>(6/60),:];
    #Minimum Goodness of Fit 
    #best_psms = best_psms[best_psms[:,:GOF].>0.6,:]; 
    println("SIZE best_psms $ms_file_idx after filter ", size(best_psms))
    #best_psms_old = best_psms

    lock(lk) do 
        BPSMS[ms_file_idx] = best_psms;
    end
    best_psms = nothing
end
println("Finished peak integration in ", integration_time, " seconds")
println("Train target-decoy model at precursor level...")
best_psms = vcat(values(BPSMS)...);

#CSV.write("/Users/n.t.wamsley/TEST_DATA/best_psms_100123_mod.csv", best_psms)

#best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_100123.csv"));
transform!(best_psms, AsTable(:) => ByRow(psm -> prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers) => :accession_numbers)
transform!(best_psms, AsTable(:) => ByRow(psm -> MS_TABLE_PATHS[psm[:ms_file_idx]]) => :file_path)



CSV.write("/Users/n.t.wamsley/TEST_DATA/TEST_ALL_unscored.csv", best_psms)

for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:peak_area] = best_psms[i,:peak_area]/base_peak_intensity
end


for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:H] = best_psms[i,:H]/base_peak_intensity
end

for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:weight] = best_psms[i,:weight]/base_peak_intensity
end

features = [:hyperscore,:total_ions,:intensity_explained,
            :poisson,:spectral_contrast,:entropy_sim,
            :RT_error,:scribe_score,:RT,:charge,
            :city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn,:err_norm,:error]
append!(features, [:peak_area,:GOF,:asymmetry,:points_above_FWHM_01,:points_above_FWHM,:H,:ρ,:FWHM, :FWHM_01, :y_ladder, :b_ladder, :n_obs,:peak_area_ms1,:prec_mz,:sequence_length,:spectrum_peaks,
:weight_sum,:hyperscore_sum,:entropy_sum,:scribe_sum,:ions_sum,:data_points,:ratio_sum,:base_width]);

#Train Model 
best_psms[:,:q_value] .= 0.0
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
println("unique ", length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx])))
CSV.write("/Users/n.t.wamsley/TEST_DATA/TEST_ALL_SCORED.csv", best_psms)
value_counts(df, col) = combine(groupby(df, col), nrow)
println(value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:ms_file_idx]))
#best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_100523_5ppm_15.csv"))
#CSV.write("/Users/n.t.wamsley/TEST_DATA/best_psms_100423_demult_adjusted.csv", PIONEER)
#=
for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:peak_area] = best_psms[i,:peak_area]/base_peak_intensity
end


for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:H] = best_psms[i,:H]/base_peak_intensity
end

for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:H] = best_psms[i,:weight]/base_peak_intensity
end
CSV.write("/Users/n.t.wamsley/TEST_DATA/best_psms_100523_5ppm_15.csv", best_psms)

CSV.write("/Users/n.t.wamsley/TEST_DATA/best_psms_100423_H_adjusted.csv", PIONEER)
=#
#=
#Model Features 
best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_FIRSTTEST.csv"));
features = [:hyperscore,:total_ions,:intensity_explained,
            :poisson,:spectral_contrast,:entropy_sim,
            :RT_error,:scribe_score,:RT,:charge,
            :city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn,:err_norm,:error]
append!(features, [:peak_area,:GOF,:asymmetry,:points_above_FWHM_01,:H,:ρ]);
#Train Model 
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 1.0, 
                        min_child_weight = 10, 
                        gamma = 10, 
                        subsample = 0.5, 
                        n_folds = 2, 
                        num_round = 200, 
                        max_depth = 10, 
                        eta = 0.0375, 
                        max_train_size = size(best_psms)[1]);

getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:]
=#


#=
for i in ProgressBar(size(best_psms)[1])
    prec_id = best_psms[i,:precursor_idx]
    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
    best_psms[i,:peak_area] = best_psms[i,:peak_area]/base_peak_intensity
end

#for i in ProgressBar(size(best_psms)[1])
#    prec_id = best_psms[i,:precursor_idx]
3    base_peak_intensity = prosit_lib["precursors"][prec_id].base_peak_intensity
#    best_psms[i,:peak_area_ms1] = best_psms[i,:peak_area_ms1]/base_peak_intensity
#end

test = unstack(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false), :],[:precursor_idx,:sequence,:charge],:ms_file_idx,:peak_area)
test = dropmissing(test)
transform!(test, AsTable(["1","2","3"]) => ByRow(x->mean(log10.(collect(x)))) => :mean)

plot(test[:,:mean], log10.(test[:,"1"]).-test[:,:mean], seriestype = :scatter)

plot(test[:,:mean], log10.(test[:,"2"]).-test[:,:mean], seriestype = :scatter, alpha = 0.01)
hline!([0.0])
hline!([mean(log10.(test[:,"2"]).-test[:,:mean])])

test[:,"2"] = test[:,"2"].*(medians[1]/medians[2])
test[:,"3"] = test[:,"3"].*(medians[1]/medians[3])

transform!(test, AsTable(["1","2","3"]) => ByRow(x->100*std(x)/mean(x)) => :CV)

test_ms1 = unstack(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false), :],[:precursor_idx,:sequence,:charge],:ms_file_idx,:peak_area_ms1)
test_ms1 = dropmissing(test_ms1)

transform!(test, AsTable(["1","2","3"]) => ByRow(x->100*std(x)/mean(x)) => :CV)

medians = [median(test[:,"1"]),median(test[:,"2"]),median(test[:,"3"])]


plot(log10.(test[:,"1"]), log10.(test[:,"2"]), seriestype=:scatter, alpha = 0.1)

histogram2d(log10.(test[:,"1"]), log10.(test[:,"2"]))
plot!(LinRange(0, 8, 100), LinRange(0, 8, 100))
plot!(LinRange(0, 8, 100), collect(LinRange(0, 8, 100)).+0.5, color = :black)
plot!(LinRange(0, 8, 100), collect(LinRange(0, 8, 100)).-0.5, color = :black)

transform!(test, AsTable(:) => ByRow(x -> std([x[c] for c in [:1,:2,:3]])/mean([x[c] for c in [:1,:2,:3]])) => :CV)

transform!(test, AsTable(["1","2","3"]) => ByRow(x->100*std(x)/mean(x)) => :CV)

transform(test, ["1","2","3"] => std)
best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_FIRSTTEST.csv"))

unique(best_psms[(best_psms[:,:q_value].<=0.05) .& (best_psms[:,:decoy].==false),[:precursor_idx]])

unique(best_psms[(best_psms[:,:q_value].<=0.05) .& (best_psms[:,:decoy].==false) .& (best_psms[:,:ms_file_idx].==1),[:precursor_idx]])
unique(best_psms[(best_psms[:,:q_value].<=0.05) .& (best_psms[:,:decoy].==false) .& (best_psms[:,:ms_file_idx].==1),[:sequence, :charge]])

unique(best_psms[(best_psms[:,:q_value].<=0.05) .& (best_psms[:,:decoy].==false) .& (best_psms[:,:ms_file_idx].==2),[:precursor_idx]])

unique(best_psms[(best_psms[:,:q_value].<=0.05) .& (best_psms[:,:decoy].==false) .& (best_psms[:,:ms_file_idx].==3),[:precursor_idx]])

#=
bins = LinRange(0, 1, 100)
histogram(best_psms[:,:prob][best_psms[:,:decoy]], alpha = 0.5, bins = bins, normalize = :probability)
histogram!(best_psms[:,:prob][best_psms[:,:decoy].==false], alpha = 0.5, bins = bins, normalize = :probability)
value_counts(df, col) = combine(groupby(df, col), nrow)
bins = LinRange(0, 50, 100)
histogram(best_psms[:,:FWHM][(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01)]*60, bins = bins)



bins = LinRange(0, 1, 100)
histogram(best_psms[:,:entropy_sim][(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01)], alpha = 0.5, bins = bins, normalize = :probability)
histogram!(best_psms[:,:entropy_sim][(best_psms[:,:decoy])], alpha = 0.5, bins = bins, normalize = :probability)

bins = LinRange(0, 15, 100)
histogram(best_psms[:,:scribe_score][(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01)], alpha = 0.5, bins = bins, normalize = :probability)
histogram!(best_psms[:,:scribe_score][(best_psms[:,:decoy])], alpha = 0.5, bins = bins, normalize = :probability)

bins = LinRange(0.5, 1.0, 100)
histogram(best_psms[:,:spectral_contrast][(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01)], alpha = 0.5, bins = bins, normalize = :probability)
histogram!(best_psms[:,:spectral_contrast][(best_psms[:,:decoy])], alpha = 0.5, bins = bins, normalize = :probability)

bins = LinRange(0.0, 30, 100)
histogram(best_psms[:,:RT_error][(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01)], alpha = 0.5, bins = bins, normalize = :probability)
histogram!(best_psms[:,:RT_error][(best_psms[:,:decoy])], alpha = 0.5, bins = bins, normalize = :probability)
=#
=#