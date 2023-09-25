
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
using Interpolations, XGBoost, SavitzkyGolay, NumericalIntegration, ExpectationMaximization
##########
#Parse Arguments 
##########
#Example Usage 
#julia --threads 24 ./src/Routines/LibrarySearch/routine.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/Projects/PROSIT/TEST_DATA/MOUSE_TEST /Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423 -s true 
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
    expected_matches = Int64(params["expected_matches"]),
    frag_ppm_err = Float64(params["frag_ppm_err"]),
    frag_tol_quantile = Float32(params["frag_tol_quantile"]),
    init_frag_tol = Float64(params["frag_tol_presearch"]),
    nnls_max_iter = Int64(nnls_params["nnls_max_iter"]),
    max_peaks = typeof(params["max_peaks"]) == Bool ? params["max_peaks"] : Int64(params["max_peaks"]),
    max_peak_width = Float64(params["max_peak_width"]),
    min_frag_count = Int64(params["min_frag_count"]),
    min_frag_count_presearch = Int64(params["min_frag_count_presearch"]),
    min_matched_ratio = Float32(params["min_matched_ratio"]),
    min_matched_ratio_presearch = Float32(params["min_matched_ratio_presearch"]),
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
    λ = Float32(nnls_params["lambda"]),
    γ = Float32(nnls_params["gamma"]),
);
#Search parameters
first_search_params = Dict(
    :collect_frag_errs => true,
    :expected_matches => params_[:expected_matches],
    :frag_ppm_err => params_[:frag_ppm_err],
    :fragment_tolerance => params_[:init_frag_tol],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count_presearch],
    :min_matched_ratio => params_[:min_matched_ratio_presearch],
    :min_spectral_contrast => params_[:min_spectral_contrast_presearch],
    :nmf_tol => params_[:nnls_tol],
    :precursor_tolerance => params_[:prec_tolerance],
    :quadrupole_isolation_width => params_[:quadrupole_isolation_width],
    :regularize => params_[:regularize],
    :rt_bounds => params_[:rt_bounds],
    :rt_tol => params_[:rt_tol],
    :sample_rate => params_[:sample_rate],
    :topN => params_[:topN_presearch],
    :λ => params_[:λ],
    :γ => params_[:γ]
);
main_search_params = Dict(
    :expected_matches => params_[:expected_matches],
    :frag_tol_quantile => params_[:frag_tol_quantile],
    :max_iter => params_[:nnls_max_iter],
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
    :topN => params_[:topN],
    :λ => params_[:λ],
    :γ => params_[:γ]
);
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
                                                                                            "entropySimilarity.jl"]];
                                                                                        
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
                                                                                    ]];


                                                                                                                                 
#Files needed for PSM scoring
[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","LibraryXTandem.jl"]]

[include(joinpath(pwd(), "src", "Routines", "PRM","IS-PRM",jl_file)) for jl_file in ["getScanPairs.jl"]]

##########
#Load Spectral Library
#Need to find a way to speed this up.  
frags_list_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"frags.*\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_list_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors.*\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
frag_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"prosit.*\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];

println("Loading spectral libraries into main memory...")
spec_load_time = @timed begin
    #const frag_list = load_object(frags_list_path)
    #const precursors_list = load_object(precursors_list_path)
    #const frag_index = load_object(frag_index_path)
    @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/frags_mouse_detailed_33NCEcorrected_start1.jld2" frags_mouse_detailed_33NCEcorrected_start1
    @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1
    @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/prosit_mouse_33NCEcorrected_start1_5ppm_15irt.jld2" prosit_mouse_33NCEcorrected_start1_5ppm_15irt 
end
frag_index = prosit_mouse_33NCEcorrected_start1_5ppm_15irt
precursors_list = precursors_mouse_detailed_33NCEcorrected_start1
frag_list = frags_mouse_detailed_33NCEcorrected_start1
println("Loaded spectral libraries in ", spec_load_time.time, " seconds")

###########
#Load RAW File
#MS_TABLE_PATHS = ["/Users/n.t.wamsley/Projects/PROSIT/TEST_DATA/MOUSE_TEST/MA5171_MOC1_DMSO_R01_PZ_DIA.arrow"],
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
    :frag_tol_quantile => 0.975,
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
    :rt_bounds => (0.0, 200.0),
    :rt_tol => 20.0,
    :sample_rate => 1.0,
    :topN => 1000,
    :λ => zero(Float32),
    :γ => zero(Float32)
)
integrate_ms1_params = Dict(
        :expected_matches => 1000000,
        :frag_tol_quantile => 0.975,
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
Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS[1:1]))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    rtPSMs, all_matches =  firstSearch(
                                        MS_TABLE,
                                        frag_index,  
                                        frag_list, 
                                        x->x, #RT to iRT map'
                                        UInt32(ms_file_idx), #MS_FILE_IDX
                                        Laplace(zero(Float64), 10.0),
                                        first_search_params,
                                        scan_range = (0, length(MS_TABLE[:masses]))
                                        );


    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end

    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors_list[psm[:precursor_idx]]))) => :iRT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT);
    transform!(rtPSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors_list[psm[:precursor_idx]])) => :decoy);
    rtPSMs = rtPSMs[rtPSMs[:,:decoy].==false,:];

    best_precursors = Set(rtPSMs[:,:precursor_idx]);
    best_matches = [match for match in all_matches if match.prec_id ∈ best_precursors];
    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];
    #Model fragment errors with a mixture model of a uniform and laplace distribution 
    PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
    frag_err_dist = estimateErrorDistribution(frag_ppm_errs, Laplace{Float64}, 0.0, 3.0, 30.0,
                    f_out = PLOT_PATH );

    RT_to_iRT_map = KDEmapping(rtPSMs[:,:RT], rtPSMs[:,:iRT], n = 50, bandwidth = 2.0);
    plotRTAlign(rtPSMs[:,:RT], rtPSMs[:,:iRT], RT_to_iRT_map, 
                f_out = PLOT_PATH);

    lock(lk) do 
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
PSMs_dict = Dict{Int64, DataFrame}()
main_search_time = @timed Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS[1:1]))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)    
        PSMs = mainLibrarySearch(
                                                MS_TABLE,
                                                frag_index,  
                                                frag_list, 
                                                RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
                                                UInt32(1), #MS_FILE_IDX
                                                frag_err_dist_dict[ms_file_idx],
                                                main_search_params,
                                                #scan_range = (101357, 111357),
                                                scan_range = (0, length(MS_TABLE[:masses]))
                                            );
        PSMs = PSMs[PSMs[:,:weight].>1000.0,:];
        refinePSMs!(PSMs, MS_TABLE, precursors_list);
        lock(lk) do 
            PSMs_dict[ms_file_idx] = PSMs
        end
        #println("TEST length(prec_counts) ", length(prec_counts))
end
println("Finished main search in ", main_search_time.time, "seconds")
############
println("Getting best PSMs...")
PSMs = vcat(values(PSMs_dict)...);
###########
#Clean features
############
PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf;
PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= maximum(PSMs[(PSMs[:,:matched_ratio]).!=Inf,:matched_ratio]);
replace!(PSMs[:,:city_block], -Inf => minimum(PSMs[PSMs[:,:city_block].!=-Inf,:city_block]));
replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]));
#PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/PSMs_080423.csv"))
transform!(PSMs, AsTable(:) => ByRow(psm -> length(collect(eachmatch(r"ox", psm[:sequence])))) => [:Mox]);

############
#Target-Decoy discrimination
############

features = [:hyperscore,:total_ions,:intensity_explained,:error,
            :poisson,:spectral_contrast,:entropy_sim,
            :RT_error,:scribe_score,:y_ladder,:RT,:n_obs,:charge,
            :city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn]


println("XGBoost for 10% FDR PSMs...")
main_search_xgboost_time = @timed rankPSMs!(PSMs, features, 
                colsample_bytree = 1.0, 
                min_child_weight = 10, 
                gamma = 10, 
                #subsample = 0.25, 
                subsample = 0.5,
                n_folds = 2,
                num_round = 200, 
                eta = 0.0375, 
                max_depth = 5,
                max_train_size = size(PSMs)[1]);

getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
println("Target PSMs at 1% FDR: ", sum((PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false)))

##########
#Regroup PSMs by file id 
PSMs = groupby(PSMs,:ms_file_idx);

#PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:precursor_idx]
#########
#save psms
#########
#CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv", PSMs)
###########
#Integrate 
#best_psms_old = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/best_psms_080423.csv"))
println("Begining Integration ...")
best_psms_dict = Dict{Int64, DataFrame}()
integration_time = @timed begin
MS_TABLE_PATH = MS_TABLE_PATHS[1]
ms_file_idx = 1
Threads.@threads for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS[1:1]))

    MS_TABLE = Arrow.Table(MS_TABLE_PATH)    
    
    #Need to change this. Highest intensity given a 1% fdr threshold. 
    #best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[ms_file_idx][PSMs[ms_file_idx][:,:q_value].<=0.1,:], :precursor_idx));

    #Get Best PSM per precursor
    best_psms = combine(sdf -> getBestPSM(sdf), groupby(PSMs[PSMs[:,:q_value].<=0.25,:], [:sequence,:charge]));
    
    transform!(best_psms, AsTable(:) => ByRow(psm -> 
                precursors_list[psm[:precursor_idx]].mz
                ) => :prec_mz
                );
    
    #Need to sort RTs 
    sort!(best_psms,:RT, rev = false);
    #Build RT index of precursors to integrate
    rt_index = buildRTIndex(best_psms);

    println("Integrating MS2...")
    ms2_chroms = integrateMS2(MS_TABLE, 
                                    frag_list, 
                                    rt_index,
                                    UInt32(ms_file_idx), 
                                    frag_err_dist_dict[ms_file_idx],
                                    integrate_ms2_params, 
                                    scan_range = (0, length(MS_TABLE[:scanNumber]))
                                    #can_range = (101357, 110357)
                                    );
    
    #Integrate MS2 Chromatograms 
    #Need to speed up this part. 
    transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms2_chroms, 
                                                UInt32(psm[:precursor_idx]), 
                                                (0.1f0, 0.15f0, 0.15f0, Float32(psm[:RT]), psm[:weight]), isplot = false)) => [:peak_area,:GOF,:FWHM,:FWHM_01,:asymmetry,:points_above_FWHM,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H]);
    
    #Remove Peaks with 0 MS2 intensity or fewer than 6 points accross the peak. 
    best_psms = best_psms[(ismissing.(best_psms[:,:peak_area]).==false).&(best_psms[:,:points_above_FWHM].>=1).&(best_psms[:,:points_above_FWHM_01].>=5),:];
    best_psms[:,:RT_error] = abs.(best_psms[:,:tᵣ] .- best_psms[:,:RT_pred]);

    #Get Predicted Isotope Distributions 
    #For some reason this requires a threadlock. Need to investiage further. 
    isotopes = UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}()
    #lock(lk) do 
    isotopes = getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], best_psms[:,:charge], QRoots(4), 4);
    #end

    prec_rt_table = sort(collect(zip(best_psms[:,:RT], UInt32.(best_psms[:,:precursor_idx]))), by = x->first(x));

    println("Integrating MS1...")
    ms1_chroms = integrateMS1(MS_TABLE, 
                                    isotopes, 
                                    prec_rt_table, 
                                    UInt32(ms_file_idx), 
                                    frag_err_dist_dict[ms_file_idx], 
                                    integrate_ms1_params,
                                    scan_range = (0, length(MS_TABLE[:scanNumber])));

    #Get MS1/MS2 Chromatogram Correlations and Offsets 
    #transform!(best_psms, AsTable(:) => ByRow(psm -> getCrossCorr(ms1_chroms, ms2_chroms, UInt32(psm[:precursor_idx]))) => [:offset,:cross_cor]);

    #Integrate MS1 Chromatograms 
    transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms1_chroms, UInt32(psm[:precursor_idx]), 
                                        (0.1f0, 
                                        0.15f0, 
                                        0.15f0, 
                                        Float32(psm[:tᵣ]), 
                                        Float32(1e4)),
                                        isplot = false)) => [:peak_area_ms1,
                                        :GOF_ms1,
                                        :FWHM_ms1,
                                        :FWHM_01_ms1,
                                        :asymmetry_ms1,
                                        :points_above_FWHM_ms1,
                                        :points_above_FWHM_01_ms1,
                                        :σ_ms1,:tᵣ_ms1,:τ_ms1,:H_ms1]);
    
    #QC Metrics
    best_psms = best_psms[best_psms[:,:peak_area].>0.0,:];
    best_psms = best_psms[best_psms[:,:FWHM].>(6/60),:];
    best_psms = best_psms[best_psms[:,:GOF].>0.5,:];

    lock(lk) do 
        best_psms_dict[ms_file_idx] = best_psms;
    end

end
end
integratePrecursor(ms1_chroms, UInt32(best_psms[N,:precursor_idx]), isplot = false)
N = N + 1

histogram(log2.(abs.(best_psms[best_psms[:,:decoy],:asymmetry])), alpha = 0.5)
histogram!(log2.(abs.(best_psms[best_psms[:,:decoy].==false,:asymmetry])), alpha = 0.5)

histogram(log2.(abs.(best_psms[best_psms[:,:decoy],:peak_area_ms1])), alpha = 0.5)
histogram!(log2.(abs.(best_psms[best_psms[:,:decoy].==false,:peak_area_ms1])), alpha = 0.5)

bins = LinRange(-2, 2, 200)
histogram((best_psms[best_psms[:,:decoy],:δt]), alpha = 0.5, normalize=:pdf, bins = bins)
histogram!((best_psms[best_psms[:,:decoy].==false,:δt]), alpha = 0.5, normalize = :pdf, bins = bins)

bins = LinRange(-2, 2, 200)
histogram((best_psms[best_psms[:,:decoy],:ρ]), alpha = 0.5, normalize=:pdf, bins = bins)
histogram!((best_psms[best_psms[:,:decoy].==false,:ρ]), alpha = 0.5, normalize = :pdf, bins = bins)


#bins = LinRange(-2, 2, 200)
histogram((best_psms[best_psms[:,:decoy],:points_above_FWHM_01]), alpha = 0.5, normalize=:probability)
histogram!((best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].<=0.01),:points_above_FWHM_01]), alpha = 0.5, normalize = :probability)

bins = LinRange(-2, 2, 100)
histogram2d((best_psms[(best_psms[:,:decoy].==true) .& (ismissing.(best_psms[:,:ρ]).==false),:δt]), 
            (best_psms[(best_psms[:,:decoy].==true) .& (ismissing.(best_psms[:,:ρ]).==false),:entropy_sim]), 
alpha = 0.5, normalize = :probability, bins = bins)

histogram2d((best_psms[(best_psms[:,:decoy].==false) .& (ismissing.(best_psms[:,:ρ]).==false),:δt]), 
            (best_psms[(best_psms[:,:decoy].==false) .& (ismissing.(best_psms[:,:ρ]).==false),:entropy_sim]), 
alpha = 0.5, normalize = :probability, bins = bins)

bins = LinRange(-2, 2, 100)
histogram(log2.(best_psms[(best_psms[:,:decoy].==true),:weight]), alpha = 0.5, normalize = :probability)
histogram!(log2.(best_psms[(best_psms[:,:decoy].==false),:weight]), alpha = 0.5, normalize = :probability)
histogram!(log2.(best_psms[(best_psms[:,:decoy].==false).& (best_psms[:,:q_value].<=0.01),:weight]), alpha = 0.5, normalize = :probability)
histogram!(best_psms[(best_psms[:,:decoy].==false),:δt], alpha = 0.5, normalize = :probability, bins = bins)

histogram((best_psms[(best_psms[:,:decoy].==true),:prec_ρ]), alpha = 0.5, normalize = :probability)
#histogram!((best_psms[(best_psms[:,:decoy].==false),:entropy_sim]), alpha = 0.5, normalize = :probability)
#histogram!(log2.(best_psms[(best_psms[:,:decoy].==false),:peak_area]), alpha = 0.5, normalize = :probability)
histogram!((best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].<=0.01),:prec_ρ]), alpha = 0.5, normalize = :probability)




histogram((best_psms[(best_psms[:,:decoy].==true),:prec_ρ]), alpha = 0.5, normalize = :probability)
histogram!((best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].<=0.01),:prec_ρ]), alpha = 0.5, normalize = :probability)

bins = LinRange(0.5, 1.0, 200)
histogram((best_psms[(best_psms[:,:decoy].==true),:GOF]), alpha = 0.5, normalize = :probability, bins = bins)
histogram!((best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].<=0.01),:GOF]), alpha = 0.5, normalize = :probability, bins = bins)

bins = (LinRange(0, 30, 100), LinRange(0, 30, 100))
histogram2d(log2.(best_psms[:,:peak_area]), log2.(max.(best_psms[:,:peak_area_ms1], 0.0)), alpha = 0.5, normalize = :probability, bins = bins)

#=
bins = (LinRange(0, 8, 100), LinRange(0, 8, 100))
histogram2d(log10.(best_psms[(best_psms[:,:decoy].==true) .& (best_psms[:,:q_value].<=0.01),:peak_area]), 
log10.(max.(best_psms[(best_psms[:,:decoy].==true) .& (best_psms[:,:q_value].<=0.01),:peak_area_ms1], 0.0)), alpha = 0.5, normalize = :probability, bins = bins)
=#

histogram!((best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].<=0.01),:n_precs]), alpha = 0.5, normalize = :probability)



histogram!(best_psms[(best_psms[:,:decoy].==false),:δt], alpha = 0.5, normalize = :probability, bins = bins)


best_psms[(best_psms[:,:decoy].==true).&(best_psms[:,:δt].==0.0),:]

println("Finished peak integration in ", integration_time.time, " seconds")
println("Train target-decoy model at precursor level...")
#best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),[:precursor_idx,:q_value,:matched_ratio,:entropy_sim,:intensity]]
best_psms = vcat(values(best_psms_dict)...);
#Model Features 
features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast,:RT_error,:y_ladder,:RT,:entropy_sim,:n_obs,:charge,:city_block,:matched_ratio,:scribe_score, :missed_cleavage,:Mox,:best_rank,:topn];
append!(features, [:peak_area,:GOF,:FWHM,:FWHM_01,:asymmetry,:points_above_FWHM,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H,:ρ,:δt,:peak_area_ms1,:points_above_FWHM_ms1]);

features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast,:RT_error,:RT,:entropy_sim,:charge,:city_block,:matched_ratio,:scribe_score, :missed_cleavage,:Mox,:best_rank,:topn];
append!(features, [:peak_area,:GOF,:asymmetry,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H,:ρ,:peak_area_ms1]);

best_psms = best_psms[(best_psms[:,:GOF].>=0.5),:];
  

#Train Model 
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 1.0, 
                        min_child_weight = 10, 
                        gamma = 10, 
                        subsample = 0.5, 
                        n_folds = 5, 
                        num_round = 200, 
                        max_depth = 10, 
                        eta = 0.0375, 
                        max_train_size = size(best_psms)[1]);
getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);


println("Number of unique Precursors ", length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx])))



sum((best_psms[:,:entropy_sim].>0.4).&(best_psms[:,:decoy]))
sum((best_psms[:,:entropy_sim].>0.4).&(best_psms[:,:decoy].==false))

MS_TABLE = Arrow.Table(MS_TABLE_PATHS[1])    
    
#Need to change this. Highest intensity given a 1% fdr threshold. 
best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[ms_file_idx][PSMs[ms_file_idx][:,:q_value].<=0.1,:], :precursor_idx));

transform!(best_psms, AsTable(:) => ByRow(psm -> 
            precursors_list[psm[:precursor_idx]].mz
            ) => :prec_mz
            );

#Need to sort RTs 
sort!(best_psms,:RT, rev = false);
#Build RT index of precursors to integrate
rt_index = buildRTIndex(best_psms);

println("Integrating MS2...")
ms2_chroms = integrateMS2(MS_TABLE, 
                                frag_list, 
                                rt_index,
                                UInt32(ms_file_idx), 
                                frag_err_dist_dict[ms_file_idx],
                                integrate_ms2_params, 
                                scan_range = (0, length(MS_TABLE[:scanNumber]))
                                #can_range = (101357, 110357)
                                );

ms2_chroms_square = integrateMS2(MS_TABLE, 
                                frag_list, 
                                rt_index,
                                UInt32(ms_file_idx), 
                                frag_err_dist_dict[ms_file_idx],
                                integrate_ms2_params, 
                                scan_range = (0, length(MS_TABLE[:scanNumber]))
                                #can_range = (101357, 110357)
                                );

N = 10065
best_psms_passing = best_psms[(best_psms[:,:q_value].<0.01) .& (best_psms[:,:decoy].==false),:]
N = 10065
N = 9936
#include("src/Routines/LibrarySearch/integrateChroms.jl")
integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
#integratePrecursor(ms2_chrom_square, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
N += 1

include("src/Routines/LibrarySearch/integrateChroms.jl")
ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
#Integrate MS2 Chromatograms 
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms2_chroms, UInt32(psm[:precursor_idx]), psm[:scan_idx], isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);


include("src/Routines/LibrarySearch/SearchRAW.jl")
fragment_intensities = integrateMS2(MS_TABLE, 
                                frag_list, 
                                rt_index,
                                UInt32(ms_file_idx), 
                                frag_err_dist_dict[ms_file_idx],
                                integrate_ms2_params, 
                                scan_range = (0, length(MS_TABLE[:scanNumber]))
                                #can_range = (101357, 110357)
                                );

p = plot()
#p = plot(title = "TEST", fontfamily="helvetica")
for (color, t) in enumerate(keys(fragment_intensities))
    plot!(p, fragment_intensities[t], color = color, legend = true, label = t)
    plot!(p, fragment_intensities[t], seriestype=:scatter, color = color, label = nothing, show = true)
end

Hs, X, weights, IDtoROW = integrateMS2(MS_TABLE, 
                                frag_list, 
                                rt_index,
                                UInt32(ms_file_idx), 
                                frag_err_dist_dict[ms_file_idx],
                                integrate_ms2_params, 
                                scan_range = (0, length(MS_TABLE[:scanNumber]))
                                #can_range = (101357, 110357)
                                );

solveHuber!(Hs, Hs*weights .- X, weights, Float32(1000), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);
Hs[:,17]
X[Hs[:,17].!=0.0]

Hs = Matrix(Hs)
Hs[2568,17] = 0.2
Hs = sparse(Hs)

test = frags_mouse_detailed_33NCEcorrected_start1[9301047]

for chrom in ProgressBar(ms2_chroms[2:end])
    chrom[:,:mz] = [MS_TABLE[:precursorMZ][scan] for scan in chrom[:, :scan_idx]]
end

transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(ms2_chroms, UInt32(psm[:precursor_idx]), psm[:scan_idx], isplot = false)) => [:intensity, :count, :SN, :slope, :peak_error,:apex,:fwhm]);
   

best_psms_passing = best_psms[(best_psms[:,:q_value].<0.01) .& (best_psms[:,:decoy].==false),:]
best_psms[(best_psms[:,:q_value].<0.01) .& (best_psms[:,:decoy].==false),:][10000:10010,[:sequence,:RT,:intensity,:prob]]
N = 10000
#best_psms_passing = best_psms[(best_psms[:,:q_value].<0.01) .& (best_psms[:,:decoy].==false),:]
integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
#ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
N += 1

p = [6e5, 30.2, 0.15]
m(t, p) = p[1]*exp.((-1).*((t .- p[2]).^2)./(2*p[3]^2))
m(ts, p)
plot(squared_error[:,:rt],
squared_error[:,:weight], seriestype=:scatter,
alpha = 0.5)
plot!(huber_loss[:,:rt],
huber_loss[:,:weight], seriestype=:scatter,
alpha = 0.5)
plot!(ts, m(ts, p))


include("src/Routines/LibrarySearch/searchRAW.jl")
ms2_chroms = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    #scan_range = (0, length(MS_TABLE[:scanNumber]))
    scan_range = (40000, 50000)
#scan_range = (101357, 102357)
);

ms2_chroms_huber = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    #scan_range = (0, length(MS_TABLE[:scanNumber]))
    scan_range = (40000, 50000)
#scan_range = (101357, 102357)
);
include("src/ML/sparseNNLS.jl")
include("src/Routines/LibrarySearch/searchRAW.jl")
ms2_chroms_huber_2 = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    #scan_range = (0, length(MS_TABLE[:scanNumber]))
    scan_range = (40000, 50000)
#scan_range = (101357, 102357)
);
include("src/Routines/LibrarySearch/searchRAW.jl")
include("src/ML/sparseNNLS.jl")
@time ms2_chroms_huber_3 = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    #scan_range = (0, length(MS_TABLE[:scanNumber]))
    scan_range = (40000, 70000)
#scan_range = (101357, 102357)
);
#integratePrecursor(ms2_chroms_huber, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
#integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)

N = 10045
N = 10250

N = 10296

N = 10045

#N = 10250
#N = 10143
N = 10303
squared_error = ms2_chroms_square[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
plot(squared_error[:,:rt],
squared_error[:,:weight], seriestype=:scatter,
alpha = 0.5)
#plot!(ms2_chroms_huber[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:rt],
#ms2_chroms_huber[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight], seriestype=:scatter,
#alpha = 0.5)
#plot!(ms2_chroms_huber_2[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:rt],
#ms2_chroms_huber_2[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight], seriestype=:scatter,
#alpha = 0.5)
huber_loss = ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
plot!(huber_loss[:,:rt],
huber_loss[:,:weight], seriestype=:scatter,
alpha = 0.5)
hline!([20000])
hline!([5000])
N += 1

julia> precursors_mouse_detailed_33NCEcorrected_start1[0x0078dd68]
LibraryPrecursor{Float32}(46.953674f0, 414.7449f0, 0.7946125f0, false, 0x02, 0x002ab610, UInt32[0x00000001], "13869", "TEALPGLK", 0x00, 0x01, 0x08)

julia> precursors_mouse_detailed_33NCEcorrected_start1[8952516]
LibraryPrecursor{Float32}(71.09662f0, 414.7449f0, 0.8327717f0, false, 0x02, 0x0029f55f, UInt32[0x00000001], "1412", "SVDLPGLK", 0x00, 0x01, 0x08)
plot(ms2_chroms[(precursor_idx=UInt32(0x0078dd68),)][:,:rt],
ms2_chroms[(precursor_idx=UInt32(0x0078dd68),)][:,:weight], seriestype=:scatter,
alpha = 0.5)
plot!(ms2_chroms_huber_3[(precursor_idx=UInt32(0x0078dd68),)][:,:rt],
ms2_chroms_huber_3[(precursor_idx=UInt32(0x0078dd68),)][:,:weight], seriestype=:scatter,
alpha = 0.5)

N += 1

cor(ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight],
ms2_chroms_huber_2[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight])

N += 1

integratePrecursor(ms2_chroms_huber, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)

N += 1

PSMs[1][PSMs[1][:,:precursor_idx] .== best_psms_passing[N,:precursor_idx],[:sequence,:precursor_idx,:weight,:total_ions,:best_rank,:entropy_sim,:matched_ratio,:spectral_contrast,:scribe_score,:RT,:q_value,:prob]]
#Why not found in ms2_chroms?
10097


plot(ms2_chroms[(precursor_idx=UInt32(best_psms_passing[10162,:precursor_idx]),)][16:29,:rt], 
ms2_chroms[(precursor_idx=UInt32(best_psms_passing[10162,:precursor_idx]),)][16:29,:weight],
seriestype=:scatter)


integratePrecursor(ms2_chroms, UInt32(best_psms_passing[10097,:precursor_idx]), isplot = true)



N = 10108

ms2_chroms[(precursor_idx=UInt32(best_psms_passing[10151,:precursor_idx]),)]


ms2_chroms[(precursor_idx=UInt32(best_psms_passing[10152,:precursor_idx]),)]
integratePrecursor(ms2_chroms, UInt32(best_psms_passing[10152,:precursor_idx]), isplot = true)
best_psms_passing[10152,:scan_idx]
precursor_idx = 889551 

integratePrecursor(ms2_chroms, UInt32(best_psms_passing[10152,:precursor_idx]),(best_psms_passing[10152,:scan_idx]), isplot = true)
best_psms_passing[10152,:scan_idx]
precursor_idx = 889551 


10162 #get rid of points with only one match
10166 #Tail when there are zeros? No walkback?
N = 10172 #Inteference?
10185 #Really bad MS2 bin overlap?
10211 #worse interference?
10230 #The worst interference?
10258 #Chosses wrong integration boundary. early cuttoff
N = 10346 #missed boundary entirely
N = 10191 # boundary
N = 10392 #overlap bin. Missed boundary?
N = 10402 #Two fragments matching. probably non_specific. #prec_idx = 1773834


PSMs[1][PSMs[1][:,:precursor_idx].==6874401,:]

plot!(ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:rt], 
ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight],
seriestype=:scatter)
MS_TABLE = Arrow.Table(MS_TABLE_PATHS[1])
ms_file_idx = 1
X, Hs, IDtoROW, last_matched_col =  firstSearch(
                                    MS_TABLE,
                                    frag_index,  
                                    frag_list, 
                                    x->x, #RT to iRT map'
                                    UInt32(ms_file_idx), #MS_FILE_IDX
                                    Laplace(zero(Float64), 10.0),
                                    first_search_params,
                                    scan_range = (101357, 101357)
                                    );