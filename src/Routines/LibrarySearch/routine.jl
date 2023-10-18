
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
SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_indy4b3_ally2b1/"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATHS = MS_TABLE_PATHS[1:1]
EXPERIMENT_NAME = "TEST_2"
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
    b_min_ind = Int64(params["b_min_ind"]),
    y_min_ind = Int64(params["y_min_ind"]),
    expected_matches = Int64(params["expected_matches"]),
    frag_ppm_err = Float64(params["frag_ppm_err"]),
    frag_tol_quantile = Float32(params["frag_tol_quantile"]),
    frag_tol_presearch = Float64(params["frag_tol_presearch"]),
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
    most_intense = Bool(params["most_intense"]),
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
    :frag_tol_presearch => params_[:frag_tol_presearch],
    :max_iter => params_[:nnls_max_iter],
    :max_peaks => params_[:max_peaks],
    :min_frag_count => params_[:min_frag_count_presearch],
    :min_frag_count_index_search => params_[:min_frag_count_index_search],
    :min_matched_ratio => params_[:min_matched_ratio_presearch],
    :min_matched_ratio_index_search => params_[:min_matched_ratio_index_search],
    :min_spectral_contrast => params_[:min_spectral_contrast_presearch],
    :most_intense => params_[:most_intense],
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
    :most_intense => params_[:most_intense],
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
f_det_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_det\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
spec_load_time = @timed begin
    f_index = load(f_index_path);
    prosit_lib["f_index"] = f_index["f_index"]
    f_det = load(f_det_path)
    prosit_lib["f_det"] = f_det["f_det_new"];
    precursors = load(precursors_path)
    prosit_lib["precursors"] = precursors["precursors"];
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
                                        Laplace(zero(Float64), first_search_params[:frag_tol_presearch]),
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

        println("Finished main search for $ms_file_idx in ", sub_search_time.time, " seconds")
        #@save "/Users/n.t.wamsley/TEST_DATA/PSMs_TEST_101623.jld2" PSMs;=
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
        lock(lk) do 
            isotopes = getIsotopes(PSMs[!,:sequence], 
                                    PSMs[!,:precursor_idx], 
                                    PSMs[!,:charge], QRoots(4), 4);
        end
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
        PSMs[:,"sequence_length"] = length.(replace.(PSMs[:,"sequence"], "M(ox)" => "M"));

        for i in range(1, size(PSMs)[1])
            if ismissing(PSMs[i,:ρ])
                continue
            end
            if isnan(PSMs[i,:ρ])
                PSMs[i,:ρ] = missing
            end
        end

        lock(lk) do 
            println("ms_file_idx: $ms_file_idx has ", size(PSMs))
            BPSMS[ms_file_idx] = PSMs
        end
end
println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")

best_psms = vcat(values(BPSMS)...)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms.jld2"); best_psms)

open(joinpath(MS_DATA_DIR, "Search", "RESULTS", "meta.txt"), "w") do f
    write(f, "Main Search Time $main_search_time \n")
end

open(joinpath(MS_DATA_DIR, "Search", "PARAMS","params.json"), "w") do f
    write(f, JSON.json(JSON.parse(read(ARGS["params_json"], String))))
end

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