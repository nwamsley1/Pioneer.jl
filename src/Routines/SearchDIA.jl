function SearchDIA(params_path::String)
    #println("JLD2 version is: ", Pkg.installed()["JLD2"])
    #println("TESTTESTTEST")
    total_time = @timed begin
    #params_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/OlsenMixedPrositParams.json"
    if !isabspath(params_path)
        params_path = joinpath(@__DIR__, "../../", params_path)
    end
    params = JSON.parse(read(params_path, String));
    MS_DATA_DIR = params["ms_data_dir"];
    SPEC_LIB_DIR = params["library_folder"];
    if !isabspath(SPEC_LIB_DIR)
        SPEC_LIB_DIR =  joinpath(@__DIR__, "../../", SPEC_LIB_DIR)
    end

    #Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
    if isabspath(MS_DATA_DIR)
        MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
    else
        println("B")
        MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR)
        MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
    end

    params_ = parseParams(params)

    qc_plot_folder, rt_alignment_folder, mass_err_estimation_folder, results_folder, temp_folder = makeOutputDirectories(
        joinpath(params_[:benchmark_params]["results_folder"], "RESULTS"),
        params
    )

    first_search_psms_folder = joinpath(temp_folder, "first_search_psms")
    if !isdir(first_search_psms_folder)
        mkpath(first_search_psms_folder)
    end

    irt_indices_folder = joinpath(temp_folder, "irt_indices_folder")
    if !isdir(irt_indices_folder)
        mkpath(irt_indices_folder)
    end

    quant_psms_folder = joinpath(temp_folder, "quant_psms_folder")
    if !isdir(quant_psms_folder )
        mkpath(quant_psms_folder )
    end

    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    if !isdir( passing_psms_folder )
    mkpath( passing_psms_folder )
    end

    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    if !isdir( passing_proteins_folder )
       mkpath( passing_proteins_folder )
    end

    second_quant_folder = joinpath(temp_folder, "second_quant")
    if !isdir( second_quant_folder)
        mkpath( second_quant_folder)
    end
    ###########
    #Load Spectral Libraries
    ###########
    spec_lib = loadSpectralLibrary(SPEC_LIB_DIR);
    precursors = spec_lib["precursors"];
    unique_proteins = unique(precursors[:accession_numbers]);
    accession_number_to_id = Dict(zip(unique_proteins, range(one(UInt32), UInt32(length(unique_proteins)))));
    ###########
    #Set CV Folds 
    ###########
    pid_to_cv_fold = getCVFolds(
        collect(range(UInt32(1), UInt32(length(spec_lib["precursors"][:sequence])))),#precursor id's, 
        spec_lib["precursors"][:accession_numbers]
        );
    ###########
    #Load Pre-Allocated Data Structures. One of each for each thread. 
    ###########
    N = Threads.nthreads()
    M = 250000
    n_precursors = length(spec_lib["precursors"][:mz])
    ionMatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    ionMisses = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    all_fmatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors ) for _ in range(1, N)];
    ionTemplates = [[DetailedFrag{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    iso_splines = parseIsoXML(joinpath(@__DIR__,"../../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml"));
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
    #File Names parsing
    file_id_to_parsed_name, parsed_fnames,file_path_to_parsed_name = parseFileNames(MS_TABLE_PATHS)

    ###########
    #Tune Parameters
    println("Parameter Tuning Search...")
    RT_to_iRT_map_dict, frag_err_dist_dict, irt_errs = parameterTuningSearch(rt_alignment_folder,
                                                                            mass_err_estimation_folder,
                                                                            MS_TABLE_PATHS,
                                                                            params_,
                                                                            spec_lib,
                                                                            ionMatches,
                                                                            ionMisses,
                                                                            all_fmatches,
                                                                            IDtoCOL,
                                                                            ionTemplates,
                                                                            iso_splines,
                                                                            scored_PSMs,
                                                                            unscored_PSMs,
                                                                            spectral_scores,
                                                                            precs)

    peak_fwhms, psms_paths = firstSearch(
        first_search_psms_folder,
        RT_to_iRT_map_dict,
        frag_err_dist_dict,
        irt_errs,
        file_id_to_parsed_name,
        MS_TABLE_PATHS,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precs
    )
    #first_psms = DataFrame(Arrow.Table([fname for fname in readdir(first_search_psms_folder,join=true) if endswith(fname, ".arrow")]))
    #first_psms = first_psms[first_psms[!,:q_value].<=0.01,:]
    #value_counts(df, col) = combine(groupby(df, col), nrow)
    #psms_counts = value_counts(first_psms, :ms_file_idx)
    #CSV.write(joinpath(results_folder, "first_search_psms_counts.csv"), psms_counts)
    #first_psms = nothing
    #test_table = DataFrame(Arrow.Table(collect(values(psms_paths))))
    #test_table = test_table[test_table[!,:q_value].<=0.01,:]
    #value_counts(test_table,:ms_file_idx)
    ##########
    #Combine First Search Results
    ##########
    println("Combining First Search Results...")
    #Use high scoring psms to map library retention times (irt)
    #onto the empirical retention times and vise-versa.
    #uniform B-spline
    #println("typeof(RT_to_iRT_map_dict) ", typeof(RT_to_iRT_map_dict))
    #println("RT_to_iRT_map_dict ", RT_to_iRT_map_dict)
    irt_rt, rt_irt = mapLibraryToEmpiricalRT(
    psms_paths,#Need to change this to a list of temporrary file paths]
    RT_to_iRT_map_dict,
    min_prob = Float64(params_[:irt_mapping_params]["min_prob"])
    )
    #Summarize list of best N precursors accross all runs. 
    precursor_dict = getBestPrecursorsAccrossRuns(
            psms_paths, 
            spec_lib["precursors"][:mz],
            rt_irt, 
            max_q_val = Float32(params_[:summarize_first_search_params]["max_q_val_for_irt"]),
            max_precursors = Int64(params_[:summarize_first_search_params]["max_precursors"]));
    irt_errs = nothing
    if length(keys(peak_fwhms)) > 1
        irt_errs = getIrtErrs(
            peak_fwhms,
            precursor_dict,
            params_
        )
    else
        irt_errs = Dict(zip(keys(peak_fwhms), zeros(Float32, length(keys(peak_fwhms)))))
    end
    bin_rt_size =  params_[:summarize_first_search_params]["max_irt_bin_size"]
    prec_to_irt =  map(x->(irt = x[:best_irt], mz = x[:mz]), precursor_dict)
    rt_index_paths = makeRTIndices(
         irt_indices_folder,
         psms_paths,
         prec_to_irt,
         rt_irt,
         min_prob = params_[:summarize_first_search_params]["max_prob_to_impute"])
     
    ############
    #Estimate Optimal Huber Loss Smoothing Parameter
    delta0 = params_[:deconvolution_params]["huber_delta0"];
    delta_exp = params_[:deconvolution_params]["huber_delta_exp"];
    delta_iters = params_[:deconvolution_params]["huber_delta_iters"];
    huber_δs = Float32[delta0*(delta_exp^i) for i in range(1, delta_iters)];

    println("Optimize Huber Loss Smoothing Parameter...")
    params_[:deconvolution_params]["huber_delta"] = getHuberLossParam(
        huber_δs,
        precursor_dict,
        MS_TABLE_PATHS,
        precursors[:is_decoy],
        frag_err_dist_dict,
        rt_index_paths,
        bin_rt_size,
        rt_irt,
        irt_errs,
        chromatograms,
        file_path_to_parsed_name,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        complex_scored_PSMs,
        complex_unscored_PSMs,
        complex_spectral_scores,
        precursor_weights);

    #if params_[:output_params]["delete_temp"]
    #    rm(first_search_psms_folder,recursive=true)
    #end      
        println("huber ", params_[:deconvolution_params]["huber_delta"])
    ############
    #Quantitative Search
    println("Begining Quantitative Search...")
    ms_table_path_to_psms_path = quantSearch(
        frag_err_dist_dict,
        pid_to_cv_fold,
        prec_to_irt,
        quant_psms_folder,
        rt_index_paths,
        bin_rt_size,
        rt_irt,
        irt_errs,
        chromatograms,
        file_path_to_parsed_name,
        MS_TABLE_PATHS,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        complex_scored_PSMs,
        complex_unscored_PSMs,
        complex_spectral_scores,
        precursor_weights
        );

    println("Traning Target-Decoy Model...")
    best_psms = samplePSMsForXgboost(quant_psms_folder, params_[:xgboost_params]["max_n_samples"]);
    models = scoreTraces!(best_psms,readdir(quant_psms_folder, join=true), precursors);
    #Arrow.write("/Users/n.t.wamsley/Desktop/test_psms_quant.arrow", best_psms)
    #Wipe memory
    best_psms = nothing
    GC.gc()
    best_traces = getBestTraces(
        quant_psms_folder,
        Float32(params_[:xgboost_params]["min_best_trace_prob"]));
    #Path for merged quant psms scores 
    merged_quant_path = joinpath(temp_folder, "merged_quant.arrow")
    #Sort quant tables in descenging order of probability and remove 
    #sub-optimal isotope traces 
    sortAndFilterQuantTables(
        quant_psms_folder,
        merged_quant_path,
        best_traces
    )
    #Merge sorted tables into a single list with two columns
    #for "prob" and "target"
    mergeSortedPSMScores(
                    quant_psms_folder, 
                    merged_quant_path
                    )
    #functions to convert "prob" to q-value and posterior-error-probability "pep" 
    precursor_pep_spline = getPEPSpline(merged_quant_path, 
                                        :prob, 
                                        min_pep_points_per_bin = params_[:xgboost_params]["precursor_prob_spline_points_per_bin"], 
                                        n_spline_bins = 5)
    precursor_qval_interp = getQValueSpline(merged_quant_path, 
                                            :prob, 
                                            min_pep_points_per_bin = max(params_[:xgboost_params]["precursor_q_value_interpolation_points_per_bin"], 3))

    getPSMsPassingQVal(
                                    quant_psms_folder, 
                                    passing_psms_folder,
                                    precursor_pep_spline,
                                    precursor_qval_interp,
                                    0.01f0)
    #Delete Quant PSMs
    #if params_[:output_params]["delete_temp"]
    #    rm(quant_psms_folder,recursive=true)
    #end
    #=
        quant_psms_test = DataFrame(Arrow.Table("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/RESULTS/temp/merged_quant.arrow"))
        target_probs = quant_psms_test[quant_psms_test[!,:target],:prob]
        decoy_probs = quant_psms_test[quant_psms_test[!,:target].==false,:prob]
        tbins = LinRange(0, 1, 100)
        histogram(target_probs, normalize=:probability, bins=tbins, alpha = 0.5)
        histogram!(decoy_probs, normalize=:probability, bins=tbins, alpha = 0.5)



        quant_psms_test = DataFrame(Arrow.Table(readdir("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/RESULTS/temp/quant_psms_folder", join=true)))
        target_probs = quant_psms_test[quant_psms_test[!,:target],:max_unmatched_residual]
        decoy_probs = quant_psms_test[quant_psms_test[!,:target].==false,:max_unmatched_residual]
        tbins = LinRange(minimum(target_probs), maximum(target_probs), 100)
        histogram(target_probs, normalize=:probability, bins=tbins, alpha = 0.5)
        histogram!(decoy_probs, normalize=:probability, bins=tbins, alpha = 0.5)


        target_probs = quant_psms_test[quant_psms_test[!,:target],:fitted_manhattan_distance]
        decoy_probs = quant_psms_test[quant_psms_test[!,:target].==false,:fitted_manhattan_distance]
        tbins = LinRange(minimum(target_probs), maximum(target_probs), 100)
        histogram(target_probs, normalize=:probability, bins=tbins, alpha = 0.5)
        histogram!(decoy_probs, normalize=:probability, bins=tbins, alpha = 0.5)


        quant_psms_test = DataFrame(Arrow.Table(readdir("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/RESULTS/temp/first_search_psms", join=true)))
        quant_psms_test[!,:target] = [precursors[:target][pid] for pid in quant_psms_test[!,:precursor_idx]]
        target_probs = quant_psms_test[quant_psms_test[!,:target],:prob]
        tbins = LinRange(0, 1, 100)
        histogram(target_probs, normalize=:probability, bins=tbins, alpha = 0.5)
        histogram!(decoy_probs, normalize=:probability, bins=tbins, alpha = 0.5)

    =#
    ###########
    #Score Protein Groups
    sorted_pg_score_path = getProteinGroups(
            passing_psms_folder, 
            passing_proteins_folder,
            temp_folder,
            precursors[:accession_numbers],
            accession_number_to_id,
            precursors[:sequence])
    if params_[:output_params]["delete_temp"]
        rm(passing_proteins_folder,recursive=true)
    end
    #=
    pg_pep_spline = getPEPSpline(sorted_pg_score_path, 
                                :max_pg_score, 
                                min_pep_points_per_bin = params_[:xgboost_params]["pg_prob_spline_points_per_bin"], 
                                n_spline_bins = 5)
    =#
    pg_qval_interp = getQValueSpline(sorted_pg_score_path, 
                                    :max_pg_score, 
                                    min_pep_points_per_bin = max(params_[:xgboost_params]["pg_q_value_interpolation_points_per_bin"], 3))

    ###########
    #Re-quantify with 1% fdr precursors 
    println("Re-Quantifying Precursors at FDR Threshold...")
    test_prop = secondQuantSearch!( 
                    file_path_to_parsed_name,
                    passing_psms_folder,
                    second_quant_folder,
                    frag_err_dist_dict,
                    rt_index_paths,
                    bin_rt_size,
                    rt_irt,
                    irt_errs,
                    chromatograms,
                    MS_TABLE_PATHS,
                    params_,
                    spec_lib["precursors"],
                    accession_number_to_id,
                    spec_lib,
                    ionMatches,
                    ionMisses,
                    IDtoCOL,
                    ionTemplates,
                    iso_splines,
                    complex_scored_PSMs,
                    complex_unscored_PSMs,
                    complex_spectral_scores,
                    precursor_weights)

    if params_[:output_params]["delete_temp"]
        rm(passing_psms_folder,recursive=true)
    end
    ###########
    #Normalize Quant 
    ###########
    println("Normalization...")
    normalizeQuant(
            second_quant_folder,
            :peak_area,
            N = params_[:normalization_params]["n_rt_bins"],
            spline_n_knots = params_[:normalization_params]["spline_n_knots"]
        )

    merged_second_quant_path = joinpath(temp_folder, "joined_second_quant.arrow")
    if isfile(merged_second_quant_path)
        rm(merged_second_quant_path, force = true)
    end
    if isfile(joinpath(results_folder, "precursors_long.arrow"))
        rm(joinpath(results_folder, "precursors_long.arrow"), force=true)
    end
    mergeSortedArrowTables(
        second_quant_folder,
        joinpath(results_folder, "precursors_long.arrow"),
        (:protein_idx,:precursor_idx),
        N = 1000000
    )
    if params_[:output_params]["delete_temp"]
        rm(second_quant_folder,recursive=true)
    end
    precursors_wide_arrow = writePrecursorCSV(
        joinpath(results_folder, "precursors_long.arrow"),
        sort(collect(values(file_id_to_parsed_name))),
        false,#normalized
        write_csv = params_[:output_params]["write_csv"],
        )
     #Summarize Precursor ID's
     #value_counts(df, col) = combine(groupby(df, col), nrow)
     println("Max LFQ...")
    if isfile( joinpath(results_folder, "protein_groups_long.arrow"))
        rm( joinpath(results_folder, "protein_groups_long.arrow"), force=true)
    end
    LFQ(
        DataFrame(Arrow.Table(joinpath(results_folder, "precursors_long.arrow"))),
        joinpath(results_folder, "protein_groups_long.arrow"),
        :peak_area,
        file_id_to_parsed_name,
        0.01f0,
        pg_qval_interp,
        batch_size = 100000
    )
    protein_groups_wide_arrow = writeProteinGroupsCSV(
        joinpath(results_folder, "protein_groups_long.arrow"),
        precursors[:sequence],
        precursors[:isotopic_mods],
        precursors[:structural_mods],
        precursors[:prec_charge],
        sort(collect(values(file_id_to_parsed_name))),
        write_csv = params_[:output_params]["write_csv"],
        )

        #=
        new_path = "/Users/n.t.wamsley/Desktop/astral_test_AltimeterFirstTry"
        mkdir(new_path)
        ThreeProteomeAnalysis(
        "/Users/n.t.wamsley/Desktop/testresults/RESULTS/RESULTS",
            "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt",
            new_path
        )

                ThreeProteomeAnalysis(
        "/Users/n.t.wamsley/Desktop/testresults/RESULTS/RESULTS",
            "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt",
            "/Users/n.t.wamsley/Desktop/astral_test_out_nonnorm"
        )

        =#
    println("QC Plots")
    if isfile(joinpath(qc_plot_folder, "QC_PLOTS.pdf"))
        rm(joinpath(qc_plot_folder, "QC_PLOTS.pdf"))
    end
    qcPlots(
        precursors_wide_arrow,
        protein_groups_wide_arrow,
        params_,
        precursors,
        parsed_fnames,
        qc_plot_folder,
        file_id_to_parsed_name,
        MS_TABLE_PATHS,
        irt_rt,
        frag_err_dist_dict
    )
    if params_[:output_params]["delete_temp"]
        rm(temp_folder,recursive=true)
    end
end
    return 10
end
