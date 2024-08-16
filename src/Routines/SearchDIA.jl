function SearchDIA(params_path::String)
    #params_path = "data/example_config/LibrarySearch.json"
 
     params = JSON.parse(read(params_path, String));
     #params = JSON.parse(read("data/example_config/LibrarySearch.json", String));
     MS_DATA_DIR = params["ms_data_dir"];
     SPEC_LIB_DIR = params["library_folder"];
 
     #Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
     MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
 
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


     ###########
     #Load Spectral Libraries
     ###########
     spec_lib = loadSpectralLibrary(SPEC_LIB_DIR)
     precursors = spec_lib["precursors"]
     ###########
     #Set CV Folds 
     ###########
     pid_to_cv_fold = getCVFolds(
         collect(range(UInt32(1), UInt32(length(spec_lib["precursors"][:sequence])))),#precursor id's, 
         spec_lib["precursors"][:accession_numbers]
         )
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
     file_id_to_parsed_name, parsed_fnames = parseFileNames(MS_TABLE_PATHS)
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
     println("Parameter Tuning Search...")
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

     #old_dict = copy(PSMs_Dict)
     ##########
     #Combine First Search Results
     ##########
     println("Combining First Search Results...")
     #Use high scoring psms to map library retention times (irt)
     #onto the empirical retention times and vise-versa.
     #uniform B-spline
     irt_rt, rt_irt = mapLibraryToEmpiricalRT(
        psms_paths,#Need to change this to a list of temporrary file paths
        min_prob = Float64(params_[:irt_mapping_params]["min_prob"])
     )
     #Summarize list of best N precursors accross all runs. 
     prec_to_irt = getBestPrecursorsAccrossRuns(
             psms_paths, 
             spec_lib["precursors"][:mz],
             rt_irt, 
             max_precursors = Int64(params_[:summarize_first_search_params]["max_precursors"]))
             map(x->(var_irt = x[:var_irt], n = x[:n]), prec_to_irt)

    irt_errs = getIrtErrs(
        peak_fwhms,
        prec_to_irt,
        params_
    )

    ms_table = Arrow.Table(first(MS_TABLE_PATHS))
    bin_rt_size = min(
    median(diff(ms_table[:retentionTime][ms_table[:msOrder].==1]))*5.5,
    params_[:summarize_first_search_params]["max_irt_bin_size"]
    )
    ms_table = nothing

    prec_to_irt = map(x->(irt = x[:best_irt], mz = x[:mz]), prec_to_irt)
    rt_index_paths = makeRTIndices(
         irt_indices_folder,
         psms_paths,
         prec_to_irt,
         rt_irt,
         min_prob = params_[:summarize_first_search_params]["max_prob_to_impute"])
     
     ############
     #Quantitative Search
     println("Begining Quantitative Search...")
     quantSearch(
         frag_err_dist_dict,
         pid_to_cv_fold,
         prec_to_irt,
         quant_psms_folder,
         rt_index_paths,
         bin_rt_size,
         rt_irt,
         irt_err,
         chromatograms,
         file_id_to_parsed_name,
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
         )

     println("Traning Target-Decoy Model...")

     best_psms = samplePSMsForXgboost(quant_psms_folder, params_[:xgboost_params]["max_n_samples"])
     scoreTraces!(best_psms,readdir(quant_psms_folder, join=true), precursors)
 
     #Get protein names 
     transform!(best_psms, AsTable(:) => ByRow(psm -> 
     precursors[:accession_numbers][psm[:precursor_idx]]
     ) => :accession_numbers
     );
     traces_passing = Set(best_psms[(best_psms[!,:q_value].<=params_[:q_value]).&(best_psms[!,:target]),:precursor_idx])
     ###########
     #Score Protein Groups
     scored_proteins = scoreProteinGroups!(best_psms)
     protein_to_q_value = Dict{Tuple{UInt32, String}, Float32}()
     for i in range(1, size(scored_proteins, 1))
         protein_to_q_value[
             (
                 UInt32(scored_proteins[i,:ms_file_idx]),
                 scored_proteins[i,:accession_numbers]
             )
         ] = scored_proteins[i,:q_value]
     end
     ###########
     #Re-quantify with 1% fdr precursors 
     best_psms[!,:peak_area] = zeros(Float32, size(best_psms, 1))
     best_psms[!,:new_best_scan] = zeros(UInt32, size(best_psms, 1))
     grouped_best_psms = groupby(best_psms, :file_name)
     println("Re-Quantifying Precursors at FDR Threshold...")
     secondQuant!( 
                     grouped_best_psms,   
                     frag_err_dist_dict,
                     traces_passing,
                     rt_index_paths,
                     RT_iRT,
                     irt_err,
                     chromatograms,
                     file_id_to_parsed_name,
                     MS_TABLE_PATHS,
                     params_,
                     spec_lib["precursors"],
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
 
     ###########
     #Ungroup precursors and get the best trace for each 
     best_psms = DataFrame(grouped_best_psms)
     #Must be based on weight and not peak_area because peak_area was not calculated for decoys
     getBestTrace!(best_psms, params_[:q_value], :weight)
     filter!(x->x.best_trace, best_psms)
     #precursor level q_value 
     getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
     #filter out decoys and low-scoring targets
     filter!(x->(x.q_value<=params_[:q_value])&(x.target)&(!isnan(x.peak_area)), best_psms)
     #IDs_PER_FILE = value_counts(best_psms, [:file_name])
     best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
     ###########
     #Normalize Quant 
     ###########
     println("Normalization...")
     normalizeQuant(
             best_psms,
             :peak_area,
             N = params_[:normalization_params]["n_rt_bins"],
             spline_n_knots = params_[:normalization_params]["spline_n_knots"],
             max_q_value = params_[:normalization_params]["max_q_value"],
             min_points_above_FWHM = params_[:normalization_params]["min_points_above_FWHM"]
         )
     
     ############
     #Prep for Protein Inference 
     best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
     best_psms[!,:ms_file_idx] =  UInt32.(best_psms[!,:ms_file_idx])
     best_psms[!,:peak_area] =  allowmissing(best_psms[!,:peak_area])
     best_psms[!,:peak_area_normalized] =  allowmissing(best_psms[!,:peak_area_normalized])
 
     gbpsms = groupby(best_psms,:ms_file_idx)
     file_idx_dicts = [Dict{UInt32, @NamedTuple{prob::Float32, qvalue::Float32, peak_area::Float32}}() for _ in range(1, length(gbpsms))]
     for (file_idx, bpsms) in pairs(gbpsms)
         for i in range(1, size(bpsms, 1))
                 file_idx_dicts[file_idx[:ms_file_idx]][bpsms[i,:precursor_idx]] = (
                 prob = bpsms[i,:prob],
                 qvalue = bpsms[i,:q_value],
                 peak_area = bpsms[i,:peak_area])
         end
     end
     best_psms[!,:structural_mods] = [precursors[:structural_mods][pid] for pid in best_psms[!,:precursor_idx]]
     best_psms[!,:isotopic_mods] = [precursors[:isotopic_mods][pid] for pid in best_psms[!,:precursor_idx]]
 
     features = [:species,:accession_numbers,:sequence,:structural_mods,
                 :isotopic_mods,:charge,:precursor_idx,:target,:weight,:peak_area,:peak_area_normalized,
                 :scan_idx,:prob,:q_value,:prec_mz,:RT,:irt_obs,:irt_pred,:best_rank,:best_rank_iso,:topn,:topn_iso,:longest_y,:longest_b,:b_count,
                 :y_count,:p_count,:non_cannonical_count,:isotope_count
                 ,:matched_ratio]
     println("Writting PSMs...")
     sort!(best_psms,[:species,:accession_numbers,:sequence,:target])
     Arrow.write(joinpath(results_folder,"psms_long.arrow"),best_psms[!,features])
     wide_psms_quant = unstack(best_psms,[:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],:file_name,:peak_area_normalized)
     sort!(wide_psms_quant,[:species,:accession_numbers,:sequence,:target])
     CSV.write(joinpath(results_folder,"psms_wide.csv"),wide_psms_quant)
 
     #Summarize Precursor ID's
     value_counts(df, col) = combine(groupby(df, col), nrow)
 
     precursor_id_table = value_counts(best_psms,:file_name)
     #CSV.write("/Users/n.t.wamsley/Desktop/precursor_ids_table.csv")
     println("Max LFQ...")
     protein_quant = proteinQuant(
         best_psms,
         protein_to_q_value
     )
 
     protein_quant[!,:experiments] = UInt32.(protein_quant[!,:experiments])
     protein_quant[!,:file_name] = [file_id_to_parsed_name[ms_file_idx] for ms_file_idx in protein_quant[!,:experiments]]
 
     #Wide DataFormat 
     wide_protein_quant = unstack(protein_quant,[:species,:protein,:target],:experiments,:log2_abundance)
     sort!(wide_protein_quant,[:species,:protein,:target])
     CSV.write(joinpath(results_folder,"proteins_wide.csv"),wide_protein_quant)
     println("QC Plots")
     qcPlots(
         best_psms,
         protein_quant,
         params_,
         parsed_fnames,
         qc_plot_folder,
         file_id_to_parsed_name,
         MS_TABLE_PATHS,
         iRT_RT,
         frag_err_dist_dict
     )
end