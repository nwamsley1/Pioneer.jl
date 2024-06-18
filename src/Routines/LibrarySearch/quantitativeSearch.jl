function quantitationSearch(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    params::Any;
                    kwargs...)
    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    thread_tasks, total_peaks = partitionScansToThreads(spectra[:masses],
                                                        spectra[:retentionTime],
                                                        spectra[:centerMass],
                                                        spectra[:msOrder],

                                                        Threads.nthreads(),
                                                        1)
    
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return quantPSMs(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                kwargs[:fragment_lookup_table], 
                                kwargs[:ms_file_idx],
                                kwargs[:rt_to_irt_spline],
                                kwargs[:mass_err_model],
                                Float32(params[:deconvolution_params]["huber_delta"]),
                                Float32(params[:deconvolution_params]["lambda"]),
                                Int64(params[:deconvolution_params]["max_iter_newton"]),
                                Int64(params[:deconvolution_params]["max_iter_bisection"]),
                                Int64(params[:deconvolution_params]["max_iter_outer"]),
                                Float32(params[:deconvolution_params]["accuracy_newton"]),
                                Float32(params[:deconvolution_params]["accuracy_bisection"]),
                                Float32(params[:deconvolution_params]["max_diff"]),
                                kwargs[:ion_matches][thread_id],
                                kwargs[:ion_misses][thread_id],
                                kwargs[:id_to_col][thread_id],
                                kwargs[:ion_templates][thread_id],
                                kwargs[:iso_splines],
                                kwargs[:chromatograms][thread_id],
                                kwargs[:scored_psms][thread_id],
                                kwargs[:unscored_psms][thread_id],
                                kwargs[:spectral_scores][thread_id],
                                kwargs[:precursor_weights][thread_id],
                                (3, 0),
                                Int64(params[:quant_search_params]["min_frag_count"]),
                                Float32(params[:quant_search_params]["min_spectral_contrast"]),
                                Float32(params[:quant_search_params]["min_log2_matched_ratio"]),
                                Tuple([Int64(x) for x in params[:quant_search_params]["min_topn_of_m"]]),
                                Int64(params[:quant_search_params]["max_best_rank"]),
                                params[:quant_search_params]["n_frag_isotopes"],
                                kwargs[:rt_index], 
                                kwargs[:irt_err],
                                Set(2),
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end
BPSMS = Dict{Int64, DataFrame}()
#PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
#PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];



quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    params_[:deconvolution_params]["huber_delta"] = median(
        [quantile(x, 0.25) for x in MS_TABLE[:intensities]])*params_[:deconvolution_params]["huber_delta_prop"] 
       

        #params_[:deconvolution_params]["huber_delta"] = 100.0f0
        params_[:deconvolution_params]["lambda"] = 0.0f0
        params_[:deconvolution_params]["accuracy_bisection"] = 10.0
        params_[:deconvolution_params]["accuracy_newton"] = 10.0
       

        @time RESULT = quantitationSearch(
            MS_TABLE, 
            params_;
            precursors = prosit_lib["precursors"],
            fragment_lookup_table = library_fragment_lookup_table,
            rt_index = RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
            ms_file_idx = UInt32(ms_file_idx), 
            rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
            mass_err_model = frag_err_dist_dict[ms_file_idx],
            irt_err = irt_err,#irt_errs[ms_file_idx]/3,
            ion_matches = ionMatches,
            ion_misses = ionMisses,
            id_to_col = IDtoCOL,
            ion_templates = ionTemplates,
            iso_splines = iso_splines,
            chromatograms = chromatograms,
            scored_psms = complex_scored_PSMs,
            unscored_psms = complex_unscored_PSMs,
            spectral_scores = complex_spectral_scores,
            precursor_weights = precursor_weights,
            );

        #Format Chromatograms 
        chroms = vcat([last(x) for x in RESULT]...);
        #sort!(chroms,:rt,alg=QuickSort)
        getIsotopesCaptured!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        gchroms = groupby(chroms,[:precursor_idx])

        #Format spectral scores 
        psms = vcat([first(x) for x in RESULT]...);
        addSecondSearchColumns!(psms, 
                                        MS_TABLE, 
                                        prosit_lib["precursors"][:mz],
                                        prosit_lib["precursors"][:prec_charge], 
                                        prosit_lib["precursors"][:is_decoy],
                                        pid_to_cv_fold);
        psms[!,:charge2] = UInt8.(psms[!,:charge].==2)
        filter!(x->isnan(x.entropy_score)==false, psms)
        features = [:intercept, :charge2, :total_ions, :err_norm, 
        :scribe, :city_block, :city_block_fitted, 
        :spectral_contrast, :entropy_score, :weight]
       # if sum(psms[!,:p_count])>0
       #     push!(features, :p_count)
       # end
        getIsotopesCaptured!(psms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        psms[!,:prob] = zeros(Float32, size(psms, 1));
        scoreSecondSearchPSMs!(psms,features);
        spectral_scores = groupby(psms, [:precursor_idx,:isotopes_captured]);

        #Integrate Chromatograms and Combine with Spectral Simmilarity Scores
        precs_to_integrate = keys(spectral_scores) ∩ keys(gchroms)
        n_precs = length(precs_to_integrate)
        psms_integrated = initQuantScans(precs_to_integrate)


        integratePrecursors(spectral_scores,
                            gchroms,
                            precs_to_integrate,
                            psms_integrated,
                            λ=Float32(params_[:quant_search_params]["WH_smoothing_strength"]))
        #Remove examples where maximum point on chromatogram doesn't have a psm
        filter!(x->!ismissing(x.scan_idx), psms_integrated)
        #Add additional features 
        addPostIntegrationFeatures!(psms_integrated, 
                                    MS_TABLE, 
                                    prosit_lib["precursors"][:sequence],
                                    prosit_lib["precursors"][:structural_mods],
                                    prosit_lib["precursors"][:mz],
                                    prosit_lib["precursors"][:irt],
                                    prosit_lib["precursors"][:prec_charge],
                                    prosit_lib["precursors"][:missed_cleavages],
                                    ms_file_idx,
                                    file_id_to_parsed_name,
                                    RT_iRT,
                                    precID_to_iRT
                                    );
        psms_integrated[!,:file_name].=file_id_to_parsed_name[ms_file_idx]
        BPSMS[ms_file_idx] = psms_integrated;
end

best_psms = vcat(values(BPSMS)...)