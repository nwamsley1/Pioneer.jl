
function quantSearch(
    frag_err_dist_dict,
    pid_to_cv_fold,
    precID_to_iRT,
    quant_psms_folder,
    rt_index_paths,
    bin_rt_size,
    RT_iRT,
    irt_errs,
    quad_model_dict,
    isotope_trace_type,
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
    complex_scored_psms,
    complex_unscored_psms,
    complex_spectral_scores,
    precursor_weights
    )

    function quantSearch(
                        #Mandatory Args
                        spectra::Arrow.Table, 
                        params::Any;
                        kwargs...)
        ########
        #Each thread needs to handle a similair number of peaks. 
        #For example if there are 10,000 scans and two threads, choose n so that
        #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
        #of scans have an equal number of fragment peaks in the spectra
        thread_tasks, total_peaks = partitionScansToThreads(spectra[:mz_array],
                                                            spectra[:retentionTime],
                                                            spectra[:centerMz],
                                                            spectra[:msOrder],

                                                            Threads.nthreads(),
                                                            1)
        
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin 
                thread_id = first(thread_task)
                return secondSearch(
                                    spectra,
                                    last(thread_task), #getRange(thread_task),
                                    spec_lib["precursors"],
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
                                    kwargs[:scored_psms][thread_id],
                                    kwargs[:unscored_psms][thread_id],
                                    kwargs[:spectral_scores][thread_id],
                                    kwargs[:precursor_weights][thread_id],
                                    kwargs[:quad_transmission_model],
                                    (3, 0),
                                    Int64(params[:quant_search_params]["min_y_count"]),
                                    Int64(params[:quant_search_params]["min_frag_count"]),
                                    Float32(params[:quant_search_params]["min_spectral_contrast"]),
                                    Float32(params[:quant_search_params]["min_log2_matched_ratio"]),
                                    Tuple([Int64(x) for x in params[:quant_search_params]["min_topn_of_m"]]),
                                    Int64(params[:quant_search_params]["max_best_rank"]),
                                    params[:quant_search_params]["n_frag_isotopes"],
                                    UInt8(params[:quant_search_params]["max_frag_rank"]),
                                    kwargs[:rt_index], 
                                    kwargs[:irt_err],
                                    Set(2),
                                )
            end
        end
        psms = fetch.(tasks)
        return psms
    end
    ms_table_path_to_psms_path = Dict{String, String}()
    quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
        try
        parsed_fname = file_path_to_parsed_name[MS_TABLE_PATH]
        rt_df = DataFrame(Arrow.Table(rt_index_paths[parsed_fname]))
        rt_index = buildRtIndex(rt_df,
                                bin_rt_size = bin_rt_size)
        rt_irt = RT_iRT[parsed_fname]
        MS_TABLE = Arrow.Table(MS_TABLE_PATH);
            precursors = spec_lib["precursors"]
            psms = vcat(quantSearch(
                MS_TABLE, 
                params_;
                precursors = spec_lib["precursors"],
                fragment_lookup_table = spec_lib["f_det"],
                rt_index = rt_index,
                ms_file_idx = UInt32(ms_file_idx), 
                rt_to_irt_spline = rt_irt,
                mass_err_model = frag_err_dist_dict[ms_file_idx],
                irt_err = irt_errs[parsed_fname],#irt_errs[ms_file_idx]/3,
                ion_matches = ionMatches,
                ion_misses = ionMisses,
                id_to_col = IDtoCOL,
                ion_templates = ionTemplates,
                iso_splines = iso_splines,
                chromatograms = chromatograms,
                scored_psms = complex_scored_psms,
                unscored_psms = complex_unscored_psms,
                spectral_scores = complex_spectral_scores,
                precursor_weights = precursor_weights,
                quad_transmission_model = quad_model_dict[ms_file_idx],
                )...);
            addSecondSearchColumns!(psms, 
                                            MS_TABLE[:retentionTime],
                                            spec_lib["precursors"][:prec_charge], 
                                            spec_lib["precursors"][:is_decoy],
                                            pid_to_cv_fold);
            psms[!,:charge2] = UInt8.(psms[!,:charge].==2);
            getIsotopesCaptured!(psms,  
                                    isotope_trace_type,
                                    quad_model_dict[ms_file_idx],
                                    psms[!,:scan_idx],
                                    spec_lib["precursors"][:prec_charge], 
                                    spec_lib["precursors"][:mz], 
                                    MS_TABLE[:centerMz],
                                    MS_TABLE[:isolationWidthMz]);
            psms[!,:best_scan] = zeros(Bool, size(psms, 1));
            filter!(x->first(x.isotopes_captured)<2, psms);
            initSummaryColumns!(psms);
            for (key, gpsms) in pairs(groupby(psms, getPsmGroupbyCols(isotope_trace_type)))
                getSummaryScores!(
                    gpsms, 
                    gpsms[!,:weight],
                    gpsms[!,:gof],
                    gpsms[!,:matched_ratio],
                    gpsms[!,:fitted_manhattan_distance],
                    gpsms[!,:fitted_spectral_contrast],
                    gpsms[!,:y_count]
                );
            end
            filter!(x->x.best_scan, psms);
            addPostIntegrationFeatures!(
                psms, 
                precursors[:sequence],
                precursors[:structural_mods],
                precursors[:mz],
                precursors[:irt],
                precursors[:prec_charge],
                precursors[:missed_cleavages],
                MS_TABLE[:TIC],
                MS_TABLE[:mz_array],
                ms_file_idx,
                rt_irt,
                precID_to_iRT

            );

            temp_path = joinpath(quant_psms_folder, parsed_fname*".arrow")
            psms[!,:prob], psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))

            Arrow.write(
                temp_path,
                psms,
                )
        catch e
            throw(e)
            @warn "Quant search for $MS_TABLE_PATH failed..."
        end
    end
    return
end