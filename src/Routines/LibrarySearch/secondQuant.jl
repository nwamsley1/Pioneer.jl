function secondQuant!( 
    grouped_best_psms,   
    frag_err_dist_dict,
    traces_passing,
    RT_INDICES,
    RT_iRT,
    irt_err,
    chromatograms,
    file_id_to_parsed_name,
    MS_TABLE_PATHS,
    params_,
    precursors,
    prosit_lib,
    library_fragment_lookup_table,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_PSMs,
    complex_unscored_PSMs,
    complex_spectral_scores,
    precursor_weights)
function secondQuant(
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
            return getChromatograms(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                kwargs[:traces_passing],
                                kwargs[:fragment_lookup_table], 
                                kwargs[:ms_file_idx],
                                kwargs[:rt_to_irt_spline],
                                kwargs[:mass_err_model],
                                kwargs[:quad_transmission_func],
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
                                kwargs[:unscored_psms][thread_id],
                                kwargs[:spectral_scores][thread_id],
                                kwargs[:precursor_weights][thread_id],
                                (3, 0),
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

quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    params_[:deconvolution_params]["huber_delta"] = median(
        [quantile(x, 0.25) for x in MS_TABLE[:intensities]])*params_[:deconvolution_params]["huber_delta_prop"] 
        chroms = vcat(secondQuant(
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
            traces_passing = traces_passing,
            quad_transmission_func = QuadTransmission(1.0f0, 1000.0f0)
            )...);


        #Format Chromatograms 
        getIsotopesCaptured!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        filter!(x->first(x.isotopes_captured)<2, chroms)
        gchroms = groupby(chroms,[:precursor_idx,:isotopes_captured])
        sub_bpsms = grouped_best_psms[(file_name = file_id_to_parsed_name[ms_file_idx],)]
        integratePrecursors(
                            gchroms,
                            sub_bpsms[!,:precursor_idx],
                            sub_bpsms[!,:isotopes_captured],
                            sub_bpsms[!,:scan_idx],
                            sub_bpsms[!,:peak_area],
                            sub_bpsms[!,:new_best_scan],
                            λ=Float32(params_[:quant_search_params]["WH_smoothing_strength"])
                            )
        #BPSMS[ms_file_idx] = psms_integrated;
end
end