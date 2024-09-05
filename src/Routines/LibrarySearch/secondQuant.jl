function secondQuantSearch!( 
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
    precursors,
    accession_number_to_id,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_psms,
    complex_unscored_psms,
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
    parsed_fname = file_path_to_parsed_name[MS_TABLE_PATH]
    rt_df = DataFrame(Arrow.Table(rt_index_paths[parsed_fname]))
    rt_index = buildRtIndex(rt_df,
                            bin_rt_size = bin_rt_size)
    #Map raw file name to psms table 
    sub_bpsms_path = joinpath(passing_psms_folder, parsed_fname.*".arrow")
    #Get psms table for this raw file. 
    sub_bpsms = DataFrame(
                    Tables.columntable(#Need a modifiable deep copy
                            Arrow.Table(
                                sub_bpsms_path
                                )))
    filter!(x->x.target, sub_bpsms) #remove decoys IMPOrtANT
    sub_bpsms[!,:peak_area] = zeros(Float32, size(sub_bpsms, 1))
    sub_bpsms[!,:new_best_scan] = zeros(UInt32, size(sub_bpsms, 1))

        chroms = vcat(secondQuant(
            MS_TABLE, 
            params_;
            precursors = spec_lib["precursors"],
            fragment_lookup_table = spec_lib["f_det"],
            rt_index = rt_index,
            ms_file_idx = UInt32(ms_file_idx), 
            rt_to_irt_spline = rt_irt[parsed_fname],
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
            traces_passing = Set(sub_bpsms[!,:precursor_idx]),
            quad_transmission_func = QuadTransmission(params_[:quad_transmission]["overhang"], params_[:quad_transmission]["smoothness"])
            )...);


        #Format Chromatograms 
        getIsotopesCaptured!(chroms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        filter!(x->first(x.isotopes_captured)<2, chroms)
        gchroms = groupby(chroms,[:precursor_idx,:isotopes_captured])
        integratePrecursors(
                            gchroms,
                            sub_bpsms[!,:precursor_idx],
                            sub_bpsms[!,:isotopes_captured],
                            sub_bpsms[!,:scan_idx],
                            sub_bpsms[!,:peak_area],
                            sub_bpsms[!,:new_best_scan],
                            λ=Float32(params_[:quant_search_params]["WH_smoothing_strength"]),
                            n_pad = params_[:quant_search_params]["n_pad"],
                            max_apex_offset = params_[:quant_search_params]["max_apex_offset"]
                            )

        temp_path = joinpath(second_quant_folder, file_path_to_parsed_name[MS_TABLE_PATH]*".arrow")
        #psms[!,:prob], psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
        filter!(x->!(isnan(x.peak_area)), sub_bpsms)
        sub_bpsms[!,:accession_numbers] = [precursors[:accession_numbers][prec_idx] for prec_idx in sub_bpsms[!,:precursor_idx]]
        sub_bpsms[!,:protein_idx] = [accession_number_to_id[accession_numbers] for accession_numbers in sub_bpsms[!,:accession_numbers]]
        sub_bpsms[!,:ms_file_idx] =  UInt32.(sub_bpsms[!,:ms_file_idx])
        sort!(sub_bpsms, [:protein_idx,:precursor_idx,:ms_file_idx])
        sub_bpsms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in sub_bpsms[!,:precursor_idx]]
        sub_bpsms[!,:peak_area] =  allowmissing(sub_bpsms[!,:peak_area])
        sub_bpsms[!,:peak_area_normalized] =  allowmissing(zeros(Float32, size(sub_bpsms, 1)))
        sub_bpsms[!,:structural_mods] = allowmissing([precursors[:structural_mods][pid] for pid in sub_bpsms[!,:precursor_idx]])
        sub_bpsms[!,:isotopic_mods] = allowmissing([precursors[:isotopic_mods][pid] for pid in  sub_bpsms[!,:precursor_idx]])
        sub_bpsms[!,:charge] = allowmissing([precursors[:prec_charge][pid] for pid in  sub_bpsms[!,:precursor_idx]])
        sub_bpsms[!,:sequence] = allowmissing([precursors[:sequence][pid] for pid in  sub_bpsms[!,:precursor_idx]])
        sub_bpsms[!,:file_name] .= file_path_to_parsed_name[MS_TABLE_PATH]
        Arrow.write(
            temp_path,
            sub_bpsms,
            )
end
end