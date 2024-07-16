PSMs_Dict = Dictionary{String, DataFrame}()

main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
    params_[:first_search_params]["n_frag_isotopes"] = 1
    #params_[:first_search_params]["min_spectral_contrast"] = 0.5#acos(0.5)
    #params_[:first_search_params]["min_log2_matched_ratio"] = 0.0
    PSMs = vcat(LibrarySearch(
        MS_TABLE,
        params_;
        frag_index = prosit_lib["f_index"],
        precursors = prosit_lib["precursors"],
        fragment_lookup_table = library_fragment_lookup_table,
        rt_to_irt_spline =  RT_to_iRT_map_dict[ms_file_idx],
        ms_file_idx = UInt32(ms_file_idx),
        irt_tol = irt_errs[ms_file_idx],
        ion_matches = ionMatches,
        ion_misses = ionMisses,
        fmatches = all_fmatches,
        id_to_col = IDtoCOL,
        ion_templates = ionTemplates,
        iso_splines = iso_splines,
        scored_psms = scored_PSMs,
        unscored_psms = unscored_PSMs,
        spectral_scores = spectral_scores,
        prec_to_score = precs,
        mass_err_model = frag_err_dist_dict[ms_file_idx],
        sample_rate = Inf,
        params = params_[:first_search_params],
        quad_transmission_func = QuadTransmission(1.25f0, 5.0f0)
                        )...)

    addMainSearchColumns!(PSMs, MS_TABLE, 
                        prosit_lib["precursors"][:structural_mods],
                        prosit_lib["precursors"][:missed_cleavages],
                        prosit_lib["precursors"][:is_decoy],
                        prosit_lib["precursors"][:irt],
                        prosit_lib["precursors"][:prec_charge]);
    
    #Observed iRT estimates based on pre-search
    PSMs[!,:iRT_observed] = RT_to_iRT_map_dict[ms_file_idx].(PSMs[!,:RT])
    PSMs[!,:iRT_error] = Float16.(abs.(PSMs[!,:iRT_observed] .- PSMs[!,:iRT_predicted]))
    psms = copy(PSMs)

    PSMs[!,:i_count] = Float16.(PSMs[!,:i_count]./(PSMs[!,:y_count].+PSMs[!,:b_count].+PSMs[!,:p_count]))
    filter!(x->isinf(x.i_count)==false, PSMs)
    filter!(x->isnan(x.i_count)==false, PSMs)
    PSMs[!,:charge2] = UInt8.(PSMs[!,:charge].==2)

    column_names = [
                    :spectral_contrast,
                    :city_block,
                    :entropy_score,
                    :scribe,
                    #:combined,
                    :charge2,
                    :poisson,
                    :iRT_error,
                    :missed_cleavage,
                    :Mox,
                    :charge,
                    :TIC,
                    :y_count,
                    #:i_count,
                    :err_norm,
                    :spectrum_peak_count,
                    :intercept]
    #PSMs[!,:combined] .= zero(Float16)
    #for i in range(1, size(PSMs, 1))
    #    PSMs[i,:combined] = PSMs[i,:scribe] + PSMs[i,:city_block] + log(PSMs[i,:spectral_contrast]) - PSMs[i,:entropy_score]
    #end
    #PSMs[!,:combined] = exp.(PSMs[!,:combined])
    #PSMs[!,:scribe2] = exp.(Float32.(PSMs[!,:scribe]))
    scoreMainSearchPSMs!(PSMs,
                                column_names,
                                n_train_rounds = params_[:first_search_params]["n_train_rounds_probit"],
                                max_iter_per_round = params_[:first_search_params]["max_iter_probit"],
                                max_q_value = params_[:first_search_params]["max_q_value_probit_rescore"]);

    getProbs!(PSMs);
    getBestPSMs!(PSMs,
                    prosit_lib["precursors"][:mz],
                    max_q_value = Float64(params_[:first_search_params]["max_q_value_filter"]),
                    max_psms = Int64(params_[:first_search_params]["max_precursors_passing"])
                )

    insert!(PSMs_Dict, 
        file_id_to_parsed_name[ms_file_idx], 
        PSMs
    );
end
