PSMs_Dict = Dictionary{String, DataFrame}()
main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
    params_[:first_search_params]["n_frag_isotopes"] = 3
    @time PSMs = vcat(LibrarySearch(
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
        sample_rate = Inf,#params_[:presearch_params]["sample_rate"],
        params = params_[:first_search_params]
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
    
    column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,
                    :iRT_error,:missed_cleavage,:Mox,
                    :charge,:TIC,
                    :y_count,:err_norm,:spectrum_peak_count,:intercept]
    if sum(PSMs[!,:p_count])>0
        push!(column_names, :p_count)
    end
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
                sum(PSMs[!,:q_value].<=0.01)
                sum(PSMs[!,:q_value].<=0.1)
    insert!(PSMs_Dict, 
        file_id_to_parsed_name[ms_file_idx], 
        PSMs
    );
end