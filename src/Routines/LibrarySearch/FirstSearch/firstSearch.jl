

function firstSearch(
                    temp_folder,
                    rt_to_irt_map_dict,
                    frag_err_dist_dict,
                    irt_errs,
                    quad_model_dict,
                    nce_model_dict,
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
                    scored_psms,
                    unscored_psms,
                    spectral_scores,
                    precs)

peak_fwhms = Dictionary{String, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}}()
psms_paths = Dictionary{String, String}()
lft = spec_lib["f_det"]
main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    try
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
    lft = updateNceModel(lft, nce_model_dict[ms_file_idx])
    psms = vcat(LibrarySearch(
        MS_TABLE,
        params_;
        frag_index = spec_lib["f_index"],
        precursors = spec_lib["precursors"],
        fragment_lookup_table = lft,
        rt_to_irt_spline =  rt_to_irt_map_dict[ms_file_idx],
        ms_file_idx = UInt32(ms_file_idx),
        irt_tol = irt_errs[ms_file_idx],
        ion_matches = ionMatches,
        ion_misses = ionMisses,
        fmatches = all_fmatches,
        id_to_col = IDtoCOL,
        ion_templates = ionTemplates,
        iso_splines = iso_splines,
        scored_psms = scored_psms,
        unscored_psms = unscored_psms,
        spectral_scores = spectral_scores,
        prec_to_score = precs,
        mass_err_model = frag_err_dist_dict[ms_file_idx],#MassErrorModel(0.9f0, (12.0f0, 17.0f0)),#frag_err_dist_dict[ms_file_idx],
        sample_rate = Inf,
        params = params_[:first_search_params],
        isotope_err_bounds = params_[:isotope_err_bounds],
        quad_transmission_model = quad_model_dict[ms_file_idx]
                        )...)
    #println("size(psms) for $MS_TABLE_PATH ", size(psms, 1))
    #println("add colums time")
    #println("size(psms, 2) ", size(psms, 2))
    #println("size(psms, 1) ", size(psms, 1))
    addMainSearchColumns!(psms,
                        rt_to_irt_map_dict[ms_file_idx],
                        spec_lib["precursors"][:structural_mods],
                        spec_lib["precursors"][:missed_cleavages],
                        spec_lib["precursors"][:is_decoy],
                        spec_lib["precursors"][:irt],
                        spec_lib["precursors"][:prec_charge],
                        MS_TABLE[:retentionTime],
                        MS_TABLE[:TIC],
                        MS_TABLE[:mz_array]);
   
    #Think irt sort is for figured fwhm later? Do this after scoring?
    #println("charge time ")
    #Observed irt estimates based on pre-search
    psms[!,:irt_observed] = rt_to_irt_map_dict[ms_file_idx].(psms[!,:rt])
    psms[!,:irt_error] = Float16.(abs.(psms[!,:irt_observed] .- psms[!,:irt_predicted]))
    #psms = copy(psms)
    #psms[!,:i_count] = Float16.(psms[!,:i_count]./(psms[!,:y_count].+psms[!,:b_count].+psms[!,:p_count]))
    #filter!(x->isinf(x.i_count)==false, psms)
    #filter!(x->isnan(x.i_count)==false, psms)
    psms[!,:charge2] = UInt8.(psms[!,:charge].==2)
    #temp_path =  joinpath(temp_folder, file_id_to_parsed_name[ms_file_idx]*"full.arrow")
    ##Arrow.write(
    #    temp_path,
    #    psms
    #    )
    psms[!,:ms_file_idx] .= UInt32(ms_file_idx)
    column_names = [
                    :spectral_contrast,
                    :city_block,
                    :entropy_score,
                    :scribe,
                    :charge2,
                    :poisson,
                    :irt_error,
                    :missed_cleavage,
                    :Mox,
                    :charge,
                    :TIC,
                    :y_count,
                    :err_norm,
                    :spectrum_peak_count,
                    :intercept]
    select!(psms, vcat(column_names, [:ms_file_idx,:score,:precursor_idx,:scan_idx,:q_value,:log2_summed_intensity,:irt,:rt,:irt_predicted,
    :target]))
    try 
        #println("score time")
        scoreMainSearchpsms!(psms,
                                column_names,
                                n_train_rounds = params_[:first_search_params]["n_train_rounds_probit"],
                                max_iter_per_round = params_[:first_search_params]["max_iter_probit"],
                                max_q_value = params_[:first_search_params]["max_q_value_probit_rescore"]);
    catch e
        throw(e)
        psms[!,:score] .= one(Float32)
    end
    #temp_path =  joinpath(temp_folder, file_id_to_parsed_name[ms_file_idx]*"fullscored.arrow")
    #Arrow.write(
    #    temp_path,
    #    psms
    #    )
    #println("probs time ")
    select!(psms, 
    [:ms_file_idx,:score,:precursor_idx,:scan_idx,:q_value,:log2_summed_intensity,:irt,:rt,:irt_predicted])
    getProbs!(psms);
    #println("getbestpsms time ")
    temp_path =  joinpath(temp_folder, file_id_to_parsed_name[ms_file_idx]*"fullscoredprobs.arrow")
    #Arrow.write(
    #    temp_path,
    #    psms
    #    )
    #println("sort time ")
    sort!(psms, :irt);
    getBestPSMs!(psms,
                    spec_lib["precursors"][:mz],
                    max_q_val = Float32(params_[:summarize_first_search_params]["max_q_val_for_irt"]),
                    max_psms = Int64(params_[:first_search_params]["max_precursors_passing"])
                )
    #Arrow.write(
    #    temp_path,
    #    psms
    #    )
    #irt_diffs is the difference between the first and last psm for a precursor
    #below the `max_q_val_for_irt` threshold. Used as a proxy for peak width
    fwhms = skipmissing(psms[!,:fwhm])
    fwhm_points = 0
    for fwhm in fwhms
        if !ismissing(fwhm)
            fwhm_points += 1
        end
    end

    if fwhm_points < params_[:summarize_first_search_params]["min_inference_points"]
        #@warn "Not enough datapoints to infer peak width. n: "*string(fwhm_points)*"\n 
        #using default "*string(params_[:summarize_first_search_params]["default_irt_width"])
        #Integration width is double the fwhm + n times the fwhm standard deviation
        #as estimated from mad (median absolute deviation normalized to a robust 
        #estimate of standard deviation )
        insert!(
            peak_fwhms,
            file_id_to_parsed_name[ms_file_idx], 
            (median_fwhm = median(fwhms), mad_fwhm = mad(fwhms, normalize = true))
        )
    else
        insert!(
            peak_fwhms,
            file_id_to_parsed_name[ms_file_idx], 
            (median_fwhm = median(fwhms), mad_fwhm = mad(fwhms, normalize = true))
        )
    end
    temp_path =  joinpath(temp_folder, file_id_to_parsed_name[ms_file_idx]*".arrow")
    psms[!,:ms_file_idx] .= UInt32(ms_file_idx)
    #println("write time ")
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx,:scan_idx,:precursor_idx,:rt,:irt_predicted,:q_value,:score,:prob,:scan_count])
        )

    insert!(
        psms_paths,
        file_id_to_parsed_name[ms_file_idx],
        temp_path
    )
    catch e
        throw(e)
        @warn "First search failed for $MS_TABLE_PATH"
        continue
    end
end
return peak_fwhms, psms_paths
end

