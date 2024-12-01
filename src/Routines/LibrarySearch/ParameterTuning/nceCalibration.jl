#Clear directory of previous QC Plots
function nceTuningSearch(
                                rt_to_irt_map_dict,
                                frag_err_dist_dict,
                                irt_errs,
                                MS_TABLE_PATHS,
                                params_,
                                spec_lib,
                                nce_grid,
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
    precursors = spec_lib["precursors"]
    nce_model_dict = Dict{Int64, Any}()
    for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(enumerate(MS_TABLE_PATHS))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        RESULT =  NceScanningSearch(
                                MS_TABLE,
                                params_,
                                nce_grid;
                                frag_index = spec_lib["f_index"],
                                precursors = precursors,
                                fragment_lookup_table = spec_lib["f_det"],
                                rt_to_irt_spline =  rt_to_irt_map_dict[ms_file_idx],
                                ms_file_idx = UInt32(1),
                                irt_tol =  irt_errs[ms_file_idx],
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
                                mass_err_model = frag_err_dist_dict[ms_file_idx],
                                sample_rate = params_[:presearch_params]["sample_rate"],
                                params = params_[:presearch_params],
                                isotope_err_bounds = params_[:isotope_err_bounds],
                                quad_transmission_model = SquareQuadModel(0.0f0));

        psms = RESULT#vcat([result for result in RESULT]...)
        psms[!,:best_psms ] .= false

        addPreSearchColumns!(psms, 
                                    MS_TABLE, 
                                    spec_lib["precursors"][:is_decoy],
                                    spec_lib["precursors"][:irt],
                                    spec_lib["precursors"][:prec_charge],
                                    MS_TABLE[:retentionTime],
                                    MS_TABLE[:TIC]
        )


        spsms = DataFrame(map(pairs(groupby(psms,[:precursor_idx,:scan_idx]))) do (key, psms)
            max_arg = argmax(psms[!,:scribe])
            return psms[max_arg,:]
        end)
        scorePresearch!(spsms)
        getQvalues!(spsms[!,:prob], spsms[!,:target], spsms[!,:q_value])
        filter!(x->(x.target)&(x.q_value<=0.01), spsms)
        passing_precs = Set(spsms[!,:precursor_idx])
        filter!(x->x.precursor_idx∈passing_precs, psms)
        psms[!,:best_psms] .= false
        psms[!,:prec_mz] = [precursors[:mz][pid] for pid in psms[!,:precursor_idx]]
        for (key, psms) in pairs(groupby(psms, :precursor_idx))
            psms[argmax(psms[!,:scribe]),:best_psms] = true
        end
        filter!(x->x.best_psms, psms)
        nce_model_dict[ms_file_idx] = fit_nce_model(PiecewiseNceModel(0.0f0),
        psms[!,:prec_mz], psms[!,:nce], psms[!,:charge], NCE_MODEL_BREAKPOINT)
        #return psms,  fit_nce_model(PiecewiseNceModel(0.0f0),
        #psms[!,:prec_mz], psms[!,:nce], psms[!,:charge], NCE_MODEL_BREAKPOINT)
    end
    return nce_model_dict
end