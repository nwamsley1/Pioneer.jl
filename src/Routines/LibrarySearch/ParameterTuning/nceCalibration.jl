#Clear directory of previous QC Plots
function nceTuningSearch(
                                MS_TABLE_PATH,
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

    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    n = 0
    rtpsms = nothing
    mass_err_model = MassErrorModel(
        0.0f0,
        (
            Float32(params_[:presearch_params]["frag_tol_ppm"]), 
            Float32(params_[:presearch_params]["frag_tol_ppm"])
        )
    )  
    while n <= params_[:presearch_params]["max_presearch_iters"]
        RESULT =  NceScanningSearch(
                                MS_TABLE,
                                params_,
                                nce_grid;
                                frag_index = spec_lib["presearch_f_index"],
                                precursors = spec_lib["precursors"],
                                fragment_lookup_table = spec_lib["f_det"],
                                rt_to_irt_spline =  x->x,
                                ms_file_idx = UInt32(1),
                                irt_tol = Inf,
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
                                mass_err_model = mass_err_model,
                                sample_rate = params_[:presearch_params]["sample_rate"],
                                params = params_[:presearch_params],
                                isotope_err_bounds = params_[:isotope_err_bounds],
                                quad_transmission_model = SquareQuadModel(0.0f0));

        psms = RESULT#vcat([result for result in RESULT]...)
        println("final size(psms) ", size(psms))
        psms[!,:best_psms ] .= false
        if iszero(size(psms, 1))
            n += 1
            continue
        end
        addPreSearchColumns!(psms, 
                                    MS_TABLE, 
                                    spec_lib["precursors"][:is_decoy],
                                    spec_lib["precursors"][:irt],
                                    spec_lib["precursors"][:prec_charge],
                                    MS_TABLE[:retentionTime],
                                    MS_TABLE[:TIC]
                                )
        if rtpsms === nothing
            rtpsms = psms
        else
            rtpsms = vcat(rtpsms, psms)
        end            
        try
            #Only allow targets?
            scorePresearch!(rtpsms)
            getQvalues!(rtpsms[!,:prob], rtpsms[!,:target], rtpsms[!,:q_value])
        catch
            n += 1
            continue
        end
            
        if sum(rtpsms[!,:q_value].<=params_[:presearch_params]["max_qval"]) >= params_[:presearch_params]["min_samples"]
            filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], rtpsms)
            rtpsms[!,:best_psms] .= false
            grouped_psms = groupby(rtpsms,[:nce,:precursor_idx])
            for psms in grouped_psms
                best_idx = argmax(psms.prob)
                psms[best_idx,:best_psms] = true
            end
            filter!(x->x.best_psms, rtpsms)
            break
        end
        n += 1
    end

    rtpsms
end