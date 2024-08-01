#Clear directory of previous QC Plots
function parameterTuningSearch(rt_alignment_folder,
                                mass_err_estimation_folder,
                                MS_TABLE_PATHS,
                                params_,
                                spec_lib,
                                library_fragment_lookup_table,
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
    [rm(joinpath(rt_alignment_folder, x)) for x in readdir(rt_alignment_folder)]
    [rm(joinpath(mass_err_estimation_folder, x)) for x in readdir(mass_err_estimation_folder)]

    RT_to_iRT_map_dict = Dict{Int64, Any}()
    frag_err_dist_dict = Dict{Int64,MassErrorModel}()
    irt_errs = Dict{Int64, Float64}()
    #MS_TABLE_PATH = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_4.arrow"
    #=
    MS_TABLE_PATHS = ["/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/astral_test/20230324_OLEP08_200ng_30min_E10H50Y40_180K_2Th3p5ms_01.arrow"]
    MS_TABLE_PATH = MS_TABLE_PATHS[1]
    ms_file_idx = 1
    =#
    for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        out_fname = String(first(split(splitpath(MS_TABLE_PATH)[end],".")));
        #Randomly sample spectra to search and retain only the 
        #most probable psms as specified in "first_seach_params"
        data_points = 0
        n = 0
        rtPSMs = nothing
        mass_err_model = MassErrorModel(
            0.0f0,
            (
                Float32(params_[:presearch_params]["frag_tol_ppm"]), 
                Float32(params_[:presearch_params]["frag_tol_ppm"])
            )
        )  
        while n <= params_[:presearch_params]["max_presearch_iters"]
            RESULT =  LibrarySearch(
                                    MS_TABLE,
                                    params_;
                                    frag_index = spec_lib["presearch_f_index"],
                                    precursors = spec_lib["precursors"],
                                    fragment_lookup_table = spec_lib["f_det"],
                                    rt_to_irt_spline =  x->x,
                                    ms_file_idx = UInt32(ms_file_idx),
                                    irt_tol = Inf,
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
                                    mass_err_model = mass_err_model,
                                    sample_rate = params_[:presearch_params]["sample_rate"],
                                    params = params_[:presearch_params],
                                    quad_transmission_func = QuadTransmission(1.0f0, 10000.0f0)
                                                    );

            PSMs = vcat([result for result in RESULT]...)
            PSMs[!,:best_psms ] .= false
            addPreSearchColumns!(PSMs, 
                                        MS_TABLE, 
                                        spec_lib["precursors"][:is_decoy],
                                        spec_lib["precursors"][:irt],
                                        spec_lib["precursors"][:prec_charge]
                                    )

            if rtPSMs === nothing
                rtPSMs = PSMs
            else
                rtPSMs = vcat(rtPSMs, PSMs)
            end            

            scorePresearch!(rtPSMs)
            getQvalues!(rtPSMs[!,:prob], rtPSMs[!,:target], rtPSMs[!,:q_value])
            
            #println("TEST ",  sum(rtPSMs[!,:q_value].<=params_[:presearch_params]["max_qval"]))
            if sum(rtPSMs[!,:q_value].<=params_[:presearch_params]["max_qval"]) >= params_[:presearch_params]["min_samples"]
                filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], rtPSMs)
                rtPSMs[!,:best_psms] .= false
                grouped_psms = groupby(rtPSMs,:precursor_idx)
                for psms in grouped_psms
                    best_idx = argmax(psms.prob)
                    psms[best_idx,:best_psms] = true
                end
                filter!(x->x.best_psms, rtPSMs)
                break
            end
            n += 1
        end

        if n >= params_[:presearch_params]["max_presearch_iters"]
            min_samples = params_[:presearch_params]["min_samples"]
            max_iters = params_[:presearch_params]["max_presearch_iters"]
            @warn "Presearch did not find $min_samples precursors at the specified fdr in $max_iters iterations"
            filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], rtPSMs)
                rtPSMs[!,:best_psms] .= false
                grouped_psms = groupby(rtPSMs,:precursor_idx)
                for psms in grouped_psms
                    best_idx = argmax(psms.prob)
                    psms[best_idx,:best_psms] = true
                end
                filter!(x->x.best_psms, rtPSMs)
        end

        RT_to_iRT_map = UniformSpline( 
                                        rtPSMs[!,:iRT_predicted], 
                                        rtPSMs[!,:RT], 
                                        3, #Degree is always three
                                        5 #Spline knots fixed. currently not tunable parameter
                                        );

        rtPSMs[!,:iRT_observed] = RT_to_iRT_map.(rtPSMs[!,:RT])
        irt_MAD = mad(rtPSMs[!,:iRT_observed] .- rtPSMs[!,:iRT_predicted])
        irt_σ = irt_MAD #Robust estimate of standard deviation

        rtPSMs = rtPSMs[abs.(rtPSMs[!,:iRT_observed] .- rtPSMs[!,:iRT_predicted]) .< irt_MAD*10,:]

        RT_to_iRT_map = UniformSpline( 
            rtPSMs[!,:iRT_predicted], 
            rtPSMs[!,:RT], 
            3, #Degree is always three
            5 #Spline knots fixed. currently not tunable parameter
            );

        rtPSMs[!,:iRT_observed] = RT_to_iRT_map.(rtPSMs[!,:RT])
        irt_MAD = mad(rtPSMs[!,:iRT_observed] .- rtPSMs[!,:iRT_predicted])
        irt_σ = irt_MAD #Robust estimate of standard deviation

        irt_errs[ms_file_idx] = params_[:irt_err_sigma]*irt_σ
        RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map

        plotRTAlign(rtPSMs[:,:RT], 
                    rtPSMs[:,:iRT_predicted], 
                    RT_to_iRT_map, 
                    out_fdir = rt_alignment_folder,
                    out_fname = out_fname
                    );
                    
        matched_fragments = vcat(massErrorSearch(
            MS_TABLE,
            rtPSMs[!,:scan_idx],
            rtPSMs[!,:precursor_idx],
            spec_lib["f_det"],
            UInt32(ms_file_idx),
            mass_err_model,
            ionMatches,
            ionMisses,
            all_fmatches,
            ionTemplates
        )...);

        function _getPPM(a::T, b::T) where {T<:AbstractFloat}
            Float32((b - a)/(a/1e6))
        end
        #Get Retention Times and Target/Decoy Status 
        ####################
        #Use best_psms to estimate 
        #1) RT to iRT curve and 
        #2) mass error (ppm) distribution 

        frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments];
        mass_err_model = ModelMassErrs(
            frag_ppm_errs,
            frag_err_quantile = Float32(params_[:presearch_params]["frag_err_quantile"]),
            out_fdir = mass_err_estimation_folder,
            out_fname = out_fname
        )

        #PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
        #File name but remove file type
        frag_err_dist_dict[ms_file_idx] = mass_err_model
    end
    #Merge Quality Control PDFs 
    merge_pdfs([joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) if endswith(x, ".pdf")], 
                    joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), cleanup=true)
    merge_pdfs([joinpath(mass_err_estimation_folder, x) for x in readdir(mass_err_estimation_folder) if endswith(x, ".pdf")], 
        joinpath(mass_err_estimation_folder, "mass_err_estimation_folder.pdf"), cleanup=true)
        return RT_to_iRT_map_dict, frag_err_dist_dict, irt_errs
end
export parameterTuningSearch