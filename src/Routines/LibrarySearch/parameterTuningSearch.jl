#Clear directory of previous QC Plots
function parameterTuningSearch(rt_alignment_folder,
                                mass_err_estimation_folder,
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
    [rm(joinpath(rt_alignment_folder, x)) for x in readdir(rt_alignment_folder)]
    [rm(joinpath(mass_err_estimation_folder, x)) for x in readdir(mass_err_estimation_folder)]

    rt_to_irt_map_dict = Dict{Int64, Any}()
    frag_err_dist_dict = Dict{Int64,MassErrorModel}()
    irt_errs = Dict{Int64, Float64}()
    ms_file_idx_to_remove = Int64[]
    failed_ms_file_idx = Int64[]
    for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        out_fname = String(first(split(splitpath(MS_TABLE_PATH)[end],".")));
        #Randomly sample spectra to search and retain only the 
        #most probable psms as specified in "first_seach_params"
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
                                    scored_psms = scored_psms,
                                    unscored_psms = unscored_psms,
                                    spectral_scores = spectral_scores,
                                    prec_to_score = precs,
                                    mass_err_model = mass_err_model,
                                    sample_rate = params_[:presearch_params]["sample_rate"],
                                    params = params_[:presearch_params],
                                    quad_transmission_func = QuadTransmission(params_[:quad_transmission]["overhang"], params_[:quad_transmission]["smoothness"])
                                                    );

            psms = vcat([result for result in RESULT]...)
            psms[!,:best_psms ] .= false
            if iszero(size(psms, 1))
                #@warn "No psms in pre-search iteration $n for $MS_TABLE_PATH"
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
                grouped_psms = groupby(rtpsms,:precursor_idx)
                for psms in grouped_psms
                    best_idx = argmax(psms.prob)
                    psms[best_idx,:best_psms] = true
                end
                filter!(x->x.best_psms, rtpsms)
                break
            end
            n += 1
        end

        if n >= params_[:presearch_params]["max_presearch_iters"]
            min_samples = params_[:presearch_params]["min_samples"]
            max_iters = params_[:presearch_params]["max_presearch_iters"]
            @warn "Presearch did not find $min_samples precursors at the specified fdr in $max_iters iterations"
            if rtpsms === nothing
                 @warn "Presearch found zero psms in $max_iters iterations... Removing file $MS_TABLE_PATH from the search..."
                 push!(ms_file_idx_to_remove, ms_file_idx)
                 continue
            end
            filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], rtpsms)
            rtpsms[!,:best_psms] .= false
            grouped_psms = groupby(rtpsms,:precursor_idx)
            for psms in grouped_psms
                best_idx = argmax(psms.prob)
                psms[best_idx,:best_psms] = true
            end
            filter!(x->x.best_psms, rtpsms)
        end
        try
        rt_to_irt_map = UniformSpline( 
                                        rtpsms[!,:irt_predicted], 
                                        rtpsms[!,:rt], 
                                        3, #Degree is always three
                                        5 #Spline knots fixed. currently not tunable parameter
                                        );

        rtpsms[!,:irt_observed] = rt_to_irt_map.(rtpsms[!,:rt])
        irt_MAD = mad(rtpsms[!,:irt_observed] .- rtpsms[!,:irt_predicted])
        irt_σ = irt_MAD #Robust estimate of standard deviation
        #Remove outliers and re-fit spline 
        rtpsms = rtpsms[abs.(rtpsms[!,:irt_observed] .- rtpsms[!,:irt_predicted]) .< irt_MAD*10,:]
        rt_to_irt_map = UniformSpline( 
            rtpsms[!,:irt_predicted], 
            rtpsms[!,:rt], 
            3, #Degree is always three
            5 #Spline knots fixed. currently not tunable parameter
            );

        rtpsms[!,:irt_observed] = rt_to_irt_map.(rtpsms[!,:rt])
        irt_MAD = mad(rtpsms[!,:irt_observed] .- rtpsms[!,:irt_predicted])
        irt_σ = irt_MAD #Robust estimate of standard deviation

        irt_errs[ms_file_idx] = params_[:irt_mapping_params]["n_sigma_tol"]*irt_σ
        rt_to_irt_map_dict[ms_file_idx] = rt_to_irt_map

        plotRTAlign(rtpsms[:,:rt], 
                    rtpsms[:,:irt_predicted], 
                    rt_to_irt_map, 
                    out_fdir = rt_alignment_folder,
                    out_fname = out_fname
                    );
                    
        matched_fragments = vcat(massErrorSearch(
            MS_TABLE,
            rtpsms[!,:scan_idx],
            rtpsms[!,:precursor_idx],
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
        #1) rt to irt curve and 
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
        catch
            psms_count = size(rtpsms, 1)
            @warn "Presearch parameter estimation failed for $MS_TABLE_PATH with psms count: $psms_count, will substitute with 
            parameters from nearest successful presearch"
            push!(failed_ms_file_idx, ms_file_idx)
        end
    end
    #Merge Quality Control PDFs 
    merge_pdfs([joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) if endswith(x, ".pdf")], 
                    joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), cleanup=true)
    merge_pdfs([joinpath(mass_err_estimation_folder, x) for x in readdir(mass_err_estimation_folder) if endswith(x, ".pdf")], 
        joinpath(mass_err_estimation_folder, "mass_err_estimation_folder.pdf"), cleanup=true)
    return rt_to_irt_map_dict, frag_err_dist_dict, irt_errs, ms_file_idx_to_remove, failed_ms_file_idx
end