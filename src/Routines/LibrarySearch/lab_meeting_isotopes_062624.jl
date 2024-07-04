smoothness = Float32[1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 5.0, 5.0, 5.0, 5.0, 5.0]
width = Float32[2.0/2, 2.1/2, 2.2/2, 2.25/2, 2.5/2, 2.0/2, 2.1/2, 2.2/2, 2.25/2, 2.5/2]
ids = []


plot(
    LinRange(-2.0, 2.0, 100),
    [QuadTransmission(width[end], smoothness[1])(0.0, x) for x in LinRange(-2.0, 2.0, 100)]
)
plot!(
    LinRange(-2.0, 2.0, 100),
    [QuadTransmission(width[end], smoothness[end])(0.0, x) for x in LinRange(-2.0, 2.0, 100)]
)
plot!(
    LinRange(-2.0, 4.0, 100),
    [QuadTransmission(width[end], smoothness[end])(2.0, x) for x in LinRange(-2.0, 4.0, 100)]
)

plot!(
    LinRange(-2.0, 4.0, 100),
    [QuadTransmission(1.0, 1000.0)(0.0, x) for x in LinRange(-2.0, 4.0, 100)]
)

vline!([-width[end], width[end]])


vline!([-width[end], width[end]])


#for idx in range(1, length(smoothness))
quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    params_[:deconvolution_params]["huber_delta"] = median(
        [quantile(x, 0.25) for x in MS_TABLE[:intensities]])*params_[:deconvolution_params]["huber_delta_prop"] 
       

        #params_[:deconvolution_params]["huber_delta"] = 100.0f0
        params_[:deconvolution_params]["lambda"] = 0.0f0
        params_[:deconvolution_params]["accuracy_bisection"] = 10.0
        params_[:deconvolution_params]["accuracy_newton"] = 10.0
        params_[:quant_search_params]["n_frag_isotopes"] = 2

        psms = vcat(secondSearch(
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
            quad_transmission_func = QuadTransmission(width[end], smoothness[end])
            )...);
        addSecondSearchColumns!(psms, 
                                        MS_TABLE, 
                                        prosit_lib["precursors"][:mz],
                                        prosit_lib["precursors"][:prec_charge], 
                                        prosit_lib["precursors"][:is_decoy],
                                        pid_to_cv_fold);
        psms[!,:charge2] = UInt8.(psms[!,:charge].==2)
        filter!(x->isnan(x.entropy_score)==false, psms)
        getIsotopesCaptured!(psms, precursors[:prec_charge],precursors[:mz], MS_TABLE)
        psms[!,:best_scan] .= false
        #filter!(x->first(x.isotopes_captured)<2, psms)
        filter!(x->!isinf(x.gof), psms)
        
        filter!(x->x.topn>1, psms)
        filter!(x->x.y_count>1, psms)
        #filter!(x->!isinf(x.scribe_fitted), psms)
        
        
        initSummaryColumns!(psms)
        for (key, gpsms) in ProgressBar(pairs(groupby(psms, [:precursor_idx,:isotopes_captured])))
            
            getSummaryScores!(
                gpsms, 
                gpsms[!,:weight],
                gpsms[!,:scribe],
                gpsms[!,:matched_ratio],
                gpsms[!,:entropy_score],
                gpsms[!,:city_block_fitted],
                gpsms[!,:y_count]
            )
        end
        
        filter!(x->x.best_scan, psms)
        addPostIntegrationFeatures!(
            psms, 
            MS_TABLE,
            precursors[:sequence],
            precursors[:structural_mods],
            precursors[:mz],
            precursors[:irt],
            precursors[:prec_charge],
            precursors[:missed_cleavages],
            ms_file_idx,
            file_id_to_parsed_name,
            RT_iRT,
            precID_to_iRT

        )
        psms[!,:file_name].=file_id_to_parsed_name[ms_file_idx]
        BPSMS[ms_file_idx] = psms;
end

best_psms = vcat(values(BPSMS)...)
filter!(x->!isinf(x.scribe_fitted), best_psms)
#=
tbins = LinRange(-12, 12, 100)
fname = :spectral_contrast_corrected
histogram(best_psms[best_psms[!,:target],fname], bins = tbins, normalize = :pdf, alpha = 0.5)
histogram!(best_psms[best_psms[!,:decoy],fname], bins = tbins, normalize = :pdf, alpha = 0.5)

tbins = LinRange(-12, 12, 100)
fname = :city_block
histogram(best_psms[best_psms[!,:target],fname], bins = tbins, normalize = :pdf, alpha = 0.5)
histogram!(best_psms[best_psms[!,:decoy],fname], bins = tbins, normalize = :pdf, alpha = 0.5)

best_psms[!,:spec_diff] = best_psms[!,:scribe_corrected] .- best_psms[!,:city_block]
tbins = LinRange(0, 1, 100)
fname = :scribe_corrected
histogram(best_psms[best_psms[!,:target],fname], bins = tbins, normalize = :pdf, alpha = 0.5)
histogram!(best_psms[best_psms[!,:decoy],fname], bins = tbins, normalize = :pdf, alpha = 0.5)


best_psms[!,:spec_diff] = best_psms[!,:scribe_corrected] .- best_psms[!,:city_block]
tbins = LinRange(0, 1, 100)
fname = :spectral_contrast
histogram(best_psms[best_psms[!,:target],fname], bins = tbins, normalize = :pdf, alpha = 0.5)
histogram!(best_psms[best_psms[!,:decoy],fname], bins = tbins, normalize = :pdf, alpha = 0.5)



=#
features = [ 
    #:max_prob,
    #:median_prob,
    #:q90_prob,
    #:pg_count,
    #:pepgroup_count,
    :missed_cleavage,
    :Mox,
    :prec_mz,
    :sequence_length,
    :charge,
    :irt_pred,
    :irt_error,
    :irt_obs,
    :RT,
    :irt_diff,
    :max_y_ions,
    :y_ions_sum,
    :longest_y,
    :y_count,
    :b_count,
    :p_count,
    :isotope_count,
    :total_ions,
    :best_rank,
    :topn,
    #:max_score,
    #:mean_score,
    :gof,
    :max_city_fitted,
    :mean_city_fitted,
    :city_block,
    :city_block_fitted,
    #:entropy_score,
    :max_entropy,
    #:scribe,
    #:scribe_fitted,
    :spectral_contrast,
    :spectral_contrast_corrected,
    :scribe_corrected,
    :max_matched_ratio,
    #:max_scribe_score,
    :err_norm,
    :error,
    :matched_ratio,
    :poisson,
    :weight,
    :log2_intensity_explained,
    :tic,
    :adjusted_intensity_explained
];
filter!(x->!isinf(x.max_scribe_score), best_psms)

best_psms[!,:accession_numbers] = [precursors[:accession_numbers][pid] for pid in best_psms[!,:precursor_idx]]
best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
best_psms[!,:decoy] = best_psms[!,:target].==false;
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 0.5, 
                        colsample_bynode = 0.5,
                        min_child_weight = 5, 
                        gamma = 1, 
                        #gamma = 0,
                        subsample = 0.25, 
                        #subsample = 1.0,
                        max_depth = 10,
                        #max_depth = 4, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        iter_scheme = [100, 100, 200],
                        #iter_scheme = [200],
                        print_importance = true);



best_psms = bst[2];
best_psms[!,:prob] =Float32.(best_psms[!,:prob])

#Calculate q-values
getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
#Get best isotope trace for each precursor
#best_psms_passing = best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]),:]#List of top scoring precursors to requantify. 
const precursors_passing = Set(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]),:precursor_idx])
#getBestTrace!(best_psms, 0.01, :weight)
##Re-calculate q-values after removing inferior isotopic trace
#getQvalues!(best_psms[!,:prob].*(best_psms[!,:best_trace]), best_psms[:,:target], best_psms[!,:q_value]);
#const traces_passing = Set([(x[1], x[2]) for x in eachrow(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]).&(best_psms[!,:best_trace]),[:precursor_idx,:isotopes_captured]])])
#Get protein names 
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][:accession_numbers][psm[:precursor_idx]]
) => :accession_numbers
);

#getBestTrace!(best_psms)
value_counts(df, col) = combine(groupby(df, col), nrow);
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:], [:file_name])

end
