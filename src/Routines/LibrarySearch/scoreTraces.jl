features = [ 
    #:max_prob,
    :median_prob,
    :q90_prob,
    :iso_rank,
    :assymetry,
    :fraction_censored,
    :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :base_width_min,
    :peak_area,
    :points_above_FWHM,
    :points_above_FWHM_01,

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
    :isotope_count,
    :total_ions,
    :best_rank,
    :topn,

    :max_score,
    :mean_score,
    :max_city_fitted,
    :mean_city_fitted,
    :city_block,
    :city_block_fitted,
    :entropy_score,
    :max_entropy,
    :scribe,
    :scribe_corrected,
    :scribe_fitted,
    :spectral_contrast,
    :spectral_contrast_corrected,
    :max_matched_ratio,
    :max_scribe_score,
    
    :data_points,
    :err_norm,
    :error,
    :matched_ratio,
    :poisson,

    :weight,
    :log2_intensity_explained,
    :TIC,
    :adjusted_intensity_explained
];

#best_psms[!,features]
best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
best_psms[!,:decoy] = best_psms[!,:target].==false;
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 0.5, 
                        colsample_bynode = 0.5,
                        min_child_weight = 5, 
                        gamma = 1, 
                        subsample = 0.25, 
                        max_depth = 10, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        iter_scheme = [100, 100, 200],
                        print_importance = false);
best_psms = bst[2];
best_psms[!,:prob] =Float32.(best_psms[!,:prob]);
getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
value_counts(df, col) = combine(groupby(df, col), nrow);
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:file_path])

length(unique(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:precursor_idx]))
println("file specific ids ", sum(IDs_PER_FILE[!,:nrow]))
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers
) => :accession_numbers
);
getBestTrace!(best_psms)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0M1_020824_q99min5max25_besttrace_021324_07.jld2"); best_psms)

