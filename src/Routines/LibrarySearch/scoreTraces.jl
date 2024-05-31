features = [ 
    #:max_prob,
    :median_prob,
    :q90_prob,
    :assymetry,
    #:fraction_censored,
    :FWHM,
    :FWHM_01,
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
    :p_count,
    #:non_cannonical_count,
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
    #:scribe_corrected,
    :scribe_fitted,
    :spectral_contrast,
    #:spectral_contrast_corrected,
    :max_matched_ratio,
    :max_scribe_score,
    
    #:data_points,
    :error_norm,
    :error,
    :matched_ratio,
    :poisson,

    :weight,
    :log2_intensity_explained,
    :tic,
    :adjusted_intensity_explained
];


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
getBestTrace!(best_psms, 0.01, :trapezoid_area)
#Re-calculate q-values after removing inferior isotopic trace
getQvalues!(best_psms[!,:prob].*(best_psms[!,:best_trace]), best_psms[:,:target], best_psms[!,:q_value]);
const traces_passing = Set([(x[1], x[2]) for x in eachrow(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]).&(best_psms[!,:best_trace]),[:precursor_idx,:isotopes_captured]])])
#Get protein names 
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][:accession_numbers][psm[:precursor_idx]]
) => :accession_numbers
);

#getBestTrace!(best_psms)
value_counts(df, col) = combine(groupby(df, col), nrow);
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01).& (best_psms[:,:decoy].==false),:], [:file_name])
