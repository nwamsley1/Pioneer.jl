function scoreTraces!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors
)
    features = [ 
        :max_prob,
        :mean_prob,
        :min_prob,
        #:max_pg_score,
        #:pg_count,
        #:pepgroup_count,
        :missed_cleavage,
        :Mox,
        :prec_mz,
        :sequence_length,
        :charge,
        :irt_pred,
        :irt_error,
        #:irt_obs,
        #:RT,
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
        :best_rank_iso,
        :topn,
        :topn_iso,
        :gof,
        :max_fitted_manhattan_distance,
        :max_fitted_spectral_contrast,
        :max_matched_residual,
        :max_unmatched_residual,
        :max_gof,
        #:entropy_score,
        #:max_entropy,
        :fitted_spectral_contrast,
        :spectral_contrast,
        :max_matched_ratio,
        :err_norm,
        #:error,
        :matched_ratio,
        :poisson,
        :weight,
        :log2_intensity_explained,
        :tic,
        #:adjusted_intensity_explained
    ];
    best_psms[!,:accession_numbers] = [precursors[:accession_numbers][pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
    best_psms[!,:decoy] = best_psms[!,:target].==false;
    xgboost_time = @timed bst = rankPSMs!(
                            best_psms, 
                            file_paths,
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
                            print_importance = false);



    #best_psms = bst[2];
    #best_psms[!,:prob] =Float32.(best_psms[!,:prob])
    ##Calculate q-values
    #getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
    return bst;#best_psms
end