function scoreTraces!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors
)
    if size(best_psms, 1) > 1000000
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    features = [ 
        :max_prob,
        :mean_prob,
        :min_prob,
        :missed_cleavage,
        :Mox,
        :prec_mz,
        :sequence_length,
        :charge,
        :irt_pred,
        :irt_error,
        :irt_diff,
        :max_y_ions,
        :y_ions_sum,
        :longest_y,
        :y_count,
        :b_count,
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
        #:matched_ratio,
        :poisson,
        :weight,
        :log2_intensity_explained,
        :tic,
    ];
    best_psms[!,:accession_numbers] = [precursors[:accession_numbers][pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
    best_psms[!,:decoy] = best_psms[!,:target].==false;
    models = rankPSMs!(
                            best_psms, 
                            file_paths,
                            features,
                            colsample_bytree = 0.5, 
                            colsample_bynode = 0.5,
                            min_child_weight = 5, 
                            gamma = 1,
                            subsample = 0.25, 
                            max_depth = 10,
                            eta = 0.05, 
                            iter_scheme = [100, 100, 200],
                            print_importance = false);
    return models;#best_psms
    else
        println("plan B")
        println("size(best_psms) ", size(best_psms))
        file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
        features = [ 
            :missed_cleavage,
            :Mox,
            #:prec_mz,
            :sequence_length,
            :charge,
            :irt_error,
            :irt_diff,
            :y_count,
            #:total_ions,
            #:topn,
            #:gof,
            :max_fitted_manhattan_distance,
            #:max_fitted_spectral_contrast,
            :max_matched_residual,
            :max_unmatched_residual,
            :max_gof,
            #:max_matched_ratio,
            :err_norm,
            :weight,
            :log2_intensity_explained,
            #:tic,
        ];
        best_psms[!,:accession_numbers] = [precursors[:accession_numbers][pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
        best_psms[!,:decoy] = best_psms[!,:target].==false;
        models = rankPSMs!(
                                best_psms, 
                                file_paths,
                                features,
                                colsample_bytree = 1.0, 
                                colsample_bynode = 1.0,
                                min_child_weight = 100, 
                                gamma = 0,
                                subsample = 1.0, 
                                max_depth = 3,
                                eta = 0.01, 
                                iter_scheme = [200],
                                print_importance = false);
        return models;#best_psms
    end
end