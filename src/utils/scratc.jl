
best_psms = samplePSMsForXgboost(quant_psms_folder, 500000)

precs_to_score =   unique(best_psms[!,:precursor_idx])
prec_to_best_score = Dictionary(
    precs_to_score,
    zeros(Float32, length(precs_to_score))
)
prec_to_best_score_old = Dictionary(
    precs_to_score,
    zeros(Float32, length(precs_to_score))
)

features = [ 
    :max_prob,
    #
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
    :best_rank_iso,
    :topn,
    :topn_iso,
    :gof,
    :max_fitted_manhattan_distance,
    :max_fitted_spectral_contrast,
    :max_matched_residual,
    :max_unmatched_residual,
    :max_gof,
    :entropy_score,
    :max_entropy,
    :fitted_spectral_contrast,
    :spectral_contrast,
    :max_matched_ratio,
    :err_norm,
    :error,
    :matched_ratio,
    :poisson,
    :weight,
    :log2_intensity_explained,
    :tic,
    :adjusted_intensity_explained
];


bst = xgboost((best_psms[!,features], best_psms[!,:target]), 
num_round=20, 
#monotone_constraints = monotone_constraints,
colsample_bytree = 0.5, 
colsample_bynode = 0.5,
gamma = 1, 
max_depth=10, 
eta = 0.05, 
min_child_weight = 5, 
subsample = 0.25, 
objective="binary:logistic",
seed = rand(UInt32),
#max_bin = 128,
watchlist=(;)
)

file_paths = readdir(quant_psms_folder, join=true)

function getBestScorePerPrec!(
    prec_to_best_score_old::Dictionary{UInt32, Float32},
    prec_to_best_score_new::Dictionary{UInt32, Float32},
    file_paths::Vector{String},
    bst::Booster,
    features::Vector{Symbol}
)
    for file_path in file_paths
        arrow_table = Arrow.Table(file_path) 
        psms = DataFrame(arrow_table)#[!,features]
        for (i, prec_idx) in enumerate(arrow_table[:precursor_idx])
            psms[i,:max_prob] = prec_to_best_score_old[prec_idx]
        end
        #Predict probabilites 
        probs = XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        for (i, prec_idx) in enumerate(arrow_table[:precursor_idx])
            if haskey(prec_to_best_score_new, prec_idx)
                if prec_to_best_score_new[prec_idx] < probs[i]
                    prec_to_best_score_new[prec_idx] = probs[i]
                end
            end
        end
    end
end

prec_to_best_score = Dictionary(
    precs_to_score,
    zeros(Float32, length(precs_to_score))
)
prec_to_best_score_old = Dictionary(
    precs_to_score,
    zeros(Float32, length(precs_to_score))
)
file_paths = readdir(quant_psms_folder, join=true)

@time getBestScorePerPrec!(
    prec_to_best_score_old,
    prec_to_best_score,
    file_paths,
    bst,
    features
)



function getBestScorePerPrec!(
    prec_to_best_score::Dictionary{UInt32, Float32},
    psms::DataFrame,
    bst::Booster,
    features::Vector{Symbol}
)
        psms
        #Predict probabilites 
        psms[!,:prob] = XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        for (i, prec_idx) in enumerate(arrow_table[:precursor_idx])
            if haskey(prec_to_best_score, prec_idx)
                if prec_to_best_score[prec_idx] < probs[i]
                    prec_to_best_score[prec_idx] = probs[i]
                end
            end
        end
end


@time getBestScorePerPrec!(
    prec_to_best_score,
    file_paths,
    bst,
    features
)


for file_path in file_paths
    # Read the Arrow table
    arrow_table = Arrow.Table(file_path)
    
    # Get the number of rows
    num_rows = length(arrow_table[1])
    
    # Calculate the number of rows to sample (1/N'th of the total)
    sample_size = min(ceil(Int, (num_rows/psms_count)*max_psms), num_rows) #ceil(Int, num_rows / N)
    
    # Generate sorted random indices for sampling
    sampled_indices = sort!(sample(1:num_rows, sample_size, replace=false))
    
    # Sample the rows and convert to DataFrame
    sampled_df = DataFrame(arrow_table)[sampled_indices, :]
    
    # Append to the result DataFrame
    append!(result_df, sampled_df)
end

