
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




function find_prob_threshold_streaming(file_paths, q_value_threshold)
    total_targets = 0
    total_decoys = 0
    
    # First pass: just count targets and decoys
    for file_path in file_paths
        df = DataFrame(Arrow.Table(file_path))
        total_targets += sum(df.target)
        total_decoys += sum(.!df.target)
    end
    
    targets_above = total_targets
    decoys_above = total_decoys
    
    # Second pass: stream through probabilities in descending order
    for prob in sort(vcat([DataFrame(Arrow.Table(fp)).prob for fp in file_paths]...), rev=true)
        current_q_value = decoys_above / (targets_above + decoys_above)
        
        if current_q_value > q_value_threshold
            return prob
        end
        
        # Update counts (assuming equal likelihood of target/decoy for simplicity)
        # In a real scenario, you'd need to check the actual target value
        prob_target = total_targets / (total_targets + total_decoys)
        targets_above -= prob_target
        decoys_above -= (1 - prob_target)
    end
    
    # If we get here, all q-values are below the threshold
    return minimum(vcat([DataFrame(Arrow.Table(fp)).prob for fp in file_paths]...))
end

function samplePSMsForXgboost(quant_psms_folder, max_psms)

    file_paths = readdir(quant_psms_folder, join=true)

    psms_count = 0

    for file_path in readdir(quant_psms_folder)
        file_path = joinpath(quant_psms_folder, file_path)
        psms_count += length(Arrow.Table(file_path)[:precursor_idx])
    end

    # Initialize an empty DataFrame to store the results
    result_df = DataFrame()

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

    return result_df
end


