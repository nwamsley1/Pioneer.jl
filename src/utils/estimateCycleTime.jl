
psms_count = 0

for file_path in readdir(quant_psms_folder)
    file_path = joinpath(quant_psms_folder, file_path)
    psms_count += length(Arrow.Table(file_path)[:precursor_idx])
end

num_rows = length(arrow_table[:precursor_idx])
n = 250000
sampled_indices = sort!(StatsBase.sample(1:num_rows, n, replace=false))
DataFrame(arrow_table)[sampled_indices,:]


function sample_and_append_arrow_files(file_paths, N)
    # Initialize an empty DataFrame to store the results
    result_df = DataFrame()

    for file_path in file_paths
        # Read the Arrow table
        arrow_table = Arrow.Table(file_path)
        
        # Get the number of rows
        num_rows = length(arrow_table[1])
        
        # Calculate the number of rows to sample (1/N'th of the total)
        sample_size = ceil(Int, num_rows / N)
        
        # Generate sorted random indices for sampling
        sampled_indices = sort!(sample(1:num_rows, sample_size, replace=false))
        
        # Sample the rows and convert to DataFrame
        sampled_df = DataFrame(arrow_table)[sampled_indices, :]
        
        # Append to the result DataFrame
        append!(result_df, sampled_df)
    end

    return result_df
end
sampled_df = sample_and_append_arrow_files(readdir(quant_psms_folder, join=true), 6)

features = [ 
    :max_prob,
    :median_prob,
    :q90_prob,
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

