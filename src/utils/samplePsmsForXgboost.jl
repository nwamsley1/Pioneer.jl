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

