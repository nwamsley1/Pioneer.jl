function find_prob_threshold(
                            file_paths::Vector{String},
                            q_value_threshold::Float32
                            )

    n_psms = 0
    for file_path in file_paths
        n_psms += length(Arrow.Table(file_path)[1])
    end

    all_probs = Vector{@NamedTuple{prob::Float32, target::Bool}}(undef, n_psms)
    n = 1
    for file_path in file_paths
        df = Arrow.Table(file_path)
        for i in range(1, length(df[1]))
            all_probs[n] = (prob = df[:prob][i], target = df[:target][i])
            n += 1
        end
    end

    # Sort probabilities in descending order
    sort!(all_probs, alg = QuickSort, rev=true)

    # Second pass: Find the probability threshold
    targets_above = 0
    decoys_above = 0
    prob_for_q_value_threshold = 0.0f0
    for (prob, target) in all_probs

        targets_above += 1
        decoys_above += (1 - target)
        
        current_q_value = decoys_above / (targets_above + decoys_above)
        
        if current_q_value > q_value_threshold
           return prob#return prob  # This is the probability threshold we're looking for
        end
    end

    # If we get here, all q-values are below the threshold
    return prob_for_q_value_threshold
end


function getPSMsPassingQVal(
                            quant_psms_folder::String, 
                            passing_psms_folder::String,
                            q_val_threshold::Float32,
                            )
    
    file_paths = readdir(quant_psms_folder, join=true)

    prob_threshold = find_prob_threshold(file_paths, q_val_threshold)

    # Initialize an empty DataFrame to store the results
    result_df = DataFrame()

    for file_path in file_paths
        # Read the Arrow table
        arrow_table = Arrow.Table(file_path)
        
               
        # Sample the rows and convert to DataFrame
        sampled_df = select!(DataFrame(arrow_table),
        [
            :precursor_idx,
            :prob,
            :target,
            :irt_obs,
            :missed_cleavage,
            :ms_file_idx])[(arrow_table[:prob].>=prob_threshold), :]
        
        # Append to the result DataFrame
        Arrow.write(
            joinpath(passing_psms_folder, basename(file_path)),
            sampled_df
        )
    end

    return
end


