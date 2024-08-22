
function find_score_threshold(
    scores::Vector{@NamedTuple{score::Float32, target::Bool}},
    q_value_threshold::Float32
)
    # Second pass: Find the probability threshold
    targets_above = 0
    decoys_above = 0
    targets = 0
    decoys = 0
    for (prob, target) in scores
        targets += target
        decoys += (1 - target)
    end

    for (score, target) in reverse(scores)

        targets_above -= target
        decoys_above -= (1 - target)
        
        current_q_value = decoys_above / (targets_above + decoys_above)
        
        if current_q_value > q_value_threshold
           return score#return prob  # This is the probability threshold we're looking for
        end
    end

    return scores[end][:score]
end

function find_prob_threshold(
                            file_paths::Vector{String},
                            best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}},
                            q_value_threshold::Float32
                            )

    n_psms = 0
    for file_path in file_paths
        n_psms += length(Arrow.Table(file_path)[1])
    end

    all_probs = Vector{@NamedTuple{score::Float32, target::Bool}}(undef, n_psms)
    n = 1
    for file_path in file_paths
        df = Arrow.Table(file_path)
        for i in range(1, length(df[1]))
            key = (precursor_idx = df[:precursor_idx][i], isotopes_captured = df[:isotopes_captured][i])
            if key ∈ best_traces
                all_probs[n] = (score = df[:prob][i], target = df[:target][i])
                n += 1
            end
        end
    end
    all_probs = all_probs[1:n]
    # Sort probabilities in descending order
    sort!(all_probs, by=x->x[:score], alg = QuickSort)

    return find_score_threshold(all_probs, q_value_threshold)
end



function getPSMsPassingQVal(
                            quant_psms_folder::String, 
                            passing_psms_folder::String,
                            best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}},
                            q_val_threshold::Float32,
                            )
    
    file_paths = readdir(quant_psms_folder, join=true)

    prob_threshold = find_prob_threshold(
                                            file_paths, 
                                            best_traces,
                                            q_val_threshold)
    println("prob_threshold $prob_threshold")
    # Initialize an empty DataFrame to store the results
    result_df = DataFrame()

    for file_path in file_paths
        # Read the Arrow table
        arrow_table = Arrow.Table(file_path)
        
        passing_psms = zeros(Bool, length(arrow_table[:prob]))
        for i in range(1, length(best_traces))
            if arrow_table[:prob][i].<prob_threshold
                continue
            end
            key = (precursor_idx = arrow_table[:precursor_idx][i], isotopes_captured = arrow_table[:isotopes_captured][i])
            if key ∈ best_traces
                passing_psms[i] = true
            end
        end
        # Sample the rows and convert to DataFrame
        sampled_df = select!(DataFrame(arrow_table),
        [
            :precursor_idx,
            :prob,
            :weight,
            :target,
            :irt_obs,
            :missed_cleavage,
            :isotopes_captured,
            :scan_idx,
            :ms_file_idx])[passing_psms, :]
        
        # Append to the result DataFrame
        Arrow.write(
            joinpath(passing_psms_folder, basename(file_path)),
            sampled_df
        )
    end

    return
end


