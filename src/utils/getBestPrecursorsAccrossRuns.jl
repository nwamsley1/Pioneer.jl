function getBestPrecursorsAccrossRuns(psms_paths::Dictionary{String, String},
                         prec_mzs::AbstractVector{Float32},
                         RT_iRT = Dict{Int64, Any};
                         max_q_val::Float32 = 0.01f0,
                         max_precursors::Int = 250000)

    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{ best_prob::Float32, 
                                                        best_irt::Float32, 
                                                        min_irt::Union{Missing, Float32}, 
                                                        max_irt::Union{Missing, Float32}, 
                                                        mz::Float32}}()
    #prec_to_best_prob = zeros(Float32, n_precursors)
    for (key, psms_path) in pairs(psms_paths) #For each data frame 
        psms = Arrow.Table(psms_path)
        #One row for each precursor 
        for row in eachindex(psms[:precursor_idx])

            #Precursor data 
            precursor_idx = psms[:precursor_idx][row]
            prob = psms[:prob][row]
            mz = prec_mzs[precursor_idx]
            max_irt, min_irt = missing, missing
            irt =  RT_iRT[key](psms[:RT][row])

            
            if haskey(prec_to_best_prob, precursor_idx)
                prec_prob = prec_to_best_prob[precursor_idx] 
                best_prob = prec_prob[:best_prob]
                best_irt = prec_prob[:best_irt]
                if psms[:q_value][row] <= max_q_val
                    if irt < coalesce(prec_prob[:min_irt], typemax(Float32))
                        min_irt = irt
                    else
                        min_irt = prec_prob[:min_irt]
                    end
                    if irt > coalesce(prec_prob[:max_irt], typemin(Float32))
                        max_irt = irt
                    else
                        max_irt = prec_prob[:max_irt]
                    end
                end
                if (best_prob < prob)
                    best_prob = prob
                    best_irt = irt
                end
                prec_to_best_prob[precursor_idx] = (
                    best_prob = best_prob,
                    best_irt = best_irt,
                    min_irt = min_irt,
                    max_irt = max_irt,
                    mz = mz)
            else
                min_irt, max_irt = missing, missing
                if psms[:q_value][row] <= max_q_val
                    min_irt, max_irt = irt, irt
                end
                val = (best_prob = prob, best_irt = irt, min_irt = min_irt, max_irt = max_irt, mz = mz)
                insert!(prec_to_best_prob, precursor_idx, val)
            end
        end
    end
    # Get the top N peptide_idx's with their corresponding maximum probabilities
    #Currently this is slow
    sort!(prec_to_best_prob, alg=PartialQuickSort(1:max_precursors), rev = true);
    N = 0
    for key in collect(keys(prec_to_best_prob))
        N += 1
        if N > max_precursors
            delete!(prec_to_best_prob, key)
        end
    end



    return prec_to_best_prob #[(prob, idx) for (idx, prob) in sort(collect(top_probs), rev=true)]
end