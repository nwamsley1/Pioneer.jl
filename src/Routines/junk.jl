function getBestTrace!(psms::DataFrame,
                       max_q_value::AbstractFloat,
                       quant_col::Symbol)

    function getScore(
        abundance::AbstractVector{Float32},
        q_value::AbstractVector{Float32},
        max_q_value::AbstractFloat)
        
        score = zero(Float32)
        for i in range(1, length(abundance))
            #if q_value[i] <= max_q_value
            if q_value[i] > 0.75
                score += abundance[i]
            end
        end 
        return score
    end
    function fillBestTrace!(
        best_psm::AbstractVector{Bool})
        for i in range(1, length(best_psm))
            best_psm[i] = true
        end
    end

    @time sort!(psms,[:precursor_idx,:isotopes_captured])
    psms[!,:best_trace] .= false
    grouped_psms = groupby(psms, :precursor_idx)
    #partition(1:length(grouped_psms), chunk_size)
    #For each precursor
    for i in range(1, length(grouped_psms))

        #For all psms for the given precursor
        #group by isotope subsets 
        iso_sets = groupby(grouped_psms[i],:isotopes_captured)
        best_iso_set = nothing
        best_score = typemin(Float32)
        for (iso_set, subpsms) in pairs(iso_sets)
            #Score is sum of peak areas where the q_value was 
            #below the threshold 
            score = getScore(
                subpsms[!,quant_col],
                #subpsms[!,:q_value],
                subpsms[!,:prob],
                max_q_value
            )
            #If the score is the best encountered so far, 
            #then this it he best isotope trace so far. 
            if score > best_score
                best_score = score
                best_iso_set = iso_set
            end
        end

        if best_iso_set !== nothing
            fillBestTrace!(iso_sets[best_iso_set][!,:best_trace])
        end
    end
    filter!(x->x.best_trace, psms);
end