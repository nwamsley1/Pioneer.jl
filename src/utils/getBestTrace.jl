
function getBestTrace!(psms::DataFrame,
                       score_col::Symbol)
    function getScore(
        scores::AbstractVector{Float32}
        )
        
        score_sum = zero(Float32)
        for score in enumerate(scores)
            score_sum += log2(score)
        end 
        return score_sum
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
                subpsms[!,:q_value],
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
    #filter!(x->x.best_trace, psms);
end


function getBestTraces(
    quant_psms_folder::String
)
    psms_trace_scores = Dictionary{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()
    for file_path in readdir(quant_psms_folder, join = true)
        if splitext(file_path)[end] != ".arrow"
            continue
        end
        psms_table = Arrow.Table(file_path)
        for i in range(1, length(psms_table[1]))
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])
            row_score = log2(psms_table[:prob][i])
            if isnan(row_score)
                row_score = zero(Float32)
            end
            if haskey(psms_trace_scores, psms_key)
                score = psms_trace_scores[psms_key]
                score += row_score
                psms_trace_scores[psms_key] = score
            else
                insert!(
                    psms_trace_scores,
                    psms_key,
                    row_score
                )
            end
        end
    end

    psms_trace_df = DataFrame(
    (precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
    isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
    score = [val for val in values(psms_trace_scores)])
    );
    psms_trace_df[!,:best_trace] .= false;
    gpsms = groupby(psms_trace_df,:precursor_idx)
    for (precursor_idx, psms) in pairs(gpsms)
        psms[argmax(psms[!,:score]),:best_trace] = true
    end
    filter!(x->x.best_trace, psms_trace_df);
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
end




psms_trace_df = DataFrame(
(precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
 isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
 score = [val for val in values(psms_trace_scores)])
);
psms_trace_df[!,:best_trace] .= false;
gpsms = groupby(psms_trace_df,:precursor_idx)
for (precursor_idx, psms) in pairs(gpsms)
    psms[argmax(psms[!,:score]),:best_trace] = true
end
filter!(x->x.best_trace, psms_trace_df);
traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);

for file_path in readdir(second_quant_folder, join = true)
    if splitext(file_path)[end] != ".arrow"
        continue
    end
    psms_table = DataFrame(
        Tables.columntable(#Need a modifiable deep copy
                Arrow.Table(
                    file_path
                    )))
    psms_table[!,:to_keep] .= false

    for i in range(1, size(psms_table, 1))
        psms_key = (precursor_idx = psms_table[i,:precursor_idx],  isotopes_captured = psms_table[i,:isotopes_captured])
        if psms_key ∈traces_passing
            psms_table[i,:to_keep] = true
        end
    end
    filter!(x->x.to_keep, psms_table)
    Arrow.write(
        file_path,
        psms_table
    )
end

