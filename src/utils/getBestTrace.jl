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