function scoreProteinGroups!(bpsms::DataFrame)
    #New columns
    bpsms[!,:peptide_count] = zeros(UInt8, size(bpsms, 1))
    bpsms[!,:max_pg_score] = zeros(Float32, size(bpsms, 1))
    gbpsms = groupby(
                        best_psms[!,
                        [:ms_file_idx,:target,:accession_numbers,:sequence,:prob,:max_pg_score,:peptide_count]], 
                        [:ms_file_idx,:target,:accession_numbers]
                    )
                    
    function setScore!(
        max_pg_scores::AbstractVector{Float32},
        peptide_counts::AbstractVector{UInt8},
        precursor_score::AbstractVector{Float32},
        sequence::AbstractVector{String})

        max_pg_score = maximum(precursor_score)
        peptide_count =length(unique(sequence))
        for j in range(1, length(max_pg_scores))
            max_pg_scores[j] = max_pg_score
            peptide_counts[j] = min(typemax(UInt8), peptide_count)
        end
    end
    for (key, psms) in pairs(gbpsms)
        setScore!(
            psms[!,:max_pg_score],
            psms[!,:peptide_count],
            psms[!,:prob],
            psms[!,:sequence]
        )
    end

    scored_proteins = select!(combine(gbpsms, first), Not([:sequence,:prob]))
    filter!(x->x.peptide_count>1,scored_proteins)
    scored_proteins[!,:q_value] = zeros(Float32, size(scored_proteins, 1))
    getQvalues!(scored_proteins[!,:max_pg_score], scored_proteins[!,:target], scored_proteins[!,:q_value]);
    filter!(x->x.q_value<=0.01, scored_proteins)
    return scored_proteins
    
end