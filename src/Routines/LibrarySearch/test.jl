function filterMatchedIons!(IDtoNMatches::ArrayDict{UInt32, UInt16}, ionMatches::Vector{FragmentMatch{Float32}}, ionMisses::Vector{FragmentMath{Float32}}, nmatches::Int64, nmisses::Int64, max_rank::Int64, min_matched_ions::Int64)
    nmatches_all, nmisses_all = nmatches, nmisses

    for i in range(1, nmatches)
        match = ionMatches[i]
        prec_id = getPrecID(match)
        if match.is_isotope 
            continue
        end
            if getRank(match) <= max_rank
                if iszero(IDtoNMatches[prec_id])
                    update!(IDtoNMatches, prec_id, one(UInt16))
                else
                    IDtoNMatches.vals[prec_id] += one(UInt16)
                end
            end
    end
    nmatches, nmisses = 0, 0
    for i in range(1, nmatches_all)
        if IDtoNMatches[getPrecID(ionMatches[i])] < min_matched_ions

            continue
        else
            nmatches += 1
            ionMatches[nmatches] = ionMatches[i]
        end
    end

    for i in range(1, nmisses_all)
        if IDtoNMatches[getPrecID(ionMisses[i])] < min_matched_ions
            continue
        else
            nmisses += 1
            ionMisses[nmisses] = ionMisses[i]
        end
    end

    reset!(IDtoNMatches)

    return nmatches_all, nmisses_all, nmatches, nmisses
end

searchRAW(
    spectra,lk,pbar,
    (first(thread_task), last(thread_task)),
    frag_index,precursors,
    ion_list, iRT_to_RT_spline,ms_file_idx,err_dist,
    selectIons!,searchScan!,collect_fmatches,expected_matches,frag_ppm_err,
    fragment_tolerance,huber_δ, unmatched_penalty_factor, IonMatchType,
    ionMatches[thread_task[2]],ionMisses[thread_task[2]],all_fmatches[thread_task[2]],IDtoCOL[thread_task[2]],ionTemplates[thread_task[2]],
    iso_splines, scored_PSMs[thread_task[2]],unscored_PSMs[thread_task[2]],spectral_scores[thread_task[2]],precursor_weights[thread_task[2]],
    precs[thread_task[2]],
    IonTemplateType,isotope_dict,isotope_err_bounds,
    max_peak_width,max_peaks,min_frag_count, min_frag_count_index_search,
    min_matched_ratio,min_index_search_score,min_spectral_contrast,min_topn_of_m,
    min_weight, most_intense,n_frag_isotopes,precursor_tolerance,quadrupole_isolation_width,
    rt_bounds,rt_index, rt_tol,sample_rate,
    spec_order,topN,topN_index_search
)