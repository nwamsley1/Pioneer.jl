function quantPSMs(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    precursors::Union{Arrow.Table, Missing},
                    library_fragment_lookup::Union{LibraryFragmentLookup{Float32}, Missing},
                    ms_file_idx::UInt32,
                    rt_to_irt::UniformSpline,
                    mass_err_model::MassErrorModel,
                    δ::Float32,
                    λ::Float32,
                    max_iter_newton::Int64,
                    max_iter_bisection::Int64,
                    max_iter_outer::Int64,
                    accuracy_newton::Float32,
                    accuracy_bisection::Float32,
                    max_diff::Float32,
                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{L},
                    iso_splines::IsotopeSplineModel,
                    chromatograms::Vector{ChromObject},
                    scored_PSMs::Vector{S},
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_frag_count::Int64,
                    min_spectral_contrast::Float32,
                    min_log2_matched_ratio::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    max_best_rank::Int64,
                    n_frag_isotopes::Int64,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    irt_tol::Float64,
                    spec_order::Set{Int64}
                    ) where {T,U<:AbstractFloat, L<:LibraryIon{Float32},
                    S<:ScoredPSM{Float32, Float16},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}

    ##########
    #Initialize 
    prec_idx, ion_idx, cycle_idx, last_val = 0, 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    isotopes = zeros(Float32, 5)
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""

    rt_idx = 0
    prec_temp_size = 0
    precs_temp = Vector{UInt32}(undef, 50000)

    
    ##########
    #Iterate through spectra
    for scan_idx in thread_task
        if scan_idx == 0 
            continue
        end
        if scan_idx > length(spectra[:masses])
            continue
        end
        ###########
        #Scan Filtering
        msn = spectra[:msOrder][scan_idx] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        if (msn < 2)
            cycle_idx += 1
        end
        #cycle_idx += (msn == 1)
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        irt = rt_to_irt(spectra[:retentionTime][scan_idx])
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
        prec_mz_string_new = string(spectra[:centerMass][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:max(length(prec_mz_string_new), 6)]

        if (irt_start_new != irt_start) | (irt_stop_new != irt_stop) | (prec_mz_string_new != prec_mz_string)
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new
            #Candidate precursors and their retention time estimates have already been determined from
            #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
            #the retention time and m/z tolerance constraints
            precs_temp_size = 0
            ion_idx, prec_idx, prec_temp_size = selectRTIndexedTransitions!(
                ionTemplates,
                precs_temp,
                precs_temp_size,
                library_fragment_lookup,
                precursors[:mz],
                precursors[:prec_charge],
                precursors[:sulfur_count],
                iso_splines,
                isotopes,
                n_frag_isotopes,
                rt_index,
                irt_start,
                irt_stop,
                spectra[:centerMass][scan_idx] - spectra[:isolationWidth][scan_idx]/2.0f0,
                spectra[:centerMass][scan_idx] + spectra[:isolationWidth][scan_idx]/2.0f0,
                (
                    spectra[:lowMass][scan_idx], spectra[:highMass][scan_idx]
                ),
                isotope_err_bounds,
                10000)
        end

        #if scan_idx != 213864
        #    continue
        #end
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:masses][scan_idx], 
                                        spectra[:intensities][scan_idx], 
                                        mass_err_model,
                                        spectra[:highMass][scan_idx],
                                        UInt32(scan_idx), 
                                        ms_file_idx)
        sort!(@view(ionMatches[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches > 2 #Few matches to do not perform de-convolution 
        
            #Spectral deconvolution. Build sparse design/template matrix for regression 
            #Sparse matrix representation of templates written to Hs. 
            #IDtoCOL maps precursor ids to their corresponding columns. 
            buildDesignMatrix!(Hs, ionMatches, ionMisses, nmatches, nmisses, IDtoCOL)
            #Adjuste size of pre-allocated arrays if needed 
            if IDtoCOL.size > length(_weights_)
                new_entries = IDtoCOL.size - length(_weights_) + 1000 
                append!(_weights_, zeros(eltype(_weights_), new_entries))
                append!(spectral_scores, Vector{eltype(spectral_scores)}(undef, new_entries))
                append!(unscored_PSMs, [eltype(unscored_PSMs)() for _ in 1:new_entries]);
            end
            #Get most recently determined weights for each precursors
            #"Hot" start
            for i in range(1, IDtoCOL.size)#pairs(IDtoCOL)
                _weights_[IDtoCOL[IDtoCOL.keys[i]]] = precursor_weights[IDtoCOL.keys[i]]
            end
            #Get initial residuals
            initResiduals!(_residuals_, Hs, _weights_);
            #Spectral deconvolution. Hybrid bisection/newtowns method
            solveHuber!(Hs, _residuals_, _weights_, 
                            δ, λ, 
                            max_iter_newton, 
                            max_iter_bisection,
                            max_iter_outer,
                            accuracy_newton,
                            accuracy_bisection,
                            10.0,#Hs.n/10.0,
                            max_diff
                            );
            #return Hs, IDtoCOL, _weights_
            test_prec_id = 2930911
            if scan_idx == 213864
                p = plot(title = "scan id $scan_idx")
                reconstructed_spectra = zeros(Float32, length(spectra[:masses][scan_idx]))
                for i in range(1, nmatches)
                    peak_id = getPeakInd(ionMatches[i])
                    prec_id = getPrecID(ionMatches[i])
                    intensity = ionMatches[i].predicted_intensity
                    w = _weights_[IDtoCOL[prec_id]]
                    reconstructed_spectra[peak_id] += w*intensity
                end
                w = _weights_[IDtoCOL[test_prec_id]]
                for i in range(1, nmatches)
                    if getPrecID(ionMatches[i])==test_prec_id
                        match = ionMatches[i]
                              plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 1, alpha = 1.0, lw = 3, label = nothing)
                              plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, reconstructed_spectra[getPeakInd(match)]], color = 2, alpha = 0.5, lw = 3, label = nothing)
                            plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 3, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                for i in range(1, nmisses)
                    if getPrecID(ionMisses[i])==test_prec_id
                        match = ionMisses[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 3, alpha = 0.5, lw = 3, label = nothing)
                        annotate!(p, match.theoretical_mz, match.predicted_intensity*w, "*")
                    end
                end
                plot!(p, [ionMatches[1].theoretical_mz], [0], color = 1, label = "Experimental Spectrum")
                plot!(p, [ionMatches[1].theoretical_mz], [0], color = 2, label = "Reconstructed Spectrum")
                plot!(p, [ionMatches[1].theoretical_mz], [0], color = 3, label = "Target Spectrum")
                #=
                for i in range(1, nmisses)
                    if getPrecID(ionMisses[i])==test_prec_id
                        match = ionMisses[i]

                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, reconstructed_spectra[getPeakInd(match)]], color = 1, alpha = 0.5, lw = 3, label = nothing)
                        #plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                =#
                savefig(p, "/Users/n.t.wamsley/Desktop/test_reconstruct_"*string(test_prec_id)*".pdf")


                prec_col = IDtoCOL[test_prec_id]
                w = _weights_[prec_col]
                p = plot(title = "scan id $scan_idx, weight: $w")
                for i in range(1, nmatches)
                    if getPrecID(ionMatches[i])==test_prec_id
                        match = ionMatches[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 1, alpha = 0.5, lw = 3, label = nothing)
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                for i in range(1, nmisses)
                    if getPrecID(ionMisses[i])==test_prec_id
                        match = ionMisses[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 1, alpha = 0.5, lw = 3, label = nothing)
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                savefig(p, "/Users/n.t.wamsley/Desktop/test_huber1100_lambda1e3_213864_"*string(test_prec_id)*".pdf")
            end
            if scan_idx == 213863
               
                prec_col = IDtoCOL[test_prec_id]
                w = _weights_[prec_col]
                p = plot(title = "scan id $scan_idx, weight: $w")
                for i in range(1, nmatches)
                    if getPrecID(ionMatches[i])==test_prec_id
                        match = ionMatches[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 1, alpha = 0.5, lw = 3, label = nothing)
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                for i in range(1, nmisses)
                    if getPrecID(ionMisses[i])==test_prec_id
                        match = ionMisses[i]
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.predicted_intensity*w], color = 1, alpha = 0.5, lw = 3, label = nothing)
                        plot!(p, [match.theoretical_mz, match.theoretical_mz], [0.0, match.intensity], color = 2, alpha = 0.5, lw = 3, label = nothing)
                    end
                end
                savefig(p, "/Users/n.t.wamsley/Desktop/test_huber1100_lambda1e3_213863_"*string(test_prec_id)*".pdf")
            end
                #for i in range(Hs.colptr[prec_col], Hs.colptr[prec_col + 1] - 1)
                #    plot(p, H.x[i]
                #end
            for j in range(1, prec_temp_size)
                if !iszero(IDtoCOL[precs_temp[j]])
                    rt_idx += 1
                    chromatograms[rt_idx] = ChromObject(
                        Float16(spectra[:retentionTime][scan_idx]),
                        _weights_[IDtoCOL[precs_temp[j]]],
                        scan_idx,
                        precs_temp[j]
                    )
                else
                    rt_idx += 1
                    chromatograms[rt_idx] = ChromObject(
                        Float16(spectra[:retentionTime][scan_idx]),
                        zero(Float32),
                        scan_idx,
                        precs_temp[j]
                    )
                end
                if rt_idx + 1 > length(chromatograms)
                    growChromObjects!(chromatograms, 500000)
                end

            end
            #Record weights for each precursor
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
            end
            getDistanceMetrics(_weights_, Hs, spectral_scores)
            ##########
            #Scoring and recording data
            ScoreFragmentMatches!(unscored_PSMs,
                                IDtoCOL,
                                ionMatches, 
                                nmatches, 
                                mass_err_model,
                                last(min_topn_of_m)
                                )

            last_val = Score!(scored_PSMs, 
                unscored_PSMs,
                spectral_scores,
                _weights_,
                IDtoCOL,
                cycle_idx,
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(spectra[:intensities][scan_idx])), 
                scan_idx,
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                max_best_rank = max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000,
                )
        else
            for j in range(1, prec_temp_size)
                rt_idx += 1
                chromatograms[rt_idx] = ChromObject(
                    Float16(spectra[:retentionTime][scan_idx]),
                    zero(Float32),
                    scan_idx,
                    precs_temp[j]
                )
                if rt_idx + 1 > length(chromatograms)
                    growChromObjects!(chromatograms, 500000)
                end
            end
        end
        ##########
        #Reset pre-allocated arrays 
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        reset!(IDtoCOL);
        reset!(Hs);
    end

    return DataFrame(@view(scored_PSMs[1:last_val])), DataFrame(@view(chromatograms[1:rt_idx]))
end

@time RESULT = quantitationSearch(
    MS_TABLE, 
    params_;
    precursors = prosit_lib["precursors"],
    fragment_lookup_table = library_fragment_lookup_table,
    rt_index = RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model = frag_err_dist_dict[ms_file_idx],
    irt_err = irt_err,#irt_errs[ms_file_idx]/3,
    ion_matches = ionMatches,
    ion_misses = ionMisses,
    id_to_col = IDtoCOL,
    ion_templates = ionTemplates,
    iso_splines = iso_splines,
    chromatograms = chromatograms,
    scored_psms = complex_scored_PSMs,
    unscored_psms = complex_unscored_PSMs,
    spectral_scores = complex_spectral_scores,
    precursor_weights = precursor_weights,
    );


    histogram(frag_ppm_errs, xlim = (-20, 20))
    frag_abs_errs = [(match.theoretical_mz-match.match_mz)/match.charge for match in matched_fragments];
    histogram(frag_abs_errs, xlim = (-0.02, 0.02))
    vline!([quantile(frag_abs_errs, 0.99)], color = :black, lw = 3)
    vline!([quantile(frag_abs_errs, 0.01)], color = :black, lw = 3)
    frag_abs_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments if match.frag_charge == 1];
    histogram(frag_abs_errs, xlim = (-20, 20), alpha = 0.5, normalize = :pdf)
    frag_abs_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in matched_fragments if match.frag_charge == 2];
    histogram!(frag_abs_errs, xlim = (-20, 20), alpha = 0.5, normalize = :pdf)
    vline!([quantile(frag_abs_errs, 0.99)], color = :black, lw = 3)
    vline!([quantile(frag_abs_errs, 0.01)], color = :black, lw = 3)


    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz)for match in matched_fragments];
    histogram2d(frag_ppm_errs, [match.theoretical_mz for match in matched_fragments], xlim = (-20, 20), alpha = 1.0, normalize = :pdf)