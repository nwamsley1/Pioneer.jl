function searchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursorIon{Float32}}, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::MassErrorModel{Float32},
                    searchScan!::Union{Function, Missing},
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel{Float64},
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, UInt8}}};
                    #keyword args
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float32 = 0.0f0,
                    huber_δ::Float32 = 1000f0,
                    unmatched_penalty_factor::Float64 = 1.0,
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                    max_peak_width::Float64 = 2.0,
                    min_frag_count::Int64 = 1,
                    min_spectral_contrast::Float32  = 0f0,
                    min_log2_matched_ratio::Float32 = -Inf32,
                    min_index_search_score::UInt8 = zero(UInt8),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    filter_by_rank::Bool = false, 
                    min_weight::Float32 = zero(Float32),
                    n_frag_isotopes::Int64 = 1,
                    quadrupole_isolation_width::Float64 = 8.5,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    irt_tol::Float64 = Inf,
                    sample_rate::Float64 = Inf,
                    spec_order::Set{Int64} = Set(2)
                    ) where {T,U<:AbstractFloat, 
                                                            L<:LibraryIon{Float32}, 
                                                            S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    thread_tasks, total_peaks = partitionScansToThreads(spectra[:masses], 
                                                        Threads.nthreads(),
                                                        1)
    pbar = ProgressBar(total = total_peaks)
    lk = ReentrantLock()

    if ismissing(precs)
        precs = [missing for _ in range(1, Threads.nthreads())]
    end
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = getThreadID(thread_task)
            return searchRAW(
                                spectra,lk,pbar,
                                getRange(thread_task),
                                frag_index,precursors,
                                ion_list, rt_to_irt_spline,ms_file_idx,err_dist,
                                searchScan!,collect_fmatches,expected_matches,frag_ppm_err,
                                huber_δ, unmatched_penalty_factor,
                                ionMatches[thread_id],ionMisses[thread_id],all_fmatches[thread_id],IDtoCOL[thread_id],ionTemplates[thread_id],
                                iso_splines, scored_PSMs[thread_id],unscored_PSMs[thread_id],spectral_scores[thread_id],precursor_weights[thread_id],
                                precs[thread_id],
                                isotope_dict,
                                isotope_err_bounds,
                                max_peak_width,min_frag_count,min_spectral_contrast,
                                min_log2_matched_ratio,min_index_search_score,min_topn_of_m,min_max_ppm,filter_by_rank,
                                min_weight,n_frag_isotopes,quadrupole_isolation_width,
                                rt_bounds,rt_index, irt_tol,sample_rate,
                                spec_order
                            )
        end
    end
    return fetch.(tasks) 
end

function searchRAW(
                    spectra::Arrow.Table,
                    lk::ReentrantLock,
                    #isprecs::Bool,
                    pbar::ProgressBar,
                    thread_task::UnitRange{Int64},
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursorIon{Float32}}, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel{Float32},
                    searchScan!::Union{Function, Missing},
                    collect_fmatches::Bool,
                    expected_matches::Int64,
                    frag_ppm_err::Float32,
                    huber_δ::Float32,
                    unmatched_penalty_factor::Float64,


                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    all_fmatches::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{L},
                    iso_splines::IsotopeSplineModel{Float64},
                    scored_PSMs::Vector{S},
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    precs::Union{Missing, Counter{UInt32, UInt8}},

                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    max_peak_width::Float64,
                    min_frag_count::Int64,
                    min_spectral_contrast::Float32,
                    min_log2_matched_ratio::Float32,
                    min_index_search_score::UInt8,
                    min_topn_of_m::Tuple{Int64, Int64},
                    min_max_ppm::Tuple{Float32, Float32},
                    filter_by_rank::Bool,
                    min_weight::Float32,
                    n_frag_isotopes::Int64,
                    quadrupole_isolation_width::Float64,
                    rt_bounds::Tuple{Float64, Float64},
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    irt_tol::Float64,
                    sample_rate::Float64,
                    spec_order::Set{Int64},
                    ) where {T,U<:AbstractFloat, L<:LibraryIon{Float32},
                    S<:ScoredPSM{Float32, Float16},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}

    thread_peaks = 0
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    prec_idx, ion_idx, cycle_idx, last_val = 0, 0, 0, 0
    prec_ids = [zero(UInt32) for _ in range(1, expected_matches)]
    Hs = SparseArray(5000);
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    isotopes = zeros(Float64, n_frag_isotopes)
    ##########
    #Iterate through spectra
    for i in thread_task
        thread_peaks += length(spectra[:masses][i])
        if thread_peaks > 100000
            lock(lk) do 
                update(pbar, thread_peaks)
            end
            thread_peaks = 0
        end

        ###########
        #Scan Filtering
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        if msn == 1
            cycle_idx += 1
        end
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

        #first(rand(1)) <= sample_rate ? nothing : continue #coin flip. Usefull for random sampling of scans. 
        iRT_low, iRT_high = getRTWindow(rt_to_irt_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, irt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        if !ismissing(precs)
            #searchScan! is the MsFragger style fragment index search
            #Loads `precs` with scores for each potential precursor matching the spectrum
            searchScan!(
                        precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        spectra[:masses][i], spectra[:intensities][i],
                        iRT_low, iRT_high,
                        frag_ppm_err,
                        mass_err_model.err_qantiles[3],
                        min_max_ppm,
                        spectra[:precursorMZ][i],
                        Float32(quadrupole_isolation_width/2.0),
                        isotope_err_bounds,
                        min_score = min_index_search_score,
                        )
            #return precs
            #println("scan_idx $i precs.matches ", precs.matches)
            #println("precs.matches ", precs.matches)
            #println("precs.size ", precs.size)
            #return precs
            #For candidate precursors exceeding the score threshold in the fragment index search...
            #Get a sorted list by m/z of fragment ion templates. The spectrum will be searched for matches to these ions only.
            #Searching with p, n, and k  precursors, theoretical fragments, and peaks, 
            #searching this way is O(n + k + n*log(n))
            #searching each precursor's fragments against the spectrum individually
            #would be O(p*n + p*k) assuming each individual precursor fragment list was already sorted
            ion_idx, prec_idx = selectTransitions!(ionTemplates, 
                                                precursors,
                                                ion_list,
                                                iso_splines,
                                                isotopes,
                                                precs,
                                                Float32(rt_to_irt_spline(spectra[:retentionTime][i])),
                                                Float32(irt_tol), #rt_tol
                                                (
                                                spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                ),
                                                isotope_err_bounds = isotope_err_bounds
                                                )::Tuple{Int64, Bool}
        else
            #Candidate precursors and their retention time estimates have already been determined from
            #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
            #the retention time and m/z tolerance constraints
            ion_idx, prec_idx = selectRTIndexedTransitions!(
                                            ionTemplates,
                                            precursors,
                                            ion_list,
                                            iso_splines,
                                            isotopes,
                                            prec_ids,
                                            rt_index,
                                            spectra[:retentionTime][i],
                                            Float32(irt_tol),#Float32(max_peak_width/2),
                                            (
                                                spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                            ),
                                        isotope_err_bounds = isotope_err_bounds)
        end
        #println("ion_idx $ion_idx prec_idx $prec_idx")
        #return precs
        #If one or fewer fragment ions matched to the spectrum, don't bother
        if ion_idx < 2
            reset!(ionTemplates, ion_idx)
            if !ismissing(precs)
            reset!(precs)
            end
            continue
        end 
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:masses][i], 
                                        spectra[:intensities][i], 
                                        frag_ppm_err,
                                        mass_err_model.err_qantiles[3],
                                        min_max_ppm,
                                        UInt32(i), 
                                        ms_file_idx)
        
        #println("nmatches $nmatches nmisses $nmisses")
        nmisses_all = nmisses
        nmatches_all = nmatches    
 
        nmatches_all, nmisses_all, nmatches, nmisses = filterMatchedIons!(IDtoCOL, 
                                                                            ionMatches, ionMisses, 
                                                                            nmatches, nmisses, 
                                                                            10000, #Arbitrarily hight
                                                                            min_frag_count #Remove precursors matching fewer than this many fragments
                                                                        )
                                                                       
        if filter_by_rank
        #println("nmatches_all $nmatches_all, nmatches $nmatches")
            _, _, nmatches, nmisses = filterMatchedIons!(IDtoCOL, 
                                                        ionMatches, ionMisses, 
                                                        nmatches, nmisses, 
                                                        last(min_topn_of_m), 
                                                        first(min_topn_of_m),
                                                        )
        end
        #println("post filter nmatches $nmatches nmisses $nmisses")
        
        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches > 2 #Few matches to do not perform de-convolution 
        
            #Spectral deconvolution. Build sparse design/template matrix for regression 
            #Sparse matrix representation of templates written to Hs. 
            #IDtoCOL maps precursor ids to their corresponding columns. 
            buildDesignMatrix!(Hs, ionMatches, ionMisses, nmatches, nmisses, IDtoCOL,
            #Reduce predicted intensities for unmatched ions by this factor
                                unmatched_penalty_factor = unmatched_penalty_factor 
                                )
            #Adjuste size of pre-allocated arrays if needed 
            if IDtoCOL.size > length(_weights_)
                new_entries = IDtoCOL.size - length(_weights_) + 1000 
                append!(_weights_, zeros(eltype(_weights_), new_entries))
                append!(spectral_scores, Vector{eltype(spectral_scores)}(undef, new_entries))
                #append!(spectral_scores, Vector{eltype(spectral_scores)}(undef, (IDtoCOL.size - length(spectral_scores) + 1000 )))
                append!(unscored_PSMs, [eltype(unscored_PSMs)() for _ in 1:new_entries]);
            end
            #Get most recently determined weights for each precursors
            #"Hot" start
            for i in range(1, IDtoCOL.size)#pairs(IDtoCOL)
                _weights_[IDtoCOL[IDtoCOL.keys[i]]] = precursor_weights[IDtoCOL.keys[i]]
            end
            #Get initial residuals
            initResiduals!(_residuals_, Hs, _weights_);
            if ismissing(precs) 
                #Spectral deconvolution. Hybrid bisection/newtowns method
                solveHuber!(Hs, _residuals_, _weights_, huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);
            end
            #Remember determined weights for each precursors
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
            end

            for i in range(1, Hs.n_vals)
                if iszero(Hs.matched[i])
                    Hs.nzval[i] = unmatched_penalty_factor*Hs.nzval[i]
                end
            end

            if ismissing(isotope_dict) 
                getDistanceMetrics(_weights_, Hs, spectral_scores)
            end

            ##########
            #Scoring and recording data
            if !ismissing(scored_PSMs)
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
                    #match_count/prec_count,
                    nmatches/(nmatches + nmisses),
                    last_val,
                    Hs.n,
                    Float32(sum(spectra[:intensities][i])), 
                    i,
                    min_spectral_contrast = min_spectral_contrast,
                    min_log2_matched_ratio = min_log2_matched_ratio,
                    min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                    min_weight = min_weight,
                    min_topn = first(min_topn_of_m),
                    block_size = 500000,
                    )
            end
        end
        #Add fragment matches to all_fmatches 
        frag_err_idx = collectFragErrs(all_fmatches, ionMatches, nmatches, frag_err_idx, collect_fmatches)
    
        ##########
        #Reset pre-allocated arrays 
        reset!(ionTemplates, ion_idx)
        reset!(ionMatches, nmatches_all)
        reset!(ionMisses, nmisses_all)
        fill!(prec_ids, zero(UInt32))
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        reset!(IDtoCOL);
        reset!(Hs);
    end

    if collect_fmatches
        return DataFrame(@view(scored_PSMs[1:last_val])), @view(all_fmatches[1:frag_err_idx])
    else
        return DataFrame(@view(scored_PSMs[1:last_val]))
    end

end

function filterMatchedIons!(IDtoNMatches::ArrayDict{UInt32, UInt16}, ionMatches::Vector{FragmentMatch{Float32}}, ionMisses::Vector{FragmentMatch{Float32}}, nmatches::Int64, nmisses::Int64, max_rank::Int64, min_matched_ions::Int64)
    nmatches_all, nmisses_all = nmatches, nmisses

    @inbounds for i in range(1, nmatches_all)
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
    @inbounds for i in range(1, nmatches_all)
        if IDtoNMatches[getPrecID(ionMatches[i])] < min_matched_ions

            continue
        else
            nmatches += 1
            ionMatches[nmatches] = ionMatches[i]
        end
    end

    @inbounds for i in range(1, nmisses_all)
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

function collectFragErrs(all_fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int, collect_fmatches::Bool) where {M<:MatchIon{Float32}}
    if collect_fmatches
        for match in range(1, nmatches)
            if n < length(all_fmatches)
                all_fmatches[n] = new_fmatches[match]
                n += 1
            else
                all_fmatches = append!(all_fmatches, [M() for x in range(1, length(all_fmatches))])
            end
        end
    end
    return n
end

function getRTWindow(irt::U, irt_tol::T) where {T,U<:AbstractFloat}
    return Float32(irt - irt_tol), Float32(irt + irt_tol)
end


