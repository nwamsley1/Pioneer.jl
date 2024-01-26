function SearchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursor{Float32}}, Missing},
                    ion_list::Union{Vector{Vector{LibraryFragment{Float32}}}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::Laplace{Float64},
                    selectIons!::Function,
                    searchScan!::Union{Function, Missing},
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{LibraryFragment{Float32}}},
                    iso_splines::IsotopeSplineModel{Float64},
                    scored_PSMs::Vector{Vector{LibPSM{Float32, Float16}}},
                    unscored_PSMs::Vector{Vector{LXTandem{Float32}}},
                    spectral_scores::Vector{Vector{SpectralScores{Float16}}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, Float32}}};
                    #keyword args
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float64 = 0.0,
                    fragment_tolerance::Float64 = 20.0,
                    huber_δ::Float32 = 1000f0,
                    unmatched_penalty_factor::Float64 = 1.0
                    IonMatchType::DataType = FragmentMatch{Float32},
                    IonTemplateType::DataType = LibraryFragment{Float32},
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                    max_iter::Int = 1000,
                    max_peak_width::Float64 = 2.0,
                    max_peaks::Union{Int64,Bool} = false, 
                    min_frag_count::Int64 = 4,
                    min_frag_count_index_search::Int64 = 0,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_index_search_score::Float32 = zero(Float32),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    filter_by_rank::False, 
                    min_weight::Float32 = zero(Float32),
                    most_intense = false,
                    n_frag_isotopes::Int64 = 1,
                    #precs::Counter{UInt32, UInt8, Float32} = Counter(UInt32, UInt8, Float32, 0),
                    #_precs::Union{Counter{UInt32, Float32}, Missing} = missing,
                    precursor_tolerance::Float64 = 5.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    rt_tol::Float64 = 30.0,
                    sample_rate::Float64 = 1.0,
                    scan_range::Tuple{Int64, Int64} = (0, 0),
                    #scored_PSMs::Union{Dict{Symbol, Vector}, Missing} = missing,
                    spec_order::Set{Int64} = Set(2),
                    topN::Int64 = 20,
                    topN_index_search::Int64 = 1000) where {T,U<:AbstractFloat}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    peaks = sum(length.(spectra[:masses]))
    peaks_per_thread = peaks÷(16)#Threads.nthreads()÷2)
    thread_tasks = []
    n = 0
    start = 1
    #Incorperate sacn range here 
    for i in range(1, length(spectra[:masses]))
        n += length(spectra[:masses][i])
        if (n > peaks_per_thread) | ((i + 1) == length(spectra[:masses]))
            push!(thread_tasks, (start, length(thread_tasks) + 1, i))
            start = i + 1
            n = 0
        end
    end
    println("length(thread_tasks) ", length(thread_tasks))
    pbar = ProgressBar(total = peaks)
    lk = ReentrantLock()

    if ismissing(precs)
        precs = [missing for _ in range(1, length(thread_tasks))]
    end
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            return searchRAW(
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
                                min_matched_ratio,min_index_search_score,min_spectral_contrast,min_topn_of_m,filter_by_rank,
                                min_weight, most_intense,n_frag_isotopes,precursor_tolerance,quadrupole_isolation_width,
                                rt_bounds,rt_index, rt_tol,sample_rate,
                                spec_order,topN,topN_index_search
                            )
            #lock(lk) do 
            #    push!(out, df_out)
            #end
        end
    end
    return fetch.(tasks) 
end

function searchRAW(
                    spectra::Arrow.Table,
                    lk::ReentrantLock,
                    #isprecs::Bool,
                    pbar::ProgressBar,
                    thread_task::Tuple{Int64, Int64},
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursor{Float32}}, Missing},
                    ion_list::Union{Vector{Vector{LibraryFragment{Float32}}}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::Laplace{Float64},
                    selectIons!::Function,
                    searchScan!::Union{Function, Missing},
                    collect_fmatches::Bool,
                    expected_matches::Int64,
                    frag_ppm_err::Float64,
                    fragment_tolerance::Float64,
                    huber_δ::Float32,
                    unmatched_penalty_factor::Float64,
                    IonMatchType::DataType,


                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    all_fmatches::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{LibraryFragment{Float32}},
                    iso_splines::IsotopeSplineModel{Float64},
                    scored_PSMs::Vector{LibPSM{Float32, Float16}},
                    unscored_PSMs::Vector{LXTandem{Float32}},
                    spectral_scores::Vector{SpectralScores{Float16}},
                    precursor_weights::Vector{Float32},
                    precs::Union{Missing, Counter{UInt32, Float32}},

                    IonTemplateType::DataType,
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    max_peak_width::Float64,
                    max_peaks::Union{Int64,Bool},
                    min_frag_count::Int64,
                    min_frag_count_index_search::Int64,
                    min_matched_ratio::Float32,
                    min_index_search_score::Float32,
                    min_spectral_contrast::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    filter_by_rank::Bool,
                    min_weight::Float32,
                    most_intense::Bool,
                    n_frag_isotopes::Int64,
                    precursor_tolerance::Float64,
                    quadrupole_isolation_width::Float64,
                    rt_bounds::Tuple{Float64, Float64},
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    rt_tol::Float64,
                    sample_rate::Float64,
                    spec_order::Set{Int64},
                    topN::Int64,
                    topN_index_search::Int64,
                    ) where {T,U<:AbstractFloat}

    thread_peaks = 0
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    prec_idx, ion_idx, cycle_idx, last_val = 0, 0, 0, 0
    minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds) #only consider scans in the bounds
    ###########
    prec_ids = [zero(UInt32) for _ in range(1, expected_matches)]
    Hs = SparseArray(5000);
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    solve_time, index_search_time, prep_time, index_ions_time = zero(Float64), zero(Float64),zero(Float64),zero(Float64)
    ##########
    #Iterate through spectra
    scans_processed = 0
    isotopes = zeros(Float64, n_frag_isotopes)
    #ncols = 0
    for i in range(first(thread_task), last(thread_task))
        #if (i < 100000) | (i > 100020)
        #   continue
        #end
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
        
        #(i >= first(thread_task)) & (i <= last(thread_task)) ? nothing : continue #Skip if outside the scan range
        
        first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

        min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

        iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        if !ismissing(searchScan!) | !ismissing(frag_index)
            searchScan!(
                        precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        min_intensity, spectra[:masses][i], spectra[:intensities][i],
                        iRT_low, iRT_high,
                        #Float32(frag_ppm_err),
                        Float32(fragment_tolerance),
                        spectra[:precursorMZ][i],
                        Float32(quadrupole_isolation_width/2.0),
                        isotope_err_bounds,
                        min_frag_count = min_frag_count_index_search, 
                        min_ratio = Float32(min_index_search_score),
                        topN = topN_index_search,#topN
                        )
            #println("prec_count $prec_count match_count $match_count ")
        end
        #selectIons! 
        #Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
        if !ismissing(precs) 
            ion_idx, prec_idx = selectIons!(ionTemplates, 
                                                precursors,
                                                ion_list,
                                                iso_splines,
                                                isotopes,
                                                precs,
                                                topN,
                                                Float32(iRT_to_RT_spline(spectra[:retentionTime][i])),
                                                Float32(rt_tol), #rt_tol
                                                (
                                                spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                ),
                                                isotope_err_bounds = isotope_err_bounds
                                                )::Tuple{Int64, Bool}
            #println("ion_idx $ion_idx prec_idx $prec_idx")
        else
            ion_idx, prec_idx = selectIons!(
                                            ionTemplates,
                                            precursors,
                                            ion_list,
                                            iso_splines,
                                            isotopes,
                                            prec_ids,
                                            rt_index,
                                            spectra[:retentionTime][i],
                                            Float32(max_peak_width/2),
                                            (
                                                spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                ),
                                        isotope_err_bounds = isotope_err_bounds)
        end
        ion_idx < 2 ? continue : nothing 
        scans_processed += 1
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks(ionTemplates, #Search the spectra for these ions 
                                    ion_idx, #search ionTemplates[1:ion_idx]
                                    ionMatches, #Fill with matched ions 
                                    ionMisses, #Fill with unmatched ions 
                                    spectra[:masses][i], 
                                    spectra[:intensities][i], 
                                    count_unmatched=true, #Should we fill "ionMisses"?
                                    δs = [frag_ppm_err], #Mass offsets 
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity, #Ignore peaks below this intensity
                                    ppm = fragment_tolerance, #Fragment match tolerance in ppm
                                    most_intense = most_intense
                                    )
                                    nmisses_all = nmisses
                                    nmatches_all = nmatches
        
        nmatches_all, nmisses_all, nmatches, nmisses = filterMatchedIons!(IDtoCOL, ionMatches, ionMisses, nmatches, nmisses, 
                                                                                10000, #Arbitrarily hight
                                                                                min_frag_count #Remove precursors matching fewer than this many fragments
                                                                        )
        if filter_by_rank
        #println("nmatches_all $nmatches_all, nmatches $nmatches")
            _, _, nmatches, nmisses = filterMatchedIons(IDtoCOL, ionMatches, ionMisses, nmatches, nmisses, 
                                                        last(min_topn_of_m), 
                                                        first(min_topn_of_m),
                                                        )
        end

        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches > 2 #Few matches to do not perform de-convolution 
            

            #Spectral deconvolution. Build sparse design/template matrix for regression 
            #Sparse matrix representation of templates written to Hs. 
            #IDtoCOL maps precursor ids to their corresponding columns. 
            buildDesignMatrix!(Hs, ionMatches, ionMisses, nmatches, nmisses, IDtoCOL,
                                unmatched_penalty_factor = unmatched_penalty_factor)
            #println("Hs.n ", Hs.n)
            if IDtoCOL.size > length(_weights_)
                append!(_weights_, zeros(eltype(_weights_), IDtoCOL.size - length(_weights_) + 1000 ))
                append!(spectral_scores, Vector{SpectralScores{Float16}}(undef, IDtoCOL.size - length(spectral_scores) + 1000 ))
                append!(unscored_PSMs, [LXTandem(Float32) for _ in 1:(IDtoCOL.size - length(unscored_PSMs) + 1000)]);
            end
            #Get most recently determined weights for eahc precursors
            #"Hot" start
            for i in range(1, IDtoCOL.size)#pairs(IDtoCOL)
                _weights_[IDtoCOL[IDtoCOL.keys[i]]] = precursor_weights[IDtoCOL.keys[i]]
            end

            #Get initial residuals
            initResiduals!(_residuals_, Hs, _weights_);
            if ismissing(precs) 
                #Spectral deconvolution.
                solve_time += @elapsed solveHuber!(Hs, _residuals_, _weights_, huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);
            end
            #Remember determined weights for eahc precursors
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
            end

            for col in range(1, Hs.n)
                for k in range(Hs.colptr[col], Hs.colptr[col + 1]-1)
                    if iszero(Hs.matched[k])
                        Hs.nzval[k] = unmatched_penalty_factor*Hs.nzval[k]
                    end
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
                                    err_dist,
                                    last(min_topn_of_m)
                                    )

                last_val = Score!(scored_PSMs, 
                    unscored_PSMs,
                    spectral_scores,
                    _weights_,
                    #match_count/prec_count,
                    nmatches/(nmatches + nmisses),
                    last_val,
                    Hs.n,
                    Float32(sum(spectra[:intensities][i])), 
                    i,
                    min_spectral_contrast = min_spectral_contrast, #Remove precursors with spectral contrast lower than this ammount
                    min_matched_ratio = min_matched_ratio,
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
            unscored_PSMs[i] = LXTandem(Float32)
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

function firstSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursor{Float32}},
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{LibraryFragment{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{LibPSM{Float32, Float16}}},
    unscored_PSMs::Vector{Vector{LXTandem{Float32}}},
    spectral_scores::Vector{Vector{SpectralScores{Float16}}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, Float32}};
    scan_range = (0, 0))

    return SearchRAW(
        spectra, 
        frag_index,
        precursors, 
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        selectTransitions!,
        searchScan!,

        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        precs,

        collect_fmatches = true,
        expected_matches = params[:expected_matches],
        frag_ppm_err = params[:frag_ppm_err],
        fragment_tolerance = params[:frag_tol_presearch],
        max_iter = params[:max_iter],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_frag_count_index_search = params[:min_frag_count_index_search],
        min_matched_ratio = params[:min_matched_ratio],
        min_index_search_score = params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast],
        min_topn_of_m = params[:min_topn_of_m],

        most_intense = false,#params[:most_intense],
        #precs = Counter(UInt32, UInt8, Float32, length(ion_list)),
        #_precs = Counter(UInt32, Float32, length(ion_list)),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        scan_range = scan_range,
        #scored_PSMs = makePSMsDict(XTandem(Float32)),
        topN = params[:topN],
        topN_index_search = params[:topN_index_search],
    )
end

function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursor{Float32}},
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    fragment_tolerance::Float64,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{LibraryFragment{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{LibPSM{Float32, Float16}}},
    unscored_PSMs::Vector{Vector{LXTandem{Float32}}},
    spectral_scores::Vector{Vector{SpectralScores{Float16}}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, Float32}};
    scan_range::Tuple{Int64, Int64} = (0, 0))

    frag_ppm_err = err_dist.μ

    return SearchRAW(
        spectra, 
        frag_index, 
        precursors,
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        selectTransitions!,
        searchScan!,

        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        precs,

        expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        huber_δ = params[:huber_δ],
        isotope_err_bounds = params[:isotope_err_bounds],
        max_iter = params[:max_iter],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_frag_count_index_search = params[:min_frag_count_index_search],
        min_matched_ratio = params[:min_matched_ratio_main_search],
        min_index_search_score = params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast],
        min_topn_of_m = params[:min_topn_of_m],
        most_intense = params[:most_intense],
        n_frag_isotopes = params[:n_frag_isotopes],
        #precs = Counter(UInt32, Float32, length(ion_list)),
        #precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = 1.0,
        scan_range = scan_range,
        #scored_PSMs = makePSMsDict(XTandem(Float32)),
        topN = params[:topN],
        topN_index_search = params[:topN_index_search]
    )
end

function integrateMS2(
    #Mandatory Args
    spectra::Arrow.Table,
    precursors::Vector{LibraryPrecursor{Float32}},
    ion_list::Vector{Vector{LibraryFragment{Float32}}},
    rt_index::retentionTimeIndex{U, T},
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    fragment_tolerance::Float64,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{LibraryFragment{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{LibPSM{Float32, Float16}}},
    unscored_PSMs::Vector{Vector{LXTandem{Float32}}},
    spectral_scores::Vector{Vector{SpectralScores{Float16}}},
    precursor_weights::Vector{Vector{Float32}};
    N = 600000*10,
    scan_range =  (0, 0)) where {U,T<:AbstractFloat}

    frag_ppm_err = err_dist.μ
    #fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])

    return SearchRAW(
        spectra, 
        missing,
        precursors,
        ion_list, 
        x->x,
        ms_file_idx,
        err_dist,
        selectRTIndexedTransitions!,
        missing,
        
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        missing,
        
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        unmatched_penalty_factor = params[:unmatched_penalty_factor],
        isotope_err_bounds = params[:isotope_err_bounds],
        max_iter = params[:max_iter],
        max_peak_width = params[:max_peak_width],
        max_peaks = params[:max_peaks],
        min_topn_of_m = params[:min_topn_of_m],
        filter_by_rank = true,
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_index_search_score = 0.6f0,#params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast], 
        min_weight = params[:min_weight],
        most_intense = params[:most_intense],
        n_frag_isotopes = params[:n_frag_isotopes],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_index = rt_index,
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = Inf,
        scan_range = scan_range,
        topN = params[:topN],
    )
end

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

function integrateMS1(
    #Mandatory Args
    spectra::Arrow.Table, 
    isotope_dict::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
    prec_rt_list::Vector{Tuple{Union{Missing,T}, UInt32}} ,
    ms_file_idx::UInt32,
    err_dist::Laplace{Float64},
    params::Dict; 
    N = 600000*10,
    scan_range = (0, 0)) where {T<:AbstractFloat}

    frag_ppm_err = err_dist.μ
    fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])

    #=
    integrate_ms1_params = (
    expected_matches = 1000000,
    frag_err_dist = frag_err_dist_dict[1],
    frag_tol_quantile = 0.975,
    max_iter = 1000,
    max_peak_width = 2.0,
    max_peaks = false,
    min_frag_count = 4,
    min_matched_ratio = Float32(0.45),
    min_spectral_contrast = Float32(0.5),
    nmf_tol = Float32(100),
    precursor_tolerance = 5.0,
    quadrupole_isolation_width = 4.25,
    regularize = false,
    rt_tol = 20.0,
    sample_rate = 1.0,
    scan_range = (0, 3000000),
    topN = 100,
    λ = zero(Float32),
    γ = zero(Float32)

    integrateMS1(
    Arrow.Table(MS_TABLE_PATHS[1]), 
    prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  
    isotopes,
    prec_rt_table, 
    RT_to_iRT_map_dict[1], #RT to iRT map'
    UInt32(1), #MS_FILE_IDX
    frag_err_dist_dict[1],
    integrate_ms1_params)
    )=#
    return SearchRAW(
        spectra, 
        missing, 
        missing, #Not ion list. Instead passing "isotope_dict"
        x->x,
        ms_file_idx,
        err_dist,  #Not really sure how to estimate this yet?
        selectIsotopes!, #Ion Selection Function for MS1 integration 
        missing,
        
        chromatograms =  Dict(:precursor_idx => zeros(UInt32, N), 
        :weight => zeros(Float32, N), 
        :scan_idx => zeros(UInt32, N),
        :rt => zeros(Float32, N), 
        :frag_count => zeros(Int64, N),
        :rank => zeros(UInt8, N),
        :cycle_idx => zeros(UInt32, N)),
                expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        fragment_tolerance = fragment_tolerance,
        IonTemplateType = Isotope{Float32},
        IonMatchType = PrecursorMatch{Float32},
        isotope_dict = isotope_dict,
        max_iter = params[:max_iter],
        max_peak_width = params[:max_peak_width],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_spectral_contrast = params[:min_spectral_contrast],
        most_intense = params[:most_intense],
        nmf_tol = params[:nmf_tol],
        precursor_tolerance = params[:precursor_tolerance],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        regularize = params[:regularize],
        rt_bounds = params[:rt_bounds],
        rt_index = prec_rt_list,
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        scan_range = scan_range,
        spec_order = Set(1),
        topN = params[:topN],
        λ = params[:λ],
        γ = params[:γ]
    )
end

function fillChroms!(chroms::Dict{Symbol, Vector}, id_to_row::UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}, n::Int64, scan_idx::Int64, cycle_idx::Int64, prec_ids::Vector{UInt32}, prec_idx::Int64, frag_counts::Accumulator{UInt32,Int64}, weights::Vector{T}, retention_time::U; block_size = 100000) where {T,U<:AbstractFloat}
    function inc!(chroms::Dict{Symbol, Vector}, n::Int64, scan_idx::Int64, cycle_idx::Int64, key::UInt32, weight::AbstractFloat, rt::AbstractFloat, frag_count::Int64,rank::UInt8)
        chroms[:precursor_idx][n] = key
        chroms[:weight][n] = weight
        chroms[:rt][n] = rt
        chroms[:scan_idx][n] = scan_idx
        chroms[:frag_count][n] = frag_count
        chroms[:rank][n] = rank
        chroms[:cycle_idx][n] = UInt32(cycle_idx)
    end

    for i in range(1, prec_idx)
        key = prec_ids[i]
        
        if n > length(chroms[:precursor_idx])
            #This block of code is triggered when not enough space has been allocated for the chromatograms
            #So incrase the space allocated by a given "block" size. 
            for key in keys(chroms)
                append!(chroms[key], zeros(eltype(chroms[key]), block_size))
            end
        end
        if haskey(frag_counts, key)
    
            if haskey(id_to_row, key)
                inc!(chroms, n, scan_idx, cycle_idx, key, weights[id_to_row[key][1]], retention_time, frag_counts[key], id_to_row[key][2])
            else
                inc!(chroms, n, scan_idx, cycle_idx, key, Float32(0.0),retention_time, frag_counts[key], zero(UInt8))
            end
        else
            inc!(chroms, n, scan_idx, cycle_idx, key, Float32(0.0), retention_time, 0, zero(UInt8))
        end
        n += 1
    end    
    return n
end

function getRTWindow(irt::U, max_rt::T, min_rt::T, rt_tol::T) where {T,U<:AbstractFloat}
    #(irt < min_rt) ? irt_low = Float32(-Inf) : irt_low = Float32(irt - rt_tol)
    #(irt > max_rt) ? irt_high = Float32(Inf) : irt_high = Float32(irt + rt_tol)
    #return irt_low, irt_high
    return Float32(irt - rt_tol), Float32(irt + rt_tol)
end

function getMinIntensity(intensities::AbstractArray{Union{T, Missing}}, max_peaks::Int) where {T<:AbstractFloat}
    return intensities[sortperm(intensities, rev = true)[min(max_peaks, length(intensities))]]
end

function getMinIntensity(intensities::AbstractArray{Union{T, Missing}}, max_peaks::Bool) where {T<:AbstractFloat}
    return zero(T)
end

function reset!(fms::Vector{M}, last_non_empty::Int64) where {M<:Match}
    for i in range(1, last_non_empty)
        fms[i] = M()
    end
end

#=function reset!(fms::Vector{PrecursorMatch{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        fms[i] = PrecursorMatch{T}()
    end
end=#

function reset!(fms::Vector{Isotope{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        fms[i] = Isotope{T}()
    end
end

function reset!(lf::Vector{LibraryFragment{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        lf[i] = LibraryFragment{T}()
    end
end

function collectFragErrs(all_fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int, collect_fmatches::Bool) where {M<:Match}
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

#=

function SearchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursor{Float32}}, Missing},
                    ion_list::Union{Vector{Vector{LibraryFragment{Float32}}}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::Laplace{Float64},
                    selectIons!::Function,
                    searchScan!::Union{Function, Missing};
                    #keyword args
                    chromatograms::Union{Dict{Symbol, Vector}, Missing} = missing,
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float64 = 0.0,
                    fragment_tolerance::Float64 = 20.0,
                    huber_δ::Float32 = 1000f0,
                    IonMatchType::DataType = FragmentMatch{Float32},
                    IonTemplateType::DataType = LibraryFragment{Float32},
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    max_iter::Int = 1000,
                    max_peak_width::Float64 = 2.0,
                    max_peaks::Union{Int64,Bool} = false, 
                    min_frag_count::Int64 = 4,
                    min_frag_count_index_search::Int64 = 0,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_index_search_score::Float32 = zero(Float32),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    min_topn::Int64 = 2,
                    min_weight::Float32 = zero(Float32),
                    most_intense = false,
                    nmf_tol::Float32 = Float32(100.0),
                    #precs::Counter{UInt32, UInt8, Float32} = Counter(UInt32, UInt8, Float32, 0),
                    precs::Union{Counter{UInt32, Float32}, Missing} = missing,
                    precursor_tolerance::Float64 = 5.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    regularize::Bool = false,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    rt_tol::Float64 = 30.0,
                    sample_rate::Float64 = 1.0,
                    scan_range::Tuple{Int64, Int64} = (0, 0),
                    scored_PSMs::Union{Dict{Symbol, Vector}, Missing} = missing,
                    spec_order::Set{Int64} = Set(2),
                    topN::Int64 = 20,
                    topN_index_search::Int64 = 1000,
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32)) where {T,U<:AbstractFloat}
    println("max_peak_width $max_peak_width")

    function distributeScans(N::Int64, m::Int64)
        scan_idx = 0
        scans = Vector{NamedTuple{(:first_scan, :last_scan, :thread_id), Tuple{Int64, Int64, Int64}}}()
        proc = 0
        i = 1
        first_scan = 1
        last_scan = 1
        while scan_idx <= N
            if (i > m) | (scan_idx == N)
                i = 0
                proc += 1
                push!(scans, (first_scan = first_scan, 
                              last_scan = scan_idx,
                              thread_id = proc%Threads.nthreads() + 1)
                      )
                scan_idx += 1
                first_scan = scan_idx
                continue
            end
            scan_idx += 1
            i += 1
        end
        return scans
    end
    
      
   
    precs = [Counter(UInt32, Float32, length(ion_list)) for _ in range(1, Threads.nthreads())]
    println(length(precs))
    ##########
    #Initialize 

    println("min_spectral_contrast ", min_spectral_contrast) #Remove precursors with spectral contrast lower than this ammount
    println("min_matched_ratio ", min_matched_ratio) 
    println("min_frag_count ", min_frag_count) 
    println("min_weight ", min_weight) 
    println("min_topn ", min_topn) 
    ###########
    #Pre-allocate Arrays to save (lots) of time in garbage collection. 
    all_fmatches = Vector{IonMatchType}()
    collect_fmatches ? all_fmatches = [[IonMatchType() for x in range(1, expected_matches)] for _ in range(1, Threads.nthreads())] : nothing

    #These are overwritten for every searched spectrum. "expected_matches"
    #is a guess for the largest array size that would be needed for any spectrum. 
    #If the guess is too small, the arrays will simply be increased in size as needed
    #by a pre-determined block-size. 
    println("START")
    ionMatches = [[IonMatchType() for _ in range(1, expected_matches)]  for _ in range(1, Threads.nthreads())]#IonMatchType is something that inherits from the "Match" class. 
    ionMisses = [[IonMatchType() for _ in range(1, expected_matches)] for _ in range(1, Threads.nthreads())]
    ionTemplates = [[IonTemplateType() for _ in range(1, expected_matches)]  for _ in range(1, Threads.nthreads())]
    prec_ids = [[zero(UInt32) for _ in range(1, expected_matches)] for _ in range(1, Threads.nthreads())]

    IDtoCOL = nothing
    if ismissing(precs)
        IDtoCOL = [ArrayDict(UInt32, UInt16, length(precursors)) for _ in range(1, Threads.nthreads())]
    else
        IDtoCOL = [ArrayDict(UInt32, UInt16, length(first(precs).ids)) for _ in range(1, Threads.nthreads())]
    end
    #H_COLS, H_ROWS, H_VALS, H_MASK = zeros(Int64, expected_matches), zeros(Int64, expected_matches), zeros(Float32, expected_matches), zeros(Float32, expected_matches)
    scored_PSMs = [Vector{LibPSM{Float32, Float16}}(undef, 5000) for _ in range(1, Threads.nthreads())];
    unscored_PSMs = [[LXTandem(Float32) for _ in 1:5000] for _ in range(1, Threads.nthreads())];
    spectral_scores = [Vector{SpectralScores{Float16}}(undef, 5000) for _ in range(1, Threads.nthreads())];
    Hs = [SparseArray(50000) for _ in range(1, Threads.nthreads())];
    println("STOP")

    #weights
    precursor_weights = ""
    if ismissing(ion_list)
        precursor_weights = [zeros(Float32, maximum(keys(isotope_dict))) for _ in range(1, Threads.nthreads())]
    else
        precursor_weights = [zeros(Float32, length(ion_list)) for _ in range(1, Threads.nthreads())]
    end
    _weights_ = [zeros(Float32, 5000) for _ in range(1, Threads.nthreads())];
    _residuals_ = [zeros(Float32, 5000) for _ in range(1, Threads.nthreads())];
    last_vals = [0 for _ in range(1, Threads.nthreads())];
    #fragment_intensities = Dictionary{String, Vector{Tuple{Float32, Float32}}}()
    ##########
    #Iterate through spectra
    scan_batches = distributeScans(length(spectra[:masses]), 10000)
    Threads.@threads :static for i in ProgressBar(1:length(scan_batches))
        ThreadID = scan_batches[i][:thread_id] 
        scan_range = (scan_batches[i][:first_scan], scan_batches[i][:last_scan])
        msms_counts = Dict{Int64, Int64}()
        frag_err_idx = 1
        prec_idx = 0
        ion_idx = 0
        cycle_idx = 0 
        minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds)
        println("scan range ", scan_range, " assigned to thread ", ThreadID)
        for i in range(scan_batches[i][:first_scan], scan_batches[i][:last_scan])
            ###########
            #Scan Filtering
            msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
            if msn == 1
                cycle_idx += 1
            end
            msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
            msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

            (i >= first(scan_range)) & (i <= last(scan_range)) ? nothing : continue #Skip if outside the scan range
            first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

            min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

            iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

            ##########
            #Ion Template Selection
            #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  
            if !ismissing(searchScan!) | !ismissing(frag_index)
                prec_count, match_count = searchScan!(
                            precs[ThreadID], #counter which keeps track of plausible matches 
                            frag_index, 
                            min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
                            iRT_low, iRT_high,
                            Float32(fragment_tolerance), 
                            Float32(precursor_tolerance),
                            Float32(quadrupole_isolation_width/2.0),
                            min_frag_count = min_frag_count_index_search, 
                            min_ratio = Float32(min_index_search_score),
                            topN = topN_index_search,#topN
                            )
                
            end
            #selectIons! 
            #Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
            if !ismissing(precs) 
                ion_idx, prec_idx = selectIons!(ionTemplates[ThreadID], 
                                                                            precursors,
                                                                            ion_list,
                                                                            precs[ThreadID],
                                                                            topN,
                                                                            Float32(iRT_to_RT_spline(spectra[:retentionTime][i])),
                                                                            Float32(rt_tol), #rt_tol
                                                                            (
                                                                            spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                                            spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                                            )
                                                                            )::Tuple{Int64, Bool}
            else
                ion_idx, prec_idx = selectIons!(
                                                ionTemplates[ThreadID],
                                                precursors,
                                                ion_list,
                                                prec_ids[ThreadID],
                                                rt_index,
                                                spectra[:retentionTime][i],
                                                Float32(max_peak_width/2),
                                                spectra[:precursorMZ][i],
                                                Float32(quadrupole_isolation_width/2.0),
                                                (
                                                    spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                                    spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                                    )
                                                )
            end
            ion_idx < 2 ? continue : nothing 

            ##########
            #Match sorted list of plausible ions to the observed spectra
            nmatches, nmisses = matchPeaks(ionTemplates[ThreadID], #Search the spectra for these ions 
                                                                ion_idx, #search ionTemplates[1:ion_idx]
                                                                ionMatches[ThreadID], #Fill with matched ions 
                                                                ionMisses[ThreadID], #Fill with unmatched ions 
                                                                spectra[:masses][i], 
                                                                spectra[:intensities][i], 
                                                                count_unmatched=true, #Should we fill "ionMisses"?
                                                                δs = [frag_ppm_err], #Mass offsets 
                                                                scan_idx = UInt32(i),
                                                                ms_file_idx = ms_file_idx,
                                                                min_intensity = min_intensity, #Ignore peaks below this intensity
                                                                ppm = fragment_tolerance, #Fragment match tolerance in ppm
                                                                most_intense = most_intense
                                                                )

            ##########
            #Spectral Deconvolution and Distance Metrics 
            if nmatches > 2 #Few matches to do not perform de-convolution 
                
                #Spectral deconvolution. Build sparse design/template matrix for regression 
                #Sparse matrix representation of templates written to Hs. 
                #IDtoCOL maps precursor ids to their corresponding columns. 
                buildDesignMatrix!(Hs[ThreadID], ionMatches[ThreadID], ionMisses[ThreadID], nmatches, nmisses, IDtoCOL[ThreadID])

                if IDtoCOL[ThreadID].size > length(_weights_[ThreadID])
                    append!(_weights_[ThreadID], zeros(eltype(_weights_[ThreadID]), IDtoCOL[ThreadID].size - length(_weights_[ThreadID]) + 1000 ))
                    append!(spectral_scores[ThreadID], Vector{SpectralScores{Float16}}(undef, IDtoCOL[ThreadID].size - length(spectral_scores[ThreadID]) + 1000 ))
                    append!(unscored_PSMs[ThreadID], [LXTandem(Float32) for _ in 1:(IDtoCOL[ThreadID].size - length(unscored_PSMs[ThreadID]) + 1000)]);
                end
                #Get most recently determined weights for eahc precursors
                #"Hot" start
                for i in range(1, IDtoCOL[ThreadID].size)#pairs(IDtoCOL)
                    _weights_[ThreadID][IDtoCOL[ThreadID][IDtoCOL[ThreadID].keys[i]]] = precursor_weights[ThreadID][IDtoCOL[ThreadID].keys[i]]
                end

                #Get initial residuals
                initResiduals!(_residuals_[ThreadID], Hs[ThreadID], _weights_[ThreadID]);
                if ismissing(precs) 
                    #Spectral deconvolution.
                    solveHuber!(Hs[ThreadID], _residuals_[ThreadID], _weights_[ThreadID], 
                                                        huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs[ThreadID].n);
                end
                #Remember determined weights for eahc precursors
                for i in range(1, IDtoCOL[ThreadID].size)
                    precursor_weights[ThreadID][IDtoCOL[ThreadID].keys[i]] = _weights_[ThreadID][IDtoCOL[ThreadID][IDtoCOL[ThreadID].keys[i]]]# = precursor_weights[id]
                end

                if ismissing(isotope_dict) 
                    getDistanceMetrics(_weights_[ThreadID], Hs[ThreadID], spectral_scores[ThreadID])
                end

                ##########
                #Scoring and recording data
                if !ismissing(scored_PSMs)

                    ScoreFragmentMatches!(unscored_PSMs[ThreadID],
                                        IDtoCOL[ThreadID],
                                        ionMatches[ThreadID], 
                                        nmatches, 
                                        err_dist)

                    last_vals[ThreadID] = Score!(scored_PSMs[ThreadID], 
                        unscored_PSMs[ThreadID],
                        spectral_scores[ThreadID],
                        _weights_[ThreadID],
                        #match_count/prec_count,
                        nmatches/(nmatches + nmisses),
                        last_vals[ThreadID],
                        Hs[ThreadID].n,
                        Float32(sum(spectra[:intensities][i])), 
                        i,
                        min_spectral_contrast = min_spectral_contrast, #Remove precursors with spectral contrast lower than this ammount
                        min_matched_ratio = min_matched_ratio,
                        min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                        min_weight = min_weight,
                        min_topn = min_topn,
                        block_size = 500000,
                        )
                end
            end
            #Add fragment matches to all_fmatches 
            frag_err_idx = collectFragErrs(all_fmatches[ThreadID], ionMatches[ThreadID], nmatches, frag_err_idx, collect_fmatches)
            
            ##########
            #Update Chromatograms 
            #=
            if !ismissing(chromatograms)
                frag_counts = counter(UInt32) #Basically a Dictionary that counts the number of matched ions (values) for each precursor (keys)
                for match_idx in range(1, nmatches) #fragmentMatches
                    DataStructures.inc!(frag_counts, ionMatches[match_idx].prec_id)
                end
                #Add precursor templates with their weights and retention times to the chromatogram table 
                chrom_idx = fillChroms!(chromatograms, IDtoCOL, chrom_idx, i, cycle_idx, prec_ids, prec_idx, frag_counts, weights, spectra[:retentionTime][i])
            end
            =#
            ##########
            #Reset pre-allocated arrays 
            reset!(ionTemplates[ThreadID], ion_idx)
            reset!(ionMatches[ThreadID], nmatches)
            reset!(ionMisses[ThreadID], nmisses)
            fill!(prec_ids[ThreadID], zero(UInt32))
            for i in range(1, Hs[ThreadID].n)
                unscored_PSMs[ThreadID][i] = LXTandem(Float32)
            end
            reset!(IDtoCOL[ThreadID]);
            reset!(Hs[ThreadID]);
        end
    end
    #return fragment_intensities
    ############
    #Return Chromatograms and Score/Feature Table
    #return all_fmatches
    if collect_fmatches
        #return DataFrame(scored_PSMs), all_fmatches
        return  DataFrame(@view(scored_PSMs[1][1:last_vals[1]])), all_fmatches[1]
    else
        return DataFrame(@view(scored_PSMs[1][1:last_vals[1]]))
    end
    #=
        if ismissing(chromatograms)
            #return all_fmatches
            #return DataFrame(scored_PSMs)
            DataFrame(@view(scored_PSMs[1:last_val]))
        elseif ismissing(scored_PSMs)
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return groupby(DataFrame(chromatograms), :precursor_idx)
        else
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return DataFrame(scored_PSMs), groupby(DataFrame(chromatograms), :precursor_idx)
        end
    =#
end

function SearchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    ion_list::Union{Vector{Vector{LibraryFragment{Float32}}}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::Laplace{Float64},
                    selectIons!::Function,
                    searchScan!::Union{Function, Missing};
                    #keyword args
                    b_min_ind::Int64 = 3,
                    y_min_ind::Int64 = 4,
                    chromatograms::Union{Dict{Symbol, Vector}, Missing} = missing,
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float64 = 0.0,
                    fragment_tolerance::Float64 = 20.0,
                    huber_δ::Float32 = 1000f0,
                    IonMatchType::DataType = FragmentMatch{Float32},
                    IonTemplateType::DataType = LibraryFragment{Float32},
                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    max_iter::Int = 1000,
                    max_peak_width::Float64 = 2.0,
                    max_peaks::Union{Int64,Bool} = false, 
                    min_frag_count::Int64 = 4,
                    min_frag_count_index_search::Int64 = 0,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_matched_ratio_index_search::Float32 = zero(Float32),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    most_intense = false,
                    nmf_tol::Float32 = Float32(100.0),
                    precs::Counter{UInt32, UInt8, Float32} = Counter(UInt32, UInt8, Float32, 0),
                    precursor_tolerance::Float64 = 5.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    regularize::Bool = false,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    rt_tol::Float64 = 30.0,
                    sample_rate::Float64 = 1.0,
                    scan_range::Tuple{Int64, Int64} = (0, 0),
                    scored_PSMs::Union{Dict{Symbol, Vector}, Missing} = missing,
                    spec_order::Set{Int64} = Set(2),
                    topN::Int64 = 20,
                    topN_index_search::Int64 = 1000,
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32)) where {T,U<:AbstractFloat}
    println("max_peak_width $max_peak_width")

    ##########
    #Initialize 
    weights = Float32[]
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    chrom_idx = 1
    prec_idx = 0
    ion_idx = 0
    cycle_idx = 0
    minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds)

    fragment_tolerance = Float64(8.1)
    ###########
    #Pre-allocate Arrays to save (lots) of time in garbage collection. 
    all_fmatches = Vector{IonMatchType}()
    collect_fmatches ? all_fmatches = [IonMatchType() for x in range(1, expected_matches)] : nothing

    #These are overwritten for every searched spectrum. "expected_matches"
    #is a guess for the largest array size that would be needed for any spectrum. 
    #If the guess is too small, the arrays will simply be increased in size as needed
    #by a pre-determined block-size. 
    ionMatches = [IonMatchType() for _ in range(1, expected_matches)] #IonMatchType is something that inherits from the "Match" class. 
    ionMisses = [IonMatchType() for _ in range(1, expected_matches)]
    ionTemplates = [IonTemplateType() for _ in range(1, expected_matches)] 
    prec_ids = [zero(UInt32) for _ in range(1, expected_matches)]
    H_COLS, H_ROWS, H_VALS, H_MASK = zeros(Int64, expected_matches), zeros(Int64, expected_matches), zeros(Float32, expected_matches), zeros(Float32, expected_matches)
    

    #weights
    precursor_weights = ""
    if ismissing(ion_list)
        precursor_weights = zeros(Float32, maximum(keys(isotope_dict)))
    else
        precursor_weights = zeros(Float32, length(ion_list))
    end

    #fragment_intensities = Dictionary{String, Vector{Tuple{Float32, Float32}}}()
    solve_time = 0.0
    index_search_time = 0.0
    prep_time = 0.0
    index_ions_time = 0.0
    ##########
    #Iterate through spectra
    #for i in ProgressBar(range(first(scan_range), last(scan_range)))
    for i in range(first(scan_range), last(scan_range))
    #for i in range(1, size(spectra[:masses])[1])

        ###########
        #Scan Filtering
        #(i%10000) == 0 ? println(i) : nothing
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        if msn == 1
            cycle_idx += 1
        end
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type

        (i >= first(scan_range)) & (i <= last(scan_range)) ? nothing : continue #Skip if outside the scan range
        first(rand(1)) <= sample_rate ? nothing : continue #dice-roll. Usefull for random sampling of scans. 

        min_intensity = getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

        iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  
        if !ismissing(searchScan!) | !ismissing(frag_index)
            index_search_time += @elapsed prec_count, match_count = searchScan!(
                        precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
                        iRT_low, iRT_high,
                        Float32(fragment_tolerance), 
                        Float32(precursor_tolerance),
                        Float32(quadrupole_isolation_width/2.0),
                        min_frag_count = min_frag_count_index_search, 
                        min_ratio = Float32(min_matched_ratio_index_search),
                        topN = topN_index_search,#topN
                        )
            
        end
        #selectIons! 
        #Get a sorted list by m/z of ion templates (fills ionTemplates). The spectrum will be searched for matches to these ions only.
        if ismissing(isotope_dict) 
            
            index_ions_time += @elapsed ion_idx, prec_idx = selectIons!(ionTemplates, 
                                               ion_list,
                                               precs,
                                               topN,
                                               prec_ids,
                                               rt_index,
                                               spectra[:retentionTime][i],
                                               Float32(max_peak_width/2.0), #rt_tol
                                               spectra[:precursorMZ][i], #prec_mz
                                               Float32(quadrupole_isolation_width/2.0) #prec_tol
                                               )::Tuple{Int64, Bool}
        else
            ion_idx, prec_idx = selectIons!(
                                            ionTemplates,
                                            rt_index,
                                            isotope_dict,
                                            prec_ids,
                                            spectra[:retentionTime][i],
                                            Float32(max_peak_width/2.0)
                                            )
        end
        ion_idx < 2 ? continue : nothing 

        ##########
        #Match sorted list of plausible ions to the observed spectra
        prep_time += @elapsed nmatches, nmisses = matchPeaks(ionTemplates, #Search the spectra for these ions 
                                    ion_idx, #search ionTemplates[1:ion_idx]
                                    ionMatches, #Fill with matched ions 
                                    ionMisses, #Fill with unmatched ions 
                                    spectra[:masses][i], 
                                    spectra[:intensities][i], 
                                    count_unmatched=true, #Should we fill "ionMisses"?
                                    δs = [frag_ppm_err], #Mass offsets 
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity, #Ignore peaks below this intensity
                                    ppm = fragment_tolerance, #Fragment match tolerance in ppm
                                    most_intense = most_intense
                                    )
        #if prec_idx
        #    println("nmatches $nmatches, nmisses $nmisses, scan_idx $i")
        #    println("nmatches $nmatches, nmisses $nmisses, scan_idx $i")
        #end

        ##########
        #Spectral Deconvolution and Distance Metrics 
        if nmatches < 2 #Few matches to do not perform de-convolution 
            IDtoCOL = UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()
        else #Spectral deconvolution. Build sparse design/template matrix for nnls regression 
            #if prec_idx
            #    println("PASSED scan_idx $i")
            #    println("PASSED scan_idx $i")
            #end
            #IDtoCOL_weights = UnorderedDictionary{UInt32, UInt32}()
            #prep_time += @elapsed begin

                #Sparse matrix representation of templates written to Hs. 
                #Hs mask is true for ions that will be used to score and false otherwise (may exclude y1 ions for example)
                #IDtoCOL maps precursor ids to their corresponding columns. 
                #X, Hs, Hs_mask, IDtoCOL, last_matched_col = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, H_COLS, H_ROWS, H_VALS, H_MASK, y_min_ind = y_min_ind, b_min_ind = b_min_ind)
            X, Hs, IDtoCOL, last_matched_col = buildDesignMatrix(ionMatches, ionMisses, nmatches, nmisses, H_COLS, H_ROWS, H_VALS)
            return X, Hs, IDtoCOL, last_matched_col, ionMatches, ionMisses, nmatches, nmisses
            #if (Hs.n > 5000)
            #    println("Hs.n ", Hs.n)
            #    println("Hs.n ", Hs.n)
            #    println("scan_id $i")
            #    println("scan_id $i")
            #end
            #println("Hs.n ", Hs.n)
            #println("Hs.m ", Hs.m)
            
            weights = Vector{Float32}(undef, Hs.n)
            #end

            for (id, row) in pairs(IDtoCOL)
                weights[first(row)] = precursor_weights[id]
            end

            solve_time += @elapsed solveHuber!(Hs, Hs*weights .- X, weights, huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);

            for (id, row) in pairs(IDtoCOL)
                precursor_weights[id] = weights[first(row)]# = precursor_weights[id]
            end

            #return IDtoCOL, weights, Hs, X, r, last_matched_col
            #return Hs
            if ismissing(isotope_dict) 
                scores = getDistanceMetrics(X, weights, Hs, last_matched_col)
            end
            #for (id, row) in pairs(IDtoCOL_weights)
            #    precursor_weights[id] = weights[row]# = precursor_weights[id]
            #end
            #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            ##########
            #Scoring and recording data
            if !ismissing(scored_PSMs)
                unscored_PSMs = UnorderedDictionary{UInt32, XTandem{Float32}}()

                ScoreFragmentMatches!(unscored_PSMs, ionMatches, nmatches, err_dist)
                #Score unscored_PSMs and write them to scored_PSMs
                #return scores
                #if prec_idx
                #    println("SCORED scan_idx $i")
                #    println("SCORED scan_idx $i")
                #    println("unscored_PSMs ", unscored_PSMs)
                #end
                Score!(scored_PSMs, 
                        unscored_PSMs, 
                        length(spectra[:intensities][i]), 
                        Float64(sum(spectra[:intensities][i])), 
                        match_count/prec_count, 
                        scores, #Named Tuple of spectrum simmilarity/distance measures 
                        weights, #Coefficients for each precursor in the spectral deconvolution
                        IDtoCOL,
                        scan_idx = i,
                        min_spectral_contrast = min_spectral_contrast, #Remove precursors with spectral contrast lower than this ammount
                        min_matched_ratio = min_matched_ratio,
                        min_frag_count = min_frag_count #Remove precursors with fewer fragments 
                        )
            end
        end
        #Add fragment matches to all_fmatches 
        frag_err_idx = collectFragErrs(all_fmatches, ionMatches, nmatches, frag_err_idx, collect_fmatches)
        
        ##########
        #Update Chromatograms 
        if !ismissing(chromatograms)
            frag_counts = counter(UInt32) #Basically a Dictionary that counts the number of matched ions (values) for each precursor (keys)
            for match_idx in range(1, nmatches) #fragmentMatches
                DataStructures.inc!(frag_counts, ionMatches[match_idx].prec_id)
            end
            #Add precursor templates with their weights and retention times to the chromatogram table 
            chrom_idx = fillChroms!(chromatograms, IDtoCOL, chrom_idx, i, cycle_idx, prec_ids, prec_idx, frag_counts, weights, spectra[:retentionTime][i])
        end

        ##########
        #Reset pre-allocated arrays 
        reset!(ionTemplates, ion_idx)
        reset!(ionMatches, nmatches), reset!(ionMisses, nmisses)
        fill!(prec_ids, zero(UInt32))

    end

    println("solve_time $solve_time")
    println("prep_time $prep_time")
    println("index_search_time $index_search_time")
    println("index_ions_time $index_ions_time")
    #return fragment_intensities
    ############
    #Return Chromatograms and Score/Feature Table
    #return all_fmatches
    if collect_fmatches
        return DataFrame(scored_PSMs), all_fmatches
    else
        if ismissing(chromatograms)
            #return all_fmatches
            return DataFrame(scored_PSMs)
        elseif ismissing(scored_PSMs)
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return groupby(DataFrame(chromatograms), :precursor_idx)
        else
            chromatograms = DataFrame(chromatograms)
            sort!(chromatograms, [:precursor_idx,:rt], alg=QuickSort);
            return DataFrame(scored_PSMs), groupby(DataFrame(chromatograms), :precursor_idx)
        end
    end

end
=#