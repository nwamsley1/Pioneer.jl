function SearchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursorIon{Float32}}, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    err_dist::MassErrorModel{Float32},
                    selectIons!::Function,
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
                    max_peaks::Union{Int64,Bool} = false, 
                    min_frag_count::Int64 = 4,
                    min_frag_count_index_search::Int64 = 0,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_index_search_score::UInt8 = zero(UInt8),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    filter_by_rank::Bool = false, 
                    min_weight::Float32 = zero(Float32),
                    n_frag_isotopes::Int64 = 1,
                    #precs::Counter{UInt32, UInt8, Float32} = Counter(UInt32, UInt8, Float32, 0),
                    #_precs::Union{Counter{UInt32, Float32}, Missing} = missing,
                    quadrupole_isolation_width::Float64 = 8.5,
                    rt_bounds::Tuple{Float64, Float64} = (0.0, 0.0),
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    rt_tol::Float64 = 30.0,
                    sample_rate::Float64 = 1.0,
                    #scored_PSMs::Union{Dict{Symbol, Vector}, Missing} = missing,
                    spec_order::Set{Int64} = Set(2),
                    topN::Int64 = 20,
                    topN_index_search::Int64 = 1000) where {T,U<:AbstractFloat, 
                                                            L<:LibraryIon{Float32}, 
                                                            S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    peaks = sum(length.(spectra[:masses]))
    peaks_per_thread = peaks÷(Threads.nthreads())#Threads.nthreads()÷2)
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
                                huber_δ, unmatched_penalty_factor,
                                ionMatches[thread_task[2]],ionMisses[thread_task[2]],all_fmatches[thread_task[2]],IDtoCOL[thread_task[2]],ionTemplates[thread_task[2]],
                                iso_splines, scored_PSMs[thread_task[2]],unscored_PSMs[thread_task[2]],spectral_scores[thread_task[2]],precursor_weights[thread_task[2]],
                                precs[thread_task[2]],
                                isotope_dict,
                                isotope_err_bounds,
                                max_peak_width,max_peaks,min_frag_count, min_frag_count_index_search,
                                min_matched_ratio,min_index_search_score,min_spectral_contrast,min_topn_of_m,min_max_ppm,filter_by_rank,
                                min_weight,n_frag_isotopes,quadrupole_isolation_width,
                                rt_bounds,rt_index, rt_tol,sample_rate,
                                spec_order,topN,topN_index_search
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
                    thread_task::Tuple{Int64, Int64},
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Vector{LibraryPrecursorIon{Float32}}, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    iRT_to_RT_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel{Float32},
                    selectIons!::Function,
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
                    max_peaks::Union{Int64,Bool},
                    min_frag_count::Int64,
                    min_frag_count_index_search::Int64,
                    min_matched_ratio::Float32,
                    min_index_search_score::UInt8,
                    min_spectral_contrast::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    min_max_ppm::Tuple{Float32, Float32},
                    filter_by_rank::Bool,
                    min_weight::Float32,
                    n_frag_isotopes::Int64,
                    quadrupole_isolation_width::Float64,
                    rt_bounds::Tuple{Float64, Float64},
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    rt_tol::Float64,
                    sample_rate::Float64,
                    spec_order::Set{Int64},
                    topN::Int64,
                    topN_index_search::Int64,
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
    minimum_rt, maximum_rt = first(rt_bounds), last(rt_bounds) #only consider scans in the bounds
    ###########
    prec_ids = [zero(UInt32) for _ in range(1, expected_matches)]
    Hs = SparseArray(5000);
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    ##########
    #Iterate through spectra
    scans_processed = 0
    isotopes = zeros(Float64, n_frag_isotopes)
    #ncols = 0
    for i in range(first(thread_task), last(thread_task))
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

        min_intensity = 0.0# getMinIntensity(spectra[:intensities][i], max_peaks) #Ignore peaks in the spectrum below this minimum intensity

        iRT_low, iRT_high = getRTWindow(iRT_to_RT_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, maximum_rt, minimum_rt, rt_tol) #Convert RT to expected iRT window

        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        if !ismissing(searchScan!) | !ismissing(frag_index)
            searchScan!(
                        precs, #counter which keeps track of plausible matches 
                        frag_index, 
                        min_intensity, 
                        spectra[:masses][i], spectra[:intensities][i],
                        iRT_low, iRT_high,
                        frag_ppm_err,
                        mass_err_model.err_qantiles[3],
                        min_max_ppm,
                        spectra[:precursorMZ][i],
                        Float32(quadrupole_isolation_width/2.0),
                        isotope_err_bounds,
                        min_frag_count = min_frag_count_index_search, 
                        min_score = min_index_search_score,
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
            _, _, nmatches, nmisses = filterMatchedIons!(IDtoCOL, ionMatches, ionMisses, nmatches, nmisses, 
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
                solveHuber!(Hs, _residuals_, _weights_, huber_δ, max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n);
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
                                    mass_err_model,
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

function firstSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}}
    ) where {S<:ScoredPSM{Float32, Float16},
    Q<:UnscoredPSM{Float32},
    R<:SpectralScores{Float16}}#where {S<:ScoredPSM{Float32, Float16}, LibraryIon{Float32}}

    xs = 0:1
    A = collect(xs)
    err_dist = MassErrorModel((zero(Float32), zero(Float32), zero(Float32)), zero(Float32))
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
        frag_ppm_err = Float32(params[:frag_ppm_err]),
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_frag_count_index_search = params[:min_frag_count_index_search],
        min_matched_ratio = params[:min_matched_ratio],
        min_index_search_score = params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast],
        min_topn_of_m = params[:min_topn_of_m],
        min_max_ppm = (40.0f0, 40.0f0),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = params[:sample_rate],
        #scored_PSMs = makePSMsDict(XTandem(Float32)),
        topN = params[:topN],
        topN_index_search = params[:topN_index_search],
    )
end

function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    err_dist::MassErrorModel{Float32},
    fragment_tolerance::Float64,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}};
    scan_range::Tuple{Int64, Int64} = (0, 0)) where {S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}#where {S<:ScoredPSM{Float32, Float16}, LibraryIon{Float32}}

    frag_ppm_err = getLocation(err_dist)
    
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
        isotope_err_bounds = params[:isotope_err_bounds],
        max_peaks = params[:max_peaks],
        min_frag_count = params[:min_frag_count],
        min_frag_count_index_search = params[:min_frag_count_index_search],
        min_matched_ratio = params[:min_matched_ratio_main_search],
        min_index_search_score = params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast],
        min_topn_of_m = params[:min_topn_of_m],
        min_max_ppm = (10.0f0, 20.0f0),
        n_frag_isotopes = params[:n_frag_isotopes],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = 1.0,
        topN = params[:topN],
        topN_index_search = params[:topN_index_search]
    )
end

function integrateMS2_(
    #Mandatory Args
    spectra::Arrow.Table,
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    rt_index::retentionTimeIndex{Float32, Float32},
    ms_file_idx::UInt32,
    err_dist::MassErrorModel{Float32},
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}}
    ) where {S<:ScoredPSM{Float32, Float16},
                                            Q<:UnscoredPSM{Float32},
                                            R<:SpectralScores{Float16}}
    frag_ppm_err = Float32(getLocation(err_dist))
    println("frag_ppm_err $frag_ppm_err")
    
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
        
        frag_ppm_err = Float32(frag_ppm_err),
        unmatched_penalty_factor = params[:unmatched_penalty_factor],
        isotope_err_bounds = params[:isotope_err_bounds],
        max_peak_width = params[:max_peak_width],
        max_peaks = params[:max_peaks],
        min_topn_of_m = params[:min_topn_of_m],
        filter_by_rank = true,
        huber_δ = params[:huber_δ],
        min_frag_count = params[:min_frag_count],
        min_matched_ratio = params[:min_matched_ratio],
        min_index_search_score = zero(UInt8),#params[:min_index_search_score],
        min_spectral_contrast = params[:min_spectral_contrast], 
        min_weight = params[:min_weight],
        min_max_ppm = (5.0f0, 25.0f0),
        n_frag_isotopes = params[:n_frag_isotopes],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_index = rt_index,
        rt_bounds = params[:rt_bounds],
        rt_tol = params[:rt_tol],
        sample_rate = Inf,
        topN = params[:topN],
    )
end

function filterMatchedIons!(IDtoNMatches::ArrayDict{UInt32, UInt16}, ionMatches::Vector{FragmentMatch{Float32}}, ionMisses::Vector{FragmentMatch{Float32}}, nmatches::Int64, nmisses::Int64, max_rank::Int64, min_matched_ions::Int64)
    nmatches_all, nmisses_all = nmatches, nmisses

    for i in range(1, nmatches_all)
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

function reset!(lf::Vector{L}, last_non_empty::Int64) where {L<:LibraryFragmentIon{Float32}}
    for i in range(1, last_non_empty)
        lf[i] = L()
    end
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

function reset!(fms::Vector{M}, last_non_empty::Int64) where {M<:MatchIon{Float32}}
    for i in range(1, last_non_empty)
        fms[i] = M()
    end
end
