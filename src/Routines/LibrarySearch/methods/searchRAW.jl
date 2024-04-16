function searchRAW(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Arrow.Table, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel{Float32},
                    searchScan!::Union{Function, Missing},
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel{Float32},
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, UInt8}}};
                    #keyword args
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float32 = 0.0f0,

                    δ::Float32 = 10000f0,
                    λ::Float32 = 0f0,
                    max_iter_newton::Int64 = 100,
                    max_iter_bisection::Int64 = 100,
                    max_iter_outer::Int64 = 100,
                    accuracy_newton::Float32 = 100f0,
                    accuracy_bisection::Float32 = 100000f0,
                    max_diff::Float32 = 0.01f0,

                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                    min_frag_count::Int64 = 1,
                    min_spectral_contrast::Float32  = 0f0,
                    min_log2_matched_ratio::Float32 = -Inf32,
                    min_index_search_score::UInt8 = zero(UInt8),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    filter_by_rank::Bool = false, 
                    filter_by_count::Bool = true,
                    max_best_rank::Int64 = one(Int64),
                    n_frag_isotopes::Int64 = 1,
                    quadrupole_isolation_width::Float64 = 8.5,
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
    @time thread_tasks, total_peaks = partitionScansToThreads(spectra[:masses],
                                                        spectra[:retentionTime],
                                                         spectra[:precursorMZ],
                                                         spectra[:msOrder],
                                                        Threads.nthreads(),
                                                        1)

    if ismissing(precs)
        precs = [missing for _ in range(1, Threads.nthreads())]
    end
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return searchFragmentIndex(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                frag_index,
                                scan_to_prec_idx,
                                rt_to_irt_spline,
                                mass_err_model,
                                searchScan!,
                                frag_ppm_err,
                                precs[thread_id],
                                isotope_err_bounds,
                                min_index_search_score,
                                min_max_ppm,
                                quadrupole_isolation_width,
                                irt_tol,
                                sample_rate,
                                spec_order
                            )
        end
    end
    precursors_passed_scoring = fetch.(tasks)
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                ion_list, 
                                rt_to_irt_spline,
                                ms_file_idx,
                                mass_err_model,
                                collect_fmatches,
                                frag_ppm_err,
                                δ,
                                λ,
                                max_iter_newton,
                                max_iter_bisection,
                                max_iter_outer,
                                accuracy_newton,
                                accuracy_bisection,
                                max_diff,
                                ionMatches[thread_id],
                                ionMisses[thread_id],
                                all_fmatches[thread_id],
                                IDtoCOL[thread_id],
                                ionTemplates[thread_id],
                                scored_PSMs[thread_id],
                                unscored_PSMs[thread_id],
                                spectral_scores[thread_id], 
                                precursor_weights[thread_id],
                                precs[thread_id],
                                isotope_dict,
                                isotope_err_bounds,
                                min_frag_count,
                                min_spectral_contrast,
                                min_log2_matched_ratio,
                                min_topn_of_m,
                                min_max_ppm,
                                filter_by_rank,
                                filter_by_count,
                                max_best_rank,
                                n_frag_isotopes,
                                quadrupole_isolation_width,
                                irt_tol,
                                spec_order
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end

function searchFragmentIndex(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    rt_to_irt_spline::Any,
                    mass_err_model::MassErrorModel{Float32},
                    searchScan!::Union{Function, Missing},
                    frag_ppm_err::Float32,
                    precs::Union{Missing, Counter{UInt32, UInt8}},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_index_search_score::UInt8,
                    min_max_ppm::Tuple{Float32, Float32},
                    quadrupole_isolation_width::Float64,
                    irt_tol::Float64,
                    sample_rate::Float64,
                    spec_order::Set{Int64},
                    )

    thread_peaks = 0
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    cycle_idx = 0
    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1
    ##########
    #Iterate through spectra
    for i in thread_task
        if i == 0 
            continue
        end
        if i > length(spectra[:masses])
            continue
        end
        thread_peaks += length(spectra[:masses][i])
        if thread_peaks > 100000
            #lock(lk) do 
            #    update(pbar, thread_peaks)
            #end
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
        
        #if (i < 100000) | (i > 100000)
        #    continue
        #end
        first(rand(1)) <= sample_rate ? nothing : continue #coin flip. Usefull for random sampling of scans. 
        iRT_low, iRT_high = getRTWindow(rt_to_irt_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, irt_tol) #Convert RT to expected iRT window
        
        #Get correct rt_bin_idx 
        while getHigh(getRTBin(frag_index, rt_bin_idx)) < iRT_low
            rt_bin_idx += 1
            if rt_bin_idx >length(getRTBins(frag_index))
                rt_bin_idx = length(getRTBins(frag_index))
                break
            end 
        end

        #if i != 100000
        #    continue
        #end
        #println("rt_bin_idx $rt_bin_idx iRT_low $iRT_low getRTBin(frag_index, rt_bin_idx) ", getRTBin(frag_index, rt_bin_idx))
        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        #if !ismissing(precs)
            #searchScan! is the MsFragger style fragment index search
            #Loads `precs` with scores for each potential precursor matching the spectrum
            #println("length(spectra[:masses][i]) ", length(spectra[:masses][i]))
            #println("spec_order", spec_order)
            searchScan!(
                        precs, #counter which keeps track of plausible matches 
                        getRTBins(frag_index),
                        getFragBins(frag_index),
                        getFragments(frag_index), 
                        spectra[:masses][i], spectra[:intensities][i],
                        rt_bin_idx, 
                        iRT_high,
                        frag_ppm_err,
                        mass_err_model,
                        min_max_ppm,
                        spectra[:precursorMZ][i],
                        Float32(quadrupole_isolation_width/2.0),
                        isotope_err_bounds
                        )
            
            match_count, prec_count = filterPrecursorMatches!(precs, min_index_search_score)
            
            if getID(precs, 1)>0

                start_idx = prec_id + 1
                #stop_idx = start_idx
                n = 1
                while n <= precs.matches
                    prec_id += 1
                    if prec_id > length(precursors_passed_scoring)
                        append!(precursors_passed_scoring, 
                                Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring))
                                )
                    end
                    precursors_passed_scoring[prec_id] = getID(precs, n)
                    n += 1
                end
                scan_to_prec_idx[i] = start_idx:prec_id#stop_idx
            else
                scan_to_prec_idx[i] = missing
            end
            #if !ismissing(precs)
            reset!(precs)
            #end
            continue

    end

    return precursors_passed_scoring[1:prec_id]

end

function getPSMS(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    precursors::Union{Arrow.Table, Missing},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    precursors_passed_scoring::Vector{UInt32},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel{Float32},
                    collect_fmatches::Bool,
                    frag_ppm_err::Float32,
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
                    all_fmatches::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{L},
                    scored_PSMs::Vector{S},
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    precs::Union{Missing, Counter{UInt32, UInt8}},

                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_frag_count::Int64,
                    min_spectral_contrast::Float32,
                    min_log2_matched_ratio::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    min_max_ppm::Tuple{Float32, Float32},
                    filter_by_rank::Bool,
                    filter_by_count::Bool,
                    max_best_rank::Int64,
                    n_frag_isotopes::Int64,
                    quadrupole_isolation_width::Float64,
                    irt_tol::Float64,
                    spec_order::Set{Int64},
                    ) where {L<:LibraryIon{Float32},
                    S<:ScoredPSM{Float32, Float16},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    frag_err_idx = 1
    prec_idx, ion_idx, cycle_idx, last_val = 0, 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    #prec_id = 0
    #precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    ##########
    #Iterate through spectra
    for i in thread_task
        if i == 0 
            continue
        end
        if i > length(spectra[:masses])
            continue
        end
        if ismissing(scan_to_prec_idx[i])
            continue
        end

        ###########
        #Scan Filtering
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        if msn == 1
            cycle_idx += 1
        end
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type
        
        #if (i < 100000) | (i > 100000)
        #    continue
        #end
        iRT_low, iRT_high = getRTWindow(rt_to_irt_spline(spectra[:retentionTime][i])::Union{Float64,Float32}, irt_tol) #Convert RT to expected iRT window
        
        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        #Candidate precursors and their retention time estimates have already been determined from
        #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
        #the retention time and m/z tolerance constraints
        ion_idx, prec_idx = selectTransitions!(
                                        ionTemplates,
                                        scan_to_prec_idx[i],
                                        precursors_passed_scoring,
                                        precursors[:mz],
                                        precursors[:prec_charge],
                                        precursors[:irt],
                                        ion_list,
                                        Float32(rt_to_irt_spline(spectra[:retentionTime][i])),
                                        Float32(irt_tol), #rt_tol
                                        (
                                        spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                        spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0)
                                        ),
                                        isotope_err_bounds = isotope_err_bounds
                                        )
        if ion_idx < 2
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
                                        mass_err_model,
                                        min_max_ppm,
                                        UInt32(i), 
                                        ms_file_idx)
        #println("nmatches $nmatches nmisses $nmisses")
        nmisses_all = nmisses
        nmatches_all = nmatches    
 
        if filter_by_count
            nmatches_all, nmisses_all, nmatches, nmisses = filterMatchedIons!(IDtoCOL, 
                                                                                ionMatches, 
                                                                                ionMisses, 
                                                                                nmatches, 
                                                                                nmisses,
                                                                                min_frag_count #Remove precursors matching fewer than this many fragments
                                                                            )
        end                                                       
        if filter_by_rank
        #println("nmatches_all $nmatches_all, nmatches $nmatches")
            _, _, nmatches, nmisses = filterMatchedIonsTop!(IDtoCOL, 
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
                solveHuber!(Hs, _residuals_, _weights_, δ, λ, 
                                max_iter_newton, 
                                max_iter_bisection,
                                max_iter_outer,
                                accuracy_newton,
                                accuracy_bisection,
                                Hs.n,
                                max_diff
                                );
            end
            #Remember determined weights for each precursors
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
            end

            if ismissing(isotope_dict) 
                getDistanceMetrics(_weights_, Hs, spectral_scores)
            end

            ##########
            #Scoring and recording data
            #if false==true#!ismissing(scored_PSMs)
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
                    max_best_rank = max_best_rank,
                    min_topn = first(min_topn_of_m),
                    block_size = 500000,
                    )
            #end
        end
        #Add fragment matches to all_fmatches 
        frag_err_idx = collectFragErrs(all_fmatches, ionMatches, nmatches, frag_err_idx, collect_fmatches)
    
        ##########
        #Reset pre-allocated arrays 
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        reset!(IDtoCOL);
        reset!(Hs);
    end
    #println("prec_idx $prec_idx")
    #println("scan_to_prec_idx $scan_to_prec_idx")

    if collect_fmatches
        return DataFrame(@view(scored_PSMs[1:last_val])), @view(all_fmatches[1:frag_err_idx])
    else
        return DataFrame(@view(scored_PSMs[1:last_val]))
    end
end

function quantPSMs(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    precursors::Union{Arrow.Table, Missing},
                    library_fragment_lookup::Union{LibraryFragmentLookup{Float32}, Missing},
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel{Float32},
                    frag_ppm_err::Float32,
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
                    iso_splines::IsotopeSplineModel{Float32},
                    scored_PSMs::Vector{S},
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_frag_count::Int64,
                    min_spectral_contrast::Float32,
                    min_log2_matched_ratio::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    min_max_ppm::Tuple{Float32, Float32},
                    max_best_rank::Int64,
                    n_frag_isotopes::Int64,
                    quadrupole_isolation_width::Float64,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    irt_tol::Float64,
                    spec_order::Set{Int64},
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
    isotopes = zeros(Float32, n_frag_isotopes)
    rt_start, rt_stop = 1, 1
    #test_n = 0
    #n_templates = zeros(Int64, length(thread_task))
    n = 0
    ##########
    #Iterate through spectra
    for i in thread_task
        n += 1
        if i == 0 
            continue
        end
        if i > length(spectra[:masses])
            continue
        end
        ###########
        #Scan Filtering
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        cycle_idx += (msn == 1)
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        #a = (getHigh(getRTBin(rt_index, rt_bin_high)) + irt_tol < spectra[:retentionTime][i]) 
        #b = (getLow(getRTBin(rt_index, rt_bin_high)) - irt_tol > spectra[:retentionTime][i]) 
        rt = spectra[:retentionTime][i]
        rt_start_new = max(searchsortedfirst(rt_index.rt_bins, rt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
        rt_stop_new = min(searchsortedlast(rt_index.rt_bins, rt + irt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    
        
        if (rt_start_new != rt_start) | rt_stop_new != rt_stop
            rt_start = rt_start_new
            rt_stop = rt_stop_new
            #Candidate precursors and their retention time estimates have already been determined from
            #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
            #the retention time and m/z tolerance constraints
            #test_n += 1
            ion_idx, prec_idx = selectRTIndexedTransitions!(
                                            ionTemplates,
                                            library_fragment_lookup,
                                            precursors[:mz],
                                            precursors[:prec_charge],
                                            precursors[:sulfur_count],
                                            iso_splines,
                                            isotopes,
                                            rt_index,
                                            rt_start,
                                            rt_stop,
                                            spectra[:precursorMZ][i] - Float32(quadrupole_isolation_width/2.0),
                                            spectra[:precursorMZ][i] + Float32(quadrupole_isolation_width/2.0),
                                            isotope_err_bounds,
                                            10000)
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
                                        mass_err_model,
                                        min_max_ppm,
                                        UInt32(i), 
                                        ms_file_idx)
        #println([x for x in ionMisses[1:nmisses] if (x.prec_id==0x00a2e2fa)])
        #println("any zeros yet ?", any([iszero(getPredictedIntenisty(x)) for x in ionMisses[1:nmisses]]))
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
                            Hs.n,
                            max_diff
                            );
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
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(spectra[:intensities][i])), 
                i,
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                max_best_rank = max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000,
                )
        end

        ##########
        #Reset pre-allocated arrays 
        #n_templates[n] = Hs.n
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        #reset!(ionTemplates, ion_idx)
        reset!(IDtoCOL);
        reset!(Hs);
    end
    #println("describe(n_templates) ", describe(n_templates))
    #println("test_n $test_n")
    return DataFrame(@view(scored_PSMs[1:last_val]))
end

function filterMatchedIonsTop!(IDtoNMatches::ArrayDict{UInt32, UInt16}, ionMatches::Vector{FragmentMatch{Float32}}, ionMisses::Vector{FragmentMatch{Float32}}, nmatches::Int64, nmisses::Int64, max_rank::Int64, min_matched_ions::Int64)
    nmatches_all, nmisses_all = nmatches, nmisses

    @inbounds for i in range(1, nmatches_all)
        match = ionMatches[i]
        prec_id = getPrecID(match)
        if match.is_isotope 
            continue
        end
            if (getRank(match) <= max_rank) .& (match.is_isotope==false)
                #if getIonType(match) == 'y' #| (filter_y==false)
                    if iszero(IDtoNMatches[prec_id])
                        update!(IDtoNMatches, prec_id, one(UInt16))
                    else
                        IDtoNMatches.vals[prec_id] += one(UInt16)
                    end
                #end
            end
    end
    nmatches, nmisses = 0, 0
    @inbounds for i in range(1, nmatches_all)
        if IDtoNMatches[getPrecID(ionMatches[i])] <= min_matched_ions

            continue
        else
            nmatches += 1
            ionMatches[nmatches] = ionMatches[i]
        end
    end

    @inbounds for i in range(1, nmisses_all)
        if IDtoNMatches[getPrecID(ionMisses[i])] <= min_matched_ions
            continue
        else
            nmisses += 1
            ionMisses[nmisses] = ionMisses[i]
        end
    end

    reset!(IDtoNMatches)

    return nmatches_all, nmisses_all, nmatches, nmisses
end

function filterMatchedIons!(IDtoNMatches::ArrayDict{UInt32, UInt16}, 
                                ionMatches::Vector{FragmentMatch{Float32}}, 
                                ionMisses::Vector{FragmentMatch{Float32}}, 
                                nmatches::Int64, 
                                nmisses::Int64, 
                                min_matched_ions::Int64
                                )
    nmatches_all, nmisses_all = nmatches, nmisses

    @inbounds for i in range(1, nmatches_all)
        match = ionMatches[i]
        prec_id = getPrecID(match)
        if match.is_isotope 
            continue
        end
            #if getRank(match) <= max_rank
                #if ((getIonType(match) == 'y') | (filter_on_y==false)) & ((match.is_isotope==false) | (filter_on_mono == false))
                    if iszero(IDtoNMatches[prec_id])
                        update!(IDtoNMatches, prec_id, one(UInt16))
                    else
                        IDtoNMatches.vals[prec_id] += one(UInt16)
                    end
            #end
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


function buildSubsetLibraryFragmentLookupTable!(
                                            precursors_passed_scoring::Vector{Vector{UInt32}},
                                            lookup_table::LibraryFragmentLookup{Float32},
                                            n_precursors::Int64
                                            )
    #prec_id_conversion = zeros(UInt32, n_precursors)
    @time begin
        precursors = sort(unique(vcat(precursors_passed_scoring...)))
        prec_frag_ranges = Vector{UnitRange{UInt32}}(undef, n_precursors)
    end
    fragment_count = 0
    n = 0
    @time for prec_idx in precursors
        n += 1
        prec_size = length(getPrecFragRange(lookup_table, prec_idx))
        prec_frag_ranges[prec_idx] = UInt32(fragment_count + 1):UInt32(fragment_count + prec_size)
        fragment_count += prec_size
    end
    println("fragment_count $fragment_count")
    fragments = Vector{DetailedFrag{Float32}}(undef, fragment_count)
    n = 0
    @time for prec_idx in precursors
        for i in getPrecFragRange(lookup_table, prec_idx)
            n+=1
            fragments[n] = lookup_table.frags[i]
        end
    end

    return LibraryFragmentLookup(fragments, prec_frag_ranges)
end

