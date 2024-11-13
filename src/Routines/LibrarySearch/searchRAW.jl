function searchFragmentIndex(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    rt_to_irt_spline::Any,
                    mass_err_model::MassErrorModel,
                    searchScan!::Union{Function, Missing},
                    precs::Union{Missing, Counter{UInt32, UInt8}},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    quad_transmission_model::QuadTransmissionModel,
                    min_index_search_score::UInt8,
                    irt_tol::Float64,
                    sample_rate::Float64,
                    spec_order::Set{Int64},
                    )
    thread_peaks = 0
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1
    ##########
    #Iterate through spectra
    for i in thread_task
        if i == 0 
            continue
        end
        if i > length(spectra[:mz_array])
            continue
        end
        thread_peaks += length(spectra[:mz_array][i])
        if thread_peaks > 100000

            thread_peaks = 0
        end
        ###########
        #Scan Filtering
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type
        
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
        ##########
        #Ion Template Selection
        #searchScan! is the MsFragger style fragment index search
        #Loads `precs` with scores for each potential precursor matching the spectrum
        quad_transmission_function = getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][i], spectra[:isolationWidthMz][i])
        searchScan!(
                    precs, #counter which keeps track of plausible matches 
                    getRTBins(frag_index),
                    getFragBins(frag_index),
                    getFragments(frag_index), 
                    spectra[:mz_array][i], spectra[:intensity_array][i],
                    rt_bin_idx, 
                    iRT_high,
                    mass_err_model,
                    quad_transmission_function,
                    isotope_err_bounds
                    )
        
        match_count, prec_count = filterPrecursorMatches!(precs, min_index_search_score)
        
        if getID(precs, 1)>0
            start_idx = prec_id + 1
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
        reset!(precs)
        continue
    end
    return precursors_passed_scoring[1:prec_id]
end
function getPSMS(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},
                    precursors::Union{Arrow.Table, Missing},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    precursors_passed_scoring::Vector{UInt32},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    quad_transmission_model::QuadTransmissionModel,
                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    IDtoCOL::ArrayDict{UInt32, UInt16},
                    ionTemplates::Vector{L},
                    iso_splines::IsotopeSplineModel,
                    scored_PSMs::Vector{<:ScoredPSM},
                    unscored_PSMs::Vector{<:UnscoredPSM},
                    spectral_scores::Vector{<:SpectralScores},
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_frag_count::Int64,
                    min_spectral_contrast::AbstractFloat,
                    min_log2_matched_ratio::AbstractFloat,
                    min_topn_of_m::Tuple{Int64, Int64},
                    max_best_rank::Int64,
                    n_frag_isotopes::Int64,
                    irt_tol::AbstractFloat,
                    spec_order::Set{Int64},
                    ) where {L<:LibraryIon{Float32}}
    ##########
    #Initialize 
    msms_counts = Dict{Int64, Int64}()
    prec_idx, ion_idx, last_val = 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    isotopes = zeros(Float32, 5)#n_frag_isotopes);
    precursor_transmission = zeros(Float32, 5)
    ##########
    #Iterate through spectra
    for i in thread_task
        if i == 0 
            continue
        end
        if i > length(spectra[:mz_array])
            continue
        end
        if ismissing(scan_to_prec_idx[i])
            continue
        end
        ###########
        #Scan Filtering
        msn = spectra[:msOrder][i] #An integer 1, 2, 3.. for MS1, MS2, MS3 ...
        msn ∈ spec_order ? nothing : continue #Skip scans outside spec order. (Skips non-MS2 scans is spec_order = Set(2))
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 #Update counter for each MSN scan type
        ##########
        #Ion Template Selection
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
                                        precursors[:sulfur_count],
                                        iso_splines,
                                        getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][i], spectra[:isolationWidthMz][i]),  
                                        precursor_transmission,
                                        isotopes,
                                        n_frag_isotopes,
                                        ion_list,
                                        Float32(rt_to_irt_spline(spectra[:retentionTime][i])),
                                        Float32(irt_tol),                                  
                                        ( 
                                            spectra[:lowMz][i], spectra[:highMz][i]
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
                                        spectra[:mz_array][i], 
                                        spectra[:intensity_array][i], 
                                        mass_err_model,
                                        spectra[:highMz][i],
                                        UInt32(i), 
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
            if IDtoCOL.size > length(spectral_scores)
                new_entries = IDtoCOL.size - length(spectral_scores) + 1000 
                append!(spectral_scores, Vector{eltype(spectral_scores)}(undef, new_entries))
                append!(unscored_PSMs, [eltype(unscored_PSMs)() for _ in 1:new_entries]);
            end 

            getDistanceMetrics(Hs, spectral_scores)

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
                IDtoCOL,
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(spectra[:intensity_array][i])), 
                i,
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                max_best_rank = max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000,
                )
        end
        #Add fragment matches to all_fmatches 
        ##########
        #Reset pre-allocated arrays 
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        reset!(IDtoCOL);
        reset!(Hs);
    end
    return DataFrame(@view(scored_PSMs[1:last_val]))
end
function getMassErrors(
                    spectra::Arrow.Table,
                    library_fragment_lookup::LibraryFragmentLookup{Float32},
                    thread_task::UnitRange{Int64},
                    scan_idxs::Vector{UInt32},
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    precursors_passed_scoring::Vector{UInt32},
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    ionMatches::Vector{FragmentMatch{Float32}},
                    ionMisses::Vector{FragmentMatch{Float32}},
                    all_fmatches::Vector{FragmentMatch{Float32}},
                    ionTemplates::Vector{L};
                    max_rank::Int64 = 5
                    ) where {L<:LibraryIon{Float32}}
    ##########
    #Initialize 
    frag_err_idx = 0
    prec_idx, ion_idx = 0, 0
    #prec_id = 0
    #precursors_passed_scoring = Vector{UInt32}(undef, 250000)                                                                                                      
    ##########
    #Iterate through spectra
    for n in thread_task
        i = scan_idxs[n]
        if i == 0 
            continue
        end
        if i > length(spectra[:mz_array])
            continue
        end
        if ismissing(scan_to_prec_idx[i])
            continue
        end
        ##########
        #Ion Template Selection
        #SearchScan! applies a fragment-index based search (MS-Fragger) to quickly identify high-probability candidates to explain the spectrum  

        #Candidate precursors and their retention time estimates have already been determined from
        #A previous serach and are incoded in the `rt_index`. Add candidate precursors that fall within
        #the retention time and m/z tolerance constraints
        ion_idx, prec_idx = selectTransitions!(
                                        ionTemplates,
                                        library_fragment_lookup,
                                        scan_to_prec_idx[i],
                                        precursors_passed_scoring,
                                        max_rank = max_rank
                                        )
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:mz_array][i], 
                                        spectra[:intensity_array][i], 
                                        mass_err_model,
                                        spectra[:highMz][i],
                                        UInt32(i), 
                                        ms_file_idx)  
        #sort!(@view(ionMatches[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        #Add fragment matches to all_fmatches 
        frag_err_idx = collectFragErrs(all_fmatches, ionMatches, nmatches, frag_err_idx)
    end
    return @view(all_fmatches[1:frag_err_idx])
end


function huberTuningSearch(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    prec_set::Set{Tuple{UInt32,UInt32}},
                    scan_idxs::Set{UInt32},
                    precursors::Union{Arrow.Table, Missing},
                    library_fragment_lookup::Union{LibraryFragmentLookup{Float32}, Missing},
                    ms_file_idx::UInt32,
                    rt_to_irt::UniformSpline,
                    mass_err_model::MassErrorModel,
                    δs::Vector{Float32},
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
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    quad_transmission_model::QuadTransmissionModel,
                    isotope_err_bounds::Tuple{Int64, Int64},
                    n_frag_isotopes::Int64,
                    rt_index::retentionTimeIndex{Float32, Float32},
                    irt_tol::Float32,
                    spec_order::Set{Int64}
                    ) where {L<:LibraryIon{Float32},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}

    ##########
    #Initialize 
    prec_idx, ion_idx, cycle_idx = 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);

    tuning_results = Dict(
        :precursor_idx => UInt32[],
        :scan_idx => UInt32[],
        :weight => Float32[],
        :huber_δ => Float32[]
    )

    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""

    rt_idx = 0
    prec_temp_size = 0
    precs_temp = Vector{UInt32}(undef, 50000)

    
    ##########
    #Iterate through spectra
    for scan_idx in thread_task
        if scan_idx ∉ scan_idxs
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
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]
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
                getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),  
                precursor_transmission,
                isotopes,
                n_frag_isotopes,
                rt_index,
                irt_start,
                irt_stop,
                (
                    spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]
                ),
                isotope_err_bounds,
                10000)
        end
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:mz_array][scan_idx], 
                                        spectra[:intensity_array][scan_idx], 
                                        mass_err_model,
                                        spectra[:highMz][scan_idx],
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

            for δ in δs
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
                #fill!(_residuals_, zero(Float32))
                #fill!(_weights_, zero(Float32))
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
                #Record weights for each precursor
                for i in range(1, IDtoCOL.size)
                    precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
                    pid = IDtoCOL.keys[i]
                    if (pid,UInt32(scan_idx)) ∈ prec_set
                        weight = _weights_[IDtoCOL[IDtoCOL.keys[i]]]
                        push!(tuning_results[:precursor_idx],pid)
                        push!(tuning_results[:weight],weight)
                        push!(tuning_results[:huber_δ], δ)
                        push!(tuning_results[:scan_idx], scan_idx)
                    end
                end
            end
            reset!(IDtoCOL);
        end
        ##########
        #Reset pre-allocated arrays 
        reset!(IDtoCOL);
        reset!(Hs);
    end
    return DataFrame(tuning_results)
end

function QuadTransmissionSearch(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},
                    scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}},
                    scan_idxs::Set{UInt32},
                    precursors::Union{Arrow.Table, Missing},
                    library_fragment_lookup::Union{LibraryFragmentLookup{Float32}, Missing},
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    δ::Float32,
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
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    spec_order::Set{Int64}
                    ) where {L<:LibraryIon{Float32},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}

    ##########
    #Initialize 
    prec_idx, ion_idx, cycle_idx = 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);

    tuning_results = Dict(
        :precursor_idx => UInt32[],
        :scan_idx => UInt32[],
        :weight => Float32[],
        :iso_idx => UInt8[],
        :center_mz => Float32[],
        :n_matches => UInt8[]
    )

    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    ##########
    #Iterate through spectra
    for scan_idx in thread_task
        if scan_idx ∉ scan_idxs
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
        ion_idx, prec_idx = selectTransitionsForQuadEstimation!(
            scan_idx_to_prec_idx[scan_idx],
            ionTemplates,
            library_fragment_lookup,
            precursors[:mz],
            precursors[:prec_charge],
            precursors[:sulfur_count],
            iso_splines,
            precursor_transmission,
            isotopes,
            (
                spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]
            ),
            10000)
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:mz_array][scan_idx], 
                                        spectra[:intensity_array][scan_idx], 
                                        mass_err_model,
                                        spectra[:highMz][scan_idx],
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
                            δ, 0.0f0, 
                            max_iter_newton, 
                            max_iter_bisection,
                            max_iter_outer,
                            accuracy_newton,
                            accuracy_bisection,
                            10.0,
                            max_diff
                            );
            #Record weights for each precursor
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
                id = IDtoCOL.keys[i]
                isotope_idx = UInt32(((id - 1) % 3) + 1)
                pid = UInt32(((id - 1) ÷ 3) + 1)
                    colid = IDtoCOL[IDtoCOL.keys[i]]
                    weight = _weights_[colid]
                    n_matches = 0
                    for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1)
                        n_matches += Hs.matched[j]
                    end
                    push!(tuning_results[:precursor_idx],pid)
                    push!(tuning_results[:weight],weight)
                    push!(tuning_results[:iso_idx], isotope_idx)
                    push!(tuning_results[:scan_idx], scan_idx)
                    push!(tuning_results[:center_mz], spectra[:centerMz][scan_idx])
                    push!(tuning_results[:n_matches], n_matches)                end
            end
        ##########
        #Reset pre-allocated arrays 
        reset!(IDtoCOL);
        reset!(Hs);
    end
    return DataFrame(tuning_results)
end
function secondSearch(
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
                    scored_PSMs::Vector{S},
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    quad_transmission_model::QuadTransmissionModel,
                    isotope_err_bounds::Tuple{Int64, Int64},
                    min_y_count::Int64,
                    min_frag_count::Int64,
                    min_spectral_contrast::Float32,
                    min_log2_matched_ratio::Float32,
                    min_topn_of_m::Tuple{Int64, Int64},
                    max_best_rank::Int64,
                    n_frag_isotopes::Int64,
                    rt_index::retentionTimeIndex{Float32, Float32},
                    irt_tol::Float32,
                    spec_order::Set{Int64}
                    ) where {L<:LibraryIon{Float32},
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
    precursor_transmission = zeros(Float32, 5)
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
        if scan_idx > length(spectra[:mz_array])
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
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]

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
                getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),  
                precursor_transmission,
                isotopes,
                n_frag_isotopes,
                rt_index,
                irt_start,
                irt_stop, 
                (
                    spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]
                ),
                isotope_err_bounds,
                10000)
        end
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:mz_array][scan_idx], 
                                        spectra[:intensity_array][scan_idx], 
                                        mass_err_model,
                                        spectra[:highMz][scan_idx],
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
            #Record weights for each precursor
            for i in range(1, IDtoCOL.size)
                precursor_weights[IDtoCOL.keys[i]] = _weights_[IDtoCOL[IDtoCOL.keys[i]]]# = precursor_weights[id]
            end
            getDistanceMetrics(_weights_, _residuals_, Hs, spectral_scores)
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
                Float32(sum(spectra[:intensity_array][scan_idx])), 
                scan_idx,
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_y_count = min_y_count,
                min_frag_count = min_frag_count, #Remove precursors with fewer fragments 
                max_best_rank = max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000,
                )
        end
        ##########
        #Reset pre-allocated arrays 
        for i in range(1, Hs.n)
            unscored_PSMs[i] = eltype(unscored_PSMs)()
        end
        reset!(IDtoCOL);
        reset!(Hs);
    end
    return DataFrame(@view(scored_PSMs[1:last_val]))
end

function getChromatograms(
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},#UnitRange{Int64},
                    precursors::Union{Arrow.Table, Missing},
                    precursors_passing::Set{UInt32},
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
                    unscored_PSMs::Vector{Q},
                    spectral_scores::Vector{R},
                    precursor_weights::Vector{Float32},
                    quad_transmission_model::QuadTransmissionModel,
                    isotope_err_bounds::Tuple{Int64, Int64},
                    n_frag_isotopes::Int64,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing},
                    irt_tol::Float32,
                    spec_order::Set{Int64}
                    ) where {T,U<:AbstractFloat, L<:LibraryIon{Float32},
                    Q<:UnscoredPSM{Float32},
                    R<:SpectralScores{Float16}}

    ##########
    #Initialize 
    prec_idx, ion_idx, cycle_idx = 0, 0, 0
    Hs = SparseArray(UInt32(5000));
    _weights_ = zeros(Float32, 5000);
    _residuals_ = zeros(Float32, 5000);
    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    reached_max_iters = 0
    rt_idx = 0
    prec_temp_size = 0
    precs_temp = Vector{UInt32}(undef, 50000)
    test_vals = Float32[]
    ##########
    #Iterate through spectra
    for scan_idx in thread_task
        if scan_idx == 0 
            continue
        end
        if scan_idx > length(spectra[:mz_array])
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
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
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
                precursors_passing,
                library_fragment_lookup,
                precursors[:mz],
                precursors[:prec_charge],
                precursors[:sulfur_count],
                iso_splines,
                getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),  
                precursor_transmission,
                isotopes,
                n_frag_isotopes,
                rt_index,
                irt_start,
                irt_stop,
                (
                    spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]
                ),
                isotope_err_bounds,
                10000)
        end
        ##########
        #Match sorted list of plausible ions to the observed spectra
        nmatches, nmisses = matchPeaks!(ionMatches, 
                                        ionMisses, 
                                        ionTemplates, 
                                        ion_idx, 
                                        spectra[:mz_array][scan_idx], 
                                        spectra[:intensity_array][scan_idx], 
                                        mass_err_model,
                                        spectra[:highMz][scan_idx],
                                        UInt32(scan_idx), 
                                        ms_file_idx)
        sort!(@view(ionMatches[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        #if scan_idx == 152308
        #   return ionMatches[1:nmatches], ionMisses[1:nmisses], ionTemplates[1:ion_idx]
        #end
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
            #_δ_ = (sort(abs.(_residuals_[1:Hs.m]), rev=true)[min(Hs.n, Hs.m)])
            #_δ_ = Float32(100000.0f0)#Float32(iqr(abs.(_residuals_[1:Hs.m]))*1.5)
            #Spectral deconvolution. Hybrid bisection/newtowns method
            #_δ_ = 100000.0f0
            #for _index_ in range(1, 2)
                solveHuber!(Hs, _residuals_, _weights_, 
                                δ, 
                                #_δ_,
                                λ, 
                                max_iter_newton, 
                                max_iter_bisection,
                                max_iter_outer,
                                accuracy_newton,
                                accuracy_bisection,
                                10.0,#Hs.n/10.0,
                                max_diff
                                );
                #_δ_ = Float32(quantile(abs.(_residuals_[1:Hs.m]), 0.9))
            #    _δ_ = maximum(_weights_[1:Hs.n])/100.0f0
                #_δ_ = 10.0f0
            #end
                #_δ_ = Float32(iqr(abs.(_residuals_[1:Hs.m]))*3)
                            #Spectral deconvolution. Hybrid bisection/newtowns method
            #push!(test_vals, _δ_)
            #end
            
            #=
            iters = 0
            for _δ_ in [1500.0f0]
                iters = solveHuber!(Hs, _residuals_, _weights_, 
                _δ_, λ, 
                max_iter_newton, 
                max_iter_bisection,
                max_iter_outer,
                accuracy_newton,
                accuracy_bisection,
                10.0,#Hs.n/10.0,
                max_diff
                );
                if iters >= max_iter_outer
                    reached_max_iters += 1
                end
            end
            =#
            #push!(reached_max_iters, iters)
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
            ##########
            #Scoring and recording data
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
    return DataFrame(@view(chromatograms[1:rt_idx]))
end

function massErrorSearch(
                    #Mandatory Args
                    spectra::Arrow.Table,
                    scan_idxs::Vector{UInt32},
                    precursors_passed_scoring::Vector{UInt32},
                    library_fragment_lookup::LibraryFragmentLookup{Float32},
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionTemplates::Vector{Vector{L}}
                    ) where {L<:LibraryIon{Float32}}

    #Sort scans and precursors 
    sorted_scan_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_scan_indices]
    precursors_passed_scoring = precursors_passed_scoring[sorted_scan_indices]
    function getScanToPrecIdx(scan_idxs::Vector{UInt32}, n_scans::Int64)
        scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_scans)
        start_idx, stop_idx = 1, 1
        for i in range(1, length(scan_idxs))
            stop_idx = i 
            if scan_idxs[start_idx] == scan_idxs[stop_idx]
                scan_to_prec_idx[scan_idxs[i]] = start_idx:stop_idx
            else
                scan_to_prec_idx[scan_idxs[i]] = i:i
                start_idx = i
            end
        end
        return scan_to_prec_idx
    end

    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra[:mz_array]))
  
    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    
    #Build thread tasks 
    thread_task_size = length(scan_idxs)÷Threads.nthreads()
    thread_tasks = []
    start_idx, stop_idx =1, min(thread_task_size - 1, length(scan_idxs))
    for i in range(1, Threads.nthreads())
        push!(thread_tasks, (i, start_idx:stop_idx))
        start_idx = stop_idx + 1
        stop_idx =  min(start_idx + thread_task_size - 1, length(scan_idxs))
    end

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getMassErrors(
                                spectra,
                                library_fragment_lookup,
                                last(thread_task),
                                scan_idxs,
                                scan_to_prec_idx,
                                precursors_passed_scoring,
                                ms_file_idx,
                                mass_err_model,
                                ionMatches[thread_id],
                                ionMisses[thread_id],
                                all_fmatches[thread_id],
                                ionTemplates[thread_id]
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end
function LibrarySearch(
    spectra::Arrow.Table,
    params::NamedTuple;
    kwargs...)
    #println("TEST $sample_rate masserrmodel ", kwargs[:mass_err_model])
    thread_tasks, total_peaks = partitionScansToThreads(
                                                            spectra[:mz_array],
                                                            spectra[:retentionTime],
                                                            spectra[:centerMz],
                                                            spectra[:msOrder],
                                                            Threads.nthreads(),
                                                            1
                                                        )

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))
    #println("start frag index search...")
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return searchFragmentIndex(
                                spectra,
                                last(thread_task),
                                kwargs[:frag_index],
                                scan_to_prec_idx,
                                kwargs[:rt_to_irt_spline],
                                kwargs[:mass_err_model],
                                searchScan!,
                                kwargs[:prec_to_score][thread_id],
                                Tuple([Int64(x) for x in kwargs[:isotope_err_bounds]]),
                                kwargs[:quad_transmission_model],
                                UInt8(kwargs[:params]["min_index_search_score"]),
                                Float64(kwargs[:irt_tol]),
                                kwargs[:sample_rate],
                                Set(2)
                            )
        end
    end
    precursors_passed_scoring = fetch.(tasks)
    #println("Finished frag index search...")
    #println("start psms thing...")
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                kwargs[:precursors],
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                kwargs[:fragment_lookup_table], 
                                kwargs[:rt_to_irt_spline],
                                kwargs[:ms_file_idx],
                                kwargs[:mass_err_model],
                                kwargs[:quad_transmission_model],
                                kwargs[:ion_matches][thread_id],
                                kwargs[:ion_misses][thread_id],
                                kwargs[:id_to_col][thread_id],
                                kwargs[:ion_templates][thread_id],
                                kwargs[:iso_splines],
                                kwargs[:scored_psms][thread_id],
                                kwargs[:unscored_psms][thread_id],
                                kwargs[:spectral_scores][thread_id],
                                Tuple([Int64(x) for x in kwargs[:isotope_err_bounds]]),
                                kwargs[:params]["min_frag_count"],
                                Float32(kwargs[:params]["min_spectral_contrast"]),
                                Float32(kwargs[:params]["min_log2_matched_ratio"]),
                                Tuple([Int64(x) for x in kwargs[:params]["min_topn_of_m"]]),
                                kwargs[:params]["max_best_rank"],
                                Int64(kwargs[:params]["n_frag_isotopes"]),
                                Float32(kwargs[:irt_tol]),
                                Set(2)
                            )
        end
    end
    tasks_out = fetch.(tasks)
    return tasks_out
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

function collectFragErrs(all_fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int) where {M<:MatchIon{Float32}}
    for match in range(1, nmatches)
        if n < length(all_fmatches)
            n += 1
            all_fmatches[n] = new_fmatches[match]
        else
            all_fmatches = append!(all_fmatches, [M() for x in range(1, length(all_fmatches))])
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
    #println("fragment_count $fragment_count")
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

