function searchFragmentIndex(
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    frag_index::FragmentIndex{Float32},
                    spectra::Arrow.Table,
                    thread_task::Vector{Int64},
                    search_data::S,
                    params::P,
                    qtm::Q,
                    mem::M,
                    rt_to_irt_spline::Any,
                    irt_tol::AbstractFloat,
                    ) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1
    for scan_idx in thread_task
        # Skip invalid indices
        (scan_idx <= 0 || scan_idx > length(spectra[:mz_array])) && continue
        (spectra[:msOrder][scan_idx] ∉ getSpecOrder(params) || rand() > getSampleRate(params)) && continue

        # Update RT bin index based on iRT window
        irt_lo, irt_hi = getRTWindow(rt_to_irt_spline(spectra[:retentionTime][scan_idx]), irt_tol)

        #=
        while rt_bin_idx <= length(getRTBins(frag_index)) && getHigh(getRTBin(frag_index, rt_bin_idx)) < irt_lo
            rt_bin_idx += 1
        end
        rt_bin_idx = min(rt_bin_idx, length(getRTBins(frag_index)))
        =#
        while getHigh(getRTBin(frag_index, rt_bin_idx)) < irt_lo
            rt_bin_idx += 1
            if rt_bin_idx >length(getRTBins(frag_index))
                rt_bin_idx = length(getRTBins(frag_index))
                break
            end 
        end
        # Fragment index search for matching precursors
        searchScan!(
            getPrecursorScores(search_data),
            getRTBins(frag_index),
            getFragBins(frag_index),
            getFragments(frag_index),
            spectra[:mz_array][scan_idx],
            spectra[:intensity_array][scan_idx],
            rt_bin_idx,
            irt_hi,
            mem,
            getQuadTransmissionFunction(qtm, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),
            getIsotopeErrBounds(params)
        )
        
        # Filter precursor matches based on score
        match_count, prec_count = filterPrecursorMatches!(getPrecursorScores(search_data), getMinIndexSearchScore(params))
        #=
        if getID(getPrecursorScores(search_data), 1) > 0
            start_idx = prec_id + 1
            for n in 1:getPrecursorScores(search_data).matches
                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring, similar(precursors_passed_scoring))
                end
                precursors_passed_scoring[prec_id] = getID(getPrecursorScores(search_data), n)
            end
            scan_to_prec_idx[scan_idx] = start_idx:prec_id
        else
            scan_to_prec_idx[scan_idx] = missing
        end
        
        reset!(getPrecursorScores(search_data))
        =#
        if getID(getPrecursorScores(search_data), 1)>0
            start_idx = prec_id + 1
            n = 1
            while n <= getPrecursorScores(search_data).matches
                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring, 
                            Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring))
                            )
                end
                precursors_passed_scoring[prec_id] = getID(getPrecursorScores(search_data), n)
                n += 1
            end
            scan_to_prec_idx[scan_idx] = start_idx:prec_id#stop_idx
        else
            scan_to_prec_idx[scan_idx] = missing
        end

        reset!(getPrecursorScores(search_data))
    end

    return precursors_passed_scoring[1:prec_id]
end

function getPSMS(
    ms_file_idx::UInt32,
    spectra::Arrow.Table,
    thread_task::Vector{Int64},
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup,
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}}, 
    precursors_passed_scoring::Vector{UInt32},
    search_data::S,
    params::P,
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    irt_tol::AbstractFloat) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:SearchParameters}

    msms_counts = Dict{Int64, Int64}()
    last_val = 0
    Hs = SparseArray(UInt32(5000))
    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    for scan_idx in thread_task
        (scan_idx == 0 || scan_idx > length(spectra[:mz_array])) && continue
        ismissing(scan_to_prec_idx[scan_idx]) && continue

        # Scan Filtering
        msn = spectra[:msOrder][scan_idx]
        msn ∉ getSpecOrder(params) && continue
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1

        # Ion Template Selection
        ion_idx = selectTransitions!(
            getIonTemplates(search_data), 
            StandardTransitionSelection(), 
            getPrecEstimation(params),
            ion_list,
            scan_to_prec_idx[scan_idx], precursors_passed_scoring,
            precursors[:mz], 
            precursors[:prec_charge], 
            precursors[:sulfur_count],
            precursors[:irt],
            getIsoSplines(search_data),
            getQuadTransmissionFunction(qtm, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),
            precursor_transmission, isotopes, getNFragIsotopes(params),
            getMaxFragRank(params),
            Float32(rt_to_irt_spline(spectra[:retentionTime][scan_idx])),
            Float32(irt_tol), 
            (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        ion_idx < 2 && continue


        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data), 
            getIonMisses(search_data), 
            getIonTemplates(search_data), 
            ion_idx, 
            spectra[:mz_array][scan_idx], 
            spectra[:intensity_array][scan_idx], 
            mem,
            spectra[:highMz][scan_idx],
            UInt32(scan_idx), 
            ms_file_idx
        )
        
        sort!(@view(getIonMatches(search_data)[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        # Process matches
        if nmatches > 2
            buildDesignMatrix!(Hs, getIonMatches(search_data), getIonMisses(search_data), 
                             nmatches, nmisses, getIdToCol(search_data))

            if getIdToCol(search_data).size > length(getSpectralScores(search_data))
                new_entries = getIdToCol(search_data).size - length(getSpectralScores(search_data)) + 1000 
                append!(getSpectralScores(search_data), Vector{eltype(getSpectralScores(search_data))}(undef, new_entries))
                append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
            end 

            getDistanceMetrics(Hs, getSpectralScores(search_data))

            ScoreFragmentMatches!(
                getUnscoredPsms(search_data),
                getIdToCol(search_data),
                getIonMatches(search_data), 
                nmatches, 
                mem,
                last(getMinTopNofM(params))
            )

            last_val = Score!(
                getScoredPsms(search_data), 
                getUnscoredPsms(search_data),
                getSpectralScores(search_data),
                getIdToCol(search_data),
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(spectra[:intensity_array][scan_idx])), 
                scan_idx,
                min_spectral_contrast = getMinSpectralContrast(params),
                min_log2_matched_ratio = getMinLog2MatchedRatio(params),
                min_frag_count = getMinFragCount(params),
                max_best_rank = getMaxBestRank(params),
                min_topn = first(getMinTopNofM(params)),
                block_size = 500000
            )
        end
        # Reset arrays
        for scan_idx in range(1, Hs.n)
            getUnscoredPsms(search_data)[scan_idx] = eltype(getUnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end
    return DataFrame(@view(getScoredPsms(search_data)[1:last_val]))
end

function library_search(spectra::Arrow.Table, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:ParameterTuningSearchParameters}
    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getPresearchFragmentIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getIRTTol(search_parameters),
                )...)
end

function library_search(spectra::Arrow.Table, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:SearchParameters}
    
    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getFragmentIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getIrtErrs(search_context)[getParsedFileName(search_context, ms_file_idx)]
                )...)
end

function library_search(
    spectra::Arrow.Table,
    search_context::SearchContext,
    search_parameters::P,
    ms_file_idx::Int64) where {P<:NceTuningSearchParameters}

    return LibrarySearchNceTuning(
        spectra,
        UInt32(ms_file_idx),
        getFragmentIndex(getSpecLib(search_context)),
        getSpecLib(search_context),
        getSearchData(search_context),
        getQuadTransmissionModel(search_context, ms_file_idx),
        getMassErrorModel(search_context, ms_file_idx),
        getRtIrtModel(search_context, ms_file_idx),
        search_parameters,
        search_parameters.nce_grid,
        getIrtErrs(search_context)[getParsedFileName(search_context, ms_file_idx)]
    )
end

function LibrarySearch(
    spectra::Arrow.Table,
    ms_file_idx::UInt32,
    fragment_index::FragmentIndex{Float32},
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel, 
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures, 
        P<:FragmentIndexSearchParameters}
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
                return searchFragmentIndex(
                        scan_to_prec_idx,
                        fragment_index,
                        spectra,
                        last(thread_task),
                        search_data[thread_id],
                        params,
                        qtm,
                        mem,
                        rt_to_irt_spline,
                        irt_tol,
                    )
        end
    end
    
    precursors_passed_scoring = fetch.(tasks)

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                ms_file_idx,
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                getPrecursors(spec_lib),
                                getFragmentLookupTable(spec_lib),
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                search_data[thread_id],
                                params,
                                qtm,
                                mem,
                                rt_to_irt_spline,
                                irt_tol)
        end
    end

    return fetch.(tasks)
end

function LibrarySearchNceTuning(
    spectra::Arrow.Table,
    ms_file_idx::UInt32,
    fragment_index::FragmentIndex{Float32},
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    nce_grid::AbstractVector{<:AbstractFloat},
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel, 
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures, 
        P<:FragmentIndexSearchParameters}

    thread_tasks = partition_scans(spectra, Threads.nthreads())
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))

    # Do fragment index search once
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return searchFragmentIndex(
                scan_to_prec_idx,
                fragment_index,
                spectra,
                last(thread_task),
                search_data[thread_id],
                params,
                qtm,
                mem,
                rt_to_irt_spline,
                irt_tol
            )
        end
    end
    
    precursors_passed_scoring = fetch.(tasks)

    # For each NCE value, run getPSMS using the same fragment index results
    all_results = map(nce_grid) do nce
        # Update NCE model in fragment lookup table

        setNceModel!(
            getFragmentLookupTable(spec_lib),
            PiecewiseNceModel(nce)
        )

        # Run getPSMS with updated NCE model
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin 
                thread_id = first(thread_task)
                psms = getPSMS(
                    ms_file_idx,
                    spectra,
                    last(thread_task),
                    getPrecursors(spec_lib),
                    getFragmentLookupTable(spec_lib),
                    scan_to_prec_idx,
                    precursors_passed_scoring[thread_id],
                    search_data[thread_id],
                    params,
                    qtm,
                    mem,
                    rt_to_irt_spline,
                    irt_tol
                )
                if !isempty(psms)
                    psms[!, :nce] .= nce
                end
                psms
            end
        end
        
        vcat(fetch.(tasks)...)
    end

    return vcat(all_results...)
end

function massErrorSearch(
    spectra::Arrow.Table,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P) where {
        M<:MassErrorModel, 
        S<:SearchDataStructures, 
        P<:SearchParameters
        }

    function getScanToPrecIdx(scan_idxs::Vector{UInt32}, n_scans::Int64)
        scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_scans)
        start_idx = stop_idx = 1
        
        for i in 1:length(scan_idxs)
            stop_idx = i
            if scan_idxs[start_idx] == scan_idxs[stop_idx]
                scan_to_prec_idx[scan_idxs[i]] = start_idx:stop_idx
            else
                scan_to_prec_idx[scan_idxs[i]] = i:i
                start_idx = i
            end
        end
        scan_to_prec_idx
    end

    function collectFragErrs(fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int) where {M<:MatchIon{Float32}}
        for match in range(1, nmatches)
            if n < length(fmatches)
                n += 1
                fmatches[n] = new_fmatches[match]
            else
                fmatches = append!(fmatches, [M() for x in range(1, length(fmatches))])
            end
        end
        return n
    end

    # Sort scans and setup thread tasks
    sorted_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_indices]
    precursor_idxs = precursor_idxs[sorted_indices]
    
    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra[:mz_array]))
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    # Process mass errors in parallel
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            frag_err_idx = 0
            
            for scan_idx in last(thread_task)
                (scan_idx == 0 || scan_idx > length(spectra[:mz_array])) && continue
                ismissing(scan_to_prec_idx[scan_idx]) && continue

                # Select transitions for mass error estimation
                ion_idx = selectTransitions!(
                    getIonTemplates(search_data[thread_id]),
                    MassErrEstimationStrategy(),
                    getPrecEstimation(params),
                    getFragmentLookupTable(spec_lib),
                    scan_to_prec_idx[scan_idx],
                    precursor_idxs,
                    max_rank = 5#getMaxBestRank(params)
                )

                # Match peaks and collect errors
                nmatches, nmisses = matchPeaks!(
                    getIonMatches(search_data[thread_id]),
                    getIonMisses(search_data[thread_id]),
                    getIonTemplates(search_data[thread_id]),
                    ion_idx,
                    spectra[:mz_array][scan_idx],
                    spectra[:intensity_array][scan_idx],
                    mem,
                    spectra[:highMz][scan_idx],
                    UInt32(scan_idx),
                    ms_file_idx
                )

                frag_err_idx = collectFragErrs(
                    getMassErrMatches(search_data[thread_id]),
                    getIonMatches(search_data[thread_id]),
                    nmatches,
                    frag_err_idx
                )
            end
            
            @view(getMassErrMatches(search_data[thread_id])[1:frag_err_idx])
        end
    end
    fetch.(tasks)
end
function getMassErrors(
                    spectra::Arrow.Table,
                    library_fragment_lookup::LibraryFragmentLookup,
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
        ion_idx = selectTransitions!(ionTemplates, MassErrEstimationStrategy(), FullPrecCapture(),
                                        library_fragment_lookup, scan_to_prec_idx[i],
                                        precursors_passed_scoring, max_rank = max_rank
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
                    library_fragment_lookup::Union{LibraryFragmentLookup, Missing},
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
                    max_frag_rank::UInt8,
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
            ion_idx = selectTransitions!(
                ionTemplates, RTIndexedTransitionSelection(), PartialPrecCapture(), library_fragment_lookup,
                precs_temp, precursors[:mz], precursors[:prec_charge],
                precursors[:sulfur_count], iso_splines,
                getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]), 
                precursor_transmission,
                isotopes, n_frag_isotopes, max_frag_rank, rt_index,
                irt_start, irt_stop, 
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                precursors_passing = nothing,
                isotope_err_bounds = isotope_err_bounds,
                block_size = 10000
            )

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
                    library_fragment_lookup::Union{LibraryFragmentLookup, Missing},
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

        ion_idx = selectTransitions!(
            ionTemplates, QuadEstimationTransitionSelection(), PartialPrecCapture(), library_fragment_lookup,
            scan_idx_to_prec_idx[scan_idx], precursors[:mz],precursors[:prec_charge],
            precursors[:sulfur_count], iso_splines,
            precursor_transmission, isotopes,
            (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
            block_size = 10000
        )
        #=
        ion_idx = selectTransitions!(
            ionTemplates, RTIndexedTransitionSelection(), PartialPrecCapture(), library_fragment_lookup,
            precs_temp, precursors[:mz], precursors[:prec_charge],
            precursors[:sulfur_count], iso_splines, getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),
            precursor_transmission,
            isotopes, n_frag_isotopes, max_frag_rank, rt_index,
            irt_start, irt_stop, 
            (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
            precursors_passing = precursors_passing,
            isotope_err_bounds = isotope_err_bounds,
            block_size = 10000
        )
        =#

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
                    push!(tuning_results[:n_matches], n_matches)                
            end
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
                    library_fragment_lookup::Union{LibraryFragmentLookup, Missing},
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
                    max_frag_rank::UInt8,
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
    #prec_temp_size = 0
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
            #=
            precs_temp_size = 0
            ion_idx, _, prec_temp_size = selectTransitions!(
                ionTemplates, RTIndexedTransitionSelection(), library_fragment_lookup,
                precs_temp, precs_temp_size, precursors[:mz], precursors[:prec_charge],
                precursors[:sulfur_count], iso_splines, quad_transmission_func, precursor_transmission,
                isotopes, n_frag_isotopes, max_frag_rank, rt_index,
                irt_start, irt_stop, 
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                precursors_passing = nothing,
                isotope_err_bounds = isotope_err_bounds,
                block_size = 10000
            )
            =#
            ion_idx = selectTransitions!(
                ionTemplates, RTIndexedTransitionSelection(), PartialPrecCapture(), library_fragment_lookup,
                precs_temp, precursors[:mz], precursors[:prec_charge],
                precursors[:sulfur_count], iso_splines,
                getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]), 
                precursor_transmission,
                isotopes, n_frag_isotopes, max_frag_rank, rt_index,
                irt_start, irt_stop, 
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                precursors_passing = nothing,
                isotope_err_bounds = isotope_err_bounds,
                block_size = 10000
            )

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
                    library_fragment_lookup::Union{LibraryFragmentLookup, Missing},
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
                    isotope_trace_type::IsotopeTraceType,
                    quad_transmission_model::QuadTransmissionModel,
                    isotope_err_bounds::Tuple{Int64, Int64},
                    n_frag_isotopes::Int64,
                    max_frag_rank::UInt8,
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
    quad_transmission_function = nothing
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
            quad_transmission_function = getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx])
            ion_idx = selectTransitions!(
                ionTemplates, RTIndexedTransitionSelection(), PartialPrecCapture(), library_fragment_lookup,
                precs_temp, precursors[:mz], precursors[:prec_charge],
                precursors[:sulfur_count], iso_splines, getQuadTransmissionFunction(quad_transmission_model, spectra[:centerMz][scan_idx], spectra[:isolationWidthMz][scan_idx]),
                precursor_transmission,
                isotopes, n_frag_isotopes, max_frag_rank, rt_index,
                irt_start, irt_stop, 
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                precursors_passing = precursors_passing,
                isotope_err_bounds = isotope_err_bounds,
                block_size = 10000
            )
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
                    library_fragment_lookup::LibraryFragmentLookup,
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


function NceScanningSearch(
    spectra::Arrow.Table,
    params::NamedTuple,
    nce_grid::LinRange{Float32, Int64};
    kwargs...)

    thread_tasks, total_peaks = partitionScansToThreads(
                                                            spectra[:mz_array],
                                                            spectra[:retentionTime],
                                                            spectra[:centerMz],
                                                            spectra[:msOrder],
                                                            Threads.nthreads(),
                                                            1
                                                        )

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))

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
    
    all_results = []
    flt = kwargs[:fragment_lookup_table] 
    prec_estimation = Bool(kwargs[:params]["abreviate_precursor_calc"]) ? PartialPrecCapture() : FullPrecCapture()
    for _nce_ in nce_grid
        flt = updateNceModel(
            flt,
                PiecewiseNceModel(_nce_))
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin
                #kwargs[:fragment_lookup_table] = _nce_
                thread_id = first(thread_task)
                return getPSMS(
                    spectra,
                    last(thread_task),
                    kwargs[:precursors],
                    scan_to_prec_idx,
                    precursors_passed_scoring[thread_id],
                    flt,
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
                    UInt8(kwargs[:params]["max_frag_rank"]),
                    prec_estimation,
                    Float32(kwargs[:irt_tol]),
                    Set(2)
                )
            end
        end
        
        # Fetch results and add NCE column
        results_for_nce = fetch.(tasks)
        tasks_out = vcat(results_for_nce...)
        tasks_out[!, :nce] .= _nce_
        
        # Add to results array
        push!(all_results, tasks_out)
    end
    final_results = vcat(all_results...)
    return final_results

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


