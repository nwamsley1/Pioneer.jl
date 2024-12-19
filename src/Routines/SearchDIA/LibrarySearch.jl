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
        #if scan_idx % 50 != 0
        #    continue
        #end
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
    precursors::BasicLibraryPrecursors,
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
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data), 
            StandardTransitionSelection(), 
            getPrecEstimation(params),
            ion_list,
            scan_to_prec_idx[scan_idx], precursors_passed_scoring,
            getMz(precursors),#[:mz], 
            getCharge(precursors),#[:prec_charge], 
            getSulfurCount(precursors),#[:sulfur_count],
            getIrt(precursors),#[:irt],
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
                    getIrtErrors(search_context)[ms_file_idx]
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
        getIrtErrors(search_context)[ms_file_idx]
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
                else
                    psms = missing
                end

                return psms
            end
        end
        
        vcat(skipmissing(fetch.(tasks))...)
    end

    return vcat(all_results...)
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

