# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

function searchFragmentIndex(
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    frag_index::FragmentIndex{Float32},
                    spectra::MassSpecData,
                    thread_task::Vector{Int64},
                    search_data::S,
                    params::P,
                    qtm::Q,
                    mem::M,
                    rt_to_irt_spline::Any,
                    irt_tol::AbstractFloat
                    ) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1
    for scan_idx in thread_task
        #if scan_idx % 50 != 0
        #    continue
        #end
        # Skip invalid indices
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        getMsOrder(spectra, scan_idx) ∉ getSpecOrder(params) && continue

        # Update RT bin index based on iRT window
        irt_lo, irt_hi = getRTWindow(rt_to_irt_spline(getRetentionTime(spectra, scan_idx)), irt_tol)
        while rt_bin_idx < length(getRTBins(frag_index)) && getHigh(getRTBin(frag_index, rt_bin_idx)) < irt_lo
            rt_bin_idx += 1
        end
        while rt_bin_idx > 1 && getLow(getRTBin(frag_index, rt_bin_idx)) > irt_lo
            rt_bin_idx -= 1
        end
        
        # Fragment index search for matching precursors
        searchScan!(
            getPrecursorScores(search_data),
            getRTBins(frag_index),
            getFragBins(frag_index),
            getFragments(frag_index),
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            rt_bin_idx,
            irt_hi,
            mem,
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            getIsotopeErrBounds(params)
        )

        # Filter precursor matches based on score
        match_count, prec_count = filterPrecursorMatches!(getPrecursorScores(search_data), getMinIndexSearchScore(params))

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
    spectra::MassSpecData,
    thread_task::Vector{Int64},
    precursors::LibraryPrecursors,
    ion_list::LibraryFragmentLookup,
    nce_model::NceModel{Float32},
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
    total_skipped_topn = 0
    Hs = SparseArray(UInt32(5000))
    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    for scan_idx in thread_task
        (scan_idx == 0 || scan_idx > length(spectra)) && continue
        ismissing(scan_to_prec_idx[scan_idx]) && continue

        # Scan Filtering
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1

        # Ion Template Selection
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            StandardTransitionSelection(),
            getPrecEstimation(params),
            ion_list,
            nce_model,
            scan_to_prec_idx[scan_idx], precursors_passed_scoring,
            getMz(precursors),
            getCharge(precursors),
            getSulfurCount(precursors),
            getIrt(precursors),
            getIsoSplines(search_data),
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            precursor_transmission, isotopes, getNFragIsotopes(params),
            getMaxFragRank(params),
            Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx))),
            Float32(irt_tol),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        ion_idx < 2 && continue


        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx)
        )
        
        sort!(@view(getIonMatches(search_data)[1:nmatches]), alg=QuickSort, lt=ion_match_lt)
        # Process matches
        if nmatches > 2
            buildDesignMatrix!(Hs, getIonMatches(search_data), getIonMisses(search_data), 
                             nmatches, nmisses, getIdToCol(search_data))

            if getIdToCol(search_data).size > length(getSpectralScores(search_data))
                new_entries = getIdToCol(search_data).size - length(getSpectralScores(search_data)) + 1000 
                append!(getSpectralScores(search_data), Vector{eltype(getSpectralScores(search_data))}(undef, new_entries))
                append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
            end 

            getDistanceMetrics(Hs, getSpectralScores(search_data); params.relative_improvement_threshold)

            ScoreFragmentMatches!(
                getUnscoredPsms(search_data),
                getIdToCol(search_data),
                getIonMatches(search_data), 
                nmatches, 
                mem,
                last(getMinTopNofM(params))
            )

            score_result = Score!(
                getScoredPsms(search_data),
                getUnscoredPsms(search_data),
                getSpectralScores(search_data),
                getIdToCol(search_data),
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(getIntensityArray(spectra, scan_idx))),
                scan_idx,
                min_spectral_contrast = getMinSpectralContrast(params),
                min_log2_matched_ratio = getMinLog2MatchedRatio(params),
                min_frag_count = getMinFragCount(params),
                max_best_rank = getMaxBestRank(params),
                min_topn = first(getMinTopNofM(params)),
                block_size = 500000
            )
            last_val = score_result.last_val
            total_skipped_topn += score_result.skipped_topn
        end
        # Reset arrays
        for scan_idx in range(1, Hs.n)
            getUnscoredPsms(search_data)[scan_idx] = eltype(getUnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end
    @debug "FirstPass filter summary: kept=$last_val, skipped_topn=$total_skipped_topn"
    return DataFrame(@view(getScoredPsms(search_data)[1:last_val]))
end

function library_search(spectra::MassSpecData, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:ParameterTuningSearchParameters}
    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getPresearchPartitionedIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getNceModel(search_context, ms_file_idx),
                    getIRTTol(search_parameters),
                )...)
end

function library_search(spectra::MassSpecData, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:SearchParameters}

    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getPartitionedIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getNceModel(search_context, ms_file_idx),
                    getIrtErrors(search_context)[ms_file_idx]
                )...)
end

function library_search(
    spectra::MassSpecData,
    search_context::SearchContext,
    search_parameters::P,
    ms_file_idx::Int64) where {P<:NceTuningSearchParameters}

    return LibrarySearchNceTuning(
        spectra,
        UInt32(ms_file_idx),
        getPartitionedIndex(getSpecLib(search_context)),
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
    spectra::MassSpecData,
    ms_file_idx::UInt32,
    partitioned_index::LocalPartitionedFragmentIndex{Float32},
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    nce_model::NceModel{Float32},
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel,
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures,
        P<:FragmentIndexSearchParameters}
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    n_threads = Threads.nthreads()
    precursor_mzs = getMz(getPrecursors(spec_lib))

    # Collect valid MS2 scan indices
    all_scan_idxs = Int[]
    for tt in thread_tasks
        append!(all_scan_idxs, last(tt))
    end
    filter!(si -> si > 0 && si <= length(spectra) && getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)

    # Single call — handles threading internally
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    precursors_passed_scoring = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    # getPSMS phase: pass the single flat vector
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            return getPSMS(
                ms_file_idx,
                spectra,
                last(thread_task),
                getPrecursors(spec_lib),
                getFragmentLookupTable(spec_lib),
                nce_model,
                scan_to_prec_idx,
                precursors_passed_scoring,
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
    spectra::MassSpecData,
    ms_file_idx::UInt32,
    partitioned_index::LocalPartitionedFragmentIndex{Float32},
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

    # Check for valid inputs
    if isempty(nce_grid)
        @user_warn "LibrarySearchNceTuning: Empty NCE grid provided"
        return DataFrame()
    end

    @debug_l2 "LibrarySearchNceTuning: Starting with $(length(spectra)) scans, NCE grid: $nce_grid"

    thread_tasks = partition_scans(spectra, Threads.nthreads())
    n_threads = Threads.nthreads()
    precursor_mzs = getMz(getPrecursors(spec_lib))

    @debug_l2 "LibrarySearchNceTuning: Created $(length(thread_tasks)) thread tasks"
    for (i, task) in enumerate(thread_tasks)
        task_range = last(task)
        n_scans = length(task_range)
        @debug_l2 "  Thread $i: $n_scans scans"
    end

    # Collect valid MS2 scan indices
    all_scan_idxs = Int[]
    for tt in thread_tasks
        append!(all_scan_idxs, last(tt))
    end
    filter!(si -> si > 0 && si <= length(spectra) && getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)

    # Single fragment index search
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    precursors_passed_scoring = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    # For each NCE value, run getPSMS using the same fragment index results
    all_results = map(nce_grid) do nce
        # Create NCE model for this grid point
        nce_model = PiecewiseNceModel(nce)

        # Run getPSMS with this NCE model
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin
                thread_id = first(thread_task)
                try
                    psms = getPSMS(
                        ms_file_idx,
                        spectra,
                        last(thread_task),
                        getPrecursors(spec_lib),
                        getFragmentLookupTable(spec_lib),
                        nce_model,
                        scan_to_prec_idx,
                        precursors_passed_scoring,
                        search_data[thread_id],
                        params,
                        qtm,
                        mem,
                        rt_to_irt_spline,
                        irt_tol
                    )
                    if !isempty(psms)
                        psms[!, :nce] .= nce
                        return psms
                    else
                        return missing
                    end
                catch e
                    @user_warn "PSM search failed on thread $thread_id for NCE $nce: $e"
                    return missing
                end
            end
        end
        # Ensure we always return a DataFrame, not Vector{Any}
        fetched = collect(skipmissing(fetch.(tasks)))
        isempty(fetched) ? DataFrame() : vcat(fetched...)
    end

    # filter out empty DFs (which are actually Vectors instead of DataFrames)
    nonempty_dfs = filter(df -> df isa DataFrame, all_results)
    # Ensure we always return a DataFrame, not Vector{Any}
    return isempty(nonempty_dfs) ? DataFrame() : vcat(nonempty_dfs...)
end

"""
    fragment_index_search_only(spectra, search_context, params, ms_file_idx)

Run only Phase 1 of LibrarySearch (fragment index search) without getPSMS scoring.
Returns (scan_to_prec_idx, precursors_passed) for direct use in second pass.
"""
function fragment_index_search_only(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:FragmentIndexSearchParameters}
    spec_lib = getSpecLib(search_context)
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt_spline = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = haskey(getIrtErrors(search_context), ms_file_idx) ? getIrtErrors(search_context)[ms_file_idx] : Float32(Inf)
    precursor_mzs = getMz(getPrecursors(spec_lib))
    n_threads = Threads.nthreads()
    partitioned_index = getPartitionedIndex(spec_lib)

    thread_tasks = partition_scans(spectra, n_threads)

    # Collect valid MS2 scan indices
    all_scan_idxs = Int[]
    for tt in thread_tasks
        append!(all_scan_idxs, last(tt))
    end
    filter!(si -> si > 0 && si <= length(spectra) && getMsOrder(spectra, si) ∈ getSpecOrder(params), all_scan_idxs)

    # Single call — handles threading internally
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    precursors_passed = searchFragmentIndexPartitionMajorHinted(
        scan_to_prec_idx, partitioned_index, spectra, all_scan_idxs,
        n_threads, params, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    return scan_to_prec_idx, precursors_passed
end

"""
    write_fragment_index_matches(scan_to_prec_idx, precursors_passed, output_path)

Flatten (scan_idx, precursor_idx) pairs from fragment index search and write to Arrow.
"""
function write_fragment_index_matches(
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    output_path::String
)
    scan_idxs = UInt32[]
    prec_idxs = UInt32[]
    for (scan_idx, range) in enumerate(scan_to_prec_idx)
        ismissing(range) && continue
        for i in range
            push!(scan_idxs, UInt32(scan_idx))
            push!(prec_idxs, precursors_passed[i])
        end
    end
    Arrow.write(output_path, (scan_idx=scan_idxs, precursor_idx=prec_idxs))
end

"""
    load_fragment_index_matches(path, n_scans)

Read Arrow file and reconstruct scan_to_prec_idx + precursors_passed vectors.
"""
function load_fragment_index_matches(path::String, n_scans::Int)
    tbl = Arrow.Table(path)
    scan_idxs = tbl[:scan_idx]
    prec_idxs = tbl[:precursor_idx]

    precursors_passed = Vector{UInt32}(prec_idxs)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(missing, n_scans)

    n = length(scan_idxs)
    n == 0 && return scan_to_prec_idx, precursors_passed

    i = 1
    while i <= n
        current_scan = scan_idxs[i]
        start_idx = i
        while i <= n && scan_idxs[i] == current_scan
            i += 1
        end
        scan_to_prec_idx[current_scan] = start_idx:(i - 1)
    end

    return scan_to_prec_idx, precursors_passed
end

function getRTWindow(irt::U, irt_tol::T) where {T,U<:AbstractFloat}
    return Float32(irt - irt_tol), Float32(irt + irt_tol)
end
