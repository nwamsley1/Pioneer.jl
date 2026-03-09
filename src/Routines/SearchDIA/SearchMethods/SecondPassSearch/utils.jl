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

#==========================================================
Core Search Funtions
==========================================================#
"""
    perform_second_pass_search(spectra::MassSpecData, rt_index::retentionTimeIndex,
                             search_context::SearchContext, params::SecondPassSearchParameters,
                             ms_file_idx::Int64) -> DataFrame

Execute second pass search across MS/MS data.

# Arguments
- `spectra`: MS/MS spectral data
- `rt_index`: Retention time index for efficient searching
- `search_context`: Search context with libraries and models
- `params`: Search parameters
- `ms_file_idx`: MS file index

# Process
1. Partitions scans across threads
2. Processes each scan batch in parallel
3. Combines results into single DataFrame
"""
function perform_second_pass_search(
    spectra::MassSpecData,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            return process_scans!(
                last(thread_task),
                spectra,
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx,
                MS2CHROM()
            )
        end
    end
    
    # Collect results with detailed error logging per task
    results = Vector{DataFrame}(undef, length(tasks))
    for (i, t) in enumerate(tasks)
        try
            results[i] = fetch(t)
        catch e
            bt = catch_backtrace()
            @user_error "SecondPassSearch task $(i) failed while fetching results (MS2CHROM)"
            @user_error sprint(showerror, e, bt)
            rethrow(e)
        end
    end
    return vcat(results...)
end

function perform_second_pass_search(
    spectra::MassSpecData,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    precursors_passing::Set{UInt32},
    isotopes_dict::UnorderedDictionary{UInt32, Vector{Isotope{T}}},
    ::MS1CHROM
) where {T<:AbstractFloat}
    thread_tasks = partition_scans(spectra, Threads.nthreads(), ms_order_select = 1)
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            return process_scans!(
                last(thread_task),
                spectra,
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx,
                precursors_passing,
                isotopes_dict,
                MS1CHROM()
            )
        end
    end
    # Collect results with detailed error logging per task
    results = Vector{DataFrame}(undef, length(tasks))
    for (i, t) in enumerate(tasks)
        try
            results[i] = fetch(t)
        catch e
            bt = catch_backtrace()
            @user_error "SecondPassSearch task $(i) failed while fetching results (MS1CHROM)"
            @user_error sprint(showerror, e, bt)
            rethrow(e)
        end
    end
    return vcat(results...)
end

"""
    perform_second_pass_search(spectra, scan_to_prec_idx, precursors_passed,
                             search_context, params, ms_file_idx, ::MS2CHROM)

Execute second pass search using pre-computed fragment index matches
instead of RT index. Each scan has its own precursor list.
"""
function perform_second_pass_search(
    spectra::MassSpecData,
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            return process_scans_fragindex!(
                last(thread_task),
                spectra,
                scan_to_prec_idx,
                precursors_passed,
                search_context,
                search_data,
                params,
                ms_file_idx
            )
        end
    end

    results = Vector{DataFrame}(undef, length(tasks))
    for (i, t) in enumerate(tasks)
        try
            results[i] = fetch(t)
        catch e
            bt = catch_backtrace()
            @user_error "SecondPassSearch (fragindex) task $(i) failed"
            @user_error sprint(showerror, e, bt)
            rethrow(e)
        end
    end
    return vcat(results...)
end

"""
    process_scans_fragindex!(scan_range, spectra, scan_to_prec_idx, precursors_passed,
                            search_context, search_data, params, ms_file_idx)

Process scans using fragment index matches. Each scan gets its own precursor list
from scan_to_prec_idx, eliminating the RT-based caching optimization.
"""
function process_scans_fragindex!(
    scan_range::Vector{Int64},
    spectra::MassSpecData,
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed::Vector{UInt32},
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    last_val = 0
    cycle_idx = 0
    total_skipped_weight = 0
    total_skipped_frag_count = 0

    nce_model = getNceModel(search_context, ms_file_idx)
    precursors = getPrecursors(getSpecLib(search_context))

    for scan_idx in scan_range
        ((scan_idx < 1) || scan_idx > length(spectra)) && continue
        msn = getMsOrder(spectra, scan_idx)
        if msn < 2
            cycle_idx += 1
        end
        msn ∉ params.spec_order && continue

        # Skip scans with no fragment index matches
        ismissing(scan_to_prec_idx[scan_idx]) && continue

        # Select transitions using StandardTransitionSelection with explicit precursor list
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            StandardTransitionSelection(),
            params.prec_estimation,
            getFragmentLookupTable(getSpecLib(search_context)),
            nce_model,
            scan_to_prec_idx[scan_idx],
            precursors_passed,
            getMz(precursors),
            getCharge(precursors),
            getSulfurCount(precursors),
            getIrt(precursors),
            getIsoSplines(search_data),
            getQuadTransmissionFunction(
                getQuadTransmissionModel(search_context, ms_file_idx),
                getCenterMz(spectra, scan_idx),
                getIsolationWidthMz(spectra, scan_idx)
            ),
            getPrecursorTransmission(search_data),
            getIsotopes(search_data),
            params.n_frag_isotopes,
            params.max_frag_rank,
            Float32(0.0),       # iRT value (not used for filtering)
            Float32(1e10),      # iRT_tol = effectively infinite (no RT filtering)
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = params.isotope_err_bounds,
            block_size = 10000
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
            getMassErrorModel(search_context, ms_file_idx),
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )

        nmatches ≤ 2 && continue

        # Build design matrix
        buildDesignMatrix!(
            Hs,
            getIonMatches(search_data),
            getIonMisses(search_data),
            nmatches,
            nmisses,
            getIdToCol(search_data)
        )

        # Handle weights
        if getIdToCol(search_data).size > length(weights)
            resize_arrays!(search_data, weights)
        end

        initialize_weights!(search_data, weights, precursor_weights)

        # Solve deconvolution
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            getHuberDelta(search_context),
            params.lambda,
            params.max_iter_newton,
            params.max_iter_bisection,
            params.max_iter_outer,
            search_context.deconvolution_stop_tolerance[],
            search_context.deconvolution_stop_tolerance[],
            params.max_diff,
            params.reg_type
        )

        # Update precursor weights
        update_precursor_weights!(search_data, weights, precursor_weights)

        # Score PSMs
        getDistanceMetrics(weights, residuals, Hs, getComplexSpectralScores(search_data))

        ScoreFragmentMatches!(
            getComplexUnscoredPsms(search_data),
            getIdToCol(search_data),
            getIonMatches(search_data),
            nmatches,
            getMassErrorModel(search_context, ms_file_idx),
            last(params.min_topn_of_m)
        )

        score_result = Score!(
            getComplexScoredPsms(search_data),
            getComplexUnscoredPsms(search_data),
            getComplexSpectralScores(search_data),
            weights,
            getIdToCol(search_data),
            cycle_idx,
            nmatches/(nmatches + nmisses),
            last_val,
            Hs.n,
            Float32(sum(getIntensityArray(spectra, scan_idx))),
            scan_idx;
            min_spectral_contrast = params.min_spectral_contrast,
            min_log2_matched_ratio = params.min_log2_matched_ratio,
            min_y_count = params.min_y_count,
            min_frag_count = params.min_frag_count,
            max_best_rank = params.max_best_rank,
            min_topn = first(params.min_topn_of_m),
            block_size = 500000
        )
        last_val = score_result.last_val
        total_skipped_weight += score_result.skipped_weight
        total_skipped_frag_count += score_result.skipped_frag_count

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    @info "SecondPass (fragindex) filter summary: kept=$last_val, skipped_weight=$total_skipped_weight, skipped_frag_count=$total_skipped_frag_count"
    return DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val]))
end

"""
    process_scans!(scan_range::Vector{Int64}, spectra::MassSpecData,
                  rt_index::retentionTimeIndex, search_context::SearchContext,
                  search_data::SearchDataStructures, params::SecondPassSearchParameters,
                  ms_file_idx::Int64) -> DataFrame

Process a batch of scans with RT bin caching.

# Process
1. Tracks RT bins and m/z windows for efficient transition selection
2. For each scan:
   - Selects transitions based on RT window
   - Matches peaks
   - Builds design matrix
   - Performs deconvolution
   - Scores PSMs
"""
function process_scans!(
    scan_range::Vector{Int64},
    spectra::MassSpecData,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    last_val = 0
    total_skipped_weight = 0
    total_skipped_frag_count = 0

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    ion_idx = 0
    cycle_idx = 0

    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    nce_model = getNceModel(search_context, ms_file_idx)

    for scan_idx in scan_range
        ((scan_idx < 1) || scan_idx > length(spectra)) && continue
        msn = getMsOrder(spectra, scan_idx)
        if msn < 2
            cycle_idx += 1
        end
        msn ∉ params.spec_order && continue

        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Check for m/z change
        prec_mz_string_new = string(getCenterMz(spectra, scan_idx))
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]

        # Update transitions if window changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || 
           (prec_mz_string_new != prec_mz_string)
            
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new

            ion_idx, _ = selectTransitions!(
                getIonTemplates(search_data),
                RTIndexedTransitionSelection(),
                params.prec_estimation,
                getFragmentLookupTable(getSpecLib(search_context)),
                nce_model,
                getPrecIds(search_data),
                getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
                getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
                getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
                getIsoSplines(search_data),
                getQuadTransmissionFunction(
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getCenterMz(spectra, scan_idx),
                    getIsolationWidthMz(spectra, scan_idx)
                ),
                getPrecursorTransmission(search_data),
                getIsotopes(search_data),
                params.n_frag_isotopes,
                params.max_frag_rank,
                rt_index,
                irt_start,
                irt_stop,
                (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
                block_size = 10000,
                min_fraction_transmitted = params.min_fraction_transmitted
            )
        end

        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            getMassErrorModel(search_context, ms_file_idx),
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )

        nmatches ≤ 2 && continue

        # Process scan
        buildDesignMatrix!(
            Hs,
            getIonMatches(search_data),
            getIonMisses(search_data),
            nmatches,
            nmisses,
            getIdToCol(search_data)
        )

        # Handle weights
        if getIdToCol(search_data).size > length(weights)
            resize_arrays!(search_data, weights)
        end

        initialize_weights!(search_data, weights, precursor_weights)
        
        # Solve deconvolution problem
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            getHuberDelta(search_context),
            params.lambda,
            params.max_iter_newton,
            params.max_iter_bisection,
            params.max_iter_outer,
            search_context.deconvolution_stop_tolerance[],#params.accuracy_newton,
            search_context.deconvolution_stop_tolerance[],#params.accuracy_bisection,
            params.max_diff,
            params.reg_type
        )

        # Update precursor weights
        update_precursor_weights!(search_data, weights, precursor_weights)

        # Score PSMs
        getDistanceMetrics(weights, residuals, Hs, getComplexSpectralScores(search_data))
        
        ScoreFragmentMatches!(
            getComplexUnscoredPsms(search_data),
            getIdToCol(search_data),
            getIonMatches(search_data),
            nmatches,
            getMassErrorModel(search_context, ms_file_idx),
            last(params.min_topn_of_m)
        )

        score_result = Score!(
            getComplexScoredPsms(search_data),
            getComplexUnscoredPsms(search_data),
            getComplexSpectralScores(search_data),
            weights,
            getIdToCol(search_data),
            cycle_idx,
            nmatches/(nmatches + nmisses),
            last_val,
            Hs.n,
            Float32(sum(getIntensityArray(spectra, scan_idx))),
            scan_idx;
            min_spectral_contrast = params.min_spectral_contrast,
            min_log2_matched_ratio = params.min_log2_matched_ratio,
            min_y_count = params.min_y_count,
            min_frag_count = params.min_frag_count,
            max_best_rank = params.max_best_rank,
            min_topn = first(params.min_topn_of_m),
            block_size = 500000
        )
        last_val = score_result.last_val
        total_skipped_weight += score_result.skipped_weight
        total_skipped_frag_count += score_result.skipped_frag_count

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    @info "SecondPass (RT-indexed) filter summary: kept=$last_val, skipped_weight=$total_skipped_weight, skipped_frag_count=$total_skipped_frag_count"
    return DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val]))
end

"""
    build_chromatograms(spectra::MassSpecData, scan_range::Vector{Int64},
                       precursors_passing::Set{UInt32}, rt_index::retentionTimeIndex,
                       search_context::SearchContext, search_data::SearchDataStructures,
                       params::IntegrateChromatogramSearchParameters,
                       ms_file_idx::Int64) -> DataFrame

Build chromatograms for a range of scans with RT bin caching.

# Process
1. Tracks RT bins for efficient transition selection
2. Selects transitions based on RT windows
3. Matches peaks and performs deconvolution
4. Records chromatogram points with weights
"""
function process_scans!(
    scan_range::Vector{Int64},
    spectra::MassSpecData,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    precursors_passing::Set{UInt32},
    isotopes_dict::UnorderedDictionary{UInt32, Vector{Isotope{T}}},
    ::MS1CHROM
) where {T<:AbstractFloat}
    #######
    # Initialize working arrays
    #mem = MassErrorModel(
    #    getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
    #    (6.0f0, 6.0f0)
    #)
    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    ion_templates = Vector{Isotope{Float32}}(undef, 100000)
    ion_matches = [PrecursorMatch{Float32}() for _ in range(1, 10000)]
    ion_misses = [PrecursorMatch{Float32}() for _ in range(1, 10000)]
    precursors = getPrecursors(getSpecLib(search_context))
    pair_ids = getPairIdx(precursors)
    pair_id_dict = Dictionary{
        UInt32, #pair_id
        UInt8 #Number of matches to spectrum
    }()

    # NEW: Create m/z grouping map for MS1
    mz_grouping = MzGroupingMap(UInt32(100000))  # 5 decimal place precision

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    ion_idx = 0
    last_val = 0
    # Note: prec_ids scratch array not required in MS1 path
    irt_tol = getIrtErrors(search_context)[ms_file_idx]

    for scan_idx in scan_range
        empty!(pair_id_dict)
        ((scan_idx < 1) || (scan_idx > length(spectra))) && continue

        # Process MS1 scans
        if getMsOrder(spectra, scan_idx) != 1
            continue
        end

        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
        irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Update transitions if window changed
        # reset per-scan counters
        ion_idx = 0
        for rt_bin_idx in irt_start:irt_stop #All retention time bins covering the current scan 
            precs = rt_index.rt_bins[rt_bin_idx].prec #Precursor idxs for the current retention time bin
            for i in 1:length(precs)
                prec_idx = first(precs[i])
                if prec_idx in precursors_passing #If the precursor passed the first search ('precursors_passing' is a set)
                    pair_id = pair_ids[prec_idx] #Each target precursor has a complement decoy with the same 'pair_id' and vice versa
                    if !haskey(pair_id_dict, pair_id)
                        insert!(pair_id_dict, pair_id, zero(UInt8))
                    else
                        #If another precursor (the respective target or decoy complement) in this pair has already been added. 
                        continue
                    end
                    # Track only in local structures for MS1 path (no scratch needed)


                    for iso in isotopes_dict[prec_idx] #Add the isotopes for the precursor to match to the scan 
                        ion_idx += 1
                        if ion_idx > length(ion_templates)
                            append!(ion_templates, Vector{Isotope{Float32}}(undef, length(ion_templates)))
                        end
                        ion_templates[ion_idx] = iso
                    end
                end
            end
        end
        
        #Sort the precursor isotopes by m/z
        sort!(@view(ion_templates[1:ion_idx]), by = x->(getMZ(x)), alg=PartialQuickSort(1:ion_idx))
        # Match peaks
        # Ensure match/miss buffers are large enough (cannot exceed ion_idx)
        if ion_idx > length(ion_matches)
            append!(ion_matches, [PrecursorMatch{Float32}() for _ in 1:max(ion_idx - length(ion_matches), length(ion_matches))])
        end
        if ion_idx > length(ion_misses)
            append!(ion_misses, [PrecursorMatch{Float32}() for _ in 1:max(ion_idx - length(ion_misses), length(ion_misses))])
        end
        nmatches, nmisses = matchPeaks!(
            ion_matches,
            ion_misses,
            ion_templates,
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )

        #Which precursors matched isotopes
        for i in range(1, nmatches)
            pair_id = pair_ids[getPrecID(ion_matches[i])]
            n_match = pair_id_dict[pair_id]
            if getIsoIdx(ion_matches[i]) <= 2
                n_match += 1
            end
            pair_id_dict[pair_id] = n_match 
        end

        #Removed matched fragments for precursors that did not match sufficiently many isotopes 
        new_nmatches = 0
        for i in range(1, nmatches)
            pair_id = pair_ids[getPrecID(ion_matches[i])]
            n_match = pair_id_dict[pair_id]
            if n_match >= 2
                new_nmatches += 1
                ion_matches[new_nmatches] = ion_matches[i]
            end
        end

        #Removed unmatched fragments for precursors that did not match sufficiently many isotopes 
        new_nmisses = 0
        for i in range(1, nmisses)
            pair_id = pair_ids[getPrecID(ion_misses[i])]
            n_match = pair_id_dict[pair_id]
            if n_match >= 2
                new_nmisses += 1
                ion_misses[new_nmisses] = ion_misses[i]
            end
        end

        nmatches = new_nmatches
        nmisses = new_nmisses 
        
        sort!(@view(ion_matches[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        # Process matches
        if nmatches > 2

            # Reset grouping for this scan
            reset!(mz_grouping)

            # Use MS1-specific design matrix construction with m/z grouping
            buildDesignMatrixMS1!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                mz_grouping,
                precursors
            )

            # m/z grouping completed for this scan

            # Populate id_to_col mapping from m/z grouping so MS1 scoring has valid indices
            id_to_col = getIdToCol(search_data)
            reset!(id_to_col)
            for (mz_group, col) in mz_grouping.mz_to_col
                if haskey(mz_grouping.mz_group_to_precids, mz_group)
                    for prec_id in mz_grouping.mz_group_to_precids[mz_group]
                        update!(id_to_col, prec_id, col)
                    end
                end
            end

            # Ensure arrays are sized for the number of grouped columns
            if mz_grouping.current_col > length(weights)
                new_entries = Int(mz_grouping.current_col) - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(getMs1SpectralScores(search_data), length(getMs1SpectralScores(search_data)) + new_entries)
                append!(getMs1UnscoredPsms(search_data), [eltype(getMs1UnscoredPsms(search_data))() for _ in 1:new_entries])
            end

            # Initialize active group weights to zero (simple baseline)
            @inbounds @fastmath for col in 1:Int(mz_grouping.current_col)
                weights[col] = 0.0f0
            end

            # Solve deconvolution
            initResiduals!(residuals, Hs, weights)
            solveHuber!(
                Hs,
                residuals,
                weights,
                params.ms1_huber_delta,
                params.ms1_lambda,
                params.max_iter_newton,
                params.max_iter_bisection,
                params.max_iter_outer,
                search_context.deconvolution_stop_tolerance[],#params.accuracy_newton,
                search_context.deconvolution_stop_tolerance[],#params.accuracy_bisection,
                params.max_diff,
                params.ms1_reg_type
            )

            # NEW: Distribute grouped coefficients back to individual precursors
            distribute_ms1_coefficients!(
                precursor_weights,  # Array indexed by precursor ID
                weights,            # Array indexed by column number (group coefficients)
                mz_grouping
            )

            # Update precursor weights - already handled by distribute_ms1_coefficients!
            # No need for additional update

            getDistanceMetrics(weights, residuals, Hs, getMs1SpectralScores(search_data))

            ScoreFragmentMatches!(
                getMs1UnscoredPsms(search_data),
                getIdToCol(search_data),
                ion_matches,
                nmatches,
                mem,
                last(params.min_topn_of_m)
                )
        end

        last_val = Score!(
            getMs1ScoredPsms(search_data),
            getMs1UnscoredPsms(search_data),
            getMs1SpectralScores(search_data),
            weights,
            getIdToCol(search_data),
            zero(Int64),
            last_val,
            Hs.n,
            scan_idx;
            block_size = 500000
        )
        # Reset arrays
        for i in 1:Hs.n
            getMs1UnscoredPsms(search_data)[i] = eltype(getMs1UnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end

    return DataFrame(@view(getMs1ScoredPsms(search_data)[1:last_val]))
end


#==========================================================
Temporarry array management functions
==========================================================#
"""
    resize_arrays!(search_data::SearchDataStructures, weights::Vector{Float32})

Resize working arrays when needed for larger deconvolution problems.

Expands:
- weights array
- spectral scores
- unscored PSMs array
"""
function resize_arrays!(search_data::SearchDataStructures, weights::Vector{Float32})
    new_entries = getIdToCol(search_data).size - length(weights) + 1000
    resize!(weights, length(weights) + new_entries)
    resize!(getComplexSpectralScores(search_data), length(getComplexSpectralScores(search_data)) + new_entries)
    append!(getComplexUnscoredPsms(search_data), [eltype(getComplexUnscoredPsms(search_data))() for _ in 1:new_entries])
end

"""
    initialize_weights!(search_data::SearchDataStructures, weights::Vector{Float32},
                       precursor_weights::Vector{Float32})

Initialize weights for deconvolution from precursor weights.
"""
function initialize_weights!(
    search_data::SearchDataStructures,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32}
)
    for i in 1:getIdToCol(search_data).size
        weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
            precursor_weights[getIdToCol(search_data).keys[i]]
    end
end

"""
    update_precursor_weights!(search_data::SearchDataStructures, weights::Vector{Float32},
                            precursor_weights::Vector{Float32})

Update precursor weights after deconvolution solution.
"""
function update_precursor_weights!(
    search_data::SearchDataStructures,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32}
)
    for i in 1:getIdToCol(search_data).size
        id = getIdToCol(search_data).keys[i]
        colid = getIdToCol(search_data)[id]
        precursor_weights[id] = weights[colid]
    end
end


"""
    reset_arrays!(search_data::SearchDataStructures, Hs::SparseArray)

Reset arrays between scan processing.

Clears:
- Unscored PSMs
- ID to column mapping
- Sparse matrix
"""
function reset_arrays!(search_data::SearchDataStructures, Hs::SparseArray)
    for i in 1:Hs.n
        getComplexUnscoredPsms(search_data)[i] = eltype(getComplexUnscoredPsms(search_data))()
    end
    reset!(getIdToCol(search_data))
    reset!(Hs)
end

#==========================================================
PSMs processing
==========================================================#
"""
    add_second_search_columns!(psms::DataFrame, scan_retention_time::AbstractVector{Float32},
                             prec_charge::AbstractVector{UInt8}, prec_is_decoy::AbstractVector{Bool},
                             precursors::LibraryPrecursors)

Add essential columns to PSM DataFrame for second pass analysis.

# Added Columns
- Retention time
- Charge state
- Target/decoy status
- Ion counts
- Error metrics
- CV fold assignments
"""
function add_second_search_columns!(psms::DataFrame, 
                        scan_retention_time::AbstractVector{Float32},
                        prec_charge::AbstractVector{UInt8},
                        prec_is_decoy::AbstractVector{Bool},
                        precursors::LibraryPrecursors,
                        #prec_id_to_cv_fold::Dictionary{UInt32, UInt8})
)
    
    ###########################
    #Correct Weights by base-peak intensity
    filter!(x->x.weight>0.0, psms);
    ###########################
    #Allocate new columns
   
    #Threads.@threads for i in ProgressBar(range(1, size(psms)[1]))
    N = size(psms, 1);
    decoys = zeros(Bool, N);
    rt = zeros(Float32, N);
    #TIC = zeros(Float16, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    prec_charges = zeros(UInt8, N);
    #prec_mzs = zeros(Float32, N);
    cv_fold = zeros(UInt8, N);
    scan_idxs::Vector{UInt32} = psms[!,:scan_idx]
    prec_idxs::Vector{UInt32} = psms[!,:precursor_idx]
    y_count::Vector{UInt8} = psms[!,:y_count]
    b_count::Vector{UInt8} = psms[!,:b_count]
    isotope_count::Vector{UInt8} = psms[!,:isotope_count]
    error::Vector{Float32} = psms[!,:error]
    #psms[!,:total_ions]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    #scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = prec_idxs[i]
                scan_idx = scan_idxs[i]

                decoys[i] = prec_is_decoy[prec_idx];
                targets[i] = decoys[i] == false
                rt[i] = Float32(scan_retention_time[scan_idx]);
                prec_charges[i] = prec_charge[prec_idx];
                total_ions[i] = UInt16(y_count[i] + b_count[i] + isotope_count[i]);
                err_norm[i] = Float16(min((2^error[i])/max(total_ions[i], one(UInt16)), Float32(6e4)))
                if isinf(matched_ratio[i])
                    matched_ratio[i] = Float16(60000)*sign(matched_ratio[i])
                end
                cv_fold[i] = getCvFold(precursors, prec_idx)#prec_id_to_cv_fold[prec_idx]
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:rt] = rt
    #psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = prec_charges
    #psms[!,:prec_mz] = prec_mzs
    psms[!,:cv_fold] = cv_fold
    psms[!,:charge2] = Vector{UInt8}(psms[!, :charge] .== 2)
    #######
    sort!(psms,:rt); #Sorting before grouping is critical. 
    return nothing
end

"""
    get_isotopes_captured!(chroms::DataFrame, isotope_trace_type::IsotopeTraceType,
                          quad_transmission_model::QuadTransmissionModel, ...)

Determine which isotopes are captured in each isolation window.

Uses quadrupole transmission model and scan parameters to calculate
isotope coverage.
"""
function get_isotopes_captured!(chroms::DataFrame, 
                                isotope_trace_type::IsotopeTraceType,
                                quad_transmission_model::QuadTransmissionModel,
                                search_data::Vector{SimpleLibrarySearch{IsotopeSplineModel{40, Float32}}},
                                scan_idx::AbstractVector{UInt32},
                                prec_charge::AbstractArray{UInt8},
                                prec_mz::AbstractArray{Float32},
                                sulfur_count::AbstractArray{UInt8},
                                centerMz::AbstractVector{Union{Missing, Float32}},
                                isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    #sum(MS2_CHROMS.weight.!=0.0)
    isotopes_captured = Vector{Tuple{Int8, Int8}}(undef, size(chroms, 1))
    precursor_fraction_transmitted = Vector{Float32}(undef, size(chroms, 1))
    
    tasks_per_thread = 5
    chunk_size = max(1, size(chroms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(chroms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            thread_id = (first(chunk) % Threads.nthreads()) + 1
            iso_splines = getIsoSplines(search_data[thread_id])

            for i in chunk
                
                prec_id = chroms[i,:precursor_idx]
                mz = prec_mz[prec_id]
                charge = prec_charge[prec_id]
                sulfur = sulfur_count[prec_id]
                scan_id = scan_idx[i]
                scan_mz = coalesce(centerMz[scan_id], zero(Float32))::Float32
                window_width = coalesce(isolationWidthMz[scan_id], zero(Float32))::Float32

                #Needs to be based on the scan definition and not the fitted model
                #because the isotopes_captured annotation must be consistent between runs 
                low_mz, high_mz = Float32(scan_mz - window_width/2), Float32(scan_mz + window_width/2)
                isotopes = getPrecursorIsotopeSet(mz, charge, low_mz, high_mz)
                ## get precursor transmission precent
                precursor_fraction_transmitted[i] = getPrecursorFractionTransmitted!(
                    iso_splines,
                    (1,5),
                    getQuadTransmissionFunction(quad_transmission_model, scan_mz, window_width),
                    mz,
                    charge,
                    sulfur)

                isotopes_captured[i] = isotopes
            end
        end
    end
    fetch.(tasks)
    chroms[!,:isotopes_captured] = isotopes_captured
    chroms[!,:precursor_fraction_transmitted] = precursor_fraction_transmitted
    return nothing
end


"""
    add_features!(psms::DataFrame, search_context::SearchContext, ...)

Add feature columns to PSMs for scoring and analysis.

# Added Features
- RT and iRT metrics
- Sequence properties
- Intensity metrics
- Spectrum characteristics
"""
function add_features!(psms::DataFrame, 
                        search_context::SearchContext,
                                    tic::AbstractVector{Float32},
                                    masses::AbstractArray,
                                    ms_file_idx::Integer,
                                    rt_to_irt_interp::RtConversionModel,
                                    prec_id_to_irt::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}}
                                    )

    precursor_sequence = getSequence(getPrecursors(getSpecLib(search_context)))#[:sequence],
    structural_mods = getStructuralMods(getPrecursors(getSpecLib(search_context)))#[:structural_mods],
    prec_mz = getMz(getPrecursors(getSpecLib(search_context)))#[:mz],
    prec_irt = getIrt(getPrecursors(getSpecLib(search_context)))#[:irt],
    prec_charge = getCharge(getPrecursors(getSpecLib(search_context)))#[:prec_charge],
    entrap_group_ids = getEntrapmentGroupId(getPrecursors(getSpecLib(search_context)))
    precursor_missed_cleavage = getMissedCleavages(getPrecursors(getSpecLib(search_context)))#[:missed_cleavages],
    precursor_pair_idxs = getPairIdx(getPrecursors(getSpecLib(search_context)))
    #filter!(x -> x.best_scan, psms);
    filter!(x->x.weight>0, psms);
    #filter!(x->x.data_points>0, psms)
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    irt_diff = zeros(Float32, N)
    irt_obs = zeros(Float32, N)
    ms1_irt_diff = zeros(Float32, N)
    irt_pred = zeros(Float32, N)
    irt_error = zeros(Float32, N)
    pair_idxs = zeros(UInt32, N)
    entrap_group_id = zeros(UInt8, N)
    #psms[!,:missed_cleavage] .= zero(UInt8);
    #psms[!,:sequence] .= "";
    missed_cleavage = zeros(UInt8, N);
    #sequence = Vector{String}(undef, N);
    #stripped_sequence = Vector{String}(undef, N);
    adjusted_intensity_explained = zeros(Float16, N);
    #log2_base_peak_intensity = zeros(Float16, N);
    prec_charges = zeros(UInt8, N)
    sequence_length = zeros(UInt8, N);
    #b_y_overlap = zeros(Bool, N);
    spectrum_peak_count = zeros(Float16, N);
    prec_mzs = zeros(Float32, size(psms, 1));
    Mox = zeros(UInt8, N);
    TIC = zeros(Float16, N);

    # Amino acid count features (20 standard amino acids)
    aa_counts = Dict{Char, Vector{UInt8}}()
    for aa in "ACDEFGHIKLMNPQRSTVWY"
        aa_counts[aa] = zeros(UInt8, N)
    end

    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    #masses = MS_TABLE[:mz_array]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    #longest_y::Vector{UInt8} = psms[!,:longest_y]
    #longest_b::Vector{UInt8} = psms[!,:longest_b]
    rt::Vector{Float32} = psms[!,:rt]
    ms1_rt::Vector{Float32} = psms[!,:rt_ms1]
    ms1_missing::Vector{Bool} = psms[!,:ms1_features_missing]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}
    #precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    function countMOX(seq::Missing)
        return zero(UInt8)
    end


    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin 
            for i in chunk
                prec_idx = precursor_idx[i]
                entrap_group_id[i] = entrap_group_ids[prec_idx]
                irt_obs[i] = rt_to_irt_interp(rt[i])
                irt_pred[i] = getPredIrt(search_context, prec_idx)#prec_irt[prec_idx]
                #irt_diff[i] = abs(irt_obs[i] - first(prec_id_to_irt[prec_idx]))
                irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_irt)
                if !ms1_missing[i]
                    ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx))
                else
                    ms1_irt_diff[i] = 0f0
                end
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])

                missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
                #sequence[i] = precursor_sequence[prec_idx]
                stripped_seq = replace(precursor_sequence[prec_idx], r"\(.*?\)" => "")
                sequence_length[i] = length(stripped_seq)
                for ch in stripped_seq
                    if haskey(aa_counts, ch)
                        aa_counts[ch][i] += one(UInt8)
                    end
                end
                Mox[i] = countMOX(structural_mods[prec_idx])::UInt8
                #sequence_length[i] = length(stripped_sequence[i])
                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                prec_charges[i] = prec_charge[prec_idx]
                #b_y_overlap[i] = ((sequence_length[i] - longest_y[i])>longest_b[i]) &  (longest_b[i] > 0) & (longest_y[i] > 0);
                pair_idxs[i] = extract_pair_idx(precursor_pair_idxs, prec_idx)
                spectrum_peak_count[i] = length(masses[scan_idx[i]])
         
                prec_mzs[i] = prec_mz[prec_idx];
            end
        end
    end
    fetch.(tasks)

    # Build has_spectral_pair: true if the paired target/decoy precursor also has a PSM in this file
    prec_idx_set = Set(precursor_idx)
    has_spectral_pair = BitVector(undef, N)
    for i in 1:N
        has_spectral_pair[i] = pair_idxs[i] != zero(UInt32) && pair_idxs[i] in prec_idx_set
    end
    n_paired = count(has_spectral_pair)
    pct_paired = round(100*n_paired/N, digits=1)
    @info "Spectral pair coverage: $n_paired / $N PSMs ($pct_paired%) have paired target/decoy in file\n"
    psms[!, :has_spectral_pair] = has_spectral_pair

    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_diff] = irt_diff
    psms[!,:irt_error] = irt_error
    psms[!,:ms1_irt_diff] = ms1_irt_diff
    psms[!,:missed_cleavage] = missed_cleavage
    #psms[!,:sequence] = sequence
    #psms[!,:stripped_sequence] = stripped_sequence
    psms[!,:Mox] = Mox
    psms[!,:sequence_length] = sequence_length

    psms[!,:tic] = TIC
    psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
    psms[!,:charge] = prec_charges
    
    #psms[!,:b_y_overlap] = b_y_overlap
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:pair_id] = pair_idxs
    psms[!,:prec_mz] = prec_mzs
    psms[!,:entrapment_group_id] = entrap_group_id
    psms[!,:ms_file_idx] .= ms_file_idx

    # Amino acid count columns
    for aa in "ACDEFGHIKLMNPQRSTVWY"
        psms[!, Symbol("aa_", aa)] = aa_counts[aa]
    end
    return nothing
end

"""
    add_precursor_ms2_features!(psms, spectra, search_context, ms_file_idx)

Check if the monoisotopic precursor m/z appears as a peak in each MS2 scan.
Adds column `log2_prec_ms2_intensity` (0.0 when not found).
"""
function add_precursor_ms2_features!(psms::DataFrame, spectra::MassSpecData, search_context::SearchContext, ms_file_idx::Int64)
    prec_mzs = getMz(getPrecursors(getSpecLib(search_context)))
    mass_err_model = getMassErrorModel(search_context, ms_file_idx)
    mz_arrays = getMzArrays(spectra)
    int_arrays = getIntensityArrays(spectra)

    N = nrow(psms)
    log2_intensity = zeros(Float32, N)

    precursor_idx = psms[!, :precursor_idx]
    scan_idx = psms[!, :scan_idx]

    # getMzBounds assumes empirical m/z is already corrected (see MassErrorModel.jl).
    # We must apply getCorrectedMz to each raw peak before comparing, matching matchPeaks!.
    # Use a wider binary-search window to account for raw-vs-corrected offset.
    offset = abs(getMassCorrection(mass_err_model))

    for i in 1:N
        prec_mz = prec_mzs[precursor_idx[i]]
        mz_arr = mz_arrays[scan_idx[i]]
        int_arr = int_arrays[scan_idx[i]]

        low, high = getMzBounds(mass_err_model, prec_mz)
        # Widen binary search to cover raw (uncorrected) m/z values
        ppm_offset = offset * (prec_mz / 1f6)
        lo_idx = searchsortedfirst(mz_arr, low - ppm_offset)
        hi_idx = searchsortedlast(mz_arr, high + ppm_offset)

        if lo_idx <= hi_idx && lo_idx <= length(mz_arr)
            best_int = Float32(0)
            for j in lo_idx:hi_idx
                corrected_mz = getCorrectedMz(mass_err_model, mz_arr[j])
                if corrected_mz >= low && corrected_mz <= high
                    peak_int = Float32(coalesce(int_arr[j], 0f0))
                    best_int = max(best_int, peak_int)
                end
            end
            if best_int > 0
                log2_intensity[i] = log2(best_int)
            end
        end
    end

    psms[!, :log2_prec_ms2_intensity] = log2_intensity

    n_found = count(>(0), log2_intensity)
    @info "  Precursor MS2 features: $(n_found)/$(N) ($(round(100*n_found/N, digits=1))%) precursor m/z peaks found in MS2 scans"
    if n_found > 0
        found_mask = log2_intensity .> 0
        vals = log2_intensity[found_mask]
        @info "    log2_prec_ms2_intensity range: $(round(minimum(vals), digits=1)) to $(round(maximum(vals), digits=1)), mean=$(round(sum(vals)/length(vals), digits=1))"
    end
end

#==========================================================
Summary Statistics
==========================================================#
"""
    init_summary_columns!(psms::DataFrame)

Initialize columns for summary statistics across PSM groups.

# Added Columns
- Maximum entropy, goodness-of-fit
- Peak intensity metrics
- Ion coverage statistics
"""
function init_summary_columns!(
    psms::DataFrame,
    )

    new_cols = [
        (:max_entropy,              Float16)
        (:max_gof,         Float16)
        (:max_fitted_manhattan_distance,          Float16)
        (:max_fitted_spectral_contrast,         Float16)
        (:max_scribe, Float16)
        (:y_ions_sum,               UInt16)
        (:max_y_ions,               UInt16)
        (:max_matched_ratio,        Float16)
        (:num_scans,        UInt16)
        (:smoothness,        Float32)
        (:weights,        Vector{Float32})
        (:irts,         Vector{Float32})
        ];

        N = size(psms, 1)
        for (col_name, col_type) in new_cols
            if col_type <: AbstractVector
                # col_type is something like Vector{Float16};
                # create an array of length N, each element is an empty vector of that subtype.
                psms[!, col_name] = [col_type() for _ in 1:N]
            else
                # scalar numeric type, use zeros as before
                psms[!, col_name] = zeros(col_type, N)
            end
        end
        return psms
end

#==========================================================
DIA-NN Recovery Diagnostics
==========================================================#

const DIANN_PERFILE_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/DIA-NN_results/per_file_precursor_indices"
const DIANN_GLOBAL_ARROW = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/DIA-NN_results/diann_precursor_indices_3P_rank_based.arrow"

"""
    load_diann_reference(file_name::String) -> (per_file::Set{UInt32}, global_set::Set{UInt32})

Load DIA-NN reference precursor indices for diagnostics.
Returns empty sets if files are not found (non-fatal).
"""
function load_diann_reference(file_name::String)
    per_file = Set{UInt32}()
    global_set = Set{UInt32}()

    # Per-file reference
    pf_path = joinpath(DIANN_PERFILE_DIR, "$(file_name).arrow")
    if isfile(pf_path)
        t = Arrow.Table(pf_path)
        per_file = Set{UInt32}(t.precursor_idx)
        @info "  DIA-NN per-file ($file_name): $(length(per_file)) unique precursor_idx values\n"
    end

    # Global reference
    if isfile(DIANN_GLOBAL_ARROW)
        t = Arrow.Table(DIANN_GLOBAL_ARROW)
        global_set = Set{UInt32}(t.precursor_idx)
        @info "  DIA-NN global: $(length(global_set)) unique precursor_idx values\n"
    end

    return per_file, global_set
end

"""
    log_diann_recovery(label, psms, diann_file, diann_global)

Log how many DIA-NN precursors have at least one surviving target PSM.
"""
function log_diann_recovery(label::String, psms::DataFrame,
                            diann_file::Set{UInt32}, diann_global::Set{UInt32})
    (isempty(diann_file) && isempty(diann_global)) && return

    # Unique target precursor_idxs with surviving PSMs
    target_mask = psms[!, :target]
    surviving_precs = Set{UInt32}()
    prec_col = psms[!, :precursor_idx]
    for i in 1:nrow(psms)
        if target_mask[i]
            push!(surviving_precs, prec_col[i])
        end
    end

    if !isempty(diann_file)
        n_recovered = length(intersect(surviving_precs, diann_file))
        n_total = length(diann_file)
        pct = round(100.0 * n_recovered / max(1, n_total), digits=1)
        @info "  [$label] DIA-NN file recovery: $n_recovered / $n_total ($pct%)\n"
    end
    if !isempty(diann_global)
        n_recovered = length(intersect(surviving_precs, diann_global))
        n_total = length(diann_global)
        pct = round(100.0 * n_recovered / max(1, n_total), digits=1)
        @info "  [$label] DIA-NN global recovery: $n_recovered / $n_total ($pct%)\n"
    end
end

#==========================================================
Iterative Prescore Filter (DIA-NN-style diff-cov + probit)
==========================================================#

"""
Per-scan features for the iterative prescore filter.
These are all computed within a single scan (no cross-scan summaries).
"""
const ITERATIVE_PRESCORE_FEATURES = [
    :fitted_manhattan_distance,   # per-scan spectral distance
    :max_matched_residual,        # per-scan: largest matched fragment residual
    :gof,                         # per-scan goodness-of-fit
    :max_unmatched_residual,      # per-scan: largest unmatched fragment residual
    :poisson,                     # per-scan Poisson score
    :irt_error,                   # per-scan: |observed iRT - predicted iRT|
    :y_count,                     # per-scan: number of y-ions matched
    :scribe,                      # per-scan: scribe score
    :err_norm,                    # per-scan: normalized mass error
    :spectral_contrast,           # per-scan: cosine similarity
]

const LGBM_RECOVERY_FEATURES = [
    # Core spectral quality
    :fitted_manhattan_distance, :max_matched_residual, :gof,
    :max_unmatched_residual, :poisson, :spectral_contrast, :err_norm,
    # Scores / weights
    :scribe, :weight, :log2_intensity_explained,
    # Ion counts
    :y_count, :b_count, :isotope_count, :total_ions,
    # RT
    :irt_error, :irt_diff,
    # Peptide properties
    :charge, :sequence_length, :missed_cleavage, :Mox, :prec_mz,
    # Other
    :tic, :best_rank, :matched_ratio, :spectrum_peak_count,
    # Amino acid composition
    :aa_H, :aa_P, :aa_L,
    # Precursor in MS2
    :log2_prec_ms2_intensity,
]

"""
    ensure_ms1_stub_columns!(psms::DataFrame)

Add placeholder MS1 columns needed by add_features!() during iterative
prescore filtering, before the real MS1 join in process_search_results!.
"""
function ensure_ms1_stub_columns!(psms::DataFrame)
    n = nrow(psms)
    if !hasproperty(psms, :rt_ms1)
        psms[!, :rt_ms1] = fill(Float32(-1), n)
    end
    if !hasproperty(psms, :ms1_features_missing)
        psms[!, :ms1_features_missing] = trues(n)
    end
    if !hasproperty(psms, :ms1_ms2_rt_diff)
        psms[!, :ms1_ms2_rt_diff] = fill(Float32(-1), n)
    end
end

"""
    sanitize_prescore_features!(psms::DataFrame, features::Vector{Symbol})

Replace Inf/NaN values with 0.0 in prescore feature columns.
Float16 columns (e.g., err_norm) can overflow to Inf, which would
crash the covariance matrix solve in train_difference_covariance!.
"""
function sanitize_prescore_features!(psms::DataFrame, features::Vector{Symbol})
    n_fixed = 0
    for f in features
        hasproperty(psms, f) || continue
        col = psms[!, f]
        for i in eachindex(col)
            if !isfinite(col[i])
                col[i] = zero(eltype(col))
                n_fixed += 1
            end
        end
    end
    if n_fixed > 0
        @info "  Sanitized $n_fixed non-finite feature values (Inf/NaN → 0)\n"
    end
end

"""
    get_best_psm_per_precursor(psms::DataFrame) -> DataFrame

Reduce to one PSM per precursor_idx, keeping the row with highest `weight`.
Returns a copy (does not modify input).
"""
function get_best_psm_per_precursor(psms::DataFrame)
    sorted = sort(psms, [:precursor_idx, order(:weight, rev=true)])
    best = combine(groupby(sorted, :precursor_idx), first)
    return best
end

"""
    get_best_psm_per_precursor_by_score(psms::DataFrame, score_col::Symbol) -> DataFrame

Reduce to one PSM per precursor_idx, keeping the row with highest `score_col`.
Returns a copy (does not modify input).
"""
function get_best_psm_per_precursor_by_score(psms::DataFrame, score_col::Symbol)
    sorted = sort(psms, [:precursor_idx, order(score_col, rev=true)])
    best = combine(groupby(sorted, :precursor_idx), first)
    return best
end

"""
    form_target_decoy_pairs(best_psms::DataFrame) -> Vector{Tuple{Int,Int}}

Find matched target-decoy pairs from best-PSM-per-precursor DataFrame.
Each pair is (target_row_idx, decoy_row_idx).
Uses the `pair_id` column which holds the partner's precursor_idx.
"""
function form_target_decoy_pairs(best_psms::DataFrame)
    prec_to_row = Dict{UInt32, Int}()
    for i in 1:nrow(best_psms)
        prec_to_row[best_psms[i, :precursor_idx]] = i
    end

    pairs = Vector{Tuple{Int,Int}}()
    visited = falses(nrow(best_psms))

    for i in 1:nrow(best_psms)
        visited[i] && continue
        partner_pid = best_psms[i, :pair_id]
        partner_pid == zero(UInt32) && continue

        partner_row = get(prec_to_row, partner_pid, 0)
        partner_row == 0 && continue
        visited[partner_row] && continue

        is_target_i = best_psms[i, :target]
        is_target_j = best_psms[partner_row, :target]
        (is_target_i == is_target_j) && continue

        t_row = is_target_i ? i : partner_row
        d_row = is_target_i ? partner_row : i

        push!(pairs, (t_row, d_row))
        visited[i] = true
        visited[partner_row] = true
    end

    return pairs
end

"""
    train_difference_covariance!(best_psms, pairs, features; regularization) -> Union{Vector{Float64}, Nothing}

DIA-NN-style difference-covariance classifier.
  delta_k = features_target[k] - features_decoy[k]
  w = (cov(delta) + eps*I) \\ mean(delta)
"""
function train_difference_covariance!(
    best_psms::DataFrame,
    pairs::Vector{Tuple{Int,Int}},
    features::Vector{Symbol};
    regularization::Float64 = 1e-9
)
    N = length(pairs)
    p = length(features)

    if N < 20
        return nothing
    end

    # Mean feature differences (target - decoy)
    ds_mean = zeros(Float64, p)
    for (t_row, d_row) in pairs
        for j in 1:p
            ds_mean[j] += Float64(best_psms[t_row, features[j]]) -
                          Float64(best_psms[d_row, features[j]])
        end
    end
    ds_mean ./= N

    # Covariance of differences
    A = zeros(Float64, p, p)
    for (t_row, d_row) in pairs
        for i in 1:p
            di = (Float64(best_psms[t_row, features[i]]) -
                  Float64(best_psms[d_row, features[i]])) - ds_mean[i]
            for j in i:p
                dj = (Float64(best_psms[t_row, features[j]]) -
                      Float64(best_psms[d_row, features[j]])) - ds_mean[j]
                A[i,j] += di * dj
            end
        end
    end
    for i in 1:p, j in i:p
        A[i,j] /= (N - 1)
        A[j,i] = A[i,j]
    end
    for i in 1:p
        A[i,i] += regularization
    end

    w = A \ ds_mean

    @info "    Diff-cov weights ($N pairs, $p features):\n"
    for (fname, wval) in zip(features, w)
        @info "      $(rpad(fname, 35)) $(round(wval, digits=6))\n"
    end

    return w
end

"""
    apply_linear_scores!(psms, w, features)

Score all PSMs: score_i = w' * features_i. Stores in psms[!, :score].
"""
function apply_linear_scores!(psms::DataFrame, w::Vector{Float64}, features::Vector{Symbol})
    n = nrow(psms)
    scores = zeros(Float32, n)
    for i in 1:n
        s = 0.0
        for (j, f) in enumerate(features)
            s += w[j] * Float64(psms[i, f])
        end
        scores[i] = Float32(s)
    end
    psms[!, :score] = scores
end

"""
    train_probit_on_features!(psms, train_mask, features; max_iter) -> Vector{Float64}

Train probit model on the masked subset. Returns coefficient vector.
"""
function train_probit_on_features!(
    psms::DataFrame,
    train_mask::BitVector,
    features::Vector{Symbol};
    max_iter::Int = 20
)
    X_train = DataFrame(Matrix{Float64}(psms[train_mask, features]), features)
    targets_train = Vector{Bool}(psms[train_mask, :target])
    n_train = sum(train_mask)

    chunk_size = max(1, n_train ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n_train, chunk_size)

    beta = zeros(Float64, length(features))
    beta = ProbitRegression(beta, X_train, targets_train, data_chunks; max_iter=max_iter)

    @info "    Probit coefficients:\n"
    for (fname, coef) in zip(features, beta)
        @info "      $(rpad(fname, 35)) $(round(coef, digits=6))\n"
    end

    return beta
end

"""
    apply_probit_scores!(psms, beta, features)

Apply probit probability scores to all PSMs. Stores in psms[!, :probit_score].
"""
function apply_probit_scores!(
    psms::DataFrame,
    beta::Vector{Float64},
    features::Vector{Symbol}
)
    n = nrow(psms)
    X_all = DataFrame(Matrix{Float64}(psms[!, features]), features)
    scores = zeros(Float32, n)

    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n, chunk_size)

    ModelPredictProbs!(scores, X_all, beta, data_chunks)
    psms[!, :probit_score] = scores
end

"""
    rerun_search_with_filter(spectra, search_context, params, ms_file_idx, surviving_pairs)

Rebuild scan_to_prec_idx with only surviving (precursor_idx, scan_idx) pairs,
then re-run perform_second_pass_search. This re-solves the Huber deconvolution
with reduced precursor competition per scan.
"""
function rerun_search_with_filter(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    surviving_pairs::Set{Tuple{UInt32, UInt32}}
)
    frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
    scan_to_prec_idx_orig, precursors_passed_orig = load_fragment_index_matches(
        frag_match_path, length(spectra)
    )
    n_original = length(precursors_passed_orig)

    new_precursors = UInt32[]
    new_scan_to_prec = Vector{Union{Missing, UnitRange{Int64}}}(missing, length(spectra))

    for scan_idx in 1:length(spectra)
        range = scan_to_prec_idx_orig[scan_idx]
        ismissing(range) && continue

        start_new = length(new_precursors) + 1
        for idx in range
            pid = precursors_passed_orig[idx]
            if (pid, UInt32(scan_idx)) in surviving_pairs
                push!(new_precursors, pid)
            end
        end
        end_new = length(new_precursors)

        if end_new >= start_new
            new_scan_to_prec[scan_idx] = start_new:end_new
        end
    end

    n_filtered = length(new_precursors)
    pct = round(100.0 * n_filtered / max(1, n_original), digits=1)
    @info "    Re-search input: $n_filtered / $n_original entries ($pct%)\n"

    return perform_second_pass_search(
        spectra, new_scan_to_prec, new_precursors,
        search_context, params, ms_file_idx, MS2CHROM()
    )
end

"""
    rerun_search_with_precursor_filter(spectra, search_context, params, ms_file_idx, passing_precs)

Like `rerun_search_with_filter`, but filters by precursor identity only (not per-scan).
Keeps all scans for precursors in `passing_precs`, then re-solves deconvolution
with reduced competition.
"""
function rerun_search_with_precursor_filter(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    passing_precs::Set{UInt32}
)
    frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
    scan_to_prec_idx_orig, precursors_passed_orig = load_fragment_index_matches(
        frag_match_path, length(spectra)
    )
    n_original = length(precursors_passed_orig)

    new_precursors = UInt32[]
    new_scan_to_prec = Vector{Union{Missing, UnitRange{Int64}}}(missing, length(spectra))

    for scan_idx in 1:length(spectra)
        range = scan_to_prec_idx_orig[scan_idx]
        ismissing(range) && continue

        start_new = length(new_precursors) + 1
        for idx in range
            pid = precursors_passed_orig[idx]
            if pid in passing_precs
                push!(new_precursors, pid)
            end
        end
        end_new = length(new_precursors)

        if end_new >= start_new
            new_scan_to_prec[scan_idx] = start_new:end_new
        end
    end

    n_filtered = length(new_precursors)
    pct = round(100.0 * n_filtered / max(1, n_original), digits=1)
    @info "    Global re-search input: $n_filtered / $n_original entries ($pct%)\n"

    return perform_second_pass_search(
        spectra, new_scan_to_prec, new_precursors,
        search_context, params, ms_file_idx, MS2CHROM()
    )
end

"""
    iterative_prescore_filter!(psms, search_context, spectra, params, ms_file_idx)

Single-round prescore filter with low-score retrain:
  1. Two-pass probit on best-per-precursor PSMs (like FirstPassSearch)
  2. Retrain a second probit on q>0.10 subset to recover borderline precursors
  3. Exclude precursors failing both models at q>0.10
  4. Re-run search with surviving (precursor_idx, scan_idx) pairs

Returns fresh PSMs from the re-search (raw columns only,
ready for process_search_results!).
"""
function iterative_prescore_filter!(
    psms::DataFrame,
    search_context::SearchContext,
    spectra::MassSpecData,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    n_initial = nrow(psms)
    @info "=== Prescore Filter: starting with $n_initial PSMs ===\n"

    # Load DIA-NN reference sets for recovery diagnostics
    file_name = getParsedFileName(search_context, ms_file_idx)
    diann_file, diann_global = load_diann_reference(file_name)
    if !isempty(diann_file) || !isempty(diann_global)
        @info "  DIA-NN reference loaded: $(length(diann_file)) per-file, $(length(diann_global)) global precursors\n"
    end

    # ── Step 1: Add columns + quality filter ──────────────────────
    add_second_search_columns!(psms,
        getRetentionTimes(spectra),
        getCharge(getPrecursors(getSpecLib(search_context))),
        getIsDecoy(getPrecursors(getSpecLib(search_context))),
        getPrecursors(getSpecLib(search_context))
    )

    get_isotopes_captured!(psms,
        getIsotopeTraceType(params),
        getQuadTransmissionModel(search_context, ms_file_idx),
        getSearchData(search_context),
        psms[!, :scan_idx],
        getCharge(getPrecursors(getSpecLib(search_context))),
        getMz(getPrecursors(getSpecLib(search_context))),
        getSulfurCount(getPrecursors(getSpecLib(search_context))),
        getCenterMzs(spectra),
        getIsolationWidthMzs(spectra)
    )

    filter!(row -> row.precursor_fraction_transmitted >= params.min_fraction_transmitted, psms)
    filter!(row -> row.weight > 1e-6, psms)
    @info "  $(nrow(psms)) PSMs after quality filters\n"

    # ── Step 2: Add ML features ───────────────────────────────────
    ensure_ms1_stub_columns!(psms)
    add_features!(psms, search_context,
        getTICs(spectra), getMzArrays(spectra),
        ms_file_idx,
        getRtIrtModel(search_context, ms_file_idx),
        getPrecursorDict(search_context)
    )
    add_precursor_ms2_features!(psms, spectra, search_context, ms_file_idx)

    sanitize_prescore_features!(psms, collect(ITERATIVE_PRESCORE_FEATURES))
    features = filter(f -> hasproperty(psms, f), collect(ITERATIVE_PRESCORE_FEATURES))
    @info "  Using $(length(features)) features: $(join(features, ", "))\n"

    # ── Step 3a: Two-pass probit on ALL PSMs ─────────────────────
    targets_all = psms[!, :target]
    n_all = nrow(psms)
    @info "  All-PSM probit: $n_all PSMs ($(count(targets_all)) targets)\n"

    # Pass 1: seed from scribe q≤0.01 on all PSMs + all decoys
    scribe_q_all = zeros(Float64, n_all)
    get_qvalues!(psms[!, :scribe], targets_all, scribe_q_all)
    pass1_mask_all = BitVector(((scribe_q_all .<= 0.01) .& targets_all) .| .!targets_all)
    n_p1_t_all = count((scribe_q_all .<= 0.01) .& targets_all)
    n_p1_d_all = count(.!targets_all)
    @info "  All-PSM probit pass 1: $n_p1_t_all targets (scribe q≤0.01) + $n_p1_d_all decoys\n"

    beta_all = train_probit_on_features!(psms, pass1_mask_all, features)
    apply_probit_scores!(psms, beta_all, features)

    # Pass 2: retrain on probit q≤0.01 + all decoys
    q_all_p1 = zeros(Float64, n_all)
    get_qvalues!(psms[!, :probit_score], targets_all, q_all_p1)
    n_1pct_all = count((q_all_p1 .<= 0.01) .& targets_all)
    @info "  All-PSM probit pass 1: $n_1pct_all targets @ 1% FDR\n"

    pass2_mask_all = BitVector(((q_all_p1 .<= 0.01) .& targets_all) .| .!targets_all)
    beta_all = train_probit_on_features!(psms, pass2_mask_all, features)
    apply_probit_scores!(psms, beta_all, features)

    # ── Step 3b: Best PSM per precursor (by probit_score) ──────
    best_psms = get_best_psm_per_precursor_by_score(psms, :probit_score)
    @info "  $(nrow(best_psms)) unique precursors from $(nrow(psms)) PSMs (selected by probit_score)\n"

    # ── Step 4: LightGBM on best PSMs ──────────────────────────
    targets_col = best_psms[!, :target]

    lgbm_main_features = filter(f -> hasproperty(best_psms, f), collect(LGBM_RECOVERY_FEATURES))
    @info "  LightGBM main model: $(nrow(best_psms)) PSMs, $(length(lgbm_main_features)) features\n"

    lgbm_main_feature_df = best_psms[!, lgbm_main_features]

    classifier_main = build_lightgbm_classifier(
        num_iterations = 200,
        max_depth = 5,
        num_leaves = 31,
        min_data_in_leaf = 50,
        feature_fraction = 0.8,
        bagging_fraction = 0.8,
        bagging_freq = 1,
        is_unbalance = true,
    )
    lgbm_main_model = fit_lightgbm_model(classifier_main, lgbm_main_feature_df, targets_col)

    lgbm_main_scores = lightgbm_predict(lgbm_main_model, lgbm_main_feature_df)
    q_main = zeros(Float64, nrow(best_psms))
    get_qvalues!(lgbm_main_scores, targets_col, q_main)

    n_pass_1 = count((q_main .<= 0.01) .& targets_col)
    n_pass_5 = count((q_main .<= 0.05) .& targets_col)
    n_pass_10 = count((q_main .<= 0.10) .& targets_col)
    @info "  Main LightGBM: $n_pass_1 @ 1%, $n_pass_5 @ 5%, $n_pass_10 @ 10% FDR target precursors\n"

    # Feature importances
    imp_main = importance(lgbm_main_model)
    if imp_main !== nothing
        sorted_imp_main = sort(imp_main, by = x -> -x[2])
        @info "  LightGBM main model feature importances ($(length(sorted_imp_main)) features):"
        for (fname, gain) in sorted_imp_main
            @info "    $fname: $(round(gain, digits=1))"
        end
        @info ""
    end

    # DIA-NN diagnostics
    if !isempty(diann_file) || !isempty(diann_global)
        for q_thresh in [0.01, 0.05, 0.10]
            lgbm_precs = Set{UInt32}(best_psms[i, :precursor_idx]
                                      for i in 1:nrow(best_psms) if (q_main[i] <= q_thresh) && targets_col[i])
            n_file = isempty(diann_file) ? 0 : length(intersect(lgbm_precs, diann_file))
            n_global = isempty(diann_global) ? 0 : length(intersect(lgbm_precs, diann_global))
            @info "    LightGBM main q≤$q_thresh: $(length(lgbm_precs)) targets, DIA-NN file=$n_file global=$n_global\n"
        end
    end

    # ── Step 4b: Save per-file LightGBM scores for cross-run aggregation ──
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
    mkpath(prescore_dir)
    score_df = DataFrame(
        precursor_idx = best_psms[!, :precursor_idx],
        lgbm_prob = Float32.(lgbm_main_scores),
        target = best_psms[!, :target]
    )
    writeArrow(joinpath(prescore_dir, "$(file_name).arrow"), score_df)

    # ── Step 5: Collect passing precursors ────────────────────────
    passing_precs = Set{UInt32}(best_psms[i, :precursor_idx]
                                 for i in 1:nrow(best_psms) if q_main[i] <= params.prescore_qvalue_threshold)

    n_passing = length(passing_precs)
    @info "  $n_passing precursors survive (main model)\n"

    # ── Step 7: Collect surviving (precursor_idx, scan_idx) pairs ─
    # From the ORIGINAL all-PSMs DataFrame (not just best-per-precursor)
    surviving_pairs = Set{Tuple{UInt32, UInt32}}()
    prec_col = psms[!, :precursor_idx]
    scan_col = psms[!, :scan_idx]
    for i in 1:nrow(psms)
        if prec_col[i] in passing_precs
            push!(surviving_pairs, (prec_col[i], scan_col[i]))
        end
    end
    n_surviving_precs = length(Set(p[1] for p in surviving_pairs))
    surviving_prec_set = Set(p[1] for p in surviving_pairs)
    n_diann_file_present = isempty(diann_file) ? 0 : length(intersect(surviving_prec_set, diann_file))
    n_diann_file_absent  = isempty(diann_file) ? 0 : length(setdiff(diann_file, surviving_prec_set))
    n_diann_global_present = isempty(diann_global) ? 0 : length(intersect(surviving_prec_set, diann_global))
    n_diann_global_absent  = isempty(diann_global) ? 0 : length(setdiff(diann_global, surviving_prec_set))
    @info "  $(length(surviving_pairs)) surviving (prec, scan) pairs across $n_surviving_precs precursors\n"
    @info "  DIA-NN file: $n_diann_file_present present, $n_diann_file_absent absent | global: $n_diann_global_present present, $n_diann_global_absent absent\n"

    # ── Step 8: Re-run search with filtered pairs ─────────────────
    @info "  Re-running search with filtered pairs...\n"
    psms = rerun_search_with_filter(
        spectra, search_context, params, ms_file_idx, surviving_pairs
    )

    n_final = nrow(psms)
    pct_reduction = round(100.0 * (1 - n_final / n_initial), digits=1)
    @info "=== Prescore Filter complete: $n_initial -> $n_final PSMs ($pct_reduction% reduction) ===\n"

    # Final recovery check — psms here are raw from re-search, need target column
    if !isempty(diann_file) || !isempty(diann_global)
        is_decoy = getIsDecoy(getPrecursors(getSpecLib(search_context)))
        prec_col_final = psms[!, :precursor_idx]
        tmp_target = BitVector(undef, nrow(psms))
        for i in 1:nrow(psms)
            tmp_target[i] = !is_decoy[prec_col_final[i]]
        end
        psms[!, :target] = tmp_target
        log_diann_recovery("Final PSMs (to ScoringSearch)", psms, diann_file, diann_global)
        select!(psms, Not(:target))
    end

    return psms
end

#==========================================================
Prescore (Best Scan Selection)
==========================================================#

const PRESCORE_FEATURES = [
    :fitted_spectral_contrast,
    :err_norm,
    :log2_intensity_explained,
    :gof,
    :scribe,
    :weight,
]

"""
    train_and_apply_prescore!(psms::DataFrame)

Train a lightweight probit model on per-scan features to produce a preliminary
quality score (`prescore`) for each scan. Then for each precursor group, mark
the scan with the highest prescore as `best_scan = true`.

This replaces the default `argmax(weight)` apex selection with a multi-feature
quality-based selection.
"""
function train_and_apply_prescore!(
    psms::DataFrame;
    features::Vector{Symbol} = PRESCORE_FEATURES,
    n_train_rounds::Int = 2,
    max_iter::Int = 20
)
    n = nrow(psms)
    if n < 100
        @user_warn "Too few PSMs ($n) for probit prescore — falling back to weight-based apex selection"
        return false
    end

    @user_info "Training probit prescore model on $n scans with features: $(join(features, ", "))"

    # Prepare labels
    targets = psms[!, :target]
    n_targets = sum(targets)
    n_decoys = n - n_targets
    @user_info "  Prescore training data: $n_targets targets, $n_decoys decoys"

    # Prepare feature matrix (convert to Float64 for probit)
    X = DataFrame(Matrix{Float64}(psms[!, features]), features)



    # Partition for parallel IRLS
    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n, chunk_size)

    # Train probit model iteratively
    β = zeros(Float64, length(features))
    scores = zeros(Float32, n)

    for iround in 1:n_train_rounds
        if iround > 1
            # Compute q-values for training data selection
            q_vals = zeros(Float64, n)
            get_qvalues!(scores, targets, q_vals)

            # Filter training data: targets at q <= 0.05 + all decoys
            train_mask = (q_vals .<= 0.05 .&& targets) .|| .!targets
            n_train = sum(train_mask)
            X_train = X[train_mask, :]
            targets_train = targets[train_mask]
            chunks_train = Iterators.partition(1:n_train,
                                   max(1, n_train ÷ (10 * Threads.nthreads())))
            β = ProbitRegression(β, X_train, targets_train, chunks_train; max_iter=max_iter)
        else
            β = ProbitRegression(β, X, targets, data_chunks; max_iter=max_iter)
        end

        # Score all scans
        ModelPredict!(scores, X, β, data_chunks)
    end

    # Log probit coefficients (feature importances)
    @user_info "  Probit prescore coefficients:"
    for (fname, coef) in zip(features, β)
        @user_info "    $(rpad(fname, 30)) $(round(coef, digits=4))"
    end

    # Final scoring with probabilities
    ModelPredictProbs!(scores, X, β, data_chunks)

    # Store prescores
    psms[!, :prescore] = scores

    # Mark best_scan based on highest prescore per precursor group
    psms[!, :best_scan] .= false
    n_changed = 0
    n_groups = 0
    for gpsms in groupby(psms, :precursor_idx)
        n_groups += 1
        probit_best = argmax(gpsms[!, :prescore])
        weight_best = argmax(gpsms[!, :weight])
        if probit_best != weight_best
            n_changed += 1
        end
        gpsms[probit_best, :best_scan] = true
    end

    pct_changed = round(100*n_changed/max(n_groups,1), digits=1)
    @user_info "  Prescore selection changed apex in $n_changed / $n_groups precursor groups ($pct_changed%)"

    return true
end

"""
    get_summary_scores!(psms::SubDataFrame, weight::AbstractVector{Float32}, ...)

Calculate summary statistics for a group of related PSMs.

Computes maximum values and sums around apex scan for various
scoring metrics.
"""
function get_summary_scores!(
                            psms::SubDataFrame,
                            weight::AbstractVector{Float32},
                            gof::AbstractVector{Float16},
                            matched_ratio::AbstractVector{Float16},
                            #entropy::AbstractVector{Float16},
                            fitted_manhattan_distance::AbstractVector{Float16},
                            fitted_spectral_contrast::AbstractVector{Float16},
                            scribe::AbstractVector{Float16},
                            y_count::AbstractVector{UInt8},
                            rt_to_irt_interp::RtConversionModel
                        )

    max_gof = -100.0
    max_matched_ratio = -100.0
   # max_entropy = -100.0
    max_fitted_manhattan_distance = -100.0
    max_fitted_spectral_contrast= -100
    max_scribe = -100
    count = 0
    y_ions_sum = 0
    max_y_ions = 0
    smoothness = 0.0f0

    # Use probit-selected best_scan if available, otherwise fall back to max weight
    best_scan_col = psms[!, :best_scan]
    if any(best_scan_col)
        apex_scan = findfirst(best_scan_col)
    else
        apex_scan = argmax(psms[!, :weight])
    end
    #Need to make sure there is not a big gap.
    start = max(1, apex_scan - 2)
    stop = min(length(weight), apex_scan + 2)

    @inbounds @fastmath for i in range(start, stop)
        if gof[i]>max_gof
            max_gof =gof[i]
        end

        if matched_ratio[i]>max_matched_ratio
            max_matched_ratio = matched_ratio[i]
        end

        #if entropy[i]>max_entropy
        #    max_entropy=entropy[i]
        #end

        if fitted_manhattan_distance[i]>max_fitted_manhattan_distance
            max_fitted_manhattan_distance = fitted_manhattan_distance[i]
        end

        if fitted_spectral_contrast[i]>max_fitted_spectral_contrast
            max_fitted_spectral_contrast = fitted_spectral_contrast[i]
        end

        if scribe[i]>max_scribe
            max_scribe = scribe[i]
        end
    
        y_ions_sum += y_count[i]
        if y_count[i] > max_y_ions
            max_y_ions = y_count[i]
        end

        count += 1
    end    

    irts = rt_to_irt_interp.(psms.rt)
    
    @inbounds @fastmath for i in range(1, length(weight))
        if length(weight) == 1
            smoothness = (-2*weight[i] / weight[apex_scan])^2
        else
            if (i == 1)
                smoothness += (((weight[i+1] - weight[i]) / (psms.rt[i+1] - psms.rt[i]) + (-weight[i]) / (psms.rt[i+1] - psms.rt[i])) / weight[apex_scan]) ^2
            elseif (i > 1) & (i < length(weight))
                smoothness += (((weight[i-1] - weight[i]) / (psms.rt[i] - psms.rt[i-1]) + (weight[i+1]-weight[i]) / (psms.rt[i+1] - psms.rt[i])) / weight[apex_scan]) ^2
            elseif (i == length(weight))
                smoothness += (((weight[i-1] - weight[i]) / (psms.rt[i] - psms.rt[i-1]) + (-weight[i]) / (psms.rt[i] - psms.rt[i-1])) / weight[apex_scan]) ^2
            end
        end
    end

   

    psms.max_gof[apex_scan] = max_gof
    psms.max_matched_ratio[apex_scan] = max_matched_ratio
   # psms.max_entropy[apex_scan] = max_entropy
    psms.max_fitted_manhattan_distance[apex_scan] = max_fitted_manhattan_distance
    psms.max_fitted_spectral_contrast[apex_scan] = max_fitted_spectral_contrast
    psms.max_scribe[apex_scan] = max_scribe
    psms.y_ions_sum[apex_scan] = y_ions_sum
    psms.max_y_ions[apex_scan] = max_y_ions
    psms.num_scans[apex_scan] = length(weight)
    psms.smoothness[apex_scan] = smoothness
    psms.weights[apex_scan] = weight
    psms.irts[apex_scan] = irts
    psms.best_scan[apex_scan] = true

end

function parseMs1Psms(
    psms::DataFrame,
    spectra::MassSpecData,
    ms2_rt_lookup::Dict{UInt32, Float32}
)
    if !hasproperty(psms, :precursor_idx) || (size(psms, 1) == 0)
        return DataFrame()
    end

    # Add RT column
    rts = zeros(Float32, size(psms, 1))
    for i in range(1, size(psms, 1))
        scan_idx = psms[i, :scan_idx]
        rts[i] = getRetentionTime(spectra, scan_idx)
    end
    psms[!, :rt] = rts

    # Pre-allocate new columns for max intensity features
    psms[!, :rt_max_intensity] = zeros(Float32, size(psms, 1))
    psms[!, :rt_diff_max_intensity] = zeros(Float32, size(psms, 1))

    # Group by precursor and apply hybrid selection
    return combine(groupby(psms, :precursor_idx)) do group
        # Find max intensity scan
        max_intensity_idx = argmax(group[!, :weight])
        max_intensity_rt = group[max_intensity_idx, :rt]

        # Find RT-closest scan to MS2 apex
        precursor_idx = group[1, :precursor_idx]
        if haskey(ms2_rt_lookup, precursor_idx)
            ms2_rt = ms2_rt_lookup[precursor_idx]
            rt_diffs = abs.(group[!, :rt] .- ms2_rt)
            closest_rt_idx = argmin(rt_diffs)

            # Start with RT-closest scan features
            result_row = group[closest_rt_idx:closest_rt_idx, :]

            # Set max intensity RT features for this row
            result_row[1, :rt_max_intensity] = max_intensity_rt
            result_row[1, :rt_diff_max_intensity] = abs(max_intensity_rt - ms2_rt)

            return result_row
        else
            # Fallback to max intensity if no MS2 match
            fallback_row = group[max_intensity_idx:max_intensity_idx, :]
            # Set fallback values for new features
            fallback_row[1, :rt_max_intensity] = max_intensity_rt
            fallback_row[1, :rt_diff_max_intensity] = Float32(0.0)  # No reference RT
            return fallback_row
        end
    end
end

#==========================================================
Global Cross-Run Prescore Aggregation
==========================================================#

const GLOBAL_PRESCORE_QVALUE_THRESHOLD = 0.15f0

"""
Log-odds average of top-N probabilities, converted back to probability space.
Same approach as FirstPassSearch's _logodds_combine.
"""
function _prescore_logodds_combine(probs::Vector{Float32}, top_n::Int)::Float32
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    sorted = sort(probs; rev=true)
    selected = @view sorted[1:n]
    eps = 1f-6
    lo = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(lo) / n
    return 1.0f0 / (1 + exp(-avg))
end

"""
    aggregate_prescore_globally!(search_context::SearchContext) -> Set{UInt32}

Load per-file LightGBM prescore arrow files, aggregate via log-odds averaging,
compute global q-values, and return the set of precursor indices passing at
q ≤ GLOBAL_PRESCORE_QVALUE_THRESHOLD.
"""
function aggregate_prescore_globally!(search_context::SearchContext)
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")

    # Collect per-file best probs per precursor
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    prec_is_target = Dictionary{UInt32, Bool}()
    n_valid_files = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        file_name = getParsedFileName(ms_data, ms_file_idx)
        score_path = joinpath(prescore_dir, "$(file_name).arrow")
        !isfile(score_path) && continue

        scores = Arrow.Table(score_path)
        n_valid_files += 1

        for i in eachindex(scores[:precursor_idx])
            pid = scores[:precursor_idx][i]
            p = scores[:lgbm_prob][i]
            if haskey(prec_probs_by_run, pid)
                push!(prec_probs_by_run[pid], p)
            else
                insert!(prec_probs_by_run, pid, Float32[p])
            end
            if !haskey(prec_is_target, pid)
                insert!(prec_is_target, pid, scores[:target][i])
            end
        end
    end

    # Aggregate via log-odds
    sqrt_n = max(1, floor(Int, sqrt(n_valid_files)))
    n_unique = length(prec_probs_by_run)
    global_prec_idxs = Vector{UInt32}(undef, n_unique)
    global_probs = Vector{Float32}(undef, n_unique)
    global_targets = Vector{Bool}(undef, n_unique)

    for (i, (pid, probs)) in enumerate(pairs(prec_probs_by_run))
        global_prec_idxs[i] = pid
        global_probs[i] = _prescore_logodds_combine(probs, sqrt_n)
        global_targets[i] = prec_is_target[pid]
    end

    # Q-values on global scores
    global_qvals = Vector{Float32}(undef, n_unique)
    get_qvalues!(global_probs, global_targets, global_qvals; doSort=true)

    # Diagnostics
    n_targets = count(global_targets)
    n_1pct = count(i -> global_qvals[i] <= 0.01f0 && global_targets[i], eachindex(global_qvals))
    n_5pct = count(i -> global_qvals[i] <= 0.05f0 && global_targets[i], eachindex(global_qvals))
    n_10pct = count(i -> global_qvals[i] <= 0.10f0 && global_targets[i], eachindex(global_qvals))
    n_15pct = count(i -> global_qvals[i] <= 0.15f0 && global_targets[i], eachindex(global_qvals))
    @info "Global prescore: $n_unique precursors ($n_targets targets) from $n_valid_files files (top_n=$sqrt_n)"
    @info "  Targets at q≤0.01: $n_1pct, q≤0.05: $n_5pct, q≤0.10: $n_10pct, q≤0.15: $n_15pct"

    # Build passing set
    passing = Set{UInt32}()
    for i in eachindex(global_prec_idxs)
        if global_qvals[i] <= GLOBAL_PRESCORE_QVALUE_THRESHOLD
            push!(passing, global_prec_idxs[i])
        end
    end

    @info "Global prescore filter: $(length(passing)) precursors pass at q≤$(GLOBAL_PRESCORE_QVALUE_THRESHOLD)"

    # ── DIA-NN diagnostics ──────────────────────────────────────
    # Build lookup: precursor_idx → global_qval for all scored precursors
    prec_to_global_qval = Dictionary{UInt32, Float32}()
    for i in eachindex(global_prec_idxs)
        insert!(prec_to_global_qval, global_prec_idxs[i], global_qvals[i])
    end

    # Global DIA-NN reference
    if isfile(DIANN_GLOBAL_ARROW)
        diann_global = Set{UInt32}(Arrow.Table(DIANN_GLOBAL_ARROW).precursor_idx)
        n_diann = length(diann_global)
        n_pass = count(pid -> pid in passing, diann_global)
        n_scored_fail = count(pid -> haskey(prec_to_global_qval, pid) && !(pid in passing), diann_global)
        n_not_scored = n_diann - n_pass - n_scored_fail
        @info "DIA-NN global ($n_diann): pass=$n_pass, scored-but-filtered=$n_scored_fail, not-scored=$n_not_scored"
    end

    # Per-file DIA-NN references
    if isdir(DIANN_PERFILE_DIR)
        for ms_file_idx in 1:n_files
            getFailedIndicator(ms_data, ms_file_idx) && continue
            file_name = getParsedFileName(ms_data, ms_file_idx)
            pf_path = joinpath(DIANN_PERFILE_DIR, "$(file_name).arrow")
            !isfile(pf_path) && continue

            diann_file = Set{UInt32}(Arrow.Table(pf_path).precursor_idx)
            n_diann_f = length(diann_file)
            n_pass_f = count(pid -> pid in passing, diann_file)
            n_scored_fail_f = count(pid -> haskey(prec_to_global_qval, pid) && !(pid in passing), diann_file)
            n_not_scored_f = n_diann_f - n_pass_f - n_scored_fail_f
            @info "  DIA-NN file $file_name ($n_diann_f): pass=$n_pass_f, scored-but-filtered=$n_scored_fail_f, not-scored=$n_not_scored_f"
        end
    end

    return passing
end

"""
    filter_arrow_files_to_passing!(search_context::SearchContext, passing_precs::Set{UInt32})

Read each fold-split second pass arrow file, filter rows to keep only precursors
in the passing set, and rewrite in-place.
"""
function filter_arrow_files_to_passing!(search_context::SearchContext, passing_precs::Set{UInt32})
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    n_removed_total = 0
    n_kept_total = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        base_path = getSecondPassPsms(ms_data, ms_file_idx)
        isempty(base_path) && continue

        for fold in UInt8[0, 1]
            fold_path = getSecondPassPsmsFold(ms_data, ms_file_idx, fold)
            !isfile(fold_path) && continue

            tbl = DataFrame(Tables.columntable(Arrow.Table(fold_path)))
            n_before = nrow(tbl)
            filter!(row -> row.precursor_idx in passing_precs, tbl)
            n_after = nrow(tbl)
            n_removed_total += (n_before - n_after)
            n_kept_total += n_after

            writeArrow(fold_path, tbl)  # Overwrite in-place
        end
    end

    @info "Global filter applied: kept $n_kept_total PSMs, removed $n_removed_total across all files"
end

"""
    rerun_globally_filtered!(search_context, params, passing_precs)

Re-run deconvolution for each file with only globally-passing precursors,
then recompute all PSM features and overwrite fold-split arrow files.
This reduces precursor competition in the design matrix, producing better
weights and features for the final LightGBM model in ScoringSearch.
"""
function rerun_globally_filtered!(
    search_context::SearchContext,
    params::P,
    passing_precs::Set{UInt32}
) where {P<:SecondPassSearchParameters}
    msdr = getMassSpecData(search_context)
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    n_rerun = 0
    n_psms_total = 0

    @info "=== Global deconvolution re-run: $(length(passing_precs)) passing precursors ==="

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        base_path = getSecondPassPsms(ms_data, ms_file_idx)
        isempty(base_path) && continue

        file_name = getParsedFileName(search_context, ms_file_idx)
        @info "  Re-running file $ms_file_idx ($file_name)..."

        try
            # Load spectra for this file
            spectra = getMSData(msdr, ms_file_idx)

            # Re-run deconvolution with only globally-passing precursors
            psms = rerun_search_with_precursor_filter(
                spectra, search_context, params, ms_file_idx, passing_precs
            )

            if nrow(psms) == 0
                @info "    No PSMs after global re-search for file $ms_file_idx, keeping existing arrow files"
                continue
            end

            # Recompute all features on fresh PSMs
            ms1_psms = DataFrame()  # MS1 scoring disabled in bypass mode
            psms = compute_psm_features!(psms, ms1_psms, params, search_context, ms_file_idx, spectra)

            if nrow(psms) == 0
                @info "    No PSMs after feature computation for file $ms_file_idx, keeping existing arrow files"
                continue
            end

            # Overwrite fold-split arrow files
            base_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
            for fold in UInt8[0, 1]
                fold_path = joinpath(base_dir, "$(file_name)_fold$(fold).arrow")
                fold_mask = psms.cv_fold .== fold
                if any(fold_mask)
                    writeArrow(fold_path, psms[fold_mask, :])
                elseif isfile(fold_path)
                    # Remove stale fold file if no PSMs for this fold
                    rm(fold_path)
                end
            end

            n_rerun += 1
            n_psms_total += nrow(psms)
            @info "    File $ms_file_idx: $(nrow(psms)) PSMs after global re-deconvolution"

        catch e
            @user_warn "Global re-deconvolution failed for file $ms_file_idx ($file_name), keeping existing arrow files. Error: $(sprint(showerror, e))"
        end
    end

    @info "=== Global re-run complete: $n_rerun files, $n_psms_total total PSMs ==="
end
