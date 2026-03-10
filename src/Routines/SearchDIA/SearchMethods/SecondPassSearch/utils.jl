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
    ::MS2CHROM;
    n_frag_isotopes::Int64 = params.n_frag_isotopes,
    max_frag_rank::UInt8 = params.max_frag_rank,
    min_frag_count::Int64 = params.min_frag_count
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
                ms_file_idx;
                n_frag_isotopes = n_frag_isotopes,
                max_frag_rank = max_frag_rank,
                min_frag_count = min_frag_count
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
    ms_file_idx::Int64;
    n_frag_isotopes::Int64 = params.n_frag_isotopes,
    max_frag_rank::UInt8 = params.max_frag_rank,
    min_frag_count::Int64 = params.min_frag_count
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
    total_skipped_matched_ratio = 0

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
            n_frag_isotopes,
            max_frag_rank,
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
            min_frag_count = min_frag_count,
            max_best_rank = params.max_best_rank,
            min_topn = first(params.min_topn_of_m),
            block_size = 500000
        )
        last_val = score_result.last_val
        total_skipped_weight += score_result.skipped_weight
        total_skipped_frag_count += score_result.skipped_frag_count
        total_skipped_matched_ratio += score_result.skipped_matched_ratio

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    @info "SecondPass (fragindex) filter summary: kept=$last_val, skipped_weight=$total_skipped_weight, skipped_frag_count=$total_skipped_frag_count, skipped_matched_ratio=$total_skipped_matched_ratio"
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
    total_skipped_matched_ratio = 0

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
        total_skipped_matched_ratio += score_result.skipped_matched_ratio

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    @info "SecondPass (RT-indexed) filter summary: kept=$last_val, skipped_weight=$total_skipped_weight, skipped_frag_count=$total_skipped_frag_count, skipped_matched_ratio=$total_skipped_matched_ratio"
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
    n_before_weight = nrow(psms)
    filter!(x->x.weight>0.0, psms);
    n_removed_weight = n_before_weight - nrow(psms)
    @info "  Weight > 0 filter: removed $n_removed_weight/$n_before_weight PSMs"
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
                                    prec_id_to_irt::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}};
                                    prescore_only::Bool=false
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
    # weight>0 already enforced upstream in add_second_search_columns!
    #filter!(x->x.data_points>0, psms)
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    irt_obs = zeros(Float32, N)
    irt_pred = zeros(Float32, N)
    irt_error = zeros(Float32, N)
    missed_cleavage = zeros(UInt8, N);
    spectrum_peak_count = zeros(Float16, N);
    Mox = zeros(UInt8, N);

    # Columns only needed for Phase 2 (full feature set)
    if !prescore_only
        irt_diff = zeros(Float32, N)
        ms1_irt_diff = zeros(Float32, N)
        pair_idxs = zeros(UInt32, N)
        entrap_group_id = zeros(UInt8, N)
        adjusted_intensity_explained = zeros(Float16, N);
        prec_charges = zeros(UInt8, N)
        sequence_length = zeros(UInt8, N);
        prec_mzs = zeros(Float32, N);
        TIC = zeros(Float16, N);

        # Amino acid count features (20 standard amino acids)
        aa_counts = Dict{Char, Vector{UInt8}}()
        for aa in "ACDEFGHIKLMNPQRSTVWY"
            aa_counts[aa] = zeros(UInt8, N)
        end
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

    # Pre-compute per-precursor string features to avoid redundant work across ~3 PSMs/precursor
    AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
    unique_precs = unique(precursor_idx)
    _mox_cache = Dict{UInt32, UInt8}()
    sizehint!(_mox_cache, length(unique_precs))
    if !prescore_only
        _seq_length_cache = Dict{UInt32, UInt8}()
        _aa_cache = Dict{UInt32, NTuple{20, UInt8}}()
        sizehint!(_seq_length_cache, length(unique_precs))
        sizehint!(_aa_cache, length(unique_precs))
        for pid in unique_precs
            stripped = replace(precursor_sequence[pid], r"\(.*?\)" => "")
            _seq_length_cache[pid] = UInt8(length(stripped))
            _mox_cache[pid] = countMOX(structural_mods[pid])
            _aa_cache[pid] = ntuple(Val(20)) do j
                UInt8(count(==(AA_ORDER[j]), stripped))
            end
        end
    else
        for pid in unique_precs
            _mox_cache[pid] = countMOX(structural_mods[pid])
        end
    end

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin 
            for i in chunk
                prec_idx = precursor_idx[i]
                irt_obs[i] = rt_to_irt_interp(rt[i])
                irt_pred[i] = getPredIrt(search_context, prec_idx)
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])
                missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
                Mox[i] = _mox_cache[prec_idx]
                spectrum_peak_count[i] = length(masses[scan_idx[i]])

                if !prescore_only
                    entrap_group_id[i] = entrap_group_ids[prec_idx]
                    irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_irt)
                    if !ms1_missing[i]
                        ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx))
                    else
                        ms1_irt_diff[i] = 0f0
                    end
                    sequence_length[i] = _seq_length_cache[prec_idx]
                    aa_tuple = _aa_cache[prec_idx]
                    @inbounds for (j, aa) in enumerate(AA_ORDER)
                        aa_counts[aa][i] = aa_tuple[j]
                    end
                    TIC[i] = Float16(log2(tic[scan_idx[i]]))
                    adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                    prec_charges[i] = prec_charge[prec_idx]
                    pair_idxs[i] = extract_pair_idx(precursor_pair_idxs, prec_idx)
                    prec_mzs[i] = prec_mz[prec_idx];
                end
            end
        end
    end
    fetch.(tasks)

    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_error] = irt_error
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:spectrum_peak_count] = spectrum_peak_count

    if !prescore_only
        psms[!,:irt_diff] = irt_diff
        psms[!,:ms1_irt_diff] = ms1_irt_diff
        psms[!,:sequence_length] = sequence_length
        psms[!,:tic] = TIC
        psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
        psms[!,:charge] = prec_charges
        psms[!,:pair_id] = pair_idxs
        psms[!,:prec_mz] = prec_mzs
        psms[!,:entrapment_group_id] = entrap_group_id
        psms[!,:ms_file_idx] .= ms_file_idx

        # Amino acid count columns
        for aa in "ACDEFGHIKLMNPQRSTVWY"
            psms[!, Symbol("aa_", aa)] = aa_counts[aa]
        end
    end
    return nothing
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
        @debug "  DIA-NN per-file ($file_name): $(length(per_file)) unique precursor_idx values"
    end

    # Global reference
    if isfile(DIANN_GLOBAL_ARROW)
        t = Arrow.Table(DIANN_GLOBAL_ARROW)
        global_set = Set{UInt32}(t.precursor_idx)
        @debug "  DIA-NN global: $(length(global_set)) unique precursor_idx values"
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
        @debug "  [$label] DIA-NN file recovery: $n_recovered / $n_total ($pct%)"
    end
    if !isempty(diann_global)
        n_recovered = length(intersect(surviving_precs, diann_global))
        n_total = length(diann_global)
        pct = round(100.0 * n_recovered / max(1, n_total), digits=1)
        @debug "  [$label] DIA-NN global recovery: $n_recovered / $n_total ($pct%)"
    end
end

#==========================================================
LightGBM Feature Set
==========================================================#

# Lean feature set for Phase 1 prescore LightGBM (fast per-file ranking)
const PRESCORE_FEATURES = [
    :fitted_manhattan_distance, :scribe, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :max_unmatched_residual, :max_matched_residual, :Mox, :spectrum_peak_count,
]

# Full feature set used in Phase 2 (ScoringSearch gets these via fold Arrow files)
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
    _sanitize_column!(col::Vector{<:AbstractFloat}) -> Int

Replace Inf/NaN values with zero in a float column. Returns count of values fixed.
Typed inner function so Julia compiles a tight loop per concrete element type.
"""
function _sanitize_column!(col::Vector{T}) where {T<:AbstractFloat}
    n = 0
    @inbounds for i in eachindex(col)
        if !isfinite(col[i])
            col[i] = zero(T)
            n += 1
        end
    end
    return n
end

# No-op for integer/bool columns (always finite)
_sanitize_column!(::AbstractVector{<:Union{Integer, Bool}}) = 0

"""
    sanitize_prescore_features!(psms::DataFrame, features::Vector{Symbol})

Replace Inf/NaN values with 0.0 in feature columns.
Float16 columns (e.g., err_norm) can overflow to Inf, which would
crash downstream ML models. Integer/Bool columns are skipped (always finite).
"""
function sanitize_prescore_features!(psms::DataFrame, features::Vector{Symbol})
    n_fixed = 0
    for f in features
        hasproperty(psms, f) || continue
        n_fixed += _sanitize_column!(psms[!, f])
    end
    if n_fixed > 0
        @info "  Sanitized $n_fixed non-finite feature values (Inf/NaN → 0)\n"
    end
end

"""
    prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra) -> DataFrame

Reusable feature computation pipeline for second pass PSMs. Called in both
Phase 1 (per-file prescore) and Phase 2 (after global re-deconvolution).

Steps:
1. add_second_search_columns! — RT, charge, target, cv_fold, err_norm, total_ions
2. get_isotopes_captured! — precursor_fraction_transmitted
3. Filter by fraction_transmitted (weight>0 enforced in step 1)
4. ensure_ms1_stub_columns! — stubs needed by add_features!
5. add_features! — irt_error, irt_diff, tic, prec_mz, sequence_length, etc.
6. sanitize_prescore_features! — replace Inf/NaN with 0
"""
function prepare_psm_features!(
    psms::DataFrame,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData;
    prescore_only::Bool=false
) where {P<:SecondPassSearchParameters}
    t0 = time()

    # 1. Add basic search columns (RT, charge, target/decoy status, cv_fold)
    add_second_search_columns!(psms,
        getRetentionTimes(spectra),
        getCharge(getPrecursors(getSpecLib(search_context))),
        getIsDecoy(getPrecursors(getSpecLib(search_context))),
        getPrecursors(getSpecLib(search_context))
    )
    t1 = time()

    # 2. Determine which precursor isotopes are captured in each scan's isolation window
    get_isotopes_captured!(
        psms,
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
    t2 = time()

    # 3. Filter by fraction_transmitted (weight>0 already enforced in add_second_search_columns!)
    to_remove = findall(psms[!, :precursor_fraction_transmitted] .< params.min_fraction_transmitted)
    deleteat!(psms, to_remove)
    t3 = time()

    # 4. Add MS1 stub columns (needed by add_features!)
    ensure_ms1_stub_columns!(psms)
    t4 = time()

    # 5. Add ML features (irt_error, irt_diff, tic, prec_mz, sequence_length, etc.)
    add_features!(
        psms,
        search_context,
        getTICs(spectra),
        getMzArrays(spectra),
        ms_file_idx,
        getRtIrtModel(search_context, ms_file_idx),
        getPrecursorDict(search_context);
        prescore_only=prescore_only
    )
    t5 = time()

    # 6. Sanitize features (Inf/NaN → 0)
    features_to_sanitize = prescore_only ? PRESCORE_FEATURES : LGBM_RECOVERY_FEATURES
    sanitize_prescore_features!(psms, collect(features_to_sanitize))
    t6 = time()

    r = s -> round(s, digits=3)
    @info "  prepare_psm_features! ($(nrow(psms)) PSMs, prescore_only=$prescore_only): " *
          "columns=$(r(t1-t0))s, isotopes=$(r(t2-t1))s, filter=$(r(t3-t2))s, " *
          "ms1_stubs=$(r(t4-t3))s, add_features=$(r(t5-t4))s, sanitize=$(r(t6-t5))s, " *
          "total=$(r(t6-t0))s"

    return psms
end

"""
    train_lgbm_and_select_best(psms; features) -> (best_psms, scores, q_values, timings)

Train LightGBM on ALL PSMs (all scans), predict scores, select best scan
per precursor by LightGBM score, compute q-values, and log diagnostics.

Returns:
- best_psms: DataFrame with one row per precursor (best by LightGBM score)
- scores: Vector{Float32} of LightGBM probabilities for best_psms
- q_values: Vector{Float64} of q-values for best_psms
"""
function train_lgbm_and_select_best(
    psms::DataFrame;
    features::Vector{Symbol} = collect(PRESCORE_FEATURES)
)
    t0 = time()

    # Filter out PSMs with fewer than 3 matched fragments (b + y ions)
    n_before = nrow(psms)
    precs_before = length(unique(psms[!, :precursor_idx]))
    frag_mask = (psms[!, :y_count] .+ psms[!, :b_count]) .>= UInt8(4)
    psms = psms[frag_mask, :]
    n_removed = n_before - nrow(psms)
    precs_removed = precs_before - length(unique(psms[!, :precursor_idx]))
    @info "  Min fragment filter (≥4 b+y): removed $n_removed/$n_before PSMs, $precs_removed/$precs_before unique precursors"

    # Filter to available features
    available_features = filter(f -> hasproperty(psms, f), features)
    targets_col = psms[!, :target]

    # Build feature matrix for LightGBM
    feature_df = psms[!, available_features]
    t_matrix = time()

    # Train LightGBM classifier
    classifier = build_lightgbm_classifier(
        num_iterations = 50,
        learning_rate = 0.2,
        max_depth = 5,
        num_leaves = 15,
        min_data_in_leaf = 200,
        feature_fraction = 0.5,
        bagging_fraction = 0.5,
        bagging_freq = 1,
        is_unbalance = true,
    )
    # Subsample for training if > 1M PSMs
    max_train = 1_000_000
    n_total = nrow(feature_df)
    if n_total > max_train
        train_idx = randperm(n_total)[1:max_train]
        @info "  LightGBM training: subsampled $max_train / $n_total PSMs"
        model = fit_lightgbm_model(classifier, feature_df[train_idx, :], targets_col[train_idx])
    else
        model = fit_lightgbm_model(classifier, feature_df, targets_col)
    end
    t_train = time()

    # Predict on ALL PSMs
    all_scores = lightgbm_predict(model, feature_df)
    t_predict = time()

    # Add scores to psms for best-per-precursor selection
    psms[!, :lgbm_score] = Float32.(all_scores)

    # Select best scan per precursor by LightGBM score
    best_psms = get_best_psm_per_precursor_by_score(psms, :lgbm_score)

    # Extract scores for best PSMs
    scores = best_psms[!, :lgbm_score]

    # Compute q-values
    best_targets = best_psms[!, :target]
    q_values = zeros(Float64, nrow(best_psms))
    get_qvalues!(scores, best_targets, q_values)
    t_select = time()

    # Feature importances (debug only)
    imp = importance(model)
    if imp !== nothing
        sorted_imp = sort(imp, by = x -> -x[2])
        @debug "  LightGBM feature importances ($(length(sorted_imp)) features):"
        for (fname, gain) in sorted_imp
            @debug "    $fname: $(round(gain, digits=1))"
        end
    end

    # Clean up temporary column
    select!(psms, Not(:lgbm_score))

    timings = (
        matrix = t_matrix - t0,
        train = t_train - t_matrix,
        predict = t_predict - t_train,
        select = t_select - t_predict,
    )

    return best_psms, Vector{Float32}(scores), q_values, timings
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
    scores = psms[!, score_col]
    prec_ids = psms[!, :precursor_idx]

    # Single pass: track best row index per precursor
    best_idx = Dictionary{UInt32, Int}()
    best_score = Dictionary{UInt32, Float32}()
    for i in eachindex(prec_ids)
        pid = prec_ids[i]
        s = scores[i]
        if !haskey(best_idx, pid) || s > best_score[pid]
            if haskey(best_idx, pid)
                best_idx[pid] = i
                best_score[pid] = s
            else
                insert!(best_idx, pid, i)
                insert!(best_score, pid, s)
            end
        end
    end

    # Collect best rows — sort indices for cache-friendly DataFrame access
    row_indices = collect(values(best_idx))
    sort!(row_indices)
    return psms[row_indices, :]
end

"""
    get_best_psm_per_precursor_by_prescore_scan(psms::DataFrame, best_scans::Dictionary{UInt32, UInt32}) -> DataFrame

Select one PSM per precursor using the Phase 1 LightGBM-selected scan_idx.
Falls back to highest deconvolution weight if the preferred scan isn't present.
`best_scans` maps precursor_idx → preferred scan_idx for this file.
"""
function get_best_psm_per_precursor_by_prescore_scan(
    psms::DataFrame,
    best_scans::Dictionary{UInt32, UInt32}
)
    # Sort by weight descending as fallback ordering
    sorted = sort(psms, [:precursor_idx, order(:weight, rev=true)])

    best = combine(groupby(sorted, :precursor_idx)) do group
        pid = first(group.precursor_idx)
        if haskey(best_scans, pid)
            target_scan = best_scans[pid]
            match_idx = findfirst(==(target_scan), group.scan_idx)
            if match_idx !== nothing
                return group[match_idx:match_idx, :]
            end
        end
        # Fallback: highest weight (first row, since sorted desc)
        return group[1:1, :]
    end

    return best
end

"""
    add_multi_scan_aggregates!(best_psms, all_psms)

Compute per-precursor aggregate features from all scans and join them onto
the best-scan DataFrame. Gives ScoringSearch multi-scan signal.
"""
function add_multi_scan_aggregates!(best_psms::DataFrame, all_psms::DataFrame)
    aggs = combine(groupby(all_psms, :precursor_idx),
        nrow => :num_scans,
        :y_count => maximum => :max_y_ions,
        :y_count => sum => :y_ions_sum,
        :gof => maximum => :max_gof,
        :fitted_manhattan_distance => maximum => :max_fitted_manhattan_distance,
        :fitted_spectral_contrast => maximum => :max_fitted_spectral_contrast,
        :matched_ratio => maximum => :max_matched_ratio,
        :scribe => maximum => :max_scribe,
        :weight => maximum => :max_weight,
    )
    aggs[!, :num_scans] = UInt16.(aggs[!, :num_scans])
    aggs[!, :y_ions_sum] = UInt16.(aggs[!, :y_ions_sum])
    leftjoin!(best_psms, aggs, on = :precursor_idx)
    return best_psms
end

#==========================================================
Global Cross-Run Prescore Aggregation
==========================================================#

const DEFAULT_GLOBAL_PRESCORE_QVALUE_THRESHOLD = 0.05f0

"""
    aggregate_prescore_globally!(search_context::SearchContext, qvalue_threshold::Float32, aggregation::PrescoreAggregationStrategy) -> (Set{UInt32}, Dictionary{UInt32, Dictionary{Int, UInt32}})

Load per-file LightGBM prescore arrow files, calibrate per-file scores via
`calibrate_file_scores(aggregation, ...)`, aggregate via `combine_scores(aggregation, ...)`,
compute global q-values, and return (1) the set of precursor indices passing at
q ≤ qvalue_threshold and (2) per-precursor Phase 1 best scan lookup mapping
precursor_idx → (file_idx → scan_idx).
"""
function aggregate_prescore_globally!(search_context::SearchContext,
                                      qvalue_threshold::Float32=DEFAULT_GLOBAL_PRESCORE_QVALUE_THRESHOLD,
                                      aggregation::PrescoreAggregationStrategy=PEPCalibratedAggregation())
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")

    # Collect per-file best probs per precursor + iRT tracking
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    prec_is_target = Dictionary{UInt32, Bool}()
    prec_best_prob = Dictionary{UInt32, Float32}()    # best lgbm_prob seen so far
    prec_best_irt = Dictionary{UInt32, Float32}()     # iRT of globally best-scoring scan
    prec_min_irt = Dictionary{UInt32, Float32}()       # min iRT across files
    prec_max_irt = Dictionary{UInt32, Float32}()       # max iRT across files
    prec_peak_widths = Dictionary{UInt32, Vector{Float32}}()  # within-file iRT peak widths
    prec_irt_by_file = Dictionary{UInt32, Vector{Pair{Int,Float32}}}()  # (file_idx => irt) per observation
    prec_best_scan = Dictionary{UInt32, Dictionary{Int, UInt32}}()  # pid → (file_idx → scan_idx)
    n_valid_files = 0

    @info "Prescore aggregation strategy: $(typeof(aggregation))"

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        file_name = getParsedFileName(ms_data, ms_file_idx)
        score_path = joinpath(prescore_dir, "$(file_name).arrow")
        !isfile(score_path) && continue

        scores = Arrow.Table(score_path)
        n_valid_files += 1

        has_irt = :irt_obs in Tables.columnnames(scores)
        has_peak_width = :irt_peak_width in Tables.columnnames(scores)
        has_scan = :scan_idx in Tables.columnnames(scores)

        # Calibrate per-file scores using the aggregation strategy
        file_probs = Float32.(scores[:lgbm_prob])
        file_targets = Bool.(scores[:target])
        calibrated_probs = calibrate_file_scores(aggregation, file_probs, file_targets)

        for i in eachindex(scores[:precursor_idx])
            pid = scores[:precursor_idx][i]
            p = calibrated_probs[i]
            if haskey(prec_probs_by_run, pid)
                push!(prec_probs_by_run[pid], p)
            else
                insert!(prec_probs_by_run, pid, Float32[p])
            end
            if !haskey(prec_is_target, pid)
                insert!(prec_is_target, pid, file_targets[i])
            end

            # Track iRT range per precursor
            if has_irt
                irt = Float32(scores[:irt_obs][i])
                if haskey(prec_min_irt, pid)
                    prec_min_irt[pid] = min(prec_min_irt[pid], irt)
                    prec_max_irt[pid] = max(prec_max_irt[pid], irt)
                    if p > prec_best_prob[pid]
                        prec_best_prob[pid] = p
                        prec_best_irt[pid] = irt
                    end
                else
                    insert!(prec_min_irt, pid, irt)
                    insert!(prec_max_irt, pid, irt)
                    insert!(prec_best_prob, pid, p)
                    insert!(prec_best_irt, pid, irt)
                end
                if haskey(prec_irt_by_file, pid)
                    push!(prec_irt_by_file[pid], ms_file_idx => irt)
                else
                    insert!(prec_irt_by_file, pid, [ms_file_idx => irt])
                end
            end

            # Track within-file iRT peak width per precursor
            if has_peak_width
                pw = Float32(scores[:irt_peak_width][i])
                if haskey(prec_peak_widths, pid)
                    push!(prec_peak_widths[pid], pw)
                else
                    insert!(prec_peak_widths, pid, Float32[pw])
                end
            end

            # Track Phase 1 best scan per precursor per file
            if has_scan
                sid = UInt32(scores[:scan_idx][i])
                if !haskey(prec_best_scan, pid)
                    insert!(prec_best_scan, pid, Dictionary{Int, UInt32}())
                end
                insert!(prec_best_scan[pid], ms_file_idx, sid)
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
        global_probs[i] = combine_scores(aggregation, probs, sqrt_n)
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
        if global_qvals[i] <= qvalue_threshold
            push!(passing, global_prec_idxs[i])
        end
    end

    @info "Global prescore filter: $(length(passing)) precursors pass at q≤$(qvalue_threshold)"

    #= ── iRT spread diagnostics for 1% FDR targets ──────────────
    if !isempty(prec_min_irt)
        # Collect iRT diffs for target precursors at 1% global FDR
        irt_diffs = Float32[]
        for i in eachindex(global_prec_idxs)
            pid = global_prec_idxs[i]
            global_qvals[i] <= 0.01f0 && global_targets[i] && haskey(prec_min_irt, pid) || continue
            push!(irt_diffs, prec_max_irt[pid] - prec_min_irt[pid])
        end

        if !isempty(irt_diffs)
            sorted_diffs = sort(irt_diffs)
            n_d = length(sorted_diffs)
            med = sorted_diffs[max(1, div(n_d + 1, 2))]
            avg = sum(sorted_diffs) / n_d
            sd = sqrt(sum((d - avg)^2 for d in sorted_diffs) / max(1, n_d - 1))
            p25 = sorted_diffs[max(1, ceil(Int, 0.25 * n_d))]
            p75 = sorted_diffs[max(1, ceil(Int, 0.75 * n_d))]
            p90 = sorted_diffs[max(1, ceil(Int, 0.90 * n_d))]
            p95 = sorted_diffs[max(1, ceil(Int, 0.95 * n_d))]
            @info "iRT spread (max-min) for $(n_d) target precursors at q≤0.01:"
            @info "  median=$(round(med; digits=2)), mean=$(round(avg; digits=2)), std=$(round(sd; digits=2))"
            @info "  p25=$(round(p25; digits=2)), p75=$(round(p75; digits=2)), p90=$(round(p90; digits=2)), p95=$(round(p95; digits=2))"

            # Save histogram
            try
                results_dir = getDataOutDir(search_context)
                p_hist = histogram(irt_diffs;
                    xlabel = "iRT spread (max - min across files)",
                    ylabel = "Count",
                    title = "Prescore iRT spread ($(n_d) targets, q≤1%)",
                    legend = false,
                    bins = min(100, max(20, div(n_d, 50)))
                )
                savefig(p_hist, joinpath(results_dir, "prescore_irt_spread.pdf"))
                @info "Saved prescore_irt_spread.pdf to $results_dir"
            catch e
                @warn "Failed to save iRT spread histogram: $e"
            end

            # ── iRT outlier detection (per-file observations) ──
            med_spread = sorted_diffs[max(1, div(length(sorted_diffs) + 1, 2))]
            abs_devs = sort([abs(d - med_spread) for d in sorted_diffs])
            mad_raw = abs_devs[max(1, div(length(abs_devs) + 1, 2))]
            mad_normalized = mad_raw * 1.4826f0
            irt_outlier_threshold = 3.0f0 * mad_normalized

            n_outlier = 0
            n_total_obs = 0
            for i in eachindex(global_prec_idxs)
                pid = global_prec_idxs[i]
                global_qvals[i] <= 0.01f0 && global_targets[i] && haskey(prec_irt_by_file, pid) || continue
                best_irt = prec_best_irt[pid]
                for (file_idx, file_irt) in prec_irt_by_file[pid]
                    n_total_obs += 1
                    if abs(file_irt - best_irt) > irt_outlier_threshold
                        n_outlier += 1
                    end
                end
            end

            pct = round(100.0 * n_outlier / max(1, n_total_obs); digits=2)
            @info "iRT outlier detection (q≤0.01 targets): MAD=$(round(mad_raw; digits=4)), threshold=$(round(irt_outlier_threshold; digits=4))"
            @info "  $n_outlier / $n_total_obs observations ($pct%) exceed 3σ from best-scoring iRT (diagnostic only, not filtering)"
        end
    end

    # ── Within-file iRT peak width diagnostics for 1% FDR targets ──
    if !isempty(prec_peak_widths)
        # Collect all per-file peak widths for target precursors at 1% global FDR
        pw_all = Float32[]
        for i in eachindex(global_prec_idxs)
            pid = global_prec_idxs[i]
            global_qvals[i] <= 0.01f0 && global_targets[i] && haskey(prec_peak_widths, pid) || continue
            append!(pw_all, prec_peak_widths[pid])
        end

        if !isempty(pw_all)
            sorted_pw = sort(pw_all)
            n_pw = length(sorted_pw)
            med = sorted_pw[max(1, div(n_pw + 1, 2))]
            avg = sum(sorted_pw) / n_pw
            sd = sqrt(sum((w - avg)^2 for w in sorted_pw) / max(1, n_pw - 1))
            p25 = sorted_pw[max(1, ceil(Int, 0.25 * n_pw))]
            p75 = sorted_pw[max(1, ceil(Int, 0.75 * n_pw))]
            p90 = sorted_pw[max(1, ceil(Int, 0.90 * n_pw))]
            p95 = sorted_pw[max(1, ceil(Int, 0.95 * n_pw))]
            @info "iRT peak width (within-file) for $(n_pw) observations from target precursors at q≤0.01:"
            @info "  median=$(round(med; digits=2)), mean=$(round(avg; digits=2)), std=$(round(sd; digits=2))"
            @info "  p25=$(round(p25; digits=2)), p75=$(round(p75; digits=2)), p90=$(round(p90; digits=2)), p95=$(round(p95; digits=2))"

            # Save histogram
            try
                results_dir = getDataOutDir(search_context)
                p_hist = histogram(pw_all;
                    xlabel = "iRT peak width (within-file)",
                    ylabel = "Count",
                    title = "Prescore iRT peak width ($(n_pw) obs, q≤1%)",
                    legend = false,
                    bins = min(100, max(20, div(n_pw, 50)))
                )
                savefig(p_hist, joinpath(results_dir, "prescore_irt_peak_width.pdf"))
                @info "Saved prescore_irt_peak_width.pdf to $results_dir"
            catch e
                @warn "Failed to save iRT peak width histogram: $e"
            end
        end
    end
    =#

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
        @debug "DIA-NN global ($n_diann): pass=$n_pass, scored-but-filtered=$n_scored_fail, not-scored=$n_not_scored"
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
            @debug "  DIA-NN file $file_name ($n_diann_f): pass=$n_pass_f, scored-but-filtered=$n_scored_fail_f, not-scored=$n_not_scored_f"
        end
    end

    return passing, prec_best_scan
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

