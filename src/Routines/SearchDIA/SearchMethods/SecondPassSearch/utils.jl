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
    thread_results = Vector{@NamedTuple{psms::DataFrame}}(undef, length(tasks))
    for (i, t) in enumerate(tasks)
        try
            thread_results[i] = fetch(t)
        catch e
            bt = catch_backtrace()
            @user_error "SecondPassSearch task $(i) failed while fetching results (MS2CHROM)"
            @user_error sprint(showerror, e, bt)
            rethrow(e)
        end
    end
    return (psms = vcat([r.psms for r in thread_results]...),)
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
    min_frag_count::Int64 = params.min_frag_count,
    min_spectral_contrast::Float32 = params.min_spectral_contrast,
    min_log2_matched_ratio::Float32 = params.min_log2_matched_ratio,
    min_topn_of_m::Tuple{Int64, Int64} = params.min_topn_of_m,
    dynamic_range::Float32 = params.dynamic_range,
    first_pass::Bool = false
)
    # Use single thread for first_pass to enable clean profiling
    n_threads = first_pass ? 1 : Threads.nthreads()
    thread_tasks = partition_scans(spectra, n_threads)

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
                min_frag_count = min_frag_count,
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_topn_of_m = min_topn_of_m,
                dynamic_range = dynamic_range,
                first_pass = first_pass
            )
        end
    end

    thread_results = Vector{@NamedTuple{psms::DataFrame}}(undef, length(tasks))
    for (i, t) in enumerate(tasks)
        try
            thread_results[i] = fetch(t)
        catch e
            bt = catch_backtrace()
            @user_error "SecondPassSearch (fragindex) task $(i) failed"
            @user_error sprint(showerror, e, bt)
            rethrow(e)
        end
    end

    return (psms = vcat([r.psms for r in thread_results]...),)
end

"""
    filter_low_quality_precursors!(weights, Hs, min_frag_count; dynamic_range=1e-3f0)

Zero out weights for precursors that will inevitably fail downstream filters,
before expensive `getDistanceMetrics` computation. Applies two filters:
1. Dynamic range: weights < `dynamic_range * max_weight` are zeroed
2. Isotope-aware fragment count: requires at least `min_frag_count` monoisotopic OR
   `min_frag_count` M+1 matched fragments. This catches PSMs early that would otherwise
   pass total fragment count but fail the later b+y monoisotopic filter.
"""
function filter_low_quality_precursors!(
    weights::Vector{Float32},
    Hs::SparseArray,
    min_frag_count::Int;
    dynamic_range::Float32 = Float32(1e-3)
)
    # Find max weight in this scan
    max_weight = zero(Float32)
    @inbounds for col in 1:Hs.n
        max_weight = max(max_weight, weights[col])
    end
    weight_threshold = dynamic_range * max_weight

    @inbounds for col in 1:Hs.n
        # Dynamic range filter
        if weights[col] < weight_threshold
            weights[col] = zero(Float32)
            continue
        end
        # Isotope-aware fragment count filter
        mono_count = 0
        m1_count = 0
        for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            if Hs.matched[i]
                iso = Hs.isotope[i]
                if iso == 0x00
                    mono_count += 1
                elseif iso == 0x01
                    m1_count += 1
                end
            end
        end
        if mono_count < min_frag_count && m1_count < min_frag_count
            weights[col] = zero(Float32)
        end
    end

    return nothing
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
    min_frag_count::Int64 = params.min_frag_count,
    min_spectral_contrast::Float32 = params.min_spectral_contrast,
    min_log2_matched_ratio::Float32 = params.min_log2_matched_ratio,
    min_topn_of_m::Tuple{Int64, Int64} = params.min_topn_of_m,
    dynamic_range::Float32 = params.dynamic_range,
    first_pass::Bool = false
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    last_val = 0
    cycle_idx = 0

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
            getHighMz(spectra, scan_idx)
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

        converged, _ = solveOLS!(
            Hs,
            residuals,
            weights,
            colnorm2,
            params.max_iter_outer,
            params.max_diff
        )
        if !converged
            reset_arrays!(search_data, Hs)
            continue
        end

        # Update precursor weights and score PSMs
        update_precursor_weights!(search_data, weights, precursor_weights)

        # Filter low-quality precursors before expensive spectral scoring
        filter_low_quality_precursors!(
            weights, Hs, min_frag_count; dynamic_range = dynamic_range
        )

        ScoreFragmentMatches!(
            getComplexUnscoredPsms(search_data),
            getIdToCol(search_data),
            getIonMatches(search_data),
            nmatches,
            getMassErrorModel(search_context, ms_file_idx),
            last(min_topn_of_m)
        )

        if first_pass
            getDistanceMetrics(weights, residuals, Hs, getFirstPassSpectralScores(search_data))
            score_result = Score!(
                getFirstPassScoredPsms(search_data),
                getComplexUnscoredPsms(search_data),
                getFirstPassSpectralScores(search_data),
                weights,
                getIdToCol(search_data),
                cycle_idx,
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(getIntensityArray(spectra, scan_idx))),
                scan_idx;
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_y_count = params.min_y_count,
                min_frag_count = min_frag_count,
                max_best_rank = params.max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000
            )
        else
            getDistanceMetrics(weights, residuals, Hs, getComplexSpectralScores(search_data))
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
                min_spectral_contrast = min_spectral_contrast,
                min_log2_matched_ratio = min_log2_matched_ratio,
                min_y_count = params.min_y_count,
                min_frag_count = min_frag_count,
                max_best_rank = params.max_best_rank,
                min_topn = first(min_topn_of_m),
                block_size = 500000
            )
        end
        last_val = score_result.last_val

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    if first_pass
        return (psms = DataFrame(@view(getFirstPassScoredPsms(search_data)[1:last_val])),)
    else
        return (psms = DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val])),)
    end
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
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    last_val = 0

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
            getHighMz(spectra, scan_idx)
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
        converged, _ = solveOLS!(
            Hs,
            residuals,
            weights,
            colnorm2,
            params.max_iter_outer,
            params.max_diff
        )
        if !converged
            reset_arrays!(search_data, Hs)
            continue
        end

        # Update precursor weights
        update_precursor_weights!(search_data, weights, precursor_weights)

        # Filter low-quality precursors before expensive spectral scoring
        filter_low_quality_precursors!(
            weights, Hs, params.min_frag_count; dynamic_range = params.dynamic_range
        )

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

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    return (psms = DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val])),)
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
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    ion_templates = Vector{Isotope{Float32}}(undef, 100000)
    ion_matches = [PrecursorMatch{Float32}() for _ in range(1, 10000)]
    ion_misses = [UnmatchedIon() for _ in range(1, 10000)]
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
            append!(ion_misses, [UnmatchedIon() for _ in 1:max(ion_idx - length(ion_misses), length(ion_misses))])
        end
        nmatches, nmisses = matchPeaks!(
            ion_matches,
            ion_misses,
            ion_templates,
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx)
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
        
        sort!(@view(ion_matches[1:nmatches]), alg=QuickSort, lt=ion_match_lt)
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
                resize!(colnorm2, length(colnorm2) + new_entries)
                resize!(getMs1SpectralScores(search_data), length(getMs1SpectralScores(search_data)) + new_entries)
                append!(getMs1UnscoredPsms(search_data), [eltype(getMs1UnscoredPsms(search_data))() for _ in 1:new_entries])
            end

            # Initialize active group weights to zero (simple baseline)
            @inbounds @fastmath for col in 1:Int(mz_grouping.current_col)
                weights[col] = 0.0f0
            end

            # Solve deconvolution
            initResiduals!(residuals, Hs, weights)
            solveOLS!(
                Hs,
                residuals,
                weights,
                colnorm2,
                params.max_iter_outer,
                params.max_diff
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
    resize!(getColNorm2(search_data), length(getColNorm2(search_data)) + new_entries)
    resize!(getComplexSpectralScores(search_data), length(getComplexSpectralScores(search_data)) + new_entries)
    resize!(getFirstPassSpectralScores(search_data), length(getFirstPassSpectralScores(search_data)) + new_entries)
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
                        precursors::LibraryPrecursors;
                        prescore_only::Bool=false
)
    # Allocate new columns
   
    #Threads.@threads for i in ProgressBar(range(1, size(psms)[1]))
    N = size(psms, 1);
    decoys = zeros(Bool, N);
    rt = zeros(Float32, N);
    #TIC = zeros(Float16, N);
    total_ions::Vector{UInt16} = psms[!,:total_ions]
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    prec_charges = zeros(UInt8, N);
    #prec_mzs = zeros(Float32, N);
    cv_fold = zeros(UInt8, N);
    scan_idxs::Vector{UInt32} = psms[!,:scan_idx]
    prec_idxs::Vector{UInt32} = psms[!,:precursor_idx]
    error::Vector{Float32} = psms[!,:error]

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
                err_norm[i] = Float16(min((2^min(error[i], 15f0))/max(total_ions[i], one(UInt16)), Float32(6e4)))
                cv_fold[i] = getCvFold(precursors, prec_idx)#prec_id_to_cv_fold[prec_idx]
            end
        end
    end
    fetch.(tasks)
    psms[!,:decoy] = decoys
    psms[!,:rt] = rt
    #psms[!,:TIC] = TIC
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = prec_charges
    #psms[!,:prec_mz] = prec_mzs
    psms[!,:cv_fold] = cv_fold
    psms[!,:charge2] = Vector{UInt8}(psms[!, :charge] .== 2)
    #######
    if !prescore_only
        sort!(psms,:rt); #Sorting before grouping is critical.
    end
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
    get_fraction_transmitted!(chroms, quad_transmission_model, search_data, ...)

Compute only `precursor_fraction_transmitted` without isotope set calculation.
Used for CombineTraces mode where `isotopes_captured` is not needed.
"""
function get_fraction_transmitted!(chroms::DataFrame,
                                   quad_transmission_model::QuadTransmissionModel,
                                   search_data::Vector{SimpleLibrarySearch{IsotopeSplineModel{40, Float32}}},
                                   scan_idx::AbstractVector{UInt32},
                                   prec_charge::AbstractArray{UInt8},
                                   prec_mz::AbstractArray{Float32},
                                   sulfur_count::AbstractArray{UInt8},
                                   centerMz::AbstractVector{Union{Missing, Float32}},
                                   isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    precursor_fraction_transmitted = Vector{Float32}(undef, size(chroms, 1))

    # Cache column vector once (avoid 25M DataFrame row-indexing lookups)
    prec_idx_col = chroms[!, :precursor_idx]

    tasks_per_thread = 5
    chunk_size = max(1, size(chroms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(chroms, 1), chunk_size)

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            thread_id = (first(chunk) % Threads.nthreads()) + 1
            iso_splines = getIsoSplines(search_data[thread_id])
            # Pre-allocate transmission buffer per task (eliminates 25M×zeros(Float32,5) allocations)
            precursor_transmission = zeros(Float32, 5)

            for i in chunk
                prec_id = prec_idx_col[i]
                mz = prec_mz[prec_id]
                charge = prec_charge[prec_id]
                sulfur = sulfur_count[prec_id]
                scan_id = scan_idx[i]
                scan_mz = coalesce(centerMz[scan_id], zero(Float32))::Float32
                window_width = coalesce(isolationWidthMz[scan_id], zero(Float32))::Float32

                precursor_fraction_transmitted[i] = getPrecursorFractionTransmitted!(
                    precursor_transmission,
                    iso_splines,
                    (1,5),
                    getQuadTransmissionFunction(quad_transmission_model, scan_mz, window_width),
                    mz,
                    charge,
                    sulfur)
            end
        end
    end
    fetch.(tasks)
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
                                    rt_to_irt_interp::RtConversionModel;
                                    prescore_only::Bool=false
                                    )

    precursors_lib = getPrecursors(getSpecLib(search_context))
    structural_mods = getStructuralMods(precursors_lib)
    prec_mz = getMz(precursors_lib)
    prec_irt = getIrt(precursors_lib)
    prec_charge = getCharge(precursors_lib)
    entrap_group_ids = getEntrapmentGroupId(precursors_lib)
    precursor_missed_cleavage = getMissedCleavages(precursors_lib)
    precursor_pair_idxs = getPairIdx(precursors_lib)
    prec_length = getLength(precursors_lib)
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
    sequence_length = zeros(UInt8, N);

    # Columns only needed for Phase 2 (full feature set)
    if !prescore_only
        irt_diff = zeros(Float32, N)
        pair_idxs = zeros(UInt32, N)
        entrap_group_id = zeros(UInt8, N)
        adjusted_intensity_explained = zeros(Float16, N);
        prec_charges = zeros(UInt8, N)
        prec_mzs = zeros(Float32, N);
        TIC = zeros(Float16, N);
    end

    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    rt::Vector{Float32} = psms[!,:rt]
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}

    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    function countMOX(seq::Missing)
        return zero(UInt8)
    end

    # Lazy-populate Mox cache: Vector-backed for O(1) lookup, computed on first access per precursor
    n_lib = length(prec_irt)
    _mox_vals = Vector{UInt8}(undef, n_lib)
    _mox_computed = falses(n_lib)

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size)

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = precursor_idx[i]
                irt_obs[i] = rt_to_irt_interp(rt[i])
                irt_pred[i] = prec_irt[prec_idx]
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])
                missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
                # Lazy Mox: compute once per precursor (benign race — same value)
                if !_mox_computed[prec_idx]
                    _mox_vals[prec_idx] = countMOX(structural_mods[prec_idx])
                    _mox_computed[prec_idx] = true
                end
                Mox[i] = _mox_vals[prec_idx]
                spectrum_peak_count[i] = length(masses[scan_idx[i]])
                sequence_length[i] = prec_length[prec_idx]

                if !prescore_only
                    entrap_group_id[i] = entrap_group_ids[prec_idx]
                    irt_diff[i] = abs(irt_obs[i] - prec_irt[prec_idx])
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
    psms[!,:sequence_length] = sequence_length

    if !prescore_only
        psms[!,:irt_diff] = irt_diff
        psms[!,:tic] = TIC
        psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
        psms[!,:charge] = prec_charges
        psms[!,:pair_id] = pair_idxs
        psms[!,:prec_mz] = prec_mzs
        psms[!,:entrapment_group_id] = entrap_group_id
        psms[!,:ms_file_idx] .= ms_file_idx
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
Typed inner function to collect target precursor indices, avoiding dynamic dispatch
on Arrow/DataFrame column types.
"""
function _collect_target_precs!(out::Set{UInt32}, targets::AbstractVector{Bool}, precs::AbstractVector{UInt32})
    @inbounds for i in eachindex(targets, precs)
        if targets[i]
            push!(out, precs[i])
        end
    end
    return out
end

"""
    log_diann_recovery(label, psms, diann_file, diann_global)

Log how many DIA-NN precursors have at least one surviving target PSM.
"""
function log_diann_recovery(label::String, psms::DataFrame,
                            diann_file::Set{UInt32}, diann_global::Set{UInt32})
    (isempty(diann_file) && isempty(diann_global)) && return

    # Unique target precursor_idxs with surviving PSMs
    surviving_precs = Set{UInt32}()
    _collect_target_precs!(surviving_precs, psms[!, :target], psms[!, :precursor_idx])

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
    :fitted_manhattan_distance, :fitted_spectral_contrast, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :max_unmatched_residual, :max_matched_residual, :Mox, :spectrum_peak_count,
    :sequence_length,
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
]

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
2. get_isotopes_captured! — precursor_fraction_transmitted (skipped for prescore)
3. Filter by fraction_transmitted (skipped for prescore)
4. add_features! — irt_error, irt_diff, tic, prec_mz, sequence_length, etc.
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
        getPrecursors(getSpecLib(search_context));
        prescore_only=prescore_only
    )
    t1 = time()

    if prescore_only
        # Skip isotope computation and fraction_transmitted filter for prescore path —
        # these PSMs get filtered in SecondPassSearch anyway, and no PRESCORE_FEATURES use isotope info.
        t2 = t1
        t3 = t1
    else
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
    end

    # 4. Add ML features (irt_error, irt_diff, tic, prec_mz, sequence_length, etc.)
    add_features!(
        psms,
        search_context,
        getTICs(spectra),
        getMzArrays(spectra),
        ms_file_idx,
        getRtIrtModel(search_context, ms_file_idx);
        prescore_only=prescore_only
    )
    t4 = time()

    r = s -> round(s, digits=3)
    @info "  prepare_psm_features! ($(nrow(psms)) PSMs, prescore_only=$prescore_only): " *
          "columns=$(r(t1-t0))s, isotopes=$(r(t2-t1))s, filter=$(r(t3-t2))s, " *
          "add_features=$(r(t4-t3))s, " *
          "total=$(r(t4-t0))s"

    return psms
end

"""
    train_lgbm_and_select_best(psms; features) -> (best_psms, scores, timings)

Train LightGBM on ALL PSMs (all scans), predict scores, select best scan
per precursor by LightGBM score, and log diagnostics.

Returns:
- best_psms: DataFrame with one row per precursor (best by LightGBM score)
- scores: Vector{Float32} of LightGBM probabilities for best_psms
- timings: NamedTuple with timing breakdowns
"""
function train_lgbm_and_select_best(
    psms::DataFrame;
    features::Vector{Symbol} = collect(PRESCORE_FEATURES),
)
    t0 = time()

    # Filter to available features
    available_features = filter(f -> hasproperty(psms, f), features)
    targets_col = psms[!, :target]

    # Build feature matrix ONCE (e.g. 13M×14 Float32)
    X_all = feature_matrix(psms, available_features)
    n_total = size(X_all, 1)
    t_matrix = time()

    # Two-fold cross-validation using existing cv_fold column
    cv_fold = psms[!, :cv_fold]
    idx0 = findall(cv_fold .== 0)
    idx1 = findall(cv_fold .== 1)
    all_scores = Vector{Float64}(undef, n_total)
    max_train = max(1_000_000, n_total ÷ 4)
    last_classifier = nothing

    for (train_idx, test_idx) in [(idx1, idx0), (idx0, idx1)]
        classifier = build_lightgbm_classifier(
            num_iterations = 50,
            learning_rate = 0.2,
            max_depth = 3,
            num_leaves = 10,
            min_data_in_leaf = 5,
            feature_fraction = 0.5,
            bagging_fraction = 0.5,
            bagging_freq = 1,
            is_unbalance = false,
            max_bin = 1023,
        )
        # Subsample training set if > 10M PSMs
        n_train_available = length(train_idx)
        if n_train_available > max_train
            train_idx = train_idx[randperm(n_train_available)[1:max_train]]
            @info "  LightGBM CV fold: subsampled $max_train / $n_train_available training PSMs"
        end
        X_train = X_all[train_idx, :]
        y_train = _prepare_labels(targets_col[train_idx])

        # Check for degenerate case (all same label)
        unique_labels = unique(y_train)
        if length(unique_labels) == 1
            all_scores[test_idx] .= (unique_labels[1] == 0 ? 0.0 : 1.0)
        else
            LightGBM.fit!(classifier, X_train, y_train; verbosity = -1)
            raw = LightGBM.predict(classifier, X_all[test_idx, :])
            all_scores[test_idx] .= ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
            last_classifier = classifier
        end
    end
    model = if last_classifier !== nothing
        LightGBMModel(last_classifier, available_features, nothing)
    else
        LightGBMModel(nothing, available_features, 0.0f0)
    end
    t_train_cv = time()

    # Add scores to psms for best-per-precursor selection
    psms[!, :lgbm_score] = Float32.(all_scores)

    # Select best scan per precursor by LightGBM score
    psms = select_best_per_precursor!(psms, :lgbm_score)

    # Extract scores for best PSMs
    scores = psms[!, :lgbm_score]
    t_best = time()

    # Feature importances
    imp = importance(model)
    if imp !== nothing
        sorted_imp = sort(imp, by = x -> -x[2])
        @user_info "FirstPass LightGBM Feature Importances (gain):"
        for (fname, gain) in sorted_imp
            @user_info "  $(rpad(fname, 40)) $(round(gain, digits=2))"
        end
    end

    timings = (
        matrix = t_matrix - t0,
        train_cv = t_train_cv - t_matrix,
        best = t_best - t_train_cv,
    )

    return psms, Vector{Float32}(scores), timings
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
    select_best_per_precursor!(psms::DataFrame, score_col::Symbol) -> DataFrame

Keeps one row per precursor_idx (highest `score_col`). Also computes `irt_range`
(max - min observed iRT across all scans) per precursor if `:irt_obs` exists.
Returns the filtered DataFrame with the added `irt_range` column.
"""
function select_best_per_precursor!(psms::DataFrame, score_col::Symbol)
    scores = psms[!, score_col]::Vector{Float32}
    prec_ids = psms[!, :precursor_idx]::Vector{UInt32}
    has_irt = hasproperty(psms, :irt_obs)
    irt_obs = has_irt ? psms[!, :irt_obs]::Vector{Float32} : nothing
    n = nrow(psms)

    # Single pass: track best row index, best score, and iRT range per precursor
    best_idx = Dict{UInt32, Int}()
    best_score = Dict{UInt32, Float32}()
    sizehint!(best_idx, n)
    sizehint!(best_score, n)

    min_irt = has_irt ? Dict{UInt32, Float32}() : nothing
    max_irt = has_irt ? Dict{UInt32, Float32}() : nothing
    if has_irt
        sizehint!(min_irt, n)
        sizehint!(max_irt, n)
    end

    @inbounds for i in 1:n
        pid = prec_ids[i]
        s = scores[i]
        if !haskey(best_idx, pid)
            best_idx[pid] = i
            best_score[pid] = s
            if has_irt
                irt = irt_obs[i]
                min_irt[pid] = irt
                max_irt[pid] = irt
            end
        else
            if s > best_score[pid]
                best_idx[pid] = i
                best_score[pid] = s
            end
            if has_irt
                irt = irt_obs[i]
                min_irt[pid] = min(min_irt[pid], irt)
                max_irt[pid] = max(max_irt[pid], irt)
            end
        end
    end

    # Build keep mask and select rows
    keep = falses(n)
    @inbounds for idx in values(best_idx)
        keep[idx] = true
    end

    result = psms[keep, :]

    # Add irt_range column to the reduced table
    if has_irt
        result_pids = result[!, :precursor_idx]::Vector{UInt32}
        irt_range_col = Vector{Float32}(undef, nrow(result))
        @inbounds for i in eachindex(result_pids)
            pid = result_pids[i]
            irt_range_col[i] = max_irt[pid] - min_irt[pid]
        end
        result[!, :irt_range] = irt_range_col
    end

    return result
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

"""
Typed inner function for the per-file aggregation loop in `aggregate_prescore_globally!`.
Accepts Arrow columns as `AbstractVector` so the compiler specializes on concrete column
types at each call site, eliminating dynamic dispatch inside the hot loop.

Accumulates per-file probabilities and tracks the iRT range from the file with the
highest lgbm_prob per precursor (for FWHM estimation).
"""
function _aggregate_file_scores!(
        prec_idx_col::AbstractVector{UInt32},
        calibrated_probs::AbstractVector{Float32},
        irt_range_col::AbstractVector{Float32},
        prec_probs_by_run::Dictionary{UInt32, Vector{Float32}},
        prec_best_prob::Dictionary{UInt32, Float32},
        prec_best_irt_range::Dictionary{UInt32, Float32})

    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        p = calibrated_probs[i]
        if haskey(prec_probs_by_run, pid)
            push!(prec_probs_by_run[pid], p)
            if p > prec_best_prob[pid]
                prec_best_prob[pid] = p
                prec_best_irt_range[pid] = irt_range_col[i]
            end
        else
            insert!(prec_probs_by_run, pid, Float32[p])
            insert!(prec_best_prob, pid, p)
            insert!(prec_best_irt_range, pid, irt_range_col[i])
        end
    end
    return nothing
end

const DEFAULT_GLOBAL_PRESCORE_QVALUE_THRESHOLD = 0.05f0

"""
    aggregate_prescore_globally!(search_context, qvalue_threshold; fold_suffix="") -> (Set{UInt32}, Float32)

Load per-file LightGBM prescore arrow files (with optional `fold_suffix` appended to
filenames), aggregate raw LightGBM probabilities via log-odds averaging of top-√n values,
compute global q-values, and return (1) the set of precursor indices passing at
q ≤ qvalue_threshold and (2) the median iRT range (FWHM proxy) from high-confidence
target precursors (q ≤ 0.1%).

Uses raw log-odds aggregation (no PEP calibration). Simulation showed that isotonic
regression PEP calibration introduces systematic anti-conservative FDR bias because
each sample's own label influences its own PEP estimate (label leakage). Raw log-odds
of CV-scored LightGBM probabilities is perfectly calibrated.
"""
function aggregate_prescore_globally!(search_context::SearchContext,
                                      qvalue_threshold::Float32=DEFAULT_GLOBAL_PRESCORE_QVALUE_THRESHOLD;
                                      fold_suffix::String="")
    r(t) = round(t; digits=2)
    t_total_start = time()

    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")

    # Collect per-file probs and iRT ranges per precursor
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    prec_best_prob = Dictionary{UInt32, Float32}()
    prec_best_irt_range = Dictionary{UInt32, Float32}()
    n_valid_files = 0

    t_reads = 0.0
    t_loop = 0.0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        file_name = getParsedFileName(ms_data, ms_file_idx)
        score_path = joinpath(prescore_dir, "$(file_name)$(fold_suffix).arrow")
        !isfile(score_path) && continue

        t_r = time()
        scores = Arrow.Table(score_path)
        t_reads += time() - t_r
        n_valid_files += 1

        t_l = time()
        _aggregate_file_scores!(
            scores[:precursor_idx], scores[:lgbm_prob], scores[:irt_range],
            prec_probs_by_run, prec_best_prob, prec_best_irt_range
        )
        t_loop += time() - t_l
    end

    # Aggregate via log-odds; look up target/decoy from spectral library
    t_qval_start = time()
    is_decoy = getIsDecoy(getPrecursors(getSpecLib(search_context)))
    sqrt_n = max(1, floor(Int, sqrt(n_valid_files)))
    n_unique = length(prec_probs_by_run)
    global_prec_idxs = Vector{UInt32}(undef, n_unique)
    global_probs = Vector{Float32}(undef, n_unique)
    global_targets = Vector{Bool}(undef, n_unique)

    for (i, (pid, probs)) in enumerate(pairs(prec_probs_by_run))
        global_prec_idxs[i] = pid
        global_probs[i] = _logodds_combine(probs, sqrt_n, 0.1f0)
        global_targets[i] = !is_decoy[pid]
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

    # Estimate FWHM from iRT ranges of high-confidence targets (q ≤ 0.1%)
    fwhm_qval_cutoff = 0.001f0
    hc_irt_ranges = Float32[]
    for i in eachindex(global_prec_idxs)
        if global_qvals[i] <= fwhm_qval_cutoff && global_targets[i]
            pid = global_prec_idxs[i]
            haskey(prec_best_irt_range, pid) && push!(hc_irt_ranges, prec_best_irt_range[pid])
        end
    end
    if !isempty(hc_irt_ranges)
        median_irt_range = Float32(median(hc_irt_ranges))
        n_hc = length(hc_irt_ranges)
        @info "  FWHM estimate: median_irt_range=$(round(median_irt_range, digits=4)) from $n_hc targets at q≤$(fwhm_qval_cutoff)"
    else
        median_irt_range = 0.2f0
        @info "  FWHM estimate: no targets at q≤$(fwhm_qval_cutoff), using default $(median_irt_range)"
    end

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

    t_qval = time() - t_qval_start
    t_total = time() - t_total_start
    @info "Prescore aggregation: file_reads=$(r(t_reads))s, loop=$(r(t_loop))s, combination+qvalues=$(r(t_qval))s, total=$(r(t_total))s"

    return passing, median_irt_range
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

