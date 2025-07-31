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
    
    # Collect results
    results = fetch.(tasks)
    return vcat([r for r in results]...)
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
    
    return vcat(fetch.(tasks)...)
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
    

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    ion_idx = 0
    cycle_idx = 0

    irt_tol = getIrtErrors(search_context)[ms_file_idx]

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
                block_size = 10000
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

        last_val = Score!(
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

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
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
    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    ion_idx = 0
    last_val = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]

    for scan_idx in scan_range
        empty!(pair_id_dict)
        ((scan_idx<1) | (scan_idx > length(spectra))) && continue

        # Process MS1 scans
        if getMsOrder(spectra, scan_idx) != 1
            continue
        end

        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
        irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Update transitions if window changed
        prec_temp_size = 0
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
                    prec_temp_size += 1 #Precursors in the scan 
                    if prec_temp_size > length(precs_temp) #Adjust size of placeholder if necessary 
                        append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
                    end
                    precs_temp[prec_temp_size] = prec_idx


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
            buildDesignMatrix!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                getIdToCol(search_data)
            )

            # Handle array resizing
            if getIdToCol(search_data).size > length(weights)
                new_entries = getIdToCol(search_data).size - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
                append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])

                resize!(getMs1SpectralScores(search_data), length(getMs1SpectralScores(search_data)) + new_entries)
                append!(getMs1UnscoredPsms(search_data), [eltype(getMs1UnscoredPsms(search_data))() for _ in 1:new_entries])
            end

            # Initialize weights
            for i in 1:getIdToCol(search_data).size
                weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
                    precursor_weights[getIdToCol(search_data).keys[i]]
            end

            # Solve deconvolution
            initResiduals!(residuals, Hs, weights)
            solveHuber!(
                Hs,
                residuals,
                weights,
                Float32(1e8),#getHuberDelta(search_context),
                Float32(0.001),#params.lambda,
                params.max_iter_newton,
                params.max_iter_bisection,
                1000,#params.max_iter_outer,
                search_context.deconvolution_stop_tolerance[],#params.accuracy_newton,
                search_context.deconvolution_stop_tolerance[],#params.accuracy_bisection,
                params.max_diff,
                NoNorm()
            )
            # Update precursor weights
            for i in 1:getIdToCol(search_data).size
                precursor_weights[getIdToCol(search_data).keys[i]] = 
                    weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]]
            end

            getDistanceMetrics(weights, residuals, Hs, getMs1SpectralScores(search_data))

            ScoreFragmentMatches!(
                getMs1UnscoredPsms(search_data),
                getIdToCol(search_data),
                ion_matches,
                nmatches,
                getMassErrorModel(search_context, ms_file_idx),
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
                err_norm[i] = min(Float16((2^error[i])/(total_ions[i])), 6e4)
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
                sequence_length[i] = length(replace(precursor_sequence[prec_idx], r"\(.*?\)" => ""))#replace.(sequence[i], "M(ox)" => "M");
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
    return nothing
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

    apex_scan = argmax(psms[!,:weight])
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
    spectra::MassSpecData
)
    if !hasproperty(psms, :precursor_idx) || (size(psms, 1) == 0)
        return DataFrame()
    end
    rts = zeros(Float32, size(psms, 1))
    for i in range(1, size(psms, 1))
        scan_idx = psms[i, :scan_idx]
        rts[i] = getRetentionTime(spectra, scan_idx)
    end
    psms[!,:rt] = rts

    #Group the psms by precursor_idx and return apex scan 
    return combine(groupby(psms,:precursor_idx)) do group
        # Find the row with the maximum value
        idx = argmax(group[!, :weight])
        # Return that row as a 1-row DataFrame
        return group[idx, :]
    end
end
