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

"""
    integrate_precursors(chromatograms::DataFrame, isotope_trace_type::IsotopeTraceType,
                        precursor_idx::AbstractVector{UInt32}, isotopes_captured::AbstractVector{Tuple{Int8, Int8}},
                        apex_scan_idx::AbstractVector{UInt32}, peak_area::AbstractVector{Float32},
                        new_best_scan::AbstractVector{UInt32}, ms_file_idx::Int64;
                        λ::Float32=1.0f0, n_pad::Int64=20, max_apex_offset::Int64=2)

Integrate chromatographic peaks for multiple precursors in parallel.

# Arguments
- `chromatograms`: DataFrame containing chromatogram data
- `min_fraction_transmitted`: Min precursor isotope distribution transmitted
- `isotope_trace_type`: Strategy for handling isotope traces
- `precursor_idx`: Vector of precursor indices
- `isotopes_captured`: Isotopes captured for each precursor
- `apex_scan_idx`: Initial apex scan indices
- `peak_area`: Vector to store integrated areas
- `new_best_scan`: Vector to store updated apex scans
- `ms_file_idx`: MS file index
- `λ`: Smoothing parameter
- `n_pad`: Number of padding points
- `max_apex_offset`: Maximum allowed apex shift

Smooths and integrates chromatograms in parallel batches
"""
function integrate_precursors(chromatograms::DataFrame,
                             isotope_trace_type::IsotopeTraceType,
                             min_fraction_transmitted::Float32,
                             precursor_idx::AbstractVector{UInt32},
                             isotopes_captured::AbstractVector{Tuple{Int8, Int8}},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32},
                             points_integrated::AbstractVector{UInt32},
                             precursor_fraction_transmitted_traces::AbstractVector{String},
                             isotopes_captured_traces::AbstractVector{String},
                             ms_file_idx::Int64; 
                             λ::Float32 = 1.0f0,
                             n_pad::Int64 = 20,
                             max_apex_offset::Int64 = 2,
                             test_print::Bool = false
                             )
    chromatogram_keys = [:precursor_idx]
    if seperateTraces(isotope_trace_type)
        chromatogram_keys = [:precursor_idx,:isotopes_captured]
    end
    grouped_chroms = groupby(chromatograms, chromatogram_keys)
    dtype = Float32
    #Maximal size of a chromatogram
    N = maximum(size(c,1) for c in grouped_chroms) + (2*n_pad)
    group_keys = keys(grouped_chroms)

    if SINGLE_THREAD_INTEGRATE[]
        # Single-threaded path for accurate allocation profiling
        # (@allocated is unreliable in multi-threaded code due to gc_live_bytes() contamination)
        b = zeros(Float32, N)
        u2 = zeros(Float32, length(b))
        ws = WHWorkspace(N)
        state = Chromatogram(zeros(dtype, N), zeros(dtype, N), N)
        _alloc_lookup = Int64(0)
        _alloc_prep = Int64(0)
        _alloc_strings = Int64(0)
        _alloc_integrate = Int64(0)
        _n_precs = 0

        for i in 1:length(precursor_idx)
            prec_id = precursor_idx[i]
            iso_set = isotopes_captured[i]
            apex_scan = apex_scan_idx[i]

            _alloc_lookup += @allocated begin
            if seperateTraces(isotope_trace_type)
                (precursor_idx = prec_id, isotopes_captured = iso_set) ∉ group_keys ? continue : nothing
                chrom = grouped_chroms[(precursor_idx = prec_id, isotopes_captured = iso_set)]
            else
                (precursor_idx = prec_id,) ∉ group_keys ? continue : nothing
                chrom = grouped_chroms[(precursor_idx = prec_id,)]
            end
            end # @allocated lookup

            _alloc_prep += @allocated begin
            avg_cycle_time = (chrom.rt[end] - chrom.rt[1]) /  length(chrom.rt)
            first_pos = findfirst(x->x>0.0, chrom[!,:intensity])
            last_pos = findlast(x->x>0.0, chrom[!,:intensity])
            isnothing(first_pos) ? continue : nothing
            chrom = view(chrom, first_pos:last_pos, :)
            chrom_scan_idx = chrom[!, :scan_idx]::AbstractVector{UInt32}
            if test_print == false
                apex_scan = findfirst(x->x==apex_scan, chrom_scan_idx)
                isnothing(apex_scan) ? continue : nothing
            else
                min_diff = typemax(Int64)
                nearest_idx = i
                for i in range(1, size(chrom, 1))
                    if abs(chrom[i,:scan_idx] - apex_scan) < min_diff
                        min_diff = abs(chrom[i,:scan_idx] - apex_scan)
                        nearest_idx = i
                    end
                end
                for i in range(max(1, nearest_idx - 5), min(size(chrom, 1), nearest_idx + 5))
                    if chrom[i,:intensity] > chrom[nearest_idx,:intensity]
                        nearest_idx = i
                    end
                end
                apex_scan = nearest_idx
            end
            end # @allocated prep

            if !ismissing(precursor_fraction_transmitted_traces)
                _alloc_strings += @allocated begin
                precursor_fraction_transmitted_traces[i], isotopes_captured_traces[i] = get_isolated_isotopes_strings(chrom[!, :precursor_fraction_transmitted], chrom[!, :isotopes_captured])
                end # @allocated strings
            end

            _alloc_integrate += @allocated begin
            peak_area[i], new_best_scan[i], points_integrated[i] = integrate_chrom(
                            chrom,
                            apex_scan,
                            b,
                            u2,
                            ws,
                            state,
                            avg_cycle_time,
                            λ,
                            n_pad = n_pad,
                            max_apex_offset = max_apex_offset,
                            isplot = false
                            );
            end # @allocated integrate

            _n_precs += 1
            reset!(state)
        end
        @info "  integrate_precursors SINGLE-THREAD ($_n_precs precs): lookup=$(round(_alloc_lookup/1e9,digits=3))GB prep=$(round(_alloc_prep/1e9,digits=3))GB strings=$(round(_alloc_strings/1e9,digits=3))GB integrate=$(round(_alloc_integrate/1e9,digits=3))GB total=$(round((_alloc_lookup+_alloc_prep+_alloc_strings+_alloc_integrate)/1e9,digits=3))GB"
        return nothing
    end

    thread_tasks = partitionThreadTasks(length(precursor_idx), 10, Threads.nthreads())
    tasks = map(thread_tasks) do chunk
        Threads.@spawn begin
            #chromdf = DataFrame()
            b = zeros(Float32, N);
            u2 = zeros(Float32, length(b));
            ws = WHWorkspace(N)

            state = Chromatogram(
                zeros(dtype, N), #t
                zeros(dtype, N), #data
                N #max index
                )

            for i in chunk
                prec_id = precursor_idx[i]
                iso_set = isotopes_captured[i]
                apex_scan = apex_scan_idx[i]

                if seperateTraces(isotope_trace_type)
                    (precursor_idx = prec_id, isotopes_captured = iso_set) ∉ group_keys ? continue : nothing
                    chrom = grouped_chroms[(precursor_idx = prec_id, isotopes_captured = iso_set)]
                else
                    (precursor_idx = prec_id,) ∉ group_keys ? continue : nothing
                    chrom = grouped_chroms[(precursor_idx = prec_id,)]
                end

                avg_cycle_time = (chrom.rt[end] - chrom.rt[1]) /  length(chrom.rt)
                first_pos = findfirst(x->x>0.0, chrom[!,:intensity])
                last_pos = findlast(x->x>0.0, chrom[!,:intensity])
                isnothing(first_pos) ? continue : nothing
                chrom = view(chrom, first_pos:last_pos, :)
                chrom_scan_idx = chrom[!, :scan_idx]::AbstractVector{UInt32}
                if test_print == false
                    apex_scan = findfirst(x->x==apex_scan, chrom_scan_idx)
                    isnothing(apex_scan) ? continue : nothing
                else
                    min_diff = typemax(Int64)
                    nearest_idx = i
                    for i in range(1, size(chrom, 1))
                        if abs(chrom[i,:scan_idx] - apex_scan) < min_diff
                            min_diff = abs(chrom[i,:scan_idx] - apex_scan)
                            nearest_idx = i
                        end
                    end
                    for i in range(max(1, nearest_idx - 5), min(size(chrom, 1), nearest_idx + 5))
                        if chrom[i,:intensity] > chrom[nearest_idx,:intensity]
                            nearest_idx = i
                        end
                    end
                    apex_scan = nearest_idx
                end

                if !ismissing(precursor_fraction_transmitted_traces)
                    precursor_fraction_transmitted_traces[i], isotopes_captured_traces[i] = get_isolated_isotopes_strings(chrom[!, :precursor_fraction_transmitted], chrom[!, :isotopes_captured])
                end

                peak_area[i], new_best_scan[i], points_integrated[i] = integrate_chrom(
                                chrom,
                                apex_scan,
                                b,
                                u2,
                                ws,
                                state,
                                avg_cycle_time,
                                λ,
                                n_pad = n_pad,
                                max_apex_offset = max_apex_offset,
                                isplot = false
                                );

                reset!(state)
            end
            return
        end
    end
    return fetch.(tasks)
end

#==========================================================
Chromatogram Building Functions
==========================================================#
"""
    extract_chromatograms(spectra::MassSpecData, passing_psms::DataFrame,
                         rt_index::retentionTimeIndex, search_context::SearchContext,
                         params::IntegrateChromatogramSearchParameters,
                         ms_file_idx::Int64) -> DataFrame

Extract chromatograms for all passing PSMs.

Uses parallel processing to build chromatograms across scan ranges.
Set `SINGLE_THREAD_CHROM_EXTRACT[] = true` to run single-threaded for profiling.
"""
const SINGLE_THREAD_CHROM_EXTRACT = Ref(false)
const SINGLE_THREAD_INTEGRATE = Ref(false)

function extract_chromatograms(
    spectra::MassSpecData,
    passing_psms::DataFrame,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    chrom_type::CHROMATOGRAM
)
    if typeof(chrom_type)==typeof(MS2CHROM())
        ms_order_select = 2
    else
        ms_order_select = 1
    end

    thread_tasks = partition_scans(spectra, Threads.nthreads(), ms_order_select = ms_order_select)

    # Pre-create Sets for each thread to avoid concurrent DataFrame access
    # This eliminates race conditions when multiple threads call passing_psms[!, :precursor_idx]
    precursor_sets = [Set(passing_psms[!, :precursor_idx]) for _ in 1:Threads.nthreads()]

    if SINGLE_THREAD_CHROM_EXTRACT[]
        # Single-threaded: run all partitions sequentially, each with its own search_data
        results = Vector{DataFrame}(undef, length(thread_tasks))
        for (idx, thread_task) in enumerate(thread_tasks)
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            results[idx] = build_chromatograms(
                spectra,
                last(thread_task),
                precursor_sets[thread_id],
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx,
                chrom_type
            )
        end
        return vcat(results...)
    else
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin
                thread_id = first(thread_task)
                search_data = getSearchData(search_context)[thread_id]

                return build_chromatograms(
                    spectra,
                    last(thread_task),
                    precursor_sets[thread_id],
                    rt_index,
                    search_context,
                    search_data,
                    params,
                    ms_file_idx,
                    chrom_type
                )
            end
        end
        thread_dfs = fetch.(tasks)
        return vcat(thread_dfs...)
    end
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
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    # Initialize working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    # Pre-allocate chromatograms with better size estimate (~100 points per scan average)
    estimated_points = length(scan_range) * 100
    chromatograms = Vector{MS2ChromObject}(undef, max(estimated_points, 10000))

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz = zero(Float32)
    ion_idx = 0
    rt_idx = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    nce_model = getNceModel(search_context, ms_file_idx)

    # Cache accessors that are constant for the entire function
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    mass_error_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)
    frag_lookup = getFragmentLookupTable(spec_lib)
    ion_templates = getIonTemplates(search_data)
    ion_matches = getIonMatches(search_data)
    ion_misses = getIonMisses(search_data)
    id_to_col = getIdToCol(search_data)
    spectral_scores = getSpectralScores(search_data)
    unscored_psms = getUnscoredPsms(search_data)

    i = 1
    for scan_idx in scan_range
        ((scan_idx<1) | (scan_idx > length(spectra))) && continue
        # Process MS2 scans
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && continue

        # Calculate RT window
        irt = rt_irt_model(getRetentionTime(spectra, scan_idx))
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Check for m/z change (numeric comparison replaces string allocation)
        prec_mz_new = getCenterMz(spectra, scan_idx)

        # Update transitions if window changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_new != prec_mz)
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz = prec_mz_new
            prec_temp_size = 0
            quad_func = getQuadTransmissionFunction(
                quad_model,
                getCenterMz(spectra, scan_idx),
                getIsolationWidthMz(spectra, scan_idx)
            )

            ion_idx, prec_temp_size = selectTransitions!(
                ion_templates,
                RTIndexedTransitionSelection(),
                params.prec_estimation,
                frag_lookup,
                nce_model,
                precs_temp,
                prec_mz_arr,
                prec_charge_arr,
                prec_sulfur_arr,
                getIsoSplines(search_data),
                quad_func,
                getPrecursorTransmission(search_data),
                getIsotopes(search_data),
                params.n_frag_isotopes,
                params.max_frag_rank,
                rt_index,
                irt_start,
                irt_stop,
                (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
                precursors_passing = precursors_passing,
                isotope_err_bounds = params.isotope_err_bounds,
                block_size = 10000,
                min_fraction_transmitted = params.min_fraction_transmitted
            )
        end

        # Match peaks
        nmatches, nmisses = matchPeaks!(
            ion_matches,
            ion_misses,
            ion_templates,
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mass_error_model,
            getHighMz(spectra, scan_idx)
        )

        sort!(@view(ion_matches[1:nmatches]), alg=QuickSort, lt=ion_match_lt)

        # Process matches
        if nmatches > 2
            i += 1
            # Build design matrix + deconvolution
            buildDesignMatrix!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                id_to_col
            )

            # Handle array resizing
            if id_to_col.size > length(weights)
                new_entries = id_to_col.size - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
                resize!(spectral_scores, length(spectral_scores) + new_entries)
                # Avoid list comprehension allocation - use direct resize and loop
                old_length = length(unscored_psms)
                resize!(unscored_psms, old_length + new_entries)
                for i in (old_length + 1):length(unscored_psms)
                    unscored_psms[i] = eltype(unscored_psms)()
                end
            end

            # Initialize weights
            for i in 1:id_to_col.size
                weights[id_to_col[id_to_col.keys[i]]] =
                    precursor_weights[id_to_col.keys[i]]
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

            # Record chromatogram points with weights
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end

                col = id_to_col[precs_temp[j]]
                if !iszero(col)
                    chromatograms[rt_idx] = MS2ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        weights[col],
                        scan_idx,
                        precs_temp[j]
                    )
                else
                    chromatograms[rt_idx] = MS2ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        zero(Float32),
                        scan_idx,
                        precs_temp[j]
                    )
                end
            end

            # Update precursor weights
            for i in 1:id_to_col.size
                precursor_weights[id_to_col.keys[i]] =
                    weights[id_to_col[id_to_col.keys[i]]]
            end

        else
            # Record zero intensity points for non-matching precursors
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end

                chromatograms[rt_idx] = MS2ChromObject(
                    Float32(getRetentionTime(spectra, scan_idx)),
                    zero(Float32),
                    scan_idx,
                    precs_temp[j]
                )
            end
        end

        # Reset arrays
        for i in 1:Hs.n
            unscored_psms[i] = eltype(unscored_psms)()
        end
        reset!(id_to_col)
        reset!(Hs)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
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
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS1CHROM
)
    # Initialize working arrays
    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    chromatograms = Vector{MS1ChromObject}(undef, 500000)  # Initial size
    ion_templates = Vector{Isotope{Float32}}(undef, 100000)
    ion_matches = [PrecursorMatch{Float32}() for _ in range(1, 10000)]
    ion_misses = [UnmatchedIon() for _ in range(1, 10000)]
    precursors = getPrecursors(getSpecLib(search_context))
    seqs = [getSequence(precursors)[pid] for pid in precursors_passing]
    pids = [pid for pid in precursors_passing]
    pcharge = [getCharge(precursors)[pid] for pid in precursors_passing]
    pmz = [getMz(precursors)[pid] for pid in precursors_passing]
    isotopes_dict = getIsotopes(seqs, pmz, pids, pcharge, QRoots(5), 5)

    # NEW: Create m/z grouping map for MS1
    mz_grouping = MzGroupingMap(UInt32(100000))  # 5 decimal place precision

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    ion_idx = 0
    rt_idx = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    id_to_col = getIdToCol(search_data)
    spectral_scores = getSpectralScores(search_data)
    unscored_psms = getUnscoredPsms(search_data)
    i = 1
    for scan_idx in scan_range

        ((scan_idx<1) | (scan_idx > length(spectra))) && continue
        # Process MS1 scans
        #msn = getMsOrder(spectra, scan_idx)
        #msn ∉ one(UInt8) && continue
        if getMsOrder(spectra, scan_idx) != 1
            continue
        end
        iso_count = Dictionary{UInt32, @NamedTuple{matched_mono::Bool, iso_count::UInt8}}()
        # Calculate RT window
        irt = rt_irt_model(getRetentionTime(spectra, scan_idx))
        irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Update transitions if window changed
        prec_temp_size = 0
        ion_idx = 0
        for rt_bin_idx in irt_start:irt_stop
            precs = rt_index.rt_bins[rt_bin_idx].prec
            for i in 1:length(precs)
                prec_idx = first(precs[i])
                if prec_idx in precursors_passing
                    prec_temp_size += 1
                    if prec_temp_size > length(precs_temp)
                        append!(precs_temp, Vector{UInt32}(undef, 1000))
                    end
                    precs_temp[prec_temp_size] = prec_idx
                    for iso in isotopes_dict[prec_idx]
                        ion_idx += 1
                        ion_templates[ion_idx] = iso
                    end
                end
            end
        end
        #Probably more efficient way to do this later 
        for i in range(1, ion_idx)
            _ion_ = ion_templates[i]
            pid = getPrecID(_ion_)
            if haskey(iso_count, pid)
                matched_mono = false
                if iso_count[pid].matched_mono
                    matched_mono = true
                elseif getIsoIdx(_ion_)==one(UInt8)
                    matched_mono = true
                end
                iso_count[pid] = (matched_mono = matched_mono, iso_count = iso_count[pid].iso_count + one(UInt8))
            else
                insert!(iso_count, 
                pid,
                (matched_mono = getIsoIdx(_ion_)==one(UInt8), iso_count = one(UInt8))
            )
            end
            iso_count
        end
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
            getHighMz(spectra, scan_idx)
        )

        #nmisses -= 1
        sort!(@view(ion_matches[1:nmatches]), alg=QuickSort, lt=ion_match_lt)
        #println("nmatches $nmatches nmisses $nmisses")
        # Process matches
        if nmatches > 2
            i += 1

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

            # Handle array resizing
            if id_to_col.size > length(weights)
                new_entries = id_to_col.size - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
                resize!(spectral_scores, length(spectral_scores) + new_entries)
                # Avoid list comprehension allocation - use direct resize and loop
                old_length = length(unscored_psms)
                resize!(unscored_psms, old_length + new_entries)
                for i in (old_length + 1):length(unscored_psms)
                    unscored_psms[i] = eltype(unscored_psms)()
                end
            end

            # Initialize weights
            for i in 1:id_to_col.size
                weights[id_to_col[id_to_col.keys[i]]] =
                    precursor_weights[id_to_col.keys[i]]
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

            # Record chromatogram points with weights
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end

                # Get weight from precursor_weights array (now contains distributed group coefficients)
                prec_id = precs_temp[j]
                weight = prec_id <= length(precursor_weights) ? precursor_weights[prec_id] : 0.0f0

                if weight > 0.0f0
                    chromatograms[rt_idx] = MS1ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        weight,
                        iso_count[prec_id].matched_mono,
                        iso_count[prec_id].iso_count,
                        scan_idx,
                        prec_id
                    )
                else
                    chromatograms[rt_idx] = MS1ChromObject(
                        Float32(getRetentionTime(spectra, scan_idx)),
                        zero(Float32),
                        false,
                        UInt8(0),
                        scan_idx,
                        precs_temp[j]
                    )
                end
            end

            # Update precursor weights - already handled by distribute_ms1_coefficients!
            # No need to update here since distribution already populated precursor_weights

        else
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)  # Exponential growth
                end
                chromatograms[rt_idx] = MS1ChromObject(
                    Float32(getRetentionTime(spectra, scan_idx)),
                    zero(Float32),
                    false,
                    UInt8(0),
                    scan_idx,
                    precs_temp[j]
                )
            end
        end

        # Reset arrays
        for i in 1:Hs.n
            unscored_psms[i] = eltype(unscored_psms)()
        end
        reset!(id_to_col)
        reset!(Hs)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
end

#==========================================================
Chromatogram Building Functions
==========================================================#
"""
    process_final_psms!(psms::DataFrame, search_context::SearchContext,
                       parsed_fname::String, ms_file_idx::Int64)

Process final PSMs after integration.

# Added Columns
- Protein information (accession numbers, indices)
- Peak areas and normalization
- Sequence information
- Modification details
- File metadata
"""
function process_final_psms!(
    psms::DataFrame,
    search_context::SearchContext,
    parsed_fname::String,
    ms_file_idx::Int64
)
    # Remove invalid peak areas
    filter!(row -> !isnan(row.peak_area::Float32), psms)
    filter!(row -> row.peak_area::Float32 > 0.0, psms)
    # Add columns
    precursors = getPrecursors(getSpecLib(search_context))
    n = size(psms, 1)
    accession_numbers = Vector{String}(undef, n)
    ms_file_idxs = Vector{UInt16}(undef, n)
    species = Vector{String}(undef, n)
    peak_area = Vector{Union{Missing, Float32}}(undef, n)
    peak_area_normalized = Vector{Union{Missing, Float32}}(undef, n)
    structural_mods = Vector{Union{Missing, String}}(undef, n)
    isotopic_mods = Vector{Union{Missing, String}}(undef, n)
    charge = Vector{UInt8}(undef, n)
    sequence = Vector{String}(undef, n)
    file_name = Vector{String}(undef, n)
    psms_precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}

    for i in range(1, n)
        pid = psms_precursor_idx[i]
        accession_numbers[i] = getAccessionNumbers(precursors)[pid]
    end
    psms[!, :accession_numbers] = accession_numbers
    
    #Critical to sort correctly for the batch-wise MaxLFQ algorithm. 
    #Otherwise could have different precursors for the same protein group 
    #split between two batches. 
    fast_df_sort!(psms, [:inferred_protein_group, :precursor_idx])

    parsed_fname = getFileIdToName(getMSData(search_context), ms_file_idx)
    for i in range(1, n)
        pid = psms_precursor_idx[i]
        ms_file_idxs[i] = UInt32(ms_file_idx)
        species[i] = join(sort(unique(split(coalesce(getProteomeIdentifiers(precursors)[pid], ""),';'))),';')
        peak_area[i] = psms[i,:peak_area]
        peak_area_normalized[i] = zero(Float32)
        structural_mods[i] = getStructuralMods(precursors)[pid]
        isotopic_mods[i] = getIsotopicMods(precursors)[pid]
        charge[i] = getCharge(precursors)[pid]
        sequence[i] = getSequence(precursors)[pid]
        file_name[i] = parsed_fname
    end

    psms[!,:ms_file_idx] = ms_file_idxs
    psms[!,:species] = species
    psms[!,:peak_area] = peak_area
    psms[!,:peak_area_normalized] = peak_area_normalized
    psms[!,:structural_mods] = structural_mods
    psms[!,:isotopic_mods] = isotopic_mods
    psms[!,:charge] = charge 
    psms[!,:sequence] = sequence
    psms[!,:file_name] = file_name
    
    return nothing
end




"""
    unique_pair_strings(pt::AbstractVector,
                        iso::AbstractVector) -> (String, String)

`pt`  - precursor_fraction_transmitted  (e.g. Vector{Float32})  
`iso` - isotopes_captured                (e.g. Vector{Tuple{Int8,Int8}})

Returns a pair of semicolon-separated strings that list **unique** `(pt,iso)`
pairs in the *order of their first appearance*.
"""
function get_isolated_isotopes_strings(pt::AbstractVector, iso::AbstractVector)
    @assert length(pt) == length(iso)

    seen   = Dict{Tuple{Any,Any},Bool}()   # tracks already‑emitted pairs
    uniq_pt  = String[]                    # keeps order
    uniq_iso = String[]

    sig2 = x -> @sprintf("%.2g", x)

    for (p, t) in zip(pt, iso)
        key = (p,t)
        haskey(seen, key) && continue      # skip duplicates

        push!(uniq_pt,  sig2.(p))         # or @sprintf("%.3f", p) for fixed format
        push!(uniq_iso, string(t))         # Tuple prints as "(0,1)" etc.
        seen[key] = true
    end

    order = reverse(sortperm(uniq_pt))

    return join(uniq_pt[order], ';'), join(uniq_iso[order], ';')
end