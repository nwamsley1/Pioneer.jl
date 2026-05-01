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
    build_chrom_index(prec_col::AbstractVector{UInt32}, n_rows::Int)

Build a lookup from precursor ID → row range in a sorted chromatogram table.

The chromatogram DataFrame is sorted by `(precursor_idx, rt)`, so all rows for
a given precursor are contiguous. This function does a single pass to record
each precursor's `start:stop` row range and tracks the longest chromatogram
(used to size scratch buffers).

Returns `(chrom_index, max_chrom_len)`.
"""
function build_chrom_index(prec_col::AbstractVector{UInt32}, n_rows::Int)
    if n_rows == 0
        return Dict{UInt32, UnitRange{Int}}(), 0
    end
    max_len = 0
    chrom_index = Dict{UInt32, UnitRange{Int}}()
    start_idx = 1
    current_key = prec_col[1]
    for row_idx in 2:n_rows
        next_key = prec_col[row_idx]
        if next_key != current_key
            chrom_index[current_key] = start_idx:(row_idx - 1)
            max_len = max(max_len, row_idx - start_idx)
            start_idx = row_idx
            current_key = next_key
        end
    end
    chrom_index[current_key] = start_idx:n_rows
    max_len = max(max_len, n_rows - start_idx + 1)
    return chrom_index, max_len
end

"""
    find_nearest_scan(scan_indices::AbstractVector{UInt32}, target_scan::UInt32)

Return the 1-based position within `scan_indices` whose value is closest to
`target_scan`. Used to map the PSM's apex scan index (global scan number) into
the local chromatogram coordinate after trimming zero-intensity edges.
"""
function find_nearest_scan(scan_indices::AbstractVector{UInt32}, target_scan::UInt32)
    nearest_idx = 1
    min_diff = typemax(Int64)
    for j in eachindex(scan_indices)
        d = abs(Int64(scan_indices[j]) - Int64(target_scan))
        if d < min_diff
            min_diff = d
            nearest_idx = j
        end
    end
    return nearest_idx
end

"""
    integrate_precursors(chromatograms, min_fraction_transmitted, precursor_idx,
                         apex_scan_idx, peak_area, new_best_scan, points_integrated;
                         λ=1.0f0)

Integrate chromatographic peaks for multiple precursors in parallel.

For each precursor, looks up its chromatogram rows via `build_chrom_index`,
trims leading/trailing zero-intensity scans, maps the PSM apex into the trimmed
coordinate, then calls `integrate_chrom` (Whittaker-Henderson smooth → boundary
detection → trapezoidal integration). Results are written into the output vectors
`peak_area`, `new_best_scan`, and `points_integrated`.
"""
function integrate_precursors(chromatograms::DataFrame,
                             min_fraction_transmitted::Float32,
                             precursor_idx::AbstractVector{UInt32},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32},
                             points_integrated::AbstractVector{UInt32};
                             λ::Float32 = 1.0f0
                             )
    n_pad = Int64(0)
    prec_col = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    rt_all = chromatograms[!, :rt]::AbstractVector{Float32}
    scan_idx_all = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    intensity_all = chromatograms[!, :intensity]::AbstractVector{Float32}
    fraction_all = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}

    chrom_index, max_chrom_len = build_chrom_index(prec_col, nrow(chromatograms))
    N = max_chrom_len + (2*n_pad)

    # Partition work into chunks, then distribute chunks to exactly nthreads() tasks
    # (one task per thread — prevents workspace corruption from task switching at yield points)
    _n_threads = Threads.nthreads()
    all_chunks = collect(partitionThreadTasks(length(precursor_idx), 10, _n_threads))
    n_tasks = min(_n_threads, length(all_chunks))

    task_ranges = [Int[] for _ in 1:n_tasks]
    for (chunk_idx, _) in enumerate(all_chunks)
        push!(task_ranges[((chunk_idx - 1) % n_tasks) + 1], chunk_idx)
    end

    n_thread_slots = Threads.maxthreadid()
    ws_by_thread = [WHWorkspace(N) for _ in 1:n_thread_slots]
    state_by_thread = [Chromatogram(zeros(Float32, N), zeros(Float32, N), 0) for _ in 1:n_thread_slots]

    function run_integration_batch!(batch_id::Int, ws::WHWorkspace, state::Chromatogram)
        for chunk_idx in task_ranges[batch_id]
            chunk = all_chunks[chunk_idx]
            for i in chunk
                prec_id = precursor_idx[i]
                apex_scan = apex_scan_idx[i]

                if !haskey(chrom_index, prec_id)
                    continue
                end
                chrom_range = chrom_index[prec_id]

                avg_cycle_time = (rt_all[last(chrom_range)] - rt_all[first(chrom_range)]) / length(chrom_range)

                # Skip entirely zero-valued chromatograms
                if isnothing(findfirst(x -> x > 0.0f0, @view intensity_all[chrom_range]))
                    continue
                end

                # Map global apex scan index to local position in chromatogram
                apex_scan = find_nearest_scan(
                    @view(scan_idx_all[chrom_range]), apex_scan
                )

                peak_area[i], new_best_scan[i], points_integrated[i] = integrate_chrom(
                    @view(rt_all[chrom_range]),
                    @view(scan_idx_all[chrom_range]),
                    @view(intensity_all[chrom_range]),
                    @view(fraction_all[chrom_range]),
                    apex_scan,
                    ws,
                    state,
                    avg_cycle_time,
                    λ,
                    min_fraction_transmitted = min_fraction_transmitted,
                    n_pad = n_pad
                )
                reset!(state)
            end
        end
        return nothing
    end

    Threads.@threads :static for batch_id in 1:n_tasks
        tid = Threads.threadid()
        run_integration_batch!(batch_id, ws_by_thread[tid], state_by_thread[tid])
    end

    # Clamp NaN and negative values to zero for downstream processing
    for i in eachindex(precursor_idx)
        if isnan(peak_area[i]) || peak_area[i] < 0f0
            peak_area[i] = 0f0
        end
    end

    return nothing
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
"""

# Quad-window helpers (originally in selectTransitions/rtIndexTransitionSelection.jl;
# kept here as the sole remaining caller after the fused migration).

function getQuadrupoleBounds(quad_func::QuadTransmissionFunction)
    return getPrecMinBound(quad_func), getPrecMaxBound(quad_func)
end

function getPrecursorWindowRange(
    precs, min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{I, I}
) where {I<:Integer}
    start = searchsortedfirst(precs, by = x->last(x),
        min_prec_mz - first(isotope_err_bounds)*NEUTRON/2)
    stop = searchsortedlast(precs, by = x->last(x),
        max_prec_mz + last(isotope_err_bounds)*NEUTRON/2)
    return start, stop
end

function withinQuadrupoleBounds(
    prec_mz::Float32, prec_charge::UInt8,
    min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{I, I}
) where {I<:Integer}
    mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
    mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge
    return mz_low ≤ prec_mz ≤ mz_high
end

"""
    collect_rt_window_precursors!(precs_temp, rt_index, rt_start_idx, rt_stop_idx,
                                   precursors_passing, prec_mzs, prec_charges,
                                   prec_sulfur_counts, iso_splines, quad_func,
                                   precursor_transmission, isotope_err_bounds,
                                   min_fraction_transmitted, precursor_rt_map,
                                   scan_rt, rt_binned_tol, rt_tol_fallback) -> Int

Walk the RT-bin range, applying quad-window, precursors_passing allowlist,
per-precursor RT, isotope_err_bounds, and min_fraction_transmitted filters.
Writes the surviving precursor ids into `precs_temp[1:n]` and returns `n`.

Extracted from the classic `RTIndexedTransitionSelection` path so
`run_fused!(FusedRTIndexed(), ...)` can consume a pre-filtered precursor
list without repeating the filtering (skipped via `check_prec_filters`).
"""
function collect_rt_window_precursors!(
    precs_temp::Vector{UInt32},
    rt_index::retentionTimeIndex,
    rt_start_idx::Int64, rt_stop_idx::Int64,
    precursors_passing::Union{Set{UInt32}, Nothing},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotope_err_bounds::Tuple{I, I},
    min_fraction_transmitted::Float32,
    precursor_rt_map::Union{Dict{UInt32, Float32}, Nothing},
    scan_rt::Float32,
    rt_binned_tol::Union{RTBinnedTolerance, Nothing},
    rt_tol_fallback::Float32) where {I<:Integer}

    min_prec_mz, max_prec_mz = getQuadrupoleBounds(quad_transmission_func)
    size = 0
    has_rt_filter = precursor_rt_map !== nothing

    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec
        start, stop = getPrecursorWindowRange(precs, min_prec_mz, max_prec_mz, isotope_err_bounds)

        for i in start:stop
            prec_idx = first(precs[i])
            (!isnothing(precursors_passing) && prec_idx ∉ precursors_passing) && continue

            if has_rt_filter
                prec_rt_val = get(precursor_rt_map, prec_idx, NaN32)
                if !isnan(prec_rt_val)
                    prec_tol = rt_binned_tol !== nothing ? get_rt_tol(rt_binned_tol, prec_rt_val) : rt_tol_fallback
                    abs(scan_rt - prec_rt_val) > prec_tol && continue
                end
            end

            prec_charge = prec_charges[prec_idx]
            prec_mz = prec_mzs[prec_idx]
            prec_sulfur_count = prec_sulfur_counts[prec_idx]

            !withinQuadrupoleBounds(prec_mz, prec_charge, min_prec_mz, max_prec_mz, isotope_err_bounds) && continue

            fraction_transmitted = getPrecursorFractionTransmitted!(
                precursor_transmission, iso_splines, (1, 5),
                quad_transmission_func, prec_mz, prec_charge, prec_sulfur_count)
            fraction_transmitted < min_fraction_transmitted && continue

            prec_isotope_set = getPrecursorIsotopeSet(prec_mz, prec_charge, quad_transmission_func)
            last(prec_isotope_set) < 0 && continue

            size += 1
            if size > length(precs_temp)
                append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
            end
            precs_temp[size] = prec_idx
        end
    end
    return size
end

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

    # Build precursor RT map for per-precursor symmetric window filtering
    _pids = passing_psms[!, :precursor_idx]::Vector{UInt32}
    _rts = passing_psms[!, :rt]::Vector{Float32}
    precursor_rt_map = Dict{UInt32, Float32}()
    sizehint!(precursor_rt_map, length(_pids))
    for i in eachindex(_pids)
        precursor_rt_map[_pids[i]] = _rts[i]
    end

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            return build_chromatograms(
                spectra,
                last(thread_task),
                precursor_sets[thread_id],
                precursor_rt_map,
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
    precursor_rt_map::Dict{UInt32, Float32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS2CHROM
)
    # Fused-kernel working arrays.
    Hs = getHsFused(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    fused_scratch = getFusedScratch(search_data)
    corr_mz = getScanCorrectedMz(search_data)
    obs_low = getScanObsLow(search_data)
    obs_high = getScanObsHigh(search_data)
    isotopes_buf = getIsotopes(search_data)
    prec_trans_buf = getPrecursorTransmission(search_data)
    id_to_col = getIdToCol(search_data)
    unscored_psms = getComplexUnscoredPsms(search_data)

    # Pre-allocate chromatograms with better size estimate (~100 points per scan average)
    estimated_points = length(scan_range) * 100
    chromatograms = Vector{MS2ChromObject}(undef, max(estimated_points, 10000))

    # RT bin tracking state — RT index is in empirical RT space
    rt_bin_start, rt_bin_stop = 1, 1
    prec_mz = zero(Float32)
    rt_idx = 0
    precs_temp = getPrecIds(search_data)
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    has_rt_tol = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol = has_rt_tol ? getRtTolerance(search_context, ms_file_idx) : nothing
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    nce_model = getNceModel(search_context, ms_file_idx)
    mass_error_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)
    prec_irt_arr = getIrt(precursors)
    frag_lookup = getFragmentLookupTable(spec_lib)

    kind = FusedRTIndexed(params.prec_estimation, UInt8(params.max_frag_rank))

    for scan_idx in scan_range
        ((scan_idx<1) | (scan_idx > length(spectra))) && continue
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && continue

        rt = getRetentionTime(spectra, scan_idx)

        if rt_binned_tol !== nothing
            rt_tol_local = get_rt_tol(rt_binned_tol, Float32(rt))
        else
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end

        rt_bin_start_new = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local, lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop_new = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        prec_mz_new = getCenterMz(spectra, scan_idx)
        if (rt_bin_start_new != rt_bin_start) || (rt_bin_stop_new != rt_bin_stop) || (prec_mz_new != prec_mz)
            rt_bin_start = rt_bin_start_new
            rt_bin_stop = rt_bin_stop_new
            prec_mz = prec_mz_new
        end

        quad_func = getQuadTransmissionFunction(
            quad_model,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx)
        )

        # 1. Enumerate precursors in the RT window (+ quad + per-prec RT +
        #    min_fraction_transmitted + allowlist filters). Writes to
        #    precs_temp[1:prec_temp_size].
        prec_temp_size = collect_rt_window_precursors!(
            precs_temp, rt_index, rt_bin_start, rt_bin_stop,
            precursors_passing, prec_mz_arr, prec_charge_arr, prec_sulfur_arr,
            getIsoSplines(search_data), quad_func, prec_trans_buf,
            params.isotope_err_bounds, params.min_fraction_transmitted,
            precursor_rt_map, Float32(rt), rt_binned_tol, rt_tol_local)

        if prec_temp_size == 0
            reset!(id_to_col); reset!(Hs)
            continue
        end

        # 2. Pre-correct peak m/z for the scan.
        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                          mass_error_model, scan_mz, scan_int)

        # 3. Fused match+build. iRT + iso_err_bounds filters already done
        #    upstream — skipped by FusedRTIndexed's check_prec_filters = false.
        nmatches, nmisses = run_fused!(
            kind,
            Hs, unscored_psms, id_to_col, fused_scratch,
            corr_mz, obs_low, obs_high, peak_mz_len,
            isotopes_buf, prec_trans_buf,
            frag_lookup, nce_model,
            precs_temp, 1:prec_temp_size,
            prec_mz_arr, prec_charge_arr, prec_sulfur_arr, prec_irt_arr,
            getIsoSplines(search_data), quad_func, mass_error_model,
            scan_int, 0f0, Float32(Inf),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx)),
            params.n_frag_isotopes,
            params.isotope_err_bounds)

        if nmatches > 2
            # Resize weight buffers for new columns.
            n_active_cols = n_active(id_to_col)
            if n_active_cols > length(weights)
                new_entries = n_active_cols - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
            end

            initialize_weights!(id_to_col, weights, precursor_weights)

            solve_deconvolution!(
                params.deconvolution_solver,
                Hs, residuals, weights, colnorm2,
                getMu(search_data), getObserved(search_data),
                params.max_iter_outer, params.max_diff)

            # Record chromatogram points (weighted for matches, zero otherwise).
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)
                end
                col = id_to_col[precs_temp[j]]
                chromatograms[rt_idx] = MS2ChromObject(
                    Float32(getRetentionTime(spectra, scan_idx)),
                    iszero(col) ? zero(Float32) : weights[col],
                    scan_idx,
                    precs_temp[j])
            end

            update_precursor_weights!(id_to_col, weights, precursor_weights)
        else
            # Too few matches → zero-weight points for every enumerated precursor.
            for j in 1:prec_temp_size
                rt_idx += 1
                if rt_idx + 1 > length(chromatograms)
                    resize!(chromatograms, length(chromatograms) * 2)
                end
                chromatograms[rt_idx] = MS2ChromObject(
                    Float32(getRetentionTime(spectra, scan_idx)),
                    zero(Float32), scan_idx, precs_temp[j])
            end
        end

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
end

# MS1 chromatogram extraction is currently unwired (the ms1_quant knob
# was deleted from the public schema). The body is preserved in a block-
# comment so it can be revived — but it uses classic ion_matches /
# buildDesignMatrix! / PrecursorMatch machinery that the fused cleanup
# PR will delete. If MS1 quant is re-enabled in the future, this needs
# a fused port (a FusedMs1 kind + an MS1-specific collector, or a
# parallel lightweight kernel that doesn't touch SparseArrayFused).
#=
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

    # RT bin tracking state — RT index is in empirical RT space
    rt_bin_start, rt_bin_stop = 1, 1
    ion_idx = 0
    rt_idx = 0
    precs_temp = getPrecIds(search_data)  # Use search_data's prec_ids
    prec_temp_size = 0
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    # RT-adaptive tolerance for MS1
    has_rt_tol_ms1 = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol_ms1 = has_rt_tol_ms1 ? getRtTolerance(search_context, ms_file_idx) : nothing
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
        # Calculate RT window — RT index is in empirical RT space
        rt = getRetentionTime(spectra, scan_idx)

        if rt_binned_tol_ms1 !== nothing
            # RT-adaptive: look up RT tolerance directly in RT space
            rt_tol_local = get_rt_tol(rt_binned_tol_ms1, Float32(rt))
        else
            # Fallback: convert global iRT tolerance to RT space via local slope
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end

        rt_bin_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local, lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Update transitions if window changed
        prec_temp_size = 0
        ion_idx = 0
        for rt_bin_idx in rt_bin_start:rt_bin_stop
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

            initialize_weights!(id_to_col, weights, precursor_weights)

            # Solve deconvolution
            solve_deconvolution!(
                params.deconvolution_solver,
                Hs, residuals, weights, colnorm2,
                getMu(search_data), getObserved(search_data),
                params.max_iter_outer, params.max_diff
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

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    return DataFrame(@view(chromatograms[1:rt_idx]))
end
=#

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
    ms_file_idxs = Vector{UInt32}(undef, n)
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

    # No sort here. MaxLFQSearch sorts the merged PSMs by :inferred_protein_group
    # before its chunked-merge so chunk boundaries align with protein-group
    # boundaries. ProteinInferenceSearch (which runs between this method and
    # MaxLFQ) is what populates :inferred_protein_group, so a sort by that
    # column at this stage isn't possible anyway.

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
