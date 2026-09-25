# In-memory filtered and sampled view of mass spectrometry data.
#
# FilteredMassSpecData copies a subset of scans from a source MassSpecData
# into Julia arrays, optionally applying top-N peak filtering. Scan selection
# uses RT-binned priority ordering to ensure uniform gradient coverage.
#
# Used by ParameterTuningSearch and NceTuningSearch to work on a small
# representative subset of scans for fast calibration.

using Random

const SCAN_PRIORITY_N_RT_BINS = 30

mutable struct FilteredMassSpecData{T<:AbstractFloat} <: MassSpecData
    # Peak data arrays
    mz_arrays::Vector{Vector{T}}
    intensity_arrays::Vector{Vector{T}}

    # Metadata arrays (one per sampled scan)
    scan_headers::Vector{String}
    scan_numbers::Vector{Int32}
    base_peak_mzs::Vector{T}
    base_peak_intensities::Vector{T}
    injection_times::Vector{T}
    retention_times::Vector{T}
    precursor_mzs::Vector{T}
    isolation_widths::Vector{T}
    precursor_charges::Vector{Int8}
    ms_orders::Vector{UInt8}
    cycle_idxs::Vector{UInt32}
    center_mzs::Vector{T}
    TICs::Vector{T}
    low_mzs::Vector{T}
    high_mzs::Vector{T}
    isolation_width_mzs::Vector{T}

    # Index mapping: filtered index → original scan index
    original_scan_indices::Vector{UInt32}

    # Reference to source data
    original_data::MassSpecData

    # Filtering configuration
    topn::Union{Nothing, Int}
    min_intensity::T
    max_scans::Int
    target_ms_order::Union{Nothing, UInt8}
    rng::MersenneTwister

    # Scan selection state
    scan_priority_order::Vector{Int32}
    n_scans_sampled::Int32
    rt_bin_assignments::Vector{Int8}
    n_rt_bins::Int
    sampled_scans::Set{Int}
    total_ms2_scans::Int

    n::Int
end

# ============================================================================
# Scan Selection Helpers
# ============================================================================

"""
Compute RT bins for all scans. Returns (bin_assignments, rt_min, rt_max, bin_width).
"""
function compute_rt_bins(spectra::MassSpecData, n_bins::Int = 15)
    rt_values = getRetentionTimes(spectra)
    isempty(rt_values) && return Int8[], 0.0, 0.0, 0.0

    rt_min, rt_max = extrema(rt_values)
    rt_min == rt_max && return fill(Int8(1), length(rt_values)), rt_min, rt_max, 0.0

    bin_width = (rt_max - rt_min) / n_bins
    rt_bin_assignments = Vector{Int8}(undef, length(rt_values))
    for (idx, rt) in enumerate(rt_values)
        rt_bin_assignments[idx] = Int8(min(max(1, ceil(Int, (rt - rt_min) / bin_width)), n_bins))
    end
    return rt_bin_assignments, rt_min, rt_max, bin_width
end

"""
Sort scans within RT bins by TIC (descending). Returns (sorted_scans, bin_starts, bin_ends).
"""
function sort_scans_by_peak_density(
    spectra::MassSpecData,
    target_ms_order::Union{Nothing, UInt8},
    rt_bin_assignments::Vector{Int8},
    n_bins::Int = 15
)
    ms_orders = getMsOrders(spectra)
    target_scan_indices = if target_ms_order !== nothing
        [Int32(i) for i in 1:length(spectra) if ms_orders[i] == target_ms_order]
    else
        Int32.(1:length(spectra))
    end

    n_target_scans = length(target_scan_indices)
    n_target_scans == 0 && return Int32[], Int32[], Int32[]

    all_sorted_scans = Vector{Int32}(undef, n_target_scans)

    # Count and compute bin boundaries
    bin_counts = zeros(Int32, n_bins)
    for scan_idx in target_scan_indices
        bin_counts[rt_bin_assignments[scan_idx]] += 1
    end
    bin_starts = Vector{Int32}(undef, n_bins)
    bin_ends = Vector{Int32}(undef, n_bins)
    cumsum = 0
    for i in 1:n_bins
        bin_starts[i] = cumsum + 1
        cumsum += bin_counts[i]
        bin_ends[i] = cumsum
    end

    # Distribute scans into bins
    fill!(bin_counts, 0)
    for scan_idx in target_scan_indices
        bin_idx = rt_bin_assignments[scan_idx]
        all_sorted_scans[bin_starts[bin_idx] + bin_counts[bin_idx]] = scan_idx
        bin_counts[bin_idx] += 1
    end

    # Sort each bin by TIC descending
    tics = getTICs(spectra)
    for i in 1:n_bins
        if bin_starts[i] <= bin_ends[i]
            sort!(view(all_sorted_scans, bin_starts[i]:bin_ends[i]), by=idx -> tics[idx], rev=true)
        end
    end

    return all_sorted_scans, bin_starts, bin_ends
end

"""
Create interleaved priority order: round-robin one scan from each bin.
"""
function create_priority_order(
    all_sorted_scans::Vector{Int32},
    bin_starts::Vector{Int32},
    bin_ends::Vector{Int32}
)
    n_bins = length(bin_starts)
    n_total_scans = length(all_sorted_scans)
    n_total_scans == 0 && return Int32[]

    priority_order = Vector{Int32}(undef, n_total_scans)
    bin_positions = copy(bin_starts)
    write_idx = 1
    scans_remaining = n_total_scans

    while scans_remaining > 0
        for bin_idx in 1:n_bins
            if bin_positions[bin_idx] <= bin_ends[bin_idx]
                priority_order[write_idx] = all_sorted_scans[bin_positions[bin_idx]]
                bin_positions[bin_idx] += 1
                write_idx += 1
                scans_remaining -= 1
                scans_remaining == 0 && break
            end
        end
    end
    return priority_order
end

"""
    get_ms2_scan_priority_order(spectra::MassSpecData) -> Vector{Int32}

Return MS2 scan indices ordered by priority: round-robin across RT bins,
sorted by TIC within each bin (highest TIC first). This is the same ordering
used by FilteredMassSpecData for scan selection, but without allocating copies.
"""
function get_ms2_scan_priority_order(spectra::MassSpecData)
    rt_bin_assignments, _, _, _ = compute_rt_bins(spectra, SCAN_PRIORITY_N_RT_BINS)
    all_sorted_scans, bin_starts, bin_ends = sort_scans_by_peak_density(
        spectra, UInt8(2), rt_bin_assignments, SCAN_PRIORITY_N_RT_BINS)
    return create_priority_order(all_sorted_scans, bin_starts, bin_ends)
end

# ============================================================================
# Constructor
# ============================================================================

"""
    FilteredMassSpecData(original; max_scans=2500, topn=nothing, kwargs...)

Create a filtered/sampled view using RT-binned priority scan selection.
"""
function FilteredMassSpecData(
    original::MassSpecData;
    max_scans::Int = 2500,
    topn::Union{Nothing, Int} = nothing,
    min_intensity::Union{Nothing, AbstractFloat} = nothing,
    target_ms_order::Union{Nothing, UInt8} = UInt8(2),
    n_rt_bins::Int = 15,
    seed::Union{Nothing, Int} = nothing
)
    T = Float32
    rng = seed === nothing ? MersenneTwister() : MersenneTwister(seed)
    min_intensity_typed = min_intensity === nothing ? zero(T) : T(min_intensity)

    # Scan selection
    rt_bin_assignments, _, _, _ = compute_rt_bins(original, n_rt_bins)
    all_sorted_scans, bin_starts, bin_ends = sort_scans_by_peak_density(original, target_ms_order, rt_bin_assignments, n_rt_bins)
    scan_priority_order = create_priority_order(all_sorted_scans, bin_starts, bin_ends)
    total_ms2_scans = length(scan_priority_order)

    n_to_sample = min(max_scans, total_ms2_scans)
    scan_indices_to_sample = if n_to_sample == 0
        UInt32[]
    elseif n_to_sample == total_ms2_scans
        convert(Vector{UInt32}, scan_priority_order)
    else
        UInt32.(scan_priority_order[1:n_to_sample])
    end
    n_sampled = length(scan_indices_to_sample)

    # Pre-allocate arrays
    mz_arrays = Vector{Vector{T}}(undef, n_sampled)
    intensity_arrays = Vector{Vector{T}}(undef, n_sampled)
    scan_headers = Vector{String}(undef, n_sampled)
    scan_numbers = Vector{Int32}(undef, n_sampled)
    base_peak_mzs = Vector{T}(undef, n_sampled)
    base_peak_intensities = Vector{T}(undef, n_sampled)
    injection_times = Vector{T}(undef, n_sampled)
    retention_times = Vector{T}(undef, n_sampled)
    precursor_mzs = Vector{T}(undef, n_sampled)
    isolation_widths = Vector{T}(undef, n_sampled)
    precursor_charges = Vector{Int8}(undef, n_sampled)
    ms_orders = Vector{UInt8}(undef, n_sampled)
    cycle_idxs = Vector{UInt32}(undef, n_sampled)
    center_mzs = Vector{T}(undef, n_sampled)
    TICs = Vector{T}(undef, n_sampled)
    low_mzs = Vector{T}(undef, n_sampled)
    high_mzs = Vector{T}(undef, n_sampled)
    isolation_width_mzs = Vector{T}(undef, n_sampled)

    indices_buffer = topn !== nothing ? Vector{Int}(undef, 10000) : Int[]

    # Copy scan data
    for (i, oi) in enumerate(scan_indices_to_sample)
        mz_array = getMzArray(original, oi)
        int_array = getIntensityArray(original, oi)

        if topn !== nothing && length(mz_array) > topn
            mz_arrays[i], intensity_arrays[i] = filterTopNPeaks(mz_array, int_array, topn, indices_buffer, min_intensity_typed)
        else
            mz_arrays[i] = [T(x) for x in mz_array if !ismissing(x)]
            intensity_arrays[i] = [T(x) for x in int_array if !ismissing(x)]
        end

        scan_headers[i] = getScanHeader(original, oi)
        scan_numbers[i] = getScanNumber(original, oi)
        base_peak_mzs[i] = T(getBasePeakMz(original, oi))
        base_peak_intensities[i] = T(getBasePeakIntensity(original, oi))
        injection_times[i] = zero(T)
        retention_times[i] = T(getRetentionTime(original, oi))
        precursor_mzs[i] = T(coalesce(getCenterMz(original, oi), zero(T)))
        isolation_widths[i] = T(coalesce(getIsolationWidthMz(original, oi), zero(T)))
        precursor_charges[i] = Int8(0)
        ms_orders[i] = getMsOrder(original, oi)
        cycle_idxs[i] = getCycleIdx(original, oi)
        center_mzs[i] = T(coalesce(getCenterMz(original, oi), zero(T)))
        TICs[i] = T(getTIC(original, oi))
        low_mzs[i] = T(getLowMz(original, oi))
        high_mzs[i] = T(getHighMz(original, oi))
        isolation_width_mzs[i] = T(coalesce(getIsolationWidthMz(original, oi), zero(T)))
    end

    return FilteredMassSpecData{T}(
        mz_arrays, intensity_arrays,
        scan_headers, scan_numbers, base_peak_mzs, base_peak_intensities,
        injection_times, retention_times, precursor_mzs, isolation_widths,
        precursor_charges, ms_orders, cycle_idxs, center_mzs, TICs,
        low_mzs, high_mzs, isolation_width_mzs,
        scan_indices_to_sample, original, topn, min_intensity_typed,
        max_scans, target_ms_order, rng, scan_priority_order,
        Int32(n_sampled), rt_bin_assignments, n_rt_bins,
        Set(scan_indices_to_sample), total_ms2_scans, n_sampled
    )
end

"""
Filter to keep only top N peaks by intensity, maintaining m/z order.
"""
function filterTopNPeaks(
    mz_array::AbstractArray,
    intensity_array::AbstractArray,
    topn::Int,
    indices_buffer::Vector{Int},
    min_intensity::T
) where {T<:AbstractFloat}
    n_peaks = length(mz_array)

    valid_count = 0
    for i in 1:n_peaks
        if !ismissing(intensity_array[i]) && intensity_array[i] >= min_intensity
            valid_count += 1
            if valid_count > length(indices_buffer)
                resize!(indices_buffer, max(valid_count, 2 * length(indices_buffer)))
            end
            indices_buffer[valid_count] = i
        end
    end

    n_to_keep = min(topn, valid_count)

    if n_to_keep == valid_count
        mz_filtered = Vector{T}(undef, n_to_keep)
        intensity_filtered = Vector{T}(undef, n_to_keep)
        for i in 1:n_to_keep
            idx = indices_buffer[i]
            mz_filtered[i] = T(mz_array[idx])
            intensity_filtered[i] = T(intensity_array[idx])
        end
        perm = sortperm(mz_filtered)
        return mz_filtered[perm], intensity_filtered[perm]
    end

    partialsortperm!(
        view(indices_buffer, 1:valid_count),
        [intensity_array[indices_buffer[i]] for i in 1:valid_count],
        1:n_to_keep, rev=true
    )

    mz_filtered = Vector{T}(undef, n_to_keep)
    intensity_filtered = Vector{T}(undef, n_to_keep)
    for i in 1:n_to_keep
        idx = indices_buffer[i]
        mz_filtered[i] = T(mz_array[idx])
        intensity_filtered[i] = T(intensity_array[idx])
    end
    perm = sortperm(mz_filtered)
    return mz_filtered[perm], intensity_filtered[perm]
end

# ============================================================================
# MassSpecData Interface — Singular Getters
# ============================================================================

Base.length(ms_data::FilteredMassSpecData) = ms_data.n
getMzArray(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = convert(Vector{Union{Missing, T}}, ms_data.mz_arrays[scan_idx])
getIntensityArray(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = convert(Vector{Union{Missing, T}}, ms_data.intensity_arrays[scan_idx])
getScanHeader(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.scan_headers[scan_idx]
getScanNumber(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.scan_numbers[scan_idx]
getBasePeakMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.base_peak_mzs[scan_idx]
getBasePeakIntensity(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.base_peak_intensities[scan_idx]
getRetentionTime(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.retention_times[scan_idx]
getPrecursorMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.precursor_mzs[scan_idx]
getIsolationWidthMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.isolation_width_mzs[scan_idx]
getMsOrder(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.ms_orders[scan_idx]
getCycleIdx(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.cycle_idxs[scan_idx]
getCenterMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.center_mzs[scan_idx]
getTIC(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.TICs[scan_idx]
getLowMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.low_mzs[scan_idx]
getHighMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.high_mzs[scan_idx]

# ============================================================================
# MassSpecData Interface — Plural/Batch Getters
# ============================================================================

getMzArrays(ms_data::FilteredMassSpecData) = ms_data.mz_arrays
getIntensityArrays(ms_data::FilteredMassSpecData) = ms_data.intensity_arrays
getRetentionTimes(ms_data::FilteredMassSpecData) = ms_data.retention_times
getTICs(ms_data::FilteredMassSpecData) = ms_data.TICs
getCenterMzs(ms_data::FilteredMassSpecData{T}) where T = convert(Vector{Union{Missing, T}}, ms_data.center_mzs)
getIsolationWidthMzs(ms_data::FilteredMassSpecData) = ms_data.isolation_width_mzs
getMsOrders(ms_data::FilteredMassSpecData) = ms_data.ms_orders
getCycleIdxs(ms_data::FilteredMassSpecData) = ms_data.cycle_idxs

# ============================================================================
# Index Mapping
# ============================================================================

getOriginalScanIndex(ms_data::FilteredMassSpecData, filtered_idx::Integer)::UInt32 = ms_data.original_scan_indices[filtered_idx]
getOriginalScanIndices(ms_data::FilteredMassSpecData)::Vector{UInt32} = ms_data.original_scan_indices

# ============================================================================
# Append (progressive sampling)
# ============================================================================

"""
Append additional scans from the priority order. Returns count of scans added.
"""
function Base.append!(
    filtered::FilteredMassSpecData{T};
    max_additional_scans::Int = 2500
) where {T<:AbstractFloat}
    original = filtered.original_data
    n_available = length(filtered.scan_priority_order) - filtered.n_scans_sampled
    n_to_add = min(max_additional_scans, n_available)
    n_to_add <= 0 && return 0

    start_idx = filtered.n_scans_sampled + 1
    end_idx = filtered.n_scans_sampled + n_to_add
    new_scan_indices = UInt32.(filtered.scan_priority_order[start_idx:end_idx])
    filtered.n_scans_sampled += Int32(n_to_add)
    for idx in new_scan_indices
        push!(filtered.sampled_scans, idx)
    end

    indices_buffer = filtered.topn !== nothing ? Vector{Int}(undef, 10000) : Int[]

    for oi in new_scan_indices
        mz_array = getMzArray(original, oi)
        int_array = getIntensityArray(original, oi)

        if filtered.topn !== nothing && length(mz_array) > filtered.topn
            mz_f, int_f = filterTopNPeaks(mz_array, int_array, filtered.topn, indices_buffer, filtered.min_intensity)
            push!(filtered.mz_arrays, mz_f)
            push!(filtered.intensity_arrays, int_f)
        else
            push!(filtered.mz_arrays, [T(x) for x in mz_array if !ismissing(x)])
            push!(filtered.intensity_arrays, [T(x) for x in int_array if !ismissing(x)])
        end

        push!(filtered.scan_headers, getScanHeader(original, oi))
        push!(filtered.scan_numbers, getScanNumber(original, oi))
        push!(filtered.base_peak_mzs, T(getBasePeakMz(original, oi)))
        push!(filtered.base_peak_intensities, T(getBasePeakIntensity(original, oi)))
        push!(filtered.injection_times, zero(T))
        push!(filtered.retention_times, T(getRetentionTime(original, oi)))
        push!(filtered.precursor_mzs, T(coalesce(getCenterMz(original, oi), zero(T))))
        push!(filtered.isolation_widths, T(coalesce(getIsolationWidthMz(original, oi), zero(T))))
        push!(filtered.precursor_charges, Int8(0))
        push!(filtered.ms_orders, getMsOrder(original, oi))
        push!(filtered.cycle_idxs, getCycleIdx(original, oi))
        push!(filtered.center_mzs, T(coalesce(getCenterMz(original, oi), zero(T))))
        push!(filtered.TICs, T(getTIC(original, oi)))
        push!(filtered.low_mzs, T(getLowMz(original, oi)))
        push!(filtered.high_mzs, T(getHighMz(original, oi)))
        push!(filtered.isolation_width_mzs, T(coalesce(getIsolationWidthMz(original, oi), zero(T))))
        push!(filtered.original_scan_indices, oi)
    end

    filtered.n += length(new_scan_indices)
    return length(new_scan_indices)
end
