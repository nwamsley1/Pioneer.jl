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

using Random

"""
In-memory filtered and sampled view of mass spectrometry data.
Uses intelligent scan selection based on RT binning and peak density.
"""
mutable struct FilteredMassSpecData{T<:AbstractFloat} <: MassSpecData
    # Peak data arrays - store filtered/sampled data
    mz_arrays::Vector{Vector{T}}
    intensity_arrays::Vector{Vector{T}}
    
    # Metadata arrays - one entry per sampled scan
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
    center_mzs::Vector{T}
    TICs::Vector{T}
    
    # Additional scan metadata
    low_mzs::Vector{T}
    high_mzs::Vector{T}
    isolation_width_mzs::Vector{T}
    
    # Index mapping - critical for correct PSM reporting
    original_scan_indices::Vector{UInt32}  # Maps filtered index -> original scan index
    
    # Reference to source data
    original_data::MassSpecData
    
    # Filtering configuration
    topn::Union{Nothing, Int}           # Nothing = no peak filtering
    min_intensity::T                    # Minimum intensity threshold
    max_scans::Int                      # Maximum number of scans to sample
    target_ms_order::Union{Nothing, UInt8}  # Which MS order to include
    rng::MersenneTwister               # RNG for reproducible sampling
    
    # Intelligent scan selection
    scan_priority_order::Vector{Int32}    # Pre-computed scan ordering by RT bin and peak density
    n_scans_sampled::Int32                # How many scans from priority order are included
    rt_bin_assignments::Vector{Int8}      # Which RT bin each scan belongs to
    n_rt_bins::Int                        # Number of RT bins used
    
    # Sampling state
    sampled_scans::Set{Int}            # Track which scans have been sampled
    total_ms2_scans::Int                # Total MS2 scans in original data
    
    # Current size
    n::Int                             # Number of sampled scans
end

# ============================================================================
# Intelligent Scan Selection Functions
# ============================================================================

"""
Compute RT bins for all scans in the data.
Returns bin assignments, RT range, and bin width.
"""
function compute_rt_bins(spectra::MassSpecData, n_bins::Int = 15)
    rt_values = getRetentionTimes(spectra)
    
    # Handle empty case
    if isempty(rt_values)
        return Int8[], 0.0, 0.0, 0.0
    end
    
    rt_min, rt_max = extrema(rt_values)
    
    # Handle edge case of single RT value
    if rt_min == rt_max
        return fill(Int8(1), length(rt_values)), rt_min, rt_max, 0.0
    end
    
    bin_width = (rt_max - rt_min) / n_bins
    
    # Assign each scan to a bin
    rt_bin_assignments = Vector{Int8}(undef, length(rt_values))
    for (idx, rt) in enumerate(rt_values)
        # Ensure bin_idx is 1-based (ceil(0.0) returns 0, so we need max(1, ...))
        bin_idx = min(max(1, ceil(Int, (rt - rt_min) / bin_width)), n_bins)
        rt_bin_assignments[idx] = Int8(bin_idx)
    end
    
    return rt_bin_assignments, rt_min, rt_max, bin_width
end

"""
Sort scans within RT bins by peak density (number of peaks).
Returns sorted array and bin boundaries.
"""
function sort_scans_by_peak_density(
    spectra::MassSpecData,
    target_ms_order::Union{Nothing, UInt8},
    rt_bin_assignments::Vector{Int8},
    n_bins::Int = 15
)
    # Get MS orders and filter for target
    ms_orders = getMsOrders(spectra)
    
    # Collect eligible scan indices
    if target_ms_order !== nothing
        target_scan_indices = [Int32(i) for i in 1:length(spectra) 
                               if ms_orders[i] == target_ms_order]
    else
        target_scan_indices = Int32.(1:length(spectra))
    end
    
    n_target_scans = length(target_scan_indices)
    
    # Handle empty case
    if n_target_scans == 0
        return Int32[], Int32[], Int32[]
    end
    
    # Pre-allocate single array for all scans
    all_sorted_scans = Vector{Int32}(undef, n_target_scans)
    
    # Count scans per bin
    bin_counts = zeros(Int32, n_bins)
    for scan_idx in target_scan_indices
        bin_idx = rt_bin_assignments[scan_idx]
        bin_counts[bin_idx] += 1
    end
    
    # Calculate bin boundaries
    bin_starts = Vector{Int32}(undef, n_bins)
    bin_ends = Vector{Int32}(undef, n_bins)
    cumsum = 0
    for i in 1:n_bins
        bin_starts[i] = cumsum + 1
        cumsum += bin_counts[i]
        bin_ends[i] = cumsum
    end
    
    # Reset bin_counts to use as write positions
    fill!(bin_counts, 0)
    
    # Write scans to their bin positions
    for scan_idx in target_scan_indices
        bin_idx = rt_bin_assignments[scan_idx]
        write_pos = bin_starts[bin_idx] + bin_counts[bin_idx]
        all_sorted_scans[write_pos] = scan_idx
        bin_counts[bin_idx] += 1
    end
    
    # Get peak counts (use TIC as proxy for peak density)
    tics = getTICs(spectra)
    
    # In-place sort each bin by peak count (descending)
    for i in 1:n_bins
        if bin_starts[i] <= bin_ends[i]
            # Sort this bin's slice in-place by TIC (descending)
            bin_slice = view(all_sorted_scans, bin_starts[i]:bin_ends[i])
            sort!(bin_slice, by=idx -> tics[idx], rev=true)
        end
    end
    
    return all_sorted_scans, bin_starts, bin_ends
end

"""
Create interleaved priority order from sorted bins.
Takes one scan from each bin in round-robin fashion.
"""
function create_priority_order(
    all_sorted_scans::Vector{Int32},
    bin_starts::Vector{Int32},
    bin_ends::Vector{Int32}
)
    n_bins = length(bin_starts)
    n_total_scans = length(all_sorted_scans)
    
    # Handle empty case
    if n_total_scans == 0
        return Int32[]
    end
    
    # Pre-allocate priority order array
    priority_order = Vector{Int32}(undef, n_total_scans)
    
    # Track current position in each bin
    bin_positions = copy(bin_starts)
    
    # Interleave: take one scan from each bin in round-robin fashion
    write_idx = 1
    scans_remaining = n_total_scans
    
    while scans_remaining > 0
        for bin_idx in 1:n_bins
            if bin_positions[bin_idx] <= bin_ends[bin_idx]
                # Write next scan from this bin
                priority_order[write_idx] = all_sorted_scans[bin_positions[bin_idx]]
                bin_positions[bin_idx] += 1
                write_idx += 1
                scans_remaining -= 1
                
                if scans_remaining == 0
                    break
                end
            end
        end
    end
    
    return priority_order
end

"""
    FilteredMassSpecData(original; max_scans, topn, kwargs...)

Create a filtered and sampled view of mass spectrometry data using intelligent scan selection.

# Arguments
- `original::MassSpecData`: Original data to filter/sample

# Keyword Arguments
- `max_scans::Int`: Maximum number of scans to sample (default: 2500)
- `topn::Union{Nothing, Int}`: Number of top peaks to keep per scan (default: nothing = no filtering)
- `min_intensity`: Minimum intensity threshold for peaks
- `target_ms_order`: MS order to filter for (default: UInt8(2) for MS2)
- `n_rt_bins::Int`: Number of RT bins for intelligent selection (default: 15)
- `seed`: Random seed for fallback sampling if needed
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
    # Determine float type from original data
    T = Float32  # Default, but could inspect original data type
    
    # Initialize RNG for fallback
    rng = seed === nothing ? MersenneTwister() : MersenneTwister(seed)
    
    # Convert min_intensity to correct type
    min_intensity_typed = min_intensity === nothing ? zero(T) : T(min_intensity)
    
    # Phase 1: Intelligent scan selection
    # Step 1: Compute RT bins
    rt_bin_assignments, rt_min, rt_max, bin_width = compute_rt_bins(original, n_rt_bins)
    
    # Step 2: Sort scans within bins by peak density
    all_sorted_scans, bin_starts, bin_ends = sort_scans_by_peak_density(
        original, target_ms_order, rt_bin_assignments, n_rt_bins
    )
    
    # Step 3: Create interleaved priority order
    scan_priority_order = create_priority_order(all_sorted_scans, bin_starts, bin_ends)
    
    total_ms2_scans = length(scan_priority_order)
    
    # Step 4: Sample initial scans following priority order
    n_to_sample = min(max_scans, total_ms2_scans)
    if n_to_sample == 0
        # Handle empty case - no scans to sample
        scan_indices_to_sample = UInt32[]
    elseif n_to_sample == total_ms2_scans
        scan_indices_to_sample = convert(Vector{UInt32}, scan_priority_order)
    else
        # Take first n_to_sample from priority order
        scan_indices_to_sample = UInt32.(scan_priority_order[1:n_to_sample])
    end
    
    n_sampled = length(scan_indices_to_sample)
    
    # Phase 2: Pre-allocate all arrays
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
    center_mzs = Vector{T}(undef, n_sampled)
    TICs = Vector{T}(undef, n_sampled)
    low_mzs = Vector{T}(undef, n_sampled)
    high_mzs = Vector{T}(undef, n_sampled)
    isolation_width_mzs = Vector{T}(undef, n_sampled)
    
    # Pre-allocate working buffer for topN filtering
    indices_buffer = topn !== nothing ? Vector{Int}(undef, 10000) : Int[]
    
    # Phase 3: Process each sampled scan
    for (filtered_idx, original_idx) in enumerate(scan_indices_to_sample)
        # Get peak data
        mz_array = getMzArray(original, original_idx)
        intensity_array = getIntensityArray(original, original_idx)
        
        # Apply filtering if requested
        if topn !== nothing && length(mz_array) > topn
            mz_filtered, intensity_filtered = filterTopNPeaks(
                mz_array, 
                intensity_array, 
                topn, 
                indices_buffer,
                min_intensity_typed
            )
            mz_arrays[filtered_idx] = mz_filtered
            intensity_arrays[filtered_idx] = intensity_filtered
        else
            # No filtering - convert Arrow arrays to Julia arrays
            mz_arrays[filtered_idx] = [T(x) for x in mz_array if !ismissing(x)]
            intensity_arrays[filtered_idx] = [T(x) for x in intensity_array if !ismissing(x)]
        end
        
        # Copy all metadata
        scan_headers[filtered_idx] = getScanHeader(original, original_idx)
        scan_numbers[filtered_idx] = getScanNumber(original, original_idx)
        base_peak_mzs[filtered_idx] = T(getBasePeakMz(original, original_idx))
        base_peak_intensities[filtered_idx] = T(getBasePeakIntensity(original, original_idx))
        injection_times[filtered_idx] = zero(T)  # Default to zero - injection time not available in BasicNonIonMobilityMassSpecData
        retention_times[filtered_idx] = T(getRetentionTime(original, original_idx))
        precursor_mzs[filtered_idx] = T(getPrecursorMz(original, original_idx))
        isolation_widths[filtered_idx] = T(getIsolationWidthMz(original, original_idx))
        precursor_charges[filtered_idx] = Int8(0)  # Default charge - not available in scan-level data
        ms_orders[filtered_idx] = getMsOrder(original, original_idx)
        center_mzs[filtered_idx] = T(getCenterMz(original, original_idx))
        TICs[filtered_idx] = T(getTIC(original, original_idx))
        low_mzs[filtered_idx] = T(getLowMz(original, original_idx))
        high_mzs[filtered_idx] = T(getHighMz(original, original_idx))
        isolation_width_mzs[filtered_idx] = T(getIsolationWidthMz(original, original_idx))
    end
    
    # Create the struct
    return FilteredMassSpecData{T}(
        mz_arrays, intensity_arrays,
        scan_headers, scan_numbers,
        base_peak_mzs, base_peak_intensities,
        injection_times, retention_times,
        precursor_mzs, isolation_widths,
        precursor_charges, ms_orders,
        center_mzs, TICs,
        low_mzs, high_mzs, isolation_width_mzs,
        scan_indices_to_sample,  # original_scan_indices
        original,                # reference
        topn, min_intensity_typed, max_scans, target_ms_order, rng,
        scan_priority_order,     # intelligent scan order
        Int32(n_sampled),        # n_scans_sampled
        rt_bin_assignments,      # RT bin assignments
        n_rt_bins,              # number of RT bins
        Set(scan_indices_to_sample),  # sampled_scans
        total_ms2_scans,             # total_ms2_scans
        n_sampled                     # n
    )
end

"""
Filter to keep only top N peaks by intensity, maintaining m/z order.
Returns filtered arrays, not views, for memory locality.
"""
function filterTopNPeaks(
    mz_array::AbstractArray,
    intensity_array::AbstractArray,
    topn::Int,
    indices_buffer::Vector{Int},
    min_intensity::T
) where {T<:AbstractFloat}
    n_peaks = length(mz_array)
    
    # Count valid peaks (non-missing, above threshold)
    valid_count = 0
    for i in 1:n_peaks
        if !ismissing(intensity_array[i]) && intensity_array[i] >= min_intensity
            valid_count += 1
            indices_buffer[valid_count] = i
        end
    end
    
    # Determine how many to keep
    n_to_keep = min(topn, valid_count)
    
    # If keeping all valid peaks, no need to sort by intensity
    if n_to_keep == valid_count
        mz_filtered = Vector{T}(undef, n_to_keep)
        intensity_filtered = Vector{T}(undef, n_to_keep)
        
        for i in 1:n_to_keep
            idx = indices_buffer[i]
            mz_filtered[i] = T(mz_array[idx])
            intensity_filtered[i] = T(intensity_array[idx])
        end
        
        # Sort by m/z
        perm = sortperm(mz_filtered)
        return mz_filtered[perm], intensity_filtered[perm]
    end
    
    # Need to select top N by intensity
    # Partial sort to find top N
    partialsortperm!(
        view(indices_buffer, 1:valid_count),
        [intensity_array[indices_buffer[i]] for i in 1:valid_count],
        1:n_to_keep,
        rev=true
    )
    
    # Extract top N
    mz_filtered = Vector{T}(undef, n_to_keep)
    intensity_filtered = Vector{T}(undef, n_to_keep)
    
    for i in 1:n_to_keep
        idx = indices_buffer[i]
        mz_filtered[i] = T(mz_array[idx])
        intensity_filtered[i] = T(intensity_array[idx])
    end
    
    # Sort by m/z for correct peak matching
    perm = sortperm(mz_filtered)
    return mz_filtered[perm], intensity_filtered[perm]
end

# Implement MassSpecData interface
Base.length(ms_data::FilteredMassSpecData) = ms_data.n

# Peak data access - return as arrays that can contain Missing for compatibility
getMzArray(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = convert(Vector{Union{Missing, T}}, ms_data.mz_arrays[scan_idx])
getIntensityArray(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = convert(Vector{Union{Missing, T}}, ms_data.intensity_arrays[scan_idx])

# Metadata access - implement ALL methods from MassSpecData interface
getScanHeader(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.scan_headers[scan_idx]
getScanNumber(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.scan_numbers[scan_idx]
getBasePeakMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.base_peak_mzs[scan_idx]
getBasePeakIntensity(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.base_peak_intensities[scan_idx]
getInjectionTime(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.injection_times[scan_idx]
getRetentionTime(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.retention_times[scan_idx]
getPrecursorMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.precursor_mzs[scan_idx]
getIsolationWidthMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.isolation_width_mzs[scan_idx]
getPrecursorCharge(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.precursor_charges[scan_idx]
getMsOrder(ms_data::FilteredMassSpecData, scan_idx::Integer) = ms_data.ms_orders[scan_idx]
getCenterMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.center_mzs[scan_idx]
getTIC(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.TICs[scan_idx]
getLowMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.low_mzs[scan_idx]
getHighMz(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.high_mzs[scan_idx]
getIsolationWidth(ms_data::FilteredMassSpecData{T}, scan_idx::Integer) where T = ms_data.isolation_widths[scan_idx]

# Batch getters - return the arrays
getMzArrays(ms_data::FilteredMassSpecData{T}) where T = ms_data.mz_arrays
getIntensityArrays(ms_data::FilteredMassSpecData{T}) where T = ms_data.intensity_arrays
getScanHeaders(ms_data::FilteredMassSpecData) = ms_data.scan_headers
getScanNumbers(ms_data::FilteredMassSpecData) = ms_data.scan_numbers
getBasePeakMzs(ms_data::FilteredMassSpecData{T}) where T = ms_data.base_peak_mzs
getBasePeakIntensities(ms_data::FilteredMassSpecData{T}) where T = ms_data.base_peak_intensities
getInjectionTimes(ms_data::FilteredMassSpecData{T}) where T = ms_data.injection_times
getRetentionTimes(ms_data::FilteredMassSpecData{T}) where T = ms_data.retention_times
getPrecursorMzs(ms_data::FilteredMassSpecData{T}) where T = ms_data.precursor_mzs
getIsolationWidthMzs(ms_data::FilteredMassSpecData{T}) where T = ms_data.isolation_width_mzs
getPrecursorCharges(ms_data::FilteredMassSpecData) = ms_data.precursor_charges
getMsOrders(ms_data::FilteredMassSpecData) = ms_data.ms_orders
getCenterMzs(ms_data::FilteredMassSpecData{T}) where T = convert(Vector{Union{Missing, T}}, ms_data.center_mzs)
getTICs(ms_data::FilteredMassSpecData{T}) where T = ms_data.TICs
getLowMzs(ms_data::FilteredMassSpecData{T}) where T = ms_data.low_mzs
getHighMzs(ms_data::FilteredMassSpecData{T}) where T = ms_data.high_mzs

# Special methods for index mapping
getOriginalScanIndex(ms_data::FilteredMassSpecData, filtered_idx::Integer)::UInt32 = ms_data.original_scan_indices[filtered_idx]
getOriginalScanIndices(ms_data::FilteredMassSpecData)::Vector{UInt32} = ms_data.original_scan_indices

# Helper to check if more scans available for sampling
hasUnsampledScans(ms_data::FilteredMassSpecData) = length(ms_data.sampled_scans) < ms_data.total_ms2_scans

# Helper to get count of unsampled scans
getUnsampledCount(ms_data::FilteredMassSpecData) = ms_data.total_ms2_scans - length(ms_data.sampled_scans)

"""
Append additional sampled scans to the filtered data using intelligent priority order.
Returns the number of scans added.
"""
function Base.append!(
    filtered::FilteredMassSpecData{T};
    max_additional_scans::Int = 2500
) where {T<:AbstractFloat}
    
    original = filtered.original_data
    
    # Determine how many more scans we can add from priority order
    n_available = length(filtered.scan_priority_order) - filtered.n_scans_sampled
    n_to_add = min(max_additional_scans, n_available)
    
    # Return early if nothing to add
    if n_to_add <= 0
        return 0
    end
    
    # Get next scans from priority order
    start_idx = filtered.n_scans_sampled + 1
    end_idx = filtered.n_scans_sampled + n_to_add
    new_scan_indices = UInt32.(filtered.scan_priority_order[start_idx:end_idx])
    
    # Update sampling counter
    filtered.n_scans_sampled += Int32(n_to_add)
    
    # Add to sampled set
    for idx in new_scan_indices
        push!(filtered.sampled_scans, idx)
    end
    
    # Return early if no new scans
    isempty(new_scan_indices) && return 0
    
    # Phase 2: Process new scans
    indices_buffer = filtered.topn !== nothing ? Vector{Int}(undef, 10000) : Int[]
    
    for original_idx in new_scan_indices
        # Get peak data
        mz_array = getMzArray(original, original_idx)
        intensity_array = getIntensityArray(original, original_idx)
        
        # Apply filtering if configured
        if filtered.topn !== nothing && length(mz_array) > filtered.topn
            mz_filtered, intensity_filtered = filterTopNPeaks(
                mz_array, 
                intensity_array, 
                filtered.topn, 
                indices_buffer,
                filtered.min_intensity
            )
            push!(filtered.mz_arrays, mz_filtered)
            push!(filtered.intensity_arrays, intensity_filtered)
        else
            # No filtering - convert to Julia arrays
            push!(filtered.mz_arrays, [T(x) for x in mz_array if !ismissing(x)])
            push!(filtered.intensity_arrays, [T(x) for x in intensity_array if !ismissing(x)])
        end
        
        # Append all metadata
        push!(filtered.scan_headers, getScanHeader(original, original_idx))
        push!(filtered.scan_numbers, getScanNumber(original, original_idx))
        push!(filtered.base_peak_mzs, T(getBasePeakMz(original, original_idx)))
        push!(filtered.base_peak_intensities, T(getBasePeakIntensity(original, original_idx)))
        push!(filtered.injection_times, zero(T))  # Default to zero - injection time not available in BasicNonIonMobilityMassSpecData
        push!(filtered.retention_times, T(getRetentionTime(original, original_idx)))
        push!(filtered.precursor_mzs, T(getPrecursorMz(original, original_idx)))
        push!(filtered.isolation_widths, T(getIsolationWidthMz(original, original_idx)))
        push!(filtered.precursor_charges, Int8(0))  # Default charge - not available in scan-level data
        push!(filtered.ms_orders, getMsOrder(original, original_idx))
        push!(filtered.center_mzs, T(getCenterMz(original, original_idx)))
        push!(filtered.TICs, T(getTIC(original, original_idx)))
        push!(filtered.low_mzs, T(getLowMz(original, original_idx)))
        push!(filtered.high_mzs, T(getHighMz(original, original_idx)))
        push!(filtered.isolation_width_mzs, T(getIsolationWidthMz(original, original_idx)))
        
        # Update index mapping
        push!(filtered.original_scan_indices, original_idx)
    end
    
    # Update count
    filtered.n += length(new_scan_indices)

    return length(new_scan_indices)
end

#==========================================================
IndexedMassSpecData - View-based wrapper for memory-efficient scan selection
==========================================================#

"""
    IndexedMassSpecData{T<:MassSpecData}

A view-like wrapper around MassSpecData that only exposes selected scan indices.
Acts as a transparent proxy for library_search - all getter methods redirect
through the scan index mapping to only show the selected scans.

# Purpose
- Allows library_search to process only a subset of scans without modification
- No data duplication - only stores index mapping
- Dynamic scan selection by swapping scan_indices vector
- Perfect for progressive sampling without replacement

# Example Usage
```julia
# Create view of only specific scans
target_scans = [100, 250, 890, 1200]  # The scans we want to process
indexed_spectra = IndexedMassSpecData(original_spectra, target_scans)

# Library search now only "sees" these 4 scans as indices 1, 2, 3, 4
psms = library_search(indexed_spectra, search_context, params, ms_file_idx)
# No post-filtering needed - only target scans were processed
```
"""
struct IndexedMassSpecData{T<:MassSpecData} <: MassSpecData
    # The underlying data source (never modified)
    original_data::T

    # The scan indices we want to expose (1-based indices into original_data)
    scan_indices::Vector{Int32}

    # Cached length for performance
    n_scans::Int32

    function IndexedMassSpecData(original_data::T, scan_indices::Vector{Int32}) where T<:MassSpecData
        # Validate scan indices
        if !isempty(scan_indices)
            max_idx = maximum(scan_indices)
            min_idx = minimum(scan_indices)
            if max_idx > length(original_data) || min_idx < 1
                throw(BoundsError("Scan indices must be between 1 and $(length(original_data))"))
            end
        end

        new{T}(original_data, scan_indices, Int32(length(scan_indices)))
    end
end

# =============================================================================
# Core Interface Methods
# =============================================================================

"""
Return the number of scans visible through this indexed view.
"""
Base.length(data::IndexedMassSpecData) = Int(data.n_scans)

"""
Map virtual index (1-based index into the view) to actual index in original data.
"""
@inline function get_actual_index(data::IndexedMassSpecData, virtual_idx::Integer)
    @boundscheck if virtual_idx < 1 || virtual_idx > data.n_scans
        throw(BoundsError(data, virtual_idx))
    end
    return data.scan_indices[virtual_idx]
end

# =============================================================================
# Peak Data Access Methods
# =============================================================================

function getMzArray(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getMzArray(data.original_data, actual_idx)
end

function getIntensityArray(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIntensityArray(data.original_data, actual_idx)
end

# =============================================================================
# Scan Metadata Access Methods
# =============================================================================

function getScanHeader(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getScanHeader(data.original_data, actual_idx)
end

function getScanNumber(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getScanNumber(data.original_data, actual_idx)
end

function getRetentionTime(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getRetentionTime(data.original_data, actual_idx)
end

function getMsOrder(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getMsOrder(data.original_data, actual_idx)
end

function getBasePeakMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getBasePeakMz(data.original_data, actual_idx)
end

function getBasePeakIntensity(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getBasePeakIntensity(data.original_data, actual_idx)
end

function getTIC(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getTIC(data.original_data, actual_idx)
end

function getInjectionTime(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getInjectionTime(data.original_data, actual_idx)
end

function getPrecursorMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getPrecursorMz(data.original_data, actual_idx)
end

function getPrecursorCharge(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getPrecursorCharge(data.original_data, actual_idx)
end

function getIsolationWidthMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIsolationWidthMz(data.original_data, actual_idx)
end

function getIsolationWidth(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIsolationWidth(data.original_data, actual_idx)
end

function getCenterMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getCenterMz(data.original_data, actual_idx)
end

function getLowMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getLowMz(data.original_data, actual_idx)
end

function getHighMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getHighMz(data.original_data, actual_idx)
end

# =============================================================================
# Batch Getter Methods (for performance)
# =============================================================================

function getMzArrays(data::IndexedMassSpecData)
    return [getMzArray(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getIntensityArrays(data::IndexedMassSpecData)
    return [getIntensityArray(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getScanHeaders(data::IndexedMassSpecData)
    return [getScanHeader(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getScanNumbers(data::IndexedMassSpecData)
    return [getScanNumber(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getRetentionTimes(data::IndexedMassSpecData)
    return [getRetentionTime(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getMsOrders(data::IndexedMassSpecData)
    return [getMsOrder(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getBasePeakMzs(data::IndexedMassSpecData)
    return [getBasePeakMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getBasePeakIntensities(data::IndexedMassSpecData)
    return [getBasePeakIntensity(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getTICs(data::IndexedMassSpecData)
    return [getTIC(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getInjectionTimes(data::IndexedMassSpecData)
    return [getInjectionTime(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getPrecursorMzs(data::IndexedMassSpecData)
    return [getPrecursorMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getPrecursorCharges(data::IndexedMassSpecData)
    return [getPrecursorCharge(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getIsolationWidthMzs(data::IndexedMassSpecData)
    return [getIsolationWidthMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getCenterMzs(data::IndexedMassSpecData)
    return [getCenterMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getLowMzs(data::IndexedMassSpecData)
    return [getLowMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

function getHighMzs(data::IndexedMassSpecData)
    return [getHighMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

# =============================================================================
# Utility Methods
# =============================================================================

"""
    create_scan_mapping(data::IndexedMassSpecData)

Create a mapping from virtual scan indices (1-based in view) to actual scan indices.
Returns a dictionary for easy lookup during PSM processing.
"""
function create_scan_mapping(data::IndexedMassSpecData)
    mapping = Dict{Int32, Int32}()
    for (virtual_idx, actual_idx) in enumerate(data.scan_indices)
        mapping[Int32(virtual_idx)] = actual_idx
    end
    return mapping
end

# Index mapping methods for IndexedMassSpecData
getOriginalScanIndex(data::IndexedMassSpecData, virtual_idx::Integer)::Int32 = data.scan_indices[virtual_idx]
getOriginalScanIndices(data::IndexedMassSpecData)::Vector{Int32} = data.scan_indices