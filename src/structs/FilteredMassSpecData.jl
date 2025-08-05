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
Samples a fraction of scans and stores only the top N peaks by intensity for each.
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
    sample_rate::Float64                # Fraction of scans to sample
    target_ms_order::Union{Nothing, UInt8}  # Which MS order to include
    rng::MersenneTwister               # RNG for reproducible sampling
    
    # Sampling state
    sampled_scans::Set{Int}            # Track which scans have been sampled
    
    # Current size
    n::Int                             # Number of sampled scans
end

"""
    FilteredMassSpecData(original, sample_rate, topn; kwargs...)

Create a filtered and sampled view of mass spectrometry data.

# Arguments
- `original::MassSpecData`: Original data to filter/sample
- `sample_rate::Float64`: Fraction of scans to sample (0.0 to 1.0)
- `topn::Union{Nothing, Int}`: Number of top peaks to keep per scan (nothing = no filtering)

# Keyword Arguments
- `min_intensity`: Minimum intensity threshold for peaks
- `target_ms_order`: MS order to filter for (default: UInt8(2) for MS2)
- `seed`: Random seed for reproducible sampling
"""
function FilteredMassSpecData(
    original::MassSpecData,
    sample_rate::AbstractFloat,
    topn::Union{Nothing, Int} = nothing;
    min_intensity::Union{Nothing, AbstractFloat} = nothing,
    target_ms_order::Union{Nothing, UInt8} = UInt8(2),
    seed::Union{Nothing, Int} = nothing
)
    # Determine float type from original data
    T = Float32  # Default, but could inspect original data type
    
    # Initialize RNG
    rng = seed === nothing ? MersenneTwister() : MersenneTwister(seed)
    
    # Convert min_intensity to correct type
    min_intensity_typed = min_intensity === nothing ? zero(T) : T(min_intensity)
    
    # Phase 1: Determine which scans to sample
    scan_indices_to_sample = UInt32[]
    for scan_idx in 1:length(original)
        # Apply MS order filter
        if target_ms_order !== nothing && getMsOrder(original, scan_idx) != target_ms_order
            continue
        end
        
        # Apply sampling
        if rand(rng) <= sample_rate
            push!(scan_indices_to_sample, UInt32(scan_idx))
        end
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
        topn, min_intensity_typed, sample_rate, target_ms_order, rng,
        Set(scan_indices_to_sample),  # sampled_scans
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
hasUnsampledScans(ms_data::FilteredMassSpecData) = length(ms_data.sampled_scans) < length(ms_data.original_data)

"""
Append additional sampled scans to the filtered data.
Returns the number of scans added.
"""
function Base.append!(
    filtered::FilteredMassSpecData{T},
    additional_sample_rate::Float64;
    max_additional_scans::Union{Nothing, Int} = nothing
) where {T<:AbstractFloat}
    
    original = filtered.original_data
    new_scan_indices = UInt32[]
    
    # Phase 1: Select additional scans to sample
    for scan_idx in 1:length(original)
        # Skip if already sampled
        scan_idx in filtered.sampled_scans && continue
        
        # Apply MS order filter if it was used originally
        if filtered.target_ms_order !== nothing && getMsOrder(original, scan_idx) != filtered.target_ms_order
            continue
        end
        
        # Sample based on rate
        if rand(filtered.rng) <= additional_sample_rate
            push!(new_scan_indices, UInt32(scan_idx))
            push!(filtered.sampled_scans, scan_idx)
            
            # Check limit
            if max_additional_scans !== nothing && length(new_scan_indices) >= max_additional_scans
                break
            end
        end
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