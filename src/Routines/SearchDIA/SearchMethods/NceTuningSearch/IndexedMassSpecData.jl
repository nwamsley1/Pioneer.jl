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

# Interface Transparency
The wrapper implements the complete MassSpecData interface:
- `length(indexed_spectra)` returns `length(scan_indices)`
- `getMzArray(indexed_spectra, 1)` returns data from `original_spectra[scan_indices[1]]`
- All other getters work similarly - virtual index mapped to actual index
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

"""
Get m/z array for scan at virtual index. Maps to actual scan in original data.
"""
function getMzArray(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getMzArray(data.original_data, actual_idx)
end

"""
Get intensity array for scan at virtual index. Maps to actual scan in original data.
"""
function getIntensityArray(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIntensityArray(data.original_data, actual_idx)
end

# =============================================================================
# Scan Metadata Access Methods
# =============================================================================

"""
Get scan header for virtual index.
"""
function getScanHeader(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getScanHeader(data.original_data, actual_idx)
end

"""
Get scan number for virtual index.
"""
function getScanNumber(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getScanNumber(data.original_data, actual_idx)
end

"""
Get retention time for virtual index.
"""
function getRetentionTime(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getRetentionTime(data.original_data, actual_idx)
end

"""
Get MS order for virtual index.
"""
function getMsOrder(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getMsOrder(data.original_data, actual_idx)
end

"""
Get base peak m/z for virtual index.
"""
function getBasePeakMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getBasePeakMz(data.original_data, actual_idx)
end

"""
Get base peak intensity for virtual index.
"""
function getBasePeakIntensity(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getBasePeakIntensity(data.original_data, actual_idx)
end

"""
Get total ion current (TIC) for virtual index.
"""
function getTIC(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getTIC(data.original_data, actual_idx)
end

"""
Get injection time for virtual index.
"""
function getInjectionTime(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getInjectionTime(data.original_data, actual_idx)
end

# =============================================================================
# Precursor Information Methods
# =============================================================================

"""
Get precursor m/z for virtual index.
"""
function getPrecursorMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getPrecursorMz(data.original_data, actual_idx)
end

"""
Get precursor charge for virtual index.
"""
function getPrecursorCharge(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getPrecursorCharge(data.original_data, actual_idx)
end

"""
Get isolation width in m/z for virtual index.
"""
function getIsolationWidthMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIsolationWidthMz(data.original_data, actual_idx)
end

"""
Get isolation width for virtual index.
"""
function getIsolationWidth(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getIsolationWidth(data.original_data, actual_idx)
end

"""
Get center m/z for virtual index.
"""
function getCenterMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getCenterMz(data.original_data, actual_idx)
end

"""
Get low m/z for virtual index.
"""
function getLowMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getLowMz(data.original_data, actual_idx)
end

"""
Get high m/z for virtual index.
"""
function getHighMz(data::IndexedMassSpecData, virtual_idx::Integer)
    actual_idx = get_actual_index(data, virtual_idx)
    return getHighMz(data.original_data, actual_idx)
end

# =============================================================================
# Batch Getter Methods (for performance)
# =============================================================================

"""
Get all m/z arrays visible through this view.
"""
function getMzArrays(data::IndexedMassSpecData)
    return [getMzArray(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all intensity arrays visible through this view.
"""
function getIntensityArrays(data::IndexedMassSpecData)
    return [getIntensityArray(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all scan headers visible through this view.
"""
function getScanHeaders(data::IndexedMassSpecData)
    return [getScanHeader(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all scan numbers visible through this view.
"""
function getScanNumbers(data::IndexedMassSpecData)
    return [getScanNumber(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all retention times visible through this view.
"""
function getRetentionTimes(data::IndexedMassSpecData)
    return [getRetentionTime(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all MS orders visible through this view.
"""
function getMsOrders(data::IndexedMassSpecData)
    return [getMsOrder(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all base peak m/z values visible through this view.
"""
function getBasePeakMzs(data::IndexedMassSpecData)
    return [getBasePeakMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all base peak intensities visible through this view.
"""
function getBasePeakIntensities(data::IndexedMassSpecData)
    return [getBasePeakIntensity(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all TIC values visible through this view.
"""
function getTICs(data::IndexedMassSpecData)
    return [getTIC(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all injection times visible through this view.
"""
function getInjectionTimes(data::IndexedMassSpecData)
    return [getInjectionTime(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all precursor m/z values visible through this view.
"""
function getPrecursorMzs(data::IndexedMassSpecData)
    return [getPrecursorMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all precursor charges visible through this view.
"""
function getPrecursorCharges(data::IndexedMassSpecData)
    return [getPrecursorCharge(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all isolation widths in m/z visible through this view.
"""
function getIsolationWidthMzs(data::IndexedMassSpecData)
    return [getIsolationWidthMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all center m/z values visible through this view.
"""
function getCenterMzs(data::IndexedMassSpecData)
    return [getCenterMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all low m/z values visible through this view.
"""
function getLowMzs(data::IndexedMassSpecData)
    return [getLowMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

"""
Get all high m/z values visible through this view.
"""
function getHighMzs(data::IndexedMassSpecData)
    return [getHighMz(data.original_data, actual_idx) for actual_idx in data.scan_indices]
end

# =============================================================================
# Utility Methods
# =============================================================================

"""
    extend_indexed_view!(data::IndexedMassSpecData, additional_indices::Vector{Int32})

Extend the indexed view to include additional scan indices.
Useful for progressive sampling without replacement.

# Example
```julia
# Start with 5% of scans
indexed_data = IndexedMassSpecData(spectra, first_5_percent_indices)

# Add next 10% of scans
extend_indexed_view!(indexed_data, next_10_percent_indices)
```
"""
function extend_indexed_view!(data::IndexedMassSpecData, additional_indices::Vector{Int32})
    # Validate new indices
    if !isempty(additional_indices)
        max_idx = maximum(additional_indices)
        min_idx = minimum(additional_indices)
        if max_idx > length(data.original_data) || min_idx < 1
            throw(BoundsError("Additional scan indices must be between 1 and $(length(data.original_data))"))
        end
    end

    # Extend the scan indices
    append!(data.scan_indices, additional_indices)
    data.n_scans = Int32(length(data.scan_indices))

    return data
end

"""
    get_actual_scan_indices(data::IndexedMassSpecData)

Get the actual scan indices that this view exposes.
Useful for debugging or creating scan index mappings in PSM results.
"""
function get_actual_scan_indices(data::IndexedMassSpecData)
    return copy(data.scan_indices)
end

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