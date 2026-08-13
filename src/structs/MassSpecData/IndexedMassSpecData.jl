# View-based wrapper for memory-efficient scan selection.
#
# IndexedMassSpecData exposes a subset of scans from a source MassSpecData
# via index remapping. No data is copied — all getters delegate to the
# original data through the scan index mapping.
#
# Used by ParameterTuningSearch for progressive sampling without replacement.

"""
    IndexedMassSpecData{T<:MassSpecData}

A view-like wrapper that only exposes selected scan indices from the
underlying MassSpecData. Virtual index 1 maps to scan_indices[1], etc.
"""
struct IndexedMassSpecData{T<:MassSpecData} <: MassSpecData
    original_data::T
    scan_indices::Vector{Int32}
    n_scans::Int32

    function IndexedMassSpecData(original_data::T, scan_indices::Vector{Int32}) where T<:MassSpecData
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

Base.length(data::IndexedMassSpecData) = Int(data.n_scans)

@inline function get_actual_index(data::IndexedMassSpecData, virtual_idx::Integer)
    @boundscheck if virtual_idx < 1 || virtual_idx > data.n_scans
        throw(BoundsError(data, virtual_idx))
    end
    return data.scan_indices[virtual_idx]
end

# ============================================================================
# Singular Getters (delegate through index mapping)
# ============================================================================

getMzArray(data::IndexedMassSpecData, vi::Integer) = getMzArray(data.original_data, get_actual_index(data, vi))
getIntensityArray(data::IndexedMassSpecData, vi::Integer) = getIntensityArray(data.original_data, get_actual_index(data, vi))
getScanHeader(data::IndexedMassSpecData, vi::Integer) = getScanHeader(data.original_data, get_actual_index(data, vi))
getScanNumber(data::IndexedMassSpecData, vi::Integer) = getScanNumber(data.original_data, get_actual_index(data, vi))
getRetentionTime(data::IndexedMassSpecData, vi::Integer) = getRetentionTime(data.original_data, get_actual_index(data, vi))
getMsOrder(data::IndexedMassSpecData, vi::Integer) = getMsOrder(data.original_data, get_actual_index(data, vi))
getCycleIdx(data::IndexedMassSpecData, vi::Integer) = getCycleIdx(data.original_data, get_actual_index(data, vi))
getBasePeakMz(data::IndexedMassSpecData, vi::Integer) = getBasePeakMz(data.original_data, get_actual_index(data, vi))
getBasePeakIntensity(data::IndexedMassSpecData, vi::Integer) = getBasePeakIntensity(data.original_data, get_actual_index(data, vi))
getTIC(data::IndexedMassSpecData, vi::Integer) = getTIC(data.original_data, get_actual_index(data, vi))
getPrecursorMz(data::IndexedMassSpecData, vi::Integer) = getPrecursorMz(data.original_data, get_actual_index(data, vi))
getIsolationWidthMz(data::IndexedMassSpecData, vi::Integer) = getIsolationWidthMz(data.original_data, get_actual_index(data, vi))
getCenterMz(data::IndexedMassSpecData, vi::Integer) = getCenterMz(data.original_data, get_actual_index(data, vi))
getLowMz(data::IndexedMassSpecData, vi::Integer) = getLowMz(data.original_data, get_actual_index(data, vi))
getHighMz(data::IndexedMassSpecData, vi::Integer) = getHighMz(data.original_data, get_actual_index(data, vi))

# ============================================================================
# Plural Getters (construct arrays from mapped indices)
# ============================================================================

getMzArrays(data::IndexedMassSpecData) = [getMzArray(data.original_data, ai) for ai in data.scan_indices]
getIntensityArrays(data::IndexedMassSpecData) = [getIntensityArray(data.original_data, ai) for ai in data.scan_indices]
getRetentionTimes(data::IndexedMassSpecData) = [getRetentionTime(data.original_data, ai) for ai in data.scan_indices]
getMsOrders(data::IndexedMassSpecData) = [getMsOrder(data.original_data, ai) for ai in data.scan_indices]
getCycleIdxs(data::IndexedMassSpecData) = [getCycleIdx(data.original_data, ai) for ai in data.scan_indices]
getTICs(data::IndexedMassSpecData) = [getTIC(data.original_data, ai) for ai in data.scan_indices]
getCenterMzs(data::IndexedMassSpecData) = [getCenterMz(data.original_data, ai) for ai in data.scan_indices]
getIsolationWidthMzs(data::IndexedMassSpecData) = [getIsolationWidthMz(data.original_data, ai) for ai in data.scan_indices]

# ============================================================================
# Index Mapping Utilities
# ============================================================================

getOriginalScanIndex(data::IndexedMassSpecData, virtual_idx::Integer)::Int32 = data.scan_indices[virtual_idx]
getOriginalScanIndices(data::IndexedMassSpecData)::Vector{Int32} = data.scan_indices

"""
Create a Dict mapping virtual scan indices to actual scan indices.
"""
function create_scan_mapping(data::IndexedMassSpecData)
    mapping = Dict{Int32, Int32}()
    for (virtual_idx, actual_idx) in enumerate(data.scan_indices)
        mapping[Int32(virtual_idx)] = actual_idx
    end
    return mapping
end
