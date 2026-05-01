# Getter methods for NonIonMobilityData (Basic and Batch).
#
# Each getter dispatches on NonIonMobilityData — Julia's type inference
# handles the Arrow.Primitive vs SentinelArrays.ChainedVector difference
# automatically without explicit return-type annotations.
#
# Only getters with external callers are defined here. Internal-only getters
# (scanHeader, scanNumber, basePeak*, packetType, injectionTime, precursorCharge)
# were removed — they were only used by FilteredMassSpecData's constructor to
# copy metadata, which now accesses ms_data.data columns directly.

#==========================================================
Singular getters (per-scan access)
==========================================================#

getMzArray(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getMzArrays(ms_data)[scan_idx]
getIntensityArray(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getIntensityArrays(ms_data)[scan_idx]
getRetentionTime(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getRetentionTimes(ms_data)[scan_idx]
getLowMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getLowMzs(ms_data)[scan_idx]
getHighMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getHighMzs(ms_data)[scan_idx]
getTIC(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getTICs(ms_data)[scan_idx]
getCenterMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getCenterMzs(ms_data)[scan_idx]
getIsolationWidthMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getIsolationWidthMzs(ms_data)[scan_idx]
getMsOrder(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getMsOrders(ms_data)[scan_idx]

# Aliases used by BuildSpecLib
getPrecursorMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getCenterMz(ms_data, scan_idx)

# Used by FilteredMassSpecData constructor
getScanHeader(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = ms_data.data[:scanHeader][scan_idx]
getScanNumber(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = ms_data.data[:scanNumber][scan_idx]
getBasePeakMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = ms_data.data[:basePeakMz][scan_idx]
getBasePeakIntensity(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = ms_data.data[:basePeakIntensity][scan_idx]

#==========================================================
Plural/batch getters (full-column access)
==========================================================#

getMzArrays(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:mz_array]
getIntensityArrays(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:intensity_array]
getRetentionTimes(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:retentionTime]
getLowMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:lowMz]
getHighMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:highMz]
getTICs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:TIC]
getCenterMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:centerMz]
getIsolationWidthMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:isolationWidthMz]
getMsOrders(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:msOrder]
