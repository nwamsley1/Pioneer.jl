# Getter methods for NonIonMobilityData (Basic and Batch).
#
# Each getter dispatches on NonIonMobilityData. The column *representation* differs between the two
# (Arrow.Primitive for single-batch files, SentinelArrays.ChainedVector for multi-batch), so the
# plural getters below cannot be annotated with a single column type — but the scalar getters can and
# must be annotated with their element type. Without that, `.data` being an `Arrow.Table` (concrete
# but schema-dynamic) made every per-scan access infer as `Any` and box. Measured on a 132,210-scan
# KEAP1 file: `getRetentionTime` in a tight loop went 47.9 -> 31.9 B/scan by adding the assertions.
#
# Only getters with external callers are defined here. Internal-only getters
# (scanHeader, scanNumber, basePeak*, packetType, injectionTime, precursorCharge)
# were removed — they were only used by FilteredMassSpecData's constructor to
# copy metadata, which now accesses ms_data.data columns directly.

#==========================================================
Singular getters (per-scan access)
==========================================================#

getMzArray(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getMzArrays(ms_data)[scan_idx]::SubArray{Union{Missing,T},1,Arrow.Primitive{Union{Missing,T},Vector{T}},Tuple{UnitRange{Int64}},true}
getIntensityArray(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getIntensityArrays(ms_data)[scan_idx]::SubArray{Union{Missing,T},1,Arrow.Primitive{Union{Missing,T},Vector{T}},Tuple{UnitRange{Int64}},true}
# Scalar return types are asserted. `.data` is an `Arrow.Table`, which is concrete but
# schema-dynamic: `ms_data.data[:col]` infers as `AbstractVector`, so an unannotated
# `getXs(ms_data)[scan_idx]` returned `Any` and **boxed on every per-scan access**. A KEAP1 run
# showed `getRetentionTime` alone at ~43M boxed 8-byte allocations, and it is the single largest
# allocation site in `build_chromatograms` (which calls it per chromatogram row).
#
# Asserting the *scalar* type rather than the column type is what makes this safe: the column
# representation genuinely varies (`Arrow.Primitive` for single-batch files vs
# `SentinelArrays.ChainedVector` for multi-batch), which is what the header note above warns about,
# but the element type does not.
#
# Caveat: `T` is derived from the mz_array element type in `BasicMassSpecData`, not per column, so a
# file whose retentionTime/lowMz/highMz/TIC width differs from its mz_array width would raise a
# TypeError here. That is a loud, immediate failure rather than silent corruption, and it mirrors
# what `getMzArray`/`getIntensityArray` above have always asserted.
getRetentionTime(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getRetentionTimes(ms_data)[scan_idx]::T
getLowMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getLowMzs(ms_data)[scan_idx]::T
getHighMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getHighMzs(ms_data)[scan_idx]::T
getTIC(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getTICs(ms_data)[scan_idx]::T
# centerMz / isolationWidthMz are genuinely nullable in the Arrow schema; the two-member Union is
# still concrete enough for Julia to union-split instead of boxing.
getCenterMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getCenterMzs(ms_data)[scan_idx]::Union{Missing,T}
getIsolationWidthMz(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getIsolationWidthMzs(ms_data)[scan_idx]::Union{Missing,T}
getMsOrder(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = getMsOrders(ms_data)[scan_idx]::UInt8
# The `UInt32(...)` conversion is not enough on its own: its argument is `Any`, so the call is a
# dynamic dispatch and the result still inferred as `Any`. The assertion pins it.
getCycleIdx(ms_data::NonIonMobilityData{T}, scan_idx::Integer) where T = UInt32(getCycleIdxs(ms_data)[scan_idx])::UInt32

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
getCycleIdxs(ms_data::NonIonMobilityData{T}) where T =
    ms_data.cycle_idxs === nothing ? ms_data.data[:cycle_idx] : ms_data.cycle_idxs
