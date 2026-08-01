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

"""
    _MSDataCol{E,S}

The two representations a scalar `NonIonMobilityData` column can have: `Arrow.Primitive` for
single-batch files and `SentinelArrays.ChainedVector` for multi-batch ones. `E` is the element type,
`S` the underlying storage type (they differ only for nullable columns, where `E` is
`Union{Missing,S}`).

Annotating the plural getters with this Union is what makes the *scalar* getters allocation-free.
`.data` is an `Arrow.Table` — concrete, but schema-dynamic — so `ms_data.data[:col]` infers as
`AbstractVector`; indexing that is a dynamic dispatch whose result is boxed, which the scalar
assertions added earlier could only check, not prevent. A two-member Union is concrete enough for
Julia to union-split, so `getindex` becomes two statically-known branches and stops boxing, while a
single method still covers both file layouts.

Measured on a tight `getRetentionTime` loop: 31.9 -> 0.00 B/scan on a 132,210-scan multi-batch file
and a 19,275-scan single-batch file alike.

The constructor in `types.jl` already branches on exactly this two-way distinction
(`typeof(table[:msOrder]) == Arrow.Primitive{UInt8, Vector{UInt8}}`), so the assumption that these are
the only two representations is one the codebase already makes. If it is ever violated, this raises a
loud `TypeError` at the access site rather than silently degrading.
"""
const _MSDataCol{E,S} = Union{Arrow.Primitive{E, Vector{S}},
                              SentinelArrays.ChainedVector{E, Arrow.Primitive{E, Vector{S}}}}

# Deliberately NOT annotated: mz_array / intensity_array use the Arrow.List shape rather than
# Arrow.Primitive, and their per-access cost is SubArray construction (~60 B) rather than a scalar
# box, so they need separate investigation. Their scalar getters are already asserted above.
getMzArrays(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:mz_array]
getIntensityArrays(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:intensity_array]

getRetentionTimes(ms_data::NonIonMobilityData{T}) where T =
    ms_data.data[:retentionTime]::_MSDataCol{T,T}
getLowMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:lowMz]::_MSDataCol{T,T}
getHighMzs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:highMz]::_MSDataCol{T,T}
getTICs(ms_data::NonIonMobilityData{T}) where T = ms_data.data[:TIC]::_MSDataCol{T,T}
# centerMz / isolationWidthMz are nullable in the Arrow schema, so the element type carries Missing
# while the storage type does not.
getCenterMzs(ms_data::NonIonMobilityData{T}) where T =
    ms_data.data[:centerMz]::_MSDataCol{Union{Missing,T},T}
getIsolationWidthMzs(ms_data::NonIonMobilityData{T}) where T =
    ms_data.data[:isolationWidthMz]::_MSDataCol{Union{Missing,T},T}
getMsOrders(ms_data::NonIonMobilityData{T}) where T =
    ms_data.data[:msOrder]::_MSDataCol{UInt8,UInt8}
getCycleIdxs(ms_data::NonIonMobilityData{T}) where T =
    ms_data.cycle_idxs === nothing ? ms_data.data[:cycle_idx] : ms_data.cycle_idxs
