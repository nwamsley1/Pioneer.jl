# TdfsMassSpecData: Pioneer's view of a `.tdfs` (TimsSlices.jl) — IM-smoothed, m/z-centroided timsTOF slices
# stored as one zstd block per slice. Every slice is one Pioneer scan.
#
# Peak arrays are decoded on demand into a buffer the CALLER owns: `getPeaks!(buf, d, scan)` decodes the scan's
# block into `buf` (~5 µs for an MS2 slice, 30-50 µs for an MS1 slice) and returns views into it, valid until
# `buf` is used for a different scan. Each task owns its buffer (the search tasks carry one in their
# SearchDataStructures), so a task that yields or moves thread can never have its peaks overwritten by another
# task. Buffers are keyed on (file, scan), so one buffer can be reused across files.
# `getMzArray` / `getIntensityArray` are not defined for this type on purpose: a shared implicit buffer is exactly
# what this design avoids.
#
# The vectors are `Vector{Union{Missing,Float32}}` so the views match the `AbstractArray{Union{Missing,Float32}}`
# signatures of run_fused! / prepare_scan_peaks! / the fragment-index scorer unchanged (element types are
# invariant, so a Float32 view would not). No element is ever `missing`.
#
# m/z and intensity are the same Float32 values `TimsSlices.expand` writes into the slice Arrow
# (`Float32((a + b·bin/k)²)`, `Float32(stored / int_scale)`), so a search from the .tdfs and from its expanded
# Arrow see identical peaks.

using TimsSlices: TimsSlices, TdfsFile, SliceBuffer, BlockCodec

"""
    PeakDecodeBuffer()

Decode scratch for one task: `getPeaks!(buf, spectra, scan)` decodes `.tdfs` peaks into it. Arrow-backed data
ignores it. Grows to the largest slice decoded into it. Not thread-safe: one task uses one buffer.
"""
mutable struct PeakDecodeBuffer
    file_uid::Int                               # TdfsMassSpecData.uid of the decoded scan (0 = none)
    scan::Int                                   # scan currently decoded (0 = none)
    sb::SliceBuffer
    codec::BlockCodec
    mz::Vector{Union{Missing, Float32}}
    intensity::Vector{Union{Missing, Float32}}
end
PeakDecodeBuffer() = PeakDecodeBuffer(0, 0, SliceBuffer(), BlockCodec(), Union{Missing, Float32}[], Union{Missing, Float32}[])

# Identifies an opened file for the buffer's (file, scan) cache key.
const _TDFS_NEXT_UID = Threads.Atomic{Int}(1)

struct TdfsMassSpecData <: MassSpecData
    file::TdfsFile
    uid::Int                                    # unique per opened file (PeakDecodeBuffer cache key)
    n::Int
    retention_time::Vector{Float32}             # minutes
    low_mz::Vector{Float32}
    high_mz::Vector{Float32}
    tic::Vector{Float32}
    center_mz::Vector{Union{Missing, Float32}}
    isolation_width::Vector{Union{Missing, Float32}}
    collision_energy_ev::Vector{Float32}
    ms_order::Vector{UInt8}
    cycle_idx::Vector{UInt32}
    im_scan::Vector{UInt16}
    frame_id::Vector{Int32}
    n_peaks::Vector{Int32}
end

"""
    TdfsMassSpecData(dir::String)

Open a `<name>.tdfs` directory. The slice side table is loaded into plain columns; the block file is memory-mapped.
Peaks are decoded per scan with `getPeaks!`.
"""
function TdfsMassSpecData(dir::String)
    file = TimsSlices.open_tdfs(dir)
    sl = file.slices
    n = length(sl.frame_row)
    nan_to_missing(v) = Union{Missing, Float32}[isnan(x) ? missing : Float32(x) for x in v]
    mz_lo = Float32(file.meta["mz_lo"]); mz_hi = Float32(file.meta["mz_hi"])
    TdfsMassSpecData(
        file, Threads.atomic_add!(_TDFS_NEXT_UID, 1), n,
        Vector{Float32}(sl.retention_time), fill(mz_lo, n), fill(mz_hi, n), Vector{Float32}(sl.tic),
        nan_to_missing(sl.center_mz), nan_to_missing(sl.isolation_width),
        Vector{Float32}(sl.collision_energy_ev), Vector{UInt8}(sl.ms_order), Vector{UInt32}(sl.cycle_idx),
        Vector{UInt16}(sl.im_scan), Vector{Int32}(sl.frame_id), Vector{Int32}(sl.n_peaks),
    )
end

is_tdfs_path(path::AbstractString) = endswith(path, ".tdfs") && isdir(path)

Base.length(d::TdfsMassSpecData) = d.n

# ---- peak arrays: decode on demand into the caller's buffer ------------------------------------------------

function _tdfs_decode!(buf::PeakDecodeBuffer, d::TdfsMassSpecData, scan::Int)
    np = Int(@inbounds d.n_peaks[scan])
    TimsSlices.read_slice!(buf.sb, buf.codec, d.file, scan)
    if length(buf.mz) < np
        resize!(buf.mz, np); resize!(buf.intensity, np)
    end
    f = d.file; bins = buf.sb.bin; ints = buf.sb.intensity; mz = buf.mz; it = buf.intensity
    @inbounds for i in 1:np
        mz[i] = Float32(TimsSlices.bin_to_mz(f, bins[i]))
        it[i] = Float32(TimsSlices.stored_to_intensity(f, ints[i]))
    end
    buf.file_uid = d.uid
    buf.scan = scan
    buf
end

const _TdfsPeakView = SubArray{Union{Missing, Float32}, 1, Vector{Union{Missing, Float32}}, Tuple{UnitRange{Int64}}, true}

"""
    getPeaks!(buf::PeakDecodeBuffer, spectra::MassSpecData, scan_idx) -> (mz, intensity)

The scan's peak arrays. For `.tdfs` data they are decoded into `buf` (skipped when `buf` already holds this scan
of this file) and are views into it, valid until `buf` is used for another scan. Other data types return their
usual arrays and do not touch `buf`.
"""
function getPeaks!(buf::PeakDecodeBuffer, d::TdfsMassSpecData, scan_idx::Integer)
    scan = Int(scan_idx)
    1 <= scan <= d.n || throw(BoundsError(d, scan))
    (buf.scan == scan && buf.file_uid == d.uid) || _tdfs_decode!(buf, d, scan)
    np = Int(@inbounds d.n_peaks[scan])
    return (view(buf.mz, 1:np)::_TdfsPeakView, view(buf.intensity, 1:np)::_TdfsPeakView)
end
getPeaks!(::PeakDecodeBuffer, d::MassSpecData, scan_idx::Integer) =
    (getMzArray(d, scan_idx), getIntensityArray(d, scan_idx))
getPeaks!(buf::PeakDecodeBuffer, d::IndexedMassSpecData, vi::Integer) =
    getPeaks!(buf, d.original_data, get_actual_index(d, vi))

getMzArray(::TdfsMassSpecData, ::Integer) = error("TdfsMassSpecData decodes into a caller-owned buffer; use getPeaks!(buf, spectra, scan)")
getIntensityArray(::TdfsMassSpecData, ::Integer) = error("TdfsMassSpecData decodes into a caller-owned buffer; use getPeaks!(buf, spectra, scan)")
getMzArrays(::TdfsMassSpecData) = error("TdfsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getPeaks! per scan")
getIntensityArrays(::TdfsMassSpecData) = error("TdfsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getPeaks! per scan")

# ---- per-scan peak counts (no decode) -------------------------------------------------------------------------

getPeakCount(d::TdfsMassSpecData, scan_idx::Integer) = Int(d.n_peaks[scan_idx])
getPeakCounts(d::TdfsMassSpecData) = d.n_peaks
# Arrow-backed types: the list column's per-element length (offset differences; no peak data touched).
getPeakCount(ms_data::NonIonMobilityData, scan_idx::Integer) = length(getMzArray(ms_data, scan_idx))
getPeakCounts(ms_data::NonIonMobilityData) = Int32[length(x) for x in getMzArrays(ms_data)]
getPeakCount(data::IndexedMassSpecData, vi::Integer) = getPeakCount(data.original_data, get_actual_index(data, vi))
getPeakCounts(data::IndexedMassSpecData) = Int32[getPeakCount(data.original_data, i) for i in data.scan_indices]
getPeakCount(ms_data::FilteredMassSpecData, scan_idx::Integer) = length(ms_data.mz_arrays[scan_idx])
getPeakCounts(ms_data::FilteredMassSpecData) = Int32[length(x) for x in ms_data.mz_arrays]

# ---- scalar getters ---------------------------------------------------------------------------------------------

getRetentionTime(d::TdfsMassSpecData, i::Integer) = d.retention_time[i]::Float32
getLowMz(d::TdfsMassSpecData, i::Integer) = d.low_mz[i]::Float32
getHighMz(d::TdfsMassSpecData, i::Integer) = d.high_mz[i]::Float32
getTIC(d::TdfsMassSpecData, i::Integer) = d.tic[i]::Float32
getCenterMz(d::TdfsMassSpecData, i::Integer) = d.center_mz[i]::Union{Missing, Float32}
getIsolationWidthMz(d::TdfsMassSpecData, i::Integer) = d.isolation_width[i]::Union{Missing, Float32}
getMsOrder(d::TdfsMassSpecData, i::Integer) = d.ms_order[i]::UInt8
getCycleIdx(d::TdfsMassSpecData, i::Integer) = d.cycle_idx[i]::UInt32
getPrecursorMz(d::TdfsMassSpecData, i::Integer) = getCenterMz(d, i)
getCollisionEnergyEv(d::TdfsMassSpecData, i::Integer)::Float32 = d.collision_energy_ev[i]
getScanHeader(::TdfsMassSpecData, ::Integer) = ""
getScanNumber(::TdfsMassSpecData, i::Integer) = Int32(i)
getBasePeakMz(::TdfsMassSpecData, ::Integer) = missing
getBasePeakIntensity(::TdfsMassSpecData, ::Integer) = missing

# ---- plural getters ---------------------------------------------------------------------------------------------

getRetentionTimes(d::TdfsMassSpecData) = d.retention_time
getLowMzs(d::TdfsMassSpecData) = d.low_mz
getHighMzs(d::TdfsMassSpecData) = d.high_mz
getTICs(d::TdfsMassSpecData) = d.tic
getCenterMzs(d::TdfsMassSpecData) = d.center_mz
getIsolationWidthMzs(d::TdfsMassSpecData) = d.isolation_width
getMsOrders(d::TdfsMassSpecData) = d.ms_order
getCycleIdxs(d::TdfsMassSpecData) = d.cycle_idx
getImScans(d::TdfsMassSpecData) = d.im_scan
getFrameIds(d::TdfsMassSpecData) = d.frame_id
getCollisionEnergyEvs(d::TdfsMassSpecData) = d.collision_energy_ev
