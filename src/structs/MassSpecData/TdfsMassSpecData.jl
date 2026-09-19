# TdfsMassSpecData: Pioneer's view of a `.tdfs` (TimsSlices.jl) — IM-smoothed, m/z-centroided timsTOF slices
# stored as one zstd block per slice. Every slice is one Pioneer scan.
#
# Peak arrays are decoded on demand: `getMzArray(d, scan)` / `getIntensityArray(d, scan)` decode the scan's block
# into the calling THREAD's scratch (one slice, ~5 µs for an MS2 slice, 30-50 µs for an MS1 slice) and return views
# into it. A view is valid until the same thread fetches a DIFFERENT scan; fetching the same scan again (the usual
# mz-then-intensity pair) reuses the decoded slice. No caller in Pioneer holds a view across another scan's fetch
# (audit 2026-09-18, docs/bruker_timstof_progress.md §9); a task migrates threads only at a yield, and the hot
# loops have none between fetch and use. Keep it that way.
#
# The vectors are `Vector{Union{Missing,Float32}}` so the views match the `AbstractArray{Union{Missing,Float32}}`
# signatures of run_fused! / prepare_scan_peaks! / the fragment-index scorer unchanged (element types are
# invariant, so a Float32 view would not). No element is ever `missing`.
#
# m/z and intensity are the same Float32 values `TimsSlices.expand` writes into the slice Arrow
# (`Float32((a + b·bin/k)²)`, `Float32(stored / int_scale)`), so a search from the .tdfs and from its expanded
# Arrow see identical peaks.

using TimsSlices: TimsSlices, TdfsFile, SliceBuffer, BlockCodec

mutable struct TdfsSliceScratch
    scan::Int                                   # scan currently decoded (0 = none)
    sb::SliceBuffer
    codec::BlockCodec
    mz::Vector{Union{Missing, Float32}}
    intensity::Vector{Union{Missing, Float32}}
end
TdfsSliceScratch() = TdfsSliceScratch(0, SliceBuffer(), BlockCodec(), Union{Missing, Float32}[], Union{Missing, Float32}[])

struct TdfsMassSpecData <: MassSpecData
    file::TdfsFile
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
    scratch::Vector{TdfsSliceScratch}           # one per thread id
end

"""
    TdfsMassSpecData(dir::String)

Open a `<name>.tdfs` directory. The slice side table is loaded into plain columns; the block file is memory-mapped;
per-thread decode scratch is allocated once (it grows to the largest slice decoded on that thread).
"""
function TdfsMassSpecData(dir::String)
    file = TimsSlices.open_tdfs(dir)
    sl = file.slices
    n = length(sl.frame_row)
    nan_to_missing(v) = Union{Missing, Float32}[isnan(x) ? missing : Float32(x) for x in v]
    mz_lo = Float32(file.meta["mz_lo"]); mz_hi = Float32(file.meta["mz_hi"])
    TdfsMassSpecData(
        file, n,
        Vector{Float32}(sl.retention_time), fill(mz_lo, n), fill(mz_hi, n), Vector{Float32}(sl.tic),
        nan_to_missing(sl.center_mz), nan_to_missing(sl.isolation_width),
        Vector{Float32}(sl.collision_energy_ev), Vector{UInt8}(sl.ms_order), Vector{UInt32}(sl.cycle_idx),
        Vector{UInt16}(sl.im_scan), Vector{Int32}(sl.frame_id), Vector{Int32}(sl.n_peaks),
        [TdfsSliceScratch() for _ in 1:Threads.maxthreadid()],
    )
end

is_tdfs_path(path::AbstractString) = endswith(path, ".tdfs") && isdir(path)

Base.length(d::TdfsMassSpecData) = d.n

# ---- peak arrays: decode on demand into the thread's scratch ------------------------------------------------

@inline function _tdfs_scratch(d::TdfsMassSpecData)
    tid = Threads.threadid()
    tid <= length(d.scratch) || throw(ErrorException("TdfsMassSpecData: thread id $tid exceeds the $(length(d.scratch)) scratch slots allocated at open"))
    @inbounds d.scratch[tid]
end

function _tdfs_decode!(s::TdfsSliceScratch, d::TdfsMassSpecData, scan::Int)
    np = Int(@inbounds d.n_peaks[scan])
    TimsSlices.read_slice!(s.sb, s.codec, d.file, scan)
    if length(s.mz) < np
        resize!(s.mz, np); resize!(s.intensity, np)
    end
    f = d.file; bins = s.sb.bin; ints = s.sb.intensity; mz = s.mz; it = s.intensity
    @inbounds for i in 1:np
        mz[i] = Float32(TimsSlices.bin_to_mz(f, bins[i]))
        it[i] = Float32(TimsSlices.stored_to_intensity(f, ints[i]))
    end
    s.scan = scan
    s
end

@inline function _tdfs_slice(d::TdfsMassSpecData, scan_idx::Integer)
    scan = Int(scan_idx)
    1 <= scan <= d.n || throw(BoundsError(d, scan))
    s = _tdfs_scratch(d)
    s.scan == scan || _tdfs_decode!(s, d, scan)
    s
end

const _TdfsPeakView = SubArray{Union{Missing, Float32}, 1, Vector{Union{Missing, Float32}}, Tuple{UnitRange{Int64}}, true}

function getMzArray(d::TdfsMassSpecData, scan_idx::Integer)
    s = _tdfs_slice(d, scan_idx)
    return view(s.mz, 1:Int(@inbounds d.n_peaks[Int(scan_idx)]))::_TdfsPeakView
end
function getIntensityArray(d::TdfsMassSpecData, scan_idx::Integer)
    s = _tdfs_slice(d, scan_idx)
    return view(s.intensity, 1:Int(@inbounds d.n_peaks[Int(scan_idx)]))::_TdfsPeakView
end

getMzArrays(::TdfsMassSpecData) = error("TdfsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getMzArray per scan")
getIntensityArrays(::TdfsMassSpecData) = error("TdfsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getIntensityArray per scan")

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
