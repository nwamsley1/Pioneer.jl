# ScxsMassSpecData: Pioneer's view of a `.scxs` (SciexWiff.jl) — centroided SCIEX SWATH scans stored as one zstd
# block per scan. Every scan is one Pioneer scan; there is no ion mobility (`getImScans` stays `nothing`), so the
# search takes the same code paths as for Arrow data.
#
# The blocks are byte-identical to TimsSlices' per-slice `.tdfs` blocks, so they decode with
# `TimsSlices.decode_slice!` into the caller-owned `PeakDecodeBuffer` (see TdfsMassSpecData.jl for the buffer
# contract). What differs from `.tdfs` is the calibration, which SCIEX records per scan:
#
#     m/z = (cal_a · (bin / bin_scale / 40 − cal_b))²,   intensity = stored / int_scale
#
# These are the same Float32 values SciexWiff writes into the Arrow it produces from the same conversion, so a
# search from the .scxs and from that Arrow see identical peaks.

using Mmap: Mmap

struct ScxsMassSpecData <: MassSpecData
    dir::String
    uid::Int                                    # unique per opened file (PeakDecodeBuffer cache key)
    n::Int
    blocks::Vector{UInt8}                       # memory-mapped blocks.bin
    block_offset::Vector{Int64}
    block_size::Vector{Int32}
    n_peaks::Vector{Int32}
    cal_a::Vector{Float64}
    cal_b::Vector{Float64}
    bin_scale::Float64
    int_scale::Float64
    retention_time::Vector{Float32}             # minutes
    low_mz::Vector{Float32}
    high_mz::Vector{Float32}
    tic::Vector{Float32}
    center_mz::Vector{Union{Missing, Float32}}
    isolation_width::Vector{Union{Missing, Float32}}
    ms_order::Vector{UInt8}
    cycle_idx::Vector{UInt32}
end

const SCXS_FORMAT_VERSION = 1

"""
    ScxsMassSpecData(dir::String)

Open a `<name>.scxs` directory: the scan table is loaded into plain columns, the block file is memory-mapped.
Peaks are decoded per scan with `getPeaks!`.
"""
function ScxsMassSpecData(dir::String)
    meta = JSON.parsefile(joinpath(dir, "meta.json"))
    Int(meta["format_version"]) == SCXS_FORMAT_VERSION ||
        error("$dir: .scxs format version $(meta["format_version"]), this Pioneer reads $SCXS_FORMAT_VERSION")
    t = Arrow.Table(joinpath(dir, "scans.arrow"))
    n = length(t.record)
    blocks = open(joinpath(dir, "blocks.bin"), "r") do io
        Mmap.mmap(io, Vector{UInt8}, filesize(io))
    end
    nan_to_missing(v) = Union{Missing, Float32}[isnan(x) ? missing : Float32(x) for x in v]
    ScxsMassSpecData(
        dir, Threads.atomic_add!(_TDFS_NEXT_UID, 1), n, blocks,
        Vector{Int64}(t.block_offset), Vector{Int32}(t.block_size), Vector{Int32}(t.n_peaks),
        Vector{Float64}(t.cal_a), Vector{Float64}(t.cal_b), Float64(meta["bin_scale"]), Float64(meta["int_scale"]),
        Vector{Float32}(t.retention_time), Vector{Float32}(t.low_mz), Vector{Float32}(t.high_mz), Vector{Float32}(t.tic),
        nan_to_missing(t.center_mz), nan_to_missing(t.isolation_width), Vector{UInt8}(t.ms_order),
        Vector{UInt32}(t.cycle),
    )
end

is_scxs_path(path::AbstractString) = endswith(path, ".scxs") && isdir(path)

Base.length(d::ScxsMassSpecData) = d.n

function _scxs_decode!(buf::PeakDecodeBuffer, d::ScxsMassSpecData, scan::Int)
    np = Int(@inbounds d.n_peaks[scan])
    off = Int(@inbounds d.block_offset[scan]); bs = Int(@inbounds d.block_size[scan])
    off + bs <= length(d.blocks) || error("$(d.dir): scan $scan block runs beyond blocks.bin")
    TimsSlices.decode_slice!(buf.sb, buf.codec, view(d.blocks, off+1:off+bs), np)
    if length(buf.mz) < np
        resize!(buf.mz, np); resize!(buf.intensity, np)
    end
    a = @inbounds d.cal_a[scan]; b = @inbounds d.cal_b[scan]
    k = 1 / (40 * d.bin_scale); inv_int = 1 / d.int_scale
    bins = buf.sb.bin; ints = buf.sb.intensity; mz = buf.mz; it = buf.intensity
    @inbounds for i in 1:np
        v = a * (Float64(bins[i]) * k - b)
        mz[i] = Float32(v * v)
        it[i] = Float32(Float64(ints[i]) * inv_int)
    end
    buf.file_uid = d.uid
    buf.scan = scan
    buf
end

function getPeaks!(buf::PeakDecodeBuffer, d::ScxsMassSpecData, scan_idx::Integer)
    scan = Int(scan_idx)
    1 <= scan <= d.n || throw(BoundsError(d, scan))
    (buf.scan == scan && buf.file_uid == d.uid) || _scxs_decode!(buf, d, scan)
    np = Int(@inbounds d.n_peaks[scan])
    return (view(buf.mz, 1:np)::_TdfsPeakView, view(buf.intensity, 1:np)::_TdfsPeakView)
end

getMzArray(::ScxsMassSpecData, ::Integer) = error("ScxsMassSpecData decodes into a caller-owned buffer; use getPeaks!(buf, spectra, scan)")
getIntensityArray(::ScxsMassSpecData, ::Integer) = error("ScxsMassSpecData decodes into a caller-owned buffer; use getPeaks!(buf, spectra, scan)")
getMzArrays(::ScxsMassSpecData) = error("ScxsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getPeaks! per scan")
getIntensityArrays(::ScxsMassSpecData) = error("ScxsMassSpecData decodes peaks per scan; use getPeakCounts / getPeakCount for lengths or getPeaks! per scan")

getPeakCount(d::ScxsMassSpecData, scan_idx::Integer) = Int(d.n_peaks[scan_idx])
getPeakCounts(d::ScxsMassSpecData) = d.n_peaks

getRetentionTime(d::ScxsMassSpecData, i::Integer) = d.retention_time[i]::Float32
getLowMz(d::ScxsMassSpecData, i::Integer) = d.low_mz[i]::Float32
getHighMz(d::ScxsMassSpecData, i::Integer) = d.high_mz[i]::Float32
getTIC(d::ScxsMassSpecData, i::Integer) = d.tic[i]::Float32
getCenterMz(d::ScxsMassSpecData, i::Integer) = d.center_mz[i]::Union{Missing, Float32}
getIsolationWidthMz(d::ScxsMassSpecData, i::Integer) = d.isolation_width[i]::Union{Missing, Float32}
getMsOrder(d::ScxsMassSpecData, i::Integer) = d.ms_order[i]::UInt8
getCycleIdx(d::ScxsMassSpecData, i::Integer) = d.cycle_idx[i]::UInt32
getPrecursorMz(d::ScxsMassSpecData, i::Integer) = getCenterMz(d, i)
getCollisionEnergyEv(::ScxsMassSpecData, ::Integer)::Float32 = 0f0     # as for mzML-converted Arrow
getScanHeader(::ScxsMassSpecData, ::Integer) = ""
getScanNumber(::ScxsMassSpecData, i::Integer) = Int32(i)
getBasePeakMz(::ScxsMassSpecData, ::Integer) = missing
getBasePeakIntensity(::ScxsMassSpecData, ::Integer) = missing

getRetentionTimes(d::ScxsMassSpecData) = d.retention_time
getLowMzs(d::ScxsMassSpecData) = d.low_mz
getHighMzs(d::ScxsMassSpecData) = d.high_mz
getTICs(d::ScxsMassSpecData) = d.tic
getCenterMzs(d::ScxsMassSpecData) = d.center_mz
getIsolationWidthMzs(d::ScxsMassSpecData) = d.isolation_width
getMsOrders(d::ScxsMassSpecData) = d.ms_order
getCycleIdxs(d::ScxsMassSpecData) = d.cycle_idx
getCollisionEnergyEvs(::ScxsMassSpecData) = nothing
