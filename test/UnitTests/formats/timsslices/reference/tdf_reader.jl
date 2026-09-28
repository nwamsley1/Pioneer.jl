# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
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

# Prototype pure-Julia reader for Bruker timsTOF .d (TDF) bundles.
#
# Format references (all open-source, no Bruker SDK):
#   - OpenTimsTDF docs/format/02-tdf-bin-block-stream.md, 03-frame-payload-encoding.md, 04-calibration.md
#   - alphatims 1.0.8 bruker.py (process_frame, parse_decompressed_bruker_binary_type2)
#   - timsrust crates/timsrust-tdf/src/calibration.rs
#
# Usage:  julia --project=<env with SQLite, CodecZstd> tdf_reader.jl <run.d> [max_frames]

using SQLite, DataFrames, CodecZstd

struct TdfBundle
    dir::String
    db::SQLite.DB
    bin::IOStream
    meta::Dict{String,String}
    frames::DataFrame
    compression::Int
end

function open_tdf(dir::AbstractString)
    db = SQLite.DB(joinpath(dir, "analysis.tdf"))
    meta = Dict{String,String}()
    for r in DBInterface.execute(db, "SELECT Key, Value FROM GlobalMetadata")
        meta[string(r[:Key])] = string(r[:Value])
    end
    frames = DataFrame(DBInterface.execute(db, "SELECT * FROM Frames ORDER BY Id"))
    bin = open(joinpath(dir, "analysis.tdf_bin"), "r")
    TdfBundle(String(dir), db, bin, meta, frames, parse(Int, meta["TimsCompressionType"]))
end

# ---------------------------------------------------------------------------
# Binary block stream: seek(TimsId); u32 block_size, u32 scan_count, payload
# ---------------------------------------------------------------------------
function read_block(b::TdfBundle, tims_id::Integer)
    seek(b.bin, tims_id)
    block_size = read(b.bin, UInt32)
    scan_count = read(b.bin, UInt32)
    payload = block_size > 8 ? read(b.bin, Int(block_size) - 8) : UInt8[]
    return Int(scan_count), payload
end

"""
Decode a codec-2 (zstd + byte-transpose) frame payload.
Returns (scans, tofs, intensities) as UInt32 vectors, 0-based scan and tof indices.
"""
function decode_codec2(payload::Vector{UInt8}, num_scans::Int, num_peaks::Int)
    scans = Vector{UInt32}(undef, num_peaks)
    tofs  = Vector{UInt32}(undef, num_peaks)
    ints  = Vector{UInt32}(undef, num_peaks)
    num_peaks == 0 && return scans, tofs, ints
    inner = transcode(ZstdDecompressor, payload)
    n = length(inner) ÷ 4
    length(inner) == 4 * (num_scans + 2 * num_peaks) ||
        error("codec2 length invariant violated: $(length(inner)) != 4*($num_scans + 2*$num_peaks)")
    # byte-transposed u32 stream: four byte columns of length n
    logical = Vector{UInt32}(undef, n)
    @inbounds for i in 1:n
        logical[i] = UInt32(inner[i]) | (UInt32(inner[n+i]) << 8) |
                     (UInt32(inner[2n+i]) << 16) | (UInt32(inner[3n+i]) << 24)
    end
    Int(logical[1]) == num_scans || @warn "header scan count $(logical[1]) != NumScans $num_scans"
    pos = num_scans + 1          # 1-based index of first peak-stream word
    k = 0
    peaks_done = 0
    @inbounds for scan in 0:num_scans-1
        peaks_in_scan = scan + 1 < num_scans ? Int(logical[scan + 2] ÷ 2) : num_peaks - peaks_done
        accum = typemax(UInt32)   # so that a stored delta of 1 gives tof 0
        for _ in 1:peaks_in_scan
            accum += logical[pos]; pos += 1     # UInt32 wrapping add
            inten = logical[pos]; pos += 1
            k += 1
            scans[k] = UInt32(scan); tofs[k] = accum; ints[k] = inten
        end
        peaks_done += peaks_in_scan
    end
    k == num_peaks || error("decoded $k peaks, expected $num_peaks")
    return scans, tofs, ints
end

function read_frame(b::TdfBundle, frame_id::Integer)
    row = b.frames[frame_id, :]
    @assert row.Id == frame_id
    scan_count, payload = read_block(b, row.TimsId)
    num_scans = Int(row.NumScans); num_peaks = Int(row.NumPeaks)
    scan_count == num_scans || @warn "block scan_count $scan_count != Frames.NumScans $num_scans (frame $frame_id)"
    b.compression == 2 || error("only TimsCompressionType 2 implemented (got $(b.compression))")
    return decode_codec2(payload, num_scans, num_peaks)
end

# ---------------------------------------------------------------------------
# Calibration (open linear models; the proprietary MzCalibration polynomial is undocumented)
# ---------------------------------------------------------------------------
struct LinearMzCal
    intercept::Float64   # sqrt(mz) at tof index 0
    slope::Float64       # d sqrt(mz) / d tof
end
tof_to_mz(c::LinearMzCal, tof) = (c.intercept + c.slope * tof)^2
mz_to_tof(c::LinearMzCal, mz) = (sqrt(mz) - c.intercept) / c.slope

"Boundary variant: assumes tof 0 <-> MzAcqRangeLower and DigitizerNumSamples <-> MzAcqRangeUpper."
function boundary_mz_cal(meta::Dict{String,String})
    lo = parse(Float64, meta["MzAcqRangeLower"]); hi = parse(Float64, meta["MzAcqRangeUpper"])
    n = parse(Float64, meta["DigitizerNumSamples"])
    if get(meta, "AcquisitionSoftware", "") == "Bruker otofControl"
        lo -= 5; hi += 5
    end
    LinearMzCal(sqrt(lo), (sqrt(hi) - sqrt(lo)) / n)
end

struct LinearImCal
    intercept::Float64   # 1/K0 at scan 0 (= OneOverK0AcqRangeUpper)
    slope::Float64
end
scan_to_im(c::LinearImCal, scan) = c.intercept + c.slope * scan
function boundary_im_cal(meta::Dict{String,String}, frames::DataFrame)
    lo = parse(Float64, meta["OneOverK0AcqRangeLower"]); hi = parse(Float64, meta["OneOverK0AcqRangeUpper"])
    smax = maximum(frames.NumScans)
    LinearImCal(hi, (lo - hi) / smax)
end

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
function main(dir::AbstractString, max_frames::Int)
    b = open_tdf(dir)
    println("== ", dir)
    for k in ("InstrumentName", "AcquisitionSoftware", "AcquisitionSoftwareVersion", "SchemaVersionMajor",
              "SchemaVersionMinor", "TimsCompressionType", "MzAcqRangeLower", "MzAcqRangeUpper",
              "DigitizerNumSamples", "OneOverK0AcqRangeLower", "OneOverK0AcqRangeUpper", "MaxNumPeaksPerScan")
        haskey(b.meta, k) && println("  ", k, " = ", b.meta[k])
    end
    fr = b.frames
    println("  frames: ", nrow(fr), "  MS1: ", count(==(0), fr.MsMsType), "  DIA(9): ", count(==(9), fr.MsMsType),
            "  PASEF(8): ", count(==(8), fr.MsMsType), "  RT range s: ", extrema(fr.Time))
    println("  NumScans max: ", maximum(fr.NumScans), "  total NumPeaks: ", sum(fr.NumPeaks))
    mzcal = boundary_mz_cal(b.meta); imcal = boundary_im_cal(b.meta, fr)
    println("  boundary mz cal: sqrt-intercept=", mzcal.intercept, " slope=", mzcal.slope,
            " -> mz(0)=", tof_to_mz(mzcal, 0), " mz(N)=", tof_to_mz(mzcal, parse(Float64, b.meta["DigitizerNumSamples"])))
    println("  boundary im cal: im(0)=", scan_to_im(imcal, 0), " im(max)=", scan_to_im(imcal, maximum(fr.NumScans)))

    nfr = min(max_frames, nrow(fr))
    t0 = time(); npk = 0
    for fid in 1:nfr
        scans, tofs, ints = read_frame(b, fid)
        npk += length(tofs)
        if fid <= 3
            println("  frame $fid MsMsType=$(fr.MsMsType[fid]) NumPeaks=$(fr.NumPeaks[fid]) first5=",
                    collect(zip(scans[1:min(5,end)], tofs[1:min(5,end)], ints[1:min(5,end)])),
                    " last3=", collect(zip(scans[max(1,end-2):end], tofs[max(1,end-2):end], ints[max(1,end-2):end])))
            length(ints) > 0 && println("     sum(int)=", sum(Int, ints), " SummedIntensities=", fr.SummedIntensities[fid],
                    "  mz range=", tof_to_mz(mzcal, minimum(tofs)), "..", tof_to_mz(mzcal, maximum(tofs)))
        end
    end
    dt = time() - t0
    println("  decoded $nfr frames, $npk peaks in $(round(dt, digits=2)) s  ($(round(npk/1e6/dt, digits=1)) M peaks/s)")
    close(b.bin)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], length(ARGS) >= 2 ? parse(Int, ARGS[2]) : typemax(Int))
end
