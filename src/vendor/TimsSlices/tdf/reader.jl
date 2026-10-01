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

# A .d bundle opened for reading: SQLite tables in memory; analysis.tdf_bin read block by block with positioned
# reads into per-thread buffers (a memory map would count every touched page of the 3-7 GB file as resident).

struct TdfFile
    dir::String
    meta::Dict{String, String}
    frames::FrameTable
    dia::DiaScheme
    mz_cal::LinearMzCal
    mz_cal_resid_ppm::Vector{Float64}
    timebase_ns::Float64        # digitizer sample length = one TOF bin (MzCalibration.DigitizerTimebase)
    im_cal::LinearImCal
    ce_ramp::CeRamp
    bin::IOStream               # analysis.tdf_bin (positioned reads; the stream position is never used)
    bin_size::Int64
    compression::Int
    max_scans::Int
    max_peaks::Int
end

function open_tdf(dir::AbstractString)
    isdir(dir) || error("not a directory: $dir")
    db = SQLite.DB(joinpath(dir, "analysis.tdf"))
    meta = read_global_metadata(db)
    frames = read_frames(db)
    dia = read_dia_scheme(db, frames)
    mzcal_id = isempty(frames.mz_calibration) ? 1 : Int(frames.mz_calibration[1])
    cal, resid = regressed_mz_cal(db, mzcal_id)
    timebase = Float64(first(DBInterface.execute(db, "SELECT DigitizerTimebase FROM MzCalibration WHERE Id=$mzcal_id"))[:DigitizerTimebase])
    imcal = boundary_im_cal(meta, frames)
    ce = fit_ce_ramp(dia)
    close(db)
    compression = parse(Int, get(meta, "TimsCompressionType", "2"))
    compression == 2 || error("only TimsCompressionType 2 is supported (file has $compression)")
    bin = open(joinpath(dir, "analysis.tdf_bin"), "r")
    TdfFile(String(dir), meta, frames, dia, cal, resid, timebase, imcal, ce, bin, filesize(bin), compression,
            isempty(frames.num_scans) ? 0 : Int(maximum(frames.num_scans)),
            isempty(frames.num_peaks) ? 0 : Int(maximum(frames.num_peaks)))
end

n_frames(f::TdfFile) = length(f.frames)
is_ms1(f::TdfFile, i::Integer) = f.frames.msms_type[i] == 0
is_dia(f::TdfFile, i::Integer) = f.frames.msms_type[i] == 9
"Rows of the frames that are converted: MS1 and diaPASEF MS2 frames, in order."
valid_frames(f::TdfFile) = [i for i in 1:n_frames(f) if is_ms1(f, i) || is_dia(f, i)]
"Windows of frame row i (MS1: one window over all scans)."
function windows(f::TdfFile, i::Integer)
    if is_ms1(f, i)
        return (DiaWindow(0, f.frames.num_scans[i], 0f0, 0f0, 0f0),)
    else
        g = f.dia.frame_group[i]
        g == 0 && error("frame $(f.frames.id[i]) has MsMsType 9 but no DiaFrameMsMsInfo row")
        return f.dia.groups[g]
    end
end
"Number of TOF bins on the digitiser axis (upper bound for any tof value)."
n_bins(f::TdfFile) = parse(Int, get(f.meta, "DigitizerNumSamples", "0"))

"Read `n` bytes at byte offset `off` of the block file into `dst[1:n]` (thread-safe: positioned read)."
function pread!(dst::Vector{UInt8}, f::TdfFile, off::Integer, n::Integer)
    length(dst) < n && resize!(dst, n)
    n == 0 && return dst
    @static if Sys.isunix()
        done = 0
        GC.@preserve dst while done < n
            r = ccall(:pread, Cssize_t, (Cint, Ptr{UInt8}, Csize_t, Int64), fd(f.bin), pointer(dst) + done, n - done, off + done)
            r > 0 || error("pread of analysis.tdf_bin at offset $(off + done) failed: $(Libc.strerror(Libc.errno()))")
            done += r
        end
    else
        lock(f.bin) do
            seek(f.bin, off); unsafe_read(f.bin, pointer(dst), n)
        end
    end
    dst
end

"""
    raw_block!(dst, f, i) -> (payload view, block_size, scan_count)

Read frame row i's block (header + zstd payload) into `dst`.
"""
function raw_block!(dst::Vector{UInt8}, f::TdfFile, i::Integer)
    off = f.frames.tims_id[i]
    off + 8 <= f.bin_size || error("frame $(f.frames.id[i]): block offset $off beyond file")
    pread!(dst, f, off, 8)
    block_size = GC.@preserve dst unsafe_load(Ptr{UInt32}(pointer(dst)))
    scan_count = GC.@preserve dst unsafe_load(Ptr{UInt32}(pointer(dst) + 4))
    block_size >= 8 || error("frame $(f.frames.id[i]): block_size $block_size < 8")
    off + block_size <= f.bin_size || error("frame $(f.frames.id[i]): block runs beyond file")
    pread!(dst, f, off, Int(block_size))
    view(dst, 9:Int(block_size)), Int(block_size), Int(scan_count)
end

"Decode frame row i into `buf`."
function read_frame!(buf::FrameBuffer, f::TdfFile, i::Integer)
    payload, _, scan_count = raw_block!(buf.raw, f, i)
    ns = Int(f.frames.num_scans[i]); np = Int(f.frames.num_peaks[i])
    scan_count == ns || error("frame $(f.frames.id[i]): block scan_count $scan_count != NumScans $ns")
    decode_codec2!(buf, payload, ns, np)
end

Base.close(f::TdfFile) = close(f.bin)
