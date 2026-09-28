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

# Raw codec-2 frame block: zstd -> byte planes -> word stream -> per-scan prefix sums.
#
# Bruker's layout (TimsCompressionType 2), as decoded here. A frame block in analysis.tdf_bin is an 8-byte header
# (UInt32 block size, UInt32 scan count) followed by a zstd payload. Decompressed, the payload is 4n bytes of
# byte planes (see codec/planes.jl) holding n = NumScans + 2 * NumPeaks little-endian UInt32 words:
#
#   word 1                    NumScans
#   words 2 .. NumScans       2 x (peaks in scan s), for s = 0 .. NumScans-2 (the last scan gets the remainder)
#   then, scan by scan,       (tof_delta, intensity) pairs
#
# Within a scan the TOF bin is a running sum of the deltas, starting from -1 (0xFFFFFFFF, wrapping), so the
# first delta of every scan is its first bin + 1. Peaks within a scan are in increasing TOF order.

"Decoded peaks of one raw frame, per-thread scratch sized to the file's largest frame."
mutable struct FrameBuffer
    n_scans::Int
    n_peaks::Int
    scan_start::Vector{Int32}    # n_scans + 1 entries; peaks of 0-based scan s are scan_start[s+1]:scan_start[s+2]-1
    tof::Vector{UInt32}          # TOF bin (0-based)
    intensity::Vector{UInt32}
    raw::Vector{UInt8}           # the compressed block as read from disk
    planes::Vector{UInt8}        # zstd output (byte planes)
    words::Vector{UInt32}
    zstd::ZstdCtx
end
FrameBuffer() = FrameBuffer(0, 0, Int32[], UInt32[], UInt32[], UInt8[], UInt8[], UInt32[], ZstdCtx())

@inline scan_range(b::FrameBuffer, s::Integer) = (b.scan_start[s + 1]):(b.scan_start[s + 2] - 1)

"""
    decode_codec2!(buf, payload, n_scans, n_peaks)

Decode one TimsCompressionType-2 payload (the bytes after the 8-byte block header) into `buf`.
"""
function decode_codec2!(buf::FrameBuffer, payload::AbstractVector{UInt8}, n_scans::Integer, n_peaks::Integer)
    buf.n_scans = n_scans; buf.n_peaks = n_peaks
    length(buf.scan_start) < n_scans + 1 && resize!(buf.scan_start, n_scans + 1)
    length(buf.tof) < n_peaks && (resize!(buf.tof, n_peaks); resize!(buf.intensity, n_peaks))
    if n_peaks == 0
        fill!(view(buf.scan_start, 1:n_scans + 1), Int32(1))
        return buf
    end
    n = n_scans + 2n_peaks
    zstd_decompress!(buf.planes, buf.zstd, payload, 4n)
    untranspose!(buf.words, buf.planes, n)
    w = buf.words
    Int(w[1]) == n_scans || error("block header scan count $(w[1]) != Frames.NumScans $n_scans")
    # scan headers: word s+2 = 2 * peaks in scan s (s = 0 .. n_scans-2); the last scan takes the remainder
    ss = buf.scan_start
    ss[1] = Int32(1)
    @inbounds for s in 0:n_scans-2
        ss[s + 2] = ss[s + 1] + Int32(w[s + 2] >> 1)
    end
    ss[n_scans + 1] = Int32(n_peaks + 1)
    ss[n_scans] <= ss[n_scans + 1] || error("scan headers exceed NumPeaks ($(ss[n_scans] - 1) > $n_peaks)")
    # peak pairs: the TOF bin is the running sum of the deltas, restarting at -1 (typemax, wrapping) every scan
    tof = buf.tof; it = buf.intensity
    pos = n_scans + 1
    @inbounds for s in 0:n_scans-1
        acc = typemax(UInt32)
        for k in ss[s + 1]:ss[s + 2] - 1
            acc += w[pos]; tof[k] = acc; it[k] = w[pos + 1]; pos += 2
        end
    end
    buf
end
