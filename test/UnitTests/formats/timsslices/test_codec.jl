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

# Synthetic tests of the codec layer: planes, word stream, quantisation, zstd, raw codec-2 blocks.

"Build a raw Bruker codec-2 payload from per-scan (bin, intensity) lists (bins 0-based, sorted per scan)."
function make_codec2_payload(scans::Vector{Vector{Tuple{Int, Int}}}, ctx::TS.ZstdCtx; level = 3)
    ns = length(scans)
    words = UInt32[UInt32(ns)]
    for s in 1:ns-1; push!(words, UInt32(2 * length(scans[s]))); end
    for s in 1:ns
        acc = typemax(UInt32)
        for (b, i) in scans[s]
            push!(words, UInt32(b) - acc); push!(words, UInt32(i)); acc = UInt32(b)
        end
    end
    n = length(words)
    planes = UInt8[]; TS.transpose!(planes, words, n)
    z = UInt8[]; nb = TS.zstd_compress!(z, ctx, planes, 4n, level)
    z[1:nb], sum(length, scans)
end

@testset "planes" begin
    rng = MersenneTwister(1)
    for n in (0, 1, 7, 8, 9, 31, 32, 33, 100_000)
        w = rand(rng, UInt32, n)
        planes = UInt8[]; TS.transpose!(planes, w, n)
        @test length(planes) >= 4n
        n > 0 && @test planes[1:n] == UInt8.(w .& 0xff) && planes[3n+1:4n] == UInt8.(w .>> 24)
        back = UInt32[]; TS.untranspose!(back, planes, n)
        @test back[1:n] == w
    end
    @test_throws ArgumentError TS.untranspose!(UInt32[], UInt8[1, 2, 3], 1)
end

@testset "zstd" begin
    ctx = TS.ZstdCtx()
    for n in (0, 1, 1000, 1_000_000)
        src = rand(MersenneTwister(n), UInt8(0):UInt8(3), n)
        z = UInt8[]; nb = TS.zstd_compress!(z, ctx, src, n, 3)
        out = UInt8[]; m = TS.zstd_decompress!(out, ctx, view(z, 1:nb), n)
        @test m == n && out[1:n] == src
        n > 0 && @test TS.zstd_frame_content_size(view(z, 1:nb)) == n
    end
    z = UInt8[]; nb = TS.zstd_compress!(z, ctx, UInt8[1, 2, 3], 3, 1)
    @test_throws ErrorException TS.zstd_decompress!(UInt8[], ctx, view(z, 1:nb), 4)
end

"Random FrameSlices: n_slices slices with sorted fractional positions and positive intensities."
function random_slices(rng, n_slices; max_peaks = 50, max_bin = 640_000, int_max = 1e6)
    fs = FrameSlices()
    for j in 1:n_slices
        np = rand(rng, 0:max_peaks)
        pos = sort(rand(rng, np) .* max_bin)
        for k in 1:np
            TS.push_peak!(fs, pos[k], rand(rng) * int_max + 0.5)
        end
        TS.end_slice!(fs, j * 8, 1)   # empty slices are not recorded
    end
    fs
end

@testset "words and blocks" begin
    rng = MersenneTwister(2)
    codec = BlockCodec()
    for trial in 1:30
        fs = random_slices(rng, rand(rng, 0:20))
        scale = rand(rng, (1, 2, 8, 256))
        blk = SliceBlock(); quantize!(blk, fs, scale, 1.0)
        @test blk.n_slices == fs.n_slices
        # peaks sorted and unique within a slice, intensities > 0
        for j in 1:blk.n_slices
            r = TS.slice_range(blk, j)
            @test issorted(view(blk.bin, r); lt = <) && allunique(view(blk.bin, r))
            @test all(>(0), view(blk.intensity, r))
        end
        words = UInt32[]; n = TS.encode_words!(words, blk)
        @test n == 1 + blk.n_slices + 2 * TS.n_peaks(blk)
        back = SliceBlock(); TS.decode_words!(back, words, n)
        @test back.n_slices == blk.n_slices && back.ptr == blk.ptr && back.bin == blk.bin && back.intensity == blk.intensity
        nw, nb = encode_block!(codec, blk, 3)
        back2 = SliceBlock(); decode_block!(back2, codec, codec.zbuf[1:nb], nw)
        @test back2.ptr == blk.ptr && back2.bin == blk.bin && back2.intensity == blk.intensity
        # per-slice blocks (the on-disk unit): each slice round-trips on its own
        sb = SliceBuffer()
        for j in 1:blk.n_slices
            r = TS.slice_range(blk, j)
            nb = encode_slice!(codec, blk, j, 3)
            isempty(r) && @test nb == 0
            decode_slice!(sb, codec, codec.zbuf[1:nb], length(r))
            @test sb.n_peaks == length(r) && sb.bin[1:sb.n_peaks] == blk.bin[r] && sb.intensity[1:sb.n_peaks] == blk.intensity[r]
        end
    end
    # empty frame
    blk = SliceBlock(); quantize!(blk, FrameSlices(), 1, 1.0)
    nw, nb = encode_block!(codec, blk, 3)
    @test nw == 1
    back = SliceBlock(); decode_block!(back, codec, codec.zbuf[1:nb], nw)
    @test back.n_slices == 0 && isempty(back.bin)
    # extremes: bin near 2^20 * 256, intensity at typemax
    fs = FrameSlices(); TS.push_peak!(fs, 1_048_575.0, 4.294967295e9); TS.push_peak!(fs, 1_048_575.4, 10.0); TS.end_slice!(fs, 0, 1)
    quantize!(blk, fs, 256, 1.0)
    @test blk.bin == [UInt32(round(1_048_575.0 * 256)), UInt32(round(1_048_575.4 * 256))] && blk.intensity[1] == typemax(UInt32)
    nw, nb = encode_block!(codec, blk, 3); decode_block!(back, codec, codec.zbuf[1:nb], nw)
    @test back.bin == blk.bin && back.intensity == blk.intensity
    # bad streams
    @test_throws ArgumentError TS.decode_words!(back, UInt32[3, 1], 2)
    @test_throws ArgumentError TS.decode_words!(back, UInt32[1, 2, 5, 7], 4)
end

@testset "quantize: merge, drop, scale" begin
    fs = FrameSlices()
    TS.push_peak!(fs, 100.2, 5.0); TS.push_peak!(fs, 100.4, 7.0); TS.push_peak!(fs, 101.6, 0.3); TS.push_peak!(fs, 200.0, 2.0)
    TS.end_slice!(fs, 0, 1)
    TS.push_peak!(fs, 50.0, 1.0); TS.end_slice!(fs, 8, 1)
    blk = SliceBlock(); quantize!(blk, fs, 1, 1.0)
    @test blk.n_slices == 2
    # 100.2 and 100.4 merge into bin 100 (5 + 7); 101.6 -> bin 102 with intensity round(0.3) = 0 -> dropped
    @test blk.bin == UInt32[100, 200, 50] && blk.intensity == UInt32[12, 2, 1]
    @test TS.slice_range(blk, 1) == 1:2 && TS.slice_range(blk, 2) == 3:3
    # int_scale keeps the small one; bin_scale 2 separates the pair
    quantize!(blk, fs, 2, 10.0)
    @test blk.bin == UInt32[200, 201, 203, 400, 100] && blk.intensity == UInt32[50, 70, 3, 20, 10]
    # a slice that loses all its peaks stays as an empty slice
    fs2 = FrameSlices(); TS.push_peak!(fs2, 1.0, 0.1); TS.end_slice!(fs2, 0, 1); TS.push_peak!(fs2, 2.0, 1.0); TS.end_slice!(fs2, 8, 1)
    quantize!(blk, fs2, 1, 1.0)
    @test blk.n_slices == 2 && blk.ptr == Int32[1, 1, 2]
    # a slightly negative position (footprint below bin 0 at the range edge) clamps to bin 0
    fs4 = FrameSlices(); TS.push_peak!(fs4, -0.05, 3.0); TS.push_peak!(fs4, 0.6, 4.0); TS.end_slice!(fs4, 0, 1)
    quantize!(blk, fs4, 256, 1.0)
    @test blk.bin == UInt32[0, 154] && blk.intensity == UInt32[3, 4]
    quantize!(blk, fs4, 1, 1.0)
    @test blk.bin == UInt32[0, 1]
    # unsorted input raises
    fs3 = FrameSlices(); TS.push_peak!(fs3, 5.0, 1.0); TS.push_peak!(fs3, 4.0, 1.0); TS.end_slice!(fs3, 0, 1)
    @test_throws ArgumentError quantize!(blk, fs3, 1, 1.0)
end

@testset "FrameSlices bookkeeping" begin
    fs = FrameSlices()
    @test !TS.end_slice!(fs, 0, 1)          # nothing pushed -> not recorded
    TS.push_peak!(fs, 1.0, 1.0); TS.discard_slice!(fs)
    @test TS.n_peaks(fs) == 0 && !TS.end_slice!(fs, 0, 1)
    TS.push_peak!(fs, 1.0, 1.0); @test TS.end_slice!(fs, 16, 2)
    @test fs.n_slices == 1 && fs.scan == [16] && fs.window == [2] && fs.ptr == Int32[1, 2]
    TS.reset!(fs); @test fs.n_slices == 0 && fs.ptr == Int32[1]
end

@testset "raw codec-2 decode" begin
    ctx = TS.ZstdCtx()
    scans = [[(0, 10), (5, 20), (700_000, 1)], Tuple{Int, Int}[], [(97_526, 38), (98_219, 114), (111_040, 33)], [(3, 4)]]
    payload, np = make_codec2_payload(scans, ctx)
    buf = FrameBuffer()
    TS.decode_codec2!(buf, payload, length(scans), np)
    @test buf.n_scans == 4 && buf.n_peaks == 7
    @test buf.scan_start[1:5] == Int32[1, 4, 4, 7, 8]
    @test buf.tof[1:7] == UInt32[0, 5, 700_000, 97_526, 98_219, 111_040, 3]
    @test buf.intensity[1:7] == UInt32[10, 20, 1, 38, 114, 33, 4]
    @test TS.scan_range(buf, 1) == 4:3 && TS.scan_range(buf, 2) == 4:6
    # a delta of 1 encodes bin 0 (accumulator wraps from 0xFFFFFFFF): covered by scan 1's first peak above
    # wrong peak count -> length invariant fails
    @test_throws ErrorException TS.decode_codec2!(buf, payload, length(scans), np + 1)
    # empty frame
    TS.decode_codec2!(buf, UInt8[], 3, 0)
    @test buf.n_peaks == 0 && buf.scan_start[1:4] == Int32[1, 1, 1, 1]
end
