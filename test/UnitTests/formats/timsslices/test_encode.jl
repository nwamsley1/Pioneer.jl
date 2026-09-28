# encode_codec2 (tdf/block.jl) is the inverse of decode_codec2!: synthetic frames round-trip exactly.

"Decode a block produced by encode_codec2 (header + payload) with the reader's decoder."
function decode_block_bytes(block::Vector{UInt8})
    n_scans = Int(reinterpret(UInt32, block[5:8])[1])
    @test Int(reinterpret(UInt32, block[1:4])[1]) == length(block)
    buf = TS.FrameBuffer()
    npk = 0                                   # NumPeaks comes from the Frames table in a real file; recover it here
    if length(block) > 8
        planes = UInt8[]; TS.zstd_decompress!(planes, TS.ZstdCtx(), view(block, 9:length(block)), 0)
        words = UInt32[]; TS.untranspose!(words, planes, length(planes) ÷ 4)
        npk = (length(planes) ÷ 4 - n_scans) ÷ 2
    end
    TS.decode_codec2!(buf, view(block, 9:length(block)), n_scans, npk)
    buf
end

@testset "encode_codec2 round trip" begin
    rng = Random.MersenneTwister(42)
    for trial in 1:25
        n_scans = rand(rng, 1:60)
        counts = [rand(rng) < 0.3 ? 0 : rand(rng, 1:40) for _ in 1:n_scans]   # empty scans included
        scan_start = cumsum(vcat(1, counts))
        tof = UInt32[]; it = UInt32[]
        for c in counts
            bins = sort!(unique(rand(rng, UInt32(0):UInt32(700_000), c)))    # includes bin 0 (delta wraps from -1)
            append!(tof, bins); append!(it, rand(rng, UInt32(1):UInt32(60_000), length(bins)))
        end
        scan_start = cumsum(vcat(1, [length(unique(x)) for x in [tof[scan_start[s]:scan_start[s+1]-1] for s in 1:n_scans]]))
        block = TS.encode_codec2(tof, it, scan_start, n_scans)
        buf = decode_block_bytes(block)
        @test buf.n_scans == n_scans && buf.n_peaks == length(tof)
        @test buf.scan_start[1:n_scans + 1] == scan_start
        @test buf.tof[1:length(tof)] == tof && buf.intensity[1:length(it)] == it
    end
    # a frame without peaks is the 8-byte header alone
    empty = TS.encode_codec2(UInt32[], UInt32[], fill(1, 11), 10)
    @test length(empty) == 8 && reinterpret(UInt32, empty) == UInt32[8, 10]
    buf = decode_block_bytes(empty)
    @test buf.n_peaks == 0 && all(==(1), buf.scan_start[1:11])
end
