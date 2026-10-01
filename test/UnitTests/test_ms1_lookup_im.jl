# build_scan_to_ms1: nearest MS1 scan by RT on plain files; nearest MS1 frame at the nearest IM scan on
# ion-mobility packet / slice files.

using Test
using Arrow
using DataFrames
using Pioneer
using Pioneer: build_scan_to_ms1, BasicMassSpecData

function _write_ms(path, df)
    Arrow.write(path, df); BasicMassSpecData(path)
end
_base(n) = DataFrame(mz_array = [Union{Missing,Float32}[500f0] for _ in 1:n], intensity_array = [Union{Missing,Float32}[1f0] for _ in 1:n],
                     scanHeader = fill("", n), scanNumber = Int32.(1:n), packetType = zeros(Int32, n),
                     lowMz = fill(100f0, n), highMz = fill(1700f0, n), TIC = ones(Float32, n),
                     collisionEnergyField = Vector{Union{Missing,Float32}}(fill(30f0, n)), collisionEnergyEvField = zeros(Float32, n))

@testset "build_scan_to_ms1 — plain file: nearest MS1 by RT" begin
    mktempdir() do d
        n = 7
        df = _base(n)
        df.retentionTime = Float32[1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0]
        df.msOrder = UInt8[1, 2, 2, 1, 2, 2, 1]
        df.centerMz = Vector{Union{Missing,Float32}}([missing, 500f0, 600f0, missing, 500f0, 600f0, missing])
        df.isolationWidthMz = Vector{Union{Missing,Float32}}([missing, 25f0, 25f0, missing, 25f0, 25f0, missing])
        df.cycle_idx = Int32[1, 1, 1, 2, 2, 2, 3]
        ms = _write_ms(joinpath(d, "plain.arrow"), df)
        m = build_scan_to_ms1(ms)
        @test m == Int32[1, 1, 1, 4, 4, 4, 7]
        # no MS1 at all -> zeros
        df2 = copy(df); df2.msOrder .= 0x02; df2.centerMz .= 500f0; df2.isolationWidthMz .= 25f0
        @test all(iszero, build_scan_to_ms1(_write_ms(joinpath(d, "noms1.arrow"), df2)))
    end
end

@testset "build_scan_to_ms1 — packet file: nearest MS1 frame at the nearest IM scan" begin
    mktempdir() do d
        # frame 1 = MS1 with slices at scans 0, 8, 16, 24; frame 2 = MS2 window rows at scans 5, 13, 21;
        # frame 3 = MS1 slices 0, 8, 16, 24 (next cycle); frame 4 = MS2 rows 3, 27
        scans   = UInt16[0, 8, 16, 24,  5, 13, 21,  0, 8, 16, 24,  3, 27]
        frames  = Int32[1, 1, 1, 1,  2, 2, 2,  3, 3, 3, 3,  4, 4]
        orders  = UInt8[1, 1, 1, 1,  2, 2, 2,  1, 1, 1, 1,  2, 2]
        n = length(scans)
        df = _base(n)
        df.retentionTime = Float32[0.10, 0.101, 0.102, 0.103,  0.20, 0.201, 0.202,  1.10, 1.101, 1.102, 1.103,  1.20, 1.201]
        df.msOrder = orders
        df.centerMz = Vector{Union{Missing,Float32}}([o == 1 ? missing : 612.5f0 for o in orders])
        df.isolationWidthMz = Vector{Union{Missing,Float32}}([o == 1 ? missing : 25f0 for o in orders])
        df.cycle_idx = Int32[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
        df.frameId = frames; df.imScan = scans
        ms = _write_ms(joinpath(d, "packets.arrow"), df)
        m = build_scan_to_ms1(ms)
        # MS2 rows of frame 2 (RT 0.20) -> MS1 frame 1 (RT 0.10 is nearer than 1.10): scan 5 -> slice 8 (row 2),
        # 13 -> 16 (row 3), 21 -> 24 (row 4); frame 4 (RT 1.20) -> MS1 frame 3: scan 3 -> 0 (row 8), 27 -> 24 (row 11)
        @test m[5:7] == Int32[2, 3, 4]
        @test m[12:13] == Int32[8, 11]
        # MS1 rows map to themselves (nearest frame, same scan)
        @test m[1:4] == Int32[1, 2, 3, 4] && m[8:11] == Int32[8, 9, 10, 11]
        # the old nearest-RT rule would have sent every frame-2 row to row 4 (last slice of frame 1)
        @test !(all(m[5:7] .== 4))
    end
end
