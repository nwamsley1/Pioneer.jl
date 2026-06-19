using Test
using Pioneer

@testset "library search candidate filters" begin
    @testset "scan candidate filter rewrites per-scan ranges" begin
        scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(missing, 6)
        scan_to_prec_idx[2] = 1:3
        scan_to_prec_idx[5] = 4:5
        scan_to_prec_idx[6] = 6:6
        precursors_passed = UInt32[10, 20, 30, 40, 50, 60]
        center_mzs = Float32[0, 500, 0, 0, 600, 700]
        isolation_widths = fill(4.0f0, 6)

        filtered = Pioneer._filter_scan_candidates_by_window!(
            scan_to_prec_idx,
            precursors_passed,
            Int[2, 5, 6],
            center_mzs,
            isolation_widths,
            60,
            (center_mz, isolation_width) -> center_mz,
            (pid, center_mz) -> (center_mz == 500.0f0 && pid != UInt32(20)) ||
                                (center_mz == 600.0f0 && pid == UInt32(50)),
        )

        @test filtered == UInt32[10, 30, 50]
        @test scan_to_prec_idx[2] == 1:2
        @test scan_to_prec_idx[5] == 3:3
        @test ismissing(scan_to_prec_idx[6])
        @test ismissing(scan_to_prec_idx[1])
    end

    @testset "window candidate filter caches decisions per precursor window" begin
        scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(missing, 5)
        scan_to_prec_idx[2] = 1:2
        scan_to_prec_idx[3] = 3:4
        scan_to_prec_idx[5] = 5:6
        precursors_passed = UInt32[10, 20, 10, 20, 10, 20]
        center_mzs = Float32[0, 500, 500, 0, 600]
        isolation_widths = fill(4.0f0, 5)
        calls_500 = [Threads.Atomic{Int}(0) for _ in 1:20]
        calls_600 = [Threads.Atomic{Int}(0) for _ in 1:20]

        filtered = Pioneer._filter_scan_candidates_by_window!(
            scan_to_prec_idx,
            precursors_passed,
            Int[2, 3, 5],
            center_mzs,
            isolation_widths,
            20,
            (center_mz, isolation_width) -> center_mz,
            function (pid, center_mz)
                calls = center_mz == 500.0f0 ? calls_500 : calls_600
                Threads.atomic_add!(calls[Int(pid)], 1)
                return !(center_mz == 500.0f0 && pid == UInt32(20))
            end,
        )

        @test filtered == UInt32[10, 10, 10, 20]
        @test scan_to_prec_idx[2] == 1:1
        @test scan_to_prec_idx[3] == 2:2
        @test scan_to_prec_idx[5] == 3:4
        @test calls_500[10][] == 1
        @test calls_500[20][] == 1
        @test calls_600[10][] == 1
        @test calls_600[20][] == 1
    end
end
