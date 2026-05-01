# Unit tests for retentionTimeIndex and buildRtIndex.
#
# Tests cover:
# - Struct construction and accessors
# - buildRtIndex binning logic (single bin, multiple bins, edge cases)
# - Precursor m/z sorting within bins
# - Empty input handling
# - DataFrame convenience method
#
# Run standalone: julia --project=. test/UnitTests/test_retention_time_index.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

@testset "retentionTimeIndex" begin

    # ── Struct basics ────────────────────────────────────────────────────
    @testset "rtIndexBin construction and accessors" begin
        precs = [(UInt32(1), 500.0f0), (UInt32(2), 600.0f0)]
        bin = Pioneer.rtIndexBin(10.0f0, 20.0f0, precs)

        @test Pioneer.getLow(bin) == 10.0f0
        @test Pioneer.getHigh(bin) == 20.0f0
        @test length(bin.prec) == 2
        @test bin.prec[1] == (UInt32(1), 500.0f0)
        @test bin.prec[2] == (UInt32(2), 600.0f0)
    end

    @testset "retentionTimeIndex empty construction" begin
        idx = Pioneer.retentionTimeIndex(Float32, Float32)
        @test isempty(idx.rt_bins)
    end

    # ── buildRtIndex ─────────────────────────────────────────────────────
    @testset "buildRtIndex: empty input" begin
        idx = Pioneer.buildRtIndex(Float32[], Float32[], UInt32[], 0.1)
        @test isempty(idx.rt_bins)
    end

    @testset "buildRtIndex: single precursor" begin
        rts = Float32[5.0]
        mzs = Float32[500.0]
        ids = UInt32[1]
        idx = Pioneer.buildRtIndex(rts, mzs, ids, 0.1)

        @test length(idx.rt_bins) == 1
        @test Pioneer.getLow(idx.rt_bins[1]) == 5.0f0
        @test Pioneer.getHigh(idx.rt_bins[1]) == 5.0f0
        @test idx.rt_bins[1].prec[1] == (UInt32(1), 500.0f0)
    end

    @testset "buildRtIndex: all precursors in one bin" begin
        # 5 precursors spanning 0.05 RT — smaller than bin_size=0.1
        rts = Float32[1.00, 1.01, 1.02, 1.03, 1.04]
        mzs = Float32[800.0, 400.0, 600.0, 200.0, 1000.0]
        ids = UInt32[10, 20, 30, 40, 50]
        idx = Pioneer.buildRtIndex(rts, mzs, ids, 0.1)

        @test length(idx.rt_bins) == 1
        @test Pioneer.getLow(idx.rt_bins[1]) == 1.00f0
        @test Pioneer.getHigh(idx.rt_bins[1]) == 1.04f0
        # Precursors should be sorted by m/z within the bin
        bin_mzs = [last(p) for p in idx.rt_bins[1].prec]
        @test issorted(bin_mzs)
        # Verify all precursors present
        bin_ids = Set([first(p) for p in idx.rt_bins[1].prec])
        @test bin_ids == Set(UInt32[10, 20, 30, 40, 50])
    end

    @testset "buildRtIndex: multiple bins" begin
        # 6 precursors spanning 3.0 RT with bin_size=1.0 → should get ~3 bins
        rts = Float32[1.0, 1.3, 1.7, 2.5, 3.0, 3.8]
        mzs = Float32[500.0, 300.0, 700.0, 400.0, 600.0, 200.0]
        ids = UInt32[1, 2, 3, 4, 5, 6]
        idx = Pioneer.buildRtIndex(rts, mzs, ids, 1.0)

        @test length(idx.rt_bins) >= 2  # at least 2 bins

        # Each bin should have precursors sorted by m/z
        for bin in idx.rt_bins
            bin_mzs = [last(p) for p in bin.prec]
            @test issorted(bin_mzs)
        end

        # All precursors should be present across all bins
        all_ids = Set{UInt32}()
        for bin in idx.rt_bins
            for p in bin.prec
                push!(all_ids, first(p))
            end
        end
        @test all_ids == Set(UInt32[1, 2, 3, 4, 5, 6])

        # Bins should be in RT order
        bin_lows = [Pioneer.getLow(b) for b in idx.rt_bins]
        @test issorted(bin_lows)
    end

    @testset "buildRtIndex: bin boundaries are correct" begin
        # Precise test: 4 precursors at RT 0.0, 0.05, 0.15, 0.20 with bin_size=0.1
        # Expected: bin1=[0.0, 0.05], bin2=[0.15, 0.20]
        rts = Float32[0.0, 0.05, 0.15, 0.20]
        mzs = Float32[100.0, 200.0, 300.0, 400.0]
        ids = UInt32[1, 2, 3, 4]
        idx = Pioneer.buildRtIndex(rts, mzs, ids, 0.1)

        @test length(idx.rt_bins) == 2
        @test Pioneer.getLow(idx.rt_bins[1]) == 0.0f0
        @test Pioneer.getHigh(idx.rt_bins[1]) == 0.05f0
        @test length(idx.rt_bins[1].prec) == 2
        @test Pioneer.getLow(idx.rt_bins[2]) == 0.15f0
        @test Pioneer.getHigh(idx.rt_bins[2]) == 0.20f0
        @test length(idx.rt_bins[2].prec) == 2
    end

    @testset "buildRtIndex: DataFrame convenience" begin
        df = DataFrame(
            irt = Float32[1.0, 2.0, 3.0],
            prec_mz = Float32[500.0, 600.0, 700.0],
            precursor_idx = UInt32[10, 20, 30]
        )
        idx = Pioneer.buildRtIndex(df; bin_rt_size = 5.0)

        # All in one bin since span (2.0) < bin_size (5.0)
        @test length(idx.rt_bins) == 1
        @test length(idx.rt_bins[1].prec) == 3
    end

    @testset "buildRtIndex: two identical RTs" begin
        rts = Float32[5.0, 5.0]
        mzs = Float32[300.0, 600.0]
        ids = UInt32[1, 2]
        idx = Pioneer.buildRtIndex(rts, mzs, ids, 0.1)

        @test length(idx.rt_bins) == 1
        # Both precursors in same bin, sorted by m/z
        @test idx.rt_bins[1].prec[1] == (UInt32(1), 300.0f0)
        @test idx.rt_bins[1].prec[2] == (UInt32(2), 600.0f0)
    end
end
