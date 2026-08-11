@testset "Quant normalization RT binning" begin

    ob = Pioneer._occupancy_bins

    @testset "_occupancy_bins invariants" begin
        # Every bin holds at least min_occupancy, sizes differ by at most one, and
        # the bins tile 1:n exactly. The old scheme failed the last of these:
        # bin_stop = bin_idx * bin_size dropped the final n mod max_bins rows.
        for (n, m, mx) in [(1000, 100, 100), (8705, 100, 100), (250, 100, 100),
                           (10_000, 100, 100), (999, 100, 100), (12_345, 100, 50)]
            bins = ob(n, m, mx)
            @test !isempty(bins)
            @test all(length(b) >= m for b in bins)
            @test maximum(length, bins) - minimum(length, bins) <= 1
            @test first(first(bins)) == 1
            @test last(last(bins)) == n
            @test sum(length, bins) == n
            for i in 2:length(bins)          # contiguous, no gaps or overlaps
                @test first(bins[i]) == last(bins[i-1]) + 1
            end
            @test length(bins) <= mx
        end
    end

    @testset "bin count is occupancy-limited" begin
        # 8,705 PSMs at >=100 each gives 87 bins, not the 100 the old code used.
        @test length(ob(8705, 100, 100)) == 87
        # Plenty of data: capped by max_bins instead.
        @test length(ob(50_000, 100, 100)) == 100
        # Exactly enough for k bins.
        @test length(ob(1000, 100, 100)) == 10
        @test length(ob(1099, 100, 100)) == 10
        @test length(ob(1100, 100, 100)) == 11
    end

    @testset "degenerate inputs collapse rather than fragment" begin
        # Between m and 2m there is only enough for one bin.
        @test ob(100, 100, 100) == [1:100]
        @test ob(199, 100, 100) == [1:199]
        @test ob(200, 100, 100) == [1:100, 101:200]

        # Below m there is not even one bin: the caller must skip the file.
        @test isempty(ob(99, 100, 100))
        @test isempty(ob(14, 100, 100))     # the yeast lowsignal file
        @test isempty(ob(1, 100, 100))
        @test isempty(ob(0, 100, 100))

        # Defensive: nonsense parameters return empty rather than throwing.
        @test isempty(ob(1000, 0, 100))
        @test isempty(ob(1000, 100, 0))
        @test isempty(ob(-5, 100, 100))
    end

    @testset "the skip threshold excludes the square-system crash" begin
        # UniformSpline's design matrix has n_knots+3 columns, the last of which is
        # structurally always zero. At exactly that many bins the system is square,
        # `\` dispatches to lu, and it throws. _min_bins_for_spline must therefore
        # be strictly greater than the coefficient count.
        n_knots = 7
        n_coeffs = n_knots + 3
        @test Pioneer._min_bins_for_spline(n_knots) > n_coeffs
        @test Pioneer._min_bins_for_spline(n_knots) == n_coeffs + 1

        # The crash this pins: reproduce it directly so a future change to the
        # basis that removes the empty column makes this test fail loudly rather
        # than leaving a threshold nobody can justify.
        t_sq = collect(range(0.0, 1.0, length = n_coeffs))
        u_sq = collect(range(1.0, 2.0, length = n_coeffs))
        @test_throws LinearAlgebra.SingularException Pioneer.UniformSpline(u_sq, t_sq, 3, n_knots)

        # One bin either side is fine, which is why the threshold is +1 and not
        # something larger.
        for k in (n_coeffs - 1, n_coeffs + 1)
            t = collect(range(0.0, 1.0, length = k))
            u = collect(range(1.0, 2.0, length = k))
            @test Pioneer.UniformSpline(u, t, 3, n_knots) isa Pioneer.UniformSpline
        end
    end

    @testset "occupancy floor keeps fits off the crash point" begin
        # With min_occupancy = 100 the bin count is n ÷ 100, so a file lands on
        # exactly n_coeffs bins only for n in 1000:1099 -- and those are skipped,
        # because 10 < _min_bins_for_spline(7) = 11.
        n_coeffs = 10
        min_bins = Pioneer._min_bins_for_spline(7)
        for n in 1000:25:1099
            @test length(ob(n, 100, 100)) == n_coeffs
            @test length(ob(n, 100, 100)) < min_bins      # therefore skipped
        end
        # And no n at all produces a bin count that is both >= min_bins and equal
        # to n_coeffs, which is the property that makes the crash unreachable.
        for n in 100:50:20_000
            k = length(ob(n, 100, 100))
            @test !(k == n_coeffs && k >= min_bins)
        end
    end
end
