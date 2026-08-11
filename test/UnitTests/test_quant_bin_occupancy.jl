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

    @testset "coefficient count and the bin threshold" begin
        # n_knots knots bound n_knots-1 spans; a cubic basis over those spans has
        # (n_knots-1)+3 functions. The design matrix used to allocate one more,
        # which could never be filled -- that column is what made every fit rank
        # deficient and made the square case throw.
        for nk in 4:12
            @test Pioneer._n_spline_coeffs(nk) == nk + 2
            @test Pioneer._min_bins_for_spline(nk) == Pioneer._n_spline_coeffs(nk) + 1
        end

        # The design matrix is now full rank with no empty column, at every knot
        # count. This is the property the fit depends on; if it regresses, the
        # square case starts throwing again.
        for nk in 4:12
            t = collect(range(0.0, 1.0, length = 60))
            knots = collect(LinRange(0.0, 1.0, nk))
            bw = 1.0 / (nk - 1)
            X = Pioneer._build_numeric_design_matrix(t, knots, bw)
            @test size(X, 2) == Pioneer._n_spline_coeffs(nk)
            @test rank(X) == size(X, 2)
            @test !any(j -> all(iszero, @view X[:, j]), 1:size(X, 2))
        end
    end

    @testset "square and near-square systems no longer throw" begin
        # Before the _n_spline_coeffs fix, exactly n_coeffs points gave a square
        # matrix with a zero column, `\` dispatched to lu, and it raised
        # SingularException -- the failure that started this work. Fitting at,
        # below and above that count must all succeed now.
        for nk in 4:12
            nc = Pioneer._n_spline_coeffs(nk)
            for k in (nc - 1, nc, nc + 1)
                k < 4 && continue
                t = collect(range(0.0, 1.0, length = k))
                u = collect(range(1.0, 2.0, length = k))
                @test Pioneer.UniformSpline(u, t, 3, nk) isa Pioneer.UniformSpline
            end
        end
    end

    @testset "fitted splines stay finite over the fitted domain" begin
        # Guards against a rank fix that silently produces NaN/Inf coefficients on
        # awkward inputs rather than throwing.
        rng_t = [sort(rand(n)) .* 25.0 for n in (15, 40, 120, 400)]
        for t in rng_t, nk in (4, 7, 10)
            length(t) <= Pioneer._n_spline_coeffs(nk) && continue
            u = [sin(x) * 3 + 10 for x in t]
            s = Pioneer.UniformSpline(u, t, 3, nk)
            grid = collect(range(minimum(t), maximum(t), length = 51))
            @test all(isfinite, [s(g) for g in grid])
        end
    end

    @testset "occupancy floor keeps fits over-determined" begin
        # With min_occupancy = 100 the bin count is n ÷ 100. A file landing on
        # exactly n_coeffs bins would interpolate rather than smooth, so the floor
        # excludes it: n_coeffs = 9 bins needs 900-999 PSMs, below min_bins = 10.
        n_coeffs = Pioneer._n_spline_coeffs(7)
        min_bins = Pioneer._min_bins_for_spline(7)
        for n in 900:25:999
            @test length(ob(n, 100, 100)) == n_coeffs
            @test length(ob(n, 100, 100)) < min_bins      # therefore skipped
        end
        # No PSM count yields a bin count that is both accepted and exactly the
        # coefficient count, so every accepted fit is over-determined.
        for n in 100:50:20_000
            k = length(ob(n, 100, 100))
            @test !(k == n_coeffs && k >= min_bins)
        end
    end
end
