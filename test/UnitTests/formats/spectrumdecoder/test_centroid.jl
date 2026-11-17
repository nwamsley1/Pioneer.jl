# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Profile centroider (src/vendor/SpectrumDecoder/centroid.jl), on synthetic profiles.

gauss(mz, x, fwhm, h) = h * exp(-(mz - x)^2 / (2 * (fwhm / 2.35482)^2))
grid(lo, hi, step) = collect(lo:step:hi)
ppm(a, b) = (a - b) / b * 1e6

@testset "centroid" begin
    @testset "exact Gaussian recovered: ion-trap and Orbitrap-like sampling" begin
        for (step, fw, x) in ((0.125, 0.5, 500.0371), (0.0025, 0.01, 800.01234), (0.0005, 0.002, 400.2468))
            mz = grid(x - 2, x + 2, step)
            y = Float32[gauss(m, x, fw, 1e5) + 1 for m in mz]
            c = centroid!(Centroid[], mz, y; threshold = 10)
            @test length(c) == 1
            @test abs(ppm(c[1].mz, x)) < 0.01
            @test isapprox(c[1].fwhm, fw; rtol = 1e-3)
            @test isapprox(c[1].height, 1e5; rtol = 1e-3)
            @test isapprox(c[1].area, 1e5 * fw / 2.35482 * sqrt(2pi); rtol = 1e-3)
            @test c[1].resid < 1e-3
        end
    end

    @testset "apex between grid points (sub-sample position)" begin
        for frac in (0.0, 0.1, 0.25, 0.5, 0.73, 0.99)
            x = 600.0 + frac * 0.125
            mz = grid(598.0, 602.0, 0.125); y = Float32[gauss(m, x, 0.5, 5e4) for m in mz]
            c = centroid!(Centroid[], mz, y; threshold = 1)
            @test length(c) == 1 && abs(ppm(c[1].mz, x)) < 0.05
        end
    end

    @testset "several separated peaks, all found, ordered" begin
        xs = [300.1, 310.55, 450.0, 999.9]; hs = [1e3, 5e5, 2e4, 7e4]
        mz = grid(299.0, 1001.0, 0.125)
        y = Float32[sum(gauss(m, x, 0.5, h) for (x, h) in zip(xs, hs)) for m in mz]
        c = centroid!(Centroid[], mz, y; threshold = 10)
        @test length(c) == 4
        @test all(abs(ppm(c[k].mz, xs[k])) < 0.5 for k in 1:4)
        @test issorted([p.mz for p in c])
    end

    @testset "noise: m/z stays within a fraction of the sampling step" begin
        x = 700.4321; mz = grid(698.0, 703.0, 0.125)
        for seed in 1:20
            r = Float32[1 + 0.02 * (2 * ((seed * 7919 + j * 104729) % 1000) / 1000 - 1) for j in eachindex(mz)]   # deterministic ±2 %
            y = Float32[gauss(m, x, 0.5, 1e5) for m in mz] .* r .+ 5
            c = centroid!(Centroid[], mz, y; threshold = 1e3)
            @test length(c) == 1 && abs(c[1].mz - x) < 0.01
        end
    end

    @testset "threshold gates small maxima" begin
        mz = grid(499.0, 503.0, 0.125)
        y = Float32[gauss(m, 500.0, 0.5, 1e5) + gauss(m, 502.0, 0.5, 50) for m in mz]
        @test length(centroid!(Centroid[], mz, y; threshold = 10)) == 2
        @test length(centroid!(Centroid[], mz, y; threshold = 100)) == 1
        # default gate: snr × median of positive intensities
        @test length(centroid!(Centroid[], mz, y)) >= 1
    end

    @testset "reduced profile: separate runs (gaps) are never joined" begin
        # two runs of points with a gap wider than gap_factor × spacing; each run is one peak
        mz1 = grid(500.0, 500.1, 0.0025); mz2 = grid(500.3, 500.4, 0.0025)
        mz = vcat(mz1, mz2)
        y = Float32[gauss(m, 500.05, 0.01, 1e5) + gauss(m, 500.35, 0.01, 2e5) for m in mz]
        c = centroid!(Centroid[], mz, y; threshold = 1)
        @test length(c) == 2
        @test abs(ppm(c[1].mz, 500.05)) < 0.1 && abs(ppm(c[2].mz, 500.35)) < 0.1
    end

    @testset "edge cases" begin
        @test isempty(centroid!(Centroid[], Float64[], Float32[]))
        @test isempty(centroid!(Centroid[], [1.0, 2.0], Float32[1, 2]))                # fewer than 3 points
        @test isempty(centroid!(Centroid[], grid(1.0, 2.0, 0.1), fill(5f0, 11); threshold = 0))   # flat: no strict maximum
        @test isempty(centroid!(Centroid[], grid(1.0, 2.0, 0.1), zeros(Float32, 11); threshold = 0))
        # a peak at the very first / last point is not a candidate (needs a neighbour on both sides)
        mz = grid(500.0, 501.0, 0.125); y = Float32[gauss(m, 500.0, 0.5, 1e5) for m in mz]
        @test isempty(centroid!(Centroid[], mz, y; threshold = 1))
        # a single-point spike: falls back to the apex point (no fit)
        y = zeros(Float32, 9); y[5] = 100; mz = grid(1.0, 2.0, 0.125)
        c = centroid!(Centroid[], mz, y; threshold = 1)
        @test length(c) == 1 && c[1].mz == mz[5] && c[1].height == 100f0 && isnan(c[1].area)
        # output buffer is reused, not appended to
        out = Centroid[]; centroid!(out, mz, y; threshold = 1); centroid!(out, mz, y; threshold = 1)
        @test length(out) == 1
    end

    @testset "non-peak-shaped points fall back to the apex point" begin
        mz = [1.0, 1.1, 1.2, 1.3, 1.4]; y = Float32[1, 100, 101, 100, 1]   # plateau: parabola fit is degenerate-ish
        c = centroid!(Centroid[], mz, y; threshold = 0)
        @test length(c) == 1 && 1.1 <= c[1].mz <= 1.3 && c[1].height >= 100
    end

    @testset "median gate helper" begin
        @test _median_positive(Float32[0, 0, 3, 1, 2]) == 2f0
        @test _median_positive(Float32[0, 0]) == 0f0
        @test _median_positive!(Float32[], Float32[5, 1, 9, 3, 7, 2, 8, 4]; stride = 2) == 7f0   # points 1,3,5,7: 5, 9, 7, 8 -> lower median 7
    end

    @testset "allocation-free inner loop when the output buffer is warm" begin
        mz = grid(300.0, 1100.0, 0.125); y = Float32[sum(gauss(m, x, 0.5, 1e4) for x in 310:7.3:1090) + 3 for m in mz]
        out = Centroid[]; centroid!(out, mz, y; threshold = 10)
        @test (@allocated centroid!(out, mz, y; threshold = 10)) == 0
    end
end
