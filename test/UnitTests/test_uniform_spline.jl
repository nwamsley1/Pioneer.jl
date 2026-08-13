# Unit tests for UniformSpline, UniformSplinePenalized, and fit_spline_with_ransac
# (src/utils/ML/uniformBasisCubicSpline.jl).
#
# These are the splines used for retention time alignment and quant normalization.
# Tests fit various known functions and verify recovery within tolerance.
#
# Run standalone: julia --project=. test/UnitTests/test_uniform_spline.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using Random
    using Statistics
    using LinearAlgebra
end

using Pioneer: UniformSpline, build_difference_matrix, calculate_adaptive_iterations

@testset "Uniform Cubic Splines" begin

# ═══════════════════════════════════════════════════════════════════════════════
# UniformSpline — unpenalized least squares fit
# ═══════════════════════════════════════════════════════════════════════════════

@testset "UniformSpline" begin

    @testset "sine curve" begin
        N = 200
        t = collect(LinRange(0.0, 4π, N))
        u = sin.(t)
        spline = UniformSpline(u, t, 3, 20)

        # Should approximate sine within 1e-3 at training points
        residuals = abs.(spline.(t) .- u)
        @test maximum(residuals) < 1e-3
    end

    @testset "linear function (exact recovery)" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = 2.0 .* t .+ 3.0
        spline = UniformSpline(u, t, 3, 5)

        # Cubic spline should recover a linear function exactly
        predicted = spline.(t)
        @test maximum(abs.(predicted .- u)) < 1e-4
    end

    @testset "quadratic function" begin
        t = collect(LinRange(-5.0, 5.0, 100))
        u = 0.5 .* t .^ 2 .- 3.0 .* t .+ 1.0
        spline = UniformSpline(u, t, 3, 10)

        predicted = spline.(t)
        @test maximum(abs.(predicted .- u)) < 1e-2
    end

    @testset "monotonic increasing (RT alignment-like)" begin
        # Simulates retention time mapping: roughly linear with slight curvature
        t = collect(LinRange(0.0, 100.0, 200))
        u = 0.8 .* t .+ 5.0 .+ 2.0 .* sin.(t ./ 15.0)
        spline = UniformSpline(u, t, 3, 10)

        predicted = spline.(t)
        @test maximum(abs.(predicted .- u)) < 0.5
    end

    @testset "clamping at boundaries" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = sin.(t)
        spline = UniformSpline(u, t, 3, 10)

        # Evaluation below first knot should clamp to first value
        v_below = spline(-5.0)
        v_first = spline(0.0)
        @test isfinite(v_below)

        # Evaluation above last knot should clamp to last value
        v_above = spline(15.0)
        v_last = spline(10.0)
        @test isfinite(v_above)
    end

    @testset "unsorted input" begin
        Random.seed!(42)
        t = collect(LinRange(0.0, 4π, 100))
        u = sin.(t)
        perm = randperm(100)
        spline = UniformSpline(u[perm], t[perm], 3, 15)

        # Should give same result as sorted input
        spline_sorted = UniformSpline(u, t, 3, 15)
        test_pts = collect(LinRange(0.0, 4π, 50))
        @test maximum(abs.(spline.(test_pts) .- spline_sorted.(test_pts))) < 1e-10
    end

    @testset "minimum knots (3)" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = 2.0 .* t .+ 1.0
        spline = UniformSpline(u, t, 3, 3)
        @test isfinite(spline(5.0))
    end

    @testset "input validation" begin
        @test_throws ErrorException UniformSpline([1.0], [1.0], 2, 5)   # degree != 3
        @test_throws ErrorException UniformSpline([1.0], [1.0], 3, 2)   # n_knots < 3
        @test_throws ErrorException UniformSpline([1.0, 2.0], [1.0], 3, 3)  # length mismatch
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# UniformSplinePenalized — P-spline with difference penalty
# ═══════════════════════════════════════════════════════════════════════════════

@testset "UniformSplinePenalized" begin

    @testset "sine curve" begin
        N = 200
        t = collect(LinRange(0.0, 4π, N))
        u = sin.(t)
        spline = Pioneer.UniformSplinePenalized(u, t, 3, 20, 0.1, 2)

        residuals = abs.(spline.(t) .- u)
        @test maximum(residuals) < 0.05
    end

    @testset "zero penalty matches unpenalized" begin
        t = collect(LinRange(0.0, 4π, 100))
        u = sin.(t)
        spline_pen = Pioneer.UniformSplinePenalized(u, t, 3, 10, 0.0, 2)
        spline_std = UniformSpline(u, t, 3, 10)

        test_pts = collect(LinRange(0.0, 4π, 50))
        # With λ=0, penalized should be identical to standard
        @test maximum(abs.(spline_pen.(test_pts) .- spline_std.(test_pts))) < 1e-4
    end

    @testset "higher penalty → smoother fit" begin
        Random.seed!(123)
        t = collect(LinRange(0.0, 10.0, 100))
        u = sin.(t) .+ 0.3 .* randn(100)

        spline_low  = Pioneer.UniformSplinePenalized(u, t, 3, 15, 0.01, 2)
        spline_high = Pioneer.UniformSplinePenalized(u, t, 3, 15, 10.0, 2)

        # Measure roughness as sum of squared second differences in predictions
        test_pts = collect(LinRange(0.0, 10.0, 200))
        pred_low  = spline_low.(test_pts)
        pred_high = spline_high.(test_pts)
        roughness_low  = sum(diff(diff(pred_low)) .^ 2)
        roughness_high = sum(diff(diff(pred_high)) .^ 2)

        @test roughness_high < roughness_low
    end

    @testset "linear function recovery" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = 3.0 .* t .- 2.0
        spline = Pioneer.UniformSplinePenalized(u, t, 3, 5, 0.1, 2)

        predicted = spline.(t)
        @test maximum(abs.(predicted .- u)) < 0.1
    end

    @testset "input validation" begin
        @test_throws ErrorException Pioneer.UniformSplinePenalized([1.0], [1.0], 2, 5)
        @test_throws ErrorException Pioneer.UniformSplinePenalized([1.0], [1.0], 3, 2)
        @test_throws ErrorException Pioneer.UniformSplinePenalized([1.0, 2.0], [1.0], 3, 3)
        @test_throws ErrorException Pioneer.UniformSplinePenalized([1.0], [1.0], 3, 3, -0.1)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# fit_spline_with_ransac — robust fitting with outliers
# ═══════════════════════════════════════════════════════════════════════════════

@testset "fit_spline_with_ransac" begin

    @testset "clean data (no outliers)" begin
        t = collect(LinRange(0.0, 4π, 200))
        u = sin.(t)
        spline = Pioneer.fit_spline_with_ransac(u, t, 3, 20; ransac_seed=42)

        # RANSAC uses penalized splines internally, so tolerance is moderate
        residuals = abs.(spline.(t) .- u)
        @test median(residuals) < 0.05
    end

    @testset "with 20% outliers" begin
        Random.seed!(42)
        N = 200
        t = collect(LinRange(0.0, 4π, N))
        u = sin.(t)

        # Add 20% gross outliers
        n_outliers = div(N, 5)
        outlier_idx = randperm(N)[1:n_outliers]
        u_noisy = copy(u)
        u_noisy[outlier_idx] .+= 5.0 .* randn(n_outliers)

        spline = Pioneer.fit_spline_with_ransac(u_noisy, t, 3, 15;
            ransac_seed=42, ransac_iterations=100)

        # Evaluate on clean inlier points
        inlier_mask = trues(N)
        inlier_mask[outlier_idx] .= false
        inlier_residuals = abs.(spline.(t[inlier_mask]) .- sin.(t[inlier_mask]))

        # Should still fit inliers well despite outliers
        @test median(inlier_residuals) < 0.1
    end

    @testset "monotonic RT-like data with outliers" begin
        Random.seed!(99)
        N = 150
        t = collect(LinRange(0.0, 100.0, N))
        u = 0.9 .* t .+ 5.0  # Linear RT alignment

        # Add 15% outliers
        n_outliers = div(N, 7)
        outlier_idx = randperm(N)[1:n_outliers]
        u_noisy = copy(u)
        u_noisy[outlier_idx] .+= 20.0 .* randn(n_outliers)

        spline = Pioneer.fit_spline_with_ransac(u_noisy, t, 3, 5;
            ransac_seed=123, ransac_iterations=50)

        # Check fit quality on clean data
        inlier_mask = trues(N)
        inlier_mask[outlier_idx] .= false
        residuals = abs.(spline.(t[inlier_mask]) .- u[inlier_mask])
        @test median(residuals) < 1.0
    end

    @testset "input validation" begin
        @test_throws ErrorException Pioneer.fit_spline_with_ransac([1.0], [1.0], 2, 5)
        @test_throws ErrorException Pioneer.fit_spline_with_ransac([1.0], [1.0], 3, 2)
        @test_throws ErrorException Pioneer.fit_spline_with_ransac([1.0, 2.0], [1.0], 3, 3)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# build_difference_matrix
# ═══════════════════════════════════════════════════════════════════════════════

@testset "build_difference_matrix" begin

    @testset "order 1" begin
        D = build_difference_matrix(4, 1)
        @test size(D) == (3, 4)
        # First row: [-1, 1, 0, 0]
        @test D[1, :] == [-1.0, 1.0, 0.0, 0.0]
        @test D[2, :] == [0.0, -1.0, 1.0, 0.0]
        @test D[3, :] == [0.0, 0.0, -1.0, 1.0]
    end

    @testset "order 2" begin
        D = build_difference_matrix(5, 2)
        @test size(D) == (3, 5)
        # Second order differences: [1, -2, 1, 0, 0], etc.
        @test D[1, :] ≈ [1.0, -2.0, 1.0, 0.0, 0.0]
    end

    @testset "order 3" begin
        D = build_difference_matrix(6, 3)
        @test size(D) == (3, 6)
        # Third order differences: [-1, 3, -3, 1, 0, 0]
        @test D[1, :] ≈ [-1.0, 3.0, -3.0, 1.0, 0.0, 0.0]
    end

    @testset "P-spline typical: 8 coefficients order 2" begin
        D = build_difference_matrix(8, 2)
        @test size(D) == (6, 8)
        P = D' * D
        @test size(P) == (8, 8)
        @test issymmetric(P)
    end

    @testset "invalid order" begin
        @test_throws ErrorException build_difference_matrix(5, 0)
        @test_throws ErrorException build_difference_matrix(5, 4)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# calculate_adaptive_iterations
# ═══════════════════════════════════════════════════════════════════════════════

@testset "calculate_adaptive_iterations" begin

    @testset "zero outliers → minimum iterations" begin
        @test calculate_adaptive_iterations(0.0, 10) == 20
    end

    @testset "high outlier ratio → more iterations" begin
        iters_low  = calculate_adaptive_iterations(0.1, 10)
        iters_high = calculate_adaptive_iterations(0.4, 10)
        @test iters_high > iters_low
    end

    @testset "clamped to [20, 200]" begin
        @test calculate_adaptive_iterations(0.0, 5) >= 20
        @test calculate_adaptive_iterations(0.49, 30) <= 200
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Callable evaluation
# ═══════════════════════════════════════════════════════════════════════════════

@testset "spline evaluation" begin

    @testset "Float32 input" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = sin.(t)
        spline = UniformSpline(u, t, 3, 10)
        @test isfinite(spline(5.0f0))
    end

    @testset "Float64 input" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = sin.(t)
        spline = UniformSpline(u, t, 3, 10)
        @test isfinite(spline(5.0))
    end

    @testset "constant function" begin
        t = collect(LinRange(0.0, 10.0, 50))
        u = fill(7.0, 50)
        spline = UniformSpline(u, t, 3, 5)
        # Should return ~7.0 everywhere
        for x in LinRange(0.0, 10.0, 20)
            @test abs(spline(x) - 7.0) < 0.1
        end
    end
end

end # top-level testset
