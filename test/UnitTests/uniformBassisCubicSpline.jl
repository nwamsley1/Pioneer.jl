# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
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

@testset "UniformSpline Tests" begin
    @testset "Basic sine wave approximation" begin
        # Example data: Approximate a sine curve
        N = 200
        t = collect(LinRange(0.0, 4*π, N))
        u = sin.(t)

        test_spline = UniformSpline(u, t, 3, 20)
        @test maximum(abs.(test_spline.(t) .- u)) < 1e-3

        # Test with minimum knots
        test_spline_min = UniformSpline(u, t, 3, 3)
        @test test_spline_min isa UniformSpline
    end

    @testset "Zero coefficient handling - UniformSpline" begin
        # Test case 1: Nearly constant data (should produce zero higher-order coefficients)
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.0001f0 .* randn(Float32, N)  # Nearly constant

        # Use many knots to force edge polynomials to have zero coefficients
        spline = UniformSpline(u, t, 3, 25)
        @test spline isa UniformSpline

        # Verify spline evaluates correctly at all input points
        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
        @test length(predictions) == N
    end

    @testset "Zero coefficient handling - Perfectly linear data" begin
        # Test case 2: Perfectly linear data (zero curvature)
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 2.0f0 .+ 3.0f0 .* t  # y = 2 + 3x (perfectly linear)

        # Many knots with linear data should produce polynomials with zero cubic terms
        spline = UniformSpline(u, t, 3, 30)
        @test spline isa UniformSpline

        # Verify linear relationship is preserved
        predictions = [spline(ti) for ti in t]
        @test maximum(abs.(predictions .- u)) < 1e-4
    end

    @testset "Zero coefficient handling - Boundary polynomials" begin
        # Test case 3: Data that goes flat at boundaries
        N = 150
        t = collect(LinRange(0.0f0, 10.0f0, N))
        # Flat at edges, curved in middle
        u = [ti < 1.0f0 ? 1.0f0 : (ti > 9.0f0 ? 3.0f0 : 1.0f0 + sin(ti)) for ti in t]

        # High knot count to test boundary handling
        spline = UniformSpline(u, t, 3, 40)
        @test spline isa UniformSpline

        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
    end
end

@testset "UniformSplinePenalized Tests" begin
    @testset "Basic sine wave with penalty" begin
        N = 200
        t = collect(LinRange(0.0, 4*π, N))
        u = sin.(t)

        # Test with low penalty
        spline_low = UniformSplinePenalized(u, t, 3, 20, 0.01, 2)
        @test maximum(abs.(spline_low.(t) .- u)) < 1e-2

        # Test with medium penalty
        spline_med = UniformSplinePenalized(u, t, 3, 20, 0.1, 2)
        @test maximum(abs.(spline_med.(t) .- u)) < 0.1

        # Test with high penalty (should be smoother, less accurate)
        spline_high = UniformSplinePenalized(u, t, 3, 20, 1.0, 2)
        @test maximum(abs.(spline_high.(t) .- u)) < 0.5
    end

    @testset "Zero penalty equals standard spline" begin
        N = 100
        t = collect(LinRange(0.0, 10.0, N))
        u = sin.(t) .+ cos.(t ./ 2)

        spline_standard = UniformSpline(u, t, 3, 10)
        spline_no_penalty = UniformSplinePenalized(u, t, 3, 10, 0.0, 2)

        # Should give nearly identical results
        test_points = LinRange(0.0, 10.0, 50)
        for tp in test_points
            @test abs(spline_standard(tp) - spline_no_penalty(tp)) < 1e-6
        end
    end

    @testset "Zero coefficient handling - Nearly constant data with high penalty" begin
        # This is the case that triggers the bug!
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.0001f0 .* randn(Float32, N)  # Nearly constant

        # High penalty with many knots should produce zero coefficients at boundaries
        spline = UniformSplinePenalized(u, t, 3, 25, 1.0f0, 2)
        @test spline isa UniformSpline

        # Verify evaluation works
        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
        @test length(predictions) == N
    end

    @testset "Zero coefficient handling - Monotonic filtered data" begin
        # Simulate the monotonic enforcement scenario from make_spline_monotonic
        N = 5821  # From actual failure case
        t = collect(LinRange(5.2123957f0, 112.212006f0, N))

        # Create smooth data with some noise
        u_base = LinRange(-0.05f0, 27.19f0, N)
        u = collect(u_base) .+ 0.1f0 .* randn(Float32, N)

        # Apply monotonic filtering (cummax-like smoothing)
        for i in 2:N
            if u[i] < u[i-1]
                u[i] = u[i-1]
            end
        end

        # This combination (smooth data + high penalty + many knots) triggers the bug
        spline = UniformSplinePenalized(u, t, 3, 59, 1.0f0, 2)
        @test spline isa UniformSpline

        # Verify evaluation at boundaries (where zero coefficients are likely)
        @test isfinite(spline(t[1]))
        @test isfinite(spline(t[end]))
        @test isfinite(spline(median(t)))
    end

    @testset "Zero coefficient handling - Different penalty orders" begin
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.01f0 .* randn(Float32, N)

        # Test all supported penalty orders
        for order in [1, 2, 3]
            spline = UniformSplinePenalized(u, t, 3, 20, 0.5f0, order)
            @test spline isa UniformSpline
            @test all(isfinite.([spline(ti) for ti in t]))
        end
    end
end

@testset "fit_spline_with_ransac Tests" begin
    @testset "RANSAC with outliers" begin
        N = 200
        t = collect(LinRange(0.0, 10.0, N))
        u = sin.(t)

        # Add outliers
        u_noisy = copy(u)
        outlier_indices = rand(1:N, 20)
        u_noisy[outlier_indices] .+= randn(20) .* 2.0

        # RANSAC should handle outliers better than standard fitting
        spline_ransac = fit_spline_with_ransac(u_noisy, t, 3, 10, 0.1, 2)
        @test spline_ransac isa UniformSpline

        # Predictions should be reasonable
        predictions = [spline_ransac(ti) for ti in t]
        @test all(isfinite.(predictions))
    end

    @testset "RANSAC zero coefficient handling" begin
        # Test RANSAC with data that produces zero coefficients
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.1f0 .* randn(Float32, N)

        spline = fit_spline_with_ransac(u, t, 3, 25, 1.0f0, 2,
                                        ransac_iterations=20,
                                        ransac_sample_size=30)
        @test spline isa UniformSpline
        @test all(isfinite.([spline(ti) for ti in t]))
    end
end
