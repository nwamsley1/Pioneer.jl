# Unit tests for Whittaker-Henderson smoother (src/utils/ML/wittakerHendersonSmoothing.jl).
#
# Tests the pentadiagonal LDL' smoother used for chromatogram smoothing
# during peak integration. Validates that noisy signals are recovered,
# smoothing strength scales with λ, and edge cases are handled.
#
# Run standalone: julia --project=. test/UnitTests/test_whittaker_henderson.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using Random
    using Statistics
end

using Pioneer: WHWorkspace, whitsmddw!

@testset "Whittaker-Henderson Smoother" begin

# ═══════════════════════════════════════════════════════════════════════════════
# WHWorkspace construction
# ═══════════════════════════════════════════════════════════════════════════════

@testset "WHWorkspace" begin
    ws = WHWorkspace(100)
    @test ws.n_max == 100
    @test length(ws.d0) == 100
    @test length(ws.d1) == 99
    @test length(ws.d2) == 98
    @test length(ws.z) == 100
    @test all(ws.z .== 0.0f0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Smooth constant signal — should return the constant
# ═══════════════════════════════════════════════════════════════════════════════

@testset "constant signal" begin
    n = 50
    ws = WHWorkspace(n)
    x = Float32.(collect(1:n))
    y = fill(5.0f0, n)
    w = ones(Float32, n)

    z = whitsmddw!(ws, x, y, w, n, 1.0f0)

    @test length(z) == n
    @test all(abs.(z .- 5.0f0) .< 0.01f0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Smooth noisy sine — should recover smooth curve
# ═══════════════════════════════════════════════════════════════════════════════

@testset "noisy Gaussian peak recovery" begin
    Random.seed!(42)
    n = 200
    ws = WHWorkspace(n)
    x = Float32.(collect(LinRange(0.0, 100.0, n)))
    y_clean = Float32.(exp.(-(x .- 50.0f0) .^ 2 ./ 200.0f0))
    y_noisy = y_clean .+ 0.05f0 .* abs.(Float32.(randn(n)))
    w = ones(Float32, n)

    z = whitsmddw!(ws, x, y_noisy, w, n, 100.0f0)

    # Smoothed result should be closer to clean signal than noisy input
    error_noisy = mean((y_noisy .- y_clean) .^ 2)
    error_smooth = mean((z .- y_clean) .^ 2)
    @test error_smooth < error_noisy
end

# ═══════════════════════════════════════════════════════════════════════════════
# λ controls smoothness
# ═══════════════════════════════════════════════════════════════════════════════

@testset "higher λ → smoother output" begin
    Random.seed!(123)
    n = 100
    ws = WHWorkspace(n)
    x = Float32.(collect(1:n))
    y = Float32.(sin.(x ./ 10.0f0) .+ 0.5f0 .* randn(Float32, n))
    w = ones(Float32, n)

    z_low = collect(whitsmddw!(ws, x, y, w, n, 1.0f0))
    z_high = collect(whitsmddw!(ws, x, y, w, n, Float32(1e6)))

    # Measure roughness as sum of squared second differences
    roughness(v) = sum(diff(diff(v)) .^ 2)
    @test roughness(z_high) < roughness(z_low)
end

# ═══════════════════════════════════════════════════════════════════════════════
# λ = 0 → no smoothing (returns weighted input)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "zero lambda returns input" begin
    n = 30
    ws = WHWorkspace(n)
    x = Float32.(collect(1:n))
    y = Float32.(collect(1:n))  # linear
    w = ones(Float32, n)

    z = whitsmddw!(ws, x, y, w, n, 0.0f0)

    @test maximum(abs.(z .- y)) < Float32(1e-4)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Weights — zero weight ignores data point
# ═══════════════════════════════════════════════════════════════════════════════

@testset "zero weight ignores point" begin
    n = 50
    ws = WHWorkspace(n)
    x = Float32.(collect(1:n))
    y = fill(1.0f0, n)
    y[25] = 100.0f0  # outlier
    w = ones(Float32, n)
    w[25] = 0.0f0  # ignore the outlier

    z = whitsmddw!(ws, x, y, w, n, 100.0f0)

    # Smoothed value at the outlier should be close to 1.0, not 100.0
    @test z[25] < 5.0f0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Non-negative output (clamped)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "output non-negative" begin
    n = 50
    ws = WHWorkspace(n)
    x = Float32.(collect(1:n))
    y = Float32.(vcat(zeros(25), ones(25)))
    w = ones(Float32, n)

    z = whitsmddw!(ws, x, y, w, n, 10.0f0)

    @test all(z .>= 0.0f0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Non-uniform spacing (chromatogram-like)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "non-uniform spacing" begin
    Random.seed!(99)
    n = 80
    ws = WHWorkspace(n)
    # Non-uniform x: denser in the middle (like real chromatograms)
    x_raw = sort(Float32.(rand(n) .* 100.0))
    y = Float32.(exp.(-(x_raw .- 50.0f0) .^ 2 ./ 200.0f0))  # Gaussian peak
    y_noisy = y .+ 0.05f0 .* Float32.(randn(n))
    w = ones(Float32, n)

    z = whitsmddw!(ws, x_raw, y_noisy, w, n, Float32(1e3))

    @test length(z) == n
    @test all(isfinite.(z))
    # Should still find the peak near x=50
    peak_idx = argmax(collect(z))
    @test 30.0f0 < x_raw[peak_idx] < 70.0f0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Small inputs (minimum viable sizes)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "small inputs" begin
    @testset "n=3 (minimum for 2nd difference)" begin
        ws = WHWorkspace(10)
        x = Float32[1.0, 2.0, 3.0]
        y = Float32[1.0, 2.0, 3.0]
        w = Float32[1.0, 1.0, 1.0]
        z = whitsmddw!(ws, x, y, w, 3, 1.0f0)
        @test length(z) == 3
        @test all(isfinite.(z))
    end

    @testset "n=2" begin
        ws = WHWorkspace(10)
        x = Float32[1.0, 2.0]
        y = Float32[1.0, 2.0]
        w = Float32[1.0, 1.0]
        z = whitsmddw!(ws, x, y, w, 2, 1.0f0)
        @test length(z) == 2
        @test all(isfinite.(z))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Workspace reuse — same workspace for different sized problems
# ═══════════════════════════════════════════════════════════════════════════════

@testset "workspace reuse" begin
    ws = WHWorkspace(100)

    # First call with n=50
    x1 = Float32.(collect(1:50))
    y1 = Float32.(sin.(x1 ./ 5.0f0))
    w1 = ones(Float32, 50)
    z1 = collect(whitsmddw!(ws, x1, y1, w1, 50, 100.0f0))

    # Second call with n=30 (smaller)
    x2 = Float32.(collect(1:30))
    y2 = Float32.(cos.(x2 ./ 5.0f0))
    w2 = ones(Float32, 30)
    z2 = collect(whitsmddw!(ws, x2, y2, w2, 30, 100.0f0))

    # Both should produce valid results
    @test length(z1) == 50
    @test length(z2) == 30
    @test all(isfinite.(z1))
    @test all(isfinite.(z2))
    # Results should be different (different signals)
    @test z1[1:30] != z2
end

end # top-level testset
