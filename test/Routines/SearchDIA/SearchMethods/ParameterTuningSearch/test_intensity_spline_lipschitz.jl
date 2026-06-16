# Regression test for the Lipschitz step cap in
# fit_intensity_mass_error.jl. The projected gradient descent loops
# inside fit_convex_bias_spline / fit_monotone_bias_spline /
# constrained spline fitting used a fixed lr=1e-4. opnorm(H)
# scales with the row count of the design matrix, so once n_pts pushes
# L = opnorm(H) past 1/lr the optimizer diverges, IRLS reweights to
# non-finite values, and the next H \ Xty throws
# `ArgumentError: matrix contains Infs or NaNs` from LAPACK.
#
# These tests build large-N synthetic inputs (300k points) that would
# trip the original divergence and assert the fit returns a finite
# UniformSpline.

using Random
using Statistics: median
using Pioneer: fit_convex_bias_spline, fit_monotone_bias_spline,
               fit_regularized_spread_spline,
               fit_intensity_mass_error_model,
               IntensityMassErrorModel,
               MassErrSample,
               SimpleMassErrorModel,
               getCorrectedMz,
               _rt_bias_da

@testset "Intensity spline Lipschitz step cap" begin
    rng = MersenneTwister(0x5eedfa11)
    n = 300_000

    @testset "fit_convex_bias_spline (n=$n)" begin
        x = collect(range(100.0, 1500.0; length=n))
        y = 1e-3 .* (x ./ 1000.0 .- 0.7).^2 .+ 5e-4 .* randn(rng, n)
        spline = fit_convex_bias_spline(x, y; convex_up=true,
                                         n_knots=10, λ=1.0,
                                         max_iter=1000, lr=1e-4)
        @test spline !== nothing
        @test all(isfinite, spline.coeffs)
    end

    @testset "fit_monotone_bias_spline (n=$n)" begin
        x = collect(range(0.0, 30.0; length=n))
        y = 1e-3 .* (x .- 15.0) .+ 5e-4 .* randn(rng, n)
        spline = fit_monotone_bias_spline(x, y; increasing=true,
                                           n_knots=10, λ=1.0,
                                           max_iter=1000, lr=1e-4)
        @test spline !== nothing
        @test all(isfinite, spline.coeffs)
    end

    @testset "fit_regularized_spread_spline" begin
        # The spread fit consumes pre-binned (centers, scales) so its design
        # matrix is small. This guards the same shared spline-system path.
        centers = collect(range(10.0, 30.0; length=200))
        scales = max.(0.5 .+ 5.0 .* exp.(-(centers .- 10.0) ./ 5.0)
                       .+ 0.05 .* randn(rng, length(centers)), 0.05)
        spline = fit_regularized_spread_spline(centers, scales;
                                                n_knots=10, λ=0.01)
        @test spline !== nothing
        @test all(isfinite, spline.coeffs)
    end
end

@testset "Intensity mass error model learns RT-dependent bias" begin
    rng = MersenneTwister(0x71cafe)
    n = 1_800
    rts = collect(range(5.0, 85.0; length=n))
    mzs = [420.0 + 7.0 * mod(i, 120) + 0.03 * randn(rng) for i in 1:n]
    log2I = [10.0 + 0.05 * mod(37 * i, 220) for i in 1:n]
    intensities = Float32.(2.0 .^ log2I)

    rt_bias = @. 0.0014 * tanh((rts - 48.0) / 7.0)
    noise = 0.00003 .* randn(rng, n)
    da_err = rt_bias .+ noise

    samples = [
        MassErrSample(Float32(mzs[i]), Float32(mzs[i] + da_err[i]), intensities[i], Float32(rts[i]))
        for i in 1:n
    ]

    model = fit_intensity_mass_error_model(
        samples,
        SimpleMassErrorModel(0.0f0, (30.0f0, 30.0f0)),
    )

    @test model isa IntensityMassErrorModel
    @test model.rt_min ≈ Float32(first(rts))
    @test model.rt_max ≈ Float32(last(rts))

    early_bias = _rt_bias_da(model, 10.0f0)
    late_bias = _rt_bias_da(model, 80.0f0)
    @test late_bias - early_bias > 0.0015f0

    obs_mz = 700.0f0
    intensity = Float32(2.0 ^ 16.0)
    @test getCorrectedMz(model, obs_mz, intensity, 80.0f0) <
          getCorrectedMz(model, obs_mz, intensity, 10.0f0)
end
