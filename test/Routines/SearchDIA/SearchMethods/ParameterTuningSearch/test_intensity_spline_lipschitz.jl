# Regression test for the Lipschitz step cap in
# fit_intensity_mass_error.jl. The projected gradient descent loops
# inside fit_convex_bias_spline / fit_monotone_bias_spline /
# fit_monotone_convex_spread_spline used a fixed lr=1e-4. opnorm(H)
# scales with the row count of the design matrix, so once n_pts pushes
# L = opnorm(H) past 1/lr the optimizer diverges, IRLS reweights to
# non-finite values, and the next H \ Xty throws
# `ArgumentError: matrix contains Infs or NaNs` from LAPACK.
#
# These tests build large-N synthetic inputs (300k points) that would
# trip the original divergence and assert the fit returns a finite
# UniformSpline.

using Random
using Statistics: median, quantile
using Distributions: Normal
import Pioneer
using Pioneer: fit_convex_bias_spline, fit_monotone_bias_spline,
               fit_binned_regularized_bias_spline,
               fit_intensity_mass_error_model,
               fit_scout_calibrated_model,
               fit_best_convex_bias,
               fit_monotone_convex_spread_spline,
               _select_regularized_mz_spread_correction,
               _mz_spread_correction,
               binned_laplace_scale,
               make_edge_point_spline_extrap,
               make_rt_range_edge_spline_extrap,
               spline_derivative,
               IntensityMassErrorModel,
               MassErrSample,
               ScoutCalibratedMassErrorModel,
               SimpleMassErrorModel,
               EMPK_CLAMP_HI,
               EMPK_MIN_K,
               EMPK_QUANTILE,
               TUNING_BIAS_EXTRAP_EDGE_POINTS,
               TUNING_CALIBRATION_BIN_SIZE,
               TUNING_MIN_SAMPLES,
               TUNING_MZ_BIAS_BIN_SIZE,
               TUNING_MZ_BIAS_KNOTS,
               TUNING_MZ_BIAS_LAMBDA,
               TUNING_RT_BIAS_BIN_SIZE,
               TUNING_RT_BIAS_KNOTS,
               TUNING_RT_BIAS_LAMBDA

@testset "Empirical k constants use matching Gaussian lower bound" begin
    @test 0.95 <= EMPK_QUANTILE < 1.0
    @test EMPK_MIN_K ≈ Float32(quantile(Normal(), (1.0 + EMPK_QUANTILE) / 2.0))
    @test EMPK_MIN_K < EMPK_CLAMP_HI
end

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

    @testset "fit_monotone_convex_spread_spline" begin
        # The spread fit consumes pre-binned (centers, scales) so its
        # design matrix is small; the Lipschitz cap still has to be
        # finite. This guards against regression of the same code path.
        centers = collect(range(10.0, 30.0; length=200))
        scales = max.(0.5 .+ 5.0 .* exp.(-(centers .- 10.0) ./ 5.0)
                       .+ 0.05 .* randn(rng, length(centers)), 0.05)
        spline = fit_monotone_convex_spread_spline(centers, scales;
                                                    n_knots=10, λ=0.01,
                                                    max_iter=1000, lr=1e-3)
        @test spline !== nothing
        @test all(isfinite, spline.coeffs)
    end
end

@testset "Parameter tuning calibration bin size" begin
    @test TUNING_CALIBRATION_BIN_SIZE == 200

    order = collect(1:400)
    vals = Float64.(1:400)
    errors = vcat(fill(1.0, 200), fill(2.0, 200))

    centers, b_vals = binned_laplace_scale(
        order, vals, errors, TUNING_CALIBRATION_BIN_SIZE)

    @test centers == [100.5, 300.5]
    @test length(b_vals) == 2
end

@testset "Parameter tuning m/z bias smoothing" begin
    @test TUNING_MIN_SAMPLES == 1200
    @test TUNING_MZ_BIAS_LAMBDA == 0.1
    @test TUNING_MZ_BIAS_BIN_SIZE == 100
    @test TUNING_MZ_BIAS_KNOTS == 24
    @test TUNING_BIAS_EXTRAP_EDGE_POINTS == 5
end

@testset "Parameter tuning RT bias tracks local transitions" begin
    @test TUNING_RT_BIAS_BIN_SIZE == 100
    @test TUNING_RT_BIAS_KNOTS == 48
    @test TUNING_RT_BIAS_LAMBDA == 0.03

    rng = MersenneTwister(0x7254)
    rt = repeat(collect(range(10.0, 80.0; length=71)), inner=100)
    truth = @. 0.0006 -
               0.0022 / (1 + exp(-(rt - 61.0) / 1.1)) +
               0.00035 * exp(-((rt - 53.0) / 2.2)^2)
    residuals = truth .+ 0.00012 .* randn(rng, length(rt))

    spline, label = fit_binned_regularized_bias_spline(
        rt, residuals;
        n_knots=TUNING_RT_BIAS_KNOTS,
        λ=TUNING_RT_BIAS_LAMBDA,
        bin_size=TUNING_RT_BIAS_BIN_SIZE)

    probes = [53.0, 58.0, 61.0, 64.0, 70.0]
    expected = @. 0.0006 -
                  0.0022 / (1 + exp(-(probes - 61.0) / 1.1)) +
                  0.00035 * exp(-((probes - 53.0) / 2.2)^2)
    fitted = [spline(p) for p in probes]

    @test occursin("binned regularized", label)
    @test median(abs.(fitted .- expected)) < 3e-5
    @test abs(spline(64.0) - expected[4]) < 3e-5
end

@testset "Spread spline tracks high-intensity tail structure" begin
    centers = collect(range(8.0, 24.0; length=81))
    truth = @. 3.8 - 2.2 / (1 + exp(-(centers - 13.0) / 0.9)) +
               0.16 * exp(-((centers - 16.6) / 0.65)^2) -
               0.20 * exp(-((centers - 19.0) / 0.8)^2) +
               0.58 / (1 + exp(-(centers - 21.4) / 0.55))
    spline = fit_monotone_convex_spread_spline(centers, truth;
        n_knots=20, λ=0.01, max_iter=1000, lr=1e-3)

    probes = [16.6, 19.0, 21.4, 23.0]
    expected = @. 3.8 - 2.2 / (1 + exp(-(probes - 13.0) / 0.9)) +
                  0.16 * exp(-((probes - 16.6) / 0.65)^2) -
                  0.20 * exp(-((probes - 19.0) / 0.8)^2) +
                  0.58 / (1 + exp(-(probes - 21.4) / 0.55))
    fitted = [spline(p) for p in probes]

    @test spline !== nothing
    @test all(isfinite, spline.coeffs)
    @test median(abs.(fitted .- expected)) < 0.035
    @test spline(23.0) - spline(19.0) > 0.35
end

function _expected_binned_bounds(x, y; bin_size)
    valid = [
        i for i in eachindex(x)
        if isfinite(Float64(x[i])) && isfinite(Float64(y[i]))
    ]
    order = valid[sortperm(Float64.(x[valid]))]
    n = length(order)
    n_bins = max(n ÷ bin_size, 1)
    centers = Float64[]
    for b in 1:n_bins
        s = (b - 1) * bin_size + 1
        e = b == n_bins ? n : b * bin_size
        push!(centers, median(Float64.(x[order[s:e]])))
    end
    return first(centers), last(centers)
end

@testset "Bias extrapolation slopes use first and last five binned medians" begin
    n = 2_000
    x = collect(range(0.0, 100.0; length=n))
    y = zeros(Float64, n)
    y[x .<= 25.0] .= 0.003 .* x[x .<= 25.0] .+ 0.01
    y[x .>= 75.0] .= -0.004 .* x[x .>= 75.0] .+ 0.5
    y[(x .> 25.0) .& (x .< 75.0)] .= @. 0.02 * sin(x[(x .> 25.0) & (x .< 75.0)])

    spline = fit_convex_bias_spline(x, @. 0.01 * (x - 50.0)^2;
        convex_up=true, n_knots=10, λ=1.0, max_iter=1000)
    lo, hi = quantile(x, [0.02, 0.98])
    extrap = make_edge_point_spline_extrap(
        spline, x, y, lo, hi; n_edge_points=5, bin_size=100)
    expected_lo, expected_hi = _expected_binned_bounds(
        x[(lo .<= x) .& (x .<= hi)], y[(lo .<= x) .& (x .<= hi)];
        bin_size=100)

    @test extrap.lo ≈ expected_lo
    @test extrap.hi ≈ expected_hi
    @test extrap.lo_slope ≈ 0.003 atol=1e-5
    @test extrap.hi_slope ≈ -0.004 atol=1e-5
    @test !(extrap.lo_slope ≈ spline_derivative(spline, extrap.lo))
    @test !(extrap.hi_slope ≈ spline_derivative(spline, extrap.hi))
end

@testset "RT bias extrapolation is constant outside binned RT range" begin
    x = vcat(
        collect(range(0.0, 99.0; length=2_000)),
        collect(range(99.80, 100.00; length=500)),
    )
    y = @. 0.001 * x
    dense = x .>= 99.80
    y[dense] .= @. 0.001 * x[dense] + 0.040 * (x[dense] - 99.80)

    spline = fit_convex_bias_spline(x, @. 0.0 * x;
        convex_up=true, n_knots=10, λ=1.0, max_iter=1000)
    lo, hi = quantile(x, [0.02, 0.98])

    clustered_extrap = make_edge_point_spline_extrap(
        spline, x, y, lo, hi; n_edge_points=5, bin_size=100)
    rt_extrap = make_rt_range_edge_spline_extrap(
        spline, x, y, lo, hi; bound_bin_size=100)
    expected_lo, expected_hi = _expected_binned_bounds(
        x[(lo .<= x) .& (x .<= hi)], y[(lo .<= x) .& (x .<= hi)];
        bin_size=100)

    @test rt_extrap.lo ≈ expected_lo
    @test rt_extrap.hi ≈ expected_hi
    @test clustered_extrap.hi_slope > 0.02
    @test rt_extrap.lo_slope == 0.0
    @test rt_extrap.hi_slope == 0.0
end

@testset "Scout calibrated bias extrapolation keeps spline endpoint slopes" begin
    n = 2_000
    mz = collect(range(400.0, 1200.0; length=n))
    log2I = collect(range(10.0, 20.0; length=n))
    da_err = @. 2e-6 * (mz - 800.0)

    left = (mz .>= quantile(mz, 0.02)) .& (mz .<= quantile(mz, 0.05))
    right = (mz .>= quantile(mz, 0.95)) .& (mz .<= quantile(mz, 0.98))
    left_idx = findall(left)
    right_idx = findall(right)
    da_err[left_idx] .= 1e-3 .+ 2e-4 .* (mz[left_idx] .- mz[first(left_idx)])
    da_err[right_idx] .= -1e-3 .- 2e-4 .* (mz[right_idx] .- mz[first(right_idx)])

    samples = [
        MassErrSample(Float32(mz[i]), Float32(mz[i] + da_err[i]), Float32(2.0 ^ log2I[i]), Float32(i))
        for i in eachindex(mz)
    ]
    model = fit_scout_calibrated_model(samples; with_intensity=false)
    expected_lo, expected_hi = _expected_binned_bounds(
        Float64.(Float32.(mz)), Float64.(Float32.(mz .+ da_err)) .- Float64.(Float32.(mz));
        bin_size=TUNING_CALIBRATION_BIN_SIZE)

    @test model isa ScoutCalibratedMassErrorModel
    @test model.mz_bias_extrap.lo ≈ expected_lo
    @test model.mz_bias_extrap.hi ≈ expected_hi
    @test !(model.mz_bias_extrap.lo ≈ Float32(quantile(mz, 0.02)))
    @test !(model.mz_bias_extrap.hi ≈ Float32(quantile(mz, 0.98)))
    @test model.mz_bias_extrap.lo_slope ≈
          spline_derivative(model.mz_bias_spline, model.mz_bias_extrap.lo)
    @test model.mz_bias_extrap.hi_slope ≈
          spline_derivative(model.mz_bias_spline, model.mz_bias_extrap.hi)
end

@testset "Intensity mass error bias extrapolation uses boundary binned-median slopes" begin
    rng = MersenneTwister(0x199)
    n = 1_200
    mz = collect(range(400.0, 1600.0; length=n))
    log2I = collect(range(8.0, 22.0; length=n))
    intensities = Float32.(2.0 .^ log2I)
    da_bias = @. 1.5e-6 * (mz - 900.0) + 2.0e-4 * tanh((log2I - 16.0) / 2.0)
    da_err = da_bias .+ 4.0e-5 .* randn(rng, n)

    samples = [
        MassErrSample(Float32(mz[i]), Float32(mz[i] + da_err[i]), intensities[i], Float32(i))
        for i in eachindex(mz)
    ]
    old_model = SimpleMassErrorModel(0.0f0, (20.0f0, 20.0f0))
    model = fit_intensity_mass_error_model(samples, old_model)

    @test model isa IntensityMassErrorModel

    fit_mz = Float64.(Float32.(mz))
    fit_obs_mz = Float64.(Float32.(mz .+ da_err))
    fit_da_err = fit_obs_mz .- fit_mz
    fit_log2I = log2.(Float64.(intensities))

    expected_mz_extrap = make_edge_point_spline_extrap(
        model.mz_bias_spline, fit_mz, fit_da_err,
        Float32(minimum(fit_mz)), Float32(maximum(fit_mz)))
    mz_residuals = fit_da_err .- [Float64(model.mz_bias_spline(Float32(m))) for m in fit_mz]
    expected_int_extrap = make_edge_point_spline_extrap(
        model.intensity_bias_spline, fit_log2I, mz_residuals,
        Float32(minimum(fit_log2I)), Float32(maximum(fit_log2I)))

    @test model.mz_bias_extrap.lo ≈ expected_mz_extrap.lo
    @test model.mz_bias_extrap.hi ≈ expected_mz_extrap.hi
    @test model.int_bias_extrap.lo ≈ expected_int_extrap.lo
    @test model.int_bias_extrap.hi ≈ expected_int_extrap.hi
    @test !(model.mz_bias_extrap.lo ≈ Float32(quantile(mz, 0.02)))
    @test !(model.mz_bias_extrap.hi ≈ Float32(quantile(mz, 0.98)))
    @test model.mz_bias_extrap.lo_slope ≈ expected_mz_extrap.lo_slope atol=1e-8
    @test model.mz_bias_extrap.hi_slope ≈ expected_mz_extrap.hi_slope atol=1e-8
    @test model.int_bias_extrap.lo_slope ≈ expected_int_extrap.lo_slope atol=1e-7
    @test model.int_bias_extrap.hi_slope ≈ expected_int_extrap.hi_slope atol=1e-7
    @test model.spread_extrap.lo ≈ Float32(quantile(log2I, 0.05))
    @test model.spread_extrap.hi ≈ model.spread_spline.last
    @test model.spread_extrap.lo < model.spread_extrap.hi
    @test model.spread_extrap.lo_slope == 0.0f0
    @test model.spread_extrap.hi_slope == 0.0f0
end

@testset "Binned regularized m/z bias spline captures local edge structure" begin
    rng = MersenneTwister(0x51506c)
    x = repeat(collect(range(250.0, 1450.0; length=37)), inner=120)
    truth = @. 0.0025 * exp(-((x - 325.0) / 95.0)^2) -
               0.0018 * exp(-((x - 610.0) / 130.0)^2) +
               0.0008 * exp(-((x - 1180.0) / 180.0)^2)
    y = truth .+ 0.00018 .* randn(rng, length(x))

    binned_spline, binned_label = fit_binned_regularized_bias_spline(
        x, y; n_knots=12, λ=0.1, bin_size=120)
    convex_spline, _ = fit_best_convex_bias(x, y; n_knots=12, λ=0.1)

    probe = collect(260.0:40.0:1420.0)
    expected = @. 0.0025 * exp(-((probe - 325.0) / 95.0)^2) -
                  0.0018 * exp(-((probe - 610.0) / 130.0)^2) +
                  0.0008 * exp(-((probe - 1180.0) / 180.0)^2)

    binned_err = median(abs.(expected .- [binned_spline(p) for p in probe]))
    convex_err = median(abs.(expected .- [convex_spline(p) for p in probe]))

    @test occursin("binned regularized", binned_label)
    @test all(isfinite, binned_spline.coeffs)
    @test binned_err < 0.6 * convex_err
    @test binned_spline(325.0) > binned_spline(610.0)
end

@testset "Robust linear m/z spread correction ignores outlier bins" begin
    corr_centers = collect(300.0:100.0:1300.0)
    corr_vals = @. 0.45 + 0.00075 * corr_centers
    corr_vals[3] = 4.8
    corr_vals[9] = 0.18

    spline, extrap, α, β, γ, label =
        _select_regularized_mz_spread_correction(corr_centers, corr_vals)

    @test occursin("robust linear", label)
    @test all(isfinite, spline.coeffs)
    @test -0.1f0 < α < 0.9f0
    @test 0.0004f0 < β < 0.0012f0
    @test γ == 0.0f0

    low = _mz_spread_correction(spline, extrap, 300.0f0)
    mid = _mz_spread_correction(spline, extrap, 800.0f0)
    high = _mz_spread_correction(spline, extrap, 1300.0f0)
    outlier_hi = _mz_spread_correction(spline, extrap, 500.0f0)
    outlier_lo = _mz_spread_correction(spline, extrap, 1100.0f0)

    @test high > low + 0.55f0
    @test abs(mid - Float32(0.45 + 0.00075 * 800.0)) < 0.08f0
    @test outlier_hi < 2.0f0
    @test outlier_lo > 1.0f0
end
