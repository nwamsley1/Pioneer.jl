# Tests for src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl.
#
# Targets the three exported functions:
# - fit_linear_irt_model: robust (Tukey-bisquare) linear regression for RT→iRT
# - make_spline_monotonic: enforces monotonicity on a fitted spline via
#   bidirectional cumulative-max from the median
# - fit_irt_model: dispatches identity / linear / spline based on PSM count

using DataFrames
using Random

using Pioneer: fit_linear_irt_model, make_spline_monotonic, fit_irt_model
using Pioneer: LinearRtConversionModel, IdentityModel,
               InterpolationRtConversionModel, RtConversionModel
using Pioneer: UniformSpline

@testset "rt_alignment_utils" begin

    @testset "fit_linear_irt_model: clean linear recovers slope/intercept" begin
        rt = Float32.(collect(1.0:0.5:50.0))
        # iRT = 1.5 * RT + 10 (clean, no noise)
        irt = Float32.(1.5f0 .* rt .+ 10.0f0)
        psms = DataFrame(rt = rt, irt_predicted = irt)
        model, residual_std, n_params = fit_linear_irt_model(psms)
        @test model isa LinearRtConversionModel
        @test n_params == 2
        # Apply model — should match irt at sample points within tight tolerance
        for r in rt[1:5:end]
            @test model(r) ≈ 1.5f0 * r + 10.0f0 atol=1f-2
        end
        # Clean data → near-zero residual std
        @test residual_std < 1f-2
    end

    @testset "fit_linear_irt_model: degenerate (all RTs identical)" begin
        rt = fill(20.0f0, 10)
        irt = Float32.(rand(MersenneTwister(7), 10) .* 100)
        psms = DataFrame(rt = rt, irt_predicted = irt)
        model, residual_std, n_params = fit_linear_irt_model(psms)
        @test model isa LinearRtConversionModel
        @test n_params == 2
        # Slope must be 0 in the degenerate path; intercept ≈ mean(irt)
        @test model(0.0f0) ≈ Float32(mean(irt)) atol=1f-3
        # Residual std equals std(irt) per the degenerate-branch contract
        @test residual_std ≈ Float32(std(irt)) atol=1f-3
    end

    @testset "fit_linear_irt_model: robust to outliers (5%)" begin
        rng = MersenneTwister(2026)
        rt = Float32.(collect(LinRange(1.0, 50.0, 200)))
        irt = Float32.(2.0f0 .* rt .+ 5.0f0 .+ 0.5f0 .* randn(rng, length(rt)))
        # Inject ~5% wild outliers
        outlier_idx = sort!(randperm(rng, length(rt))[1:10])
        irt[outlier_idx] .+= Float32.(rand(rng, [-200.0, 200.0], length(outlier_idx)))
        psms = DataFrame(rt = rt, irt_predicted = irt)

        model, residual_std, _ = fit_linear_irt_model(psms)
        @test model isa LinearRtConversionModel
        # Robust slope should still be close to 2.0 despite outliers; OLS would
        # be pulled noticeably off.
        # Recover slope from two well-spaced predictions.
        slope_est = (model(40.0f0) - model(10.0f0)) / 30.0f0
        @test slope_est ≈ 2.0f0 atol=0.1f0
    end

    @testset "make_spline_monotonic: enforces monotonicity on non-monotonic input" begin
        # Build a spline whose raw output is non-monotonic by design (sinusoid
        # over [0, 2π]).
        rt = Float32.(collect(LinRange(0.0, 2π, 60)))
        irt = Float32.(sin.(rt))
        spline = UniformSpline(irt, rt, 3, 10)
        @test minimum(diff(Float32[spline(r) for r in rt])) < 0   # confirm non-monotonic

        mono = make_spline_monotonic(spline, rt, irt)
        # Sample on the same grid; monotonic non-decreasing
        sampled = Float32[mono(r) for r in rt]
        @test all(diff(sampled) .≥ -1f-6)
    end

    @testset "make_spline_monotonic: already-monotonic spline preserved" begin
        rt = Float32.(collect(LinRange(1.0, 50.0, 80)))
        irt = Float32.(0.5f0 .* rt .+ 2.0f0)   # linear, already monotonic
        spline = UniformSpline(irt, rt, 3, 8)

        mono = make_spline_monotonic(spline, rt, irt)
        sampled = Float32[mono(r) for r in rt]
        @test all(diff(sampled) .≥ -1f-6)
        # Output close to the original linear function on the interior
        for r in rt[10:10:end-10]
            @test mono(r) ≈ 0.5f0 * r + 2.0f0 atol=0.5f0
        end
    end

    @testset "fit_irt_model: insufficient data → IdentityModel" begin
        # min_psms default is 30; supply 10
        rt = Float32.(collect(1.0:1.0:10.0))
        irt = Float32.(rt .+ 1.0)
        psms = DataFrame(rt = rt, irt_predicted = irt)
        model, rt_out, irt_out, mad_out = fit_irt_model(psms; min_psms = 30)
        @test model isa IdentityModel
        @test isempty(rt_out)
        @test isempty(irt_out)
        @test mad_out == 0.0f0
    end

    @testset "fit_irt_model: small-dataset linear branch" begin
        # 20 PSMs (≥min_psms=10 but <30) → linear branch
        rt = Float32.(collect(1.0:1.0:20.0))
        irt = Float32.(2.0f0 .* rt .+ 3.0f0)
        psms = DataFrame(rt = rt, irt_predicted = irt)
        model, rt_out, irt_out, mad_out = fit_irt_model(psms; min_psms = 10)
        @test model isa LinearRtConversionModel
        @test length(rt_out)  == 20
        @test length(irt_out) == 20
        @test mad_out ≥ 0.0f0
        # Linear model should recover the slope on clean data
        @test (model(15.0f0) - model(5.0f0)) / 10.0f0 ≈ 2.0f0 atol=0.1f0
    end

    @testset "fit_irt_model: many-PSM branch returns a callable RtConversionModel" begin
        rng = MersenneTwister(99)
        rt = Float32.(collect(LinRange(1.0, 60.0, 200)))
        irt = Float32.(0.8f0 .* rt .+ 4.0f0 .+ 0.2f0 .* randn(rng, length(rt)))
        psms = DataFrame(rt = rt, irt_predicted = irt)
        model, rt_out, irt_out, mad_out = fit_irt_model(psms;
            min_psms = 30, ransac_threshold = 1000)

        @test model isa RtConversionModel
        # rt_out and irt_out should match the input length (after outlier removal,
        # most points retained for clean data); test for non-empty.
        @test length(rt_out) > 50
        @test length(irt_out) == length(rt_out)
        @test mad_out ≥ 0.0f0
        # Sanity: applying the model recovers the linear trend within ~5 iRT units
        for r in rt[20:30:180]
            @test abs(model(r) - (0.8f0 * r + 4.0f0)) ≤ 5.0f0
        end
    end

end
