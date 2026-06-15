# Tests for SimpleMassErrorModel, LinearDaMassErrorModel, and
# LinearBiasPpmTolMassErrorModel in src/structs/MassErrorModel.jl.
#
# Targets the getter / corrected-mz / bounds API including:
# - asymmetric Simple model (left ≠ right) and the reverse-swap semantics
#   of getMzBoundsReverse
# - degenerate zero-bias / zero-tolerance models
# - 3-arg forwarding (intensity-aware overload that ignores intensity)
# - getCorrectedMzAndBounds internal consistency
# - the Linear* models' Laplace log-density fallbacks

using Pioneer: SimpleMassErrorModel, LinearDaMassErrorModel,
               LinearBiasPpmTolMassErrorModel, IntensityMassErrorModel,
               UniformSpline, SplineExtrap,
               getRightTol, getLeftTol, getMassOffset, getMassCorrection,
               getCorrectedMz, getMzBounds, getMzBoundsReverse,
               getCorrectedMzAndBounds, laplace_log_density,
               get_default_top3_ll,
               _model_tolerance_heatmap_grid
import Pioneer
using StaticArrays: SVector

_test_const_spline_f32(value::Float32; first::Float32 = 0f0, last::Float32 = 1f0) =
    UniformSpline{12, Float32}(
        SVector{12, Float32}(
            value, 0f0, 0f0, 0f0,
            value, 0f0, 0f0, 0f0,
            value, 0f0, 0f0, 0f0,
        ),
        3,
        first,
        last,
        (last - first) / 2f0,
    )

_test_linear_spline_f32(intercept::Float32, slope::Float32; first::Float32 = 0f0, last::Float32 = 1f0) =
    let bin_width = (last - first) / 2f0
        UniformSpline{12, Float32}(
            SVector{12, Float32}(
                intercept, slope * bin_width, 0f0, 0f0,
                intercept + slope * bin_width, slope * bin_width, 0f0, 0f0,
                intercept + slope * (last - first), slope * bin_width, 0f0, 0f0,
            ),
            3,
            first,
            last,
            bin_width,
        )
    end

_test_extrap_f32(spline::UniformSpline{N, Float32}) where N =
    SplineExtrap{Float32}(
        spline.first,
        spline.last,
        spline(spline.first),
        0f0,
        spline(spline.last),
        0f0,
    )

@testset "MassErrorModel" begin

    @testset "SimpleMassErrorModel" begin
        @testset "symmetric ±10 ppm, zero bias" begin
            mem = SimpleMassErrorModel(0.0f0, (10.0f0, 10.0f0))
            @test getLeftTol(mem) == 10.0f0
            @test getRightTol(mem) == 10.0f0
            @test getMassOffset(mem) == 0.0f0
            @test getMassCorrection(mem) == 0.0f0
            @test getCorrectedMz(mem, 1000.0f0) == 1000.0f0
            # Zero bias → 3-arg matches 2-arg
            @test getCorrectedMz(mem, 1000.0f0, 0.5f0) == getCorrectedMz(mem, 1000.0f0)

            low, high = getMzBounds(mem, 1000.0f0)
            @test low ≈ 1000.0f0 - 0.01f0   # 10 ppm at 1000 = 0.01
            @test high ≈ 1000.0f0 + 0.01f0
            # Symmetric → forward and reverse identical
            @test getMzBoundsReverse(mem, 1000.0f0) == getMzBounds(mem, 1000.0f0)
        end

        @testset "asymmetric tolerance: reverse swaps left/right" begin
            # left = 20 ppm (empirical can be 20 ppm BELOW theoretical)
            # right = 5 ppm (empirical can be 5 ppm ABOVE theoretical)
            mem = SimpleMassErrorModel(0.0f0, (20.0f0, 5.0f0))
            @test getLeftTol(mem) == 20.0f0
            @test getRightTol(mem) == 5.0f0

            # Forward: empirical bounds for theoretical = 1000.
            # Empirical ∈ (theoretical - left*ppm, theoretical + right*ppm)
            low, high = getMzBounds(mem, 1000.0f0)
            @test low  ≈ 1000.0f0 - 0.020f0   # 20 ppm
            @test high ≈ 1000.0f0 + 0.005f0   # 5 ppm

            # Reverse: theoretical bounds for empirical = 1000.
            # Theoretical ∈ (empirical - right*ppm, empirical + left*ppm)
            rlow, rhigh = getMzBoundsReverse(mem, 1000.0f0)
            @test rlow  ≈ 1000.0f0 - 0.005f0   # right swap
            @test rhigh ≈ 1000.0f0 + 0.020f0   # left swap

            # If we feed the high end of the forward bound back through reverse,
            # we should bracket the original theoretical (≤ on the lower edge —
            # the round-trip can land exactly at 1000.0 in Float32).
            rlow2, rhigh2 = getMzBoundsReverse(mem, high)
            @test rlow2 ≤ 1000.0f0 ≤ rhigh2
        end

        @testset "non-zero bias correction" begin
            # +5 ppm bias: corrected = mz - 5 * mz/1e6
            mem = SimpleMassErrorModel(5.0f0, (10.0f0, 10.0f0))
            @test getMassOffset(mem) == 5.0f0
            corrected = getCorrectedMz(mem, 1000.0f0)
            @test corrected ≈ 1000.0f0 - 0.005f0
            @test corrected < 1000.0f0  # positive bias drags corrected down

            # Negative bias: corrected drifts up
            mem_neg = SimpleMassErrorModel(-7.5f0, (10.0f0, 10.0f0))
            corrected_neg = getCorrectedMz(mem_neg, 1000.0f0)
            @test corrected_neg ≈ 1000.0f0 + 0.0075f0
            @test corrected_neg > 1000.0f0
        end

        @testset "zero-tolerance: bounds collapse to point" begin
            mem = SimpleMassErrorModel(0.0f0, (0.0f0, 0.0f0))
            low, high = getMzBounds(mem, 500.0f0)
            @test low == high == 500.0f0
            rlow, rhigh = getMzBoundsReverse(mem, 500.0f0)
            @test rlow == rhigh == 500.0f0
        end

        @testset "very small and very large mz" begin
            mem = SimpleMassErrorModel(0.0f0, (10.0f0, 10.0f0))
            for mz in (Float32(50), Float32(100), Float32(2000))
                low, high = getMzBounds(mem, mz)
                @test low < mz < high
                # Compare against the function's exact computation (avoid
                # subtraction rounding near magnitudes ≫ tolerance).
                ppm = 10.0f0 * (mz / 1f6)
                @test high == Float32(mz + ppm)
                @test low  == Float32(mz - ppm)
            end
        end

        @testset "getCorrectedMzAndBounds: internal consistency" begin
            mem = SimpleMassErrorModel(3.0f0, (15.0f0, 8.0f0))
            corrected, low, high = getCorrectedMzAndBounds(mem, 800.0f0, 0.0f0)
            @test corrected == getCorrectedMz(mem, 800.0f0)
            expected_low, expected_high = getMzBoundsReverse(mem, corrected)
            @test low == expected_low
            @test high == expected_high
        end

        @testset "3-arg forwarding ignores intensity argument" begin
            mem = SimpleMassErrorModel(2.5f0, (12.0f0, 9.0f0))
            mz = 750.0f0
            @test getCorrectedMz(mem, mz, 0.1f0) == getCorrectedMz(mem, mz, 1.0f6)
            @test getCorrectedMz(mem, mz, 0.0f0) == getCorrectedMz(mem, mz)
            @test getMzBoundsReverse(mem, mz, 5.0f0) == getMzBoundsReverse(mem, mz)
        end
    end

    @testset "LinearDaMassErrorModel" begin
        @testset "zero coefficients = identity correction" begin
            m = LinearDaMassErrorModel(0.0f0, 0.0f0, 0.005f0)
            @test getCorrectedMz(m, 1000.0f0) == 1000.0f0
            @test getMassOffset(m) == 0.0f0
            @test getMassCorrection(m) == 0.0f0
            @test getRightTol(m) == 0.005f0
            @test getLeftTol(m) == 0.005f0
            low, high = getMzBounds(m, 1000.0f0)
            @test low  == 1000.0f0 - 0.005f0
            @test high == 1000.0f0 + 0.005f0
            # LinearDa is symmetric → forward == reverse
            @test getMzBoundsReverse(m, 1000.0f0) == getMzBounds(m, 1000.0f0)
        end

        @testset "intercept-only bias" begin
            m = LinearDaMassErrorModel(0.01f0, 0.0f0, 0.001f0)
            # corrected = mz - (0.01 + 0)
            @test getCorrectedMz(m, 1000.0f0) ≈ 999.99f0
        end

        @testset "intercept + slope: full linear correction" begin
            # bias_da(mz) = 0.001 + 0.000002 * mz
            m = LinearDaMassErrorModel(0.001f0, 2.0f-6, 0.002f0)
            mz = 500.0f0
            expected_bias = 0.001f0 + 2.0f-6 * mz
            @test getCorrectedMz(m, mz) ≈ mz - expected_bias

            # Verify slope contribution scales with mz
            high_mz_bias = 0.001f0 + 2.0f-6 * 2000.0f0
            @test getCorrectedMz(m, 2000.0f0) ≈ 2000.0f0 - high_mz_bias
            # Difference between mz=500 and mz=2000 corrected should reflect slope
            @test (2000.0f0 - getCorrectedMz(m, 2000.0f0)) -
                  (500.0f0 - getCorrectedMz(m, 500.0f0)) ≈ 2.0f-6 * (2000.0f0 - 500.0f0) atol=1f-5
        end

        @testset "tolerance is constant Da regardless of mz" begin
            m = LinearDaMassErrorModel(0.0f0, 0.0f0, 0.003f0)
            for mz in Float32[150, 500, 1000, 2000]
                low, high = getMzBounds(m, mz)
                # Bounds match the function's Float32 computation exactly;
                # don't try to reconstruct the tolerance via subtraction (rounding).
                @test high == Float32(mz + 0.003f0)
                @test low  == Float32(mz - 0.003f0)
            end
        end

        @testset "getCorrectedMzAndBounds: corrected propagates to bounds" begin
            m = LinearDaMassErrorModel(0.005f0, 1.0f-6, 0.002f0)
            corrected, low, high = getCorrectedMzAndBounds(m, 1500.0f0, 0.0f0)
            @test corrected == getCorrectedMz(m, 1500.0f0)
            expected_low, expected_high = getMzBoundsReverse(m, corrected)
            @test low == expected_low
            @test high == expected_high
            # (Width = 2 × tolerance_da is implied by the equality checks above;
            #  computing high - low directly loses precision in Float32 at the
            #  ~1500 mass scale.)
        end

        @testset "Laplace log-density and top3 fallbacks return zero" begin
            m = LinearDaMassErrorModel(0.0f0, 0.0f0, 0.001f0)
            @test laplace_log_density(m, 1000.0f0, 0.5f0, 0.001f0) == 0.0f0
            @test get_default_top3_ll(m) == 0.0f0
        end

        @testset "3-arg forwarding ignores intensity" begin
            m = LinearDaMassErrorModel(0.001f0, 0.0f0, 0.002f0)
            @test getCorrectedMz(m, 1000.0f0, 0.5f0) == getCorrectedMz(m, 1000.0f0)
            @test getMzBoundsReverse(m, 1000.0f0, 0.5f0) == getMzBoundsReverse(m, 1000.0f0)
        end
    end

    @testset "LinearBiasPpmTolMassErrorModel" begin
        @testset "zero coefficients, ±10 ppm tolerance" begin
            m = LinearBiasPpmTolMassErrorModel(0.0f0, 0.0f0, 10.0f0)
            @test getCorrectedMz(m, 1000.0f0) == 1000.0f0
            @test getMassOffset(m) == 0.0f0
            @test getMassCorrection(m) == 0.0f0
            @test getRightTol(m) == 10.0f0
            @test getLeftTol(m) == 10.0f0

            low, high = getMzBounds(m, 1000.0f0)
            @test low  ≈ 1000.0f0 - 0.010f0
            @test high ≈ 1000.0f0 + 0.010f0
            # Symmetric → forward == reverse
            @test getMzBoundsReverse(m, 1000.0f0) == getMzBounds(m, 1000.0f0)
        end

        @testset "linear bias correction with ppm tolerance" begin
            m = LinearBiasPpmTolMassErrorModel(0.002f0, 5.0f-7, 8.0f0)
            # Correction matches LinearDa formula
            @test getCorrectedMz(m, 1000.0f0) ≈ 1000.0f0 - (0.002f0 + 5.0f-7 * 1000.0f0)
            # Tolerance scales with mz (ppm semantics). Compare against the
            # function's exact Float32 expression to avoid subtraction rounding.
            for mz in Float32[200, 500, 1500]
                low, high = getMzBounds(m, mz)
                hw = 8.0f0 * (mz / 1f6)
                @test high == Float32(mz + hw)
                @test low  == Float32(mz - hw)
            end
        end

        @testset "getCorrectedMzAndBounds: bounds based on corrected mass" begin
            m = LinearBiasPpmTolMassErrorModel(0.001f0, 1.0f-6, 12.0f0)
            corrected, low, high = getCorrectedMzAndBounds(m, 1000.0f0, 0.0f0)
            @test corrected == getCorrectedMz(m, 1000.0f0)
            expected_low, expected_high = getMzBoundsReverse(m, corrected)
            @test low == expected_low
            @test high == expected_high
            # The bounds are ±12 ppm of CORRECTED mass, not of input.
            hw = 12.0f0 * (corrected / 1f6)
            @test high == Float32(corrected + hw)
            @test low  == Float32(corrected - hw)
        end

        @testset "Laplace log-density and top3 fallbacks return zero" begin
            m = LinearBiasPpmTolMassErrorModel(0.0f0, 0.0f0, 10.0f0)
            @test laplace_log_density(m, 1000.0f0, 0.5f0, 0.001f0) == 0.0f0
            @test get_default_top3_ll(m) == 0.0f0
        end

        @testset "3-arg forwarding ignores intensity" begin
            m = LinearBiasPpmTolMassErrorModel(0.0005f0, 0.0f0, 6.0f0)
            @test getCorrectedMz(m, 1000.0f0, 0.5f0) == getCorrectedMz(m, 1000.0f0)
            @test getMzBoundsReverse(m, 1000.0f0, 0.5f0) == getMzBoundsReverse(m, 1000.0f0)
        end
    end

    @testset "IntensityMassErrorModel RT bias" begin
        mz_bias = _test_const_spline_f32(0f0; first=100f0, last=1000f0)
        intensity_bias = _test_const_spline_f32(0f0; first=0f0, last=20f0)
        spread = _test_const_spline_f32(1f0; first=0f0, last=20f0)
        mz_spread = _test_const_spline_f32(1f0; first=100f0, last=1000f0)
        rt_bias = _test_linear_spline_f32(0f0, 0.0004f0; first=10f0, last=20f0)

        mem = IntensityMassErrorModel(
            mz_bias,
            intensity_bias,
            spread,
            rt_bias,
            _test_extrap_f32(mz_bias),
            _test_extrap_f32(intensity_bias),
            _test_extrap_f32(spread),
            _test_extrap_f32(rt_bias),
            1.96f0,
            0.02f0,
            mz_spread,
            _test_extrap_f32(mz_spread),
            1f0,
            0f0,
            0f0,
            0f0,
            0.02f0,
            10f0,
            20f0,
        )

        mz = 500f0
        intensity = 1024f0
        corrected_no_rt, _, _ = getCorrectedMzAndBounds(mem, mz, intensity)
        corrected_mid_rt, low_mid_rt, high_mid_rt = getCorrectedMzAndBounds(mem, mz, intensity, 15f0)
        corrected_hi_rt, _, _ = getCorrectedMzAndBounds(mem, mz, intensity, 20f0)

        @test corrected_no_rt == mz
        @test corrected_mid_rt ≈ mz - 0.002f0 atol=2f-5
        @test corrected_hi_rt ≈ mz - 0.004f0 atol=2f-5
        @test low_mid_rt < corrected_mid_rt < high_mid_rt
    end

    @testset "IntensityMassErrorModel tolerance heatmap grid uses fitted spread surface" begin
        mz_bias = _test_const_spline_f32(0f0; first=100f0, last=1000f0)
        intensity_bias = _test_const_spline_f32(0f0; first=0f0, last=20f0)
        spread = _test_const_spline_f32(2f0; first=0f0, last=20f0)
        mz_spread = _test_const_spline_f32(3f0; first=100f0, last=1000f0)
        rt_bias = _test_const_spline_f32(0f0; first=10f0, last=20f0)

        mem = IntensityMassErrorModel(
            mz_bias,
            intensity_bias,
            spread,
            rt_bias,
            _test_extrap_f32(mz_bias),
            _test_extrap_f32(intensity_bias),
            _test_extrap_f32(spread),
            _test_extrap_f32(rt_bias),
            2f0,
            0.02f0,
            mz_spread,
            _test_extrap_f32(mz_spread),
            1f0,
            0f0,
            0f0,
            0f0,
            0.02f0,
            10f0,
            20f0,
        )

        mz_grid, log2I_grid, tolerance_mda =
            _model_tolerance_heatmap_grid(mem, Float32[200, 800], Float64[4, 16];
                                          n_mz_bins=2, n_intensity_bins=2)

        @test mz_grid == [200.0, 800.0]
        @test log2I_grid == [4.0, 16.0]
        @test size(tolerance_mda) == (2, 2)
        @test all(tolerance_mda .≈ 2.0 * 2.0 * (log(2.0) / 0.6744897501960817) * 3.0)
    end

    @testset "IntensityMassErrorModel fits RT bias with m/z-like flexibility" begin
        n = 5000
        samples = Vector{Pioneer.MassErrSample}(undef, n)
        for i in 1:n
            frac = Float32((i - 1) / (n - 1))
            mz = 300f0 + 700f0 * Float32(((i * 37) % n) / (n - 1))
            log2I = 8f0 + 8f0 * Float32(((i * 53) % n) / (n - 1))
            rt = 5f0 + 30f0 * frac
            rt_bias = 0.003f0 * sin(8f0 * Float32(pi) * frac)
            samples[i] = Pioneer.MassErrSample(mz, mz + rt_bias, 2f0 ^ log2I, rt)
        end

        old_model = SimpleMassErrorModel(0f0, (20f0, 20f0))
        model = Pioneer.fit_intensity_mass_error_model(samples, old_model; k=1.96f0)

        @test model isa IntensityMassErrorModel
        @test length(model.rt_bias_spline.coeffs) >= length(model.mz_bias_spline.coeffs)
        @test model.rt_bias_spline.first >= 4.9f0
        @test model.rt_bias_spline.last <= 35.1f0
        @test model.rt_bias_extrap.lo_slope == 0f0
        @test model.rt_bias_extrap.hi_slope == 0f0
        @test 4.9f0 <= model.rt_bias_extrap.lo < model.rt_bias_extrap.hi <= 35.1f0
    end

end
