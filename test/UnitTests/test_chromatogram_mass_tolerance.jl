using Test
using StaticArrays
using Pioneer

function tolerance_test_spline(value::Float32, lo::Float32, hi::Float32)
    coefficients = SVector{8, Float32}(value, 0, 0, 0, value, 0, 0, 0)
    return Pioneer.UniformSpline{8, Float32}(coefficients, 3, lo, hi, hi - lo)
end

@testset "Chromatogram mass tolerance" begin
    mz_bias = tolerance_test_spline(0.002f0, 100f0, 2000f0)
    int_bias = tolerance_test_spline(-0.001f0, 0f0, 24f0)
    spread = tolerance_test_spline(1.5f0, 0f0, 24f0)
    rt_bias = tolerance_test_spline(0.001f0, 0f0, 60f0)
    mz_spread = tolerance_test_spline(1f0, 100f0, 2000f0)
    extrap(s) = Pioneer.make_spline_extrap(s, s.first, s.last)
    calibrated = Pioneer.IntensityMassErrorModel(
        mz_bias, int_bias, spread, rt_bias,
        extrap(mz_bias), extrap(int_bias), extrap(spread), extrap(rt_bias),
        2f0, 0.0001f0, mz_spread, extrap(mz_spread), 1f0, 0f0, 0f0,
        -10f0, 0.003f0, 0f0, 60f0,
    )
    models = (
        calibrated,
        Pioneer.SimpleMassErrorModel(3f0, (20f0, 5f0)),
        Pioneer.LinearDaMassErrorModel(0.002f0, 1f-6, 0.006f0),
        Pioneer.LinearBiasPpmTolMassErrorModel(0.002f0, 1f-6, 10f0),
        Pioneer.ScoutCalibratedMassErrorModel(mz_bias, extrap(mz_bias),
            int_bias, extrap(int_bias), true, 0.006f0),
    )

    for model in models
        original_fields = ntuple(i -> getfield(model, i), fieldcount(typeof(model)))
        widened = @inferred Pioneer.chromatogram_mass_error_model(model)
        @test typeof(widened) === typeof(model)
        for mz in (100f0, 500f0, 1492.7925f0, 2000f0),
            intensity in (1f0, 1000f0, 582083.75f0), rt in (0f0, 35.49328f0, 60f0)
            center, lo, hi = Pioneer.getCorrectedMzAndBounds(model, mz, intensity, rt)
            corrected, low, high = Pioneer.getCorrectedMzAndBounds(widened, mz, intensity, rt)
            @test center === corrected
            @test low <= lo <= center <= hi <= high
            # Widths are rounded to Float32 after scaling, rather than expanding
            # an already rounded interval. Allow one m/z ULP for that rounding.
            @test abs((center - low) - 1.50f0 * (center - lo)) <= eps(center)
            @test abs((high - center) - 1.50f0 * (hi - center)) <= eps(center)
        end
        @test original_fields == ntuple(i -> getfield(model, i), fieldcount(typeof(model)))
    end

    widened = Pioneer.chromatogram_mass_error_model(calibrated)
    for field in fieldnames(typeof(calibrated))
        expected = field in (:k, :conservative_tol_da) ?
            getfield(calibrated, field) * 1.50f0 : getfield(calibrated, field)
        @test isequal(getfield(widened, field), expected)
    end
    @test Pioneer.laplace_log_density(widened, 500f0, 500.003f0, 1000f0) ==
          Pioneer.laplace_log_density(calibrated, 500f0, 500.003f0, 1000f0)

    @testset "Newly accepted fragment survives the coarse lookup" begin
        mz, intensity, rt = 500f0, 1000f0, 35f0
        center, _, hi = Pioneer.getCorrectedMzAndBounds(calibrated, mz, intensity, rt)
        target = center + 1.375f0 * (hi - center)
        @test target > hi
        corrected, lows, highs = Float32[], Float32[], Float32[]
        width = Ref(0f0)
        n = Pioneer.prepare_scan_peaks!(corrected, lows, highs, widened,
            Float32[mz], Float32[intensity], rt, width)
        @test lows[1] <= target <= highs[1]
        # The legacy conservative bound can be narrower than this valid match.
        @test abs(center - target) > Pioneer.getRightTol(widened)
        low, high = Pioneer.scan_match_window(target, width[])
        start = Pioneer.bsearch_hybrid(corrected, low, 1, n)
        @test Pioneer.scan_for_nearest_in_window(corrected, lows, highs,
            start, n, target, high)[1] == 1
        # MainSearch still uses the original model and rejects the same fragment.
        Pioneer.prepare_scan_peaks!(corrected, lows, highs, calibrated,
            Float32[mz], Float32[intensity], rt, width)
        @test target > highs[1]
    end
end
