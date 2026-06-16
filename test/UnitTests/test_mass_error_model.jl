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
               LinearBiasPpmTolMassErrorModel,
               getRightTol, getLeftTol, getMassOffset, getMassCorrection,
               getCorrectedMz, getMzBounds, getMzBoundsReverse,
               getCorrectedMzAndBounds, laplace_log_density,
               get_default_top3_ll, MassErrSample

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

        @testset "4-arg forwarding ignores scan retention time" begin
            mem = SimpleMassErrorModel(2.5f0, (12.0f0, 9.0f0))
            mz = 750.0f0
            intensity = 42.0f0
            scan_rt = 31.5f0
            @test getCorrectedMz(mem, mz, intensity, scan_rt) == getCorrectedMz(mem, mz, intensity)
            @test getMzBoundsReverse(mem, mz, intensity, scan_rt) == getMzBoundsReverse(mem, mz, intensity)
            @test getCorrectedMzAndBounds(mem, mz, intensity, scan_rt) ==
                  getCorrectedMzAndBounds(mem, mz, intensity)
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

        @testset "4-arg forwarding ignores scan retention time" begin
            m = LinearDaMassErrorModel(0.001f0, 0.0f0, 0.002f0)
            mz = 1000.0f0
            intensity = 0.5f0
            scan_rt = 18.25f0
            @test getCorrectedMz(m, mz, intensity, scan_rt) == getCorrectedMz(m, mz, intensity)
            @test getMzBoundsReverse(m, mz, intensity, scan_rt) == getMzBoundsReverse(m, mz, intensity)
            @test getCorrectedMzAndBounds(m, mz, intensity, scan_rt) ==
                  getCorrectedMzAndBounds(m, mz, intensity)
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

        @testset "4-arg forwarding ignores scan retention time" begin
            m = LinearBiasPpmTolMassErrorModel(0.0005f0, 0.0f0, 6.0f0)
            mz = 1000.0f0
            intensity = 0.5f0
            scan_rt = 18.25f0
            @test getCorrectedMz(m, mz, intensity, scan_rt) == getCorrectedMz(m, mz, intensity)
            @test getMzBoundsReverse(m, mz, intensity, scan_rt) == getMzBoundsReverse(m, mz, intensity)
            @test getCorrectedMzAndBounds(m, mz, intensity, scan_rt) ==
                  getCorrectedMzAndBounds(m, mz, intensity)
        end
    end

    @testset "MassErrSample stores scan retention time" begin
        sample = MassErrSample(500.0f0, 500.002f0, 1234.0f0, 42.5f0)
        @test sample.theoretical_mz == 500.0f0
        @test sample.observed_mz == 500.002f0
        @test sample.intensity == 1234.0f0
        @test sample.rt == 42.5f0
    end

end
