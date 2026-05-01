# Lightweight tests for the small quad transmission model files and the
# IsotopeTraceType abstraction. These are mostly trivial getter / dispatch
# wrappers, but they ship at 0% coverage in the develop snapshot. The tests
# below execute every public function in each file at least once and check
# the documented invariants (NoQuad always returns 1.0; GeneralGauss decays
# from 1.0 at center to <1.0 outside).

using Pioneer: NoQuadModel, NoQuadFunction,
               GeneralGaussModel, GeneralGaussFunction, GeneralGaussQuadParams,
               getQuadTransmissionFunction, getQuadTransmissionBounds,
               getPrecMinBound, getPrecMaxBound,
               CombineTraces, SeperateTraces, seperateTraces, getPsmGroupbyCols

@testset "Quad models — basic dispatch coverage" begin

    @testset "NoQuadModel" begin
        m = NoQuadModel(1.0f0)   # 1 Da overhang on each side
        center = 500.0f0
        width  = 4.0f0

        bounds_low, bounds_high = getQuadTransmissionBounds(m, center, width)
        @test bounds_low  ≈ 498.0f0
        @test bounds_high ≈ 502.0f0
        # Bounds are NOT influenced by overhang
        @test (bounds_high - bounds_low) ≈ width

        f = getQuadTransmissionFunction(m, center, width)
        # The transmission FUNCTION incorporates the overhang on its mz limits
        @test getPrecMinBound(f) ≈ 497.0f0   # 498 - 1
        @test getPrecMaxBound(f) ≈ 503.0f0   # 502 + 1
        # NoQuadFunction returns 1.0 for any input mz, regardless of position
        for mz in (300.0f0, center, 700.0f0, -10.0f0)
            @test f(mz) === one(typeof(mz))
        end

        @testset "Float64 type stability" begin
            m64 = NoQuadModel(1.0)
            f64 = getQuadTransmissionFunction(m64, 500.0, 4.0)
            @test f64(500.0) === 1.0
            @test getPrecMinBound(f64) ≈ 497.0
        end
    end

    @testset "GeneralGaussModel" begin
        m = GeneralGaussModel(2.0f0, 0.0f0)   # b=2, no overhang
        center = 500.0f0
        width  = 4.0f0

        bounds_low, bounds_high = getQuadTransmissionBounds(m, center, width)
        @test bounds_low  ≈ 498.0f0
        @test bounds_high ≈ 502.0f0

        f = getQuadTransmissionFunction(m, center, width)
        @test getPrecMinBound(f) ≈ 498.0f0   # GG does NOT add overhang to bounds
        @test getPrecMaxBound(f) ≈ 502.0f0

        @testset "transmission peaks at center, decays outward" begin
            t_center = f(center)
            @test t_center ≈ 1.0f0
            t_off1   = f(center + 0.5f0)
            t_off2   = f(center + 1.5f0)
            t_far    = f(center + 5.0f0)
            # Strictly monotonic decay with |Δmz|
            @test 1.0f0 ≥ t_center > t_off1 > t_off2 > t_far ≥ 0.0f0
            # Mirror symmetry about center
            @test f(center - 1.5f0) ≈ f(center + 1.5f0)
        end

        @testset "overhang shifts the transmission curve outward" begin
            m_oh = GeneralGaussModel(2.0f0, 1.0f0)   # 1 Da overhang
            f_oh = getQuadTransmissionFunction(m_oh, center, width)
            # At the same off-center mz, more overhang → higher transmission
            for delta in (1.0f0, 2.0f0, 3.0f0)
                @test f_oh(center + delta) > f(center + delta)
            end
        end

        @testset "GeneralGaussQuadParams direct call (covers the abs path)" begin
            params = GeneralGaussQuadParams(2.0f0, 2.0f0, 0.0f0)
            # Symmetric: same value for ±x
            @test params(1.5f0) ≈ params(-1.5f0)
            @test params(0.0f0) ≈ 1.0f0
        end
    end

    @testset "IsotopeTraceType" begin
        @testset "CombineTraces" begin
            ct = CombineTraces(0.25f0)
            @test seperateTraces(ct) == false
            @test getPsmGroupbyCols(ct) == [:precursor_idx]
            @test ct.min_ft == 0.25f0
        end

        @testset "SeperateTraces" begin
            st = SeperateTraces()
            @test seperateTraces(st) == true
            @test getPsmGroupbyCols(st) == [:precursor_idx, :isotopes_captured]
        end
    end

end
