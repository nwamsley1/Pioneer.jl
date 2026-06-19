# Unit tests for run_fused! per-precursor pre-match filter helpers.
#
# Target functions (src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl):
#   passes_irt_filter
#   quad_window_with_iso_bounds
#   passes_prec_mz_filter
#
# Run: julia --project=. test/UnitTests/test_fused_prec_filters.jl

using Test
using Pioneer

const passes_irt_filter           = Pioneer.passes_irt_filter
const quad_window_with_iso_bounds = Pioneer.quad_window_with_iso_bounds
const passes_prec_mz_filter       = Pioneer.passes_prec_mz_filter
const iso_mz_for                  = Pioneer.iso_mz_for
const in_frag_mz_window           = Pioneer.in_frag_mz_window
const conservative_half_width     = Pioneer.conservative_half_width
const match_window                = Pioneer.match_window
const compute_ppm_err             = Pioneer.compute_ppm_err
const push_match!                 = Pioneer.push_match!
const push_miss!                  = Pioneer.push_miss!
const FusedScratch                = Pioneer.FusedScratch
# Use the module-qualified Float64 constant below so this file works both
# standalone and when included after test_setup.jl.
const C13_C12_MASS_DIFF_F32 = Pioneer.C13_C12_MASS_DIFF_F32

# A minimal QuadTransmissionFunction stand-in: GeneralGaussModel exists in
# Pioneer and supports getPrecMinBound / getPrecMaxBound. Using a fixed-bounds
# wrapper keeps the test deterministic regardless of model details.
struct FixedQuadFunc <: Pioneer.QuadTransmissionFunction
    lo::Float32
    hi::Float32
end
Pioneer.getPrecMinBound(q::FixedQuadFunc) = q.lo
Pioneer.getPrecMaxBound(q::FixedQuadFunc) = q.hi

@testset "fused per-precursor filters" begin

    @testset "passes_irt_filter" begin
        @test  passes_irt_filter(50.0f0, 50.0f0, 5.0f0)   # exact center
        @test  passes_irt_filter(45.0f0, 50.0f0, 5.0f0)   # at -tol boundary
        @test  passes_irt_filter(55.0f0, 50.0f0, 5.0f0)   # at +tol boundary
        @test !passes_irt_filter(44.999f0, 50.0f0, 5.0f0)
        @test !passes_irt_filter(55.001f0, 50.0f0, 5.0f0)
    end

    @testset "quad_window_with_iso_bounds" begin
        # Quad bounds [500, 510], iso bounds (-1, 1), charge 2.
        qfunc = FixedQuadFunc(500.0f0, 510.0f0)
        low, high = quad_window_with_iso_bounds(qfunc, 2, (-1, 1))
        @test low  === 500.0f0 - Float32(Pioneer.C13_C12_MASS_DIFF * (-1) / 2)
        @test high === 510.0f0 + Float32(Pioneer.C13_C12_MASS_DIFF *  1  / 2)

        # Charge dependence: the C13-C12 m/z offset is smaller at higher charge,
        # so the offset relative to literal quad bounds shrinks with charge.
        for charge in (1, 2, 3, 4)
            l, h = quad_window_with_iso_bounds(qfunc, charge, (-1, 1))
            @test l === 500.0f0 - Float32(Pioneer.C13_C12_MASS_DIFF * (-1) / charge)
            @test h === 510.0f0 + Float32(Pioneer.C13_C12_MASS_DIFF *   1  / charge)
        end

        # Asymmetric bounds (-2, 0) → first/last applied independently.
        low_a, high_a = quad_window_with_iso_bounds(qfunc, 2, (-2, 0))
        @test low_a  === 500.0f0 - Float32(Pioneer.C13_C12_MASS_DIFF * (-2) / 2)
        @test high_a === 510.0f0 + Float32(Pioneer.C13_C12_MASS_DIFF *   0  / 2)
    end

    @testset "passes_prec_mz_filter" begin
        @test  passes_prec_mz_filter(505.0f0, 500.0f0, 510.0f0)
        @test  passes_prec_mz_filter(500.0f0, 500.0f0, 510.0f0)   # boundary
        @test  passes_prec_mz_filter(510.0f0, 500.0f0, 510.0f0)   # boundary
        @test !passes_prec_mz_filter(499.9f0, 500.0f0, 510.0f0)
        @test !passes_prec_mz_filter(510.1f0, 500.0f0, 510.0f0)
    end

    @testset "iso_mz_for" begin
        # iso_idx=0 returns monoisotopic m/z exactly.
        @test iso_mz_for(500.0f0, 0, 1f0)        === 500.0f0
        @test iso_mz_for(800.123f0, 0, 0.5f0)    === 800.123f0

        # +1 isotope at charge 1 / 2 / 3.
        @test iso_mz_for(500.0f0, 1, 1f0)        ≈ 500.0f0 + C13_C12_MASS_DIFF_F32
        @test iso_mz_for(500.0f0, 1, 1f0/2f0)    ≈ 500.0f0 + C13_C12_MASS_DIFF_F32/2f0
        @test iso_mz_for(500.0f0, 1, 1f0/3f0)    ≈ 500.0f0 + C13_C12_MASS_DIFF_F32/3f0

        # +2 isotope at charge 2 → exactly one C13-C12 spacing above monoisotopic.
        @test iso_mz_for(500.0f0, 2, 1f0/2f0)    ≈ 500.0f0 + C13_C12_MASS_DIFF_F32

        # bit-identical to the original expression.
        for frag_mz in (123.45f0, 800.5f0, 1500.001f0),
            iso_idx in 0:3,
            charge in (1, 2, 3, 4)
            inv = 1f0 / Float32(charge)
            ref = Float32(frag_mz + Float32(iso_idx) * C13_C12_MASS_DIFF_F32 * inv)
            @test iso_mz_for(frag_mz, iso_idx, inv) === ref
        end
    end

    @testset "in_frag_mz_window" begin
        @test  in_frag_mz_window(500.0f0, 100.0f0, 2000.0f0)
        @test  in_frag_mz_window(100.0f0, 100.0f0, 2000.0f0)   # boundary
        @test  in_frag_mz_window(2000.0f0, 100.0f0, 2000.0f0)  # boundary
        @test !in_frag_mz_window(99.999f0, 100.0f0, 2000.0f0)
        @test !in_frag_mz_window(2000.001f0, 100.0f0, 2000.0f0)
    end

    @testset "conservative_half_width" begin
        # SimpleMassErrorModel(mass_offset, (left_tol, right_tol)) — both tols
        # in ppm. getMzBounds returns (mass - left_tol*ppm, mass + right_tol*ppm)
        # roughly; the exposed `right` side is what conservative_half_width uses.
        mem = Pioneer.SimpleMassErrorModel(0.0f0, (10.0f0, 10.0f0))   # ±10 ppm
        @test conservative_half_width(mem, 500.0f0)  ≈ 500.0f0 * 10f0 * 1f-6  atol=1f-3
        @test conservative_half_width(mem, 1000.0f0) ≈ 1000.0f0 * 10f0 * 1f-6 atol=1f-3
        # bit-identical to the inlined original (`hi - frag_mz` after getMzBounds).
        for mz in (123.45f0, 800.0f0, 1500.0f0)
            _, hi = Pioneer.getMzBounds(mem, mz)
            @test conservative_half_width(mem, mz) === hi - mz
        end
    end

    @testset "match_window" begin
        @test match_window(500.0f0, 0.005f0) === (499.995f0, 500.005f0)
        @test match_window(1000.0f0, 0.0f0)   === (1000.0f0, 1000.0f0)  # zero tol
        # By construction: lo === iso_mz - hw, hi === iso_mz + hw. Avoid
        # asserting `iso_mz - lo === hw` — Float32 round-trip can be sub-ULP off.
        for iso_mz in (100.0f0, 500.0f0, 1500.0f0), hw in (0.001f0, 0.01f0, 0.1f0)
            lo, hi = match_window(iso_mz, hw)
            @test lo === iso_mz - hw
            @test hi === iso_mz + hw
        end
    end

    @testset "compute_ppm_err" begin
        # Theoretical 500.000, observed 500.005 → ≈-10 ppm (observed higher).
        # 500.005f0 isn't exactly representable, so allow a 0.05 ppm tol.
        @test compute_ppm_err(500.000f0, 500.005f0) ≈ -10.0f0  atol=0.05
        # Observed lower → positive ppm.
        @test compute_ppm_err(500.005f0, 500.000f0) ≈  10.0f0  atol=0.05
        # Zero error.
        @test compute_ppm_err(800.0f0, 800.0f0) === 0.0f0
        # Bit-identical to inline form for representative values.
        for theo in (123.45f0, 500.0f0, 800.5f0, 1500.001f0)
            for delta_ppm in (-20f0, -5f0, 0f0, 5f0, 20f0)
                obs = theo - theo * delta_ppm * 1f-6
                @test compute_ppm_err(theo, obs) === (theo - obs) / (theo * 1f-6)
            end
        end
    end

    @testset "push_match! / push_miss!" begin
        # Start with capacity = 4 so we can exercise grow.
        s = FusedScratch(4)
        @test s.n == 0 && s.miss_n == 0

        push_match!(s, 7, 1.5f0, 1000.0f0, 0, 3)
        @test s.n == 1
        @test s.row[1]     === UInt32(7)
        @test s.nzval[1]   === 1.5f0
        @test s.x[1]       === 1000.0f0
        @test s.isotope[1] === UInt8(0)
        @test s.rank[1]    === UInt8(3)

        push_miss!(s, 0.25f0, 2, 4)
        @test s.miss_n == 1
        @test s.miss_nzval[1]   === 0.25f0
        @test s.miss_isotope[1] === UInt8(2)
        @test s.miss_rank[1]    === UInt8(4)

        # Force grow on the match buffer (capacity 4 → push 5).
        for k in 2:5
            push_match!(s, k * 10, Float32(k), Float32(k * 100), k - 1, k)
        end
        @test s.n == 5
        @test length(s.row) >= 5
        @test s.row[5]     === UInt32(50)
        @test s.isotope[5] === UInt8(4)
        @test s.rank[5]    === UInt8(5)

        # Independent buffers — pushing matches doesn't change miss state.
        @test s.miss_n == 1

        # Force grow on the miss buffer.
        for k in 2:6
            push_miss!(s, Float32(k) * 0.1f0, k, k + 1)
        end
        @test s.miss_n == 6
        @test length(s.miss_nzval) >= 6
        @test s.miss_nzval[6]   === 0.6f0
        @test s.miss_isotope[6] === UInt8(6)
        @test s.miss_rank[6]    === UInt8(7)
    end

    @testset "bit-identical to inline original" begin
        # Reproduce the original arithmetic verbatim and confirm the helper
        # matches it exactly. This guards against accidental refactors of
        # the float ordering (Float32(a*b/c) vs Float32(a*b)*Float32(1/c)).
        qfunc = FixedQuadFunc(498.7f0, 511.3f0)
        for prec_charge in (1, 2, 3, 4), bounds in ((-1,1), (-2,0), (0,2), (-3,3))
            ref_lo = Pioneer.getPrecMinBound(qfunc) -
                     Float32(Pioneer.C13_C12_MASS_DIFF * first(bounds) / prec_charge)
            ref_hi = Pioneer.getPrecMaxBound(qfunc) +
                     Float32(Pioneer.C13_C12_MASS_DIFF * last(bounds)  / prec_charge)
            lo, hi = quad_window_with_iso_bounds(qfunc, prec_charge, bounds)
            @test lo === ref_lo
            @test hi === ref_hi
        end
    end

end
