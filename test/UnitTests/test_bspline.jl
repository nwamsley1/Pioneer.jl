# Unit tests for B-spline evaluation (src/utils/ML/libraryBSpline.jl).
#
# Validates prepared, unrolled de Boor evaluation against a recursive
# Cox-de Boor reference oracle.
#
# Run standalone: julia --project=. test/UnitTests/test_bspline.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

# ═══════════════════════════════════════════════════════════════════════════════
# Reference: recursive Cox-de Boor (the original implementation before de Boor)
# ═══════════════════════════════════════════════════════════════════════════════

function B_recursive(x::T, k::Int, i::Int, t::NTuple{N,T}) where {N,T<:AbstractFloat}
    if k == 0
        return T(t[i] <= x < t[i+1])
    end
    c1 = if t[i+k] == t[i]
        zero(T)
    else
        ((x - t[i]) / (t[i+k] - t[i])) * B_recursive(x, k-1, i, t)
    end
    c2 = if t[i+k+1] == t[i+1]
        zero(T)
    else
        ((t[i+k+1] - x) / (t[i+k+1] - t[i+1])) * B_recursive(x, k-1, i+1, t)
    end
    return c1 + c2
end

function splevl_recursive(x::T, knots::NTuple{N,T}, c::NTuple{M,T}, k::Int) where {M,N,T<:AbstractFloat}
    n = length(knots) - k - 1
    v = zero(T)
    for i in 1:n
        v += c[i] * B_recursive(x, k, i, knots)
    end
    return v
end

# Test data matching actual codebase types
const TEST_KNOTS = (6.0f0, 13.0f0, 20.0f0, 27.0f0, 34.0f0, 41.0f0, 48.0f0, 55.0f0)
const TEST_COEFFS = (1.6181915f-6, 7.382022f-6, 7.887343f-5, 0.00023642876f0)
const TEST_DEGREE = 3

eval_spline(x, knots, coefficients) = Pioneer.splevl_prepared(
    coefficients, Pioneer.prepare_spline_fractions(x, knots))

@testset "B-Spline (de Boor)" begin

# ═══════════════════════════════════════════════════════════════════════════════
# Interior points — mid-span evaluation
# ═══════════════════════════════════════════════════════════════════════════════

@testset "interior points" begin
    for x in Float32[10.0, 16.5, 23.5, 30.5, 37.5, 44.5, 51.5]
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        observed = eval_spline(x, TEST_KNOTS, TEST_COEFFS)
        @test observed ≈ v_ref atol=eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Exact knot positions
# ═══════════════════════════════════════════════════════════════════════════════

@testset "knot positions" begin
    for x in Float32.(TEST_KNOTS)
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        observed = eval_spline(x, TEST_KNOTS, TEST_COEFFS)
        @test observed ≈ v_ref atol=eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Boundary / outside range
# ═══════════════════════════════════════════════════════════════════════════════

@testset "boundaries and outside range" begin
    for x in Float32[5.0, 5.99, 6.01, 54.99, 55.0, 55.01, 60.0]
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        observed = eval_spline(x, TEST_KNOTS, TEST_COEFFS)
        @test observed ≈ v_ref atol=eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Points near knots (±eps)
# ═══════════════════════════════════════════════════════════════════════════════

@testset "near-knot points" begin
    for ki in TEST_KNOTS
        for offset in [eps(Float32), -eps(Float32), 10*eps(Float32)]
            x = Float32(ki) + offset
            v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
            observed = eval_spline(x, TEST_KNOTS, TEST_COEFFS)
            # Allow 2 ULP: different evaluation order causes ≤1.5 ULP rounding
            @test observed ≈ v_ref atol=2*eps(Float32)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Dense sweep across full range
# ═══════════════════════════════════════════════════════════════════════════════

@testset "dense sweep (200 points)" begin
    max_diff = 0.0f0
    for i in 0:200
        x = Float32(5.0 + 51.0 * i / 200)
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        observed = eval_spline(x, TEST_KNOTS, TEST_COEFFS)
        max_diff = max(max_diff, abs(v_ref - observed))
    end
    @test max_diff <= eps(Float32)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Various coefficient patterns
# ═══════════════════════════════════════════════════════════════════════════════

@testset "coefficient patterns" begin
    coeff_sets = [
        (1.0f0, 0.0f0, 0.0f0, 0.0f0),   # only first basis
        (0.0f0, 0.0f0, 0.0f0, 1.0f0),   # only last basis
        (1.0f0, 1.0f0, 1.0f0, 1.0f0),   # uniform
        (0.0f0, 1.0f0, 1.0f0, 0.0f0),   # middle only
        (1.0f-8, 1.0f0, 1.0f-8, 1.0f0), # large dynamic range
    ]
    for coeffs in coeff_sets
        max_diff = 0.0f0
        for i in 0:100
            x = Float32(5.0 + 50.0 * i / 100)
            v_ref = splevl_recursive(x, TEST_KNOTS, coeffs, TEST_DEGREE)
            observed = eval_spline(x, TEST_KNOTS, coeffs)
            max_diff = max(max_diff, abs(v_ref - observed))
        end
        @test max_diff <= 2 * eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Zero outside domain
# ═══════════════════════════════════════════════════════════════════════════════

@testset "zero outside domain" begin
    @test eval_spline(0.0f0, TEST_KNOTS, TEST_COEFFS) == 0.0f0
    @test eval_spline(100.0f0, TEST_KNOTS, TEST_COEFFS) == 0.0f0
    @test eval_spline(-10.0f0, TEST_KNOTS, TEST_COEFFS) == 0.0f0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Non-negative for non-negative coefficients
# ═══════════════════════════════════════════════════════════════════════════════

@testset "non-negative output" begin
    for i in 0:100
        x = Float32(5.0 + 50.0 * i / 100)
        @test eval_spline(x, TEST_KNOTS, TEST_COEFFS) >= 0.0f0
    end
end

@testset "prepared fragment spline evaluation" begin
    coeff_sets = (
        TEST_COEFFS,
        (1.0f0, 0.0f0, 0.0f0, 0.0f0),
        (0.0f0, 1.0f0, 1.0f0, 0.0f0),
        (1.0f-8, 1.0f0, 1.0f-8, 1.0f0),
    )
    for x in range(6f0, prevfloat(55f0), length=201)
        fractions = Pioneer.prepare_spline_fractions(x, TEST_KNOTS)
        for coeffs in coeff_sets
            expected = splevl_recursive(x, TEST_KNOTS, coeffs, TEST_DEGREE)
            @test Pioneer.splevl_prepared(coeffs, fractions) ≈ expected atol=2eps(Float32)
        end
    end

    outside = Pioneer.prepare_spline_fractions(100f0, TEST_KNOTS)
    @test Pioneer.splevl_prepared(TEST_COEFFS, outside) == 0f0

    repeated_knots = (10f0, 10f0, 10f0, 10f0, 50f0, 50f0, 50f0, 50f0)
    for x in (10f0, 20f0, 30f0, prevfloat(50f0)), coeffs in coeff_sets
        expected = splevl_recursive(x, repeated_knots, coeffs, TEST_DEGREE)
        prepared = Pioneer.prepare_spline_fractions(x, repeated_knots)
        @test Pioneer.splevl_prepared(coeffs, prepared) ≈ expected atol=2eps(Float32)
    end
end

@testset "prepared fragment intensity models" begin
    lookup = Pioneer.SplineFragmentLookup(
        Pioneer.SplineCompactFrag{4,Float32}[], UInt64[1], TEST_KNOTS)

    function test_model(model, mzs_and_charges)
        prepared_model = Pioneer.prepare_fragment_intensity_model(lookup, model)
        for (mz, charge) in mzs_and_charges
            prepared = Pioneer.getSplineData(lookup, prepared_model, charge, mz)
            expected = Pioneer.prepare_spline_fractions(model(mz, charge), TEST_KNOTS)
            @test isequal(prepared, expected)
        end
    end

    constant_model = Pioneer.PiecewiseNceModel(30f0)
    test_model(constant_model, ((400f0, UInt8(2)), (900f0, UInt8(4))))

    dynamic_model = Pioneer.PiecewiseNceModel(500f0, 0.01f0, 20f0, 25f0, 1f0)
    test_model(dynamic_model, ((400f0, UInt8(2)), (600f0, UInt8(3))))

    binned_model = Pioneer.BinnedMedianNceModel{Float32}(
        Float32[20, 22, 30, 32],
        (0x01, 0x03, 0x00, 0x00, 0x00, 0x00),
        (0x02, 0x02, 0x00, 0x00, 0x00, 0x00),
        (300f0, 300f0, 0f0, 0f0, 0f0, 0f0),
        (200f0, 200f0, 0f0, 0f0, 0f0, 0f0),
        27f0,
    )
    test_model(binned_model, (
        (250f0, UInt8(1)),
        (550f0, UInt8(1)),
        (250f0, UInt8(2)),
        (550f0, UInt8(2)),
        (400f0, UInt8(3)),
    ))

    empty_binned_model = Pioneer.BinnedMedianNceModel{Float32}(
        Float32[],
        ntuple(_ -> 0x00, 6),
        ntuple(_ -> 0x00, 6),
        ntuple(_ -> 0f0, 6),
        ntuple(_ -> 0f0, 6),
        27f0,
    )
    test_model(empty_binned_model, ((400f0, UInt8(2)),))
end

end # top-level testset
