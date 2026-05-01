# Unit tests for B-spline evaluation (src/utils/ML/libraryBSpline.jl).
#
# Validates the iterative de Boor algorithm (splevl) against the original
# recursive Cox-de Boor implementation that it replaced (22x speedup).
# The recursive version serves as a reference oracle.
#
# Run standalone: julia --project=. test/UnitTests/test_bspline.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using StaticArrays
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

@testset "B-Spline (de Boor)" begin

# ═══════════════════════════════════════════════════════════════════════════════
# Interior points — mid-span evaluation
# ═══════════════════════════════════════════════════════════════════════════════

@testset "interior points" begin
    for x in Float32[10.0, 16.5, 23.5, 30.5, 37.5, 44.5, 51.5]
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        v_deb = Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        @test v_deb ≈ v_ref atol=eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Exact knot positions
# ═══════════════════════════════════════════════════════════════════════════════

@testset "knot positions" begin
    for x in Float32.(TEST_KNOTS)
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        v_deb = Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        @test v_deb ≈ v_ref atol=eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Boundary / outside range
# ═══════════════════════════════════════════════════════════════════════════════

@testset "boundaries and outside range" begin
    for x in Float32[5.0, 5.99, 6.01, 54.99, 55.0, 55.01, 60.0]
        v_ref = splevl_recursive(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        v_deb = Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        @test v_deb ≈ v_ref atol=eps(Float32)
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
            v_deb = Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
            # Allow 2 ULP: different evaluation order causes ≤1.5 ULP rounding
            @test v_deb ≈ v_ref atol=2*eps(Float32)
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
        v_deb = Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE)
        max_diff = max(max_diff, abs(v_ref - v_deb))
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
            v_deb = Pioneer.splevl(x, TEST_KNOTS, coeffs, TEST_DEGREE)
            max_diff = max(max_diff, abs(v_ref - v_deb))
        end
        @test max_diff <= 2 * eps(Float32)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Zero outside domain
# ═══════════════════════════════════════════════════════════════════════════════

@testset "zero outside domain" begin
    @test Pioneer.splevl(0.0f0, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE) == 0.0f0
    @test Pioneer.splevl(100.0f0, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE) == 0.0f0
    @test Pioneer.splevl(-10.0f0, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE) == 0.0f0
end

# ═══════════════════════════════════════════════════════════════════════════════
# Non-negative for non-negative coefficients
# ═══════════════════════════════════════════════════════════════════════════════

@testset "non-negative output" begin
    for i in 0:100
        x = Float32(5.0 + 50.0 * i / 100)
        @test Pioneer.splevl(x, TEST_KNOTS, TEST_COEFFS, TEST_DEGREE) >= 0.0f0
    end
end

end # top-level testset
