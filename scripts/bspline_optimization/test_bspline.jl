#=
B-Spline Optimization: Iterative de Boor vs Recursive Cox-de Boor
=================================================================
Tests correctness and benchmarks performance of the iterative de Boor algorithm
as a drop-in replacement for the recursive B() + splevl() in libraryBSpline.jl.

Usage:
    julia --project=. scripts/bspline_optimization/test_bspline.jl
=#

# BenchmarkTools is in the global env; ensure it's on LOAD_PATH
push!(LOAD_PATH, "@v#.$(VERSION.major).$(VERSION.minor)")
push!(LOAD_PATH, "@stdlib")

using BenchmarkTools
using StaticArrays
using FastGaussQuadrature
using Printf

# ============================================================
# Original recursive implementation (from libraryBSpline.jl)
# ============================================================

function B_recursive(x::T, k::Int, i::Int, t::NTuple{N,T}) where {N,T<:AbstractFloat}
    if k == 0
        return T(t[i] ≤ x < t[i+1])
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

function splint_recursive(knots::NTuple{N,T}, c::NTuple{M,T}, d::Int,
                          gqx::SVector{20,T}, gqw::SVector{20,T}) where {M,N,T<:AbstractFloat}
    i_eval = zero(T)
    for i in 1:20
        i_eval += splevl_recursive(gqx[i], knots, c, d) * gqw[i]
    end
    return i_eval
end

# ============================================================
# Iterative de Boor implementation
# ============================================================
#
# The de Boor algorithm evaluates S(x) = Σ c_i * B_i^k(x) in O(k²) operations
# by finding the knot span containing x and iteratively combining k+1 coefficients,
# instead of evaluating each B_i^k separately via O(2^k) recursive calls.
#
# For k=3, N=8, M=4: 6 multiply-adds vs 32 recursive calls.
#
# Key insight for out-of-range spans: coefficients outside [1,M] are zero-padded,
# and the corresponding knot accesses are clamped. Since zero coefficients propagate
# through the triangular iterations, the clamped knot values are only ever multiplied
# by zero, producing correct results for all x.

"""
    splevl_deboor(x, knots, c, k)

Evaluate B-spline at point x using the iterative de Boor algorithm.
Drop-in replacement for splevl(). Assumes k=3 (cubic).

Algorithm: For knot span j where knots[j] ≤ x < knots[j+1],
initialize d[1..4] = c[j-3..j] (zero-padded), then iterate:
  for r = 1:3, for p = 4:-1:(r+1):
    α = (x - t[j+p-4]) / (t[j+p-r] - t[j+p-4])
    d[p] = (1-α)*d[p-1] + α*d[p]
Result = d[4].
"""
@inline function splevl_deboor(x::T, knots::NTuple{N,T}, c::NTuple{M,T}, k::Int) where {M,N,T<:AbstractFloat}
    # Find knot span: j where knots[j] ≤ x < knots[j+1]
    j = 0
    @inbounds for idx in 1:(N-1)
        if knots[idx] ≤ x < knots[idx+1]
            j = idx
            break
        end
    end

    # x outside all knot spans → all basis functions are zero
    j == 0 && return zero(T)

    # Safe access helpers:
    # - Coefficients: zero outside [1, M]
    # - Knots: clamp to [1, N] (only matters when multiplied by zero)
    @inline _getc(i) = (1 ≤ i ≤ M) ? @inbounds(c[i]) : zero(T)
    @inline _getk(i) = @inbounds knots[clamp(i, 1, N)]

    # Initialize d[1..4] with active coefficients
    d1 = _getc(j - 3)
    d2 = _getc(j - 2)
    d3 = _getc(j - 1)
    d4 = _getc(j)

    # ---- r = 1: three updates ----
    denom = _getk(j + 3) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    denom = _getk(j + 2) - _getk(j - 1)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 1)) / denom)
    d3 = (one(T) - α) * d2 + α * d3

    denom = _getk(j + 1) - _getk(j - 2)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 2)) / denom)
    d2 = (one(T) - α) * d1 + α * d2

    # ---- r = 2: two updates ----
    denom = _getk(j + 2) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    denom = _getk(j + 1) - _getk(j - 1)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 1)) / denom)
    d3 = (one(T) - α) * d2 + α * d3

    # ---- r = 3: one update ----
    denom = _getk(j + 1) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    return d4
end

"""
    splint_deboor(knots, c, d, gqx, gqw)

Integrate B-spline using Gauss-Legendre quadrature with the de Boor evaluator.
Drop-in replacement for splint().
"""
function splint_deboor(knots::NTuple{N,T}, c::NTuple{M,T}, d::Int,
                       gqx::SVector{20,T}, gqw::SVector{20,T}) where {M,N,T<:AbstractFloat}
    i_eval = zero(T)
    @inbounds for i in 1:20
        i_eval += splevl_deboor(gqx[i], knots, c, d) * gqw[i]
    end
    return i_eval
end

# ============================================================
# Test data (matches actual codebase types)
# ============================================================

const KNOTS = (6.0f0, 13.0f0, 20.0f0, 27.0f0, 34.0f0, 41.0f0, 48.0f0, 55.0f0)
const COEFFS = (1.6181915f-6, 7.382022f-6, 7.887343f-5, 0.00023642876f0)
const DEGREE = 3

function make_quadrature(::Type{T}, x0, x1) where T
    x, w = gausslegendre(20)
    ws = (x1 - x0) / 2
    for i in 1:20
        x[i] = x0 + (x[i] + 1) * (x1 - x0) / 2
        w[i] = ws * w[i]
    end
    return SVector{20}(T.(x)), SVector{20}(T.(w))
end

# ============================================================
# Correctness tests
# ============================================================

function test_correctness()
    println("=" ^ 60)
    println("CORRECTNESS TESTS")
    println("=" ^ 60)

    knots = KNOTS
    c = COEFFS
    k = DEGREE
    all_pass = true

    # --- Test 1: Interior points across all knot spans ---
    println("\n--- Test 1: Interior points (mid-span) ---")
    test_points = Float32[10.0, 16.5, 23.5, 30.5, 37.5, 44.5, 51.5]
    max_diff = 0.0f0
    for x in test_points
        v_rec = splevl_recursive(x, knots, c, k)
        v_deb = splevl_deboor(x, knots, c, k)
        diff = abs(v_rec - v_deb)
        max_diff = max(max_diff, diff)
        if diff > eps(Float32)
            @printf("  FAIL x=%.1f: recursive=%.10e, deboor=%.10e, diff=%.2e\n", x, v_rec, v_deb, diff)
            all_pass = false
        end
    end
    @printf("  Max diff: %.2e (threshold: %.2e) %s\n",
            max_diff, eps(Float32), max_diff ≤ eps(Float32) ? "PASS" : "FAIL")

    # --- Test 2: Exact knot positions ---
    println("\n--- Test 2: Exact knot positions ---")
    max_diff = 0.0f0
    for x in Float32.(knots)
        v_rec = splevl_recursive(x, knots, c, k)
        v_deb = splevl_deboor(x, knots, c, k)
        diff = abs(v_rec - v_deb)
        max_diff = max(max_diff, diff)
        if diff > eps(Float32)
            @printf("  FAIL x=%.1f: recursive=%.10e, deboor=%.10e, diff=%.2e\n", x, v_rec, v_deb, diff)
            all_pass = false
        end
    end
    @printf("  Max diff: %.2e %s\n", max_diff, max_diff ≤ eps(Float32) ? "PASS" : "FAIL")

    # --- Test 3: Boundary and endpoint behavior ---
    println("\n--- Test 3: Boundaries (first/last knot, outside range) ---")
    boundary_points = Float32[5.0, 5.99, 6.0, 6.01, 54.99, 55.0, 55.01, 60.0]
    max_diff = 0.0f0
    for x in boundary_points
        v_rec = splevl_recursive(x, knots, c, k)
        v_deb = splevl_deboor(x, knots, c, k)
        diff = abs(v_rec - v_deb)
        max_diff = max(max_diff, diff)
        if diff > eps(Float32)
            @printf("  FAIL x=%.2f: recursive=%.10e, deboor=%.10e, diff=%.2e\n", x, v_rec, v_deb, diff)
            all_pass = false
        end
    end
    @printf("  Max diff: %.2e %s\n", max_diff, max_diff ≤ eps(Float32) ? "PASS" : "FAIL")

    # --- Test 4: Near-knot points (within eps) ---
    println("\n--- Test 4: Points near knots (±eps) ---")
    max_diff = 0.0f0
    for ki in knots
        for offset in [eps(Float32), -eps(Float32), 10*eps(Float32), -10*eps(Float32)]
            x = Float32(ki) + offset
            v_rec = splevl_recursive(x, knots, c, k)
            v_deb = splevl_deboor(x, knots, c, k)
            diff = abs(v_rec - v_deb)
            max_diff = max(max_diff, diff)
            if diff > 2 * eps(Float32)  # Allow slightly more tolerance near boundaries
                @printf("  FAIL x=%.10f: recursive=%.10e, deboor=%.10e, diff=%.2e\n", x, v_rec, v_deb, diff)
                all_pass = false
            end
        end
    end
    @printf("  Max diff: %.2e (threshold: %.2e) %s\n",
            max_diff, 2*eps(Float32), max_diff ≤ 2*eps(Float32) ? "PASS" : "FAIL")

    # --- Test 5: Dense sweep across full range ---
    println("\n--- Test 5: Dense sweep (1000 points, 5.0 to 56.0) ---")
    max_diff = 0.0f0
    worst_x = 0.0f0
    n_points = 1000
    for i in 0:n_points
        x = Float32(5.0 + (56.0 - 5.0) * i / n_points)
        v_rec = splevl_recursive(x, knots, c, k)
        v_deb = splevl_deboor(x, knots, c, k)
        diff = abs(v_rec - v_deb)
        if diff > max_diff
            max_diff = diff
            worst_x = x
        end
    end
    pass = max_diff ≤ eps(Float32)
    if !pass
        all_pass = false
    end
    @printf("  Max diff: %.2e at x=%.4f %s\n", max_diff, worst_x, pass ? "PASS" : "FAIL")

    # --- Test 6: Right endpoint (x == last knot) ---
    println("\n--- Test 6: Right endpoint (x == t[end]) ---")
    x = Float32(55.0)
    v_rec = splevl_recursive(x, knots, c, k)
    v_deb = splevl_deboor(x, knots, c, k)
    @printf("  x=55.0: recursive=%.10e, deboor=%.10e\n", v_rec, v_deb)
    if v_rec == v_deb
        println("  PASS (both return same value)")
    else
        diff = abs(v_rec - v_deb)
        @printf("  FAIL diff=%.2e\n", diff)
        all_pass = false
    end

    # --- Test 7: Different coefficient sets ---
    println("\n--- Test 7: Various coefficient sets ---")
    coeff_sets = [
        (1.0f0, 0.0f0, 0.0f0, 0.0f0),   # Only first basis
        (0.0f0, 0.0f0, 0.0f0, 1.0f0),   # Only last basis
        (1.0f0, 1.0f0, 1.0f0, 1.0f0),   # Uniform
        (0.0f0, 1.0f0, 1.0f0, 0.0f0),   # Middle only
        (1.0f-8, 1.0f0, 1.0f-8, 1.0f0), # Large dynamic range
    ]
    # Allow 2 ULP: different evaluation order causes ≤1.5 ULP rounding differences
    coeff_threshold = 2 * eps(Float32)
    max_diff_all = 0.0f0
    for (ci, coeffs) in enumerate(coeff_sets)
        max_diff = 0.0f0
        for i in 0:100
            x = Float32(5.0 + 50.0 * i / 100)
            v_rec = splevl_recursive(x, knots, coeffs, k)
            v_deb = splevl_deboor(x, knots, coeffs, k)
            diff = abs(v_rec - v_deb)
            max_diff = max(max_diff, diff)
        end
        max_diff_all = max(max_diff_all, max_diff)
        @printf("  Coeff set %d: max diff = %.2e %s\n", ci, max_diff,
                max_diff ≤ coeff_threshold ? "PASS" : "FAIL")
        if max_diff > coeff_threshold
            all_pass = false
        end
    end

    # --- Test 8: Integration (splint) ---
    println("\n--- Test 8: Integration (splint) ---")
    gqx, gqw = make_quadrature(Float32, 20.0, 40.0)
    v_rec = splint_recursive(knots, c, k, gqx, gqw)
    v_deb = splint_deboor(knots, c, k, gqx, gqw)
    diff = abs(v_rec - v_deb)
    @printf("  recursive: %.10e\n", v_rec)
    @printf("  deboor:    %.10e\n", v_deb)
    @printf("  diff:      %.2e\n", diff)
    # Integration accumulates 20 evaluations, so allow proportionally more tolerance
    int_threshold = 20 * eps(Float32)
    pass = diff ≤ int_threshold
    if !pass
        all_pass = false
    end
    @printf("  Threshold: %.2e %s\n", int_threshold, pass ? "PASS" : "FAIL")

    # --- Test 9: Integration with different ranges ---
    println("\n--- Test 9: Integration with different quadrature ranges ---")
    ranges = [(6.0, 55.0), (20.0, 40.0), (27.0, 34.0), (30.0, 50.0)]
    for (x0, x1) in ranges
        gqx_r, gqw_r = make_quadrature(Float32, x0, x1)
        v_rec = splint_recursive(knots, c, k, gqx_r, gqw_r)
        v_deb = splint_deboor(knots, c, k, gqx_r, gqw_r)
        diff = abs(v_rec - v_deb)
        pass = diff ≤ int_threshold
        if !pass
            all_pass = false
        end
        @printf("  [%.0f, %.0f]: diff=%.2e %s\n", x0, x1, diff, pass ? "PASS" : "FAIL")
    end

    println("\n" * "=" ^ 60)
    if all_pass
        println("ALL CORRECTNESS TESTS PASSED")
    else
        println("SOME TESTS FAILED")
    end
    println("=" ^ 60)

    return all_pass
end

# ============================================================
# Benchmarks
# ============================================================

function run_benchmarks()
    println("\n" * "=" ^ 60)
    println("BENCHMARKS")
    println("=" ^ 60)

    knots = KNOTS
    c = COEFFS
    k = DEGREE

    # --- splevl: single point evaluation ---
    x = 37.0f0  # Same as the commented example in libraryBSpline.jl
    println("\n--- splevl at x=$x ---")

    print("  recursive: ")
    b_rec = @benchmark splevl_recursive($x, $knots, $c, $k)
    show(stdout, MIME("text/plain"), b_rec)
    println()

    print("  deboor:    ")
    b_deb = @benchmark splevl_deboor($x, $knots, $c, $k)
    show(stdout, MIME("text/plain"), b_deb)
    println()

    speedup = median(b_rec).time / median(b_deb).time
    @printf("  Speedup: %.1fx\n", speedup)

    # --- splevl at different knot spans ---
    println("\n--- splevl at various x values ---")
    for x_val in Float32[10.0, 23.5, 30.5, 37.5, 44.5, 51.5]
        t_rec = @belapsed splevl_recursive($x_val, $knots, $c, $k)
        t_deb = @belapsed splevl_deboor($x_val, $knots, $c, $k)
        @printf("  x=%.1f: recursive=%.1fns, deboor=%.1fns, speedup=%.1fx\n",
                x_val, t_rec*1e9, t_deb*1e9, t_rec/t_deb)
    end

    # --- splint: integration ---
    println("\n--- splint (20-point Gauss-Legendre) ---")
    gqx, gqw = make_quadrature(Float32, 20.0, 40.0)

    print("  recursive: ")
    b_rec = @benchmark splint_recursive($knots, $c, $k, $gqx, $gqw)
    show(stdout, MIME("text/plain"), b_rec)
    println()

    print("  deboor:    ")
    b_deb = @benchmark splint_deboor($knots, $c, $k, $gqx, $gqw)
    show(stdout, MIME("text/plain"), b_deb)
    println()

    speedup = median(b_rec).time / median(b_deb).time
    @printf("  Speedup: %.1fx\n", speedup)

    # --- Batch evaluation (simulates hot loop) ---
    println("\n--- Batch: 10000 splevl calls (uniform x in [20, 45]) ---")
    xs = Float32[20.0f0 + 25.0f0 * i / 10000 for i in 0:9999]

    t_rec = @belapsed begin
        s = 0.0f0
        @inbounds for x in $xs
            s += splevl_recursive(x, $knots, $c, $k)
        end
        s
    end

    t_deb = @belapsed begin
        s = 0.0f0
        @inbounds for x in $xs
            s += splevl_deboor(x, $knots, $c, $k)
        end
        s
    end

    @printf("  recursive: %.2f μs (%.1f ns/eval)\n", t_rec*1e6, t_rec*1e9/10000)
    @printf("  deboor:    %.2f μs (%.1f ns/eval)\n", t_deb*1e6, t_deb*1e9/10000)
    @printf("  Speedup: %.1fx\n", t_rec/t_deb)

    println("\n" * "=" ^ 60)
end

# ============================================================
# Main
# ============================================================

function main()
    println("B-Spline Optimization: de Boor vs Recursive")
    println("Knots:  $KNOTS")
    println("Coeffs: $COEFFS")
    println("Degree: $DEGREE")

    passed = test_correctness()

    if passed
        run_benchmarks()
    else
        println("\nSkipping benchmarks due to test failures.")
        exit(1)
    end
end

main()
