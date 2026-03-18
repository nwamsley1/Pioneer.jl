# Unit tests for solveOLS! — non-negative coordinate descent OLS solver
#
# Target: Pioneer.solveOLS! (spectralLinearRegression.jl:305-371)
# Solves min ‖Hw − x‖² s.t. w ≥ 0
#
# Run: julia --project=. test/UnitTests/test_solveOLS.jl

using Test
using Pioneer
using LinearAlgebra

# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------
const SA            = Pioneer.SparseArray
const sortSparse!   = Pioneer.sortSparse!
const initResiduals! = Pioneer.initResiduals!
const solveOLS!     = Pioneer.solveOLS!

# ---------------------------------------------------------------------------
# Helper: build SparseArray from (col, row, nzval, x, matched, isotope) tuples
# ---------------------------------------------------------------------------
"""
Build a SparseArray from a vector of `(col, row, nzval, x, matched, isotope)` tuples.
Calls `sortSparse!` to set up colptr, then fixes sa.m in case the first sorted
entry has the max rowval (sortSparse! only checks entries 2..n_vals).
"""
function build_test_sparse(entries)
    n = length(entries)
    sa = SA(UInt32(max(n + 10, 50)))
    sa.n_vals = n
    for (i, (col, row, nz, xv, m, iso)) in enumerate(entries)
        sa.colval[i] = UInt16(col)
        sa.rowval[i] = UInt32(row)
        sa.nzval[i]  = Float32(nz)
        sa.x[i]      = Float32(xv)
        sa.matched[i] = m
        sa.isotope[i] = UInt8(iso)
    end
    sortSparse!(sa)
    # Fix sa.m: sortSparse! may miss the first entry's rowval
    for i in 1:sa.n_vals
        sa.m = max(sa.m, Int64(sa.rowval[i]))
    end
    return sa
end

"""Solve from scratch: zero-initialize weights, compute residuals, run OLS."""
function solve_fresh(H; max_iter=200, tol=Float32(1e-6))
    w = zeros(Float32, H.n)
    r = zeros(Float32, H.m)
    colnorm2 = zeros(Float32, H.n)
    initResiduals!(r, H, w)
    result = solveOLS!(H, r, w, colnorm2, max_iter, tol)
    return w, r, colnorm2, result
end

"""Convert SparseArray to dense Matrix{Float32} (via colptr iteration)."""
function to_dense_colptr(H)
    M = zeros(Float32, H.m, H.n)
    for col in 1:H.n
        for i in H.colptr[col]:(H.colptr[col+1]-1)
            M[H.rowval[i], col] += H.nzval[i]
        end
    end
    return M
end

# ═══════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "solveOLS!" begin

    # -------------------------------------------------------------------
    @testset "1.1 Identity matrix" begin
        # H = I(3), x = [10, 20, 30] → w = [10, 20, 30]
        entries = [
            (1, 1, 1.0, 10.0, true, 0),
            (2, 2, 1.0, 20.0, true, 0),
            (3, 3, 1.0, 30.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w, r, _, result = solve_fresh(H)

        @test w ≈ Float32[10, 20, 30] atol=1e-4
    end

    # -------------------------------------------------------------------
    @testset "1.2 Diagonal matrix" begin
        # H = diag([2, 3, 5]), x = [10, 15, 25] → w = [5, 5, 5]
        entries = [
            (1, 1, 2.0, 10.0, true, 0),
            (2, 2, 3.0, 15.0, true, 0),
            (3, 3, 5.0, 25.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w, _, _, _ = solve_fresh(H)

        @test w ≈ Float32[5, 5, 5] atol=1e-4
    end

    # -------------------------------------------------------------------
    @testset "1.3 Two overlapping columns" begin
        # H = [1 0; 1 1; 0 1], x = [2, 5, 3] → w = [2, 3]
        # H^T H = [2 1; 1 2], H^T x = [7, 8]
        entries = [
            (1, 1, 1.0, 2.0, true, 0),   # col 1, row 1
            (1, 2, 1.0, 5.0, true, 0),   # col 1, row 2
            (2, 2, 1.0, 5.0, true, 0),   # col 2, row 2 (same x as col 1 row 2)
            (2, 3, 1.0, 3.0, true, 0),   # col 2, row 3
        ]
        H = build_test_sparse(entries)
        w, _, _, _ = solve_fresh(H)

        @test w ≈ Float32[2, 3] atol=1e-3
    end

    # -------------------------------------------------------------------
    @testset "1.4 Non-negative clamping" begin
        # H = [1 0; 0 1; 1 1], x = [10, 1, 3]
        # Unconstrained: w = [22/3, -5/3] ≈ [7.33, -1.67]
        # NNLS clamps w₂=0: H₁=[1;0;1], H₁ᵀx=13, H₁ᵀH₁=2, w₁=6.5
        entries = [
            (1, 1, 1.0, 10.0, true, 0),  # col 1, row 1
            (2, 2, 1.0,  1.0, true, 0),  # col 2, row 2
            (1, 3, 1.0,  3.0, true, 0),  # col 1, row 3
            (2, 3, 1.0,  3.0, true, 0),  # col 2, row 3 (shared row)
        ]
        H = build_test_sparse(entries)
        w, _, _, _ = solve_fresh(H)

        @test w[1] ≈ 6.5f0 atol=0.1
        @test w[2] ≈ 0.0f0 atol=1e-4
        @test all(w .>= 0)   # verify non-negativity constraint
    end

    # -------------------------------------------------------------------
    @testset "1.5 Convergence flag" begin
        # Well-conditioned identity system converges quickly
        entries = [
            (1, 1, 1.0, 10.0, true, 0),
            (2, 2, 1.0, 20.0, true, 0),
            (3, 3, 1.0, 30.0, true, 0),
        ]
        H = build_test_sparse(entries)
        _, _, _, (converged, iters) = solve_fresh(H; max_iter=100, tol=Float32(1e-6))

        @test converged == true
        @test iters < 100
    end

    # -------------------------------------------------------------------
    @testset "1.6 Non-convergence with max_iter=1" begin
        # Overlapping system needs multiple iterations; 1 is not enough
        entries = [
            (1, 1, 1.0, 2.0, true, 0),
            (1, 2, 1.0, 5.0, true, 0),
            (2, 2, 1.0, 5.0, true, 0),
            (2, 3, 1.0, 3.0, true, 0),
        ]
        H = build_test_sparse(entries)
        _, _, _, (converged, iters) = solve_fresh(H; max_iter=1, tol=Float32(1e-10))

        @test converged == false
        @test iters == 1
    end

    # -------------------------------------------------------------------
    @testset "1.7 Single column" begin
        # H = [3; 4; 5], x = [6, 8, 10] → w = [2.0]
        entries = [
            (1, 1, 3.0, 6.0, true, 0),
            (1, 2, 4.0, 8.0, true, 0),
            (1, 3, 5.0, 10.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w, _, _, _ = solve_fresh(H)

        @test w ≈ Float32[2.0] atol=1e-4
    end

    # -------------------------------------------------------------------
    @testset "1.8 All-zero observed" begin
        # H = I(2), x = [0, 0] → w = [0, 0]
        entries = [
            (1, 1, 1.0, 0.0, true, 0),
            (2, 2, 1.0, 0.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w, _, _, _ = solve_fresh(H)

        @test w ≈ Float32[0, 0] atol=1e-6
    end

    # -------------------------------------------------------------------
    @testset "1.9 Residual consistency" begin
        # After solving, verify r ≈ Hw - x from dense matrix
        entries = [
            (1, 1, 1.0, 2.0, true, 0),
            (1, 2, 1.0, 5.0, true, 0),
            (2, 2, 1.0, 5.0, true, 0),
            (2, 3, 1.0, 3.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w, r, _, _ = solve_fresh(H)

        # Build dense version and compute expected residuals
        H_dense = to_dense_colptr(H)
        x_dense = Float32[2, 5, 3]  # x values for rows 1, 2, 3
        r_expected = H_dense * w - x_dense

        @test r[1:H.m] ≈ r_expected atol=1e-4
    end

    # -------------------------------------------------------------------
    @testset "1.10 Reference comparison (5-column problem)" begin
        # H: 7 rows × 5 columns with known non-negative solution
        # w_true = [3, 2, 4, 1, 5], x = H * w_true
        entries = [
            # col 1: rows [1, 2], nzval [2, 1]
            (1, 1, 2.0,  6.0, true, 0),   # x[1] = 2*3 = 6
            (1, 2, 1.0,  5.0, true, 0),   # x[2] = 1*3 + 1*2 = 5
            # col 2: rows [2, 3], nzval [1, 3]
            (2, 2, 1.0,  5.0, true, 0),   # x[2] = 5 (shared)
            (2, 3, 3.0, 14.0, true, 0),   # x[3] = 3*2 + 2*4 = 14
            # col 3: rows [3, 4], nzval [2, 2]
            (3, 3, 2.0, 14.0, true, 0),   # x[3] = 14 (shared)
            (3, 4, 2.0,  9.0, true, 0),   # x[4] = 2*4 + 1*1 = 9
            # col 4: rows [4, 5, 6], nzval [1, 3, 1]
            (4, 4, 1.0,  9.0, true, 0),   # x[4] = 9 (shared)
            (4, 5, 3.0,  3.0, true, 0),   # x[5] = 3*1 = 3
            (4, 6, 1.0, 11.0, true, 0),   # x[6] = 1*1 + 2*5 = 11
            # col 5: rows [6, 7], nzval [2, 4]
            (5, 6, 2.0, 11.0, true, 0),   # x[6] = 11 (shared)
            (5, 7, 4.0, 20.0, true, 0),   # x[7] = 4*5 = 20
        ]
        H = build_test_sparse(entries)
        w, _, _, (converged, _) = solve_fresh(H; max_iter=500)

        w_true = Float32[3, 2, 4, 1, 5]
        @test w ≈ w_true atol=0.05
        @test converged == true

        # Sanity check: compare against dense least-squares
        H_dense = Float64.(to_dense_colptr(H))
        x_dense = Float64[6, 5, 14, 9, 3, 11, 20]
        w_ref = H_dense \ x_dense
        @test Float32.(w_ref) ≈ w_true atol=0.01
    end

end  # solveOLS! testset

println("\nAll solveOLS! tests passed!")
