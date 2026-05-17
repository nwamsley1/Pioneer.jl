# Unit tests for src/utils/SpectralDeconvolution/solvePoissonMM.jl.
#
# The solver maximizes the Poisson log-likelihood   max Σ y_i log μ_i - μ_i
# under the model μ = A x and the constraint x_j ≥ 0. Most numerical
# packages (GLM.jl, NMF.jl) don't ship a non-negativity-constrained
# Poisson regression with an OLS-style API, so we build ground truth
# directly from KKT conditions and from closed-form solutions on
# diagonal A.
#
# At the optimum (with the loss-gradient sign convention used by
# `solvePoissonMM_fast!`,  L1_j = Σ_i a_ij (1 - y_i / μ_i) ):
#   x_j > 0  →  L1_j ≈ 0          (interior: gradient zero)
#   x_j = 0  →  L1_j ≥ 0          (boundary: increasing x_j wouldn't
#                                  decrease loss)

using Pioneer: SparseArrayFused, solvePoissonMM_fast!, initMu!,
               initObserved!, updateMu!, solve_deconvolution!,
               OLSSolver, PoissonMMSolver

# ────────────────────────────────────────────────────────────────────
# Helper: build a SparseArrayFused from row-indexed (col, row, val) triples.
# Triples must be grouped by column ascending. m is the number of rows.
function _make_sparse_fused(triples::Vector{<:Tuple{<:Integer,<:Integer,<:Real}},
                            n_rows::Integer, n_cols::Integer)
    n_vals = length(triples)
    sa = SparseArrayFused(Int64(max(n_vals, 1)))
    sa.m = Int64(n_rows)
    sa.n = Int64(n_cols)
    sa.n_vals = Int64(n_vals)
    if length(sa.colptr) < n_cols + 1
        append!(sa.colptr, zeros(Int64, n_cols + 1 - length(sa.colptr)))
    end
    @inbounds for c in 1:(n_cols + 1)
        sa.colptr[c] = Int64(0)
    end
    @inbounds for c in 1:n_cols
        sa.colptr[c] = Int64(0)
    end
    # Fill colptr by counting per-column entries.
    counts = zeros(Int64, n_cols)
    for (col, _row, _val) in triples
        counts[col] += 1
    end
    sa.colptr[1] = Int64(1)
    for c in 1:n_cols
        sa.colptr[c + 1] = sa.colptr[c] + counts[c]
    end
    # Place entries column-by-column following triples order
    cursor = zeros(Int64, n_cols)
    for c in 1:n_cols
        cursor[c] = sa.colptr[c]
    end
    @inbounds for (col, row, val) in triples
        idx = cursor[col]
        sa.rowval[idx] = Int64(row)
        sa.nzval[idx]  = Float32(val)
        cursor[col] += 1
    end
    return sa
end

# Loss gradient L1_j = Σ_i a_ij (1 - y_i / μ_i)
function _loss_gradient(sa, μ::Vector{Float32}, y::Vector{Float32}, j::Int)
    g = 0.0
    for i in sa.colptr[j]:(sa.colptr[j + 1] - 1)
        row = sa.rowval[i]
        a   = Float64(sa.nzval[i])
        μ_i = max(Float64(μ[row]), 1e-6)
        y_i = Float64(y[row])
        g += a * (1.0 - y_i / μ_i)
    end
    return g
end

@testset "solvePoissonMM" begin

    @testset "initMu!: μ = A * w with floor clamp" begin
        # 2 cols × 3 rows, dense block:
        #   col 1: row 1=2.0, row 2=1.0
        #   col 2: row 2=3.0, row 3=4.0
        sa = _make_sparse_fused([(1,1,2.0), (1,2,1.0), (2,2,3.0), (2,3,4.0)], 3, 2)
        w = Float32[2.0, 0.5]
        μ = zeros(Float32, 4)
        initMu!(μ, sa, w)
        # μ = [2*2, 1*2 + 3*0.5, 4*0.5] = [4.0, 3.5, 2.0]
        @test μ[1] ≈ 4.0f0
        @test μ[2] ≈ 3.5f0
        @test μ[3] ≈ 2.0f0
        # Floor clamp: μ unaffected here (all > 1e-6)
        @test all(μ[1:3] .≥ 1f-6)
    end

    @testset "initMu!: zero-weight rows clamped to POISSON_MU_FLOOR" begin
        sa = _make_sparse_fused([(1,1,2.0)], 2, 1)
        w = Float32[1.0]
        μ = zeros(Float32, 3)
        initMu!(μ, sa, w)
        @test μ[1] ≈ 2.0f0
        # Row 2 has no entries → unclamped μ would be 0; floor pushes to 1e-6
        @test μ[2] == 1f-6
    end

    @testset "initMu! grows μ if too short" begin
        sa = _make_sparse_fused([(1,1,2.0), (1,5,3.0)], 5, 1)
        w = Float32[1.0]
        μ = Float32[]   # empty; initMu! must grow it to length 5
        initMu!(μ, sa, w)
        @test length(μ) ≥ 5
        @test μ[1] ≈ 2.0f0
        @test μ[5] ≈ 3.0f0
    end

    @testset "initObserved!: extracts sa.x into y" begin
        sa = _make_sparse_fused([(1,1,2.0), (1,2,1.0), (2,3,4.0)], 3, 2)
        sa.x[1] = 7.0f0   # row 1 observation
        sa.x[2] = 5.0f0   # row 2 observation
        sa.x[3] = 9.0f0   # row 3 observation
        y = zeros(Float32, 3)
        initObserved!(y, sa)
        @test y[1] ≈ 7.0f0
        @test y[2] ≈ 5.0f0
        @test y[3] ≈ 9.0f0
    end

    @testset "updateMu!: μ updated by Δ × column entries" begin
        sa = _make_sparse_fused([(1,1,2.0), (1,2,1.0)], 3, 1)
        w = Float32[1.0]
        μ = zeros(Float32, 3)
        initMu!(μ, sa, w)
        μ_before = copy(μ)
        # Bump column 1 weight from 1.0 → 3.0 (Δ = +2.0)
        updateMu!(sa, μ, 1, 3.0f0, 1.0f0)
        @test μ[1] ≈ μ_before[1] + 2.0f0 * 2.0f0   # row 1 a=2 → +4
        @test μ[2] ≈ μ_before[2] + 2.0f0 * 1.0f0   # row 2 a=1 → +2
    end

    @testset "single-column closed form: x = Σy / Σa" begin
        # 4 rows, 1 col, all matching: a = [1, 1, 2, 4], y = [3, 1, 6, 4]
        # Poisson MLE for x: x = (Σ y_i) / (Σ a_i) = 14/8 = 1.75
        sa = _make_sparse_fused(
            [(1,1,1.0), (1,2,1.0), (1,3,2.0), (1,4,4.0)], 4, 1)
        y = Float32[3.0, 1.0, 6.0, 4.0]
        X = Float32[0.5]   # warm start
        μ = zeros(Float32, 4)
        initMu!(μ, sa, X)
        converged, iters = solvePoissonMM_fast!(sa, μ, y, X, Int64(50), 1f-5)
        @test converged
        @test X[1] ≈ 1.75f0 atol=5f-3
        @test X[1] ≥ 0
    end

    @testset "diagonal A: independent columns, closed form per column" begin
        # 3 cols, 6 rows. Each col has 2 rows. Expect:
        #   x_j = (y_{2j-1} + y_{2j}) / (a + b) for col j with values a, b
        sa = _make_sparse_fused([
            (1,1,2.0), (1,2,3.0),
            (2,3,1.0), (2,4,1.0),
            (3,5,4.0), (3,6,2.0),
        ], 6, 3)
        y = Float32[10.0, 5.0,   2.0, 8.0,   4.0, 4.0]
        X = ones(Float32, 3)
        μ = zeros(Float32, 6)
        initMu!(μ, sa, X)
        converged, _ = solvePoissonMM_fast!(sa, μ, y, X, Int64(50), 1f-5)
        @test converged
        # x_1 = (10 + 5)  / (2 + 3) = 3.0
        # x_2 = (2 + 8)   / (1 + 1) = 5.0
        # x_3 = (4 + 4)   / (4 + 2) = 8/6 ≈ 1.333
        @test X[1] ≈ 3.0f0    atol=5f-3
        @test X[2] ≈ 5.0f0    atol=5f-3
        @test X[3] ≈ 1.3333f0 atol=5f-3
        @test all(X .≥ 0)
    end

    @testset "y = 0 drives x → 0 (no observations to support)" begin
        sa = _make_sparse_fused([(1,1,2.0), (1,2,1.0), (2,2,3.0)], 2, 2)
        y = zeros(Float32, 2)
        X = Float32[1.0, 1.0]
        μ = zeros(Float32, 2)
        initMu!(μ, sa, X)
        # Poisson loss with y=0 is just Σ μ; minimum is at x=0.
        converged, _ = solvePoissonMM_fast!(sa, μ, y, X, Int64(100), 1f-6)
        @test converged
        # Tolerance: with relative threshold 1e-6 the solver may stop above 0
        @test X[1] ≤ 1f-3
        @test X[2] ≤ 1f-3
        @test all(X .≥ 0)
    end

    @testset "non-negativity active: column with no support pinned at 0" begin
        # 3 rows, 2 cols. Col 1 covers rows 1,2; col 2 covers row 3.
        # y has support only on rows 1,2 (y_3 = 0). The MLE wants x_2 → 0
        # since col 2 contributes nothing to the likelihood except cost.
        sa = _make_sparse_fused([
            (1,1,1.0), (1,2,1.0),
            (2,3,1.0),
        ], 3, 2)
        y = Float32[2.0, 4.0, 0.0]
        X = Float32[1.0, 1.0]   # warm start col 2 at 1.0
        μ = zeros(Float32, 3)
        initMu!(μ, sa, X)
        converged, _ = solvePoissonMM_fast!(sa, μ, y, X, Int64(100), 1f-6)
        @test converged
        # x_1 → (2 + 4) / (1 + 1) = 3.0
        @test X[1] ≈ 3.0f0 atol=5f-3
        @test X[2] ≤ 1f-3
        @test X[2] ≥ 0
    end

    @testset "KKT conditions hold at convergence (mixed-coupling problem)" begin
        # 4 rows, 3 cols, every col touches multiple rows. Ground truth not
        # closed-form; verify KKT instead.
        sa = _make_sparse_fused([
            (1,1,1.0), (1,2,2.0), (1,3,1.0),
            (2,1,2.0), (2,2,1.0), (2,4,1.0),
            (3,3,3.0), (3,4,2.0),
        ], 4, 3)
        y = Float32[3.0, 4.0, 6.0, 2.0]
        X = Float32[0.5, 0.5, 0.5]
        μ = zeros(Float32, 4)
        initMu!(μ, sa, X)
        converged, _ = solvePoissonMM_fast!(sa, μ, y, X, Int64(200), 1f-7)
        @test converged

        # Recompute μ from scratch with the converged X, then check KKT.
        μ_check = zeros(Float32, 4)
        initMu!(μ_check, sa, X)
        for j in 1:3
            grad = _loss_gradient(sa, μ_check, y, j)
            if X[j] > 1f-3
                # Interior: gradient ≈ 0
                @test abs(grad) ≤ 5f-2
            else
                # Boundary: gradient ≥ 0 (loss can only go up by increasing x_j)
                @test grad ≥ -5f-2
            end
        end
        @test all(X .≥ 0)
    end

    @testset "y-scaling kicks in for large counts and is restored on exit" begin
        # max(y) = 1000 → y_scale = 1000; y is divided then multiplied back.
        sa = _make_sparse_fused([(1,1,1.0), (1,2,1.0)], 2, 1)
        y = Float32[400.0, 600.0]
        y_orig = copy(y)
        X = Float32[1.0]
        μ = zeros(Float32, 2)
        initMu!(μ, sa, X)
        converged, _ = solvePoissonMM_fast!(sa, μ, y, X, Int64(50), 1f-5)
        @test converged
        # y must be restored to original after exit
        @test y ≈ y_orig
        # Closed-form: x = (400 + 600) / (1 + 1) = 500
        @test X[1] ≈ 500.0f0 atol=2.0f0
    end

    @testset "warm-start at optimum: converges in 1 outer iter" begin
        # Same single-column problem, but start at the closed-form optimum.
        sa = _make_sparse_fused([(1,1,2.0), (1,2,3.0)], 2, 1)
        y = Float32[6.0, 9.0]   # x_opt = 15 / 5 = 3.0
        X = Float32[3.0]
        μ = zeros(Float32, 2)
        initMu!(μ, sa, X)
        converged, iters = solvePoissonMM_fast!(sa, μ, y, X, Int64(20), 1f-5)
        @test converged
        @test iters ≤ 1
        @test X[1] ≈ 3.0f0 atol=1f-3
    end

    @testset "solve_deconvolution! dispatch" begin
        # Single-column problem, observations stored in sa.x (the dispatch path
        # initializes y from sa.x via initObserved!). r and colnorm2 are needed
        # to satisfy the OLS interface but unused by PoissonMMSolver beyond
        # the initial allocation check.
        sa = _make_sparse_fused([(1,1,1.0), (1,2,1.0), (1,3,2.0)], 3, 1)
        sa.x[1] = 2.0f0
        sa.x[2] = 4.0f0
        sa.x[3] = 6.0f0   # observations: y = [2, 4, 6/2] but initObserved
                          # uses sa.x[n] for the row, so for col 1 row 3 with
                          # value 2, sa.x[3] = 6 → y[3] = 6.

        w = Float32[1.0]
        r = zeros(Float32, 3)
        colnorm2 = Float32[1.0]
        μ = zeros(Float32, 3)
        y = zeros(Float32, 3)

        @testset "PoissonMMSolver path" begin
            X = copy(w)
            converged, _ = solve_deconvolution!(
                PoissonMMSolver(), sa, copy(r), X, copy(colnorm2),
                copy(μ), copy(y), Int64(50), 1f-5)
            @test converged
            # x = (Σ y) / (Σ a) = (2 + 4 + 6) / (1 + 1 + 2) = 3.0
            @test X[1] ≈ 3.0f0 atol=5f-3
            @test X[1] ≥ 0
        end

        @testset "OLSSolver path runs without error" begin
            # Verifies dispatch works; we don't assert OLS optimality here
            # because OLS treats the problem differently and isn't the
            # subject under test.
            X = copy(w)
            solve_deconvolution!(
                OLSSolver(), sa, copy(r), X, copy(colnorm2),
                copy(μ), copy(y), Int64(50), 1f-5)
            @test all(isfinite, X)
        end

        @testset "HuberSolver path restores robust deconvolution" begin
            @test isdefined(Pioneer, :HuberSolver)
            if isdefined(Pioneer, :HuberSolver)
                HuberSolver = getproperty(Pioneer, :HuberSolver)
                sa_huber = _make_sparse_fused(
                    [(1,1,1.0), (1,2,1.0), (1,3,1.0)], 3, 1)
                sa_huber.x[1] = 2.0f0
                sa_huber.x[2] = 4.0f0
                sa_huber.x[3] = 6.0f0

                X = Float32[1.0]
                converged, _ = solve_deconvolution!(
                    HuberSolver(1.0f0, 0.0f0, Int64(50), Int64(100),
                                1.0f-5, 1.0f-4, Pioneer.NoNorm()),
                    sa_huber, copy(r), X, copy(colnorm2), copy(μ), copy(y),
                    Int64(100), 1.0f-5)

                @test converged
                @test X[1] ≈ 4.0f0 atol=1.0f-2
                @test X[1] ≥ 0.0f0
            end
        end
    end

end
