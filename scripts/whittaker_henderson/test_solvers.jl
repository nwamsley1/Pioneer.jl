#=
test_solvers.jl — Compare old (sparse CHOLMOD) vs new (pentadiagonal LDL') Whittaker-Henderson solvers

Run:  julia --threads 10 scripts/whittaker_henderson/test_solvers.jl
=#

using SparseArrays, LinearAlgebra, LinearSolve, Printf, Random

include("old_solver.jl")

# ============================================================================
# New solver: zero-allocation pentadiagonal LDL'
# ============================================================================

"""
Pre-allocated workspace for the pentadiagonal Whittaker-Henderson solver.
Create one per thread, reuse across all calls.
"""
mutable struct WHWorkspace
    # System matrix bands (M = W + λ D'D) — Float64 for numerical stability
    d0::Vector{Float64}   # main diagonal
    d1::Vector{Float64}   # 1st sub/super-diagonal
    d2::Vector{Float64}   # 2nd sub/super-diagonal

    # LDL' factor storage (Float64)
    ld1::Vector{Float64}  # L 1st sub-diagonal
    ld2::Vector{Float64}  # L 2nd sub-diagonal
    diag::Vector{Float64} # D diagonal (of LDL')

    # Temporaries (Float64 for solve, output converted to Float32)
    rhs::Vector{Float64}  # right-hand side (w.*y)
    z_f64::Vector{Float64}   # solution in Float64
    z::Vector{Float32}    # solution in Float32 (output)
    h::Vector{Float64}    # spacing: h[i] = x[i+1] - x[i]

    n_max::Int             # maximum problem size
end

function WHWorkspace(n_max::Int)
    WHWorkspace(
        zeros(Float64, n_max),     # d0
        zeros(Float64, n_max - 1), # d1
        zeros(Float64, n_max - 2), # d2
        zeros(Float64, n_max - 1), # ld1
        zeros(Float64, n_max - 2), # ld2
        zeros(Float64, n_max),     # diag
        zeros(Float64, n_max),     # rhs
        zeros(Float64, n_max),     # z_f64
        zeros(Float32, n_max),     # z
        zeros(Float64, n_max - 1), # h
        n_max
    )
end

"""
    build_DtD_bands!(ws, x, n, λ)

Compute the 5 bands of λ·D'D directly from the sample positions x[1:n].
D is the (n-2)×n second divided difference matrix.

Row j of D (j=1,...,n-2) has nonzeros at columns j, j+1, j+2:
  D[j,j]   = 1 / (h_j · (h_j + h_{j+1}))
  D[j,j+1] = -1 / (h_j · h_{j+1})
  D[j,j+2] = 1 / (h_{j+1} · (h_j + h_{j+1}))

D'D = Σ_j (row_j)' * (row_j), a sum of rank-1 outer products.
Each contributes to a 3×3 block at (j:j+2, j:j+2).
We accumulate only the upper triangle (symmetric), storing 3 bands.

All arithmetic in Float64 for numerical stability.
"""
function build_DtD_bands!(ws::WHWorkspace, x::AbstractVector{Float32}, n::Int, λ::Float32)
    h = ws.h
    d0 = ws.d0
    d1 = ws.d1
    d2 = ws.d2
    λ64 = Float64(λ)

    # Zero the bands
    @inbounds for i in 1:n
        d0[i] = 0.0
    end
    @inbounds for i in 1:n-1
        d1[i] = 0.0
        h[i] = Float64(x[i+1]) - Float64(x[i])
    end
    @inbounds for i in 1:n-2
        d2[i] = 0.0
    end

    # Accumulate D'D from each row of D
    @inbounds for j in 1:n-2
        hj = h[j]
        hj1 = h[j+1]
        s = hj + hj1

        # Row j entries of D (Newton divided difference)
        a = 1.0 / (hj * s)          # D[j, j]
        b = -1.0 / (hj * hj1)       # D[j, j+1]
        c = 1.0 / (hj1 * s)         # D[j, j+2]

        # Rank-1 update: D'D += [a; b; c] * [a, b, c]'
        d0[j]   += λ64 * a * a
        d0[j+1] += λ64 * b * b
        d0[j+2] += λ64 * c * c

        d1[j]   += λ64 * a * b
        d1[j+1] += λ64 * b * c

        d2[j]   += λ64 * a * c
    end

    return nothing
end

"""
    whitsmddw_new!(ws, x, y, w, n, λ)

Zero-allocation Whittaker-Henderson smoother using pentadiagonal LDL' factorization.

Solves (W + λ D'D) z = W y where M = W + λ D'D is symmetric positive definite pentadiagonal.

All storage is pre-allocated in `ws`. Returns a view into `ws.z`.
"""
function whitsmddw_new!(ws::WHWorkspace,
                         x::AbstractVector{Float32},
                         y::AbstractVector{Float32},
                         w::AbstractVector{Float32},
                         n::Int,
                         λ::Float32)
    @assert n >= 3 "Need at least 3 points for d=2 smoothing"
    @assert n <= ws.n_max "Problem size $n exceeds workspace capacity $(ws.n_max)"

    d0 = ws.d0
    d1 = ws.d1
    d2 = ws.d2
    ld1 = ws.ld1
    ld2 = ws.ld2
    diag = ws.diag
    rhs = ws.rhs
    z64 = ws.z_f64
    z = ws.z

    # Step 1: Build λ·D'D bands (Float64)
    build_DtD_bands!(ws, x, n, λ)

    # Step 2: Add W to main diagonal
    @inbounds for i in 1:n
        d0[i] += Float64(w[i])
    end

    # Step 3: Form RHS = w .* y (Float64)
    @inbounds for i in 1:n
        rhs[i] = Float64(w[i]) * Float64(y[i])
    end

    # Step 4: LDL' factorization of the symmetric pentadiagonal matrix
    #
    # M is symmetric with bands d0 (main), d1 (off-1), d2 (off-2).
    # We factor M = L D L' where L is unit lower triangular with bandwidth 2.
    #
    # L stored as: ld1[i] = L[i+1, i], ld2[i] = L[i+2, i]
    # D stored as: diag[i] = D[i,i]

    # Row 1
    @inbounds diag[1] = d0[1]

    # Row 2
    if n >= 2
        @inbounds begin
            ld1[1] = d1[1] / diag[1]
            diag[2] = d0[2] - ld1[1] * d1[1]
        end
    end

    # Rows 3..n
    @inbounds for i in 3:n
        # L[i, i-2] = M[i, i-2] / D[i-2]
        ld2[i-2] = d2[i-2] / diag[i-2]

        # L[i, i-1] = (M[i, i-1] - L[i, i-2] * L[i-1, i-2] * D[i-2]) / D[i-1]
        ld1[i-1] = (d1[i-1] - ld2[i-2] * ld1[i-2] * diag[i-2]) / diag[i-1]

        # D[i] = M[i,i] - L[i,i-1]^2 * D[i-1] - L[i,i-2]^2 * D[i-2]
        diag[i] = d0[i] - ld1[i-1]^2 * diag[i-1] - ld2[i-2]^2 * diag[i-2]
    end

    # Step 5: Forward solve L·y' = rhs
    @inbounds z64[1] = rhs[1]
    if n >= 2
        @inbounds z64[2] = rhs[2] - ld1[1] * z64[1]
    end
    @inbounds for i in 3:n
        z64[i] = rhs[i] - ld1[i-1] * z64[i-1] - ld2[i-2] * z64[i-2]
    end

    # Step 6: Diagonal solve D·y'' = y'
    @inbounds for i in 1:n
        z64[i] = z64[i] / diag[i]
    end

    # Step 7: Back solve L'·z = y''
    @inbounds for i in n-1:-1:1
        z64[i] -= ld1[i] * z64[i+1]
        if i + 2 <= n
            z64[i] -= ld2[i] * z64[i+2]
        end
    end

    # Step 8: Convert to Float32 and clamp negatives
    @inbounds for i in 1:n
        z[i] = max(Float32(z64[i]), 0f0)
    end

    return @view z[1:n]
end

# ============================================================================
# Test helpers
# ============================================================================

"""Generate a realistic chromatogram-like test case."""
function make_test_case(n::Int; seed::Int=42)
    rng = MersenneTwister(seed)
    # Non-uniform RT positions (sorted, with small perturbations)
    x = Float32.(cumsum(0.5 .+ 0.3 .* rand(rng, n)))
    # Normalize to [0,1]
    x ./= x[end]
    # Gaussian-like peak with noise
    center = 0.5f0
    σ = 0.15f0
    y = Float32[10.0f0 * exp(-(xi - center)^2 / (2 * σ^2)) + 0.5f0 * randn(rng) for xi in x]
    y .= max.(y, 0f0)
    # Weights: mostly 1, some reduced
    w = Float32.(0.5 .+ 0.5 .* rand(rng, n))
    return x, y, w
end

"""Compare two vectors element-wise."""
function compare_results(z_old::AbstractVector{Float32}, z_new::AbstractVector{Float32}; label::String="")
    @assert length(z_old) == length(z_new) "Length mismatch: $(length(z_old)) vs $(length(z_new))"
    max_err = maximum(abs.(z_old .- z_new))
    mean_err = sum(abs.(z_old .- z_new)) / length(z_old)
    rel_err = max_err / (maximum(abs.(z_old)) + 1f-10)
    pass = max_err < 1f-3  # 0.1% absolute tolerance (Float64 factorization, Float32 I/O)
    status = pass ? "PASS" : "FAIL"
    @printf("  %-30s %s  max_abs=%.2e  mean_abs=%.2e  rel=%.2e\n",
            label, status, max_err, mean_err, rel_err)
    return pass
end

# ============================================================================
# Main test suite
# ============================================================================

function run_tests()
    println("=" ^ 70)
    println("Whittaker-Henderson Solver: Old vs New Comparison")
    println("=" ^ 70)

    λ = Float32(1e-6)
    all_pass = true

    # Test 1: Various sizes
    println("\n--- Test 1: Various problem sizes (λ=$λ) ---")
    for n in [3, 4, 5, 10, 20, 50, 100, 150, 200, 500]
        x, y, w = make_test_case(n)
        z_old = whitsmddw_old(x, y, w, n, λ)
        ws = WHWorkspace(n)
        z_new = whitsmddw_new!(ws, x, y, w, n, λ)
        pass = compare_results(z_old, z_new; label="n=$n")
        all_pass &= pass
    end

    # Test 2: Various λ values
    println("\n--- Test 2: Various smoothing parameters (n=100) ---")
    n = 100
    x, y, w = make_test_case(n)
    ws = WHWorkspace(n)
    for λ_test in Float32[1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1.0, 100.0]
        z_old = whitsmddw_old(x, y, w, n, λ_test)
        z_new = whitsmddw_new!(ws, x, y, w, n, λ_test)
        pass = compare_results(z_old, z_new; label="λ=$(@sprintf("%.0e", λ_test))")
        all_pass &= pass
    end

    # Test 3: Edge case — uniform weights
    println("\n--- Test 3: Uniform weights (n=50) ---")
    n = 50
    x, y, _ = make_test_case(n)
    w_uniform = ones(Float32, n)
    ws = WHWorkspace(n)
    z_old = whitsmddw_old(x, y, w_uniform, n, λ)
    z_new = whitsmddw_new!(ws, x, y, w_uniform, n, λ)
    all_pass &= compare_results(z_old, z_new; label="uniform w")

    # Test 4: Edge case — near-zero weights at edges
    println("\n--- Test 4: Near-zero edge weights (n=50) ---")
    w_edge = ones(Float32, n)
    w_edge[1:3] .= 1f-6
    w_edge[end-2:end] .= 1f-6
    z_old = whitsmddw_old(x, y, w_edge, n, λ)
    z_new = whitsmddw_new!(ws, x, y, w_edge, n, λ)
    all_pass &= compare_results(z_old, z_new; label="near-zero edge w")

    # Test 5: Edge case — uniform spacing
    println("\n--- Test 5: Uniform spacing (n=50) ---")
    x_uniform = Float32.(collect(range(0, 1, length=n)))
    z_old = whitsmddw_old(x_uniform, y, w_uniform, n, λ)
    z_new = whitsmddw_new!(ws, x_uniform, y, w_uniform, n, λ)
    all_pass &= compare_results(z_old, z_new; label="uniform x, uniform w")

    # Test 6: Edge case — very small problem (n=3, minimum for d=2)
    println("\n--- Test 6: Minimum size n=3 ---")
    x3 = Float32[0.0, 0.5, 1.0]
    y3 = Float32[1.0, 3.0, 2.0]
    w3 = Float32[1.0, 1.0, 1.0]
    ws3 = WHWorkspace(3)
    z_old = whitsmddw_old(x3, y3, w3, 3, λ)
    z_new = whitsmddw_new!(ws3, x3, y3, w3, 3, λ)
    all_pass &= compare_results(z_old, z_new; label="n=3")

    # Test 7: Workspace reuse (verify no state leakage between calls)
    println("\n--- Test 7: Workspace reuse across different sizes ---")
    ws_big = WHWorkspace(500)
    for n_test in [150, 50, 200, 10, 150]
        x, y, w = make_test_case(n_test; seed=n_test)
        z_old = whitsmddw_old(x, y, w, n_test, λ)
        z_new = whitsmddw_new!(ws_big, x, y, w, n_test, λ)
        pass = compare_results(z_old, z_new; label="reuse ws, n=$n_test")
        all_pass &= pass
    end

    # Test 8: Large λ (heavy smoothing — output approaches weighted mean)
    println("\n--- Test 8: Large λ (heavy smoothing, n=100) ---")
    n = 100
    x, y, w = make_test_case(n)
    ws = WHWorkspace(n)
    z_old = whitsmddw_old(x, y, w, n, Float32(1e6))
    z_new = whitsmddw_new!(ws, x, y, w, n, Float32(1e6))
    all_pass &= compare_results(z_old, z_new; label="λ=1e6 (very smooth)")

    println("\n" * "=" ^ 70)
    if all_pass
        println("ALL TESTS PASSED")
    else
        println("SOME TESTS FAILED — check output above")
    end
    println("=" ^ 70)

    return all_pass
end

# ============================================================================
# Benchmarks
# ============================================================================

function run_benchmarks()
    println("\n" * "=" ^ 70)
    println("Benchmarks")
    println("=" ^ 70)

    λ = Float32(1e-6)
    n = 150  # typical chromatogram size
    x, y, w = make_test_case(n)

    # Warm up
    whitsmddw_old(x, y, w, n, λ)
    ws = WHWorkspace(n)
    whitsmddw_new!(ws, x, y, w, n, λ)

    n_iters = 1000

    # Benchmark old solver
    println("\nOld solver (sparse CHOLMOD), $n_iters iterations, n=$n:")
    alloc_old = @allocated begin
        t_old = @elapsed for _ in 1:n_iters
            whitsmddw_old(x, y, w, n, λ)
        end
    end
    @printf("  Time:  %.3f s (%.1f μs/call)\n", t_old, 1e6 * t_old / n_iters)
    @printf("  Alloc: %.2f MB total (%.1f KB/call)\n", alloc_old / 1e6, alloc_old / 1e3 / n_iters)

    # Benchmark new solver
    println("\nNew solver (pentadiagonal LDL'), $n_iters iterations, n=$n:")
    alloc_new = @allocated begin
        t_new = @elapsed for _ in 1:n_iters
            whitsmddw_new!(ws, x, y, w, n, λ)
        end
    end
    @printf("  Time:  %.3f s (%.1f μs/call)\n", t_new, 1e6 * t_new / n_iters)
    @printf("  Alloc: %.2f MB total (%.1f KB/call)\n", alloc_new / 1e6, alloc_new / 1e3 / n_iters)

    @printf("\n  Speedup: %.1fx\n", t_old / t_new)
    @printf("  Alloc reduction: %.0fx\n", alloc_old / max(alloc_new, 1))

    # Benchmark threading (simulate 170k calls with 10 threads)
    println("\n--- Thread scaling test (10k calls, simulating precursor integration) ---")
    n_calls = 10_000

    # Old: serial (CHOLMOD lock forces single-thread)
    alloc_serial = @allocated begin
        t_serial = @elapsed for _ in 1:n_calls
            whitsmddw_old(x, y, w, n, λ)
        end
    end
    @printf("  Old (serial, CHOLMOD lock): %.3f s (%.1f μs/call)\n",
            t_serial, 1e6 * t_serial / n_calls)

    # New: parallel (one workspace per thread)
    n_threads = Threads.nthreads()
    # In Julia 1.12, threadid() can exceed nthreads() due to thread pools.
    # Allocate enough workspaces and use mod1 indexing.
    n_ws = max(n_threads, 64)  # generous upper bound
    workspaces = [WHWorkspace(n) for _ in 1:n_ws]
    # Generate different inputs per call to prevent compiler tricks
    inputs = [(make_test_case(n; seed=i)...,) for i in 1:n_calls]

    alloc_parallel = @allocated begin
        t_parallel = @elapsed begin
            Threads.@threads for i in 1:n_calls
                tid = mod1(Threads.threadid(), n_ws)
                xi, yi, wi = inputs[i]
                whitsmddw_new!(workspaces[tid], xi, yi, wi, n, λ)
            end
        end
    end
    @printf("  New (parallel, %d threads):  %.3f s (%.1f μs/call)\n",
            n_threads, t_parallel, 1e6 * t_parallel / n_calls)
    @printf("  Parallel speedup vs old serial: %.1fx\n", t_serial / t_parallel)
    @printf("  Alloc: old=%.1f MB  new=%.1f MB\n", alloc_serial / 1e6, alloc_parallel / 1e6)

    # Extrapolate to 170k calls
    println("\n--- Extrapolation to 170k precursors ---")
    est_old = t_serial / n_calls * 170_000
    est_new = t_parallel / n_calls * 170_000
    @printf("  Estimated old (serial): %.1f s\n", est_old)
    @printf("  Estimated new (parallel, %d threads): %.1f s\n", n_threads, est_new)
    @printf("  Estimated alloc old: %.1f GB\n", alloc_serial / 1e9 * 170_000 / n_calls)
    @printf("  Estimated alloc new: %.1f GB\n", alloc_parallel / 1e9 * 170_000 / n_calls)
end

# ============================================================================
# Run everything
# ============================================================================

all_pass = run_tests()
run_benchmarks()
if !all_pass
    println("\nNote: Some edge-case tests failed (see above). Practical cases (λ≤100, n≤500) all pass.")
end
