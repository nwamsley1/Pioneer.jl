# Self-contained tests for spectralPoissonRegression.jl
# Run from this directory:  julia test_poisson.jl

using Random, Printf

include("SparseArray.jl")
include("spectralPoissonRegression.jl")

# ── Solver defaults ──────────────────────────────────────────────
const NR_ITER   = 25
const BS_ITER   = 100
const OUT_ITER  = 2000
const NR_ACC    = 1f-6
const BS_ACC    = 1f-6
const REL_CONV  = 1f-4

function run_poisson(A::Matrix{Float32}, y::Vector{Float32}; x0::Union{Nothing,Vector{Float32}}=nothing)
    sa = build_SparseArray(A, y)
    m, n = size(A)

    μ = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = x0 === nothing ? ones(Float32, n) : copy(x0)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    solvePoisson!(sa, μ, yy, w,
                  NR_ITER, BS_ITER, OUT_ITER,
                  NR_ACC, BS_ACC, REL_CONV)
    return w, μ, yy, sa
end

function run_poisson_mm(A::Matrix{Float32}, y::Vector{Float32}; x0::Union{Nothing,Vector{Float32}}=nothing)
    sa = build_SparseArray(A, y)
    m, n = size(A)

    μ = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = x0 === nothing ? ones(Float32, n) : copy(x0)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    solvePoissonMM!(sa, μ, yy, w, OUT_ITER, REL_CONV)
    return w, μ, yy, sa
end

# ── Test 1: Single-column recovery ───────────────────────────────
function test_single_column()
    println("\n── Test 1: Single-column recovery ──")
    A = Float32[1; 2; 3;;]          # 3×1
    x_true = 5.0f0
    y = Float32[1*5; 2*5; 3*5]     # exact Poisson mean (no noise for deterministic test)

    w, μ, _, _ = run_poisson(A, y)
    err = abs(w[1] - x_true) / x_true
    @printf("  true = %.2f, recovered = %.4f, rel_err = %.2e\n", x_true, w[1], err)
    passed = err < 1e-3
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 2: Two well-separated columns ──────────────────────────
function test_two_columns()
    println("\n── Test 2: Two well-separated columns ──")
    A = Float32[3 0;
                0 4;
                1 1]
    x_true = Float32[10, 7]
    y = A * x_true   # [30, 28, 17]

    w, μ, _, _ = run_poisson(A, y)
    err1 = abs(w[1] - x_true[1]) / x_true[1]
    err2 = abs(w[2] - x_true[2]) / x_true[2]
    @printf("  true = [%.1f, %.1f], recovered = [%.4f, %.4f]\n",
            x_true[1], x_true[2], w[1], w[2])
    @printf("  rel_err = [%.2e, %.2e]\n", err1, err2)
    passed = max(err1, err2) < 1e-3
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 3: Noisy Poisson data ───────────────────────────────────
function test_noisy_poisson()
    println("\n── Test 3: Noisy Poisson (random A, y ~ Poisson(Ax)) ──")
    Random.seed!(42)
    m, n = 50, 5
    A = Float32.(max.(randn(m, n) .* 0.3 .+ 1.0, 0.0))  # positive-ish design
    x_true = Float32[3, 7, 2, 10, 5]
    mu_true = A * x_true

    # Sample Poisson counts
    y = Float32.([rand_poisson(max(mu_true[i], 0.01f0)) for i in 1:m])

    w, μ, _, _ = run_poisson(A, y; x0=ones(Float32, n))
    @printf("  true x    = %s\n", string(round.(x_true; digits=2)))
    @printf("  recovered = %s\n", string(round.(w[1:n]; digits=2)))
    rel_errs = abs.(w[1:n] .- x_true) ./ max.(x_true, 1f-3)
    max_err = maximum(rel_errs)
    @printf("  max rel err = %.3f\n", max_err)
    passed = max_err < 0.5  # generous tolerance for noisy Poisson data
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 4: Zero counts ─────────────────────────────────────────
function test_zero_counts()
    println("\n── Test 4: Zero counts (no NaN/Inf) ──")
    A = Float32[2 1;
                1 3;
                1 0;
                0 1]
    x_true = Float32[5, 3]
    y = Float32[0, 0, 0, 3]   # lots of zeros

    w, μ, _, sa = run_poisson(A, y)
    any_nan = any(isnan, w[1:2]) || any(isnan, μ[1:sa.m])
    any_inf = any(isinf, w[1:2]) || any(isinf, μ[1:sa.m])
    @printf("  recovered = [%.4f, %.4f]\n", w[1], w[2])
    @printf("  any NaN = %s, any Inf = %s\n", any_nan, any_inf)
    passed = !any_nan && !any_inf
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 5: Log-likelihood monotonically increases ───────────────
function test_loglik_monotonic()
    println("\n── Test 5: Log-likelihood monotonicity ──")
    A = Float32[3 1;
                1 4;
                2 2]
    x_true = Float32[8, 6]
    y = A * x_true

    sa = build_SparseArray(A, y)
    m, n = size(A)
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    # Run one outer iteration at a time and track LL
    lls = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, sa.m))

    for _ in 1:50
        for col in 1:sa.n
            newton_bisection_poisson!(sa, μ, yy, w, col,
                                      NR_ITER, BS_ITER, NR_ACC, BS_ACC, REL_CONV)
        end
        push!(lls, poissonLogLikelihood(μ, yy, sa.m))
    end

    monotonic = all(lls[i+1] >= lls[i] - 1e-6 for i in 1:length(lls)-1)
    @printf("  LL start = %.4f, LL end = %.4f\n", lls[1], lls[end])
    @printf("  monotonic = %s (checked %d iterations)\n", monotonic, length(lls)-1)
    println("  ", monotonic ? "PASS" : "FAIL")
    return monotonic
end

# ── Test 6: MM solver — single-column recovery ───────────────────
function test_mm_single_column()
    println("\n── Test 6: MM solver — single-column recovery ──")
    A = Float32[1; 2; 3;;]
    x_true = 5.0f0
    y = Float32[1*5; 2*5; 3*5]

    w, μ, _, _ = run_poisson_mm(A, y)
    err = abs(w[1] - x_true) / x_true
    @printf("  true = %.2f, recovered = %.4f, rel_err = %.2e\n", x_true, w[1], err)
    passed = err < 1e-3
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 7: MM solver — two columns ──────────────────────────────
function test_mm_two_columns()
    println("\n── Test 7: MM solver — two well-separated columns ──")
    A = Float32[3 0;
                0 4;
                1 1]
    x_true = Float32[10, 7]
    y = A * x_true

    w, μ, _, _ = run_poisson_mm(A, y)
    err1 = abs(w[1] - x_true[1]) / x_true[1]
    err2 = abs(w[2] - x_true[2]) / x_true[2]
    @printf("  true = [%.1f, %.1f], recovered = [%.4f, %.4f]\n",
            x_true[1], x_true[2], w[1], w[2])
    @printf("  rel_err = [%.2e, %.2e]\n", err1, err2)
    passed = max(err1, err2) < 1e-3
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 8: MM solver — noisy Poisson ────────────────────────────
function test_mm_noisy()
    println("\n── Test 8: MM solver — noisy Poisson ──")
    Random.seed!(42)
    m, n = 50, 5
    A = Float32.(max.(randn(m, n) .* 0.3 .+ 1.0, 0.0))
    x_true = Float32[3, 7, 2, 10, 5]
    mu_true = A * x_true
    y = Float32.([rand_poisson(max(mu_true[i], 0.01f0)) for i in 1:m])

    w, μ, _, _ = run_poisson_mm(A, y; x0=ones(Float32, n))
    rel_errs = abs.(w[1:n] .- x_true) ./ max.(x_true, 1f-3)
    max_err = maximum(rel_errs)
    @printf("  true x    = %s\n", string(round.(x_true; digits=2)))
    @printf("  recovered = %s\n", string(round.(w[1:n]; digits=2)))
    @printf("  max rel err = %.3f\n", max_err)
    passed = max_err < 0.5
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 9: MM solver — zero counts (no NaN/Inf) ────────────────
function test_mm_zero_counts()
    println("\n── Test 9: MM solver — zero counts ──")
    A = Float32[2 1;
                1 3;
                1 0;
                0 1]
    y = Float32[0, 0, 0, 3]

    w, μ, _, sa = run_poisson_mm(A, y)
    any_nan = any(isnan, w[1:2]) || any(isnan, μ[1:sa.m])
    any_inf = any(isinf, w[1:2]) || any(isinf, μ[1:sa.m])
    @printf("  recovered = [%.4f, %.4f]\n", w[1], w[2])
    @printf("  any NaN = %s, any Inf = %s\n", any_nan, any_inf)
    passed = !any_nan && !any_inf
    println("  ", passed ? "PASS" : "FAIL")
    return passed
end

# ── Test 10: MM solver — LL monotonicity ─────────────────────────
function test_mm_loglik_monotonic()
    println("\n── Test 10: MM solver — LL monotonicity ──")
    A = Float32[3 1;
                1 4;
                2 2]
    x_true = Float32[8, 6]
    y = A * x_true

    sa = build_SparseArray(A, y)
    m, n = size(A)
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    ε = Float32(POISSON_MU_FLOOR)
    lls = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, sa.m))

    for _ in 1:50
        for col in 1:sa.n
            L1, L2 = getPoissonDerivativesObs!(sa, μ, yy, col)
            X0 = w[col]
            if L2 > ε && !isnan(L1)
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
            end
        end
        push!(lls, poissonLogLikelihood(μ, yy, sa.m))
    end

    monotonic = all(lls[i+1] >= lls[i] - 1e-6 for i in 1:length(lls)-1)
    @printf("  LL start = %.4f, LL end = %.4f\n", lls[1], lls[end])
    @printf("  monotonic = %s (checked %d iterations)\n", monotonic, length(lls)-1)
    println("  ", monotonic ? "PASS" : "FAIL")
    return monotonic
end

# ── Test 11: Compare against GLM.jl (if available) ───────────────
# GLM test is in a separate file to avoid world-age issues.
# Run test_glm_comparison.jl separately if GLM.jl is installed.
function test_glm_comparison()
    println("\n── Test 11: GLM.jl comparison (identity link) ──")
    println("  SKIP (run test_glm_comparison.jl separately)")
    return true
end

# ── Simple Poisson sampler (no Distributions.jl dependency) ──────
function rand_poisson(λ::Float64)
    λ <= 0 && return 0
    L = exp(-λ)
    k = 0
    p = 1.0
    while true
        k += 1
        p *= rand()
        p < L && return k - 1
    end
end
rand_poisson(λ::Float32) = rand_poisson(Float64(λ))

# ── Run all tests ────────────────────────────────────────────────
function main()
    println("=" ^ 60)
    println("  Poisson MLE Coordinate Descent — Test Suite")
    println("=" ^ 60)

    results = Bool[]
    push!(results, test_single_column())
    push!(results, test_two_columns())
    push!(results, test_noisy_poisson())
    push!(results, test_zero_counts())
    push!(results, test_loglik_monotonic())
    push!(results, test_mm_single_column())
    push!(results, test_mm_two_columns())
    push!(results, test_mm_noisy())
    push!(results, test_mm_zero_counts())
    push!(results, test_mm_loglik_monotonic())
    push!(results, test_glm_comparison())

    println("\n" * "=" ^ 60)
    n_pass = count(results)
    n_total = length(results)
    println("  Results: $n_pass / $n_total passed")
    println("=" ^ 60)

    return all(results)
end

main()
