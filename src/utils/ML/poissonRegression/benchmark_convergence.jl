# Benchmark: Poisson Fisher vs Poisson Observed-Hessian vs Huber (≈OLS)
#
# Evaluates all solvers under BOTH Poisson LL and SSE so the comparison is fair.
# (RMSE alone favors Huber by construction; Poisson LL alone favors Poisson solvers.)
#
# Run:  julia benchmark_convergence.jl
# Requires: GLM, DataFrames, Distributions

using Random, Printf, Statistics, GLM, DataFrames, Distributions

include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Problem generation ───────────────────────────────────────────

function rand_poisson_manual(λ::Float64)
    if λ <= 0
        return 0
    elseif λ < 30
        L = exp(-λ)
        k = 0
        p = 1.0
        while true
            k += 1
            p *= rand()
            p < L && return k - 1
        end
    else
        while true
            x = round(Int, λ + sqrt(λ) * randn())
            x >= 0 && return x
        end
    end
end

"""
Generate a sparse deconvolution problem mimicking mass-spec data.
"""
function generate_problem(; m=2000, n=500, density=0.02, n_active=80, seed=2024)
    Random.seed!(seed)

    A = zeros(Float32, m, n)
    nnz_per_col = max(1, round(Int, m * density))
    for j in 1:n
        rows = sort(randperm(m)[1:nnz_per_col])
        for i in rows
            A[i, j] = Float32(rand() * 10 + 0.5)
        end
    end

    x_true = zeros(Float32, n)
    active_cols = sort(randperm(n)[1:n_active])
    for j in active_cols
        x_true[j] = Float32(rand() * 20 + 1)
    end

    mu_true = A * x_true
    y = Float32.([rand_poisson_manual(Float64(max(mu_true[i], 0.0f0))) for i in 1:m])

    return A, x_true, y, active_cols
end

# ── Residual init for Huber (no @turbo dependency) ───────────────

function initResiduals_plain!(r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer, T<:AbstractFloat}
    if length(r) < sa.m
        append!(r, zeros(T, sa.m - length(r)))
    end
    for i in 1:sa.m
        r[i] = zero(T)
    end
    for n in 1:sa.n_vals
        if iszero(r[sa.rowval[n]])
            r[sa.rowval[n]] = -sa.x[n]
        end
    end
    for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1] - 1)
            r[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
end

# ── Model-fair evaluation functions ──────────────────────────────

function sse_objective(r::Vector{Float32}, m::Int)
    s = 0.0
    @inbounds for i in 1:m
        s += Float64(r[i])^2
    end
    return s
end

"""
Compute μ = Aw from dense A and weight vector (Float64 output for precision).
"""
function compute_mu_dense(A::Matrix{Float32}, w)
    m, n = size(A)
    μ = zeros(Float64, m)
    for j in 1:min(n, length(w))
        wj = Float64(w[j])
        iszero(wj) && continue
        for i in 1:m
            μ[i] += Float64(A[i,j]) * wj
        end
    end
    return μ
end

"""
Poisson deviance: 2 Σ [y_i log(y_i/μ_i) - (y_i - μ_i)].
Standardized goodness-of-fit; lower is better.
"""
function poisson_deviance(μ::Vector{Float64}, y::Vector{Float32})
    dev = 0.0
    for i in eachindex(y)
        μ_i = max(μ[i], 1e-6)
        y_i = Float64(y[i])
        if y_i > 0
            dev += y_i * log(y_i / μ_i) - (y_i - μ_i)
        else
            dev += μ_i  # lim y→0: y*log(y/μ) → 0, so term = μ
        end
    end
    return 2.0 * dev
end

"""
SSE from predicted μ (Float64) and observed y.
"""
function sse_from_mu(μ::Vector{Float64}, y::Vector{Float32})
    s = 0.0
    for i in eachindex(y)
        s += (Float64(y[i]) - μ[i])^2
    end
    return s
end

"""
Poisson LL from predicted μ (Float64) and observed y.
"""
function poisson_ll_from_mu64(μ::Vector{Float64}, y::Vector{Float32})
    ll = 0.0
    for i in eachindex(y)
        μ_i = max(μ[i], 1e-6)
        y_i = Float64(y[i])
        ll += y_i * log(μ_i) - μ_i
    end
    return ll
end

"""
SSE from Poisson solver's internal Float32 μ and y vectors.
"""
function sse_from_poisson_state(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    s = 0.0
    @inbounds for i in 1:m
        s += (Float64(y[i]) - Float64(μ[i]))^2
    end
    return s
end

"""
Poisson LL from Huber residuals: μ_i = y_i + r_i  (since r = Aw - y).
"""
function poisson_ll_from_residuals(r::Vector{Float32}, y::Vector{Float32}, m::Int)
    ll = 0.0
    @inbounds for i in 1:m
        μ_i = max(Float64(y[i]) + Float64(r[i]), 1e-6)
        y_i = Float64(y[i])
        ll += y_i * log(μ_i) - μ_i
    end
    return ll
end

"""
SSE from Huber residuals (just Σ r_i², same as sse_objective but for consistency).
"""
function sse_from_residuals(r::Vector{Float32}, m::Int)
    return sse_objective(r, m)
end

# ── Instrumented solvers (track both Poisson LL and SSE) ─────────

function solve_poisson_fisher_instrumented(sa, y_vec;
        max_outer=500, nr_iter=25, bs_iter=100,
        nr_acc=1f-6, bs_acc=1f-6, rel_conv=1f-4)

    m, n = sa.m, sa.n
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    lls       = Float64[]
    sses      = Float64[]
    max_diffs = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, m))
    push!(sses, sse_from_poisson_state(μ, yy, m))

    for iter in 1:max_outer
        _diff = 0f0
        for col in 1:n
            δx = abs(newton_bisection_poisson!(sa, μ, yy, w, col,
                        nr_iter, bs_iter, nr_acc, bs_acc, rel_conv))
            if !iszero(w[col])
                rc = δx / abs(w[col])
                rc > _diff && (_diff = rc)
            end
        end
        push!(lls, poissonLogLikelihood(μ, yy, m))
        push!(sses, sse_from_poisson_state(μ, yy, m))
        push!(max_diffs, Float64(_diff))
        _diff < rel_conv && break
    end
    return w, lls, sses, max_diffs
end

function solve_poisson_observed_instrumented(sa, y_vec;
        max_outer=500, nr_iter=25, bs_iter=100,
        nr_acc=1f-6, bs_acc=1f-6, rel_conv=1f-4)

    m, n = sa.m, sa.n
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    lls       = Float64[]
    sses      = Float64[]
    max_diffs = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, m))
    push!(sses, sse_from_poisson_state(μ, yy, m))

    for iter in 1:max_outer
        _diff = 0f0
        for col in 1:n
            δx = abs(newton_bisection_poisson_obs!(sa, μ, yy, w, col,
                        nr_iter, bs_iter, nr_acc, bs_acc, rel_conv))
            if !iszero(w[col])
                rc = δx / abs(w[col])
                rc > _diff && (_diff = rc)
            end
        end
        push!(lls, poissonLogLikelihood(μ, yy, m))
        push!(sses, sse_from_poisson_state(μ, yy, m))
        push!(max_diffs, Float64(_diff))
        _diff < rel_conv && break
    end
    return w, lls, sses, max_diffs
end

function solve_poisson_mm_instrumented(sa, y_vec;
        max_outer=500, rel_conv=1f-4, max_inner=5)

    m, n = sa.m, sa.n
    μ  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initMu!(μ, sa, w)

    lls       = Float64[]
    sses      = Float64[]
    max_diffs = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, m))
    push!(sses, sse_from_poisson_state(μ, yy, m))

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0

    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0

        for col in 1:n
            X_before = w[col]

            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, yy, col)
                if L2 <= ε || isnan(L1)
                    break
                end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                if iszero(w[col]) || (!iszero(X0) && abs_step / abs(X0) < rel_conv)
                    break
                end
            end

            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col])
                rc > _diff && (_diff = rc)
            end
        end
        push!(lls, poissonLogLikelihood(μ, yy, m))
        push!(sses, sse_from_poisson_state(μ, yy, m))
        push!(max_diffs, Float64(_diff))
        _diff < rel_conv && break
    end
    return w, lls, sses, max_diffs
end

function solve_huber_instrumented(sa, y_vec;
        max_outer=500, nr_iter=25, bs_iter=100,
        nr_acc=1f-6, bs_acc=1f-6, rel_conv=1f-4,
        δ=1f10, λ=0f0)

    m, n = sa.m, sa.n
    r  = zeros(Float32, m)
    yy = zeros(Float32, m)
    w  = ones(Float32, n)

    initObserved!(yy, sa)
    initResiduals_plain!(r, sa, w)

    lls       = Float64[]
    sses      = Float64[]
    max_diffs = Float64[]
    push!(lls, poisson_ll_from_residuals(r, yy, m))
    push!(sses, sse_objective(r, m))

    reg = NoNorm()
    for iter in 1:max_outer
        _diff = 0f0
        for col in 1:n
            δx = abs(newton_bisection!(sa, r, w, col, δ, λ,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv))
            if !iszero(w[col])
                rc = δx / abs(w[col])
                rc > _diff && (_diff = rc)
            end
        end
        push!(lls, poisson_ll_from_residuals(r, yy, m))
        push!(sses, sse_objective(r, m))
        push!(max_diffs, Float64(_diff))
        _diff < rel_conv && break
    end
    return w, lls, sses, max_diffs
end

# ── GLM reference solutions ─────────────────────────────────────

function solve_glm_poisson(A::Matrix{Float32}, y::Vector{Float32})
    n = size(A, 2)
    X = Float64.(A)
    yy = Float64.(y)
    try
        model = glm(X, yy, Poisson(), IdentityLink();
                    dropcollinear=true, maxiter=500, atol=1e-8, rtol=1e-8)
        w = zeros(Float32, n)
        w .= Float32.(coef(model))
        return w, model
    catch e
        @warn "GLM Poisson failed: $e"
        return nothing, nothing
    end
end

function solve_glm_ols(A::Matrix{Float32}, y::Vector{Float32})
    n = size(A, 2)
    X = Float64.(A)
    yy = Float64.(y)
    try
        model = lm(X, yy; dropcollinear=true)
        w = zeros(Float32, n)
        w .= Float32.(coef(model))
        return w, model
    catch e
        @warn "GLM OLS failed: $e"
        return nothing, nothing
    end
end

# ── Unified comparison table ────────────────────────────────────

function unified_comparison_table(solvers, A, y, x_true, active_cols)
    n = length(x_true)
    n_active = length(active_cols)

    println("\n── Model-Fair Evaluation (all solvers × both likelihoods) ──")
    println()
    @printf("  %-18s  %12s  %12s  %14s  %6s  %4s (truth=%d)\n",
            "Solver", "Poisson_LL", "Poisson_Dev", "SSE", "Corr", "nnz", n_active)
    println("  " * "-"^78)

    for (label, w) in solvers
        w_n = Float32.(w[1:n])
        μ = compute_mu_dense(A, w_n)

        pll  = poisson_ll_from_mu64(μ, y)
        pdev = poisson_deviance(μ, y)
        sse  = sse_from_mu(μ, y)

        w_act = Float64.(w_n[active_cols])
        x_act = Float64.(x_true[active_cols])
        cor_val = cor(w_act, x_act)
        nnz = count(w_n .> 0.1f0)

        @printf("  %-18s  %12.2f  %12.2f  %14.2f  %6.3f  %4d\n",
                label, pll, pdev, sse, cor_val, nnz)
    end
end

# ── Main benchmark ───────────────────────────────────────────────

function main()
    println("=" ^ 80)
    println("  Poisson (Fisher) vs Poisson (Observed Hessian) vs Huber (≈OLS)")
    println("  m=2000 rows, n=500 cols, 2%% density, 80 active columns")
    println("=" ^ 80)

    A, x_true, y, active_cols = generate_problem()
    sa = build_SparseArray(A, y)

    mu_true = A * x_true
    @printf("\nProblem stats:\n")
    @printf("  A size: %d × %d, nnz = %d (%.1f%% dense)\n",
            size(A,1), size(A,2), count(!iszero, A), 100*count(!iszero, A)/length(A))
    @printf("  Active columns: %d, y range: [%.0f, %.0f], mean(y): %.1f\n",
            length(active_cols), minimum(y), maximum(y), mean(y))
    @printf("  True μ range: [%.1f, %.1f]\n", minimum(mu_true), maximum(mu_true))

    # ── Poisson Fisher ───────────────────────────────────────────
    println("\n── Poisson MLE (Fisher scoring) ──")
    t_fisher = @elapsed begin
        w_fisher, lls_fisher, sses_fisher, diffs_fisher = solve_poisson_fisher_instrumented(sa, y; max_outer=500)
    end
    n_fisher = length(diffs_fisher)
    @printf("  Converged in %d outer iterations, %.3f sec\n", n_fisher, t_fisher)
    @printf("  Final LL = %.2f,  Final SSE = %.2f\n", lls_fisher[end], sses_fisher[end])
    ll_mono_fisher = all(lls_fisher[i+1] >= lls_fisher[i] - 1e-4 for i in 1:length(lls_fisher)-1)
    @printf("  LL monotonic: %s\n", ll_mono_fisher)

    # ── Poisson Observed Hessian ─────────────────────────────────
    println("\n── Poisson MLE (Observed Hessian) ──")
    t_obs = @elapsed begin
        w_obs, lls_obs, sses_obs, diffs_obs = solve_poisson_observed_instrumented(sa, y; max_outer=500)
    end
    n_obs = length(diffs_obs)
    @printf("  Converged in %d outer iterations, %.3f sec\n", n_obs, t_obs)
    @printf("  Final LL = %.2f,  Final SSE = %.2f\n", lls_obs[end], sses_obs[end])
    ll_mono_obs = all(lls_obs[i+1] >= lls_obs[i] - 1e-4 for i in 1:length(lls_obs)-1)
    @printf("  LL monotonic: %s\n", ll_mono_obs)

    # ── Poisson MM (Cyclops-style) ─────────────────────────────
    println("\n── Poisson MLE (MM / Cyclops-style) ──")
    t_mm = @elapsed begin
        w_mm, lls_mm, sses_mm, diffs_mm = solve_poisson_mm_instrumented(sa, y; max_outer=500)
    end
    n_mm = length(diffs_mm)
    @printf("  Converged in %d outer iterations, %.3f sec\n", n_mm, t_mm)
    @printf("  Final LL = %.2f,  Final SSE = %.2f\n", lls_mm[end], sses_mm[end])
    ll_mono_mm = all(lls_mm[i+1] >= lls_mm[i] - 1e-4 for i in 1:length(lls_mm)-1)
    @printf("  LL monotonic: %s\n", ll_mono_mm)

    # ── Huber (≈OLS) ────────────────────────────────────────────
    println("\n── Huber (δ=1e10 ≈ OLS) ──")
    t_huber = @elapsed begin
        w_huber, lls_huber, sses_huber, diffs_huber = solve_huber_instrumented(sa, y; max_outer=500, δ=1f10)
    end
    n_huber = length(diffs_huber)
    @printf("  Converged in %d outer iterations, %.3f sec\n", n_huber, t_huber)
    @printf("  Final LL = %.2f,  Final SSE = %.2f\n", lls_huber[end], sses_huber[end])
    sse_mono = all(sses_huber[i+1] <= sses_huber[i] + 1e-2 for i in 1:length(sses_huber)-1)
    @printf("  SSE monotonic: %s\n", sse_mono)

    # ── GLM reference ────────────────────────────────────────────
    println("\n── GLM.jl reference solutions ──")

    t_glm_pois = @elapsed begin
        w_glm_poisson, _ = solve_glm_poisson(A, y)
    end
    if w_glm_poisson !== nothing
        @printf("  GLM Poisson (identity link): %.3f sec\n", t_glm_pois)
    else
        @printf("  GLM Poisson: FAILED (rank-deficient)\n")
    end

    t_glm_ols = @elapsed begin
        w_glm_ols, _ = solve_glm_ols(A, y)
    end
    if w_glm_ols !== nothing
        @printf("  GLM OLS: %.3f sec\n", t_glm_ols)
    else
        @printf("  GLM OLS: FAILED\n")
    end

    # ── Model-fair evaluation table ──────────────────────────────
    solvers = Pair{String, Vector{Float32}}[]
    push!(solvers, "Poisson-Fisher"  => w_fisher)
    push!(solvers, "Poisson-ObsHess" => w_obs)
    push!(solvers, "Poisson-MM"      => w_mm)
    push!(solvers, "Huber (≈OLS)"    => w_huber)
    if w_glm_poisson !== nothing
        push!(solvers, "GLM-Poisson" => w_glm_poisson)
    end
    if w_glm_ols !== nothing
        push!(solvers, "GLM-OLS" => w_glm_ols)
    end

    unified_comparison_table(solvers, A, y, x_true, active_cols)

    # ── Convergence trace (both objectives for all solvers) ──────
    println("\n── Convergence trace (both objectives for all solvers) ──")
    println("  iter   Fisher_LL      ObsHess_LL     MM_LL          Huber_LL       Fisher_SSE     ObsHess_SSE    MM_SSE         Huber_SSE")
    println("  " * "-"^130)
    max_show = max(length(diffs_fisher), length(diffs_obs), length(diffs_mm), length(diffs_huber))
    show_iters = sort(unique(vcat(
        collect(1:min(15, max_show)),
        collect(20:10:max_show),
        collect(max(1, max_show-4):max_show)
    )))
    for i in show_iters
        fl = i+1 <= length(lls_fisher) ? @sprintf("%14.2f", lls_fisher[i+1]) : "              "
        ol = i+1 <= length(lls_obs)    ? @sprintf("%14.2f", lls_obs[i+1])    : "              "
        ml = i+1 <= length(lls_mm)     ? @sprintf("%14.2f", lls_mm[i+1])     : "              "
        hl = i+1 <= length(lls_huber)  ? @sprintf("%14.2f", lls_huber[i+1])  : "              "
        fs = i+1 <= length(sses_fisher) ? @sprintf("%14.2f", sses_fisher[i+1]) : "              "
        os = i+1 <= length(sses_obs)    ? @sprintf("%14.2f", sses_obs[i+1])    : "              "
        ms = i+1 <= length(sses_mm)     ? @sprintf("%14.2f", sses_mm[i+1])     : "              "
        hs = i+1 <= length(sses_huber)  ? @sprintf("%14.2f", sses_huber[i+1])  : "              "
        @printf("  %4d  %s %s %s %s %s %s %s %s\n", i, fl, ol, ml, hl, fs, os, ms, hs)
    end

    # ── Weight-level convergence diffs ───────────────────────────
    println("\n── Weight convergence (max relative change per iteration) ──")
    println("  iter   Fisher_diff     ObsHess_diff    MM_diff         Huber_diff")
    println("  " * "-"^70)
    for i in show_iters
        fd = i <= length(diffs_fisher) ? @sprintf("%.4e", diffs_fisher[i]) : "  converged  "
        od = i <= length(diffs_obs)    ? @sprintf("%.4e", diffs_obs[i])    : "  converged  "
        md = i <= length(diffs_mm)     ? @sprintf("%.4e", diffs_mm[i])     : "  converged  "
        hd = i <= length(diffs_huber)  ? @sprintf("%.4e", diffs_huber[i])  : "  converged  "
        @printf("  %4d   %s  %s  %s  %s\n", i, fd, od, md, hd)
    end

    # ── Timing summary ───────────────────────────────────────────
    println("\n── Timing summary ──")
    @printf("  Poisson Fisher:   %3d iters, %.4f sec  (%.5f sec/iter)\n",
            n_fisher, t_fisher, t_fisher / max(n_fisher, 1))
    @printf("  Poisson ObsHess:  %3d iters, %.4f sec  (%.5f sec/iter)\n",
            n_obs, t_obs, t_obs / max(n_obs, 1))
    @printf("  Poisson MM:       %3d iters, %.4f sec  (%.5f sec/iter)\n",
            n_mm, t_mm, t_mm / max(n_mm, 1))
    @printf("  Huber (≈OLS):     %3d iters, %.4f sec  (%.5f sec/iter)\n",
            n_huber, t_huber, t_huber / max(n_huber, 1))
    if w_glm_poisson !== nothing
        @printf("  GLM Poisson:      %.4f sec\n", t_glm_pois)
    else
        @printf("  GLM Poisson:      FAILED\n")
    end
    if w_glm_ols !== nothing
        @printf("  GLM OLS:          %.4f sec\n", t_glm_ols)
    else
        @printf("  GLM OLS:          FAILED\n")
    end

    println("\n" * "=" ^ 80)
    println("  Done.")
    println("=" ^ 80)
end

main()
