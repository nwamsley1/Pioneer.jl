# Experiment: Can Float64 derivatives or OLS warm-start fix MM on production data?
#
# Run:  julia --project=<Pioneer root> test_mm_fixes.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Load real data ────────────────────────────────────────────────
data = deserialize("/Users/n.t.wamsley/Desktop/solveHuber_inputs_scan342335.jls")

sa = Main.SparseArray(
    data[:n_vals], data[:m], data[:n],
    Vector{Int64}(data[:rowval]), Vector{UInt16}(data[:colval]), data[:nzval],
    data[:matched], data[:isotope], data[:x], Vector{Int64}(data[:colptr])
)

δ_param  = data[:delta]
λ_param  = data[:lambda]
nr_iter  = data[:max_iter_newton]
bs_iter  = data[:max_iter_bisection]
max_outer = data[:max_iter_outer]
nr_acc   = data[:accuracy_newton]
bs_acc   = data[:accuracy_bisection]
rel_conv = data[:max_diff]

println("=" ^ 90)
println("  MM Fix Experiments: scan342335")
println("  m=$(sa.m), n=$(sa.n), rel_conv=$rel_conv")
println("=" ^ 90)

# ── Eval functions ────────────────────────────────────────────────
function sse_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += (Float64(y[i]) - Float64(μ[i]))^2; end; s
end
function poisson_ll_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    ll = 0.0
    @inbounds for i in 1:m
        μ_i = max(Float64(μ[i]), 1e-6)
        ll += Float64(y[i]) * log(μ_i) - μ_i
    end; ll
end

function initResiduals_plain!(r::Vector{T}, sa_l::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds for i in 1:sa_l.m; r[i] = zero(T); end
    @inbounds for n in 1:sa_l.n_vals
        if iszero(r[sa_l.rowval[n]]); r[sa_l.rowval[n]] = -sa_l.x[n]; end
    end
    @inbounds for col in 1:sa_l.n
        for n in sa_l.colptr[col]:(sa_l.colptr[col+1] - 1)
            r[sa_l.rowval[n]] += w[col] * sa_l.nzval[n]
        end
    end
end

# ══════════════════════════════════════════════════════════════════
# Experiment 1: Float64 derivative accumulation
# ══════════════════════════════════════════════════════════════════

"""Float64 derivatives — same math, higher precision accumulation."""
function getPoissonDerivativesObs_f64!(Hs::SparseArray{Ti, T},
                                       μ::Vector{T}, y::Vector{T},
                                       col::Int64) where {Ti<:Integer, T<:AbstractFloat}
    L1 = zero(Float64)
    L2 = zero(Float64)
    ε  = Float64(POISSON_MU_FLOOR)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row   = Hs.rowval[i]
        a_ij  = Float64(Hs.nzval[i])
        μ_i   = max(Float64(μ[row]), ε)
        y_i   = Float64(y[row])
        L1   += a_ij * (1.0 - y_i / μ_i)
        L2   += a_ij * a_ij * y_i / (μ_i * μ_i)
    end
    if L2 < ε
        @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            a_ij = Float64(Hs.nzval[i])
            μ_i  = max(Float64(μ[Hs.rowval[i]]), ε)
            L2  += a_ij * a_ij / μ_i
        end
    end
    return L1, L2
end

function solve_mm_f64_deriv(sa_l, init_weights; max_inner=5, mm_rel_conv=rel_conv, mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    ε = 1e-6  # Float64
    max_weight = Float32(0)
    iters = 0
    for iter in 1:mm_max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs_f64!(sa_l, μ, yy, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                # Compute step in Float64, apply in Float32
                w[col] = Float32(max(Float64(w[col]) - L1 / L2, 0.0))
                updateMu!(sa_l, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                if iszero(w[col]) || (!iszero(X0) && abs_step / abs(X0) < mm_rel_conv)
                    break
                end
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < mm_rel_conv && break
    end
    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ══════════════════════════════════════════════════════════════════
# Experiment 2: OLS warm-start → Poisson MM
# ══════════════════════════════════════════════════════════════════

function solve_huber_N_iters(sa_l, init_weights, n_iters)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, m)
    w = copy(init_weights)
    initResiduals_plain!(r, sa_l, w)
    reg = NoNorm()
    for iter in 1:n_iters
        for col in 1:n
            newton_bisection!(sa_l, r, w, col, δ_param, λ_param,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv)
        end
    end
    return w
end

function solve_mm_from_weights(sa_l, init_weights; max_inner=5, mm_rel_conv=rel_conv, mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)
    iters = 0
    for iter in 1:mm_max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, yy, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                if iszero(w[col]) || (!iszero(X0) && abs_step / abs(X0) < mm_rel_conv)
                    break
                end
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < mm_rel_conv && break
    end
    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ── Reference solvers ─────────────────────────────────────────────
function solve_obshess(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)
    max_weight = Float32(0)
    iters = 0
    for iter in 1:max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            δx = abs(newton_bisection_poisson_obs!(sa_l, μ, yy, w, col,
                        nr_iter, bs_iter, nr_acc, bs_acc, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

function solve_huber_full(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, m); w = copy(init_weights)
    y = zeros(Float32, m)
    @inbounds for i in 1:sa_l.n_vals
        if iszero(y[sa_l.rowval[i]]); y[sa_l.rowval[i]] = sa_l.x[i]; end
    end
    initResiduals_plain!(r, sa_l, w)
    reg = NoNorm()
    max_weight = Float32(0)
    iters = 0
    for iter in 1:max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            δx = abs(newton_bisection!(sa_l, r, w, col, δ_param, λ_param,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    # Compute μ from weights for eval
    μ = zeros(Float32, m); yy = zeros(Float32, m)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)
    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ══════════════════════════════════════════════════════════════════
#  Run experiments
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))

# ── Baselines ─────────────────────────────────────────────────────
println("\n── Baselines (cold-start) ──")
@printf("  %-30s  %5s  %8s  %16s  %14s  %5s\n",
        "Solver", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 85)

t = @elapsed w_h, it_h, sse_h, pll_h, nnz_h = solve_huber_full(sa, cold_w)
@printf("  %-30s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "Huber", it_h, t, pll_h, sse_h, nnz_h)

t = @elapsed w_o, it_o, sse_o, pll_o, nnz_o = solve_obshess(sa, cold_w)
@printf("  %-30s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "ObsHess (Newton+bisection)", it_o, t, pll_o, sse_o, nnz_o)

t = @elapsed w_mm, it_mm, sse_mm, pll_mm, nnz_mm = solve_mm_from_weights(sa, cold_w)
@printf("  %-30s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "MM K=5 (cold)", it_mm, t, pll_mm, sse_mm, nnz_mm)

# ── Experiment 1: Float64 derivatives ─────────────────────────────
println("\n── Experiment 1: Float64 derivative accumulation (cold-start) ──")
@printf("  %-30s  %5s  %8s  %16s  %14s  %5s\n",
        "Solver", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 85)

for ki in [1, 5, 10, 25]
    t = @elapsed w_f, it_f, sse_f, pll_f, nnz_f = solve_mm_f64_deriv(sa, cold_w; max_inner=ki)
    @printf("  %-30s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            "MM-f64 K=$ki", it_f, t, pll_f, sse_f, nnz_f)
end

# ── Experiment 2: OLS warm-start → MM ────────────────────────────
println("\n── Experiment 2: Huber warm-start → Poisson MM (cold initial) ──")
@printf("  %-30s  %5s  %5s  %8s  %16s  %14s  %5s\n",
        "Solver", "OLS", "MM", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 90)

for n_ols in [1, 2, 3, 5, 10, 12]
    t = @elapsed begin
        w_warm = solve_huber_N_iters(sa, cold_w, n_ols)
        w_r, it_r, sse_r, pll_r, nnz_r = solve_mm_from_weights(sa, w_warm)
    end
    @printf("  %-30s  %5d  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            "Huber($n_ols)→MM(K=5)", n_ols, it_r, t, pll_r, sse_r, nnz_r)
end

# ── Experiment 2b: OLS warm-start → MM with Float64 ──────────────
println("\n── Experiment 2b: Huber warm-start → MM-f64 ──")
@printf("  %-30s  %5s  %5s  %8s  %16s  %14s  %5s\n",
        "Solver", "OLS", "MM", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 90)

for n_ols in [1, 3, 5, 10, 12]
    t = @elapsed begin
        w_warm = solve_huber_N_iters(sa, cold_w, n_ols)
        w_r, it_r, sse_r, pll_r, nnz_r = solve_mm_f64_deriv(sa, w_warm; max_inner=5)
    end
    @printf("  %-30s  %5d  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            "Huber($n_ols)→MM-f64(K=5)", n_ols, it_r, t, pll_r, sse_r, nnz_r)
end

# ── Cross-solver weight correlations ─────────────────────────────
println("\n── Weight correlations for best configs ──")
# Best OLS→MM
w_best_warm = solve_huber_N_iters(sa, cold_w, 5)
w_best, _, _, _, _ = solve_mm_from_weights(sa, w_best_warm)
w_best64_warm = solve_huber_N_iters(sa, cold_w, 5)
w_best64, _, _, _, _ = solve_mm_f64_deriv(sa, w_best64_warm; max_inner=5)

wb = Float64.(w_best[1:sa.n])
wb64 = Float64.(w_best64[1:sa.n])
wo = Float64.(w_o[1:sa.n])
wh = Float64.(w_h[1:sa.n])

@printf("  ObsHess vs Huber:              %.6f\n", cor(wo, wh))
@printf("  ObsHess vs Huber(5)→MM:        %.6f\n", cor(wo, wb))
@printf("  ObsHess vs Huber(5)→MM-f64:    %.6f\n", cor(wo, wb64))
@printf("  Huber vs Huber(5)→MM:          %.6f\n", cor(wh, wb))
@printf("  Huber(5)→MM vs Huber(5)→MM-f64:%.6f\n", cor(wb, wb64))

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
