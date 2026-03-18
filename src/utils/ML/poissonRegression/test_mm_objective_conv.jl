# Experiment 3: Cyclops-style objective-based convergence
#
# Instead of converging on max relative weight change (which fires too early
# for MM's small steps), converge on relative change in the objective (LL):
#   |LL_new - LL_old| / (|LL_new| + 1) < tol
#
# This is exactly what Cyclops uses (Suchard et al. 2013).
#
# Run:  julia --project=<Pioneer root> test_mm_objective_conv.jl

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

nr_iter  = data[:max_iter_newton]
bs_iter  = data[:max_iter_bisection]
max_outer = data[:max_iter_outer]
nr_acc   = data[:accuracy_newton]
bs_acc   = data[:accuracy_bisection]
rel_conv = data[:max_diff]

println("=" ^ 90)
println("  Objective-Based Convergence Experiment: scan342335")
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

# ── MM with OBJECTIVE-BASED convergence (Cyclops-style) ───────────
function solve_mm_objconv(sa_l, init_weights;
                          max_inner=5, obj_tol=1e-6,
                          mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)

    lls   = Float64[]
    diffs = Float64[]

    old_ll = poissonLogLikelihood(μ, yy, m)
    push!(lls, old_ll)

    for iter in 1:mm_max_outer
        _diff = Float32(0)
        max_weight = Float32(0)
        weight_floor = iter > 5 ? Float32(0) : Float32(0)  # compute below

        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, yy, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                if iszero(w[col]) || (!iszero(X0) && abs_step / abs(X0) < 1f-4)
                    break
                end
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if !iszero(w[col])
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end

        new_ll = poissonLogLikelihood(μ, yy, m)
        push!(lls, new_ll)
        push!(diffs, Float64(_diff))

        # Cyclops convergence: relative change in objective
        obj_change = abs(new_ll - old_ll) / (abs(new_ll) + 1.0)
        if obj_change < obj_tol && iter > 1
            break
        end
        old_ll = new_ll
    end

    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, length(diffs), sse, pll, nnz, lls, diffs
end

# ── MM with WEIGHT-BASED convergence (current approach) ───────────
function solve_mm_weightconv(sa_l, init_weights;
                             max_inner=5, wt_tol=rel_conv,
                             mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)

    lls   = Float64[]
    diffs = Float64[]
    push!(lls, poissonLogLikelihood(μ, yy, m))

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
                if iszero(w[col]) || (!iszero(X0) && abs_step / abs(X0) < wt_tol)
                    break
                end
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        push!(lls, poissonLogLikelihood(μ, yy, m))
        push!(diffs, Float64(_diff))
        _diff < wt_tol && break
    end

    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, length(diffs), sse, pll, nnz, lls, diffs
end

# ── Reference: ObsHess Newton+bisection ───────────────────────────
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

# ══════════════════════════════════════════════════════════════════
#  Run experiments
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))

# ── Baselines ─────────────────────────────────────────────────────
println("\n── Baselines (cold-start) ──")
@printf("  %-35s  %5s  %8s  %16s  %14s  %5s\n",
        "Solver", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 90)

t_o = @elapsed w_o, it_o, sse_o, pll_o, nnz_o = solve_obshess(sa, cold_w)
@printf("  %-35s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "ObsHess (Newton+bisection)", it_o, t_o, pll_o, sse_o, nnz_o)

t_wt = @elapsed w_wt, it_wt, sse_wt, pll_wt, nnz_wt, _, _ = solve_mm_weightconv(sa, cold_w)
@printf("  %-35s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "MM K=5 (weight conv, rel=0.01)", it_wt, t_wt, pll_wt, sse_wt, nnz_wt)

# ── Objective-based convergence with varying tolerance ────────────
println("\n── MM K=5 with objective-based convergence (cold-start) ──")
@printf("  %-35s  %5s  %8s  %16s  %14s  %5s\n",
        "Objective tolerance", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 90)

for tol in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10]
    local t = @elapsed local w_t, it_t, sse_t, pll_t, nnz_t, _, _ = solve_mm_objconv(sa, cold_w; obj_tol=tol)
    @printf("  %-35s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            @sprintf("obj_tol=%.0e", tol), it_t, t, pll_t, sse_t, nnz_t)
end

# ── Convergence trace for best objective-conv ─────────────────────
println("\n── Convergence trace: obj_tol=1e-8 (cold-start, K=5) ──")
w_best, it_best, sse_best, pll_best, nnz_best, lls_best, diffs_best = solve_mm_objconv(sa, cold_w; obj_tol=1e-8)

@printf("  %5s  %16s  %16s  %12s\n", "Iter", "Poisson LL", "ΔLL_rel", "max_wt_diff")
println("  " * "─" ^ 55)
n_trace = length(diffs_best)
show_iters = sort(unique(vcat(
    collect(1:min(30, n_trace)),
    collect(35:5:min(60, n_trace)),
    collect(60:10:min(100, n_trace)),
    collect(100:50:min(500, n_trace)),
    collect(500:100:n_trace),
    [n_trace]
)))
for i in show_iters
    i > n_trace && continue
    ll_rel = i > 0 ? abs(lls_best[i+1] - lls_best[i]) / (abs(lls_best[i+1]) + 1.0) : 0.0
    @printf("  %5d  %16.6e  %16.6e  %12.4e\n",
            i, lls_best[i+1], ll_rel, diffs_best[i])
end

ll_mono = all(lls_best[i+1] >= lls_best[i] - 1.0 for i in 1:length(lls_best)-1)
@printf("\n  LL monotonic: %s\n", ll_mono)

# ── Also test varying inner steps with objective convergence ──────
println("\n── Varying inner steps K with obj_tol=1e-6 (cold-start) ──")
@printf("  %-25s  %5s  %8s  %16s  %14s  %5s\n",
        "Config", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 80)

for ki in [1, 3, 5, 10, 25]
    local t = @elapsed local w_t, it_t, sse_t, pll_t, nnz_t, _, _ = solve_mm_objconv(sa, cold_w; max_inner=ki, obj_tol=1e-6)
    @printf("  %-25s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            "K=$ki, obj_tol=1e-6", it_t, t, pll_t, sse_t, nnz_t)
end

# ── Final comparison ──────────────────────────────────────────────
println("\n── Final comparison vs ObsHess (cold-start) ──")
w_final, it_final, sse_final, pll_final, nnz_final, _, _ = solve_mm_objconv(sa, cold_w; max_inner=5, obj_tol=1e-6)
wo_f = Float64.(w_o[1:sa.n]); wf_f = Float64.(w_final[1:sa.n])
@printf("  ObsHess:             iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_o, pll_o, sse_o, nnz_o)
@printf("  MM(K=5,obj=1e-6):    iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_final, pll_final, sse_final, nnz_final)
@printf("  Weight correlation:  %.6f\n", cor(wo_f, wf_f))

# ── Warm-start comparison ─────────────────────────────────────────
warm_w = copy(data[:weights])
println("\n── Warm-start comparison ──")

t_ow = @elapsed w_ow, it_ow, sse_ow, pll_ow, nnz_ow = solve_obshess(sa, warm_w)
@printf("  ObsHess:             iters=%d, %.4f sec, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_ow, t_ow, pll_ow, sse_ow, nnz_ow)

t_mw = @elapsed w_mw, it_mw, sse_mw, pll_mw, nnz_mw, _, _ = solve_mm_objconv(sa, warm_w; max_inner=5, obj_tol=1e-6)
@printf("  MM(K=5,obj=1e-6):    iters=%d, %.4f sec, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_mw, t_mw, pll_mw, sse_mw, nnz_mw)

ww_oh = Float64.(w_ow[1:sa.n]); ww_mm = Float64.(w_mw[1:sa.n])
@printf("  Weight correlation:  %.6f\n", cor(ww_oh, ww_mm))

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
