## Compare Huber, Poisson-Fisher, Poisson-ObsHessian on real scan342335 data
## Run:  julia --project=<Pioneer root> test_real_data_comparison.jl

using Printf, Statistics, Serialization

# Deserialize needs Pioneer's SparseArray type in scope
using Pioneer

# Now load standalone solver implementations (they define their own local SparseArray)
# We must do this AFTER Pioneer so the local SparseArray shadows Pioneer's for the solvers
include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Load real data ────────────────────────────────────────────────
data = deserialize("/Users/n.t.wamsley/Desktop/solveHuber_inputs_scan342335.jls")

# Reconstruct using the LOCAL SparseArray (from SparseArray.jl include)
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
println("  Real Production Data: scan342335")
println("  m=$(sa.m) rows, n=$(sa.n) cols, n_vals=$(sa.n_vals)")
println("  δ=$δ_param, λ=$λ_param, max_outer=$max_outer, rel_conv=$rel_conv")
println("=" ^ 90)

# ── Residual init for Huber ───────────────────────────────────────
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

# ── Evaluation functions ──────────────────────────────────────────
function sse_from_residuals(r::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += Float64(r[i])^2; end; s
end

function poisson_ll_from_residuals(r::Vector{Float32}, y::Vector{Float32}, m::Int)
    ll = 0.0
    @inbounds for i in 1:m
        μ_i = max(Float64(y[i]) + Float64(r[i]), 1e-6)
        ll += Float64(y[i]) * log(μ_i) - μ_i
    end; ll
end

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

# ── Instrumented solvers (all use significance floor) ─────────────

function solve_huber_real(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, max(m, length(data[:residuals])))
    w = copy(init_weights)

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
    sse = sse_from_residuals(r, m)
    pll = poisson_ll_from_residuals(r, y, m)
    return w, iters, sse, pll
end

function solve_poisson_fisher_real(sa_l, init_weights)
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
            δx = abs(newton_bisection_poisson!(sa_l, μ, yy, w, col,
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
    return w, iters, sse, pll
end

function solve_poisson_obs_real(sa_l, init_weights)
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
    return w, iters, sse, pll
end

function solve_poisson_mm_real(sa_l, init_weights; mm_rel_conv=rel_conv, mm_max_outer=max_outer, max_inner=5)
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
    return w, iters, sse, pll
end

function solve_poisson_mm_real_traced(sa_l, init_weights; mm_rel_conv=rel_conv, mm_max_outer=max_outer, max_inner=5)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)

    lls  = Float64[poisson_ll_from_mu(μ, yy, m)]
    sses = Float64[sse_from_mu(μ, yy, m)]
    diffs = Float64[]

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
        push!(lls, poisson_ll_from_mu(μ, yy, m))
        push!(sses, sse_from_mu(μ, yy, m))
        push!(diffs, Float64(_diff))
        _diff < mm_rel_conv && break
    end
    nnz = count(w[1:n] .> 0.01)
    return w, length(diffs), lls, sses, diffs, nnz
end

# ── Run comparisons ──────────────────────────────────────────────

warm_w = copy(data[:weights])
cold_w = ones(Float32, length(data[:weights]))

n_nz = count(warm_w[1:sa.n] .> 0)
@printf("\nWarm-start weight stats (n=%d):\n", sa.n)
@printf("  nonzero: %d, max: %.4e, median(nonzero): %.4e\n",
        n_nz, maximum(warm_w[1:sa.n]),
        median(warm_w[1:sa.n][warm_w[1:sa.n] .> 0]))

println("\n" * "─" ^ 90)
println("  WARM-START (pre-converged weights from production)")
println("─" ^ 90)

t1 = @elapsed w_hw, it_hw, sse_hw, pll_hw = solve_huber_real(sa, warm_w)
@printf("  Huber:          %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_hw, t1, sse_hw, pll_hw)

t2 = @elapsed w_fw, it_fw, sse_fw, pll_fw = solve_poisson_fisher_real(sa, warm_w)
@printf("  Poisson-Fisher: %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_fw, t2, sse_fw, pll_fw)

t3 = @elapsed w_ow, it_ow, sse_ow, pll_ow = solve_poisson_obs_real(sa, warm_w)
@printf("  Poisson-ObsHes: %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_ow, t3, sse_ow, pll_ow)

t3b = @elapsed w_mw, it_mw, sse_mw, pll_mw = solve_poisson_mm_real(sa, warm_w)
@printf("  Poisson-MM:     %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_mw, t3b, sse_mw, pll_mw)

println("\n" * "─" ^ 90)
println("  COLD-START (ones initialization)")
println("─" ^ 90)

t4 = @elapsed w_hc, it_hc, sse_hc, pll_hc = solve_huber_real(sa, cold_w)
@printf("  Huber:          %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_hc, t4, sse_hc, pll_hc)

t5 = @elapsed w_fc, it_fc, sse_fc, pll_fc = solve_poisson_fisher_real(sa, cold_w)
@printf("  Poisson-Fisher: %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_fc, t5, sse_fc, pll_fc)

t6 = @elapsed w_oc, it_oc, sse_oc, pll_oc = solve_poisson_obs_real(sa, cold_w)
@printf("  Poisson-ObsHes: %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_oc, t6, sse_oc, pll_oc)

t6b = @elapsed w_mc, it_mc, sse_mc, pll_mc = solve_poisson_mm_real(sa, cold_w)
@printf("  Poisson-MM:     %4d iters, %.4f sec, SSE=%.6e, PoissonLL=%.6e\n", it_mc, t6b, sse_mc, pll_mc)

# ── Summary table ────────────────────────────────────────────────
println("\n" * "=" ^ 90)
println("  SUMMARY TABLE")
println("=" ^ 90)
@printf("\n  %-24s  %5s  %8s  %14s  %16s  %5s\n",
        "Solver", "Iters", "Time(s)", "SSE", "Poisson LL", "nnz")
println("  " * "─" ^ 80)

for (label, iters, t, sse, pll, w) in [
    ("Huber (warm)",          it_hw, t1, sse_hw, pll_hw, w_hw),
    ("Huber (cold)",          it_hc, t4, sse_hc, pll_hc, w_hc),
    ("Poisson-Fisher (warm)", it_fw, t2, sse_fw, pll_fw, w_fw),
    ("Poisson-Fisher (cold)", it_fc, t5, sse_fc, pll_fc, w_fc),
    ("Poisson-ObsHes (warm)", it_ow, t3, sse_ow, pll_ow, w_ow),
    ("Poisson-ObsHes (cold)", it_oc, t6, sse_oc, pll_oc, w_oc),
    ("Poisson-MM (warm)",     it_mw, t3b, sse_mw, pll_mw, w_mw),
    ("Poisson-MM (cold)",     it_mc, t6b, sse_mc, pll_mc, w_mc),
]
    nnz = count(w[1:sa.n] .> 0.01)
    @printf("  %-24s  %5d  %8.4f  %14.6e  %16.6e  %5d\n", label, iters, t, sse, pll, nnz)
end

# ── Weight correlations ──────────────────────────────────────────
println("\n" * "─" ^ 90)
println("  Cross-solver weight correlations (warm-start)")
println("─" ^ 90)
wh = Float64.(w_hw[1:sa.n]); wf = Float64.(w_fw[1:sa.n]); wo = Float64.(w_ow[1:sa.n]); wm = Float64.(w_mw[1:sa.n])
@printf("  Huber vs Fisher:  %.6f\n", cor(wh, wf))
@printf("  Huber vs ObsHes:  %.6f\n", cor(wh, wo))
@printf("  Huber vs MM:      %.6f\n", cor(wh, wm))
@printf("  Fisher vs ObsHes: %.6f\n", cor(wf, wo))
@printf("  Fisher vs MM:     %.6f\n", cor(wf, wm))
@printf("  ObsHes vs MM:     %.6f\n", cor(wo, wm))

println("\n" * "─" ^ 90)
println("  Cross-solver weight correlations (cold-start)")
println("─" ^ 90)
whc = Float64.(w_hc[1:sa.n]); wfc = Float64.(w_fc[1:sa.n]); woc = Float64.(w_oc[1:sa.n]); wmc = Float64.(w_mc[1:sa.n])
@printf("  Huber vs Fisher:  %.6f\n", cor(whc, wfc))
@printf("  Huber vs ObsHes:  %.6f\n", cor(whc, woc))
@printf("  Huber vs MM:      %.6f\n", cor(whc, wmc))
@printf("  Fisher vs ObsHes: %.6f\n", cor(wfc, woc))
@printf("  Fisher vs MM:     %.6f\n", cor(wfc, wmc))
@printf("  ObsHes vs MM:     %.6f\n", cor(woc, wmc))

println("\n" * "─" ^ 90)
println("  Warm vs cold consistency (per solver)")
println("─" ^ 90)
@printf("  Huber:          cor=%.6f\n", cor(wh, whc))
@printf("  Poisson-Fisher: cor=%.6f\n", cor(wf, wfc))
@printf("  Poisson-ObsHes: cor=%.6f\n", cor(wo, woc))
@printf("  Poisson-MM:     cor=%.6f\n", cor(wm, wmc))

# ── MM convergence investigation ──────────────────────────────────
println("\n" * "=" ^ 90)
println("  MM CONVERGENCE INVESTIGATION (multi-step)")
println("=" ^ 90)

# Test varying inner step counts (cold-start)
println("\n── MM with varying inner steps (cold-start, rel_conv=$(rel_conv)) ──")
@printf("  %-22s  %5s  %8s  %14s  %16s  %5s\n",
        "Inner steps", "Iters", "Time(s)", "SSE", "Poisson LL", "nnz")
println("  " * "─" ^ 80)

for ki in [1, 3, 5, 10, 25]
    t = @elapsed w_t, it_t, sse_t, pll_t = solve_poisson_mm_real(sa, cold_w;
                                                max_inner=ki)
    nnz_t = count(w_t[1:sa.n] .> 0.01)
    @printf("  %-22s  %5d  %8.4f  %14.6e  %16.6e  %5d\n",
            "K=$ki", it_t, t, sse_t, pll_t, nnz_t)
end

# Also warm-start
println("\n── MM with varying inner steps (warm-start, rel_conv=$(rel_conv)) ──")
@printf("  %-22s  %5s  %8s  %14s  %16s  %5s\n",
        "Inner steps", "Iters", "Time(s)", "SSE", "Poisson LL", "nnz")
println("  " * "─" ^ 80)

for ki in [1, 3, 5, 10, 25]
    t = @elapsed w_t, it_t, sse_t, pll_t = solve_poisson_mm_real(sa, warm_w;
                                                max_inner=ki)
    nnz_t = count(w_t[1:sa.n] .> 0.01)
    @printf("  %-22s  %5d  %8.4f  %14.6e  %16.6e  %5d\n",
            "K=$ki", it_t, t, sse_t, pll_t, nnz_t)
end

# Convergence trace for best MM config
println("\n── MM convergence trace (cold-start, K=5, rel_conv=1e-4) ──")
_, _, mm_lls, mm_sses, mm_diffs, mm_nnz = solve_poisson_mm_real_traced(sa, cold_w;
                                            mm_rel_conv=1f-4, mm_max_outer=5000, max_inner=5)

@printf("  %5s  %16s  %14s  %12s\n", "Iter", "Poisson LL", "SSE", "max_rel_diff")
println("  " * "─" ^ 55)
n_trace = length(mm_diffs)
show_iters = sort(unique(vcat(
    collect(1:min(20, n_trace)),
    collect(25:5:min(50, n_trace)),
    collect(50:25:min(200, n_trace)),
    collect(200:100:n_trace),
    [n_trace]
)))
for i in show_iters
    i > n_trace && continue
    @printf("  %5d  %16.6e  %14.6e  %12.4e\n",
            i, mm_lls[i+1], mm_sses[i+1], mm_diffs[i])
end

# LL monotonicity check
ll_mono = all(mm_lls[i+1] >= mm_lls[i] - 1.0 for i in 1:length(mm_lls)-1)
@printf("\n  LL monotonic: %s\n", ll_mono)
@printf("  Final nnz: %d / %d\n", mm_nnz, sa.n)

# Compare best MM vs ObsHess
println("\n── MM (K=5) vs ObsHess (cold-start) ──")
w_mm_best, it_mm_best, sse_mm_best, pll_mm_best = solve_poisson_mm_real(sa, cold_w;
                                                    mm_rel_conv=1f-4, mm_max_outer=5000, max_inner=5)
wm_best = Float64.(w_mm_best[1:sa.n])
@printf("  ObsHess:  iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_oc, pll_oc, sse_oc, count(w_oc[1:sa.n] .> 0.01))
@printf("  MM(K=5):  iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_mm_best, pll_mm_best, sse_mm_best, count(w_mm_best[1:sa.n] .> 0.01))
@printf("  ObsHess vs MM(K=5) weight correlation: %.6f\n", cor(woc, wm_best))

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
