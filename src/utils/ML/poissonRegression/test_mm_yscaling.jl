## Test MM solver with y-scaling on real scan342335 data
## Compares to Huber loss: runtime, iterations, inner steps, weight sparsity
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/test_mm_yscaling.jl

using Printf, Statistics, Serialization

# Need Pioneer's SparseArray type for deserialization
using Pioneer

# Local SparseArray + solvers shadow Pioneer's types
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
println("  MM Y-Scaling Test: scan342335")
println("  m=$(sa.m) rows, n=$(sa.n) cols, n_vals=$(sa.n_vals)")
println("  max_outer=$max_outer, rel_conv=$rel_conv")
println("=" ^ 90)

# ── Helpers ───────────────────────────────────────────────────────

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

function poisson_ll_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    ll = 0.0
    @inbounds for i in 1:m
        μ_i = max(Float64(μ[i]), 1e-6)
        ll += Float64(y[i]) * log(μ_i) - μ_i
    end; ll
end

function sse_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += (Float64(y[i]) - Float64(μ[i]))^2; end; s
end

# ── Huber solver (baseline) ──────────────────────────────────────

function solve_huber(sa_l, init_weights)
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
    return w, iters, y, r
end

# ── MM solver using the actual solvePoissonMM! with y-scaling ────

function solve_mm(sa_l, init_weights; max_inner=5, mm_max_outer=max_outer, mm_rel_conv=rel_conv)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)
    solvePoissonMM!(sa_l, μ, yy, w, Int64(mm_max_outer), mm_rel_conv;
                    max_inner_iter=Int64(max_inner))
    return w, μ, yy
end

# ── Instrumented MM solver to count actual outer iterations ──────

function solve_mm_counted(sa_l, init_weights; max_inner=5, mm_max_outer=max_outer, mm_rel_conv=rel_conv)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)

    # Y-scaling (mirror solvePoissonMM!)
    y_scale = Float32(0)
    @inbounds for i in 1:sa_l.m
        if yy[i] > y_scale; y_scale = yy[i]; end
    end
    if y_scale > Float32(1)
        @inbounds for i in 1:sa_l.m; yy[i] /= y_scale; end
        @inbounds for j in 1:sa_l.n; w[j] /= y_scale; end
        initMu!(μ, sa_l, w)
    end

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)
    iters = 0
    total_inner_steps = 0

    lls = Float64[poisson_ll_from_mu(μ, yy, m)]
    diffs = Float64[]

    for iter in 1:mm_max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-4) : Float32(0)
        max_weight = Float32(0)
        iter_inner = 0
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, yy, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
                iter_inner += 1
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
        total_inner_steps += iter_inner
        iters = iter
        push!(lls, poisson_ll_from_mu(μ, yy, m))
        push!(diffs, Float64(_diff))
        _diff < mm_rel_conv && break
    end

    # Unscale
    if y_scale > Float32(1)
        @inbounds for i in 1:sa_l.m; yy[i] *= y_scale; end
        @inbounds for j in 1:sa_l.n; w[j] *= y_scale; end
        initMu!(μ, sa_l, w)
    end

    return w, iters, total_inner_steps, μ, yy, lls, diffs
end

# ══════════════════════════════════════════════════════════════════
# Run tests
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))
warm_w = copy(data[:weights])

# Get y values for reference
yy_ref = zeros(Float32, sa.m)
initObserved!(yy_ref, sa)
y_max = maximum(yy_ref[1:sa.m])
y_min = minimum(yy_ref[i] for i in 1:sa.m if yy_ref[i] > 0)
@printf("\n  y range: min(nonzero)=%.4e, max=%.4e, dynamic range=%.1f OOM\n",
        y_min, y_max, log10(y_max/y_min))

# ── 1. Huber baseline ────────────────────────────────────────────
println("\n" * "=" ^ 90)
println("  1. HUBER BASELINE")
println("=" ^ 90)

# Warmup
solve_huber(sa, cold_w)

t_hc = @elapsed w_hc, it_hc, y_hc, r_hc = solve_huber(sa, cold_w)
μ_hc = zeros(Float32, sa.m)
@inbounds for i in 1:sa.m; μ_hc[i] = y_hc[i] + r_hc[i]; end
pll_hc = poisson_ll_from_mu(μ_hc, y_hc, sa.m)
sse_hc = sse_from_mu(μ_hc, y_hc, sa.m)

t_hw = @elapsed w_hw, it_hw, y_hw, r_hw = solve_huber(sa, warm_w)
μ_hw = zeros(Float32, sa.m)
@inbounds for i in 1:sa.m; μ_hw[i] = y_hw[i] + r_hw[i]; end
pll_hw = poisson_ll_from_mu(μ_hw, y_hw, sa.m)
sse_hw = sse_from_mu(μ_hw, y_hw, sa.m)

@printf("  Cold-start: %4d iters, %.4f sec, PoissonLL=%.6e, SSE=%.6e\n", it_hc, t_hc, pll_hc, sse_hc)
@printf("  Warm-start: %4d iters, %.4f sec, PoissonLL=%.6e, SSE=%.6e\n", it_hw, t_hw, pll_hw, sse_hw)

# ── 2. MM with y-scaling: varying inner iterations ───────────────
println("\n" * "=" ^ 90)
println("  2. MM WITH Y-SCALING: VARYING INNER STEPS")
println("=" ^ 90)

# Warmup
solve_mm(sa, cold_w; max_inner=5)

println("\n── Cold-start ──")
@printf("  %-12s  %5s  %8s  %10s  %16s  %14s\n",
        "Inner steps", "Iters", "Time(s)", "InnerTotal", "Poisson LL", "SSE")
println("  " * "─" ^ 80)

for ki in [1, 2, 3, 5, 10, 20]
    t = @elapsed w_t, it_t, inner_t, μ_t, yy_t, _, _ = solve_mm_counted(sa, cold_w; max_inner=ki)
    pll_t = poisson_ll_from_mu(μ_t, yy_t, sa.m)
    sse_t = sse_from_mu(μ_t, yy_t, sa.m)
    @printf("  K=%-9d  %5d  %8.4f  %10d  %16.6e  %14.6e\n",
            ki, it_t, t, inner_t, pll_t, sse_t)
end

println("\n── Warm-start ──")
@printf("  %-12s  %5s  %8s  %10s  %16s  %14s\n",
        "Inner steps", "Iters", "Time(s)", "InnerTotal", "Poisson LL", "SSE")
println("  " * "─" ^ 80)

for ki in [1, 2, 3, 5, 10, 20]
    t = @elapsed w_t, it_t, inner_t, μ_t, yy_t, _, _ = solve_mm_counted(sa, warm_w; max_inner=ki)
    pll_t = poisson_ll_from_mu(μ_t, yy_t, sa.m)
    sse_t = sse_from_mu(μ_t, yy_t, sa.m)
    @printf("  K=%-9d  %5d  %8.4f  %10d  %16.6e  %14.6e\n",
            ki, it_t, t, inner_t, pll_t, sse_t)
end

# ── 3. Detailed comparison: MM K=1 vs K=5 vs Huber ──────────────
println("\n" * "=" ^ 90)
println("  3. HEAD-TO-HEAD: Huber vs MM(K=1) vs MM(K=5)  [cold-start]")
println("=" ^ 90)

t_mm1 = @elapsed w_mm1, it_mm1, inner_mm1, μ_mm1, yy_mm1, lls_mm1, diffs_mm1 = solve_mm_counted(sa, cold_w; max_inner=1)
pll_mm1 = poisson_ll_from_mu(μ_mm1, yy_mm1, sa.m)
sse_mm1 = sse_from_mu(μ_mm1, yy_mm1, sa.m)

t_mm5 = @elapsed w_mm5, it_mm5, inner_mm5, μ_mm5, yy_mm5, lls_mm5, diffs_mm5 = solve_mm_counted(sa, cold_w; max_inner=5)
pll_mm5 = poisson_ll_from_mu(μ_mm5, yy_mm5, sa.m)
sse_mm5 = sse_from_mu(μ_mm5, yy_mm5, sa.m)

@printf("\n  %-20s  %5s  %8s  %10s  %16s  %14s\n",
        "Solver", "Iters", "Time(s)", "InnerTotal", "Poisson LL", "SSE")
println("  " * "─" ^ 80)
@printf("  %-20s  %5d  %8.4f  %10s  %16.6e  %14.6e\n",
        "Huber", it_hc, t_hc, "N/A", pll_hc, sse_hc)
@printf("  %-20s  %5d  %8.4f  %10d  %16.6e  %14.6e\n",
        "MM(K=1,yscale)", it_mm1, t_mm1, inner_mm1, pll_mm1, sse_mm1)
@printf("  %-20s  %5d  %8.4f  %10d  %16.6e  %14.6e\n",
        "MM(K=5,yscale)", it_mm5, t_mm5, inner_mm5, pll_mm5, sse_mm5)

# ── 4. Inner step utilization ─────────────────────────────────────
println("\n" * "=" ^ 90)
println("  4. INNER STEP UTILIZATION")
println("=" ^ 90)

for (ki, label) in [(1, "K=1"), (5, "K=5")]
    _, iters, inner_total, _, _, _, _ = solve_mm_counted(sa, cold_w; max_inner=ki)
    avg_inner = inner_total / (iters * sa.n)
    @printf("  %s: %d outer iters x %d cols = %d coord updates, %d inner steps total\n",
            label, iters, sa.n, iters * sa.n, inner_total)
    @printf("       avg inner steps per coord update: %.2f\n", avg_inner)
end

# ── 5. Weight sparsity analysis ──────────────────────────────────
println("\n" * "=" ^ 90)
println("  5. WEIGHT SPARSITY: fraction below 1e-5 * max_weight")
println("=" ^ 90)

for (label, w_vec) in [("Huber (cold)", w_hc),
                        ("MM K=1 (cold)", w_mm1),
                        ("MM K=5 (cold)", w_mm5)]
    ww = w_vec[1:sa.n]
    mw = maximum(ww)
    threshold = mw * 1e-5
    n_below = count(ww .< threshold)
    n_zero = count(iszero.(ww))
    n_total = sa.n
    @printf("  %-18s  max_w=%.4e, <%s: %d/%d (%.1f%%), ==0: %d/%d (%.1f%%)\n",
            label, mw, "1e-5*max",
            n_below, n_total, 100.0*n_below/n_total,
            n_zero, n_total, 100.0*n_zero/n_total)
end

# ── 6. Weight correlation across solvers ─────────────────────────
println("\n" * "=" ^ 90)
println("  6. WEIGHT CORRELATION (cold-start)")
println("=" ^ 90)

wh = Float64.(w_hc[1:sa.n])
wm1 = Float64.(w_mm1[1:sa.n])
wm5 = Float64.(w_mm5[1:sa.n])

@printf("  Huber vs MM(K=1):  cor = %.6f\n", cor(wh, wm1))
@printf("  Huber vs MM(K=5):  cor = %.6f\n", cor(wh, wm5))
@printf("  MM(K=1) vs MM(K=5): cor = %.6f\n", cor(wm1, wm5))

# ── 7. Convergence trace for MM K=5 ──────────────────────────────
println("\n" * "=" ^ 90)
println("  7. MM(K=5) CONVERGENCE TRACE (cold-start, first 30 iters + last 5)")
println("=" ^ 90)

@printf("  %5s  %16s  %12s\n", "Iter", "Poisson LL", "max_rel_diff")
println("  " * "─" ^ 40)
n_trace = length(diffs_mm5)
show_iters = sort(unique(vcat(
    collect(1:min(30, n_trace)),
    collect(max(1, n_trace-4):n_trace)
)))
for i in show_iters
    i > n_trace && continue
    @printf("  %5d  %16.6e  %12.4e\n", i, lls_mm5[i+1], diffs_mm5[i])
end

# LL monotonicity check (on scaled LL values)
ll_mono = all(lls_mm5[i+1] >= lls_mm5[i] - 1.0 for i in 1:length(lls_mm5)-1)
@printf("\n  LL monotonic (allowing 1.0 tolerance for Float32 noise): %s\n", ll_mono)

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
