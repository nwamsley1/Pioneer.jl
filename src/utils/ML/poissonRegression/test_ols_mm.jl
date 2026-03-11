## Test OLS v1/v2 and Poisson MM v1/v2 vs Huber on real scan342335 data
## OLS v1  = original solveOLS!  (incremental residuals only)
## OLS v2  = solveOLS_v2!        (periodic reinitResiduals! every 10 sweeps)
## PMM v1  = solvePoissonMM!     (original)
## PMM v2  = solvePoissonMM_v2!  (periodic initMu! every 10 sweeps)
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/test_ols_mm.jl

using Printf, Statistics, Serialization, BenchmarkTools, Profile, PProf

# Need Pioneer's SparseArray type for deserialization
using Pioneer

# Local SparseArray + solvers shadow Pioneer's types
include("SparseArray.jl")
include("spectralLinearRegression_reference.jl")
include("spectralPoissonRegression.jl")

# ── Load real data ────────────────────────────────────────────────
data = deserialize("/Users/nathanwamsley/Desktop/solveHuber_inputs_scan342335.jls")

sa = Main.SparseArray(
    data[:n_vals], data[:m], data[:n],
    Vector{Int64}(data[:rowval]), Vector{UInt16}(data[:colval]), data[:nzval],
    data[:matched], data[:isotope], data[:x], Vector{Int64}(data[:colptr])
)

max_outer     = data[:max_iter_outer]
max_outer_i64 = Int64(max_outer)
rel_conv      = data[:max_diff]

println("=" ^ 90)
println("  OLS v1/v2 + Poisson MM v1/v2 + Huber Test: scan342335")
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

function sse_from_residuals(r::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += Float64(r[i])^2; end; s
end

# ── OLS v1 wrapper ───────────────────────────────────────────────

function solve_ols_v1(sa_l, init_weights; ols_max_outer=max_outer, ols_rel_conv=rel_conv)
    r = zeros(Float32, max(sa_l.m, length(data[:residuals])))
    w = copy(init_weights)
    colnorm2 = Vector{Float32}(undef, sa_l.n)
    initResiduals_plain!(r, sa_l, w)
    solveOLS!(sa_l, r, w, colnorm2, Int64(ols_max_outer), ols_rel_conv)
    return w, r
end

# ── OLS v2 wrapper ───────────────────────────────────────────────

function solve_ols_v2(sa_l, init_weights; ols_max_outer=max_outer, ols_rel_conv=rel_conv)
    r = zeros(Float32, max(sa_l.m, length(data[:residuals])))
    w = copy(init_weights)
    colnorm2 = Vector{Float32}(undef, sa_l.n)
    initResiduals_plain!(r, sa_l, w)
    solveOLS_v2!(sa_l, r, w, colnorm2, Int64(ols_max_outer), ols_rel_conv)
    return w, r
end

# ── Poisson MM v1 wrapper ────────────────────────────────────────

function solve_pmm_v1(sa_l, init_weights; pmm_max_outer=max_outer, pmm_rel_conv=rel_conv)
    m = sa_l.m
    μ = zeros(Float32, m)
    y = zeros(Float32, m)
    w = copy(init_weights)
    initObserved!(y, sa_l)
    initMu!(μ, sa_l, w)
    solvePoissonMM!(sa_l, μ, y, w, Int64(pmm_max_outer), pmm_rel_conv)
    return w, μ, y
end

# ── Poisson MM v2 wrapper ────────────────────────────────────────

function solve_pmm_v2(sa_l, init_weights; pmm_max_outer=max_outer, pmm_rel_conv=rel_conv)
    m = sa_l.m
    μ = zeros(Float32, m)
    y = zeros(Float32, m)
    w = copy(init_weights)
    initObserved!(y, sa_l)
    initMu!(μ, sa_l, w)
    solvePoissonMM_v2!(sa_l, μ, y, w, Int64(pmm_max_outer), pmm_rel_conv)
    return w, μ, y
end

# ── Instrumented OLS v1 ─────────────────────────────────────────

function solve_ols_v1_counted(sa_l, init_weights; ols_max_outer=max_outer, ols_rel_conv=rel_conv)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, max(m, length(data[:residuals])))
    w = copy(init_weights)
    initResiduals_plain!(r, sa_l, w)

    colnorm2 = Vector{Float32}(undef, n)
    @inbounds for col in 1:n
        s = 0f0
        for i in sa_l.colptr[col]:(sa_l.colptr[col+1]-1)
            s += sa_l.nzval[i]^2
        end
        colnorm2[col] = s
    end

    max_weight = 0f0
    iters = 0
    sses = Float64[sse_from_residuals(r, m)]

    for iter in 1:ols_max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:n
            L2 = colnorm2[col]; iszero(L2) && continue
            L1 = 0f0
            @inbounds @fastmath for k in sa_l.colptr[col]:(sa_l.colptr[col+1]-1)
                L1 += sa_l.nzval[k] * r[sa_l.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1 / L2, 0f0)
            updateResiduals!(sa_l, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor; rc = δx / abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        push!(sses, sse_from_residuals(r, m))
        _diff < ols_rel_conv && break
    end
    return w, iters, r, sses
end

# ── Instrumented OLS v2 ─────────────────────────────────────────

function solve_ols_v2_counted(sa_l, init_weights; ols_max_outer=max_outer, ols_rel_conv=rel_conv, reinit_period=10)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, max(m, length(data[:residuals])))
    w = copy(init_weights)
    initResiduals_plain!(r, sa_l, w)

    colnorm2 = Vector{Float32}(undef, n)
    @inbounds for col in 1:n
        s = 0f0
        for i in sa_l.colptr[col]:(sa_l.colptr[col+1]-1)
            s += sa_l.nzval[i]^2
        end
        colnorm2[col] = s
    end

    max_weight = 0f0
    iters = 0
    sses = Float64[sse_from_residuals(r, m)]

    for iter in 1:ols_max_outer
        if iter > 1 && mod(iter - 1, reinit_period) == 0
            reinitResiduals!(r, sa_l, w)
        end
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:n
            L2 = colnorm2[col]; iszero(L2) && continue
            L1 = 0f0
            @inbounds @fastmath for k in sa_l.colptr[col]:(sa_l.colptr[col+1]-1)
                L1 += sa_l.nzval[k] * r[sa_l.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1 / L2, 0f0)
            updateResiduals!(sa_l, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor; rc = δx / abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        push!(sses, sse_from_residuals(r, m))
        _diff < ols_rel_conv && break
    end
    return w, iters, r, sses
end

# ── Instrumented Poisson MM v1 ───────────────────────────────────

function solve_pmm_v1_counted(sa_l, init_weights; pmm_max_outer=max_outer, pmm_rel_conv=rel_conv, max_inner=5)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); y = zeros(Float32, m)
    w = copy(init_weights)
    initObserved!(y, sa_l)
    initMu!(μ, sa_l, w)

    # Y-scaling (mirrors solvePoissonMM!)
    y_scale = maximum(y[1:m])
    if y_scale > 1f0
        @inbounds for i in 1:m; y[i] /= y_scale; end
        @inbounds for j in 1:n; w[j] /= y_scale; end
        initMu!(μ, sa_l, w)
    end

    ε = Float32(POISSON_MU_FLOOR)
    inner_tol = pmm_rel_conv
    max_weight = 0f0
    iters = 0
    lls = Float64[poissonLogLikelihood(μ, y, m)]

    for iter in 1:pmm_max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-4) : 0f0
        max_weight = 0f0
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, y, col)
                (L2 <= ε || isnan(L1)) && break
                X0 = w[col]; w[col] = max(w[col] - L1/L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                (iszero(w[col]) || (!iszero(X0) && abs_step/abs(X0) < inner_tol)) && break
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        push!(lls, poissonLogLikelihood(μ, y, m))
        _diff < pmm_rel_conv && break
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:m; y[i] *= y_scale; end
        @inbounds for j in 1:n; w[j] *= y_scale; end
        initMu!(μ, sa_l, w)
    end
    return w, iters, μ, y, lls
end

# ── Instrumented Poisson MM v2 ───────────────────────────────────

function solve_pmm_v2_counted(sa_l, init_weights; pmm_max_outer=max_outer, pmm_rel_conv=rel_conv, max_inner=5, reinit_period=10)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); y = zeros(Float32, m)
    w = copy(init_weights)
    initObserved!(y, sa_l)
    initMu!(μ, sa_l, w)

    # Y-scaling
    y_scale = maximum(y[1:m])
    if y_scale > 1f0
        @inbounds for i in 1:m; y[i] /= y_scale; end
        @inbounds for j in 1:n; w[j] /= y_scale; end
        initMu!(μ, sa_l, w)
    end

    ε = Float32(POISSON_MU_FLOOR)
    inner_tol = pmm_rel_conv
    max_weight = 0f0
    iters = 0
    lls = Float64[poissonLogLikelihood(μ, y, m)]

    for iter in 1:pmm_max_outer
        # Periodic μ recomputation
        if iter > 1 && mod(iter - 1, reinit_period) == 0
            initMu!(μ, sa_l, w)
        end

        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-4) : 0f0
        max_weight = 0f0
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, y, col)
                (L2 <= ε || isnan(L1)) && break
                X0 = w[col]; w[col] = max(w[col] - L1/L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
                abs_step = abs(w[col] - X0)
                (iszero(w[col]) || (!iszero(X0) && abs_step/abs(X0) < inner_tol)) && break
            end
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        push!(lls, poissonLogLikelihood(μ, y, m))
        _diff < pmm_rel_conv && break
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:m; y[i] *= y_scale; end
        @inbounds for j in 1:n; w[j] *= y_scale; end
        initMu!(μ, sa_l, w)
    end
    return w, iters, μ, y, lls
end

# ── Huber solver wrapper ─────────────────────────────────────────

function solve_huber(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, max(m, length(data[:residuals])))
    w = copy(init_weights)
    initResiduals_plain!(r, sa_l, w)
    reg = NoNorm()
    δ_param = data[:delta]; λ_param = data[:lambda]
    nr_iter = data[:max_iter_newton]; bs_iter = data[:max_iter_bisection]
    nr_acc  = data[:accuracy_newton]; bs_acc  = data[:accuracy_bisection]

    max_weight = 0f0; iters = 0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:n
            δx = abs(newton_bisection!(sa_l, r, w, col, δ_param, λ_param,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor; rc = δx / abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        _diff < rel_conv && break
    end
    return w, iters, r
end

# ══════════════════════════════════════════════════════════════════
# Run tests
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))
warm_w = copy(data[:weights])

# ── PMM optimized wrapper ──────────────────────────────────────
function solve_pmm_opt(sa_l, init_weights; pmm_max_outer=max_outer, pmm_rel_conv=rel_conv)
    m = sa_l.m
    μ = zeros(Float32, m)
    y = zeros(Float32, m)
    w = copy(init_weights)
    initObserved!(y, sa_l)
    initMu!(μ, sa_l, w)
    solvePoissonMM_opt!(sa_l, μ, y, w, Int64(pmm_max_outer), pmm_rel_conv)
    return w, μ, y
end

# ── Huber optimized wrapper ────────────────────────────────────
function solve_huber_opt(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    r = zeros(Float32, max(m, length(data[:residuals])))
    w = copy(init_weights)
    initResiduals_plain!(r, sa_l, w)
    reg = NoNorm()
    δ_param = data[:delta]; λ_param = data[:lambda]
    nr_iter = data[:max_iter_newton]; bs_iter = data[:max_iter_bisection]
    nr_acc  = data[:accuracy_newton]; bs_acc  = data[:accuracy_bisection]

    max_weight = 0f0; iters = 0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:n
            δx = abs(newton_bisection_opt!(sa_l, r, w, col, δ_param, λ_param,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor; rc = δx / abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        _diff < rel_conv && break
    end
    return w, iters, r
end

# ── 1. Warmup ALL solvers ────────────────────────────────────────
println("\n  Warming up all solvers...")
solve_ols_v1(sa, cold_w); solve_ols_v2(sa, cold_w)
solve_pmm_v1(sa, cold_w); solve_pmm_v2(sa, cold_w)
solve_pmm_opt(sa, cold_w)
solve_huber(sa, cold_w)
solve_huber_opt(sa, cold_w)
solve_ols_v1_counted(sa, cold_w); solve_ols_v2_counted(sa, cold_w)
solve_pmm_v1_counted(sa, cold_w); solve_pmm_v2_counted(sa, cold_w)
println("  Done.")

# ── 0. Correctness check: optimized vs original ─────────────────
println("\n" * "=" ^ 90)
println("  0. CORRECTNESS CHECK (optimized vs original)")
println("=" ^ 90)

w_pmm_orig, _, _ = solve_pmm_v1(sa, cold_w)
w_pmm_optv, _, _ = solve_pmm_opt(sa, cold_w)
pmm_maxdiff = maximum(abs.(w_pmm_orig[1:sa.n] .- w_pmm_optv[1:sa.n]))
pmm_cor = cor(Float64.(w_pmm_orig[1:sa.n]), Float64.(w_pmm_optv[1:sa.n]))
@printf("  PMM orig vs opt:   max|Δw| = %.4e   cor = %.10f\n", pmm_maxdiff, pmm_cor)

w_hub_orig, _, _ = solve_huber(sa, cold_w)
w_hub_optv, _, _ = solve_huber_opt(sa, cold_w)
hub_maxdiff = maximum(abs.(w_hub_orig[1:sa.n] .- w_hub_optv[1:sa.n]))
hub_cor = cor(Float64.(w_hub_orig[1:sa.n]), Float64.(w_hub_optv[1:sa.n]))
@printf("  Huber orig vs opt: max|Δw| = %.4e   cor = %.10f\n", hub_maxdiff, hub_cor)

# ── 2. Cold-start comparison ─────────────────────────────────────
println("\n" * "=" ^ 90)
println("  1. COLD-START COMPARISON")
println("=" ^ 90)

t_ov1_c = @elapsed w_ov1_c, it_ov1_c, r_ov1_c, sses_ov1_c = solve_ols_v1_counted(sa, cold_w)
sse_ov1_c = sse_from_residuals(r_ov1_c, sa.m)

t_ov2_c = @elapsed w_ov2_c, it_ov2_c, r_ov2_c, sses_ov2_c = solve_ols_v2_counted(sa, cold_w)
sse_ov2_c = sse_from_residuals(r_ov2_c, sa.m)

t_pv1_c = @elapsed w_pv1_c, it_pv1_c, μ_pv1_c, y_pv1_c, lls_pv1_c = solve_pmm_v1_counted(sa, cold_w)
t_pv2_c = @elapsed w_pv2_c, it_pv2_c, μ_pv2_c, y_pv2_c, lls_pv2_c = solve_pmm_v2_counted(sa, cold_w)

t_hub_c = @elapsed w_hub_c, it_hub_c, r_hub_c = solve_huber(sa, cold_w)
sse_hub_c = sse_from_residuals(r_hub_c, sa.m)

@printf("\n  %-20s  %5s  %10s  %16s\n", "Solver", "Iters", "Time(s)", "SSE / LL")
println("  " * "─" ^ 60)
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "OLS v1 (cold)", it_ov1_c, t_ov1_c, sse_ov1_c)
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "OLS v2 (cold)", it_ov2_c, t_ov2_c, sse_ov2_c)
@printf("  %-20s  %5d  %10.4f  %16.6e  (LL)\n",  "PMM v1 (cold)", it_pv1_c, t_pv1_c, lls_pv1_c[end])
@printf("  %-20s  %5d  %10.4f  %16.6e  (LL)\n",  "PMM v2 (cold)", it_pv2_c, t_pv2_c, lls_pv2_c[end])
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "Huber (cold)", it_hub_c, t_hub_c, sse_hub_c)

# ── 3. Warm-start comparison ─────────────────────────────────────
println("\n" * "=" ^ 90)
println("  2. WARM-START COMPARISON")
println("=" ^ 90)

t_ov1_w = @elapsed w_ov1_w, it_ov1_w, r_ov1_w, sses_ov1_w = solve_ols_v1_counted(sa, warm_w)
sse_ov1_w = sse_from_residuals(r_ov1_w, sa.m)

t_ov2_w = @elapsed w_ov2_w, it_ov2_w, r_ov2_w, sses_ov2_w = solve_ols_v2_counted(sa, warm_w)
sse_ov2_w = sse_from_residuals(r_ov2_w, sa.m)

t_pv1_w = @elapsed w_pv1_w, it_pv1_w, μ_pv1_w, y_pv1_w, lls_pv1_w = solve_pmm_v1_counted(sa, warm_w)
t_pv2_w = @elapsed w_pv2_w, it_pv2_w, μ_pv2_w, y_pv2_w, lls_pv2_w = solve_pmm_v2_counted(sa, warm_w)

t_hub_w = @elapsed w_hub_w, it_hub_w, r_hub_w = solve_huber(sa, warm_w)
sse_hub_w = sse_from_residuals(r_hub_w, sa.m)

@printf("\n  %-20s  %5s  %10s  %16s\n", "Solver", "Iters", "Time(s)", "SSE / LL")
println("  " * "─" ^ 60)
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "OLS v1 (warm)", it_ov1_w, t_ov1_w, sse_ov1_w)
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "OLS v2 (warm)", it_ov2_w, t_ov2_w, sse_ov2_w)
@printf("  %-20s  %5d  %10.4f  %16.6e  (LL)\n",  "PMM v1 (warm)", it_pv1_w, t_pv1_w, lls_pv1_w[end])
@printf("  %-20s  %5d  %10.4f  %16.6e  (LL)\n",  "PMM v2 (warm)", it_pv2_w, t_pv2_w, lls_pv2_w[end])
@printf("  %-20s  %5d  %10.4f  %16.6e  (SSE)\n", "Huber (warm)", it_hub_w, t_hub_w, sse_hub_w)

# ── 4. Weight correlation ────────────────────────────────────────
println("\n" * "=" ^ 90)
println("  3. WEIGHT CORRELATION")
println("=" ^ 90)

for (label, wov1, wov2, wpv1, wpv2, whub) in [
    ("Cold", w_ov1_c, w_ov2_c, w_pv1_c, w_pv2_c, w_hub_c),
    ("Warm", w_ov1_w, w_ov2_w, w_pv1_w, w_pv2_w, w_hub_w)]
    ov1 = Float64.(wov1[1:sa.n]); ov2 = Float64.(wov2[1:sa.n])
    pv1 = Float64.(wpv1[1:sa.n]); pv2 = Float64.(wpv2[1:sa.n])
    hub = Float64.(whub[1:sa.n])
    @printf("  %-6s  cor(OLSv1, Hub)=%.6f  cor(OLSv2, Hub)=%.6f  cor(OLSv1, OLSv2)=%.6f\n",
            label, cor(ov1, hub), cor(ov2, hub), cor(ov1, ov2))
    @printf("  %-6s  cor(PMMv1, Hub)=%.6f  cor(PMMv2, Hub)=%.6f  cor(PMMv1, PMMv2)=%.6f\n",
            "", cor(pv1, hub), cor(pv2, hub), cor(pv1, pv2))
    @printf("  %-6s  cor(OLSv1, PMMv1)=%.6f  max|OLSv1-PMMv1|=%.4e\n",
            "", cor(ov1, pv1), maximum(abs.(ov1 .- pv1)))
    @printf("  %-6s  max|v1-v2| OLS=%.4e  PMM=%.4e\n",
            "", maximum(abs.(ov1 .- ov2)), maximum(abs.(pv1 .- pv2)))
end

# ── 5. SSE / LL monotonicity ─────────────────────────────────────
println("\n" * "=" ^ 90)
println("  4. OBJECTIVE MONOTONICITY (cold-start)")
println("=" ^ 90)

# OLS: SSE should decrease
for (label, sses) in [("OLS v1", sses_ov1_c), ("OLS v2", sses_ov2_c)]
    violations = [i for i in 1:length(sses)-1 if sses[i+1] > sses[i] + 1.0]
    mono = isempty(violations)
    @printf("  %-10s  SSE monotonic (tol=1.0): %s", label, mono)
    !mono && @printf("  violations at iters: %s", join(violations, ", "))
    println()
end

# Poisson MM: LL should increase (we negate: -LL should decrease)
for (label, lls) in [("PMM v1", lls_pv1_c), ("PMM v2", lls_pv2_c)]
    violations = [i for i in 1:length(lls)-1 if lls[i+1] < lls[i] - 1.0]
    mono = isempty(violations)
    @printf("  %-10s  LL monotonic  (tol=1.0): %s", label, mono)
    !mono && @printf("  violations at iters: %s", join(violations, ", "))
    println()
end

# ── 6. Convergence traces ────────────────────────────────────────
println("\n" * "=" ^ 90)
println("  5. SSE CONVERGENCE TRACE — OLS (cold-start)")
println("=" ^ 90)

n_ov1 = length(sses_ov1_c) - 1; n_ov2 = length(sses_ov2_c) - 1
n_max_ols = max(n_ov1, n_ov2)
show_iters_ols = sort(unique(vcat(collect(1:min(20, n_max_ols)), collect(max(1, n_max_ols-4):n_max_ols))))

@printf("\n  %5s  %16s  %12s  │  %16s  %12s\n", "Iter", "SSE (v1)", "ΔSSE (v1)", "SSE (v2)", "ΔSSE (v2)")
println("  " * "─" ^ 78)
for i in show_iters_ols
    v1s = i <= n_ov1 ? @sprintf("%16.6e", sses_ov1_c[i+1]) : ""
    v1d = i <= n_ov1 ? @sprintf("%12.4e", sses_ov1_c[i+1] - sses_ov1_c[i]) : ""
    v2s = i <= n_ov2 ? @sprintf("%16.6e", sses_ov2_c[i+1]) : ""
    v2d = i <= n_ov2 ? @sprintf("%12.4e", sses_ov2_c[i+1] - sses_ov2_c[i]) : ""
    @printf("  %5d  %16s  %12s  │  %16s  %12s\n", i, v1s, v1d, v2s, v2d)
end

println("\n" * "=" ^ 90)
println("  5b. LL CONVERGENCE TRACE — Poisson MM (cold-start)")
println("=" ^ 90)

n_pv1 = length(lls_pv1_c) - 1; n_pv2 = length(lls_pv2_c) - 1
n_max_pmm = max(n_pv1, n_pv2)
show_iters_pmm = sort(unique(vcat(collect(1:min(20, n_max_pmm)), collect(max(1, n_max_pmm-4):n_max_pmm))))

@printf("\n  %5s  %16s  %12s  │  %16s  %12s\n", "Iter", "LL (v1)", "ΔLL (v1)", "LL (v2)", "ΔLL (v2)")
println("  " * "─" ^ 78)
for i in show_iters_pmm
    v1s = i <= n_pv1 ? @sprintf("%16.6e", lls_pv1_c[i+1]) : ""
    v1d = i <= n_pv1 ? @sprintf("%12.4e", lls_pv1_c[i+1] - lls_pv1_c[i]) : ""
    v2s = i <= n_pv2 ? @sprintf("%16.6e", lls_pv2_c[i+1]) : ""
    v2d = i <= n_pv2 ? @sprintf("%12.4e", lls_pv2_c[i+1] - lls_pv2_c[i]) : ""
    @printf("  %5d  %16s  %12s  │  %16s  %12s\n", i, v1s, v1d, v2s, v2d)
end

# ── 7. Weight sparsity ───────────────────────────────────────────
println("\n" * "=" ^ 90)
println("  6. WEIGHT SPARSITY")
println("=" ^ 90)

for (label, w_vec) in [("OLS v1 (cold)", w_ov1_c), ("OLS v2 (cold)", w_ov2_c),
                        ("PMM v1 (cold)", w_pv1_c), ("PMM v2 (cold)", w_pv2_c),
                        ("Huber (cold)", w_hub_c),
                        ("OLS v1 (warm)", w_ov1_w), ("OLS v2 (warm)", w_ov2_w),
                        ("PMM v1 (warm)", w_pv1_w), ("PMM v2 (warm)", w_pv2_w),
                        ("Huber (warm)", w_hub_w)]
    ww = w_vec[1:sa.n]
    mw = maximum(ww)
    threshold = mw * 1e-5
    n_below = count(ww .< threshold)
    n_zero = count(iszero.(ww))
    n_total = sa.n
    @printf("  %-18s  max_w=%.4e, <1e-5*max: %d/%d (%.1f%%), ==0: %d/%d (%.1f%%)\n",
            label, mw, n_below, n_total, 100.0*n_below/n_total,
            n_zero, n_total, 100.0*n_zero/n_total)
end

# ── 8. Timing sweep (BenchmarkTools) ─────────────────────────────
println("\n" * "=" ^ 90)
println("  7. TIMING SWEEP (BenchmarkTools @benchmark, evals=1)")
println("=" ^ 90)

# Pre-allocate buffers once (not timed)
_r   = Vector{Float32}(undef, max(sa.m, length(data[:residuals])))
_w   = Vector{Float32}(undef, length(data[:weights]))
_cn2 = Vector{Float32}(undef, sa.n)
_μ   = Vector{Float32}(undef, sa.m)
_y   = Vector{Float32}(undef, sa.m)

# Huber params at module scope for @benchmark
δ_bm     = data[:delta]
λ_bm     = data[:lambda]
nr_bm    = Int64(data[:max_iter_newton])
bs_bm    = Int64(data[:max_iter_bisection])
nr_acc_bm = data[:accuracy_newton]
bs_acc_bm = data[:accuracy_bisection]
reg_bm   = NoNorm()

for (label, init_w) in [("Cold", cold_w), ("Warm", warm_w)]
    println("\n  $label-start:")

    b_ov1 = @benchmark solveOLS!($sa, $_r, $_w, $_cn2, $max_outer_i64, $rel_conv) setup=begin
        copyto!($_w, $init_w); initResiduals_plain!($_r, $sa, $_w)
    end evals=1

    b_pv1 = @benchmark solvePoissonMM!($sa, $_μ, $_y, $_w, $max_outer_i64, $rel_conv) setup=begin
        copyto!($_w, $init_w); initObserved!($_y, $sa); initMu!($_μ, $sa, $_w)
    end evals=1

    b_pv1_opt = @benchmark solvePoissonMM_opt!($sa, $_μ, $_y, $_w, $max_outer_i64, $rel_conv) setup=begin
        copyto!($_w, $init_w); initObserved!($_y, $sa); initMu!($_μ, $sa, $_w)
    end evals=1

    b_hub = @benchmark solveHuber!($sa, $_r, $_w, $δ_bm, $λ_bm, $nr_bm, $bs_bm,
                                   $max_outer_i64, $nr_acc_bm, $bs_acc_bm,
                                   $rel_conv, $reg_bm) setup=begin
        copyto!($_w, $init_w); initResiduals_plain!($_r, $sa, $_w)
    end evals=1

    b_hub_opt = @benchmark solveHuber_opt!($sa, $_r, $_w, $δ_bm, $λ_bm, $nr_bm, $bs_bm,
                                            $max_outer_i64, $nr_acc_bm, $bs_acc_bm,
                                            $rel_conv, $reg_bm) setup=begin
        copyto!($_w, $init_w); initResiduals_plain!($_r, $sa, $_w)
    end evals=1

    # Extract median times (ns → seconds)
    t_ov1     = median(b_ov1).time / 1e9
    t_pv1     = median(b_pv1).time / 1e9
    t_pv1_opt = median(b_pv1_opt).time / 1e9
    t_hub     = median(b_hub).time / 1e9
    t_hub_opt = median(b_hub_opt).time / 1e9

    @printf("  %-6s  OLSv1: %.4f s  PMMv1: %.4f s  PMM_opt: %.4f s  Huber: %.4f s  Hub_opt: %.4f s\n",
            label, t_ov1, t_pv1, t_pv1_opt, t_hub, t_hub_opt)
    @printf("  %-6s  PMM optimization speedup: %.2fx  (%.1f%% faster)\n",
            "", t_pv1/t_pv1_opt, 100.0*(t_pv1-t_pv1_opt)/t_pv1)
    @printf("  %-6s  Huber optimization speedup: %.2fx  (%.1f%% faster)\n",
            "", t_hub/t_hub_opt, 100.0*(t_hub-t_hub_opt)/t_hub)
    @printf("  %-6s  vs Huber_opt → OLSv1: %.2fx  PMMv1_opt: %.2fx\n",
            "", t_hub_opt/t_ov1, t_hub_opt/t_pv1_opt)
    @printf("  %-6s  samples: OLSv1=%d  PMMv1=%d  PMM_opt=%d  Huber=%d  Hub_opt=%d\n",
            "", length(b_ov1.times), length(b_pv1.times),
            length(b_pv1_opt.times), length(b_hub.times), length(b_hub_opt.times))
end

# ── 9. CPU Profiling (flame graphs) ─────────────────────────────
println("\n" * "=" ^ 90)
println("  8. CPU PROFILING → flame graphs (1000 runs each, cold-start)")
println("=" ^ 90)

n_profile = 1000
profile_dir = joinpath(@__DIR__, "..", "..", "..", "..", "data")

# --- OLS v1 ---
println("  Profiling OLS v1...")
Profile.clear()
@profile for _ in 1:n_profile
    copyto!(_w, cold_w); initResiduals_plain!(_r, sa, _w)
    solveOLS!(sa, _r, _w, _cn2, max_outer_i64, rel_conv)
end
pprof(out=joinpath(profile_dir, "profile_ols_v1.pb.gz"), web=false)
println("    → wrote profile_ols_v1.pb.gz")

# --- PMM v1 ---
println("  Profiling PMM v1...")
Profile.clear()
@profile for _ in 1:n_profile
    copyto!(_w, cold_w); initObserved!(_y, sa); initMu!(_μ, sa, _w)
    solvePoissonMM!(sa, _μ, _y, _w, max_outer_i64, rel_conv)
end
pprof(out=joinpath(profile_dir, "profile_pmm_v1.pb.gz"), web=false)
println("    → wrote profile_pmm_v1.pb.gz")

# --- PMM v1 optimized ---
println("  Profiling PMM v1 optimized (1 div instead of 2)...")
Profile.clear()
@profile for _ in 1:n_profile
    copyto!(_w, cold_w); initObserved!(_y, sa); initMu!(_μ, sa, _w)
    solvePoissonMM_opt!(sa, _μ, _y, _w, max_outer_i64, rel_conv)
end
pprof(out=joinpath(profile_dir, "profile_pmm_v1_opt.pb.gz"), web=false)
println("    → wrote profile_pmm_v1_opt.pb.gz")

# --- Huber ---
println("  Profiling Huber...")
Profile.clear()
@profile for _ in 1:n_profile
    copyto!(_w, cold_w); initResiduals_plain!(_r, sa, _w)
    solveHuber!(sa, _r, _w, δ_bm, λ_bm, nr_bm, bs_bm,
                max_outer_i64, nr_acc_bm, bs_acc_bm, rel_conv, reg_bm)
end
pprof(out=joinpath(profile_dir, "profile_huber.pb.gz"), web=false)
println("    → wrote profile_huber.pb.gz")

# --- Huber optimized ---
println("  Profiling Huber optimized (x*x instead of x^2)...")
Profile.clear()
@profile for _ in 1:n_profile
    copyto!(_w, cold_w); initResiduals_plain!(_r, sa, _w)
    solveHuber_opt!(sa, _r, _w, δ_bm, λ_bm, nr_bm, bs_bm,
                    max_outer_i64, nr_acc_bm, bs_acc_bm, rel_conv, reg_bm)
end
pprof(out=joinpath(profile_dir, "profile_huber_opt.pb.gz"), web=false)
println("    → wrote profile_huber_opt.pb.gz")

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
