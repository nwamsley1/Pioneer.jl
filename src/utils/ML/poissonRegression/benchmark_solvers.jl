## Benchmark OLS v1, PMM v1, and Huber across all real-world problems
## in /Users/n.t.wamsley/Desktop/solveHuber_problems/
##
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/benchmark_solvers.jl

using Printf, Statistics, Serialization

# Need Pioneer's types for deserialization
using Pioneer

# Local SparseArray + all solvers
include("SparseArray.jl")
include("spectralLinearRegression_reference.jl")
include("spectralPoissonRegression.jl")

# ── Helpers ──────────────────────────────────────────────────────

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

# ── Instrumented solvers (compiled functions, not top-level) ────

function solveOLS_instrumented!(sa, r, w, cn, max_outer, rel_conv)
    @inbounds for col in 1:sa.n
        s = 0f0
        for i in sa.colptr[col]:(sa.colptr[col+1]-1)
            s += sa.nzval[i]^2
        end
        cn[col] = s
    end
    max_weight = 0f0
    iters = 0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:sa.n
            L2 = cn[col]; iszero(L2) && continue
            L1 = 0f0
            @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1)
                L1 += sa.nzval[k] * r[sa.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1/L2, 0f0)
            updateResiduals!(sa, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    return iters
end

function solvePMM_instrumented!(sa, μ, y, w, max_outer, rel_conv)
    # Y-scaling
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε_pmm = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0
    iters = 0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-4) : 0f0
        max_weight = 0f0
        for col in 1:sa.n
            L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
            if L2 > ε_pmm && !isnan(L1)
                X0 = w[col]
                w[col] = max(w[col] - L1/L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                δx = abs(w[col] - X0)
                w[col] > max_weight && (max_weight = w[col])
                if w[col] > weight_floor
                    rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
                end
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return iters
end

function solveHuber_instrumented!(sa, r, w, δ_hub, λ_hub, nr_max, bs_max, max_outer, nr_acc, bs_acc, rel_conv, reg)
    max_weight = 0f0
    iters = 0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : 0f0
        max_weight = 0f0
        for col in 1:sa.n
            δx = abs(newton_bisection!(sa, r, w, col, δ_hub, λ_hub,
                        nr_max, bs_max, nr_acc, bs_acc, reg, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    return iters
end

# ── Load all problems ────────────────────────────────────────────

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems in $problem_dir\n")

# ── Warmup on first problem ──────────────────────────────────────

println("Warming up solvers...")
d0 = deserialize(joinpath(problem_dir, files[1]))
begin
    x_row0 = d0[:x]
    rowval0_i64 = Vector{Int64}(d0[:rowval])
    if length(x_row0) == d0[:m]
        x_nz0 = Vector{Float32}(undef, d0[:n_vals])
        for k in 1:d0[:n_vals]; x_nz0[k] = x_row0[rowval0_i64[k]]; end
    else
        x_nz0 = x_row0
    end
end
sa0 = Main.SparseArray(
    d0[:n_vals], d0[:m], d0[:n],
    rowval0_i64, Vector{UInt16}(d0[:colval]), d0[:nzval],
    ones(Bool, d0[:n_vals]), zeros(UInt8, d0[:n_vals]), x_nz0,
    Vector{Int64}(d0[:colptr])
)
w0 = ones(Float32, sa0.n)
r0 = zeros(Float32, sa0.m); cn0 = zeros(Float32, sa0.n)
μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)

initResiduals_plain!(r0, sa0, w0)
solveOLS!(sa0, r0, w0, cn0, Int64(d0[:max_iter_outer]), d0[:max_diff])

w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePoissonMM!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), d0[:max_diff])

w0 .= 1f0; initResiduals_plain!(r0, sa0, w0)
solveHuber!(sa0, r0, w0, d0[:huber_delta], d0[:lambda],
            Int64(d0[:max_iter_newton]), Int64(d0[:max_iter_bisection]),
            Int64(d0[:max_iter_outer]), Float32(1e-6), Float32(1e-6),
            d0[:max_diff], NoNorm())
println("  Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct ProblemResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    n_vals::Int
    # OLS
    ols_time::Float64
    ols_iters::Int
    ols_sse::Float64
    # PMM
    pmm_time::Float64
    pmm_iters::Int
    pmm_sse::Float64
    pmm_ll::Float64
    # Huber
    hub_time::Float64
    hub_iters::Int
    hub_sse::Float64
    # Agreement
    cor_ols_hub::Float64
    cor_pmm_hub::Float64
    cor_ols_pmm::Float64
    maxdiff_ols_hub::Float64
    maxdiff_pmm_hub::Float64
end

results = ProblemResult[]

# ── Run all problems ─────────────────────────────────────────────

println("=" ^ 130)
@printf("  %-8s %4s %5s │ %10s %5s %12s │ %10s %5s %12s %14s │ %10s %5s %12s │ %8s %8s\n",
        "Scan", "Cols", "Rows",
        "OLS(s)", "Itr", "SSE",
        "PMM(s)", "Itr", "SSE", "LL",
        "Hub(s)", "Itr", "SSE",
        "cor(O,H)", "cor(P,H)")
println("  " * "─" ^ 126)

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))

    begin
        x_row = data[:x]
        rowval_i64 = Vector{Int64}(data[:rowval])
        if length(x_row) == data[:m]
            x_nz = Vector{Float32}(undef, data[:n_vals])
            for k in 1:data[:n_vals]; x_nz[k] = x_row[rowval_i64[k]]; end
        else
            x_nz = x_row
        end
    end
    sa = Main.SparseArray(
        data[:n_vals], data[:m], data[:n],
        rowval_i64, Vector{UInt16}(data[:colval]), data[:nzval],
        ones(Bool, data[:n_vals]), zeros(UInt8, data[:n_vals]), x_nz,
        Vector{Int64}(data[:colptr])
    )

    max_outer_i64 = Int64(data[:max_iter_outer])
    rel_conv      = data[:max_diff]
    δ_hub         = data[:huber_delta]
    λ_hub         = data[:lambda]
    nr_max        = Int64(data[:max_iter_newton])
    bs_max        = Int64(data[:max_iter_bisection])
    nr_acc        = haskey(data, :convergence_tol) ? data[:convergence_tol] : Float32(1e-6)
    bs_acc        = nr_acc
    reg           = NoNorm()

    cold_w = ones(Float32, sa.n)

    # ── OLS v1 (instrumented) ──
    r_ols  = zeros(Float32, sa.m)
    w_ols  = copy(cold_w)
    cn_ols = zeros(Float32, sa.n)
    initResiduals_plain!(r_ols, sa, w_ols)

    ols_iters = 0
    ols_time = @elapsed begin
        ols_iters = solveOLS_instrumented!(sa, r_ols, w_ols, cn_ols, max_outer_i64, rel_conv)
    end
    ols_sse = sse_from_residuals(r_ols, sa.m)

    # ── PMM v1 (instrumented) ──
    μ_pmm = zeros(Float32, sa.m)
    y_pmm = zeros(Float32, sa.m)
    w_pmm = copy(cold_w)
    initObserved!(y_pmm, sa)
    initMu!(μ_pmm, sa, w_pmm)

    pmm_iters = 0
    pmm_time = @elapsed begin
        pmm_iters = solvePMM_instrumented!(sa, μ_pmm, y_pmm, w_pmm, max_outer_i64, rel_conv)
    end
    # Compute SSE for PMM by building residuals
    r_pmm = zeros(Float32, sa.m)
    initResiduals_plain!(r_pmm, sa, w_pmm)
    pmm_sse = sse_from_residuals(r_pmm, sa.m)
    pmm_ll  = poissonLogLikelihood(μ_pmm, y_pmm, sa.m)

    # ── Huber (instrumented) ──
    r_hub = zeros(Float32, sa.m)
    w_hub = copy(cold_w)
    initResiduals_plain!(r_hub, sa, w_hub)

    hub_iters = 0
    hub_time = @elapsed begin
        hub_iters = solveHuber_instrumented!(sa, r_hub, w_hub, δ_hub, λ_hub,
                        nr_max, bs_max, max_outer_i64, nr_acc, bs_acc, rel_conv, reg)
    end
    hub_sse = sse_from_residuals(r_hub, sa.m)

    # ── Agreement ──
    wn = sa.n
    ov = Float64.(w_ols[1:wn])
    pv = Float64.(w_pmm[1:wn])
    hv = Float64.(w_hub[1:wn])

    # Correlation needs variance — if all zeros, set to 1.0
    c_oh = (all(iszero, ov) || all(iszero, hv)) ? NaN : cor(ov, hv)
    c_ph = (all(iszero, pv) || all(iszero, hv)) ? NaN : cor(pv, hv)
    c_op = (all(iszero, ov) || all(iszero, pv)) ? NaN : cor(ov, pv)
    md_oh = maximum(abs.(ov .- hv))
    md_ph = maximum(abs.(pv .- hv))

    push!(results, ProblemResult(
        data[:scan_idx], sa.n, sa.m, sa.n_vals,
        ols_time, ols_iters, ols_sse,
        pmm_time, pmm_iters, pmm_sse, pmm_ll,
        hub_time, hub_iters, hub_sse,
        c_oh, c_ph, c_op, md_oh, md_ph
    ))

    @printf("  %-8d %4d %5d │ %10.6f %5d %12.4e │ %10.6f %5d %12.4e %14.4e │ %10.6f %5d %12.4e │ %8.4f %8.4f\n",
            data[:scan_idx], sa.n, sa.m,
            ols_time, ols_iters, ols_sse,
            pmm_time, pmm_iters, pmm_sse, pmm_ll,
            hub_time, hub_iters, hub_sse,
            isnan(c_oh) ? -1.0 : c_oh,
            isnan(c_ph) ? -1.0 : c_ph)
end

# ── Summary statistics ───────────────────────────────────────────

println("\n" * "=" ^ 130)
println("  SUMMARY ACROSS $(length(results)) PROBLEMS")
println("=" ^ 130)

# Timing
ols_times = [r.ols_time for r in results]
pmm_times = [r.pmm_time for r in results]
hub_times = [r.hub_time for r in results]

@printf("\n  %-20s  %12s  %12s  %12s  %12s\n", "Timing", "Mean(s)", "Median(s)", "Min(s)", "Max(s)")
println("  " * "─" ^ 72)
@printf("  %-20s  %12.6f  %12.6f  %12.6f  %12.6f\n", "OLS v1",  mean(ols_times), median(ols_times), minimum(ols_times), maximum(ols_times))
@printf("  %-20s  %12.6f  %12.6f  %12.6f  %12.6f\n", "PMM v1",  mean(pmm_times), median(pmm_times), minimum(pmm_times), maximum(pmm_times))
@printf("  %-20s  %12.6f  %12.6f  %12.6f  %12.6f\n", "Huber",   mean(hub_times), median(hub_times), minimum(hub_times), maximum(hub_times))

# Speedup ratios (per-problem, then aggregate)
hub_vs_ols = [r.hub_time / max(r.ols_time, 1e-9) for r in results]
hub_vs_pmm = [r.hub_time / max(r.pmm_time, 1e-9) for r in results]
@printf("\n  Huber/OLS speedup:  mean=%.2fx  median=%.2fx\n", mean(hub_vs_ols), median(hub_vs_ols))
@printf("  Huber/PMM speedup:  mean=%.2fx  median=%.2fx\n", mean(hub_vs_pmm), median(hub_vs_pmm))

# Total time
@printf("\n  Total wall time:  OLS=%.4f s  PMM=%.4f s  Huber=%.4f s\n",
        sum(ols_times), sum(pmm_times), sum(hub_times))

# Iterations
ols_iters = [r.ols_iters for r in results]
pmm_iters = [r.pmm_iters for r in results]
hub_iters = [r.hub_iters for r in results]

@printf("\n  %-20s  %8s  %8s  %8s  %8s\n", "Iterations", "Mean", "Median", "Min", "Max")
println("  " * "─" ^ 56)
@printf("  %-20s  %8.1f  %8d  %8d  %8d\n", "OLS v1",  mean(ols_iters), Int(median(ols_iters)), minimum(ols_iters), maximum(ols_iters))
@printf("  %-20s  %8.1f  %8d  %8d  %8d\n", "PMM v1",  mean(pmm_iters), Int(median(pmm_iters)), minimum(pmm_iters), maximum(pmm_iters))
@printf("  %-20s  %8.1f  %8d  %8d  %8d\n", "Huber",   mean(hub_iters), Int(median(hub_iters)), minimum(hub_iters), maximum(hub_iters))

# SSE comparison
ols_sses = [r.ols_sse for r in results]
pmm_sses = [r.pmm_sse for r in results]
hub_sses = [r.hub_sse for r in results]

@printf("\n  %-20s  %12s  %12s  %12s\n", "SSE", "Mean", "Median", "Max")
println("  " * "─" ^ 60)
@printf("  %-20s  %12.4e  %12.4e  %12.4e\n", "OLS v1",  mean(ols_sses), median(ols_sses), maximum(ols_sses))
@printf("  %-20s  %12.4e  %12.4e  %12.4e\n", "PMM v1",  mean(pmm_sses), median(pmm_sses), maximum(pmm_sses))
@printf("  %-20s  %12.4e  %12.4e  %12.4e\n", "Huber",   mean(hub_sses), median(hub_sses), maximum(hub_sses))

# Who wins on SSE?
ols_best = count(r -> r.ols_sse <= min(r.pmm_sse, r.hub_sse) + 1.0, results)
pmm_best = count(r -> r.pmm_sse <= min(r.ols_sse, r.hub_sse) + 1.0, results)
hub_best = count(r -> r.hub_sse <= min(r.ols_sse, r.pmm_sse) + 1.0, results)
@printf("\n  Lowest SSE (within 1.0):  OLS=%d  PMM=%d  Huber=%d  (of %d)\n",
        ols_best, pmm_best, hub_best, length(results))

# Agreement with Huber
cors_oh = filter(!isnan, [r.cor_ols_hub for r in results])
cors_ph = filter(!isnan, [r.cor_pmm_hub for r in results])
cors_op = filter(!isnan, [r.cor_ols_pmm for r in results])
mds_oh = [r.maxdiff_ols_hub for r in results]
mds_ph = [r.maxdiff_pmm_hub for r in results]

@printf("\n  %-20s  %10s  %10s  %10s  %12s  %12s\n",
        "Agreement vs Huber", "cor mean", "cor med", "cor min", "maxΔw mean", "maxΔw max")
println("  " * "─" ^ 76)
@printf("  %-20s  %10.6f  %10.6f  %10.6f  %12.4e  %12.4e\n",
        "OLS vs Huber", mean(cors_oh), median(cors_oh), minimum(cors_oh), mean(mds_oh), maximum(mds_oh))
@printf("  %-20s  %10.6f  %10.6f  %10.6f  %12.4e  %12.4e\n",
        "PMM vs Huber", mean(cors_ph), median(cors_ph), minimum(cors_ph), mean(mds_ph), maximum(mds_ph))
if !isempty(cors_op)
    @printf("  %-20s  %10.6f  %10.6f  %10.6f\n",
            "OLS vs PMM", mean(cors_op), median(cors_op), minimum(cors_op))
end

# Breakdown by problem size
println("\n  BREAKDOWN BY PROBLEM SIZE")
println("  " * "─" ^ 100)
@printf("  %-14s  %5s │ %12s %12s %12s │ %8s %8s %8s │ %8s %8s\n",
        "Size bin", "Count",
        "OLS med(s)", "PMM med(s)", "Hub med(s)",
        "OLS itr", "PMM itr", "Hub itr",
        "cor(O,H)", "cor(P,H)")
println("  " * "─" ^ 100)

for (lo, hi, label) in [(1, 10, "1-10 cols"), (11, 30, "11-30 cols"),
                          (31, 100, "31-100 cols"), (101, 500, "101-500 cols")]
    subset = filter(r -> lo <= r.n_cols <= hi, results)
    isempty(subset) && continue
    n = length(subset)
    @printf("  %-14s  %5d │ %12.6f %12.6f %12.6f │ %8.1f %8.1f %8.1f │ %8.4f %8.4f\n",
            label, n,
            median([r.ols_time for r in subset]),
            median([r.pmm_time for r in subset]),
            median([r.hub_time for r in subset]),
            mean([r.ols_iters for r in subset]),
            mean([r.pmm_iters for r in subset]),
            mean([r.hub_iters for r in subset]),
            mean(filter(!isnan, [r.cor_ols_hub for r in subset])),
            mean(filter(!isnan, [r.cor_pmm_hub for r in subset])))
end

println("\n" * "=" ^ 130)
println("  Done.")
println("=" ^ 130)
