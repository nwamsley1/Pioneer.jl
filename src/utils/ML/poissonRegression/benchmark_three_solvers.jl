## Benchmark: PMM vs OLS vs Huber — Solution Quality + Timing
##
## Careful benchmarking protocol:
##   1. All buffers pre-allocated (zero allocation in timed region)
##   2. Full JIT warmup on first problem before any timing
##   3. Each problem timed N_RUNS times, median taken
##   4. State (weights, residuals, μ) reset before each timed run
##   5. GC disabled during timed region, forced between problems
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_three_solvers.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Helpers ──────────────────────────────────────────────────────

function load_sa(data)
    x_row = data[:x]
    rowval_i64 = Vector{Int64}(data[:rowval])
    if length(x_row) == data[:m]
        x_nz = Vector{Float32}(undef, data[:n_vals])
        for k in 1:data[:n_vals]; x_nz[k] = x_row[rowval_i64[k]]; end
    else
        x_nz = x_row
    end
    Main.SparseArray(
        data[:n_vals], data[:m], data[:n],
        rowval_i64, Vector{UInt16}(data[:colval]), data[:nzval],
        ones(Bool, data[:n_vals]), zeros(UInt8, data[:n_vals]), x_nz,
        Vector{Int64}(data[:colptr])
    )
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

# ── Timing harness ───────────────────────────────────────────────
#
# Each solver gets a "run!" function that:
#   (a) resets state from pre-allocated cold_w
#   (b) runs the solver
#   (c) returns nothing (so we don't time result construction)
#
# We time only the run! call, with GC disabled.

function run_pmm!(sa, μ, y, w, cold_w, max_outer, rel_conv, K)
    # Reset
    @inbounds for j in 1:sa.n; w[j] = cold_w[j]; end
    initObserved!(y, sa)
    initMu!(μ, sa, w)
    # Solve
    solvePoissonMM!(sa, μ, y, w, max_outer, rel_conv; max_inner_iter=Int64(K))
    return nothing
end

function run_ols!(sa, r, w, cold_w, colnorm2, max_outer, rel_conv)
    # Reset
    @inbounds for j in 1:sa.n; w[j] = cold_w[j]; end
    initResiduals_plain!(r, sa, w)
    # Solve
    solveOLS!(sa, r, w, colnorm2, max_outer, rel_conv)
    return nothing
end

function run_huber!(sa, r, w, cold_w, δ_hub, λ_hub, nr_max, bs_max, max_outer, nr_acc, bs_acc, rel_conv, reg)
    # Reset
    @inbounds for j in 1:sa.n; w[j] = cold_w[j]; end
    initResiduals_plain!(r, sa, w)
    # Solve
    solveHuber!(sa, r, w, δ_hub, λ_hub, nr_max, bs_max, max_outer, nr_acc, bs_acc, rel_conv, reg)
    return nothing
end

"""
    time_median(f, n_runs) → median_seconds

Time `f()` n_runs times with GC disabled, return median elapsed time.
"""
function time_median(f, n_runs)
    times = Vector{Float64}(undef, n_runs)
    for i in 1:n_runs
        GC.gc()  # clean slate before each run
        GC.enable(false)
        times[i] = @elapsed f()
        GC.enable(true)
    end
    return median(times)
end

# ── Config ───────────────────────────────────────────────────────

const N_RUNS = 11          # odd number → clean median
const K_PMM = 5
const PMM_REL_CONV = Float32(0.001)

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

# ── Warmup ───────────────────────────────────────────────────────
# Run each solver twice on the first problem to compile all code paths.

println("Warming up (JIT compilation)...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
cold_w0 = ones(Float32, sa0.n)

# PMM warmup
μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m); w0 = ones(Float32, sa0.n)
for _ in 1:2
    run_pmm!(sa0, μ0, y0, w0, cold_w0, Int64(d0[:max_iter_outer]), PMM_REL_CONV, K_PMM)
end

# OLS warmup
r0 = zeros(Float32, sa0.m); cn0 = zeros(Float32, sa0.n); w0 .= 1f0
for _ in 1:2
    run_ols!(sa0, r0, w0, cold_w0, cn0, Int64(d0[:max_iter_outer]), d0[:max_diff])
end

# Huber warmup
w0 .= 1f0
for _ in 1:2
    run_huber!(sa0, r0, w0, cold_w0, d0[:huber_delta], d0[:lambda],
               Int64(d0[:max_iter_newton]), Int64(d0[:max_iter_bisection]),
               Int64(d0[:max_iter_outer]), Float32(1e-6), Float32(1e-6),
               d0[:max_diff], NoNorm())
end

# Also warmup time_median itself
time_median(() -> nothing, 3)
println("  Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct ThreeSolverResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    n_vals::Int
    # Timing (median of N_RUNS)
    pmm_time::Float64
    ols_time::Float64
    hub_time::Float64
    # Solution quality
    pmm_ll::Float64
    ols_ll::Float64
    hub_ll::Float64
    pmm_sse::Float64
    ols_sse::Float64
    hub_sse::Float64
    # Agreement
    cor_pmm_ols::Float64
    cor_pmm_hub::Float64
    cor_ols_hub::Float64
end

results = ThreeSolverResult[]

# ── Run all problems ─────────────────────────────────────────────

println("=" ^ 150)
@printf("  %-8s %4s %5s %6s │ %10s %10s %10s │ %14s %14s %14s │ %8s %8s %8s\n",
        "Scan", "Cols", "Rows", "NNZ",
        "PMM(μs)", "OLS(μs)", "Hub(μs)",
        "PMM_LL", "OLS_LL", "Hub_LL",
        "cor(P,O)", "cor(P,H)", "cor(O,H)")
println("  " * "─" ^ 146)

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)

    max_outer  = Int64(data[:max_iter_outer])
    rel_conv   = data[:max_diff]
    δ_hub      = data[:huber_delta]
    λ_hub      = data[:lambda]
    nr_max     = Int64(data[:max_iter_newton])
    bs_max     = Int64(data[:max_iter_bisection])
    nr_acc     = haskey(data, :convergence_tol) ? data[:convergence_tol] : Float32(1e-6)
    bs_acc     = nr_acc
    reg        = NoNorm()
    cold_w     = ones(Float32, sa.n)

    # Pre-allocate all buffers (shared where possible, but separate for final solutions)
    μ_pmm  = zeros(Float32, sa.m)
    y_pmm  = zeros(Float32, sa.m)
    w_pmm  = ones(Float32, sa.n)

    r_ols  = zeros(Float32, sa.m)
    w_ols  = ones(Float32, sa.n)
    cn_ols = zeros(Float32, sa.n)

    r_hub  = zeros(Float32, sa.m)
    w_hub  = ones(Float32, sa.n)

    # ── Time PMM ──
    t_pmm = time_median(N_RUNS) do
        run_pmm!(sa, μ_pmm, y_pmm, w_pmm, cold_w, max_outer, PMM_REL_CONV, K_PMM)
    end

    # ── Time OLS ──
    t_ols = time_median(N_RUNS) do
        run_ols!(sa, r_ols, w_ols, cold_w, cn_ols, max_outer, rel_conv)
    end

    # ── Time Huber ──
    t_hub = time_median(N_RUNS) do
        run_huber!(sa, r_hub, w_hub, cold_w, δ_hub, λ_hub, nr_max, bs_max,
                   max_outer, nr_acc, bs_acc, rel_conv, reg)
    end

    # ── Solution quality (from last run, which is a valid converged solution) ──
    # PMM: LL from internal μ
    pmm_ll = poissonLogLikelihood(μ_pmm, y_pmm, sa.m)

    # OLS: build μ to compute LL, compute SSE from residuals
    μ_ols = zeros(Float32, sa.m)
    y_obs = zeros(Float32, sa.m)
    initMu!(μ_ols, sa, w_ols)
    initObserved!(y_obs, sa)
    ols_ll = poissonLogLikelihood(μ_ols, y_obs, sa.m)
    ols_sse = 0.0
    @inbounds for i in 1:sa.m; ols_sse += Float64(r_ols[i])^2; end

    # Huber: same
    μ_hub = zeros(Float32, sa.m)
    initMu!(μ_hub, sa, w_hub)
    hub_ll = poissonLogLikelihood(μ_hub, y_obs, sa.m)
    hub_sse = 0.0
    @inbounds for i in 1:sa.m; hub_sse += Float64(r_hub[i])^2; end

    # PMM SSE: build residuals
    r_pmm = zeros(Float32, sa.m)
    initResiduals_plain!(r_pmm, sa, w_pmm)
    pmm_sse = 0.0
    @inbounds for i in 1:sa.m; pmm_sse += Float64(r_pmm[i])^2; end

    # ── Correlations ──
    pv = Float64.(w_pmm[1:sa.n])
    ov = Float64.(w_ols[1:sa.n])
    hv = Float64.(w_hub[1:sa.n])

    safe_cor(a, b) = (all(iszero, a) || all(iszero, b)) ? NaN : cor(a, b)
    c_po = safe_cor(pv, ov)
    c_ph = safe_cor(pv, hv)
    c_oh = safe_cor(ov, hv)

    push!(results, ThreeSolverResult(
        data[:scan_idx], sa.n, sa.m, sa.n_vals,
        t_pmm, t_ols, t_hub,
        pmm_ll, ols_ll, hub_ll,
        pmm_sse, ols_sse, hub_sse,
        c_po, c_ph, c_oh
    ))

    @printf("  %-8d %4d %5d %6d │ %10.1f %10.1f %10.1f │ %14.4e %14.4e %14.4e │ %8.4f %8.4f %8.4f\n",
            data[:scan_idx], sa.n, sa.m, sa.n_vals,
            t_pmm * 1e6, t_ols * 1e6, t_hub * 1e6,
            pmm_ll, ols_ll, hub_ll,
            isnan(c_po) ? -1.0 : c_po,
            isnan(c_ph) ? -1.0 : c_ph,
            isnan(c_oh) ? -1.0 : c_oh)
end

# ── Summary ──────────────────────────────────────────────────────

println("\n" * "=" ^ 150)
println("  SUMMARY — $(length(results)) problems, PMM(K=$K_PMM, rel_conv=$PMM_REL_CONV), N_RUNS=$N_RUNS")
println("=" ^ 150)

# ── Timing ──
pmm_t = [r.pmm_time for r in results]
ols_t = [r.ols_time for r in results]
hub_t = [r.hub_time for r in results]

println("\n  TIMING (median of $N_RUNS runs per problem, in microseconds)")
@printf("  %-8s │ %10s %10s %10s %10s %10s\n",
        "Solver", "Mean", "Median", "P10", "P90", "Max")
println("  " * "─" ^ 65)
for (label, times) in [("PMM", pmm_t), ("OLS", ols_t), ("Huber", hub_t)]
    t_us = times .* 1e6
    @printf("  %-8s │ %10.1f %10.1f %10.1f %10.1f %10.1f\n",
            label, mean(t_us), median(t_us), quantile(t_us, 0.1), quantile(t_us, 0.9), maximum(t_us))
end

# Speedup ratios (per-problem)
println("\n  Speedup ratios (per-problem, then aggregated):")
pmm_vs_ols = pmm_t ./ max.(ols_t, 1e-12)
hub_vs_ols = hub_t ./ max.(ols_t, 1e-12)
pmm_vs_hub = pmm_t ./ max.(hub_t, 1e-12)
@printf("    PMM/OLS:   mean=%.2fx  median=%.2fx  (PMM is %.1f%% %s)\n",
        mean(pmm_vs_ols), median(pmm_vs_ols),
        abs(median(pmm_vs_ols) - 1.0) * 100,
        median(pmm_vs_ols) > 1 ? "slower" : "faster")
@printf("    Hub/OLS:   mean=%.2fx  median=%.2fx  (Huber is %.1f%% %s)\n",
        mean(hub_vs_ols), median(hub_vs_ols),
        abs(median(hub_vs_ols) - 1.0) * 100,
        median(hub_vs_ols) > 1 ? "slower" : "faster")
@printf("    PMM/Hub:   mean=%.2fx  median=%.2fx  (PMM is %.1f%% %s)\n",
        mean(pmm_vs_hub), median(pmm_vs_hub),
        abs(median(pmm_vs_hub) - 1.0) * 100,
        median(pmm_vs_hub) > 1 ? "slower" : "faster")

# Total wall time
@printf("\n  Total wall time (sum of medians): PMM=%.4fs  OLS=%.4fs  Huber=%.4fs\n",
        sum(pmm_t), sum(ols_t), sum(hub_t))

# ── Solution Quality: Poisson LL ──
println("\n  POISSON LOG-LIKELIHOOD (higher = better for Poisson objective)")
pmm_ll = [r.pmm_ll for r in results]
ols_ll = [r.ols_ll for r in results]
hub_ll = [r.hub_ll for r in results]

n_pmm_best_ll = count(i -> pmm_ll[i] >= ols_ll[i] && pmm_ll[i] >= hub_ll[i], 1:length(results))
n_ols_best_ll = count(i -> ols_ll[i] >= pmm_ll[i] && ols_ll[i] >= hub_ll[i], 1:length(results))
n_hub_best_ll = count(i -> hub_ll[i] >= pmm_ll[i] && hub_ll[i] >= ols_ll[i], 1:length(results))
@printf("    PMM has highest LL: %d/%d\n", n_pmm_best_ll, length(results))
@printf("    OLS has highest LL: %d/%d\n", n_ols_best_ll, length(results))
@printf("    Hub has highest LL: %d/%d\n", n_hub_best_ll, length(results))

# ── Solution Quality: SSE ──
println("\n  SUM OF SQUARED ERRORS (lower = better for SSE objective)")
pmm_sse = [r.pmm_sse for r in results]
ols_sse = [r.ols_sse for r in results]
hub_sse = [r.hub_sse for r in results]

@printf("  %-8s │ %14s %14s %14s\n", "Solver", "Mean", "Median", "Max")
println("  " * "─" ^ 55)
for (label, sses) in [("PMM", pmm_sse), ("OLS", ols_sse), ("Huber", hub_sse)]
    @printf("  %-8s │ %14.4e %14.4e %14.4e\n", label, mean(sses), median(sses), maximum(sses))
end

n_ols_best_sse = count(i -> ols_sse[i] <= min(pmm_sse[i], hub_sse[i]) + 1.0, 1:length(results))
n_pmm_best_sse = count(i -> pmm_sse[i] <= min(ols_sse[i], hub_sse[i]) + 1.0, 1:length(results))
@printf("\n    OLS has lowest SSE (within 1.0): %d/%d\n", n_ols_best_sse, length(results))
@printf("    PMM has lowest SSE (within 1.0): %d/%d\n", n_pmm_best_sse, length(results))

# ── Weight Correlations ──
println("\n  WEIGHT CORRELATIONS")
cors_po = filter(!isnan, [r.cor_pmm_ols for r in results])
cors_ph = filter(!isnan, [r.cor_pmm_hub for r in results])
cors_oh = filter(!isnan, [r.cor_ols_hub for r in results])

@printf("  %-12s │ %8s %8s %8s %8s\n", "Pair", "Mean", "Median", "Min", "Max")
println("  " * "─" ^ 52)
for (label, cors) in [("PMM vs OLS", cors_po), ("PMM vs Hub", cors_ph), ("OLS vs Hub", cors_oh)]
    @printf("  %-12s │ %8.4f %8.4f %8.4f %8.4f\n",
            label, mean(cors), median(cors), minimum(cors), maximum(cors))
end

# ── Breakdown by problem size ──
println("\n  TIMING BY PROBLEM SIZE (median μs)")
@printf("  %-14s %5s │ %10s %10s %10s │ %8s %8s\n",
        "Size bin", "Count", "PMM(μs)", "OLS(μs)", "Hub(μs)", "PMM/OLS", "Hub/OLS")
println("  " * "─" ^ 80)
for (lo, hi, label) in [(1, 10, "1-10 cols"), (11, 30, "11-30 cols"),
                          (31, 100, "31-100 cols"), (101, 500, "101-500 cols")]
    subset = filter(r -> lo <= r.n_cols <= hi, results)
    isempty(subset) && continue
    n = length(subset)
    mp = median([r.pmm_time for r in subset]) * 1e6
    mo = median([r.ols_time for r in subset]) * 1e6
    mh = median([r.hub_time for r in subset]) * 1e6
    @printf("  %-14s %5d │ %10.1f %10.1f %10.1f │ %8.2fx %8.2fx\n",
            label, n, mp, mo, mh, mp/max(mo, 0.01), mh/max(mo, 0.01))
end

println("\n" * "=" ^ 150)
println("  Done.")
println("=" ^ 150)
