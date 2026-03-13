## Benchmark: solvePoissonMM! (original) vs solvePoissonMM_fast! (optimized)
##
## Optimizations in _fast!:
##   1. Inlined + fused derivative/μ-update (saves 1 pass per inner iter 2..K)
##   2. Precomputed a_ij² (saves 1 multiply per nonzero in L2)
##   3. Column skipping (zero columns skipped after sweep 3)
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_pmm_optimization.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralLinearRegression_reference.jl")
include("spectralPoissonRegression.jl")

# ── Load problems ────────────────────────────────────────────────

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems in $problem_dir\n")

# ── Helper: build SparseArray from serialized data ───────────────

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

# ── Helper: run OLS for reference ────────────────────────────────

function initResiduals_plain!(r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds for i in 1:sa.m; r[i] = zero(T); end
    @inbounds for n in 1:sa.n_vals
        if iszero(r[sa.rowval[n]]); r[sa.rowval[n]] = -sa.x[n]; end
    end
    @inbounds for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1] - 1)
            r[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
end

# ── Warmup ───────────────────────────────────────────────────────

println("Warming up solvers...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)

# Warmup OLS
w0 = ones(Float32, sa0.n); r0 = zeros(Float32, sa0.m); cn0 = zeros(Float32, sa0.n)
initResiduals_plain!(r0, sa0, w0)
solveOLS!(sa0, r0, w0, cn0, Int64(d0[:max_iter_outer]), d0[:max_diff])

# Warmup original PMM
w0 .= 1f0; μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePoissonMM!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), d0[:max_diff])

# Warmup fast PMM
w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePoissonMM_fast!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), d0[:max_diff])

println("  Done.\n")

# ── Benchmark parameters ────────────────────────────────────────

const N_RUNS = 11  # runs per solver per problem, take median

# ── Results storage ──────────────────────────────────────────────

struct OptResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    n_vals::Int
    ols_time::Float64
    pmm_orig_time::Float64
    pmm_fast_time::Float64
    speedup::Float64          # orig / fast
    max_weight_diff::Float64  # max |w_orig - w_fast|
    ll_orig::Float64
    ll_fast::Float64
    ll_diff::Float64          # |ll_orig - ll_fast| / |ll_orig|
    cor_orig_fast::Float64
end

results = OptResult[]

# ── Run all problems ─────────────────────────────────────────────

println("=" ^ 120)
@printf("  %-8s %4s %5s %5s │ %10s │ %10s %10s %8s │ %12s %12s %10s\n",
        "Scan", "Cols", "Rows", "NNZ",
        "OLS(μs)",
        "Orig(μs)", "Fast(μs)", "Speedup",
        "maxΔw", "relΔLL", "cor(O,F)")
println("  " * "─" ^ 116)

GC.enable(false)

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])
    rel_conv  = data[:max_diff]

    # ── OLS reference (single run) ──
    w_ols = ones(Float32, sa.n)
    r_ols = zeros(Float32, sa.m); cn_ols = zeros(Float32, sa.n)
    initResiduals_plain!(r_ols, sa, w_ols)
    ols_time = @elapsed solveOLS!(sa, r_ols, w_ols, cn_ols, max_outer, rel_conv)

    # ── Original PMM: N_RUNS, take median ──
    orig_times = Vector{Float64}(undef, N_RUNS)
    w_orig = Vector{Float32}(undef, sa.n)
    μ_orig = zeros(Float32, sa.m)
    y_orig = zeros(Float32, sa.m)

    for run in 1:N_RUNS
        w_orig .= 1f0
        initObserved!(y_orig, sa)
        initMu!(μ_orig, sa, w_orig)
        orig_times[run] = @elapsed solvePoissonMM!(sa, μ_orig, y_orig, w_orig, max_outer, rel_conv)
    end
    # Last run's solution is in w_orig, μ_orig, y_orig
    ll_orig = poissonLogLikelihood(μ_orig, y_orig, sa.m)

    # ── Fast PMM: N_RUNS, take median ──
    fast_times = Vector{Float64}(undef, N_RUNS)
    w_fast = Vector{Float32}(undef, sa.n)
    μ_fast = zeros(Float32, sa.m)
    y_fast = zeros(Float32, sa.m)

    for run in 1:N_RUNS
        w_fast .= 1f0
        initObserved!(y_fast, sa)
        initMu!(μ_fast, sa, w_fast)
        fast_times[run] = @elapsed solvePoissonMM_fast!(sa, μ_fast, y_fast, w_fast, max_outer, rel_conv)
    end
    ll_fast = poissonLogLikelihood(μ_fast, y_fast, sa.m)

    # ── Compare solutions ──
    wn = sa.n
    ov = Float64.(w_orig[1:wn])
    fv = Float64.(w_fast[1:wn])
    max_wd = maximum(abs.(ov .- fv))
    ll_rel = abs(ll_orig) > 1e-10 ? abs(ll_orig - ll_fast) / abs(ll_orig) : 0.0
    c_of = (all(iszero, ov) || all(iszero, fv)) ? NaN : cor(ov, fv)

    orig_med = median(orig_times)
    fast_med = median(fast_times)
    spdup = orig_med / max(fast_med, 1e-12)

    push!(results, OptResult(
        data[:scan_idx], sa.n, sa.m, sa.n_vals,
        ols_time, orig_med, fast_med, spdup,
        max_wd, ll_orig, ll_fast, ll_rel, c_of
    ))

    @printf("  %-8d %4d %5d %5d │ %10.1f │ %10.1f %10.1f %7.2fx │ %12.4e %12.4e %10.6f\n",
            data[:scan_idx], sa.n, sa.m, sa.n_vals,
            ols_time * 1e6,
            orig_med * 1e6, fast_med * 1e6, spdup,
            max_wd, ll_rel,
            isnan(c_of) ? -1.0 : c_of)
end

GC.enable(true)

# ── Summary ──────────────────────────────────────────────────────

println("\n" * "=" ^ 120)
println("  SUMMARY ACROSS $(length(results)) PROBLEMS")
println("=" ^ 120)

ols_times  = [r.ols_time for r in results]
orig_times = [r.pmm_orig_time for r in results]
fast_times = [r.pmm_fast_time for r in results]
speedups   = [r.speedup for r in results]

@printf("\n  %-20s  %12s  %12s  %12s  %12s\n", "Timing (μs)", "Mean", "Median", "Min", "Max")
println("  " * "─" ^ 72)
@printf("  %-20s  %12.1f  %12.1f  %12.1f  %12.1f\n", "OLS",       mean(ols_times)*1e6, median(ols_times)*1e6, minimum(ols_times)*1e6, maximum(ols_times)*1e6)
@printf("  %-20s  %12.1f  %12.1f  %12.1f  %12.1f\n", "PMM original", mean(orig_times)*1e6, median(orig_times)*1e6, minimum(orig_times)*1e6, maximum(orig_times)*1e6)
@printf("  %-20s  %12.1f  %12.1f  %12.1f  %12.1f\n", "PMM fast",     mean(fast_times)*1e6, median(fast_times)*1e6, minimum(fast_times)*1e6, maximum(fast_times)*1e6)

@printf("\n  Speedup (orig/fast):  mean=%.2fx  median=%.2fx  min=%.2fx  max=%.2fx\n",
        mean(speedups), median(speedups), minimum(speedups), maximum(speedups))

# PMM/OLS ratio
orig_vs_ols = [r.pmm_orig_time / max(r.ols_time, 1e-12) for r in results]
fast_vs_ols = [r.pmm_fast_time / max(r.ols_time, 1e-12) for r in results]
@printf("\n  PMM/OLS ratio:  original=%.2fx  fast=%.2fx  (median)\n",
        median(orig_vs_ols), median(fast_vs_ols))

# Solution agreement
max_wds  = [r.max_weight_diff for r in results]
ll_diffs = [r.ll_diff for r in results]
cors     = filter(!isnan, [r.cor_orig_fast for r in results])

@printf("\n  Solution agreement (orig vs fast):\n")
@printf("    max|Δw|:   mean=%.4e  max=%.4e\n", mean(max_wds), maximum(max_wds))
@printf("    rel|ΔLL|:  mean=%.4e  max=%.4e\n", mean(ll_diffs), maximum(ll_diffs))
if !isempty(cors)
    @printf("    cor(w):    mean=%.8f  min=%.8f\n", mean(cors), minimum(cors))
end

# Total wall time
@printf("\n  Total wall time:  OLS=%.4f s  PMM orig=%.4f s  PMM fast=%.4f s\n",
        sum(ols_times), sum(orig_times), sum(fast_times))

println("\n" * "=" ^ 120)
println("  Done.")
println("=" ^ 120)
