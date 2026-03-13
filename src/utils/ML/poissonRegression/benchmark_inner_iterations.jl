## Benchmark: Inner Newton Iterations for PMM Solver
##
## Tests K=1,3,5,10 inner Newton steps per coordinate per outer sweep.
## For each K, reports: outer iterations, total Newton steps, timing,
## Poisson LL, and extended-run stability (EXT_UNSTABLE count).
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_inner_iterations.jl

using Printf, Statistics, Serialization

# Need Pioneer's types for deserialization
using Pioneer

# Local SparseArray + all solvers
include("SparseArray.jl")
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

# ── Instrumented PMM solver with inner iteration tracking ────────

"""
    solvePMM_inner_diag!(sa, μ, y, w, max_outer, rel_conv, max_inner)
        → (outer_iters, total_inner_steps, final_diff)

Like solvePoissonMM! but returns diagnostic counters:
- outer_iters: number of outer sweeps completed
- total_inner_steps: total Newton steps across all coordinates and sweeps
- final_diff: max relative weight change on the last sweep
"""
function solvePMM_inner_diag!(sa, μ, y, w, max_outer, rel_conv, max_inner)
    # Y-scaling
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0
    outer_iters = 0
    total_inner = 0
    final_diff = 0f0

    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0

        for col in 1:sa.n
            X_before = w[col]

            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                if L2 <= ε || isnan(L1)
                    break
                end
                total_inner += 1
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                if iszero(w[col]) || abs(w[col] - X0) / max(abs(w[col]), Float32(1e-10)) < Float32(1e-3)
                    break
                end
            end

            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / max(abs(w[col]), Float32(1e-10))
                rc > _diff && (_diff = rc)
            end
        end

        outer_iters = iter
        final_diff = _diff
        _diff < rel_conv && break
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return (outer_iters, total_inner, final_diff)
end

"""
    extended_pmm_inner_run!(sa, μ, y, w, extra_iters, rel_conv, max_inner)
        → max_rel_weight_change

Starting from a converged PMM solution, run `extra_iters` more outer sweeps
(each with `max_inner` inner Newton steps). Returns the maximum relative
weight change observed across all extra iterations.
"""
function extended_pmm_inner_run!(sa, μ, y, w, extra_iters, rel_conv, max_inner)
    # Y-scaling
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = maximum(w[1:sa.n])
    max_rel_change = 0f0

    for iter in 1:extra_iters
        weight_floor = max_weight * Float32(1e-2)
        max_weight = 0f0

        for col in 1:sa.n
            X_before = w[col]

            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                if L2 <= ε || isnan(L1)
                    break
                end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                if iszero(w[col]) || abs(w[col] - X0) / max(abs(w[col]), Float32(1e-10)) < Float32(1e-3)
                    break
                end
            end

            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / max(abs(w[col]), Float32(1e-10))
                rc > max_rel_change && (max_rel_change = rc)
            end
        end
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return max_rel_change
end

# ── Configuration ────────────────────────────────────────────────

const INNER_K_VALUES = [1, 3, 5, 10]
const PMM_REL_CONV   = Float32(0.001)
const EXTRA_ITERS    = 200
const EXT_UNSTABLE_THRESHOLD = 0.01

# ── Load problems ────────────────────────────────────────────────

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems in $problem_dir\n")

# ── Warmup ───────────────────────────────────────────────────────

println("Warming up solvers...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n)
μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)

for K in INNER_K_VALUES
    w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
    solvePMM_inner_diag!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), PMM_REL_CONV, K)
    extended_pmm_inner_run!(sa0, μ0, y0, w0, 10, PMM_REL_CONV, K)
end
println("  Done.\n")

# ── Per-K result storage ─────────────────────────────────────────

struct InnerResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    K::Int
    outer_iters::Int
    total_inner::Int
    time_s::Float64
    pmm_ll::Float64
    final_diff::Float64
    ext_max_change::Float64
    ext_unstable::Bool
end

all_results = Dict{Int, Vector{InnerResult}}()
for K in INNER_K_VALUES
    all_results[K] = InnerResult[]
end

# ── Run all problems × all K values ─────────────────────────────

println("=" ^ 170)
@printf("  %-8s %4s %5s", "Scan", "Cols", "Rows")
for K in INNER_K_VALUES
    @printf(" │ K=%-2d %4s %6s %10s %14s %10s %3s", K, "Oitr", "Nitr", "time(s)", "LL", "ext_Δw", "U")
end
println()
println("  " * "─" ^ 166)

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer_i64 = Int64(data[:max_iter_outer])

    @printf("  %-8d %4d %5d", data[:scan_idx], sa.n, sa.m)

    for K in INNER_K_VALUES
        # Fresh start
        w_k = ones(Float32, sa.n)
        μ_k = zeros(Float32, sa.m)
        y_k = zeros(Float32, sa.m)
        initObserved!(y_k, sa)
        initMu!(μ_k, sa, w_k)

        local oitr, nitr, fdiff
        t = @elapsed begin
            oitr, nitr, fdiff = solvePMM_inner_diag!(sa, μ_k, y_k, w_k,
                                                      max_outer_i64, PMM_REL_CONV, K)
        end

        ll = poissonLogLikelihood(μ_k, y_k, sa.m)

        # Extended stability test
        w_ext = copy(w_k)
        μ_ext = copy(μ_k)
        y_ext = copy(y_k)
        ext_change = extended_pmm_inner_run!(sa, μ_ext, y_ext, w_ext,
                                              EXTRA_ITERS, PMM_REL_CONV, K)
        unstable = ext_change > EXT_UNSTABLE_THRESHOLD

        push!(all_results[K], InnerResult(
            data[:scan_idx], sa.n, sa.m, K,
            oitr, nitr, t, ll, Float64(fdiff),
            Float64(ext_change), unstable
        ))

        @printf(" │      %4d %6d %10.6f %14.4e %10.2e %3s",
                oitr, nitr, t, ll, ext_change, unstable ? " U" : "  ")
    end
    println()
end

# ── Per-K Summary Tables ─────────────────────────────────────────

println("\n" * "=" ^ 170)
println("  INNER ITERATION BENCHMARK SUMMARY — $(length(files)) problems, rel_conv=$(PMM_REL_CONV)")
println("=" ^ 170)

println("\n  Per-K Summary:")
@printf("  %4s │ %8s %8s %8s │ %10s %10s %10s │ %10s %10s │ %14s %14s │ %4s/%d\n",
        "K", "Oitr_m", "Oitr_md", "Oitr_mx",
        "Nitr_m", "Nitr_md", "Nitr_mx",
        "Time_m(s)", "Time_md(s)",
        "LL_mean", "LL_median",
        "Unst", length(files))
println("  " * "─" ^ 130)

for K in INNER_K_VALUES
    res = all_results[K]
    oi = [r.outer_iters for r in res]
    ni = [r.total_inner for r in res]
    ts = [r.time_s for r in res]
    lls = [r.pmm_ll for r in res]
    n_unst = count(r -> r.ext_unstable, res)

    @printf("  %4d │ %8.1f %8d %8d │ %10.1f %10d %10d │ %10.6f %10.6f │ %14.4e %14.4e │ %4d/%d\n",
            K,
            mean(oi), Int(median(oi)), maximum(oi),
            mean(ni), Int(median(ni)), maximum(ni),
            mean(ts), median(ts),
            mean(lls), median(lls),
            n_unst, length(res))
end

# ── LL agreement across K values ─────────────────────────────────

println("\n  Poisson LL agreement across K values (per-problem max |LL_i - LL_j| / |LL_ref|):")
ll_diffs = Float64[]
for pi in 1:length(files)
    lls_k = [all_results[K][pi].pmm_ll for K in INNER_K_VALUES]
    ll_ref = maximum(abs.(lls_k))
    if ll_ref > 0
        push!(ll_diffs, (maximum(lls_k) - minimum(lls_k)) / ll_ref)
    end
end
if !isempty(ll_diffs)
    @printf("    mean=%.2e  median=%.2e  max=%.2e\n",
            mean(ll_diffs), median(ll_diffs), maximum(ll_diffs))
end

# ── Which K achieves highest LL most often? ──────────────────────

println("\n  Which K achieves highest LL (per problem)?")
k_wins = Dict(K => 0 for K in INNER_K_VALUES)
for pi in 1:length(files)
    best_ll = -Inf
    best_k = 0
    for K in INNER_K_VALUES
        ll = all_results[K][pi].pmm_ll
        if ll > best_ll
            best_ll = ll
            best_k = K
        end
    end
    k_wins[best_k] += 1
end
for K in INNER_K_VALUES
    @printf("    K=%2d:  %d/%d wins\n", K, k_wins[K], length(files))
end

# ── Detail for EXT_UNSTABLE cases per K ──────────────────────────

for K in INNER_K_VALUES
    unstable = filter(r -> r.ext_unstable, all_results[K])
    isempty(unstable) && continue
    println("\n  EXT_UNSTABLE problems for K=$K ($(length(unstable)) problems):")
    println("  " * "─" ^ 90)
    @printf("  %-8s %4s %5s │ %5s %6s %10s │ %14s │ %10s\n",
            "Scan", "Cols", "Rows", "Oitr", "Nitr", "time(s)", "LL", "ext_Δw")
    println("  " * "─" ^ 90)
    for r in unstable
        @printf("  %-8d %4d %5d │ %5d %6d %10.6f │ %14.4e │ %10.2e\n",
                r.scan_idx, r.n_cols, r.n_rows,
                r.outer_iters, r.total_inner, r.time_s,
                r.pmm_ll, r.ext_max_change)
    end
end

# ── Timing speedup from fewer outer iterations ──────────────────

println("\n  Timing relative to K=1:")
t1 = [r.time_s for r in all_results[1]]
for K in INNER_K_VALUES
    tk = [r.time_s for r in all_results[K]]
    ratios = tk ./ max.(t1, 1e-9)
    @printf("    K=%2d vs K=1:  mean=%.2fx  median=%.2fx  (total: %.4fs vs %.4fs)\n",
            K, mean(ratios), median(ratios), sum(tk), sum(t1))
end

println("\n" * "=" ^ 170)
println("  Done.")
println("=" ^ 170)
