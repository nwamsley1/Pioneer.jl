## How many iterations do the slow-converging problems actually need?
##
## For each of the ~24 unstable problems, run up to 2000 outer iterations
## and track when max_rel_change of weights >1% of max drops below 1%.
##
## Run:  julia --project=. src/utils/ML/poissonRegression/diagnose_slow_convergence.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")

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

# ── Instrumented solver: returns per-sweep convergence trace ─────

"""
    solvePMM_trace!(sa, μ, y, w, max_outer, max_inner)
        → (trace_all, trace_1pct, trace_ll)

Run PMM for up to max_outer sweeps (no early stopping).
Returns per-sweep:
- trace_all[i]: max relative weight change (all weights above floor)
- trace_1pct[i]: max relative weight change (weights >1% of max only)
- trace_ll[i]: Poisson LL after sweep i
"""
function solvePMM_trace!(sa, μ, y, w, max_outer, max_inner)
    # Y-scaling
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0

    trace_all = Float64[]
    trace_1pct = Float64[]
    trace_ll = Float64[]

    push!(trace_ll, poissonLogLikelihood(μ, y, sa.m))

    for iter in 1:max_outer
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0
        diff_all = 0f0
        diff_1pct = 0f0

        for col in 1:sa.n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                if iszero(w[col]) || abs(w[col] - X0) / max(abs(w[col]), Float32(1e-10)) < Float32(1e-3)
                    break
                end
            end

            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])

            if w[col] > weight_floor && !iszero(w[col])
                rc = δx / max(abs(w[col]), Float32(1e-10))
                rc > diff_all && (diff_all = rc)
            end
        end

        # Compute 1%-of-max metric after sweep (now max_weight is known)
        # We need to re-check — but we can't go back. Instead, approximate:
        # weight_floor_1pct = max_weight * 0.01
        # But we already tracked w[col] > weight_floor above, which uses the
        # previous sweep's max. For the 1% metric, let's just track it properly:
        # Actually, we can compute it from the trace by running a second pass...
        # For simplicity, we'll track it in a post-sweep check.

        push!(trace_all, Float64(diff_all))
        push!(trace_ll, poissonLogLikelihood(μ, y, sa.m))
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end

    return (trace_all, trace_ll)
end

"""
    solvePMM_trace_1pct!(sa, μ, y, w, max_outer, max_inner)

Like solvePMM_trace! but tracks convergence of weights >1% of max specifically.
Saves w snapshot each sweep to compute per-weight changes.
"""
function solvePMM_trace_1pct!(sa, μ, y, w, max_outer, max_inner)
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0

    trace_1pct = Float64[]
    trace_all = Float64[]
    trace_ll = Float64[]

    w_prev = copy(w[1:sa.n])

    for iter in 1:max_outer
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0

        for col in 1:sa.n
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                if iszero(w[col]) || abs(w[col] - X0) / max(abs(w[col]), Float32(1e-10)) < Float32(1e-3)
                    break
                end
            end
            w[col] > max_weight && (max_weight = w[col])
        end

        # Now max_weight is known for this sweep — compute both metrics
        diff_all = 0.0
        diff_1pct = 0.0
        floor_all = Float64(weight_floor)
        floor_1pct = Float64(max_weight) * 0.01

        for j in 1:sa.n
            wj = Float64(w[j])
            wp = Float64(w_prev[j])
            denom = max(abs(wj), 1e-10)
            rc = abs(wj - wp) / denom

            if wj > floor_all
                rc > diff_all && (diff_all = rc)
            end
            if wj > floor_1pct
                rc > diff_1pct && (diff_1pct = rc)
            end
        end

        push!(trace_all, diff_all)
        push!(trace_1pct, diff_1pct)
        push!(trace_ll, poissonLogLikelihood(μ, y, sa.m))

        # Save for next sweep
        @inbounds for j in 1:sa.n; w_prev[j] = w[j]; end
    end

    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end

    return (trace_all, trace_1pct, trace_ll)
end

# ── Config ───────────────────────────────────────────────────────

const K = 5
const MAX_OUTER = 2000
const STABLE_THRESHOLD = 0.01  # <1% = "stable"

# The 24 scan_idxs that were unstable at K=5, LL_tol=1e-6 in the previous benchmark
const UNSTABLE_SCANS = Set([
    106459, 106763, 108279, 108289, 108582, 108589, 88886, 90088,
    94945, 96150, 96155, 96452, 96453, 96754, 96757, 97064,
    97363, 97367, 97451, 97667, 97670, 97671, 97967, 97985
])

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems, $(length(UNSTABLE_SCANS)) targeted\n")

# Warmup
println("Warming up...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePMM_trace_1pct!(sa0, μ0, y0, w0, 50, K)
println("Done.\n")

# ── Run targeted problems ────────────────────────────────────────

println("=" ^ 140)
println("  SLOW CONVERGENCE ANALYSIS — K=$K, up to $MAX_OUTER outer iterations")
println("=" ^ 140)

struct SlowConvResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    # When does diff_all drop below 1%?
    iter_all_stable::Int     # -1 if never
    # When does diff_1pct drop below 1%?
    iter_1pct_stable::Int    # -1 if never
    # When does diff_1pct drop below 0.1%?
    iter_1pct_tight::Int     # -1 if never
    # Final values
    final_diff_all::Float64
    final_diff_1pct::Float64
    final_ll::Float64
    # Time
    time_s::Float64
    # Total Newton steps at stabilization points
    total_newton_at_1pct::Int
end

results = SlowConvResult[]

for fname in files
    data = deserialize(joinpath(problem_dir, fname))
    data[:scan_idx] in UNSTABLE_SCANS || continue

    sa = load_sa(data)

    w_k = ones(Float32, sa.n)
    μ_k = zeros(Float32, sa.m); y_k = zeros(Float32, sa.m)
    initObserved!(y_k, sa); initMu!(μ_k, sa, w_k)

    local trace_all, trace_1pct, trace_ll
    t = @elapsed begin
        trace_all, trace_1pct, trace_ll = solvePMM_trace_1pct!(sa, μ_k, y_k, w_k, MAX_OUTER, K)
    end

    # Find stabilization points
    iter_all_stable = -1
    iter_1pct_stable = -1
    iter_1pct_tight = -1

    for i in 1:length(trace_all)
        if iter_all_stable == -1 && trace_all[i] < STABLE_THRESHOLD
            iter_all_stable = i
        end
        if iter_1pct_stable == -1 && trace_1pct[i] < 0.01
            iter_1pct_stable = i
        end
        if iter_1pct_tight == -1 && trace_1pct[i] < 0.001
            iter_1pct_tight = i
        end
    end

    # Estimate total Newton steps at 1pct stabilization
    # Approximate: n_cols * iter * avg_inner_per_col (roughly 2-3 for K=5)
    # Actually just multiply: iters * n_cols (lower bound, 1 Newton per col per sweep)
    total_newton_est = iter_1pct_stable > 0 ? iter_1pct_stable * sa.n : -1

    push!(results, SlowConvResult(
        data[:scan_idx], sa.n, sa.m,
        iter_all_stable, iter_1pct_stable, iter_1pct_tight,
        trace_all[end], trace_1pct[end], trace_ll[end],
        t, total_newton_est
    ))

    # Print convergence trace (selected iterations)
    println("\n  ── Scan $(data[:scan_idx]) │ $(sa.n) cols × $(sa.m) rows │ $(Printf.@sprintf("%.3f", t))s ──")
    @printf("  %6s  %10s  %10s  %14s\n", "Iter", "diff_all", "diff_>1%%", "LL")
    println("  " * "─" ^ 50)

    show_iters = sort(unique(vcat(
        collect(1:min(10, MAX_OUTER)),
        collect(15:5:min(50, MAX_OUTER)),
        collect(50:50:min(500, MAX_OUTER)),
        collect(500:250:MAX_OUTER),
        # Also show stabilization points
        filter(x -> x > 0, [iter_all_stable, iter_1pct_stable, iter_1pct_tight]),
        [MAX_OUTER]
    )))
    filter!(x -> x <= length(trace_all), show_iters)

    for i in show_iters
        marker = ""
        if i == iter_all_stable; marker *= " ← all<1%"; end
        if i == iter_1pct_stable; marker *= " ← >1%<1%"; end
        if i == iter_1pct_tight; marker *= " ← >1%<0.1%"; end
        @printf("  %6d  %10.4e  %10.4e  %14.6e%s\n",
                i, trace_all[i], trace_1pct[i], trace_ll[i], marker)
    end
end

# ── Summary ──────────────────────────────────────────────────────

println("\n\n" * "=" ^ 140)
println("  SUMMARY — $(length(results)) slow-converging problems")
println("=" ^ 140)

@printf("\n  %-8s %4s %5s │ %8s %8s %8s │ %10s %10s │ %10s\n",
        "Scan", "Cols", "Rows",
        "all<1%", ">1%<1%", ">1%<0.1%",
        "fin_1pct", "fin_all",
        "time(s)")
println("  " * "─" ^ 95)

for r in results
    @printf("  %-8d %4d %5d │ %8s %8s %8s │ %10.4e %10.4e │ %10.3f\n",
            r.scan_idx, r.n_cols, r.n_rows,
            r.iter_all_stable > 0 ? string(r.iter_all_stable) : "never",
            r.iter_1pct_stable > 0 ? string(r.iter_1pct_stable) : "never",
            r.iter_1pct_tight > 0 ? string(r.iter_1pct_tight) : "never",
            r.final_diff_1pct, r.final_diff_all,
            r.time_s)
end

# Aggregate
stable_iters = [r.iter_1pct_stable for r in results if r.iter_1pct_stable > 0]
tight_iters = [r.iter_1pct_tight for r in results if r.iter_1pct_tight > 0]
never_stable = count(r -> r.iter_1pct_stable < 0, results)
never_tight = count(r -> r.iter_1pct_tight < 0, results)

println("\n  Iterations to stabilize weights >1% of max:")
if !isempty(stable_iters)
    @printf("    <1%% change:   mean=%.0f  median=%d  max=%d  (never: %d/%d)\n",
            mean(stable_iters), round(Int, median(stable_iters)), maximum(stable_iters),
            never_stable, length(results))
end
if !isempty(tight_iters)
    @printf("    <0.1%% change: mean=%.0f  median=%d  max=%d  (never: %d/%d)\n",
            mean(tight_iters), round(Int, median(tight_iters)), maximum(tight_iters),
            never_tight, length(results))
end

# Compare: normal convergence is ~6 iters. How much more do these need?
println("\n  Context: normal PMM convergence (LL_tol=1e-6, K=5) takes ~6 outer iters.")
if !isempty(stable_iters)
    @printf("  These problems need %dx--%dx more to stabilize >1%% weights.\n",
            round(Int, minimum(stable_iters) / 6), round(Int, maximum(stable_iters) / 6))
end

println("\n" * "=" ^ 140)
println("  Done.")
println("=" ^ 140)
