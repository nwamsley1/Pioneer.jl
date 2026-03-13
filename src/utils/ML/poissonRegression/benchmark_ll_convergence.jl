## Benchmark: LL-based convergence criterion × inner iterations K
##
## Tests K=1,3,5,10 × LL_tol=1e-3,1e-4,1e-6 (relative LL change per sweep).
## For each combo: outer iters, total Newton steps, timing, final LL, EXT_UNSTABLE.
## Then deep-dives on remaining unstable problems where weights >1% of max don't settle.
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_ll_convergence.jl

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

# ── PMM solver with LL-based convergence ─────────────────────────

"""
    solvePMM_ll_conv!(sa, μ, y, w, max_outer, ll_tol, max_inner)
        → (outer_iters, total_inner, final_ll)

PMM solver that converges based on relative LL change per sweep:
    |LL_new - LL_old| / (|LL_new| + 1) < ll_tol
"""
function solvePMM_ll_conv!(sa, μ, y, w, max_outer, ll_tol::Float64, max_inner)
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

    prev_ll = poissonLogLikelihood(μ, y, sa.m)

    for iter in 1:max_outer
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0

        for col in 1:sa.n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                if L2 <= ε || isnan(L1); break; end
                total_inner += 1
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa, μ, col, w[col], X0)
                if iszero(w[col]) || abs(w[col] - X0) / max(abs(w[col]), Float32(1e-10)) < Float32(1e-3)
                    break
                end
            end
            w[col] > max_weight && (max_weight = w[col])
        end

        cur_ll = poissonLogLikelihood(μ, y, sa.m)
        outer_iters = iter
        rel_ll_change = abs(cur_ll - prev_ll) / (abs(cur_ll) + 1.0)
        prev_ll = cur_ll

        if rel_ll_change < ll_tol
            break
        end
    end

    final_ll = prev_ll

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
        final_ll = poissonLogLikelihood(μ, y, sa.m)
    end

    return (outer_iters, total_inner, final_ll)
end

"""
    extended_run_ll!(sa, μ, y, w, extra_iters, max_inner)
        → (max_rel_change, max_rel_change_above1pct, ll_before, ll_after)

Run extra iterations and track:
- max_rel_change: across all weights above weight_floor
- max_rel_change_above1pct: only for weights >1% of max (the ones that matter)
"""
function extended_run_ll!(sa, μ, y, w, extra_iters, max_inner)
    w_before = copy(w[1:sa.n])
    ll_before = poissonLogLikelihood(μ, y, sa.m)

    # Y-scaling
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = maximum(w[1:sa.n])

    for iter in 1:extra_iters
        weight_floor = max_weight * Float32(1e-2)
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
    end

    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end

    ll_after = poissonLogLikelihood(μ, y, sa.m)

    # Compute per-weight changes
    max_w = maximum(w[1:sa.n])
    max_rel_all = 0.0
    max_rel_above1pct = 0.0
    for j in 1:sa.n
        wc = Float64(w_before[j])
        we = Float64(w[j])
        denom = max(abs(we), 1e-10)
        rc = abs(we - wc) / denom
        if rc > max_rel_all
            max_rel_all = rc
        end
        frac = max(wc, we) / max(Float64(max_w), 1e-10)
        if frac > 0.01 && rc > max_rel_above1pct
            max_rel_above1pct = rc
        end
    end

    return (max_rel_all, max_rel_above1pct, ll_before, ll_after)
end

# ── Config ───────────────────────────────────────────────────────

const K_VALUES = [1, 3, 5, 10]
const LL_TOLS = [1e-3, 1e-4, 1e-6]
const EXTRA_ITERS = 200
const EXT_THRESHOLD = 0.01

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

# Warmup
println("Warming up...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
for K in K_VALUES
    for tol in LL_TOLS
        w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
        solvePMM_ll_conv!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), tol, K)
        extended_run_ll!(sa0, μ0, y0, w0, 5, K)
    end
end
println("Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct LLResult
    scan_idx::Int
    n_cols::Int
    K::Int
    ll_tol::Float64
    outer_iters::Int
    total_inner::Int
    time_s::Float64
    final_ll::Float64
    ext_max_rel_all::Float64
    ext_max_rel_1pct::Float64
    ext_ll_before::Float64
    ext_ll_after::Float64
end

# Key: (K, ll_tol) → Vector{LLResult}
all_results = Dict{Tuple{Int,Float64}, Vector{LLResult}}()
for K in K_VALUES, tol in LL_TOLS
    all_results[(K, tol)] = LLResult[]
end

# ── Run all problems × K × LL_tol ───────────────────────────────

println("Running $(length(files)) problems × $(length(K_VALUES)) K × $(length(LL_TOLS)) LL_tol ...\n")

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    for K in K_VALUES
        for tol in LL_TOLS
            w_k = ones(Float32, sa.n)
            μ_k = zeros(Float32, sa.m); y_k = zeros(Float32, sa.m)
            initObserved!(y_k, sa); initMu!(μ_k, sa, w_k)

            local oitr, nitr, fll
            t = @elapsed begin
                oitr, nitr, fll = solvePMM_ll_conv!(sa, μ_k, y_k, w_k, max_outer, tol, K)
            end

            # Extended stability
            w_ext = copy(w_k); μ_ext = copy(μ_k); y_ext = copy(y_k)
            ext_all, ext_1pct, ll_bef, ll_aft = extended_run_ll!(sa, μ_ext, y_ext, w_ext, EXTRA_ITERS, K)

            push!(all_results[(K, tol)], LLResult(
                data[:scan_idx], sa.n, K, tol,
                oitr, nitr, t, fll,
                ext_all, ext_1pct, ll_bef, ll_aft
            ))
        end
    end

    if fi % 20 == 0
        println("  ... $fi / $(length(files))")
    end
end

println("\n  All problems done.\n")

# ── Summary table: K × LL_tol ────────────────────────────────────

println("=" ^ 160)
println("  LL-BASED CONVERGENCE BENCHMARK — $(length(files)) problems")
println("=" ^ 160)

println("\n  Per (K, LL_tol) Summary:")
@printf("  %3s  %8s │ %6s %6s %6s │ %8s %8s %8s │ %10s %10s │ %6s %6s\n",
        "K", "LL_tol",
        "Oi_m", "Oi_md", "Oi_mx",
        "Ni_m", "Ni_md", "Ni_mx",
        "Time_m(s)", "Time_tot",
        "Unst", "Uns1%")
println("  " * "─" ^ 110)

for K in K_VALUES
    for tol in LL_TOLS
        res = all_results[(K, tol)]
        oi = [r.outer_iters for r in res]
        ni = [r.total_inner for r in res]
        ts = [r.time_s for r in res]
        n_unst_all = count(r -> r.ext_max_rel_all > EXT_THRESHOLD, res)
        n_unst_1pct = count(r -> r.ext_max_rel_1pct > EXT_THRESHOLD, res)

        @printf("  %3d  %8.0e │ %6.1f %6d %6d │ %8.0f %8d %8d │ %10.6f %10.4f │ %3d/%d %3d/%d\n",
                K, tol,
                mean(oi), round(Int, median(oi)), maximum(oi),
                mean(ni), round(Int, median(ni)), maximum(ni),
                mean(ts), sum(ts),
                n_unst_all, length(res),
                n_unst_1pct, length(res))
    end
    println()
end

# ── LL agreement: does tighter tol actually improve LL? ──────────

println("  Final LL comparison (median across problems):")
@printf("  %3s │", "K")
for tol in LL_TOLS
    @printf("  LL_tol=%.0e  ", tol)
end
println()
println("  " * "─" ^ 60)
for K in K_VALUES
    @printf("  %3d │", K)
    for tol in LL_TOLS
        lls = [r.final_ll for r in all_results[(K, tol)]]
        @printf("  %14.6e", median(lls))
    end
    println()
end

# ── LL change from extended run ──────────────────────────────────

println("\n  LL change from +$EXTRA_ITERS extended iters (median |LL_after - LL_before| / |LL_after|):")
@printf("  %3s │", "K")
for tol in LL_TOLS
    @printf("  LL_tol=%.0e  ", tol)
end
println()
println("  " * "─" ^ 60)
for K in K_VALUES
    @printf("  %3d │", K)
    for tol in LL_TOLS
        res = all_results[(K, tol)]
        ll_changes = [abs(r.ext_ll_after - r.ext_ll_before) / (abs(r.ext_ll_after) + 1.0) for r in res]
        @printf("  %14.2e", median(ll_changes))
    end
    println()
end

# ══════════════════════════════════════════════════════════════════
# DEEP DIVE: problems where weights >1% of max still move >1%
# ══════════════════════════════════════════════════════════════════

# Use K=5, LL_tol=1e-6 as our tightest setting
println("\n" * "=" ^ 160)
println("  DEEP DIVE: Weights >1%% of max that move >1%% after +$EXTRA_ITERS iters")
println("  Settings: K=5, LL_tol=1e-6")
println("=" ^ 160)

deep_K = 5
deep_tol = 1e-6

# For the deep dive, re-run with full per-weight tracking
println()
n_deep = 0
for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    # Solve
    w_conv = ones(Float32, sa.n)
    μ = zeros(Float32, sa.m); y = zeros(Float32, sa.m)
    initObserved!(y, sa); initMu!(μ, sa, w_conv)
    solvePMM_ll_conv!(sa, μ, y, w_conv, max_outer, deep_tol, deep_K)
    ll_conv = poissonLogLikelihood(μ, y, sa.m)

    # Extended
    w_ext = copy(w_conv); μ_ext = copy(μ); y_ext = copy(y)
    ext_all, ext_1pct, _, ll_ext = extended_run_ll!(sa, μ_ext, y_ext, w_ext, EXTRA_ITERS, deep_K)

    ext_1pct > EXT_THRESHOLD || continue

    # This problem has weights >1% of max moving >1%
    global n_deep += 1
    max_w = maximum(w_ext[1:sa.n])

    println("  ── Scan $(data[:scan_idx]) │ $(sa.n) cols × $(sa.m) rows │ max_w=$(Printf.@sprintf("%.2e", max_w)) │ LL_conv=$(Printf.@sprintf("%.6e", ll_conv)) → LL_ext=$(Printf.@sprintf("%.6e", ll_ext)) ──")
    @printf("     ext_max_rel_all=%.4f  ext_max_rel_>1%%=%.4f  LL_rel_Δ=%.2e\n",
            ext_all, ext_1pct,
            abs(ll_ext - ll_conv) / (abs(ll_ext) + 1.0))

    # List all weights >1% of max that moved >1%
    movers = Tuple{Int, Float64, Float64, Float64, Float64}[]  # (col, w_conv, w_ext, rel_change, frac_of_max)
    for j in 1:sa.n
        wc = Float64(w_conv[j])
        we = Float64(w_ext[j])
        frac = max(wc, we) / max(Float64(max_w), 1e-10)
        frac > 0.01 || continue
        denom = max(abs(we), 1e-10)
        rc = abs(we - wc) / denom
        rc > 0.01 || continue
        push!(movers, (j, wc, we, rc, frac))
    end

    sort!(movers, by=x -> -x[5])  # sort by frac_of_max descending
    @printf("     %4s  %12s  %12s  %10s  %10s  %6s\n",
            "Col", "w_conv", "w_ext", "rel_change", "frac_max", "nnz_col")
    for (col, wc, we, rc, frac) in movers
        nnz_col = sa.colptr[col+1] - sa.colptr[col]
        @printf("     %4d  %12.4e  %12.4e  %10.4f  %10.4f  %6d\n",
                col, wc, we, rc, frac, nnz_col)
    end

    # Also show the top 5 weights (by magnitude) that did NOT move, for contrast
    stable = Tuple{Int, Float64, Float64, Float64}[]
    for j in 1:sa.n
        wc = Float64(w_conv[j])
        we = Float64(w_ext[j])
        denom = max(abs(we), 1e-10)
        rc = abs(we - wc) / denom
        frac = max(wc, we) / max(Float64(max_w), 1e-10)
        if frac > 0.01 && rc <= 0.01
            push!(stable, (j, we, rc, frac))
        end
    end
    sort!(stable, by=x -> -x[4])
    if !isempty(stable)
        println("     (stable weights >1% of max for contrast, top 5):")
        for (col, we, rc, frac) in stable[1:min(5, length(stable))]
            @printf("     %4d  %12s  %12.4e  %10.4f  %10.4f\n",
                    col, "", we, rc, frac)
        end
    end
    println()
end

if n_deep == 0
    println("  No problems with weights >1% of max moving >1% after extended run!")
end

println("  Total deep-dive problems: $n_deep / $(length(files))")

println("\n" * "=" ^ 160)
println("  Done.")
println("=" ^ 160)
