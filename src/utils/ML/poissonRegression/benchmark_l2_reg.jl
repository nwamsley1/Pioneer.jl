## Benchmark: L2 regularization effect on PMM convergence
##
## Tests λ_l2 = 0, 1e-6, 1e-4, 1e-2 with K=5, LL_tol=1e-6.
## L2 penalty is applied in y-scaled space: gradient += 2λw_j, hessian += 2λ.
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_l2_reg.jl

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

# ── PMM solver with L2 regularization + LL convergence ───────────

function solvePMM_l2!(sa, μ, y, w, max_outer, ll_tol::Float64, max_inner, λ_l2::Float32)
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
    two_λ = 2f0 * λ_l2

    prev_ll = poissonLogLikelihood(μ, y, sa.m)

    for iter in 1:max_outer
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0

        for col in 1:sa.n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                # Add L2 regularization: d/dw(λw²) = 2λw, d²/dw²(λw²) = 2λ
                L1 += two_λ * w[col]
                L2 += two_λ
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

function extended_run_l2!(sa, μ, y, w, extra_iters, max_inner, λ_l2::Float32)
    w_before = copy(w[1:sa.n])
    ll_before = poissonLogLikelihood(μ, y, sa.m)

    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = maximum(w[1:sa.n])
    two_λ = 2f0 * λ_l2

    for iter in 1:extra_iters
        weight_floor = max_weight * Float32(1e-2)
        max_weight = 0f0
        for col in 1:sa.n
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa, μ, y, col)
                L1 += two_λ * w[col]
                L2 += two_λ
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

    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end

    ll_after = poissonLogLikelihood(μ, y, sa.m)

    # Per-weight analysis
    max_w = maximum(w[1:sa.n])
    max_rel_all = 0.0
    max_rel_1pct = 0.0
    for j in 1:sa.n
        wc = Float64(w_before[j])
        we = Float64(w[j])
        denom = max(abs(we), 1e-10)
        rc = abs(we - wc) / denom
        if rc > max_rel_all
            max_rel_all = rc
        end
        frac = max(wc, we) / max(Float64(max_w), 1e-10)
        if frac > 0.01 && rc > max_rel_1pct
            max_rel_1pct = rc
        end
    end

    return (max_rel_all, max_rel_1pct, ll_before, ll_after)
end

# ── Config ───────────────────────────────────────────────────────

const K = 5
const LL_TOL = 1e-6
const EXTRA_ITERS = 200
const EXT_THRESHOLD = 0.01
const L2_VALUES = [0f0, 1f-6, 1f-4, 1f-2]

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

# Warmup
println("Warming up...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
for λ in L2_VALUES
    w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
    solvePMM_l2!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), LL_TOL, K, λ)
    extended_run_l2!(sa0, μ0, y0, w0, 5, K, λ)
end
println("Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct L2Result
    scan_idx::Int
    n_cols::Int
    λ_l2::Float32
    outer_iters::Int
    total_inner::Int
    time_s::Float64
    final_ll::Float64
    ext_max_rel_all::Float64
    ext_max_rel_1pct::Float64
    ext_ll_before::Float64
    ext_ll_after::Float64
end

all_results = Dict{Float32, Vector{L2Result}}()
for λ in L2_VALUES
    all_results[λ] = L2Result[]
end

# ── Run all problems × λ values ─────────────────────────────────

println("Running $(length(files)) problems × $(length(L2_VALUES)) λ values (K=$K, LL_tol=$LL_TOL)...\n")

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    for λ in L2_VALUES
        w_k = ones(Float32, sa.n)
        μ_k = zeros(Float32, sa.m); y_k = zeros(Float32, sa.m)
        initObserved!(y_k, sa); initMu!(μ_k, sa, w_k)

        local oitr, nitr, fll
        t = @elapsed begin
            oitr, nitr, fll = solvePMM_l2!(sa, μ_k, y_k, w_k, max_outer, LL_TOL, K, λ)
        end

        w_ext = copy(w_k); μ_ext = copy(μ_k); y_ext = copy(y_k)
        ext_all, ext_1pct, ll_bef, ll_aft = extended_run_l2!(sa, μ_ext, y_ext, w_ext, EXTRA_ITERS, K, λ)

        push!(all_results[λ], L2Result(
            data[:scan_idx], sa.n, λ, oitr, nitr, t, fll,
            ext_all, ext_1pct, ll_bef, ll_aft
        ))
    end

    if fi % 20 == 0
        println("  ... $fi / $(length(files))")
    end
end

println("\n  All problems done.\n")

# ── Summary table ────────────────────────────────────────────────

println("=" ^ 140)
println("  L2 REGULARIZATION BENCHMARK — $(length(files)) problems, K=$K, LL_tol=$LL_TOL")
println("=" ^ 140)

println("\n  Per-λ Summary:")
@printf("  %8s │ %6s %6s %6s │ %8s %8s │ %10s │ %14s %14s │ %6s %6s\n",
        "λ_l2",
        "Oi_m", "Oi_md", "Oi_mx",
        "Ni_m", "Ni_md",
        "Time_tot",
        "LL_mean", "LL_median",
        "Unst", "Uns1%")
println("  " * "─" ^ 120)

for λ in L2_VALUES
    res = all_results[λ]
    oi = [r.outer_iters for r in res]
    ni = [r.total_inner for r in res]
    ts = [r.time_s for r in res]
    lls = [r.final_ll for r in res]
    n_unst_all = count(r -> r.ext_max_rel_all > EXT_THRESHOLD, res)
    n_unst_1pct = count(r -> r.ext_max_rel_1pct > EXT_THRESHOLD, res)

    @printf("  %8.0e │ %6.1f %6d %6d │ %8.0f %8d │ %10.4f │ %14.4e %14.4e │ %3d/%d %3d/%d\n",
            λ,
            mean(oi), round(Int, median(oi)), maximum(oi),
            mean(ni), round(Int, median(ni)),
            sum(ts),
            mean(lls), median(lls),
            n_unst_all, length(res),
            n_unst_1pct, length(res))
end

# ── LL comparison: does L2 hurt solution quality? ────────────────

println("\n  LL comparison (does L2 hurt Poisson LL?):")
println("  Per-problem LL difference vs λ=0:")
ref_lls = [r.final_ll for r in all_results[0f0]]
for λ in L2_VALUES[2:end]
    lls = [r.final_ll for r in all_results[λ]]
    diffs = lls .- ref_lls
    rel_diffs = abs.(diffs) ./ (abs.(ref_lls) .+ 1.0)
    @printf("    λ=%.0e:  mean_rel_diff=%.2e  median=%.2e  max=%.2e  (negative=LL decreased)\n",
            λ, mean(rel_diffs), median(rel_diffs), maximum(rel_diffs))
    n_worse = count(d -> d < -1.0, diffs)  # LL decreased by >1.0
    @printf("             #problems where LL decreased by >1.0: %d/%d\n", n_worse, length(diffs))
end

# ── Extended run LL stability ────────────────────────────────────

println("\n  Extended run LL change (median |LL_after - LL_before| / |LL_after|):")
for λ in L2_VALUES
    res = all_results[λ]
    ll_changes = [abs(r.ext_ll_after - r.ext_ll_before) / (abs(r.ext_ll_after) + 1.0) for r in res]
    @printf("    λ=%.0e:  median=%.2e  max=%.2e\n", λ, median(ll_changes), maximum(ll_changes))
end

# ── Deep dive: which problems still unstable at best λ? ──────────

# Find the λ that minimizes unstable-1pct count
best_λ = L2_VALUES[1]
best_count = length(files) + 1
for λ in L2_VALUES
    n = count(r -> r.ext_max_rel_1pct > EXT_THRESHOLD, all_results[λ])
    if n < best_count
        global best_count = n
        global best_λ = λ
    end
end

println("\n  Best λ for stability: λ=$(Printf.@sprintf("%.0e", best_λ)) with $best_count/$(length(files)) unstable (>1% of max)")

# Show per-problem comparison for the 24 problems that were unstable at λ=0
println("\n  Per-problem detail for problems unstable at λ=0 (K=$K, LL_tol=$LL_TOL):")
@printf("  %-8s %4s │", "Scan", "Cols")
for λ in L2_VALUES
    @printf(" λ=%.0e: %6s %6s │", λ, "ext1%", "LL_rΔ")
end
println()
println("  " * "─" ^ (14 + length(L2_VALUES) * 24))

unstable_at_zero = findall(r -> r.ext_max_rel_1pct > EXT_THRESHOLD, all_results[0f0])
for idx in unstable_at_zero
    r0 = all_results[0f0][idx]
    @printf("  %-8d %4d │", r0.scan_idx, r0.n_cols)
    for λ in L2_VALUES
        r = all_results[λ][idx]
        ll_rel = abs(r.ext_ll_after - r.ext_ll_before) / (abs(r.ext_ll_after) + 1.0)
        marker = r.ext_max_rel_1pct > EXT_THRESHOLD ? " *" : "  "
        @printf(" %8.4f %8.2e%s│", r.ext_max_rel_1pct, ll_rel, marker)
    end
    println()
end

println("\n  (* = still unstable)")

println("\n" * "=" ^ 140)
println("  Done.")
println("=" ^ 140)
