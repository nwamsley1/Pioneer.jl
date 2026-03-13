## Diagnose EXT_UNSTABLE: are the moving weights tiny relative to max?
##
## For each problem flagged EXT_UNSTABLE (K=5), compare converged vs extended weights.
## Report per-weight: magnitude, relative to max, and relative change.

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")

# ── Helpers (same as benchmark) ──────────────────────────────────

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

function solvePMM_inner!(sa, μ, y, w, max_outer, rel_conv, max_inner)
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0
    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * Float32(1e-2) : 0f0
        max_weight = 0f0
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
            if w[col] > weight_floor
                rc = δx / max(abs(w[col]), Float32(1e-10))
                rc > _diff && (_diff = rc)
            end
        end
        _diff < rel_conv && break
    end
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
end

function run_extended!(sa, μ, y, w, extra_iters, max_inner)
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
            w[col] > max_weight && (max_weight = w[col])
        end
    end
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
end

# ── Config ───────────────────────────────────────────────────────

const K = 5
const PMM_REL_CONV = Float32(0.001)
const EXTRA_ITERS = 200

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

# Warmup
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePMM_inner!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), PMM_REL_CONV, K)
run_extended!(sa0, μ0, y0, w0, 10, K)
println("Warmup done.\n")

# ── Collect per-weight analysis for all problems ─────────────────

# Aggregate stats across all problems
all_rel_changes = Float64[]       # relative change of each "moved" weight
all_frac_of_max = Float64[]       # w_conv / max_weight for each "moved" weight
all_abs_weights = Float64[]       # absolute weight value (converged)
n_unstable = 0
n_total = 0

println("=" ^ 140)
@printf("  %-8s %4s │ %8s │ %6s %6s │ moved weights: %8s %8s %8s %8s %8s │ %8s\n",
        "Scan", "Cols", "max_w",
        "#moved", "#total",
        "relΔ_med", "relΔ_max", "frac_med", "frac_max", "w_med",
        "LL_Δrel")
println("  " * "─" ^ 136)

for fname in files
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    # Solve to convergence
    w_conv = ones(Float32, sa.n)
    μ = zeros(Float32, sa.m); y = zeros(Float32, sa.m)
    initObserved!(y, sa); initMu!(μ, sa, w_conv)
    solvePMM_inner!(sa, μ, y, w_conv, max_outer, PMM_REL_CONV, K)
    ll_conv = poissonLogLikelihood(μ, y, sa.m)

    # Save converged state, then extend
    w_ext = copy(w_conv)
    μ_ext = copy(μ); y_ext = copy(y)
    run_extended!(sa, μ_ext, y_ext, w_ext, EXTRA_ITERS, K)
    ll_ext = poissonLogLikelihood(μ_ext, y_ext, sa.m)

    # Analyze per-weight changes
    max_w = maximum(w_conv[1:sa.n])
    global n_total += 1

    moved_rel_changes = Float64[]
    moved_frac_of_max = Float64[]
    moved_abs_weights = Float64[]

    for j in 1:sa.n
        wc = Float64(w_conv[j])
        we = Float64(w_ext[j])
        denom = max(abs(we), 1e-10)
        rel_change = abs(we - wc) / denom
        if rel_change > 0.01  # >1% change
            push!(moved_rel_changes, rel_change)
            push!(moved_frac_of_max, wc / max(Float64(max_w), 1e-10))
            push!(moved_abs_weights, wc)
        end
    end

    is_unstable = !isempty(moved_rel_changes)
    if is_unstable
        global n_unstable += 1
        append!(all_rel_changes, moved_rel_changes)
        append!(all_frac_of_max, moved_frac_of_max)
        append!(all_abs_weights, moved_abs_weights)
    end

    ll_rel = abs(ll_ext - ll_conv) / max(abs(ll_conv), 1.0)

    if is_unstable
        @printf("  %-8d %4d │ %8.2e │ %6d %6d │          %8.4f %8.4f %8.4f %8.4f %8.2e │ %8.2e\n",
                data[:scan_idx], sa.n, max_w,
                length(moved_rel_changes), sa.n,
                median(moved_rel_changes), maximum(moved_rel_changes),
                median(moved_frac_of_max), maximum(moved_frac_of_max),
                median(moved_abs_weights),
                ll_rel)
    end
end

# ── Aggregate summary ────────────────────────────────────────────

println("\n" * "=" ^ 140)
println("  AGGREGATE ANALYSIS OF WEIGHTS WITH >1% RELATIVE CHANGE AFTER +$EXTRA_ITERS ITERS")
println("=" ^ 140)

@printf("\n  Unstable problems: %d / %d\n", n_unstable, n_total)
@printf("  Total moved weights across all unstable problems: %d\n\n", length(all_rel_changes))

println("  Relative change of moved weights (|w_ext - w_conv| / |w_ext|):")
@printf("    mean=%.4f  median=%.4f  p90=%.4f  max=%.4f\n",
        mean(all_rel_changes), median(all_rel_changes),
        quantile(all_rel_changes, 0.9), maximum(all_rel_changes))

println("\n  Fraction of max weight (w_conv / max_weight):")
@printf("    mean=%.6f  median=%.6f  p10=%.6f  p90=%.6f  max=%.6f\n",
        mean(all_frac_of_max), median(all_frac_of_max),
        quantile(all_frac_of_max, 0.1), quantile(all_frac_of_max, 0.9),
        maximum(all_frac_of_max))

println("\n  Absolute weight magnitude (converged):")
@printf("    mean=%.4e  median=%.4e  min=%.4e  max=%.4e\n",
        mean(all_abs_weights), median(all_abs_weights),
        minimum(all_abs_weights), maximum(all_abs_weights))

# Bucket by fraction of max
println("\n  Breakdown by fraction of max weight:")
for (lo, hi, label) in [
    (0.0, 1e-4, "< 0.01%"),
    (1e-4, 1e-3, "0.01-0.1%"),
    (1e-3, 1e-2, "0.1-1%"),
    (1e-2, 1e-1, "1-10%"),
    (1e-1, 1.0, "10-100%"),
    (1.0, Inf, "≥100% (IS max)")
]
    mask = lo .<= all_frac_of_max .< hi
    n = count(mask)
    n == 0 && continue
    rc = all_rel_changes[mask]
    @printf("    %-16s  n=%4d  relΔ: mean=%.4f  median=%.4f  max=%.4f\n",
            label, n, mean(rc), median(rc), maximum(rc))
end

println("\n" * "=" ^ 140)
println("  Done.")
println("=" ^ 140)
