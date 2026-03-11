# Experiment: Does scaling y by a constant fix MM convergence on production data?
#
# Poisson MLE:  max_w Σ [y_i log(Aw)_i - (Aw)_i]  s.t. w ≥ 0
# If y' = y/c, optimal w' = w/c.  Solve for w', then w = c * w'.
#
# Run:  julia --project=<Pioneer root> test_mm_scaling.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Load real data ────────────────────────────────────────────────
data = deserialize("/Users/nathanwamsley/Desktop/solveHuber_inputs_scan342335.jls")

sa = Main.SparseArray(
    data[:n_vals], data[:m], data[:n],
    Vector{Int64}(data[:rowval]), Vector{UInt16}(data[:colval]), data[:nzval],
    data[:matched], data[:isotope], data[:x], Vector{Int64}(data[:colptr])
)

nr_iter  = data[:max_iter_newton]
bs_iter  = data[:max_iter_bisection]
max_outer = data[:max_iter_outer]
nr_acc   = data[:accuracy_newton]
bs_acc   = data[:accuracy_bisection]
rel_conv = data[:max_diff]

println("=" ^ 90)
println("  Y-Scaling Experiment: scan342335")
println("  m=$(sa.m), n=$(sa.n), rel_conv=$rel_conv")
println("=" ^ 90)

# ── Eval functions ────────────────────────────────────────────────
function sse_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += (Float64(y[i]) - Float64(μ[i]))^2; end; s
end
function poisson_ll_from_mu(μ::Vector{Float32}, y::Vector{Float32}, m::Int)
    ll = 0.0
    @inbounds for i in 1:m
        μ_i = max(Float64(μ[i]), 1e-6)
        ll += Float64(y[i]) * log(μ_i) - μ_i
    end; ll
end

# ── MM solver (multi-step, returns per-iter trace) ────────────────
function solve_mm_traced(sa_l, init_weights, y_scale::Float32;
                         max_inner=5, mm_rel_conv=rel_conv, mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m)
    yy = zeros(Float32, m)
    w = copy(init_weights)

    initObserved!(yy, sa_l)
    # Scale observed values
    if y_scale != 1f0
        @inbounds for i in 1:m
            yy[i] /= y_scale
        end
    end
    initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)

    lls  = Float64[]
    diffs = Float64[]

    # Compute initial LL in original scale
    if y_scale != 1f0
        μ_orig = μ .* y_scale
        y_orig = yy .* y_scale
        push!(lls, poisson_ll_from_mu(μ_orig, y_orig, m))
    else
        push!(lls, poisson_ll_from_mu(μ, yy, m))
    end

    for iter in 1:mm_max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            X_before = w[col]
            for _k in 1:max_inner
                L1, L2 = getPoissonDerivativesObs!(sa_l, μ, yy, col)
                if L2 <= ε || isnan(L1); break; end
                X0 = w[col]
                w[col] = max(w[col] - L1 / L2, 0f0)
                updateMu!(sa_l, μ, col, w[col], X0)
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

        # LL in original scale for comparison
        if y_scale != 1f0
            μ_orig = μ .* y_scale
            y_orig = yy .* y_scale
            push!(lls, poisson_ll_from_mu(μ_orig, y_orig, m))
        else
            push!(lls, poisson_ll_from_mu(μ, yy, m))
        end
        push!(diffs, Float64(_diff))
        _diff < mm_rel_conv && break
    end

    # Rescale weights back to original
    if y_scale != 1f0
        @inbounds for j in 1:n
            w[j] *= y_scale
        end
    end

    # Eval in original scale
    μ_eval = zeros(Float32, m)
    y_eval = zeros(Float32, m)
    initObserved!(y_eval, sa_l)
    initMu!(μ_eval, sa_l, w)
    sse = sse_from_mu(μ_eval, y_eval, m)
    pll = poisson_ll_from_mu(μ_eval, y_eval, m)
    nnz = count(w[1:n] .> 0.01)

    return w, length(diffs), sse, pll, nnz, lls, diffs
end

# ── Reference: ObsHess Newton+bisection ───────────────────────────
function solve_obshess(sa_l, init_weights)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l); initMu!(μ, sa_l, w)
    max_weight = Float32(0)
    iters = 0
    for iter in 1:max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            δx = abs(newton_bisection_poisson_obs!(sa_l, μ, yy, w, col,
                        nr_iter, bs_iter, nr_acc, bs_acc, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    sse = sse_from_mu(μ, yy, m)
    pll = poisson_ll_from_mu(μ, yy, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ══════════════════════════════════════════════════════════════════
#  Run experiments
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))

# Get y stats for scaling choices
y_vals = zeros(Float32, sa.m)
initObserved!(y_vals, sa)
y_nz = y_vals[y_vals .> 0]
@printf("\ny statistics:\n")
@printf("  range: [%.0f, %.0f], mean: %.1f, median: %.1f\n",
        minimum(y_vals), maximum(y_vals), mean(y_vals), median(y_nz))
@printf("  nonzero: %d / %d\n", length(y_nz), sa.m)

# ── Baselines ─────────────────────────────────────────────────────
println("\n── Baselines (cold-start, no scaling) ──")

t = @elapsed w_o, it_o, sse_o, pll_o, nnz_o = solve_obshess(sa, cold_w)
@printf("  ObsHess (Newton):  %4d iters, %.4f sec, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_o, t, pll_o, sse_o, nnz_o)

t = @elapsed w_mm0, it_mm0, sse_mm0, pll_mm0, nnz_mm0, _, _ = solve_mm_traced(sa, cold_w, 1f0)
@printf("  MM K=5 (no scale): %4d iters, %.4f sec, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_mm0, t, pll_mm0, sse_mm0, nnz_mm0)

# ── Scaling experiments ───────────────────────────────────────────
println("\n── MM K=5 with y-scaling (cold-start) ──")
@printf("  %-28s  %5s  %8s  %16s  %14s  %5s  %14s\n",
        "Scale factor", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz", "w_max")
println("  " * "─" ^ 100)

scale_factors = [1f0, Float32(mean(y_nz)), Float32(median(y_nz)),
                 Float32(maximum(y_vals)), Float32(sqrt(mean(y_nz))),
                 1f3, 1f4, 1f5, 1f6]

for c in scale_factors
    t = @elapsed w_s, it_s, sse_s, pll_s, nnz_s, _, _ = solve_mm_traced(sa, cold_w, c)
    w_max = maximum(w_s[1:sa.n])
    label = if c == 1f0
        "1 (no scale)"
    elseif c == Float32(mean(y_nz))
        @sprintf("mean(y)=%.0f", c)
    elseif c == Float32(median(y_nz))
        @sprintf("median(y)=%.0f", c)
    elseif c == Float32(maximum(y_vals))
        @sprintf("max(y)=%.0f", c)
    elseif c == Float32(sqrt(mean(y_nz)))
        @sprintf("sqrt(mean)=%.1f", c)
    else
        @sprintf("%.0e", c)
    end
    @printf("  %-28s  %5d  %8.4f  %16.6e  %14.6e  %5d  %14.4e\n",
            label, it_s, t, pll_s, sse_s, nnz_s, w_max)
end

# ── Best scaling: convergence trace ──────────────────────────────
# Find best scale factor by LL
# Use max(y) as best scale factor based on results above
best_c = Float32(maximum(y_vals))

@printf("\n── Convergence trace for best scale c=%.0f (cold-start, K=5) ──\n", best_c)
w_best, it_best, sse_best, pll_best, nnz_best, lls_best, diffs_best = solve_mm_traced(sa, cold_w, best_c)

@printf("  %5s  %16s  %12s\n", "Iter", "Poisson LL", "max_rel_diff")
println("  " * "─" ^ 40)
n_trace = length(diffs_best)
show_iters = sort(unique(vcat(
    collect(1:min(20, n_trace)),
    collect(25:5:min(50, n_trace)),
    collect(50:25:min(200, n_trace)),
    collect(200:100:n_trace),
    [n_trace]
)))
for i in show_iters
    i > n_trace && continue
    @printf("  %5d  %16.6e  %12.4e\n", i, lls_best[i+1], diffs_best[i])
end

ll_mono = all(lls_best[i+1] >= lls_best[i] - 1.0 for i in 1:length(lls_best)-1)
@printf("\n  LL monotonic: %s\n", ll_mono)

# ── Compare best scaled MM vs ObsHess ────────────────────────────
println("\n── Best scaled MM vs ObsHess ──")
@printf("  ObsHess:       iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_o, pll_o, sse_o, nnz_o)
@printf("  MM(c=%.0f):  iters=%d, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        best_c, it_best, pll_best, sse_best, nnz_best)

wo_f = Float64.(w_o[1:sa.n])
wb_f = Float64.(w_best[1:sa.n])
@printf("  Weight correlation: %.6f\n", cor(wo_f, wb_f))

println("\n" * "=" ^ 90)
println("  Done.")
println("=" ^ 90)
