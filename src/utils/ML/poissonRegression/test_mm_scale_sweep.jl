# Systematic y-scaling investigation:
# 1. What target max(y') gives best MM convergence?
# 2. Does scaling help Huber/OLS too?
# 3. Theory: what determines the right scale?
#
# Run:  julia --project=<Pioneer root> test_mm_scale_sweep.jl

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")
include("spectralLinearRegression_reference.jl")

# ── Load real data ────────────────────────────────────────────────
data = deserialize("/Users/n.t.wamsley/Desktop/solveHuber_inputs_scan342335.jls")

sa = Main.SparseArray(
    data[:n_vals], data[:m], data[:n],
    Vector{Int64}(data[:rowval]), Vector{UInt16}(data[:colval]), data[:nzval],
    data[:matched], data[:isotope], data[:x], Vector{Int64}(data[:colptr])
)

δ_param  = data[:delta]
λ_param  = data[:lambda]
nr_iter  = data[:max_iter_newton]
bs_iter  = data[:max_iter_bisection]
max_outer = data[:max_iter_outer]
nr_acc   = data[:accuracy_newton]
bs_acc   = data[:accuracy_bisection]
rel_conv = data[:max_diff]

println("=" ^ 100)
println("  Y-Scaling Sweep: scan342335")
println("  m=$(sa.m), n=$(sa.n), rel_conv=$rel_conv")
println("=" ^ 100)

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

# ── Problem characterization ─────────────────────────────────────
y_vals = zeros(Float32, sa.m)
initObserved!(y_vals, sa)
y_nz = y_vals[y_vals .> 0]

# Compute μ_init = A * ones(n) to understand the scale mismatch
μ_init = zeros(Float32, sa.m)
w_ones = ones(Float32, sa.n)
initMu!(μ_init, sa, w_ones)

# Compute column densities (nnz per column)
col_nnz = [sa.colptr[j+1] - sa.colptr[j] for j in 1:sa.n]

println("\n── Problem characterization ──")
@printf("  y:       range=[%.0f, %.0f], mean=%.1f, median(nz)=%.1f, nnz=%d/%d\n",
        minimum(y_vals), maximum(y_vals), mean(y_vals),
        length(y_nz) > 0 ? median(y_nz) : 0.0, length(y_nz), sa.m)
@printf("  μ_init:  range=[%.4f, %.4f], mean=%.4f (with w=1)\n",
        minimum(μ_init[1:sa.m]), maximum(μ_init[1:sa.m]), mean(μ_init[1:sa.m]))
@printf("  y/μ_init ratio at nonzero y: mean=%.1f, max=%.1f\n",
        mean(y_nz ./ max.(μ_init[findall(y_vals .> 0)], 1f-6)),
        maximum(y_nz ./ max.(μ_init[findall(y_vals .> 0)], 1f-6)))
@printf("  Col nnz: range=[%d, %d], mean=%.1f\n",
        minimum(col_nnz), maximum(col_nnz), mean(col_nnz))

# True weight scale (from warm start)
warm_w = copy(data[:weights])
w_active = warm_w[1:sa.n]
w_nz = w_active[w_active .> 0]
@printf("  True weights: range=[%.2e, %.2e], mean=%.2e, median=%.2e\n",
        minimum(w_nz), maximum(w_nz), mean(w_nz), median(w_nz))
@printf("  Dynamic range: %.1f orders of magnitude\n",
        log10(maximum(w_nz) / max(minimum(w_nz), 1f-10)))

# ── MM solver with y-scaling ──────────────────────────────────────
function solve_mm_scaled(sa_l, init_weights, y_scale::Float32;
                         max_inner=5, mm_rel_conv=rel_conv, mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n
    μ = zeros(Float32, m); yy = zeros(Float32, m); w = copy(init_weights)
    initObserved!(yy, sa_l)
    if y_scale != 1f0
        @inbounds for i in 1:m; yy[i] /= y_scale; end
    end
    initMu!(μ, sa_l, w)

    ε = Float32(POISSON_MU_FLOOR)
    max_weight = Float32(0)
    iters = 0
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
        iters = iter
        _diff < mm_rel_conv && break
    end
    # Rescale
    if y_scale != 1f0
        @inbounds for j in 1:n; w[j] *= y_scale; end
    end
    # Eval in original scale
    μ_e = zeros(Float32, m); y_e = zeros(Float32, m)
    initObserved!(y_e, sa_l); initMu!(μ_e, sa_l, w)
    sse = sse_from_mu(μ_e, y_e, m)
    pll = poisson_ll_from_mu(μ_e, y_e, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ── Huber solver with y-scaling ───────────────────────────────────
# For Huber: scaling y by c means we also need to scale δ by c
# (the Huber transition parameter) to keep the same relative behavior.
# Residuals r = Aw - y scale by c, so δ must scale too.
function solve_huber_scaled(sa_l, init_weights, y_scale::Float32;
                            mm_max_outer=max_outer)
    m, n = sa_l.m, sa_l.n

    # Build a modified SparseArray with scaled x values
    if y_scale != 1f0
        x_scaled = copy(sa_l.x)
        @inbounds for i in eachindex(x_scaled); x_scaled[i] /= y_scale; end
        sa_s = Main.SparseArray(sa_l.n_vals, sa_l.m, sa_l.n,
                                sa_l.rowval, sa_l.colval, sa_l.nzval,
                                sa_l.matched, sa_l.isotope, x_scaled, sa_l.colptr)
    else
        sa_s = sa_l
    end

    r = zeros(Float32, m); w = copy(init_weights)
    y = zeros(Float32, m)
    @inbounds for i in 1:sa_s.n_vals
        if iszero(y[sa_s.rowval[i]]); y[sa_s.rowval[i]] = sa_s.x[i]; end
    end
    initResiduals_plain!(r, sa_s, w)

    # Scale δ proportionally
    δ_s = y_scale != 1f0 ? δ_param / y_scale : δ_param

    reg = NoNorm()
    max_weight = Float32(0)
    iters = 0
    for iter in 1:mm_max_outer
        _diff = Float32(0)
        weight_floor = iter > 5 ? max_weight * Float32(1e-7) : Float32(0)
        max_weight = Float32(0)
        for col in 1:n
            δx = abs(newton_bisection!(sa_s, r, w, col, δ_s, λ_param,
                        nr_iter, bs_iter, nr_acc, bs_acc, reg, rel_conv))
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / abs(w[col]); rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    # Rescale weights back
    if y_scale != 1f0
        @inbounds for j in 1:n; w[j] *= y_scale; end
    end
    # Eval in original scale
    μ_e = zeros(Float32, m); y_e = zeros(Float32, m)
    initObserved!(y_e, sa_l); initMu!(μ_e, sa_l, w)
    sse = sse_from_mu(μ_e, y_e, m)
    pll = poisson_ll_from_mu(μ_e, y_e, m)
    nnz = count(w[1:n] .> 0.01)
    return w, iters, sse, pll, nnz
end

# ── ObsHess reference ─────────────────────────────────────────────
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
#  Experiment 1: MM — sweep target max(y') values
# ══════════════════════════════════════════════════════════════════

cold_w = ones(Float32, length(data[:weights]))
max_y = maximum(y_vals)

println("\n" * "=" ^ 100)
println("  EXPERIMENT 1: MM convergence vs target max(y')")
println("  Scale c = max(y)/target, so y' = y/c has max(y') = target")
println("=" ^ 100)

println("\n── Reference (no scaling) ──")
t_ref = @elapsed (w_ref, it_ref, sse_ref, pll_ref, nnz_ref) = solve_obshess(sa, cold_w)
@printf("  ObsHess: %4d iters, %.4f sec, PoissonLL=%.6e, SSE=%.6e, nnz=%d\n",
        it_ref, t_ref, pll_ref, sse_ref, nnz_ref)

@printf("\n  %-18s  %12s  %5s  %8s  %16s  %14s  %5s  %12s\n",
        "target max(y')", "scale c", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz", "scaled w_max")
println("  " * "─" ^ 105)

# Sweep: what target value for max(y') works?
targets = [0.01, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0,
           5000.0, 10000.0, 50000.0, 100000.0, 500000.0, max_y]
for target in targets
    c = Float32(max_y / target)
    t = @elapsed (w_s, it_s, sse_s, pll_s, nnz_s) = solve_mm_scaled(sa, cold_w, c)
    scaled_w_max = maximum(w_s[1:sa.n]) / c  # max weight in scaled space
    label = target == max_y ? @sprintf("%.0f (no scale)", target) : @sprintf("%.2g", target)
    @printf("  %-18s  %12.1f  %5d  %8.4f  %16.6e  %14.6e  %5d  %12.2e\n",
            label, c, it_s, t, pll_s, sse_s, nnz_s, scaled_w_max)
end

# ══════════════════════════════════════════════════════════════════
#  Experiment 2: Does scaling help Huber/OLS?
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 100)
println("  EXPERIMENT 2: Does y-scaling help Huber convergence?")
println("=" ^ 100)

@printf("\n  %-18s  %12s  %5s  %8s  %16s  %14s  %5s\n",
        "target max(y')", "scale c", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 90)

huber_targets = [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, max_y]
for target in huber_targets
    c = Float32(max_y / target)
    t = @elapsed (w_s, it_s, sse_s, pll_s, nnz_s) = solve_huber_scaled(sa, cold_w, c)
    label = target == max_y ? @sprintf("%.0f (no scale)", target) : @sprintf("%.2g", target)
    @printf("  %-18s  %12.1f  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            label, c, it_s, t, pll_s, sse_s, nnz_s)
end

# ══════════════════════════════════════════════════════════════════
#  Experiment 3: Theory check — what's the critical ratio?
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 100)
println("  THEORETICAL ANALYSIS")
println("=" ^ 100)

# The observed Hessian step for coordinate j is:
#   δ_j = L1_j / L2_j
#   L1 = Σ A_ij (1 - y_i/μ_i)
#   L2 = Σ A_ij² y_i / μ_i²
#
# At initialization (w=1, μ_i = Σ_j A_ij):
#   L2 ∝ y_i / μ_i²
#   Step ∝ L1/L2 ∝ μ_i²/y_i * (1 - y_i/μ_i) ≈ -μ_i/1 when y_i >> μ_i
#
# So the step is O(μ_init) regardless of y_i. But the TARGET weight
# is w_true ~ y_i / (A_ij * n_overlapping_columns).
#
# Steps needed to reach target: w_true / step ≈ y_max / μ_init²
#
# With scaling y' = y/c: step is same O(μ_init), but target is w'/c.
# Steps needed: (w_true/c) / μ_init ≈ y_max / (c * μ_init)
#
# For convergence in O(1) sweeps: need c ≈ y_max / μ_init
# Since μ_init ≈ mean column value * nnz_per_row ≈ O(1-10):
# c ≈ max(y) works!

# Compute μ_init statistics for rows with nonzero y
nz_rows = findall(y_vals .> 0)
μ_init_at_y = μ_init[nz_rows]
y_at_nz = y_vals[nz_rows]
ratio = y_at_nz ./ max.(μ_init_at_y, 1f-6)

@printf("\n  At initialization (w=1), for rows with y>0:\n")
@printf("    μ_init:  mean=%.2f, median=%.2f, range=[%.4f, %.2f]\n",
        mean(μ_init_at_y), median(μ_init_at_y),
        minimum(μ_init_at_y), maximum(μ_init_at_y))
@printf("    y/μ:     mean=%.1f, median=%.1f, max=%.1f\n",
        mean(ratio), median(ratio), maximum(ratio))
@printf("    log10(y/μ) range: [%.1f, %.1f]\n",
        minimum(log10.(max.(ratio, 1f-6))), maximum(log10.(max.(ratio, 1f-6))))

@printf("\n  Theory: scale c ≈ max(y)/μ_init_typical ≈ %.0f / %.1f ≈ %.0f\n",
        max_y, median(μ_init_at_y), max_y / median(μ_init_at_y))
@printf("  Equivalently: set max(y') ≈ μ_init_typical ≈ %.1f\n",
        median(μ_init_at_y))
@printf("  This makes y'/μ ≈ O(1) so L2 ≈ O(A²/μ) (Fisher) and steps are well-scaled\n")

# ══════════════════════════════════════════════════════════════════
#  Experiment 4: Warm-start with best scaling
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 100)
println("  EXPERIMENT 4: Warm-start with best scales")
println("=" ^ 100)

@printf("\n  %-25s  %5s  %8s  %16s  %14s  %5s\n",
        "Config", "Iters", "Time(s)", "Poisson LL", "SSE", "nnz")
println("  " * "─" ^ 80)

t = @elapsed (w_s, it_s, sse_s, pll_s, nnz_s) = solve_obshess(sa, warm_w)
@printf("  %-25s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
        "ObsHess (no scale)", it_s, t, pll_s, sse_s, nnz_s)

for target in [1.0, 10.0, 100.0, max_y]
    c = Float32(max_y / target)
    # Warm-start weights need scaling too
    w_warm_scaled = copy(warm_w)
    @inbounds for j in 1:sa.n; w_warm_scaled[j] /= c; end
    t = @elapsed (w_s, it_s, sse_s, pll_s, nnz_s) = solve_mm_scaled(sa, w_warm_scaled, c)
    label = target == max_y ? "MM warm (no scale)" : @sprintf("MM warm (tgt=%.0f)", target)
    @printf("  %-25s  %5d  %8.4f  %16.6e  %14.6e  %5d\n",
            label, it_s, t, pll_s, sse_s, nnz_s)
end

# Weight correlation between best scaled MM and ObsHess
best_c = Float32(max_y / 1.0)
w_warm_sc = copy(warm_w)
@inbounds for j in 1:sa.n; w_warm_sc[j] /= best_c; end
w_mm_best, _, _, _, _ = solve_mm_scaled(sa, w_warm_sc, best_c)
w_obs_best, _, _, _, _ = solve_obshess(sa, warm_w)

wo = Float64.(w_obs_best[1:sa.n])
wm = Float64.(w_mm_best[1:sa.n])
@printf("\n  ObsHess vs MM(tgt=1) warm-start weight correlation: %.6f\n", cor(wo, wm))

println("\n" * "=" ^ 100)
println("  Done.")
println("=" ^ 100)
