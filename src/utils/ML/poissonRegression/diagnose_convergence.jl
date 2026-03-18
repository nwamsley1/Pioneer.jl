## Convergence Diagnostic for OLS, PMM, and Huber solvers
##
## For each of the 100 benchmark problems, this script reports:
##   A) Final _diff value and whether each solver truly converged
##   B) Extended PMM run (200 extra iters from converged solution)
##   C) Poisson log-likelihood for all three solvers
##   D) Flags: non-convergence, PMM LL < OLS LL, anti-correlation, extended-run instability
##
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/diagnose_convergence.jl

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

# ── Instrumented solvers returning (iters, final_diff) ───────────

function solveOLS_diag!(sa, r, w, cn, max_outer, rel_conv)
    @inbounds for col in 1:sa.n
        s = 0f0
        for i in sa.colptr[col]:(sa.colptr[col+1]-1)
            s += sa.nzval[i]^2
        end
        cn[col] = s
    end
    max_weight = 0f0
    iters = 0
    final_diff = 0f0
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
        final_diff = _diff
        _diff < rel_conv && break
    end
    return (iters, final_diff)
end

function solvePMM_diag!(sa, μ, y, w, max_outer, rel_conv)
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
    final_diff = 0f0
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
        final_diff = _diff
        _diff < rel_conv && break
    end
    # Unscale
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return (iters, final_diff)
end

function solveHuber_diag!(sa, r, w, δ_hub, λ_hub, nr_max, bs_max, max_outer, nr_acc, bs_acc, rel_conv, reg)
    max_weight = 0f0
    iters = 0
    final_diff = 0f0
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
        final_diff = _diff
        _diff < rel_conv && break
    end
    return (iters, final_diff)
end

"""
    extended_pmm_run!(sa, μ, y, w, extra_iters, rel_conv) → max_rel_weight_change

Starting from a converged PMM solution, run `extra_iters` more iterations.
Returns the maximum relative weight change observed, and mutates w/μ in place.
"""
function extended_pmm_run!(sa, μ, y, w, extra_iters, rel_conv)
    # Y-scaling (same as solvePMM_diag!)
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε_pmm = Float32(POISSON_MU_FLOOR)
    max_weight = maximum(w[1:sa.n])
    max_rel_change = 0f0
    for iter in 1:extra_iters
        weight_floor = max_weight * Float32(1e-4)
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
                    rc = δx / abs(w[col])
                    rc > max_rel_change && (max_rel_change = rc)
                end
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

# ── Load + build SparseArray helper ─────────────────────────────

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

# ── Load all problems ────────────────────────────────────────────

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems in $problem_dir\n")

# ── Warmup on first problem ─────────────────────────────────────

println("Warming up solvers...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n)
r0 = zeros(Float32, sa0.m); cn0 = zeros(Float32, sa0.n)
μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)

initResiduals_plain!(r0, sa0, w0)
solveOLS_diag!(sa0, r0, w0, cn0, Int64(d0[:max_iter_outer]), d0[:max_diff])

w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePMM_diag!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), Float32(0.001))

w0 .= 1f0; initResiduals_plain!(r0, sa0, w0)
solveHuber_diag!(sa0, r0, w0, d0[:huber_delta], d0[:lambda],
            Int64(d0[:max_iter_newton]), Int64(d0[:max_iter_bisection]),
            Int64(d0[:max_iter_outer]), Float32(1e-6), Float32(1e-6),
            d0[:max_diff], NoNorm())

# Warmup extended PMM
w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePMM_diag!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), Float32(0.001))
extended_pmm_run!(sa0, μ0, y0, w0, 10, Float32(0.001))
println("  Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct DiagResult
    scan_idx::Int
    n_cols::Int
    n_rows::Int
    # OLS convergence
    ols_iters::Int
    ols_final_diff::Float64
    ols_converged::Bool
    # PMM convergence
    pmm_iters::Int
    pmm_final_diff::Float64
    pmm_converged::Bool
    # Huber convergence
    hub_iters::Int
    hub_final_diff::Float64
    hub_converged::Bool
    # Poisson LL for all three
    ols_ll::Float64
    pmm_ll::Float64
    hub_ll::Float64
    # Extended PMM run
    ext_pmm_max_rel_change::Float64
    # Agreement
    cor_ols_hub::Float64
    cor_pmm_hub::Float64
    # Flags
    flag_pmm_not_converged::Bool
    flag_pmm_ll_lt_ols::Bool
    flag_pmm_anticorr::Bool
    flag_ext_pmm_unstable::Bool
end

results = DiagResult[]

# ── Run all problems ─────────────────────────────────────────────

const EXTRA_PMM_ITERS = 200
const EXT_UNSTABLE_THRESHOLD = 0.01  # max rel change > 1% = "unstable"
const PMM_REL_CONV = Float32(0.001)  # 10x tighter than default rel_conv (0.01)

println("=" ^ 160)
@printf("  %-8s %4s %5s │ %5s %10s %5s │ %5s %10s %5s │ %5s %10s %5s │ %14s %14s %14s │ %10s │ FLAGS\n",
        "Scan", "Cols", "Rows",
        "O_itr", "O_diff", "O_cv",
        "P_itr", "P_diff", "P_cv",
        "H_itr", "H_diff", "H_cv",
        "OLS_LL", "PMM_LL", "Hub_LL",
        "ext_Δw")
println("  " * "─" ^ 156)

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)

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

    # ── OLS ──
    r_ols  = zeros(Float32, sa.m)
    w_ols  = copy(cold_w)
    cn_ols = zeros(Float32, sa.n)
    initResiduals_plain!(r_ols, sa, w_ols)
    ols_iters, ols_diff = solveOLS_diag!(sa, r_ols, w_ols, cn_ols, max_outer_i64, rel_conv)
    ols_converged = ols_diff < rel_conv

    # ── PMM ──
    μ_pmm = zeros(Float32, sa.m)
    y_pmm = zeros(Float32, sa.m)
    w_pmm = copy(cold_w)
    initObserved!(y_pmm, sa)
    initMu!(μ_pmm, sa, w_pmm)
    pmm_iters, pmm_diff = solvePMM_diag!(sa, μ_pmm, y_pmm, w_pmm, max_outer_i64, PMM_REL_CONV)
    pmm_converged = pmm_diff < PMM_REL_CONV

    # ── Huber ──
    r_hub = zeros(Float32, sa.m)
    w_hub = copy(cold_w)
    initResiduals_plain!(r_hub, sa, w_hub)
    hub_iters, hub_diff = solveHuber_diag!(sa, r_hub, w_hub, δ_hub, λ_hub,
                    nr_max, bs_max, max_outer_i64, nr_acc, bs_acc, rel_conv, reg)
    hub_converged = hub_diff < rel_conv

    # ── Poisson LL for all three ──
    # OLS: build μ from OLS weights
    μ_ols = zeros(Float32, sa.m)
    initMu!(μ_ols, sa, w_ols)
    y_obs = zeros(Float32, sa.m)
    initObserved!(y_obs, sa)
    ols_ll = poissonLogLikelihood(μ_ols, y_obs, sa.m)

    # PMM: μ already computed
    pmm_ll = poissonLogLikelihood(μ_pmm, y_pmm, sa.m)

    # Huber: build μ from Huber weights
    μ_hub = zeros(Float32, sa.m)
    initMu!(μ_hub, sa, w_hub)
    hub_ll = poissonLogLikelihood(μ_hub, y_obs, sa.m)

    # ── Extended PMM run ──
    # Copy current PMM solution, then run 200 more iters
    w_pmm_ext = copy(w_pmm)
    μ_pmm_ext = copy(μ_pmm)
    y_pmm_ext = copy(y_pmm)
    ext_max_change = extended_pmm_run!(sa, μ_pmm_ext, y_pmm_ext, w_pmm_ext, EXTRA_PMM_ITERS, rel_conv)

    # ── Agreement ──
    wn = sa.n
    ov = Float64.(w_ols[1:wn])
    pv = Float64.(w_pmm[1:wn])
    hv = Float64.(w_hub[1:wn])
    c_oh = (all(iszero, ov) || all(iszero, hv)) ? NaN : cor(ov, hv)
    c_ph = (all(iszero, pv) || all(iszero, hv)) ? NaN : cor(pv, hv)

    # ── Flags ──
    flag_pmm_not_conv  = !pmm_converged
    flag_pmm_ll_lt_ols = pmm_ll < ols_ll
    flag_pmm_anticorr  = !isnan(c_ph) && c_ph < 0.0
    flag_ext_unstable  = ext_max_change > EXT_UNSTABLE_THRESHOLD

    push!(results, DiagResult(
        data[:scan_idx], sa.n, sa.m,
        ols_iters, Float64(ols_diff), ols_converged,
        pmm_iters, Float64(pmm_diff), pmm_converged,
        hub_iters, Float64(hub_diff), hub_converged,
        ols_ll, pmm_ll, hub_ll,
        Float64(ext_max_change),
        c_oh, c_ph,
        flag_pmm_not_conv, flag_pmm_ll_lt_ols, flag_pmm_anticorr, flag_ext_unstable
    ))

    # Build flag string
    flags = String[]
    flag_pmm_not_conv  && push!(flags, "PMM_NOCONV")
    flag_pmm_ll_lt_ols && push!(flags, "PMM_LL<OLS")
    flag_pmm_anticorr  && push!(flags, "PMM_ANTI")
    flag_ext_unstable  && push!(flags, "EXT_UNSTABLE")
    flag_str = isempty(flags) ? "" : join(flags, " ")

    @printf("  %-8d %4d %5d │ %5d %10.2e %5s │ %5d %10.2e %5s │ %5d %10.2e %5s │ %14.4e %14.4e %14.4e │ %10.2e │ %s\n",
            data[:scan_idx], sa.n, sa.m,
            ols_iters, ols_diff, ols_converged ? "  yes" : "** NO",
            pmm_iters, pmm_diff, pmm_converged ? "  yes" : "** NO",
            hub_iters, hub_diff, hub_converged ? "  yes" : "** NO",
            ols_ll, pmm_ll, hub_ll,
            ext_max_change,
            flag_str)
end

# ── Summary ──────────────────────────────────────────────────────

println("\n" * "=" ^ 160)
println("  CONVERGENCE DIAGNOSTIC SUMMARY — $(length(results)) problems")
println("=" ^ 160)

# Convergence status
n_ols_conv = count(r -> r.ols_converged, results)
n_pmm_conv = count(r -> r.pmm_converged, results)
n_hub_conv = count(r -> r.hub_converged, results)
println("\n  Convergence status:")
@printf("    OLS:   %d/%d converged\n", n_ols_conv, length(results))
@printf("    PMM:   %d/%d converged\n", n_pmm_conv, length(results))
@printf("    Huber: %d/%d converged\n", n_hub_conv, length(results))

# Iteration counts
println("\n  Iteration counts (converged problems only):")
for (label, get_iters, get_conv) in [
    ("OLS",   r -> r.ols_iters, r -> r.ols_converged),
    ("PMM",   r -> r.pmm_iters, r -> r.pmm_converged),
    ("Huber", r -> r.hub_iters, r -> r.hub_converged),
]
    conv = filter(get_conv, results)
    if !isempty(conv)
        its = [get_iters(r) for r in conv]
        @printf("    %-6s  mean=%.1f  median=%d  min=%d  max=%d\n",
                label, mean(its), Int(median(its)), minimum(its), maximum(its))
    end
end

# Final _diff values
println("\n  Final _diff values (all problems):")
for (label, get_diff) in [
    ("OLS",   r -> r.ols_final_diff),
    ("PMM",   r -> r.pmm_final_diff),
    ("Huber", r -> r.hub_final_diff),
]
    diffs = [get_diff(r) for r in results]
    @printf("    %-6s  mean=%.2e  median=%.2e  min=%.2e  max=%.2e\n",
            label, mean(diffs), median(diffs), minimum(diffs), maximum(diffs))
end

# Poisson LL comparison
println("\n  Poisson log-likelihood comparison:")
n_pmm_best_ll = count(r -> r.pmm_ll >= r.ols_ll && r.pmm_ll >= r.hub_ll, results)
n_ols_beats_pmm = count(r -> r.ols_ll > r.pmm_ll, results)
n_hub_beats_pmm = count(r -> r.hub_ll > r.pmm_ll, results)
@printf("    PMM has highest LL:     %d/%d\n", n_pmm_best_ll, length(results))
@printf("    OLS LL > PMM LL:        %d/%d  (convergence failure!)\n", n_ols_beats_pmm, length(results))
@printf("    Huber LL > PMM LL:      %d/%d  (convergence failure!)\n", n_hub_beats_pmm, length(results))

# LL differences where OLS > PMM
ols_beats = filter(r -> r.ols_ll > r.pmm_ll, results)
if !isempty(ols_beats)
    ll_gaps = [r.ols_ll - r.pmm_ll for r in ols_beats]
    @printf("    LL gap (OLS-PMM) when OLS wins:  mean=%.4e  max=%.4e\n",
            mean(ll_gaps), maximum(ll_gaps))
end

# Extended PMM stability
println("\n  Extended PMM run (+$(EXTRA_PMM_ITERS) iters from converged solution):")
ext_changes = [r.ext_pmm_max_rel_change for r in results]
@printf("    Max relative weight change:  mean=%.2e  median=%.2e  max=%.2e\n",
        mean(ext_changes), median(ext_changes), maximum(ext_changes))
n_ext_unstable = count(r -> r.flag_ext_pmm_unstable, results)
@printf("    Unstable (change > %.0e):     %d/%d\n", EXT_UNSTABLE_THRESHOLD, n_ext_unstable, length(results))

# Agreement
println("\n  PMM-Huber weight correlation:")
cors_ph = filter(!isnan, [r.cor_pmm_hub for r in results])
if !isempty(cors_ph)
    @printf("    mean=%.4f  median=%.4f  min=%.4f  max=%.4f\n",
            mean(cors_ph), median(cors_ph), minimum(cors_ph), maximum(cors_ph))
end
n_anti = count(r -> r.flag_pmm_anticorr, results)
@printf("    Anti-correlated (cor<0):  %d/%d\n", n_anti, length(results))

# Flag summary
println("\n  FLAG SUMMARY:")
n_f1 = count(r -> r.flag_pmm_not_converged, results)
n_f2 = count(r -> r.flag_pmm_ll_lt_ols, results)
n_f3 = count(r -> r.flag_pmm_anticorr, results)
n_f4 = count(r -> r.flag_ext_pmm_unstable, results)
n_any = count(r -> r.flag_pmm_not_converged || r.flag_pmm_ll_lt_ols || r.flag_pmm_anticorr || r.flag_ext_pmm_unstable, results)
@printf("    PMM_NOCONV    (didn't converge):     %d\n", n_f1)
@printf("    PMM_LL<OLS    (OLS has higher LL):   %d\n", n_f2)
@printf("    PMM_ANTI      (cor(PMM,Hub) < 0):    %d\n", n_f3)
@printf("    EXT_UNSTABLE  (extended run moved):  %d\n", n_f4)
@printf("    ANY FLAG:                            %d/%d\n", n_any, length(results))

# Detail flagged problems
flagged = filter(r -> r.flag_pmm_not_converged || r.flag_pmm_ll_lt_ols || r.flag_pmm_anticorr || r.flag_ext_pmm_unstable, results)
if !isempty(flagged)
    println("\n  FLAGGED PROBLEMS (detail):")
    println("  " * "─" ^ 120)
    @printf("  %-8s %4s │ %5s %10s %5s │ %14s %14s %14s │ %10s │ FLAGS\n",
            "Scan", "Cols",
            "P_itr", "P_diff", "P_cv",
            "OLS_LL", "PMM_LL", "Hub_LL",
            "ext_Δw")
    println("  " * "─" ^ 120)
    for r in flagged
        flags = String[]
        r.flag_pmm_not_converged && push!(flags, "PMM_NOCONV")
        r.flag_pmm_ll_lt_ols    && push!(flags, "PMM_LL<OLS")
        r.flag_pmm_anticorr     && push!(flags, "PMM_ANTI")
        r.flag_ext_pmm_unstable && push!(flags, "EXT_UNSTABLE")
        @printf("  %-8d %4d │ %5d %10.2e %5s │ %14.4e %14.4e %14.4e │ %10.2e │ %s\n",
                r.scan_idx, r.n_cols,
                r.pmm_iters, r.pmm_final_diff, r.pmm_converged ? "  yes" : "** NO",
                r.ols_ll, r.pmm_ll, r.hub_ll,
                r.ext_pmm_max_rel_change,
                join(flags, " "))
    end
end

# Clean problems (no flags)
n_clean = length(results) - length(flagged)
println("\n  Clean problems (no flags): $n_clean/$(length(results))")

println("\n" * "=" ^ 160)
println("  Done.")
println("=" ^ 160)
