## Diagnose OLS "false convergence" on cold-start problems
##
## Runs OLS with per-iteration diagnostics on all problems, identifies cases
## where OLS converges in ≤5 iterations, and compares with Huber to determine
## whether the solution is actually correct.
##
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/diagnose_ols_convergence.jl

using Printf, Statistics, Serialization

# Need Pioneer's types for deserialization
using Pioneer

# Local SparseArray + all solvers
include("SparseArray.jl")
include("spectralLinearRegression_reference.jl")
include("spectralPoissonRegression.jl")

# ── Helpers ──────────────────────────────────────────────────────

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

function sse_from_residuals(r::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += Float64(r[i])^2; end; s
end

# ── Instrumented OLS that logs per-iteration diagnostics ────────

struct OLSIterLog
    iter::Int
    max_diff::Float32        # max relative weight change (_diff)
    max_weight::Float32      # largest weight after this sweep
    n_zero::Int              # number of zero weights
    n_clamped_to_zero::Int   # number that went from >0 to 0 this sweep
    n_unclamped::Int         # number that went from 0 to >0 this sweep
    sse::Float64             # SSE after this sweep
    max_abs_change::Float32  # max absolute weight change
    sum_abs_change::Float64  # sum of all absolute weight changes
end

function solveOLS_instrumented!(
    sa::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, colnorm2::Vector{T},
    max_iter_outer::Int64, rel_conv::T
) where {Ti<:Integer, T<:AbstractFloat}

    @inbounds for col in 1:sa.n
        s = zero(T)
        for i in sa.colptr[col]:(sa.colptr[col+1]-1)
            s += sa.nzval[i] * sa.nzval[i]
        end
        colnorm2[col] = s
    end

    logs = OLSIterLog[]
    max_weight = T(0)
    prev_w = copy(w)

    for iter in 1:max_iter_outer
        _diff = T(0)
        weight_floor = iter > 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)
        max_abs = T(0)
        sum_abs = 0.0
        n_clamped = 0
        n_unclamped = 0

        for col in 1:sa.n
            L2 = colnorm2[col]; iszero(L2) && continue

            L1 = zero(T)
            @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1)
                L1 += sa.nzval[k] * r[sa.rowval[k]]
            end

            X0 = w[col]
            w[col] = max(w[col] - L1 / L2, zero(T))
            updateResiduals!(sa, r, col, w[col], X0)

            δx = abs(w[col] - X0)
            δx > max_abs && (max_abs = δx)
            sum_abs += Float64(δx)

            # Track clamping events
            if X0 > zero(T) && w[col] == zero(T)
                n_clamped += 1
            elseif X0 == zero(T) && w[col] > zero(T)
                n_unclamped += 1
            end

            if w[col] > max_weight
                max_weight = w[col]
            end
            if w[col] > weight_floor
                rc = δx / abs(w[col])
                rc > _diff && (_diff = rc)
            end
        end

        n_zero = count(iszero, @view(w[1:sa.n]))
        sse = sse_from_residuals(r, sa.m)
        push!(logs, OLSIterLog(iter, _diff, max_weight, n_zero, n_clamped, n_unclamped, sse, max_abs, sum_abs))

        _diff < rel_conv && break
    end

    return logs
end

# ── Load problem ────────────────────────────────────────────────

function load_problem(path::String)
    data = deserialize(path)
    rowval = Vector{Int64}(data[:rowval])
    colptr = Vector{Int64}(data[:colptr])
    colval = Vector{UInt16}(data[:colval])
    n_vals = data[:n_vals]

    sa = Main.SparseArray(
        n_vals, data[:m], data[:n],
        rowval, colval, data[:nzval],
        haskey(data, :matched) ? data[:matched] : ones(Bool, n_vals),
        haskey(data, :isotope) ? data[:isotope] : zeros(UInt8, n_vals),
        data[:x], colptr
    )

    max_outer = Int64(data[:max_iter_outer])
    rel_conv  = haskey(data, :max_diff) ? data[:max_diff] : Float32(0.01)

    δ_hub  = haskey(data, :huber_delta) ? data[:huber_delta] : (haskey(data, :delta) ? data[:delta] : Float32(1e6))
    λ_hub  = haskey(data, :lambda) ? data[:lambda] : Float32(0)
    nr_max = Int64(haskey(data, :max_iter_newton) ? data[:max_iter_newton] : 50)
    bs_max = Int64(haskey(data, :max_iter_bisection) ? data[:max_iter_bisection] : 100)
    nr_acc = haskey(data, :accuracy_newton) ? data[:accuracy_newton] :
             (haskey(data, :convergence_tol) ? data[:convergence_tol] : Float32(1e-6))
    bs_acc = haskey(data, :accuracy_bisection) ? data[:accuracy_bisection] : nr_acc

    return sa, max_outer, rel_conv, δ_hub, λ_hub, nr_max, bs_max, nr_acc, bs_acc
end

# ── Run warmup ──────────────────────────────────────────────────

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

println("Warming up...")
sa0, mo0, rc0, dh0, lh0, nr0, bs0, na0, ba0 = load_problem(joinpath(problem_dir, files[1]))
r0 = zeros(Float32, sa0.m); w0 = ones(Float32, sa0.n); cn0 = zeros(Float32, sa0.n)
initResiduals_plain!(r0, sa0, w0)
solveOLS!(sa0, r0, w0, cn0, mo0, rc0)
r0 = zeros(Float32, sa0.m); w0 = ones(Float32, sa0.n)
initResiduals_plain!(r0, sa0, w0)
solveHuber!(sa0, r0, w0, dh0, lh0, nr0, bs0, mo0, na0, ba0, rc0, NoNorm())
_ = solveOLS_instrumented!(sa0, zeros(Float32, sa0.m), ones(Float32, sa0.n), cn0, mo0, rc0)
println("  Done.\n")

# ══════════════════════════════════════════════════════════════════
# Phase 1: Run all problems, identify early-convergence cases
# ══════════════════════════════════════════════════════════════════

println("=" ^ 120)
println("  PHASE 1: IDENTIFY EARLY OLS CONVERGENCE (≤5 iterations)")
println("=" ^ 120)

struct ProblemSummary
    fname::String
    n_cols::Int
    n_rows::Int
    rel_conv::Float32
    ols_iters::Int
    ols_sse::Float64
    hub_iters::Int
    hub_sse::Float64
    cor_oh::Float64
    max_diff_oh::Float64
    ols_max_weight::Float32
    hub_max_weight::Float32
end

summaries = ProblemSummary[]
early_convergence_data = Dict{String, Vector{OLSIterLog}}()

for fname in files
    sa, max_outer, rel_conv, δ_hub, λ_hub, nr_max, bs_max, nr_acc, bs_acc = load_problem(joinpath(problem_dir, fname))
    cold_w = ones(Float32, sa.n)

    # Instrumented OLS
    r_ols = zeros(Float32, sa.m); w_ols = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r_ols, sa, w_ols)
    logs = solveOLS_instrumented!(sa, r_ols, w_ols, cn, max_outer, rel_conv)
    ols_sse = sse_from_residuals(r_ols, sa.m)

    # Huber for comparison (instrumented iteration count)
    r_hub = zeros(Float32, sa.m); w_hub = copy(cold_w)
    initResiduals_plain!(r_hub, sa, w_hub)
    hub_iters = 0
    max_weight_h = Float32(0)
    for iter in 1:max_outer
        _diff = Float32(0)
        wf = iter > 5 ? max_weight_h * Float32(1e-7) : Float32(0)
        max_weight_h = Float32(0)
        for col in 1:sa.n
            δx = abs(newton_bisection!(sa, r_hub, w_hub, col, δ_hub, λ_hub,
                        nr_max, bs_max, nr_acc, bs_acc, NoNorm(), rel_conv))
            w_hub[col] > max_weight_h && (max_weight_h = w_hub[col])
            if w_hub[col] > wf
                rc = δx / abs(w_hub[col]); rc > _diff && (_diff = rc)
            end
        end
        hub_iters = iter
        _diff < rel_conv && break
    end
    hub_sse = sse_from_residuals(r_hub, sa.m)

    # Agreement
    ov = Float64.(w_ols[1:sa.n]); hv = Float64.(w_hub[1:sa.n])
    c_oh = (all(iszero, ov) || all(iszero, hv) || length(ov) < 2) ? NaN : cor(ov, hv)
    md_oh = maximum(abs.(ov .- hv))

    push!(summaries, ProblemSummary(
        fname, sa.n, sa.m, rel_conv,
        length(logs), ols_sse,
        hub_iters, hub_sse,
        c_oh, md_oh,
        isempty(logs) ? 0f0 : logs[end].max_weight,
        max_weight_h
    ))

    if length(logs) <= 5
        early_convergence_data[fname] = logs
    end
end

# ── Print summary table ─────────────────────────────────────────

# Sort by OLS iterations ascending
sorted_idx = sortperm(summaries, by=s -> s.ols_iters)

println()
@printf("  %-36s %4s %5s %8s │ %5s %12s %8s │ %5s %12s %8s │ %8s %10s\n",
        "Problem", "Cols", "Rows", "rel_conv",
        "O_itr", "OLS_SSE", "O_maxW",
        "H_itr", "Hub_SSE", "H_maxW",
        "cor(O,H)", "maxΔw")
println("  " * "─" ^ 130)

n_early = 0
for i in sorted_idx
    global n_early
    s = summaries[i]
    marker = s.ols_iters <= 5 ? ">>>" : "   "
    s.ols_iters <= 5 && (n_early += 1)
    @printf("%s %-36s %4d %5d %8.4f │ %5d %12.4e %8.2f │ %5d %12.4e %8.2f │ %8.4f %10.4e\n",
            marker, s.fname, s.n_cols, s.n_rows, s.rel_conv,
            s.ols_iters, s.ols_sse, s.ols_max_weight,
            s.hub_iters, s.hub_sse, s.hub_max_weight,
            isnan(s.cor_oh) ? -1.0 : s.cor_oh, s.max_diff_oh)
end

println("\n  $n_early / $(length(summaries)) problems converged in ≤5 iterations (marked >>>)")

# ══════════════════════════════════════════════════════════════════
# Phase 2: Detailed per-iteration analysis of early-convergence cases
# ══════════════════════════════════════════════════════════════════

if n_early > 0
    println("\n" * "=" ^ 120)
    println("  PHASE 2: PER-ITERATION DIAGNOSTICS FOR EARLY-CONVERGENCE PROBLEMS")
    println("=" ^ 120)

    for s in sort(summaries[sorted_idx], by=s -> s.ols_iters)
        s.ols_iters > 5 && break
        logs = early_convergence_data[s.fname]

        println("\n  ┌─ $(s.fname)  ($(s.n_cols) cols, $(s.n_rows) rows, rel_conv=$(s.rel_conv))")
        println("  │")
        @printf("  │  %-5s │ %12s %10s │ %6s %8s %8s │ %12s %12s │ %12s\n",
                "Iter", "max_rel_Δ", "max_w",
                "zeros", "→0", "→>0",
                "SSE", "max_abs_Δ", "Σ|Δw|")
        println("  │  " * "─" ^ 100)

        for log in logs
            converged_marker = log.max_diff < s.rel_conv ? " ← CONVERGED" : ""
            @printf("  │  %5d │ %12.6e %10.4f │ %6d %8d %8d │ %12.4e %12.4e │ %12.4e%s\n",
                    log.iter, log.max_diff, log.max_weight,
                    log.n_zero, log.n_clamped_to_zero, log.n_unclamped,
                    log.sse, log.max_abs_change, log.sum_abs_change,
                    converged_marker)
        end

        # Show what the converged weights look like vs Huber
        println("  │")
        @printf("  │  OLS solution: SSE=%.4e  max_weight=%.4f  zeros=%d/%d\n",
                s.ols_sse, s.ols_max_weight, logs[end].n_zero, s.n_cols)
        @printf("  │  Huber soln:   SSE=%.4e  max_weight=%.4f\n", s.hub_sse, s.hub_max_weight)
        @printf("  │  cor(OLS,Hub)=%.6f  max|Δw|=%.4e\n",
                isnan(s.cor_oh) ? -1.0 : s.cor_oh, s.max_diff_oh)

        # Analyze the convergence criterion failure mode
        if length(logs) > 0 && logs[end].max_diff < s.rel_conv
            last = logs[end]
            n_significant = s.n_cols - last.n_zero
            println("  │")
            println("  │  DIAGNOSIS:")
            @printf("  │    Columns checked for convergence: %d of %d (%.1f%%)\n",
                    n_significant, s.n_cols, 100.0 * n_significant / s.n_cols)
            @printf("  │    Columns at zero (excluded):      %d of %d (%.1f%%)\n",
                    last.n_zero, s.n_cols, 100.0 * last.n_zero / s.n_cols)
            if last.sum_abs_change > 0.01
                println("  │    ⚠ Large total absolute change ($(Printf.@sprintf("%.4e", last.sum_abs_change))) despite \"convergence\"")
            end
            if s.ols_sse > 2 * s.hub_sse && s.hub_sse > 0
                @printf("  │    ⚠ OLS SSE is %.1fx worse than Huber SSE\n", s.ols_sse / s.hub_sse)
            end
            if !isnan(s.cor_oh) && s.cor_oh < 0.95
                @printf("  │    ⚠ Low correlation with Huber (%.4f) — solutions disagree\n", s.cor_oh)
            end
        end
        println("  └─")
    end
end

# ══════════════════════════════════════════════════════════════════
# Phase 3: The fix — run OLS with forced minimum iterations
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 120)
println("  PHASE 3: COMPARE FIXES ON EARLY-CONVERGENCE PROBLEMS")
println("=" ^ 120)

# Fix 1: Force minimum iterations (don't check convergence before iter N)
function solveOLS_miniter!(
    sa::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, cn::Vector{T},
    max_iter::Int64, rel_conv::T, min_iter::Int64
) where {Ti<:Integer, T<:AbstractFloat}

    @inbounds for col in 1:sa.n
        s = zero(T); for i in sa.colptr[col]:(sa.colptr[col+1]-1); s += sa.nzval[i]^2; end
        cn[col] = s
    end

    max_weight = T(0); iters = 0
    for iter in 1:max_iter
        _diff = T(0)
        wf = iter > 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)
        for col in 1:sa.n
            L2 = cn[col]; iszero(L2) && continue
            L1 = zero(T)
            @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1)
                L1 += sa.nzval[k] * r[sa.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1/L2, zero(T))
            updateResiduals!(sa, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > wf; rc = δx/abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        iter >= min_iter && _diff < rel_conv && break
    end
    return iters
end

# Fix 2: Include zero→nonzero transitions in convergence check
function solveOLS_countall!(
    sa::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, cn::Vector{T},
    max_iter::Int64, rel_conv::T
) where {Ti<:Integer, T<:AbstractFloat}

    @inbounds for col in 1:sa.n
        s = zero(T); for i in sa.colptr[col]:(sa.colptr[col+1]-1); s += sa.nzval[i]^2; end
        cn[col] = s
    end

    max_weight = T(0); iters = 0
    for iter in 1:max_iter
        _diff = T(0)
        wf = iter > 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)
        for col in 1:sa.n
            L2 = cn[col]; iszero(L2) && continue
            L1 = zero(T)
            @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1)
                L1 += sa.nzval[k] * r[sa.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1/L2, zero(T))
            updateResiduals!(sa, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])

            # KEY CHANGE: also count columns going TO or FROM zero
            if w[col] > wf || (X0 > zero(T) && w[col] == zero(T))
                if w[col] > zero(T)
                    rc = δx / abs(w[col])
                elseif X0 > zero(T)
                    # Weight went to zero — use X0 as denominator
                    rc = δx / abs(X0)
                else
                    rc = zero(T)
                end
                rc > _diff && (_diff = rc)
            end
        end
        iters = iter
        _diff < rel_conv && break
    end
    return iters
end

# Fix 3: Use absolute convergence check based on SSE change
function solveOLS_ssestop!(
    sa::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, cn::Vector{T},
    max_iter::Int64, rel_conv::T
) where {Ti<:Integer, T<:AbstractFloat}

    @inbounds for col in 1:sa.n
        s = zero(T); for i in sa.colptr[col]:(sa.colptr[col+1]-1); s += sa.nzval[i]^2; end
        cn[col] = s
    end

    max_weight = T(0); iters = 0; prev_sse = Inf
    for iter in 1:max_iter
        _diff = T(0)
        wf = iter > 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)
        for col in 1:sa.n
            L2 = cn[col]; iszero(L2) && continue
            L1 = zero(T)
            @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1)
                L1 += sa.nzval[k] * r[sa.rowval[k]]
            end
            X0 = w[col]; w[col] = max(w[col] - L1/L2, zero(T))
            updateResiduals!(sa, r, col, w[col], X0)
            δx = abs(w[col] - X0)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > wf; rc = δx/abs(w[col]); rc > _diff && (_diff = rc); end
        end
        iters = iter
        sse = sse_from_residuals(r, sa.m)
        # Require BOTH: relative weight convergence AND SSE has plateaued
        sse_rel_change = abs(sse - prev_sse) / max(sse, 1.0)
        if _diff < rel_conv && sse_rel_change < Float64(rel_conv)
            break
        end
        prev_sse = sse
    end
    return iters
end

# Run fixes only on the early-convergence problems
early_problems = filter(s -> s.ols_iters <= 5, summaries)

if !isempty(early_problems)
    println()
    @printf("  %-36s %4s │ %5s %12s │ %5s %12s │ %5s %12s │ %5s %12s │ %5s %12s │ %8s\n",
            "Problem", "Cols",
            "orig", "orig_SSE",
            "min5", "min5_SSE",
            "min10", "min10_SSE",
            "ctall", "ctall_SSE",
            "sseS", "sseS_SSE",
            "Hub_SSE")
    println("  " * "─" ^ 150)

    for s in sort(early_problems, by=s -> s.n_cols)
        sa, max_outer, rel_conv, δh, λh, nr, bs, nra, bsa = load_problem(joinpath(problem_dir, s.fname))
        cold_w = ones(Float32, sa.n)

        # Original
        r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
        initResiduals_plain!(r, sa, w)
        orig_iters = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(1))
        orig_sse = sse_from_residuals(r, sa.m)

        # Fix 1a: min 5 iterations
        r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
        initResiduals_plain!(r, sa, w)
        min5_iters = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(5))
        min5_sse = sse_from_residuals(r, sa.m)

        # Fix 1b: min 10 iterations
        r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
        initResiduals_plain!(r, sa, w)
        min10_iters = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(10))
        min10_sse = sse_from_residuals(r, sa.m)

        # Fix 2: count all transitions
        r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
        initResiduals_plain!(r, sa, w)
        ctall_iters = solveOLS_countall!(sa, r, w, cn, max_outer, rel_conv)
        ctall_sse = sse_from_residuals(r, sa.m)

        # Fix 3: SSE-based stopping
        r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
        initResiduals_plain!(r, sa, w)
        sses_iters = solveOLS_ssestop!(sa, r, w, cn, max_outer, rel_conv)
        sses_sse = sse_from_residuals(r, sa.m)

        @printf("  %-36s %4d │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e │ %12.4e\n",
                s.fname, s.n_cols,
                orig_iters, orig_sse,
                min5_iters, min5_sse,
                min10_iters, min10_sse,
                ctall_iters, ctall_sse,
                sses_iters, sses_sse,
                s.hub_sse)
    end
else
    println("\n  No early-convergence problems found — all problems converged in >5 iterations.")
end

# ══════════════════════════════════════════════════════════════════
# Phase 4: Verify fixes don't hurt well-behaved problems
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 120)
println("  PHASE 4: VERIFY FIXES ON WELL-BEHAVED PROBLEMS (≤30 cols)")
println("=" ^ 120)

small_problems = filter(s -> s.n_cols <= 30, summaries)
println()
@printf("  %-36s %4s │ %5s %12s │ %5s %12s │ %5s %12s │ %5s %12s │ %5s %12s\n",
        "Problem", "Cols",
        "orig", "orig_SSE",
        "min5", "min5_SSE",
        "min10", "min10_SSE",
        "ctall", "ctall_SSE",
        "sseS", "sseS_SSE")
println("  " * "─" ^ 150)

orig_total_iters = 0; min5_total = 0; min10_total = 0; ctall_total = 0; sses_total = 0
for s in sort(small_problems, by=s -> s.n_cols)[1:min(20, length(small_problems))]
    global orig_total_iters, min5_total, min10_total, ctall_total, sses_total
    sa, max_outer, rel_conv, _, _, _, _, _, _ = load_problem(joinpath(problem_dir, s.fname))
    cold_w = ones(Float32, sa.n)

    r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r, sa, w)
    oi = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(1))
    os = sse_from_residuals(r, sa.m); orig_total_iters += oi

    r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r, sa, w)
    m5i = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(5))
    m5s = sse_from_residuals(r, sa.m); min5_total += m5i

    r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r, sa, w)
    m10i = solveOLS_miniter!(sa, r, w, cn, max_outer, rel_conv, Int64(10))
    m10s = sse_from_residuals(r, sa.m); min10_total += m10i

    r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r, sa, w)
    cai = solveOLS_countall!(sa, r, w, cn, max_outer, rel_conv)
    cas = sse_from_residuals(r, sa.m); ctall_total += cai

    r = zeros(Float32, sa.m); w = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals_plain!(r, sa, w)
    ssi = solveOLS_ssestop!(sa, r, w, cn, max_outer, rel_conv)
    sss = sse_from_residuals(r, sa.m); sses_total += ssi

    @printf("  %-36s %4d │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e │ %5d %12.4e\n",
            s.fname, s.n_cols,
            oi, os, m5i, m5s, m10i, m10s, cai, cas, ssi, sss)
end

n_shown = min(20, length(small_problems))
@printf("\n  Total iterations (first %d small problems):\n", n_shown)
@printf("    Original:   %d\n", orig_total_iters)
@printf("    Min-5:      %d\n", min5_total)
@printf("    Min-10:     %d\n", min10_total)
@printf("    Count-all:  %d\n", ctall_total)
@printf("    SSE-stop:   %d\n", sses_total)

# ══════════════════════════════════════════════════════════════════
# Phase 5: Y-scaled OLS on ALL 100 problems
# ══════════════════════════════════════════════════════════════════

println("\n" * "=" ^ 120)
println("  PHASE 5: Y-SCALED OLS vs ORIGINAL OLS vs HUBER (all $(length(files)) problems)")
println("=" ^ 120)

include("solveOLS.jl")  # brings in solveOLS_yscaled!

# Warmup
sa0, mo0, rc0, _, _, _, _, _, _ = load_problem(joinpath(problem_dir, files[1]))
r0 = zeros(Float32, sa0.m); w0 = ones(Float32, sa0.n); cn0 = zeros(Float32, sa0.n)
initResiduals!(r0, sa0, w0)
solveOLS_yscaled!(sa0, r0, w0, cn0, mo0, rc0)

println()
@printf("  %-36s %4s │ %5s %12s %10s │ %5s %12s %10s │ %5s %12s %10s │ %8s %8s\n",
        "Problem", "Cols",
        "O_itr", "OLS_SSE", "O_maxW",
        "Ys_it", "Yscl_SSE", "Ys_maxW",
        "H_itr", "Hub_SSE", "H_maxW",
        "cor(Ys,H)", "cor(O,H)")
println("  " * "─" ^ 155)

n_ols_nan = 0; n_yscl_nan = 0; n_hub_nan = 0
n_yscl_better = 0; n_ols_better = 0; n_equal = 0
ols_total_t = 0.0; yscl_total_t = 0.0; hub_total_t = 0.0

for fname in files
    global n_ols_nan, n_yscl_nan, n_hub_nan, n_yscl_better, n_ols_better, n_equal
    global ols_total_t, yscl_total_t, hub_total_t
    sa, max_outer, rel_conv, δh, λh, nr, bs, nra, bsa = load_problem(joinpath(problem_dir, fname))
    cold_w = ones(Float32, sa.n)

    # Original OLS
    r_ols = zeros(Float32, sa.m); w_ols = copy(cold_w); cn = zeros(Float32, sa.n)
    initResiduals!(r_ols, sa, w_ols)
    ols_iters = 0; ols_t = @elapsed begin
        @inbounds for col in 1:sa.n
            s = 0f0; for i in sa.colptr[col]:(sa.colptr[col+1]-1); s += sa.nzval[i]^2; end; cn[col] = s
        end
        mw = 0f0
        for iter in 1:max_outer
            _d = 0f0; wf = iter > 5 ? mw * 1f-7 : 0f0; mw = 0f0
            for col in 1:sa.n
                L2 = cn[col]; iszero(L2) && continue
                L1 = 0f0
                @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1); L1 += sa.nzval[k]*r_ols[sa.rowval[k]]; end
                X0 = w_ols[col]; w_ols[col] = max(w_ols[col]-L1/L2, 0f0)
                updateResiduals!(sa, r_ols, col, w_ols[col], X0)
                dx = abs(w_ols[col]-X0); w_ols[col]>mw && (mw=w_ols[col])
                if w_ols[col]>wf; rc=dx/abs(w_ols[col]); rc>_d && (_d=rc); end
            end
            ols_iters = iter; _d < rel_conv && break
        end
    end
    ols_sse = sse_from_residuals(r_ols, sa.m)
    ols_total_t += ols_t

    # Y-scaled OLS
    r_ys = zeros(Float32, sa.m); w_ys = copy(cold_w); cn_ys = zeros(Float32, sa.n)
    initResiduals!(r_ys, sa, w_ys)
    ys_iters = 0; ys_t = @elapsed begin
        # Inline y-scaling to count iterations
        ys = 0f0
        @inbounds for n in 1:sa.n_vals; sa.x[n] > ys && (ys = sa.x[n]); end
        if ys > 1f0
            inv_s = 1f0/ys
            @inbounds for n in 1:sa.n_vals; sa.x[n] *= inv_s; end
            @inbounds for j in 1:sa.n; w_ys[j] *= inv_s; end
            initResiduals!(r_ys, sa, w_ys)
        end
        @inbounds for col in 1:sa.n
            s = 0f0; for i in sa.colptr[col]:(sa.colptr[col+1]-1); s += sa.nzval[i]^2; end; cn_ys[col] = s
        end
        mw = 0f0
        for iter in 1:max_outer
            _d = 0f0; wf = iter > 5 ? mw * 1f-7 : 0f0; mw = 0f0
            for col in 1:sa.n
                L2 = cn_ys[col]; iszero(L2) && continue
                L1 = 0f0
                @inbounds @fastmath for k in sa.colptr[col]:(sa.colptr[col+1]-1); L1 += sa.nzval[k]*r_ys[sa.rowval[k]]; end
                X0 = w_ys[col]; w_ys[col] = max(w_ys[col]-L1/L2, 0f0)
                updateResiduals!(sa, r_ys, col, w_ys[col], X0)
                dx = abs(w_ys[col]-X0); w_ys[col]>mw && (mw=w_ys[col])
                if w_ys[col]>wf; rc=dx/abs(w_ys[col]); rc>_d && (_d=rc); end
            end
            ys_iters = iter; _d < rel_conv && break
        end
        # Unscale
        if ys > 1f0
            @inbounds for n in 1:sa.n_vals; sa.x[n] *= ys; end
            @inbounds for j in 1:sa.n; w_ys[j] *= ys; end
        end
    end
    r_ys2 = zeros(Float32, sa.m); initResiduals_plain!(r_ys2, sa, w_ys)
    ys_sse = sse_from_residuals(r_ys2, sa.m)
    yscl_total_t += ys_t

    # Huber
    r_hub = zeros(Float32, sa.m); w_hub = copy(cold_w)
    initResiduals_plain!(r_hub, sa, w_hub)
    hub_iters = 0; hub_t = @elapsed begin
        mwh = 0f0
        for iter in 1:max_outer
            _d = 0f0; wf = iter > 5 ? mwh * 1f-7 : 0f0; mwh = 0f0
            for col in 1:sa.n
                dx = abs(newton_bisection!(sa, r_hub, w_hub, col, δh, λh, nr, bs, nra, bsa, NoNorm(), rel_conv))
                w_hub[col]>mwh && (mwh=w_hub[col])
                if w_hub[col]>wf; rc=dx/abs(w_hub[col]); rc>_d && (_d=rc); end
            end
            hub_iters = iter; _d < rel_conv && break
        end
    end
    hub_sse = sse_from_residuals(r_hub, sa.m)
    hub_total_t += hub_t

    # Agreement
    wn = sa.n
    ov = Float64.(w_ols[1:wn]); yv = Float64.(w_ys[1:wn]); hv = Float64.(w_hub[1:wn])
    c_yh = (all(iszero, yv) || all(iszero, hv) || length(yv) < 2) ? NaN : cor(yv, hv)
    c_oh = (all(iszero, ov) || all(iszero, hv) || length(ov) < 2) ? NaN : cor(ov, hv)

    isnan(ols_sse) && (n_ols_nan += 1)
    isnan(ys_sse)  && (n_yscl_nan += 1)
    isnan(hub_sse) && (n_hub_nan += 1)

    if !isnan(ys_sse) && !isnan(ols_sse)
        if ys_sse < ols_sse * 0.99; n_yscl_better += 1
        elseif ols_sse < ys_sse * 0.99; n_ols_better += 1
        else; n_equal += 1; end
    elseif !isnan(ys_sse) && isnan(ols_sse)
        n_yscl_better += 1
    elseif isnan(ys_sse) && !isnan(ols_sse)
        n_ols_better += 1
    end

    marker = (isnan(ols_sse) && !isnan(ys_sse)) ? "FIX" :
             (isnan(ols_sse) && isnan(ys_sse))   ? "NaN" : "   "
    @printf("%s %-36s %4d │ %5d %12.4e %10.2f │ %5d %12.4e %10.2f │ %5d %12.4e %10.2f │ %8.4f %8.4f\n",
            marker, fname, sa.n,
            ols_iters, ols_sse, maximum(abs.(w_ols[1:wn])),
            ys_iters, ys_sse, maximum(abs.(w_ys[1:wn])),
            hub_iters, hub_sse, maximum(abs.(w_hub[1:wn])),
            isnan(c_yh) ? -1.0 : c_yh,
            isnan(c_oh) ? -1.0 : c_oh)
end

println("\n  " * "─" ^ 80)
println("  SUMMARY:")
@printf("    NaN SSE:     OLS=%d   Y-scaled=%d   Huber=%d   (of %d)\n",
        n_ols_nan, n_yscl_nan, n_hub_nan, length(files))
@printf("    SSE winner:  Y-scaled better=%d   OLS better=%d   Equal=%d\n",
        n_yscl_better, n_ols_better, n_equal)
@printf("    Total time:  OLS=%.4f s   Y-scaled=%.4f s   Huber=%.4f s\n",
        ols_total_t, yscl_total_t, hub_total_t)

println("\n" * "=" ^ 120)
println("  Done.")
println("=" ^ 120)
