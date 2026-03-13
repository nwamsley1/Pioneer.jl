## Benchmark: weight_floor = 1e-2 vs 1e-3 vs 1e-4
##
## For each floor value, runs K=5 PMM with rel_conv=0.001 and also
## a high-iter trace (200 iters, no early stop) to verify convergence.
##
## Run:  julia --project=. src/utils/ML/poissonRegression/benchmark_weight_floor.jl

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

# ── PMM solver parameterized by weight_floor ─────────────────────

function solvePMM_floor!(sa, μ, y, w, max_outer, rel_conv, max_inner, wf_frac::Float32)
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

    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * wf_frac : 0f0
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
            δx = abs(w[col] - X_before)
            w[col] > max_weight && (max_weight = w[col])
            if w[col] > weight_floor
                rc = δx / max(abs(w[col]), Float32(1e-10))
                rc > _diff && (_diff = rc)
            end
        end
        outer_iters = iter
        _diff < rel_conv && break
    end

    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return (outer_iters, total_inner)
end

# ── High-iter trace: convergence curve without early stopping ────
# Runs inside a SINGLE y-scaling context (no double-scaling artifact)

function solvePMM_trace!(sa, μ, y, w, max_outer, max_inner, wf_frac::Float32)
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    max_weight = 0f0

    trace_diff = Float64[]

    for iter in 1:max_outer
        _diff = 0f0
        weight_floor = iter > 5 ? max_weight * wf_frac : 0f0
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
        push!(trace_diff, Float64(_diff))
    end

    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] *= y_scale; end
        @inbounds for j in 1:sa.n; w[j] *= y_scale; end
        initMu!(μ, sa, w)
    end
    return trace_diff
end

# ── Config ───────────────────────────────────────────────────────

const K = 5
const REL_CONV = Float32(0.001)
const FLOOR_VALUES = [Float32(1e-2), Float32(1e-3), Float32(1e-4)]
const TRACE_ITERS = 200  # for convergence trace (no early stop)

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Found $(length(files)) problems\n")

# Warmup
println("Warming up...")
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
for wf in FLOOR_VALUES
    w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
    solvePMM_floor!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), REL_CONV, K, wf)
    w0 .= 1f0; initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
    solvePMM_trace!(sa0, μ0, y0, w0, 20, K, wf)
end
println("Done.\n")

# ── Results storage ──────────────────────────────────────────────

struct FloorResult
    scan_idx::Int
    n_cols::Int
    wf::Float32
    outer_iters::Int
    total_inner::Int
    time_s::Float64
    final_ll::Float64
    # From trace: iter where diff first < 0.001, and final diff at iter 200
    trace_conv_iter::Int
    trace_final_diff::Float64
end

all_results = Dict{Float32, Vector{FloorResult}}()
for wf in FLOOR_VALUES
    all_results[wf] = FloorResult[]
end

# ── Run all problems ─────────────────────────────────────────────

println("Running $(length(files)) problems × $(length(FLOOR_VALUES)) floor values...\n")

for (fi, fname) in enumerate(files)
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    for wf in FLOOR_VALUES
        # Timed solve with early stopping
        w_k = ones(Float32, sa.n)
        μ_k = zeros(Float32, sa.m); y_k = zeros(Float32, sa.m)
        initObserved!(y_k, sa); initMu!(μ_k, sa, w_k)

        local oitr, nitr
        t = @elapsed begin
            oitr, nitr = solvePMM_floor!(sa, μ_k, y_k, w_k, max_outer, REL_CONV, K, wf)
        end
        ll = poissonLogLikelihood(μ_k, y_k, sa.m)

        # Trace: 200 iters, no early stop, single y-scaling context
        w_t = ones(Float32, sa.n)
        μ_t = zeros(Float32, sa.m); y_t = zeros(Float32, sa.m)
        initObserved!(y_t, sa); initMu!(μ_t, sa, w_t)
        trace = solvePMM_trace!(sa, μ_t, y_t, w_t, TRACE_ITERS, K, wf)

        # Find first iter where diff < rel_conv
        trace_conv = -1
        for i in 1:length(trace)
            if trace[i] < Float64(REL_CONV)
                trace_conv = i
                break
            end
        end

        push!(all_results[wf], FloorResult(
            data[:scan_idx], sa.n, wf, oitr, nitr, t, ll,
            trace_conv, trace[end]
        ))
    end

    if fi % 25 == 0
        println("  ... $fi / $(length(files))")
    end
end

println("\n  All problems done.\n")

# ── Summary table ────────────────────────────────────────────────

println("=" ^ 130)
println("  WEIGHT FLOOR BENCHMARK — $(length(files)) problems, K=$K, rel_conv=$REL_CONV")
println("=" ^ 130)

println("\n  Per-floor Summary (with early stopping at rel_conv=$REL_CONV):")
@printf("  %10s │ %6s %6s %6s │ %8s %8s │ %10s %10s │ %14s\n",
        "wf_frac",
        "Oi_m", "Oi_md", "Oi_mx",
        "Ni_m", "Ni_md",
        "Time_m(s)", "Time_tot",
        "LL_median")
println("  " * "─" ^ 100)

for wf in FLOOR_VALUES
    res = all_results[wf]
    oi = [r.outer_iters for r in res]
    ni = [r.total_inner for r in res]
    ts = [r.time_s for r in res]
    lls = [r.final_ll for r in res]
    @printf("  %10.0e │ %6.1f %6d %6d │ %8.0f %8d │ %10.6f %10.4f │ %14.4e\n",
            wf,
            mean(oi), round(Int, median(oi)), maximum(oi),
            mean(ni), round(Int, median(ni)),
            mean(ts), sum(ts),
            median(lls))
end

# ── Trace-based convergence (no early stopping) ─────────────────

println("\n  Convergence from 200-iter trace (no early stopping):")
@printf("  %10s │ %8s %8s %8s │ %8s │ %12s %12s\n",
        "wf_frac",
        "conv_m", "conv_md", "conv_mx",
        "#noconv",
        "fin_diff_md", "fin_diff_mx")
println("  " * "─" ^ 85)

for wf in FLOOR_VALUES
    res = all_results[wf]
    conv_iters = [r.trace_conv_iter for r in res if r.trace_conv_iter > 0]
    n_noconv = count(r -> r.trace_conv_iter < 0, res)
    fin_diffs = [r.trace_final_diff for r in res]

    if !isempty(conv_iters)
        @printf("  %10.0e │ %8.1f %8d %8d │ %8d │ %12.2e %12.2e\n",
                wf,
                mean(conv_iters), round(Int, median(conv_iters)), maximum(conv_iters),
                n_noconv,
                median(fin_diffs), maximum(fin_diffs))
    else
        @printf("  %10.0e │ %8s %8s %8s │ %8d │ %12.2e %12.2e\n",
                wf, "n/a", "n/a", "n/a", n_noconv,
                median(fin_diffs), maximum(fin_diffs))
    end
end

# ── LL comparison across floors ──────────────────────────────────

println("\n  LL agreement across floor values:")
ref_lls = [r.final_ll for r in all_results[FLOOR_VALUES[1]]]
for wf in FLOOR_VALUES[2:end]
    lls = [r.final_ll for r in all_results[wf]]
    rel_diffs = abs.(lls .- ref_lls) ./ (abs.(ref_lls) .+ 1.0)
    @printf("    floor=%.0e vs %.0e:  median_rel_diff=%.2e  max=%.2e\n",
            wf, FLOOR_VALUES[1], median(rel_diffs), maximum(rel_diffs))
end

# ── Per-problem detail: any disagreements? ───────────────────────

println("\n  Problems where iter count differs across floors:")
@printf("  %-8s %4s │", "Scan", "Cols")
for wf in FLOOR_VALUES
    @printf("  wf=%.0e", wf)
end
println("  │ notes")
println("  " * "─" ^ 70)

n_differ = 0
for pi in 1:length(files)
    iters = [all_results[wf][pi].outer_iters for wf in FLOOR_VALUES]
    if minimum(iters) != maximum(iters)
        global n_differ += 1
        r = all_results[FLOOR_VALUES[1]][pi]
        @printf("  %-8d %4d │", r.scan_idx, r.n_cols)
        for wf in FLOOR_VALUES
            @printf("  %8d", all_results[wf][pi].outer_iters)
        end
        # Check if LL differs meaningfully
        lls = [all_results[wf][pi].final_ll for wf in FLOOR_VALUES]
        ll_spread = (maximum(lls) - minimum(lls)) / (abs(maximum(lls)) + 1.0)
        @printf("  │ LL_spread=%.2e", ll_spread)
        println()
        n_differ >= 30 && (println("  ... (truncated)"); break)
    end
end
println("  Total: $n_differ / $(length(files)) problems with different iter counts")

# ── Trace detail for a few interesting problems ──────────────────

println("\n  Convergence trace comparison (first 15 iters) for 3 sample problems:")
sample_indices = [1, 50, 100]
for si in sample_indices
    si > length(files) && continue
    r = all_results[FLOOR_VALUES[1]][si]
    println("\n  ── Scan $(r.scan_idx) │ $(r.n_cols) cols ──")
    @printf("  %4s", "iter")
    for wf in FLOOR_VALUES
        @printf("  wf=%.0e", wf)
    end
    println()
    println("  " * "─" ^ 40)

    # Re-run traces for this problem to show side by side
    data = deserialize(joinpath(problem_dir, files[si]))
    sa = load_sa(data)
    traces = Dict{Float32, Vector{Float64}}()
    for wf in FLOOR_VALUES
        w_t = ones(Float32, sa.n)
        μ_t = zeros(Float32, sa.m); y_t = zeros(Float32, sa.m)
        initObserved!(y_t, sa); initMu!(μ_t, sa, w_t)
        traces[wf] = solvePMM_trace!(sa, μ_t, y_t, w_t, 30, K, wf)
    end

    for i in 1:min(30, TRACE_ITERS)
        @printf("  %4d", i)
        for wf in FLOOR_VALUES
            @printf("  %10.2e", traces[wf][i])
        end
        # Mark convergence
        markers = String[]
        for wf in FLOOR_VALUES
            if i > 1 && traces[wf][i] < Float64(REL_CONV) && traces[wf][i-1] >= Float64(REL_CONV)
                push!(markers, @sprintf("←%.0e", wf))
            end
        end
        !isempty(markers) && print("  ", join(markers, " "))
        println()
    end
end

println("\n" * "=" ^ 130)
println("  Done.")
println("=" ^ 130)
