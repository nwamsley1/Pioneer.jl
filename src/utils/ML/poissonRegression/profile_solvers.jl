## Profile OLS v1, PMM v1, and Huber on small/medium/large problems
## Generates per-solver, per-problem flame graphs (.pb.gz)
##
## Problems:
##   Small:  scan 88280  (9 cols,  98 rows)
##   Medium: scan 97068  (35 cols, 399 rows)
##   Large:  scan 342335 (449 cols, 4566 rows)
##
## Run:  julia --project=<Pioneer root> src/utils/ML/poissonRegression/profile_solvers.jl

using Printf, Statistics, Serialization, BenchmarkTools, Profile, PProf

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

# ── Load problems ────────────────────────────────────────────────

problems = [
    ("small_88280",  "/Users/n.t.wamsley/Desktop/solveHuber_problems/problem_88280_9cols.jls"),
    ("med_97068",    "/Users/n.t.wamsley/Desktop/solveHuber_problems/problem_97068_35cols.jls"),
    ("large_342335", "/Users/n.t.wamsley/Desktop/solveHuber_inputs_scan342335.jls"),
]

profile_dir = joinpath(@__DIR__, "..", "..", "..", "..", "data")

# ── Load and build SparseArrays ──────────────────────────────────

function load_problem(path::String)
    data = deserialize(path)
    # The two file formats have slightly different field types
    rowval = haskey(data, :rowval) ? Vector{Int64}(data[:rowval]) : error("no rowval")
    colptr = haskey(data, :colptr) ? Vector{Int64}(data[:colptr]) : error("no colptr")
    colval = Vector{UInt16}(data[:colval])
    n_vals = data[:n_vals]

    x_row = data[:x]
    if length(x_row) == data[:m]
        x_nz = Vector{Float32}(undef, n_vals)
        for k in 1:n_vals; x_nz[k] = x_row[rowval[k]]; end
    else
        x_nz = x_row
    end

    sa = Main.SparseArray(
        n_vals, data[:m], data[:n],
        rowval, colval, data[:nzval],
        haskey(data, :matched) ? data[:matched] : ones(Bool, n_vals),
        haskey(data, :isotope) ? data[:isotope] : zeros(UInt8, n_vals),
        x_nz, colptr
    )

    max_outer = Int64(data[:max_iter_outer])
    rel_conv  = haskey(data, :max_diff) ? data[:max_diff] : Float32(0.01)

    # Huber params
    δ_hub  = haskey(data, :huber_delta) ? data[:huber_delta] : (haskey(data, :delta) ? data[:delta] : Float32(1e6))
    λ_hub  = haskey(data, :lambda) ? data[:lambda] : Float32(0)
    nr_max = Int64(haskey(data, :max_iter_newton) ? data[:max_iter_newton] : 50)
    bs_max = Int64(haskey(data, :max_iter_bisection) ? data[:max_iter_bisection] : 100)
    nr_acc = haskey(data, :accuracy_newton) ? data[:accuracy_newton] : Float32(1e-6)
    bs_acc = haskey(data, :accuracy_bisection) ? data[:accuracy_bisection] : Float32(1e-6)

    return sa, max_outer, rel_conv, δ_hub, λ_hub, nr_max, bs_max, nr_acc, bs_acc
end

# ── Warmup ───────────────────────────────────────────────────────

println("Loading problems and warming up...")
problem_data = Dict{String, Any}()

for (label, path) in problems
    sa, max_outer, rel_conv, δ_hub, λ_hub, nr_max, bs_max, nr_acc, bs_acc = load_problem(path)
    problem_data[label] = (sa=sa, max_outer=max_outer, rel_conv=rel_conv,
                           δ_hub=δ_hub, λ_hub=λ_hub, nr_max=nr_max, bs_max=bs_max,
                           nr_acc=nr_acc, bs_acc=bs_acc)
    println("  $label: m=$(sa.m), n=$(sa.n), n_vals=$(sa.n_vals)")

    # Warmup all three solvers
    cold_w = ones(Float32, sa.n)

    r = zeros(Float32, sa.m); cn = zeros(Float32, sa.n); w = copy(cold_w)
    initResiduals_plain!(r, sa, w)
    solveOLS!(sa, r, w, cn, max_outer, rel_conv)

    μ = zeros(Float32, sa.m); y = zeros(Float32, sa.m); w = copy(cold_w)
    initObserved!(y, sa); initMu!(μ, sa, w)
    solvePoissonMM!(sa, μ, y, w, max_outer, rel_conv)

    r = zeros(Float32, sa.m); w = copy(cold_w)
    initResiduals_plain!(r, sa, w)
    solveHuber!(sa, r, w, δ_hub, λ_hub, nr_max, bs_max, max_outer, nr_acc, bs_acc, rel_conv, NoNorm())
end
println("  Done.\n")

# ══════════════════════════════════════════════════════════════════
# For each problem: benchmark timing, then profile each solver
# ══════════════════════════════════════════════════════════════════

for (label, _path) in problems
    pd = problem_data[label]
    sa = pd.sa
    max_outer = pd.max_outer
    rel_conv = pd.rel_conv
    cold_w = ones(Float32, sa.n)

    println("=" ^ 100)
    println("  PROBLEM: $label  (m=$(sa.m), n=$(sa.n), n_vals=$(sa.n_vals))")
    println("=" ^ 100)

    # Pre-allocate buffers
    _r   = Vector{Float32}(undef, sa.m)
    _w   = Vector{Float32}(undef, sa.n)
    _cn2 = Vector{Float32}(undef, sa.n)
    _μ   = Vector{Float32}(undef, sa.m)
    _y   = Vector{Float32}(undef, sa.m)

    # ── 1. Benchmark timing ──────────────────────────────────────

    println("\n  1. TIMING (@benchmark, cold-start, evals=1)")
    println("  " * "─" ^ 60)

    b_ols = @benchmark solveOLS!($sa, $_r, $_w, $_cn2, $max_outer, $rel_conv) setup=begin
        copyto!($_w, $cold_w); initResiduals_plain!($_r, $sa, $_w)
    end evals=1

    b_pmm = @benchmark solvePoissonMM!($sa, $_μ, $_y, $_w, $max_outer, $rel_conv) setup=begin
        copyto!($_w, $cold_w); initObserved!($_y, $sa); initMu!($_μ, $sa, $_w)
    end evals=1

    b_hub = @benchmark solveHuber!($sa, $_r, $_w, $(pd.δ_hub), $(pd.λ_hub), $(pd.nr_max), $(pd.bs_max),
                                    $max_outer, $(pd.nr_acc), $(pd.bs_acc), $rel_conv, $(NoNorm())) setup=begin
        copyto!($_w, $cold_w); initResiduals_plain!($_r, $sa, $_w)
    end evals=1

    t_ols = median(b_ols).time / 1e9
    t_pmm = median(b_pmm).time / 1e9
    t_hub = median(b_hub).time / 1e9

    @printf("  OLS v1:  %.6f s  (n=%d samples)\n", t_ols, length(b_ols.times))
    @printf("  PMM v1:  %.6f s  (n=%d samples)\n", t_pmm, length(b_pmm.times))
    @printf("  Huber:   %.6f s  (n=%d samples)\n", t_hub, length(b_hub.times))
    @printf("  Ratios:  Hub/OLS=%.2fx  Hub/PMM=%.2fx  PMM/OLS=%.2fx\n",
            t_hub/t_ols, t_hub/t_pmm, t_pmm/t_ols)

    # ── 2. Solution comparison ───────────────────────────────────

    println("\n  2. SOLUTION COMPARISON (cold-start)")
    println("  " * "─" ^ 60)

    # OLS
    r_ols = zeros(Float32, sa.m); w_ols = copy(cold_w); cn_ols = zeros(Float32, sa.n)
    initResiduals_plain!(r_ols, sa, w_ols)
    solveOLS!(sa, r_ols, w_ols, cn_ols, max_outer, rel_conv)
    sse_ols = sse_from_residuals(r_ols, sa.m)

    # PMM
    μ_pmm = zeros(Float32, sa.m); y_pmm = zeros(Float32, sa.m); w_pmm = copy(cold_w)
    initObserved!(y_pmm, sa); initMu!(μ_pmm, sa, w_pmm)
    solvePoissonMM!(sa, μ_pmm, y_pmm, w_pmm, max_outer, rel_conv)
    r_pmm = zeros(Float32, sa.m); initResiduals_plain!(r_pmm, sa, w_pmm)
    sse_pmm = sse_from_residuals(r_pmm, sa.m)
    ll_pmm = poissonLogLikelihood(μ_pmm, y_pmm, sa.m)

    # Huber
    r_hub = zeros(Float32, sa.m); w_hub = copy(cold_w)
    initResiduals_plain!(r_hub, sa, w_hub)
    solveHuber!(sa, r_hub, w_hub, pd.δ_hub, pd.λ_hub, pd.nr_max, pd.bs_max,
                max_outer, pd.nr_acc, pd.bs_acc, rel_conv, NoNorm())
    sse_hub = sse_from_residuals(r_hub, sa.m)

    @printf("  OLS SSE:   %16.6e\n", sse_ols)
    @printf("  PMM SSE:   %16.6e   LL: %16.6e\n", sse_pmm, ll_pmm)
    @printf("  Huber SSE: %16.6e\n", sse_hub)

    ov = Float64.(w_ols[1:sa.n]); pv = Float64.(w_pmm[1:sa.n]); hv = Float64.(w_hub[1:sa.n])
    if !all(iszero, ov) && !all(iszero, hv)
        @printf("  cor(OLS,Hub)=%.6f  cor(PMM,Hub)=%.6f  cor(OLS,PMM)=%.6f\n",
                cor(ov, hv), cor(pv, hv), cor(ov, pv))
    end
    @printf("  max|Δw| OLS-Hub=%.4e  PMM-Hub=%.4e\n",
            maximum(abs.(ov .- hv)), maximum(abs.(pv .- hv)))

    # Weight sparsity
    for (slabel, ww) in [("OLS", w_ols), ("PMM", w_pmm), ("Huber", w_hub)]
        w_slice = ww[1:sa.n]
        mw = maximum(w_slice)
        n_zero = count(iszero, w_slice)
        @printf("  %-6s  max_w=%.4e  zeros=%d/%d (%.1f%%)\n",
                slabel, mw, n_zero, sa.n, 100.0*n_zero/sa.n)
    end

    # ── 3. CPU Profiling ─────────────────────────────────────────

    # Scale profile runs inversely with problem size
    n_profile = sa.n < 20 ? 5000 : (sa.n < 100 ? 2000 : 1000)

    println("\n  3. CPU PROFILING ($n_profile runs each, cold-start)")
    println("  " * "─" ^ 60)

    # --- OLS v1 ---
    print("  Profiling OLS v1...")
    Profile.clear()
    @profile for _ in 1:n_profile
        copyto!(_w, cold_w); initResiduals_plain!(_r, sa, _w)
        solveOLS!(sa, _r, _w, _cn2, max_outer, rel_conv)
    end
    outfile = joinpath(profile_dir, "profile_$(label)_ols.pb.gz")
    pprof(out=outfile, web=false)
    println(" → $(basename(outfile))")

    # --- PMM v1 ---
    print("  Profiling PMM v1...")
    Profile.clear()
    @profile for _ in 1:n_profile
        copyto!(_w, cold_w); initObserved!(_y, sa); initMu!(_μ, sa, _w)
        solvePoissonMM!(sa, _μ, _y, _w, max_outer, rel_conv)
    end
    outfile = joinpath(profile_dir, "profile_$(label)_pmm.pb.gz")
    pprof(out=outfile, web=false)
    println(" → $(basename(outfile))")

    # --- Huber ---
    print("  Profiling Huber...")
    Profile.clear()
    @profile for _ in 1:n_profile
        copyto!(_w, cold_w); initResiduals_plain!(_r, sa, _w)
        solveHuber!(sa, _r, _w, pd.δ_hub, pd.λ_hub, pd.nr_max, pd.bs_max,
                    max_outer, pd.nr_acc, pd.bs_acc, rel_conv, NoNorm())
    end
    outfile = joinpath(profile_dir, "profile_$(label)_huber.pb.gz")
    pprof(out=outfile, web=false)
    println(" → $(basename(outfile))")

    println()
end

println("=" ^ 100)
println("  Done. Profile files written to: $profile_dir")
println("=" ^ 100)
