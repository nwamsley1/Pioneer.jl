## Test Y-scaled OLS vs Original OLS on a FEW real serialized problems
##
## Run:  julia --project=. src/utils/ML/poissonRegression/test_yscaled_real.jl

using Pioneer, Serialization, Printf

include(joinpath(@__DIR__, "SparseArray.jl"))
include(joinpath(@__DIR__, "solveOLS.jl"))

function sse_f64(r::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += Float64(r[i])^2; end; s
end

function compute_sse(sa, w)
    r = zeros(Float32, sa.m)
    @inbounds for n in 1:sa.n_vals
        if iszero(r[sa.rowval[n]]); r[sa.rowval[n]] = -sa.x[n]; end
    end
    @inbounds for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1]-1)
            r[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
    sse_f64(r, sa.m)
end

function load_sa(path)
    data = deserialize(path)
    rowval = Vector{Int64}(data[:rowval])
    n_vals = data[:n_vals]

    # data[:x] is per-row (length m), but SparseArray.x is per-nonzero (length n_vals)
    # Expand: x_nz[k] = x_row[rowval[k]]
    x_row = data[:x]
    if length(x_row) == data[:m]
        x_nz = Vector{Float32}(undef, n_vals)
        for k in 1:n_vals
            x_nz[k] = x_row[rowval[k]]
        end
    else
        x_nz = x_row  # already per-nonzero
    end

    sa = Main.SparseArray(
        n_vals, data[:m], data[:n],
        rowval, Vector{UInt16}(data[:colval]), data[:nzval],
        haskey(data, :matched) ? data[:matched] : ones(Bool, n_vals),
        haskey(data, :isotope) ? data[:isotope] : zeros(UInt8, n_vals),
        x_nz, Vector{Int64}(data[:colptr])
    )
    max_outer = Int64(data[:max_iter_outer])
    rel_conv  = haskey(data, :max_diff) ? data[:max_diff] : Float32(0.01)
    return sa, max_outer, rel_conv
end

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"

# Run ALL problems
test_files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))
println("Testing $(length(test_files)) problems\n")

# Warmup on first problem
sa0, mo0, rc0 = load_sa(joinpath(problem_dir, test_files[1]))
r0 = zeros(Float32, sa0.m); w0 = ones(Float32, sa0.n); cn0 = zeros(Float32, sa0.n)
initResiduals!(r0, sa0, w0); solveOLS!(sa0, r0, w0, cn0, mo0, rc0)
w0 .= 1f0; initResiduals!(r0, sa0, w0); solveOLS_yscaled!(sa0, r0, w0, cn0, mo0, rc0)
sa0 = nothing; GC.gc()
println("Warmup done.\n")

println("=" ^ 130)
println("  OLS vs Y-SCALED OLS on real problems (cold start w=ones)")
println("=" ^ 130)
@printf("\n  %-36s %4s %5s │ %12s %12s │ %12s %12s │ %7s\n",
        "Problem", "Cols", "Rows",
        "OLS_SSE", "OLS_maxW",
        "Yscl_SSE", "Yscl_maxW",
        "Status")
println("  " * "─" ^ 120)

n_ols_nan = 0; n_yscl_nan = 0; n_fixed = 0

for fname in test_files
    global n_ols_nan, n_yscl_nan, n_fixed
    sa, max_outer, rel_conv = load_sa(joinpath(problem_dir, fname))

    # Data stats
    max_x = maximum(sa.x[i] for i in 1:sa.n_vals)
    max_nz = maximum(sa.nzval[i] for i in 1:sa.n_vals)
    min_nz = minimum(sa.nzval[i] for i in 1:sa.n_vals if sa.nzval[i] > 0)

    # OLS (cold start w=ones)
    r = zeros(Float32, sa.m); w = ones(Float32, sa.n); cn = zeros(Float32, sa.n)
    initResiduals!(r, sa, w)
    solveOLS!(sa, r, w, cn, max_outer, rel_conv)
    ols_sse = compute_sse(sa, w)
    ols_maxw = maximum(abs.(w[1:sa.n]))

    # Y-scaled OLS (cold start w=ones)
    r2 = zeros(Float32, sa.m); w2 = ones(Float32, sa.n); cn2 = zeros(Float32, sa.n)
    initResiduals!(r2, sa, w2)
    solveOLS_yscaled!(sa, r2, w2, cn2, max_outer, rel_conv)
    ys_sse = compute_sse(sa, w2)
    ys_maxw = maximum(abs.(w2[1:sa.n]))

    o_bad = isnan(ols_sse) || isinf(ols_sse)
    y_bad = isnan(ys_sse) || isinf(ys_sse)
    o_bad && (n_ols_nan += 1)
    y_bad && (n_yscl_nan += 1)
    if o_bad && !y_bad; n_fixed += 1; end

    status = if o_bad && !y_bad; "FIXED"
    elseif o_bad && y_bad; "BOTH"
    elseif !o_bad && !y_bad; "OK"
    else; "REGRESS"
    end

    @printf("  %-36s %4d %5d │ %12.4e %12.4e │ %12.4e %12.4e │ %7s\n",
            fname, sa.n, sa.m, ols_sse, ols_maxw, ys_sse, ys_maxw, status)
    @printf("    (max_x=%.2e  nzval=[%.2e, %.2e])\n", max_x, min_nz, max_nz)

    # Free memory
    sa = nothing; r = nothing; w = nothing; r2 = nothing; w2 = nothing
    GC.gc(false)
end

println("\n  " * "─" ^ 80)
@printf("  NaN/Inf SSE:  OLS=%d   Y-scaled=%d   (of %d)\n", n_ols_nan, n_yscl_nan, length(test_files))
@printf("  FIXED:        %d problems (NaN→finite with y-scaling)\n", n_fixed)
println("=" ^ 130)
