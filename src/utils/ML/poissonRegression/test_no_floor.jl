## Quick test: what happens with weight_floor = 0 (no floor)?
## Does any problem fail to converge?

using Printf, Statistics, Serialization
using Pioneer

include("SparseArray.jl")
include("spectralPoissonRegression.jl")

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

function solvePMM_nofloor!(sa, μ, y, w, max_outer, rel_conv, max_inner)
    y_scale = maximum(y[1:sa.m])
    if y_scale > 1f0
        @inbounds for i in 1:sa.m; y[i] /= y_scale; end
        @inbounds for j in 1:sa.n; w[j] /= y_scale; end
        initMu!(μ, sa, w)
    end
    ε = Float32(POISSON_MU_FLOOR)
    outer_iters = 0

    for iter in 1:max_outer
        _diff = 0f0

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
            # NO floor — every nonzero weight participates in convergence check
            if !iszero(w[col])
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
    return outer_iters
end

const K = 5
const REL_CONV = Float32(0.001)

problem_dir = "/Users/n.t.wamsley/Desktop/solveHuber_problems"
files = sort(filter(f -> endswith(f, ".jls"), readdir(problem_dir)))

# Warmup
d0 = deserialize(joinpath(problem_dir, files[1]))
sa0 = load_sa(d0)
w0 = ones(Float32, sa0.n); μ0 = zeros(Float32, sa0.m); y0 = zeros(Float32, sa0.m)
initObserved!(y0, sa0); initMu!(μ0, sa0, w0)
solvePMM_nofloor!(sa0, μ0, y0, w0, Int64(d0[:max_iter_outer]), REL_CONV, K)

iters_nofloor = Int[]
iters_floor = Int[]
hit_max = Int[]

for fname in files
    data = deserialize(joinpath(problem_dir, fname))
    sa = load_sa(data)
    max_outer = Int64(data[:max_iter_outer])

    # No floor
    w1 = ones(Float32, sa.n); μ1 = zeros(Float32, sa.m); y1 = zeros(Float32, sa.m)
    initObserved!(y1, sa); initMu!(μ1, sa, w1)
    oi1 = solvePMM_nofloor!(sa, μ1, y1, w1, max_outer, REL_CONV, K)
    push!(iters_nofloor, oi1)
    if oi1 >= max_outer
        push!(hit_max, data[:scan_idx])
    end

    # With floor (1e-4)
    w2 = ones(Float32, sa.n); μ2 = zeros(Float32, sa.m); y2 = zeros(Float32, sa.m)
    initObserved!(y2, sa); initMu!(μ2, sa, w2)
    solvePoissonMM!(sa, μ2, y2, w2, max_outer, REL_CONV; max_inner_iter=Int64(K))
    # Count iters by checking if we hit max
    # Actually, solvePoissonMM! doesn't return iters. Just compare weights.
    push!(iters_floor, 0)  # placeholder
end

println("=" ^ 80)
println("  NO FLOOR TEST — $(length(files)) problems, K=$K, rel_conv=$REL_CONV")
println("=" ^ 80)

@printf("\n  No-floor iters: mean=%.1f  median=%d  max=%d\n",
        mean(iters_nofloor), round(Int, median(iters_nofloor)), maximum(iters_nofloor))
@printf("  Hit max_iter: %d / %d\n", length(hit_max), length(files))

if !isempty(hit_max)
    println("\n  Problems that hit max_iter (FAILED TO CONVERGE):")
    for s in hit_max
        @printf("    scan=%d\n", s)
    end
end

println("\n" * "=" ^ 80)
