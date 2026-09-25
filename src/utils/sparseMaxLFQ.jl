using Random: MersenneTwister, shuffle!, rand
using LinearAlgebra: dot, norm
using Statistics: median!, mean

const SPARSE_MAXLFQ_PARTNERS = 16
const SPARSE_MAXLFQ_SEED = 0

function _lfq_root!(parents, i)
    while parents[i] != i
        parents[i] = parents[parents[i]]
        i = parents[i]
    end
    return i
end

function _lfq_union!(parents, sizes, i, j)
    a, b = _lfq_root!(parents, i), _lfq_root!(parents, j)
    a == b && return false
    sizes[a] < sizes[b] && ((a, b) = (b, a))
    parents[b] = a
    sizes[a] += sizes[b]
    return true
end

# A spanning forest preserves the complete overlap graph's connectivity.
# Additional partners are sampled through shared precursors, without computing
# an all-pairs overlap matrix or using intensities/experimental conditions.
function _sparse_lfq_pairs(X, partners, seed)
    np, nr = size(X)
    parents, sizes = collect(1:nr), ones(Int, nr)
    run_precursors = [Int[] for _ in 1:nr]
    precursor_runs = [Int[] for _ in 1:np]
    pairs = Set{Tuple{Int,Int}}()
    rng = MersenneTwister(seed)
    for p in 1:np
        runs = precursor_runs[p]
        for r in 1:nr
            v = X[p, r]
            ismissing(v) && continue
            isfinite(v) || throw(ArgumentError("Sparse MaxLFQ expects finite log intensities or missing"))
            push!(runs, r)
            push!(run_precursors[r], p)
        end
        shuffle!(rng, runs)
        for i in 2:length(runs)
            a, b = runs[i-1], runs[i]
            _lfq_union!(parents, sizes, a, b) && push!(pairs, minmax(a, b))
        end
    end
    if partners >= nr - 1
        for a in 1:nr, b in a+1:nr
            any(p -> !ismissing(X[p, b]), run_precursors[a]) && push!(pairs, (a,b))
        end
    elseif partners > 0
        for a in 1:nr
            precursors = run_precursors[a]
            isempty(precursors) && continue
            chosen = Set{Int}()
            for _ in 1:8partners
                runs = precursor_runs[rand(rng, precursors)]
                b = rand(rng, runs)
                b == a && continue
                push!(chosen, b)
                push!(pairs, minmax(a,b))
                length(chosen) >= partners && break
            end
        end
    end
    labels = zeros(Int, nr)
    root_labels = Dict{Int,Int}()
    for r in 1:nr
        root = _lfq_root!(parents, r)
        labels[r] = get!(root_labels, root) do
            length(root_labels) + 1
        end
    end
    return sort!(collect(pairs)), labels, run_precursors
end

function _lfq_laplacian_mul!(out, x, pairs)
    fill!(out, 0)
    for (a,b) in pairs
        delta = x[a] - x[b]
        out[a] += delta
        out[b] -= delta
    end
    return out
end

# Projected Jacobi-preconditioned CG solves on the zero-mean subspace.
function _lfq_sparse_solve(pairs, rhs, degree; rtol, atol, maxiter)
    n = length(rhs)
    x = zeros(Float64, n)
    residual = copy(rhs)
    residual .-= mean(residual)
    tolerance = max(atol, rtol * norm(rhs))
    norm(residual) <= tolerance && return x, 0, norm(residual)
    z = residual ./ degree
    z .-= mean(z)
    direction = copy(z)
    product = similar(x)
    rz = dot(residual, z)
    for iteration in 1:maxiter
        _lfq_laplacian_mul!(product, direction, pairs)
        denom = dot(direction, product)
        denom > 0 || error("Sparse MaxLFQ solver lost positive curvature")
        alpha = rz / denom
        x .+= alpha .* direction
        residual .-= alpha .* product
        residual .-= mean(residual)
        if norm(residual) <= tolerance
            # Check the actual equation residual before reporting convergence.
            _lfq_laplacian_mul!(product, x, pairs)
            residual .= rhs .- product
            residual .-= mean(residual)
            if norm(residual) <= tolerance
                x .-= mean(x)
                return x, iteration, norm(residual)
            end
            z .= residual ./ degree
            z .-= mean(z)
            direction .= z
            rz = dot(residual, z)
            continue
        end
        z .= residual ./ degree
        z .-= mean(z)
        next_rz = dot(residual, z)
        direction .= z .+ (next_rz / rz) .* direction
        rz = next_rz
    end
    error("Sparse MaxLFQ did not converge within $maxiter iterations")
end

"""
    solve_sparse_maxlfq(X, run_priorities; partners=16, seed=0,
                       rtol=1e-10, atol=1e-12, maxiter=max(100, 4size(X,2)))

MaxLFQ using a sparse set of same-precursor median log ratios.
Input is a precursor-by-run matrix of finite log2 intensities or missing values.
Keeps the full overlap graph's best connected component and MaxLFQ intensity
scaling. Returns a named tuple with estimates, component_labels, edge_count,
iterations, and residual_norm. Single-run handling remains the caller's job,
as with solve_maxlfq.

A randomized spanning forest guarantees connectivity. Each run proposes at most
partners additional distinct neighbors through shared precursors, with bounded
attempts. The total edge count is at most (partners+1)*number_of_runs; partners
is not a maximum degree. Setting partners >= number_of_runs-1 uses all valid
pairs for reference comparisons. Fixed seeds are reproducible for fixed input
ordering. Graph selection depends on observation availability, not intensities.

The input matrix and its observation indexes require O(precursors*runs) and
O(observations) space respectively; no dense run-by-run matrix is allocated.
Each iterative solve step is O(edges). Convergence is checked, never silently
truncated; strict linear total runtime is not guaranteed.
"""
function solve_sparse_maxlfq(X::AbstractMatrix, run_priorities::AbstractVector;
    partners::Int=SPARSE_MAXLFQ_PARTNERS, seed::Int=SPARSE_MAXLFQ_SEED, rtol::Real=1e-10, atol::Real=1e-12,
    maxiter::Int=max(100, 4size(X,2)))
    partners >= 0 || throw(ArgumentError("partners must be nonnegative"))
    rtol > 0 && isfinite(rtol) || throw(ArgumentError("rtol must be finite and positive"))
    atol >= 0 && isfinite(atol) || throw(ArgumentError("atol must be finite and nonnegative"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    nr = size(X,2)
    length(run_priorities) == nr || throw(DimensionMismatch("Priorities must match run columns"))
    estimates = Vector{Union{Missing,Float32}}(missing, nr)
    pairs, labels, run_precursors = _sparse_lfq_pairs(X, partners, seed)
    nr == 0 && return (; estimates, component_labels=labels, edge_count=0, iterations=0, residual_norm=0.0)
    counts = zeros(Int, maximum(labels))
    priorities = zeros(Float64, length(counts))
    for r in 1:nr
        counts[labels[r]] += 1
        v = run_priorities[r]
        ismissing(v) || (priorities[labels[r]] += v)
    end
    best = argmax([(counts[i], priorities[i]) for i in eachindex(counts)])
    indices = findall(==(best), labels)
    n = length(indices)
    n < 2 && return (; estimates, component_labels=labels, edge_count=0, iterations=0, residual_norm=0.0)
    local_ids = zeros(Int, nr)
    local_ids[indices] = 1:n
    local_pairs = Tuple{Int,Int}[]
    rhs, degree = zeros(n), zeros(n)
    ratios = Float64[]
    for (a,b) in pairs
        labels[a] == best || continue
        empty!(ratios)
        for p in run_precursors[a]
            ismissing(X[p,b]) || push!(ratios, Float64(X[p,b]) - Float64(X[p,a]))
        end
        ratio = median!(ratios)
        i, j = local_ids[a], local_ids[b]
        push!(local_pairs, (i,j))
        rhs[i] -= ratio
        rhs[j] += ratio
        degree[i] += 1
        degree[j] += 1
    end
    profile, iterations, residual_norm = _lfq_sparse_solve(
        local_pairs, rhs, degree; rtol, atol, maxiter)
    max_observed = maximum(Float64(X[p,r]) for r in indices for p in run_precursors[r])
    observed_sum = sum(exp2(Float64(X[p,r])-max_observed) for r in indices for p in run_precursors[r])
    max_profile = maximum(profile)
    scale = max_observed + log2(observed_sum) - max_profile -
        log2(sum(v -> exp2(v-max_profile), profile))
    for (i,r) in enumerate(indices)
        estimates[r] = Float32(profile[i] + scale)
    end
    return (; estimates, component_labels=labels, edge_count=length(local_pairs), iterations, residual_norm)
end
