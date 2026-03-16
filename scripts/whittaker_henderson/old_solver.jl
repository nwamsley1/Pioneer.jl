#=
old_solver.jl — Standalone mock of the current whitsmddw implementation from Pioneer.jl

This extracts the exact sparse-matrix / CHOLMOD approach used in:
  src/utils/ML/wittakerHendersonSmoothing.jl

So we can benchmark it, test against alternatives, and verify correctness.
=#

using SparseArrays, LinearAlgebra, LinearSolve

const CHOLMOD_LOCK = ReentrantLock()

"""
    ddmat(x, d)

Compute the divided differencing matrix of order `d` for sample positions `x`.
Recursive construction: D_0 = I, D_d = V * diff(D_{d-1}).

Port of Paul Eilers' 2003 MATLAB function.
"""
function ddmat(x::AbstractVector{Float32}, d::Int64)
    m = length(x)
    if d == 0
        return spdiagm(0 => ones(Float64, m))
    else
        dx = x[d+1:end] .- x[1:end-d]
        V = spdiagm(0 => 1.0f0 ./ dx)
        D_lower = ddmat(x, d - 1)
        D_diff = D_lower[2:end, :] .- D_lower[1:end-1, :]
        return V * D_diff
    end
end

"""
    whitsmddw_old(x, y, w, n, λ; d=2)

Current (sparse CHOLMOD) Whittaker smoother with divided differences and weights.

Solves:  (W + λ D'D) z = W y
where W = diag(w), D = divided-difference matrix of order d.

This allocates ~5 sparse matrices + CHOLMOD factorization per call.
"""
function whitsmddw_old(x::AbstractVector{Float32},
                       y::AbstractVector{Float32},
                       w::AbstractVector{Float32},
                       n::Int,
                       λ::Float32 = 0.0002f0;
                       d::Int64 = 2)
    xs = @view x[1:n]
    ys = @view y[1:n]
    ws = @view w[1:n]

    D = ddmat(xs, d)
    W = spdiagm(0 => ws)
    M = W + λ * (D' * D)

    b = convert(Vector{Float64}, ws .* ys)
    prob = LinearProblem(Symmetric(M), b)

    sol = lock(CHOLMOD_LOCK) do
        solve(prob, CHOLMODFactorization())
    end

    z = Float32.(sol.u)
    @. z = max(z, 0f0)
    return z
end
