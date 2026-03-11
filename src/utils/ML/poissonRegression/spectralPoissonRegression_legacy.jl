# Legacy Poisson solvers — Fisher-scoring and observed-Hessian Newton+bisection variants.
#
# These functions are NOT used by the production MM solvers (solvePoissonMM! / solvePoissonMM_v2!).
# They are kept for test files (test_poisson.jl, test_glm_comparison.jl).
#
# Depends on: POISSON_MU_FLOOR, updateMu!, getPoissonDerivativesObs!
#   from spectralPoissonRegression.jl — must be included AFTER that file.

"""
    getPoissonDerivatives!(Hs, μ, y, col) → (L1, L2)

Compute gradient (L1) and Fisher information (L2) for column `col`.

    L1 = Σ_i A_ij * (1 - y_i / μ_i)
    L2 = Σ_i A_ij² / μ_i
"""
function getPoissonDerivatives!(Hs::SparseArray{Ti, T},
                                 μ::Vector{T},
                                 y::Vector{T},
                                 col::Int64) where {Ti<:Integer, T<:AbstractFloat}
    L1 = zero(Float32)
    L2 = zero(Float32)
    ε  = Float32(POISSON_MU_FLOOR)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row   = Hs.rowval[i]
        a_ij  = Hs.nzval[i]
        μ_i   = max(μ[row], ε)
        y_i   = y[row]
        L1   += a_ij * (one(Float32) - y_i / μ_i)
        L2   += a_ij * a_ij / μ_i
    end
    return L1, L2
end

"""
    getPoissonL1(Hs, μ, y, col) → L1

Gradient only (for bisection).
"""
function getPoissonL1(Hs::SparseArray{Ti, T},
                       μ::Vector{T},
                       y::Vector{T},
                       col::Int64) where {Ti<:Integer, T<:AbstractFloat}
    L1 = zero(Float32)
    ε  = Float32(POISSON_MU_FLOOR)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row   = Hs.rowval[i]
        a_ij  = Hs.nzval[i]
        μ_i   = max(μ[row], ε)
        y_i   = y[row]
        L1   += a_ij * (one(Float32) - y_i / μ_i)
    end
    return L1
end

"""
    bisection_poisson!(Hs, μ, y, X₁, col, a, b, fa, max_iter, accuracy_bisection)

Bisection root-finder on the Poisson gradient for column `col`.
"""
function bisection_poisson!(Hs::SparseArray{Ti, T},
                             μ::Vector{T},
                             y::Vector{T},
                             X₁::Vector{T},
                             col::Int64,
                             a::T, b::T, fa::Float32,
                             max_iter::Int64,
                             accuracy_bisection::T) where {Ti<:Integer, T<:AbstractFloat}
    n = 0
    c = (a + b) / 2
    # Move μ to reflect X₁[col] = c
    updateMu!(Hs, μ, col, c, X₁[col])
    X0 = X₁[col]
    X₁[col] = c
    X_init = X₁[col]
    X0 = X₁[col]

    while n < max_iter
        fc = getPoissonL1(Hs, μ, y, col)
        if sign(fc) != sign(fa)
            b = c
        else
            a = c
            fa = fc
        end
        c  = (a + b) / 2
        X0 = X₁[col]
        X₁[col] = c
        updateMu!(Hs, μ, col, X₁[col], X0)
        abs(X₁[col] - X0) < accuracy_bisection && break
        n += 1
    end
    return X₁[col] - X_init
end

"""
    newton_bisection_poisson!(Hs, μ, y, X₁, col, ...)

Newton-Raphson (Fisher scoring) with bisection fallback for column `col`.
"""
function newton_bisection_poisson!(Hs::SparseArray{Ti, T},
                                    μ::Vector{T},
                                    y::Vector{T},
                                    X₁::Vector{T},
                                    col::Int64,
                                    max_iter_newton::Int64,
                                    max_iter_bisection::Int64,
                                    accuracy_newton::T,
                                    accuracy_bisection::T,
                                    rel_tol::T = T(0.01)) where {Ti<:Integer, T<:AbstractFloat}
    n = 0
    X_init = X₁[col]
    X0     = X₁[col]
    max_l1 = typemax(T)
    max_x1 = typemax(T)

    @inbounds begin
        while n < max_iter_newton
            L1, L2 = getPoissonDerivatives!(Hs, μ, y, col)
            update_rule = L1 / L2

            if isnan(update_rule) || iszero(L2)
                n = max_iter_newton
                break
            end

            # Track positive-gradient bound for bisection
            if (sign(L1) == 1) && (L1 < max_l1)
                max_x1, max_l1 = X₁[col], L1
            end

            X0 = X₁[col]
            X₁[col] = max(X₁[col] - update_rule, zero(T))
            n += 1

            updateMu!(Hs, μ, col, X₁[col], X0)

            # Convergence check
            abs_change = abs(X₁[col] - X0)
            if !iszero(X0)
                if abs_change / abs(X0) < rel_tol
                    break
                end
            else
                if abs_change < accuracy_newton
                    break
                end
            end
        end

        # Bisection fallback
        if n == max_iter_newton
            X0 = X₁[col]
            X₁[col] = zero(T)
            updateMu!(Hs, μ, col, X₁[col], X0)
            L1 = getPoissonL1(Hs, μ, y, col)

            if sign(L1) != 1
                bisection_poisson!(Hs, μ, y, X₁, col,
                                   zero(T),
                                   min(max(max_x1, zero(Float32)), Float32(1e11)),
                                   L1,
                                   max_iter_bisection,
                                   accuracy_bisection)
            end
            return X₁[col] - X_init
        else
            return X₁[col] - X_init
        end
    end
end

"""
    solvePoisson!(Hs, μ, y, X₁, max_iter_newton, max_iter_bisection, max_iter_outer,
                  accuracy_newton, accuracy_bisection, relative_convergence_threshold)

Coordinate descent Poisson MLE solver (Fisher scoring).

- `Hs`: SparseArray (design matrix A in CSC-like format, with observed values in Hs.x)
- `μ`:  predicted intensities vector (length ≥ Hs.m), will be mutated
- `y`:  observed intensities vector (length ≥ Hs.m)
- `X₁`: weight vector (length ≥ Hs.n), will be mutated
"""
function solvePoisson!(Hs::SparseArray{Ti, T},
                        μ::Vector{T},
                        y::Vector{T},
                        X₁::Vector{T},
                        max_iter_newton::Int64,
                        max_iter_bisection::Int64,
                        max_iter_outer::Int64,
                        accuracy_newton::T,
                        accuracy_bisection::T,
                        relative_convergence_threshold::T) where {Ti<:Integer, T<:AbstractFloat}

    newton_rel_tol = relative_convergence_threshold
    max_weight = T(0)

    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            δx = abs(newton_bisection_poisson!(Hs, μ, y, X₁, col,
                                                max_iter_newton,
                                                max_iter_bisection,
                                                accuracy_newton,
                                                accuracy_bisection,
                                                newton_rel_tol))

            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            if X₁[col] > weight_floor
                rel_change = δx / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end
        if _diff < relative_convergence_threshold
            break
        end
        i += 1
    end
    return nothing
end

"""
    newton_bisection_poisson_obs!(Hs, μ, y, X₁, col, ...)

Newton-Raphson with observed Hessian and bisection fallback.
"""
function newton_bisection_poisson_obs!(Hs::SparseArray{Ti, T},
                                        μ::Vector{T},
                                        y::Vector{T},
                                        X₁::Vector{T},
                                        col::Int64,
                                        max_iter_newton::Int64,
                                        max_iter_bisection::Int64,
                                        accuracy_newton::T,
                                        accuracy_bisection::T,
                                        rel_tol::T = T(0.01)) where {Ti<:Integer, T<:AbstractFloat}
    n = 0
    X_init = X₁[col]
    X0     = X₁[col]
    max_l1 = typemax(T)
    max_x1 = typemax(T)

    @inbounds begin
        while n < max_iter_newton
            L1, L2 = getPoissonDerivativesObs!(Hs, μ, y, col)
            update_rule = L1 / L2

            if isnan(update_rule) || iszero(L2)
                n = max_iter_newton
                break
            end

            if (sign(L1) == 1) && (L1 < max_l1)
                max_x1, max_l1 = X₁[col], L1
            end

            X0 = X₁[col]
            X₁[col] = max(X₁[col] - update_rule, zero(T))
            n += 1

            updateMu!(Hs, μ, col, X₁[col], X0)

            abs_change = abs(X₁[col] - X0)
            if !iszero(X0)
                if abs_change / abs(X0) < rel_tol
                    break
                end
            else
                if abs_change < accuracy_newton
                    break
                end
            end
        end

        if n == max_iter_newton
            X0 = X₁[col]
            X₁[col] = zero(T)
            updateMu!(Hs, μ, col, X₁[col], X0)
            L1 = getPoissonL1(Hs, μ, y, col)

            if sign(L1) != 1
                bisection_poisson!(Hs, μ, y, X₁, col,
                                   zero(T),
                                   min(max(max_x1, zero(Float32)), Float32(1e11)),
                                   L1,
                                   max_iter_bisection,
                                   accuracy_bisection)
            end
            return X₁[col] - X_init
        else
            return X₁[col] - X_init
        end
    end
end

"""
    solvePoissonObs!(Hs, μ, y, X₁, ...)

Coordinate descent Poisson MLE with observed Hessian (true Newton).
"""
function solvePoissonObs!(Hs::SparseArray{Ti, T},
                           μ::Vector{T},
                           y::Vector{T},
                           X₁::Vector{T},
                           max_iter_newton::Int64,
                           max_iter_bisection::Int64,
                           max_iter_outer::Int64,
                           accuracy_newton::T,
                           accuracy_bisection::T,
                           relative_convergence_threshold::T) where {Ti<:Integer, T<:AbstractFloat}

    newton_rel_tol = relative_convergence_threshold
    max_weight = T(0)

    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            δx = abs(newton_bisection_poisson_obs!(Hs, μ, y, X₁, col,
                                                    max_iter_newton,
                                                    max_iter_bisection,
                                                    accuracy_newton,
                                                    accuracy_bisection,
                                                    newton_rel_tol))

            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            if X₁[col] > weight_floor
                rel_change = δx / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end
        if _diff < relative_convergence_threshold
            break
        end
        i += 1
    end
    return nothing
end
