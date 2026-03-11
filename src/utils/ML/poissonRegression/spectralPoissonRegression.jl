# Poisson MLE Coordinate Descent for Spectral Deconvolution
#
# Solves:  max_x  Σ_i [ y_i log(μ_i) - μ_i ]   subject to x_j ≥ 0
# where μ = Ax  (identity link Poisson GLM, no intercept)
#
# Uses Fisher-scoring coordinate descent:
#   gradient:     L1_j = Σ_i A_ij (1 - y_i / μ_i)
#   Fisher info:  L2_j = Σ_i A_ij² / μ_i
#   update:       x_j ← max(x_j - L1/L2, 0)

const POISSON_MU_FLOOR = 1f-6

"""
    updateMu!(Hs, μ, col, X1, X0)

Update predicted values μ after changing weight for column `col` from X0 to X1.
Reuses the same structure as `updateResiduals!` from the Huber solver.
"""
function updateMu!(Hs::SparseArray{Ti, T}, μ::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val  = Hs.nzval[i]
        μ[row_val] += nz_val * (X1 - X0)
    end
end

"""
    initObserved!(y, sa)

Extract observed values from `sa.x` into vector `y` (length ≥ sa.m).
Each row's observed value appears potentially multiple times in sa.x (once per
nonzero in that row); we just need one copy per row.
"""
function initObserved!(y::Vector{T}, sa::SparseArray{Ti,T}) where {Ti<:Integer, T<:AbstractFloat}
    if length(y) < sa.m
        append!(y, zeros(T, sa.m - length(y)))
    end
    @inbounds for i in 1:sa.m
        y[i] = zero(T)
    end
    @inbounds for n in 1:sa.n_vals
        row = sa.rowval[n]
        if iszero(y[row])
            y[row] = sa.x[n]
        end
    end
end

"""
    initMu!(μ, sa, w)

Compute μ = A * w from scratch, with floor clamping μ_i ≥ POISSON_MU_FLOOR.
"""
function initMu!(μ::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer, T<:AbstractFloat}
    if length(μ) < sa.m
        append!(μ, zeros(T, sa.m - length(μ)))
    end
    @inbounds for i in 1:sa.m
        μ[i] = zero(T)
    end
    @inbounds for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1] - 1)
            μ[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
    # Floor clamp
    ε = T(POISSON_MU_FLOOR)
    @inbounds for i in 1:sa.m
        if μ[i] < ε
            μ[i] = ε
        end
    end
end

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

Coordinate descent Poisson MLE solver.

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
        # After 5 iterations, ignore columns below max_weight * 1e-7 for convergence.
        # Below this floor, peaks are outside the instrument's dynamic range (~5 OOM)
        # and their oscillation is physically meaningless.
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            δx = abs(newton_bisection_poisson!(Hs, μ, y, X₁, col,
                                                max_iter_newton,
                                                max_iter_bisection,
                                                accuracy_newton,
                                                accuracy_bisection,
                                                newton_rel_tol))

            # Track max weight for significance floor
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            # Only count relative change for significant columns
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

# ══════════════════════════════════════════════════════════════════
# Observed-Hessian variant (true Newton-Raphson)
#
# The observed Hessian for the Poisson identity-link NLL is:
#   L2_obs = Σ_i A_ij² * y_i / μ_i²
#
# vs Fisher info:
#   L2_fisher = Σ_i A_ij² / μ_i
#
# The observed Hessian exactly solves each coordinate subproblem,
# guaranteeing monotonic decrease of -ℓ across outer iterations.
# ══════════════════════════════════════════════════════════════════

"""
    getPoissonDerivativesObs!(Hs, μ, y, col) → (L1, L2)

Gradient and observed Hessian for column `col`.

    L1 = Σ_i A_ij * (1 - y_i / μ_i)
    L2 = Σ_i A_ij² * y_i / μ_i²
"""
function getPoissonDerivativesObs!(Hs::SparseArray{Ti, T},
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
        L2   += a_ij * a_ij * y_i / (μ_i * μ_i)
    end
    # If L2 ≈ 0 (all y_i = 0 in this column's support), fall back to Fisher
    if L2 < ε
        @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            a_ij = Hs.nzval[i]
            μ_i  = max(μ[Hs.rowval[i]], ε)
            L2  += a_ij * a_ij / μ_i
        end
    end
    return L1, L2
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
        # After 5 iterations, ignore columns below max_weight * 1e-7 for convergence.
        # Below this floor, peaks are outside the instrument's dynamic range (~5 OOM)
        # and their oscillation is physically meaningless.
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            δx = abs(newton_bisection_poisson_obs!(Hs, μ, y, X₁, col,
                                                    max_iter_newton,
                                                    max_iter_bisection,
                                                    accuracy_newton,
                                                    accuracy_bisection,
                                                    newton_rel_tol))

            # Track max weight for significance floor
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            # Only count relative change for significant columns
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

# ══════════════════════════════════════════════════════════════════
# MM (Majorization-Minimization) variant — Cyclops-style
#
# Takes K observed-Hessian Newton steps per coordinate per sweep,
# recomputing derivatives after each step. This bridges large dynamic
# ranges that a single step cannot traverse in Float32 precision.
#
# Each step uses the observed Hessian: δ = L1/L2, w ← max(w - δ, 0).
# No bisection fallback needed — steps are bounded by non-negativity.
#
# Reference: Suchard et al. "Massive parallelization of serial
# inference algorithms for complex generalized linear models."
# ACM TOMACS 23(10), 2013.
# ══════════════════════════════════════════════════════════════════

"""
    solvePoissonMM!(Hs, μ, y, X₁, max_iter_outer,
                    relative_convergence_threshold;
                    max_inner_iter=5)

Coordinate descent Poisson MLE with MM (majorization-minimization) steps.
Takes `max_inner_iter` observed-Hessian Newton steps per coordinate per
sweep — much simpler than full Newton+bisection but handles large dynamic
ranges. Default inner iterations is 5.
"""
function solvePoissonMM!(Hs::SparseArray{Ti, T},
                          μ::Vector{T},
                          y::Vector{T},
                          X₁::Vector{T},
                          max_iter_outer::Int64,
                          relative_convergence_threshold::T;
                          max_inner_iter::Int64 = Int64(5)) where {Ti<:Integer, T<:AbstractFloat}

    # ── Y-scaling: divide y by max(y) to bring weights into tractable range ──
    # Without scaling, cold-start weights need to traverse 5-8 OOM from w=1
    # to w_true ~ O(1e4-1e8), but the observed Hessian makes steps microscopic.
    # Scaling y' = y/c means optimal w' = w/c, keeping steps well-sized.
    y_scale = T(0)
    @inbounds for i in 1:Hs.m
        if y[i] > y_scale
            y_scale = y[i]
        end
    end
    if y_scale > T(1)
        @inbounds for i in 1:Hs.m
            y[i] /= y_scale
        end
        @inbounds for j in 1:Hs.n
            X₁[j] /= y_scale
        end
        # Recompute μ = A * X₁ with scaled weights
        initMu!(μ, Hs, X₁)
    end

    # ── Core MM iteration ──
    max_weight = T(0)
    ε = T(POISSON_MU_FLOOR)
    inner_tol = relative_convergence_threshold

    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            X_before = X₁[col]

            # Multi-step MM: up to max_inner_iter observed-Hessian steps
            for _k in 1:max_inner_iter
                L1, L2 = getPoissonDerivativesObs!(Hs, μ, y, col)

                if L2 <= ε || isnan(L1)
                    break
                end

                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                updateMu!(Hs, μ, col, X₁[col], X0)

                # Early exit if this step was tiny
                abs_step = abs(X₁[col] - X0)
                if iszero(X₁[col]) || (!iszero(X0) && abs_step / abs(X0) < inner_tol)
                    break
                end
            end

            δx = abs(X₁[col] - X_before)

            # Convergence tracking with significance floor
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

    # ── Unscale: restore y and weights to original magnitude ──
    if y_scale > T(1)
        @inbounds for i in 1:Hs.m
            y[i] *= y_scale
        end
        @inbounds for j in 1:Hs.n
            X₁[j] *= y_scale
        end
        # Recompute μ = A * X₁ with unscaled weights so caller sees correct μ
        initMu!(μ, Hs, X₁)
    end

    return nothing
end

"""
    poissonLogLikelihood(μ, y, m)

Compute Poisson log-likelihood Σ_i [y_i log(μ_i) - μ_i] for diagnostic monitoring.
(Drops the -log(y_i!) constant.)
"""
function poissonLogLikelihood(μ::Vector{T}, y::Vector{T}, m::Int) where T<:AbstractFloat
    ll = zero(Float64)
    ε  = Float64(POISSON_MU_FLOOR)
    @inbounds for i in 1:m
        μ_i = max(Float64(μ[i]), ε)
        y_i = Float64(y[i])
        ll += y_i * log(μ_i) - μ_i
    end
    return ll
end
