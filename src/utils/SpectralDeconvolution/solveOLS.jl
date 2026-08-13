# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Regularization type trait — stored on parameter structs but not currently
# dispatched on (retained for backward compatibility with config files).
abstract type RegularizationType end
struct L1Norm <: RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

function solveOLS!(Hs::AbstractSparseDesignMatrix{Ti, T},
                    r::Vector{T},
                    X₁::Vector{T},
                    colnorm2::Vector{T},
                    max_iter_outer::Int64,
                    relative_convergence_threshold::T) where {Ti<:Integer,T<:AbstractFloat}

    # Precompute column norms (sum of squared entries per column)
    @inbounds for col in 1:Hs.n
        s = zero(T)
        for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            s += Hs.nzval[i]^2
        end
        colnorm2[col] = s
    end

    max_weight = T(0)
    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            cn2 = colnorm2[col]
            iszero(cn2) && continue

            # Compute gradient (X^T r for this column)
            grad = zero(T)
            @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                grad += Hs.nzval[j] * r[Hs.rowval[j]]
            end

            # Exact coordinate update: new_x = x - grad / ||col||^2, clamped >= 0
            X0 = X₁[col]
            X₁[col] = max(X0 - grad / cn2, zero(T))

            # Update residuals in-place
            δx_raw = X₁[col] - X0
            if !iszero(δx_raw)
                @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                    r[Hs.rowval[j]] += Hs.nzval[j] * δx_raw
                end
            end

            # Track max weight
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            # Convergence tracking
            if X₁[col] > weight_floor
                rel_change = abs(δx_raw) / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end

        if _diff < relative_convergence_threshold
            return (true, i)
        end
        i += 1
    end

    return (false, max_iter_outer)
end

# ─── Diagnostic counters for adaptive-LASSO solver verification ─────────────
const DIAG_ALASSO_CALLS              = Threads.Atomic{Int64}(0)
const DIAG_ALASSO_NCOLS              = Threads.Atomic{Int64}(0)
const DIAG_ALASSO_ZEROS_AFTER_OLS    = Threads.Atomic{Int64}(0)
const DIAG_ALASSO_ZEROS_AFTER_FINAL  = Threads.Atomic{Int64}(0)
const DIAG_ALASSO_LAMBDA_MAX_SUM     = Threads.Atomic{Float64}(0.0)
const DIAG_ALASSO_LAMBDA_EFF_SUM     = Threads.Atomic{Float64}(0.0)

# ─────────────────────────────────────────────────────────────────────────────
# Weighted L1 + non-negativity inner solver (single weighted-LASSO subproblem).
# Identical to solveOLS! except `grad += λ * penalty_weights[col]` is added
# before the projected coordinate update. Since w ≥ 0, |w| = w and the
# subgradient of λ·p_j·|w_j| is exactly +λ·p_j (a constant), so the projected
# update is just OLS with a shifted gradient. This is the inner subproblem for
# the iterated adaptive LASSO.
# ─────────────────────────────────────────────────────────────────────────────
function solveWeightedL1NN!(Hs::AbstractSparseDesignMatrix{Ti, T},
                             r::Vector{T},
                             X₁::Vector{T},
                             colnorm2::Vector{T},
                             penalty_weights::Vector{T},
                             λ::T,
                             max_iter_outer::Int64,
                             relative_convergence_threshold::T) where {Ti<:Integer,T<:AbstractFloat}
    max_weight = T(0)
    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)
        for col in 1:Hs.n
            cn2 = colnorm2[col]
            iszero(cn2) && continue
            grad = zero(T)
            @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                grad += Hs.nzval[j] * r[Hs.rowval[j]]
            end
            grad += λ * penalty_weights[col]    # weighted-L1 derivative for w ≥ 0
            X0 = X₁[col]
            X₁[col] = max(X0 - grad / cn2, zero(T))
            δx_raw = X₁[col] - X0
            if !iszero(δx_raw)
                @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                    r[Hs.rowval[j]] += Hs.nzval[j] * δx_raw
                end
            end
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end
            if X₁[col] > weight_floor
                rel_change = abs(δx_raw) / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end
        if _diff < relative_convergence_threshold
            return (true, i)
        end
        i += 1
    end
    return (false, max_iter_outer)
end

# ─────────────────────────────────────────────────────────────────────────────
# Simple non-iterative L1-regularized NNLS.
# Algorithm:
#   1. Compute λ_max = max_j |c_j^T y| at w=0 (glmnet convention).
#   2. λ_eff = λ_rel · λ_max.
#   3. Solve weighted-L1 NNLS with uniform penalty_weights = 1 (standard LASSO).
# Returns (converged, iters_used).
# ─────────────────────────────────────────────────────────────────────────────
function solveLasso!(Hs::AbstractSparseDesignMatrix{Ti, T},
                      r::Vector{T},
                      X₁::Vector{T},
                      colnorm2::Vector{T},
                      λ_rel::T,
                      max_iter_outer::Int64,
                      relative_convergence_threshold::T) where {Ti<:Integer,T<:AbstractFloat}
    # Precompute column norms
    @inbounds for col in 1:Hs.n
        s = zero(T)
        for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            s += Hs.nzval[j]^2
        end
        colnorm2[col] = s
    end

    # Compute λ_max at w=0
    @inbounds for j in 1:Hs.n
        X₁[j] = zero(T)
    end
    initResiduals!(r, Hs, X₁)
    Threads.atomic_add!(DIAG_ALASSO_CALLS, 1)
    λ_max = zero(T)
    @inbounds for col in 1:Hs.n
        g = zero(T)
        for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            g += Hs.nzval[j] * r[Hs.rowval[j]]
        end
        ag = abs(g)
        if ag > λ_max
            λ_max = ag
        end
    end
    λ_eff = λ_rel * λ_max

    if λ_eff == zero(T)
        # No signal — fall back to OLS
        solveOLS!(Hs, r, X₁, colnorm2, max_iter_outer, relative_convergence_threshold)
        return (true, 0)
    end

    # Single weighted-L1 NNLS solve with uniform penalty weights
    penalty_weights = ones(T, Hs.n)
    initResiduals!(r, Hs, X₁)
    converged, iters = solveWeightedL1NN!(Hs, r, X₁, colnorm2,
                                            penalty_weights, λ_eff,
                                            max_iter_outer, relative_convergence_threshold)
    n_zero_now = 0
    @inbounds for j in 1:Hs.n
        n_zero_now += (X₁[j] == zero(T))
    end
    Threads.atomic_add!(DIAG_ALASSO_ZEROS_AFTER_FINAL, n_zero_now)
    Threads.atomic_add!(DIAG_ALASSO_NCOLS, Int(Hs.n))
    Threads.atomic_add!(DIAG_ALASSO_LAMBDA_MAX_SUM, Float64(λ_max))
    Threads.atomic_add!(DIAG_ALASSO_LAMBDA_EFF_SUM, Float64(λ_eff))
    return (converged, iters)
end

# ─────────────────────────────────────────────────────────────────────────────
# Iterated Adaptive LASSO (non-negative, identity-link, OLS loss).
# Algorithm:
#   1. Run unpenalized non-negative OLS → w⁽⁰⁾.
#   2. Compute λ_max = max_j |c_j^T r| at the OLS solution; λ_eff = λ_rel · λ_max.
#   3. For t = 1..n_iters:
#        a. ε = ε_rel · max_j w_j⁽ᵗ⁻¹⁾ (relative-magnitude perturbation).
#        b. p_j = 1 / (w_j⁽ᵗ⁻¹⁾ + ε)^γ.
#        c. Solve weighted-L1 NNLS subproblem (warm-started from w⁽ᵗ⁻¹⁾).
# Returns (converged_flag, total_outer_iters_used).
# ─────────────────────────────────────────────────────────────────────────────
function solveIteratedAdaptiveLasso!(Hs::AbstractSparseDesignMatrix{Ti, T},
                                      r::Vector{T},
                                      X₁::Vector{T},
                                      colnorm2::Vector{T},
                                      λ_rel::T,
                                      γ::T,
                                      ε_rel::T,
                                      n_iters::Int,
                                      max_iter_outer::Int64,
                                      relative_convergence_threshold::T) where {Ti<:Integer,T<:AbstractFloat}
    # Precompute column norms (solveOLS! does this; do it again here for the
    # fresh entry path).
    @inbounds for col in 1:Hs.n
        s = zero(T)
        for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            s += Hs.nzval[j]^2
        end
        colnorm2[col] = s
    end

    # Step 0: compute λ_max = max_j |c_j^T y| at w=0 (glmnet convention).
    # IMPORTANT: weights enter this routine warm-started from prior scans, so
    # we must explicitly zero X₁ before computing residuals; otherwise r is
    # already small and λ_max is meaningless.
    @inbounds for j in 1:Hs.n
        X₁[j] = zero(T)
    end
    initResiduals!(r, Hs, X₁)   # now r = -y (Hs convention: r = Hw - y)
    Threads.atomic_add!(DIAG_ALASSO_CALLS, 1)
    λ_max = zero(T)
    @inbounds for col in 1:Hs.n
        g = zero(T)
        for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            g += Hs.nzval[j] * r[Hs.rowval[j]]
        end
        ag = abs(g)
        if ag > λ_max
            λ_max = ag
        end
    end
    λ_eff = λ_rel * λ_max

    # Step 1: unpenalized OLS (initial estimator)
    solveOLS!(Hs, r, X₁, colnorm2, max_iter_outer, relative_convergence_threshold)
    n_zero_after_OLS = 0
    @inbounds for j in 1:Hs.n
        n_zero_after_OLS += (X₁[j] == zero(T))
    end
    Threads.atomic_add!(DIAG_ALASSO_ZEROS_AFTER_OLS, n_zero_after_OLS)
    Threads.atomic_add!(DIAG_ALASSO_NCOLS, Int(Hs.n))

    # Skip iteration if no signal at all
    if λ_eff == zero(T)
        return (true, 0)
    end

    # Step 3: iteration
    penalty_weights = Vector{T}(undef, Hs.n)
    total_iters = 0
    for t in 1:n_iters
        max_w = zero(T)
        @inbounds for j in 1:Hs.n
            if X₁[j] > max_w
                max_w = X₁[j]
            end
        end
        ε = max(ε_rel * max_w, T(1e-10))
        @inbounds for j in 1:Hs.n
            denom = (X₁[j] + ε)^γ
            penalty_weights[j] = one(T) / denom
        end

        initResiduals!(r, Hs, X₁)
        _, iters_used = solveWeightedL1NN!(Hs, r, X₁, colnorm2,
                                             penalty_weights, λ_eff,
                                             max_iter_outer, relative_convergence_threshold)
        total_iters += iters_used
        n_zero_now = 0
        @inbounds for j in 1:Hs.n
            n_zero_now += (X₁[j] == zero(T))
        end
        if t == n_iters
            Threads.atomic_add!(DIAG_ALASSO_ZEROS_AFTER_FINAL, n_zero_now)
        end
    end
    Threads.atomic_add!(DIAG_ALASSO_LAMBDA_MAX_SUM, Float64(λ_max))
    Threads.atomic_add!(DIAG_ALASSO_LAMBDA_EFF_SUM, Float64(λ_eff))
    return (true, total_iters)
end
