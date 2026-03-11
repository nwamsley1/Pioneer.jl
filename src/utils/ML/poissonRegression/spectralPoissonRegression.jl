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
                          relative_convergence_threshold::T) where {Ti<:Integer, T<:AbstractFloat}

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
        initMu!(μ, Hs, X₁)
    end

    # ── Core MM iteration ──
    max_weight = T(0)
    ε = T(POISSON_MU_FLOOR)

    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            L1, L2 = getPoissonDerivativesObs!(Hs, μ, y, col)

            if L2 > ε && !isnan(L1)
                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                updateMu!(Hs, μ, col, X₁[col], X0)

                δx = abs(X₁[col] - X0)
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
        initMu!(μ, Hs, X₁)
    end

    return nothing
end

"""
    solvePoissonMM_v2!(Hs, μ, y, X₁, max_iter_outer,
                       relative_convergence_threshold;
                       max_inner_iter=5, reinit_period=10)

Like `solvePoissonMM!` but periodically recomputes μ = A*x from scratch
every `reinit_period` sweeps to prevent Float32 drift in `updateMu!`.
"""
function solvePoissonMM_v2!(Hs::SparseArray{Ti, T},
                             μ::Vector{T},
                             y::Vector{T},
                             X₁::Vector{T},
                             max_iter_outer::Int64,
                             relative_convergence_threshold::T;
                             reinit_period::Int64 = Int64(10)) where {Ti<:Integer, T<:AbstractFloat}

    # ── Y-scaling: divide y by max(y) to bring weights into tractable range ──
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
        initMu!(μ, Hs, X₁)
    end

    # ── Core MM iteration ──
    max_weight = T(0)
    ε = T(POISSON_MU_FLOOR)

    i = 0
    while i < max_iter_outer
        # Periodic μ recomputation to prevent Float32 drift
        if i > 0 && mod(i, reinit_period) == 0
            initMu!(μ, Hs, X₁)
        end

        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            L1, L2 = getPoissonDerivativesObs!(Hs, μ, y, col)

            if L2 > ε && !isnan(L1)
                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                updateMu!(Hs, μ, col, X₁[col], X0)

                δx = abs(X₁[col] - X0)
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
        initMu!(μ, Hs, X₁)
    end

    return nothing
end

"""
    getPoissonDerivativesObs_opt!(Hs, μ, y, col) → (L1, L2)

Optimized variant: replaces 2 divisions per nonzero with 1 division + multiplications.
`inv_μ = 1/μ_i` is computed once, then reused for both L1 and L2.
"""
function getPoissonDerivativesObs_opt!(Hs::SparseArray{Ti, T},
                                       μ::Vector{T},
                                       y::Vector{T},
                                       col::Int64) where {Ti<:Integer, T<:AbstractFloat}
    L1 = zero(Float32)
    L2 = zero(Float32)
    ε  = Float32(POISSON_MU_FLOOR)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row   = Hs.rowval[i]
        a_ij  = Hs.nzval[i]
        inv_μ = one(Float32) / max(μ[row], ε)
        y_i   = y[row]
        L1   += a_ij * (one(Float32) - y_i * inv_μ)
        L2   += a_ij * a_ij * y_i * inv_μ * inv_μ
    end
    # If L2 ≈ 0 (all y_i = 0 in this column's support), fall back to Fisher
    if L2 < ε
        @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            a_ij  = Hs.nzval[i]
            inv_μ = one(Float32) / max(μ[Hs.rowval[i]], ε)
            L2   += a_ij * a_ij * inv_μ
        end
    end
    return L1, L2
end

"""
    solvePoissonMM_opt!(Hs, μ, y, X₁, max_iter_outer,
                        relative_convergence_threshold)

Same as `solvePoissonMM!` but uses `getPoissonDerivativesObs_opt!` (1 div instead of 2).
"""
function solvePoissonMM_opt!(Hs::SparseArray{Ti, T},
                              μ::Vector{T},
                              y::Vector{T},
                              X₁::Vector{T},
                              max_iter_outer::Int64,
                              relative_convergence_threshold::T) where {Ti<:Integer, T<:AbstractFloat}

    # ── Y-scaling ──
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
        initMu!(μ, Hs, X₁)
    end

    # ── Core MM iteration ──
    max_weight = T(0)
    ε = T(POISSON_MU_FLOOR)

    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            L1, L2 = getPoissonDerivativesObs_opt!(Hs, μ, y, col)

            if L2 > ε && !isnan(L1)
                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                updateMu!(Hs, μ, col, X₁[col], X0)

                δx = abs(X₁[col] - X0)
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
        end

        if _diff < relative_convergence_threshold
            break
        end
        i += 1
    end

    # ── Unscale ──
    if y_scale > T(1)
        @inbounds for i in 1:Hs.m
            y[i] *= y_scale
        end
        @inbounds for j in 1:Hs.n
            X₁[j] *= y_scale
        end
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
