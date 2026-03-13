# Poisson MLE Coordinate Descent for Spectral Deconvolution
#
# Solves:  max_x  Σ_i [ y_i log(μ_i) - μ_i ]   subject to x_j ≥ 0
# where μ = Ax  (identity link Poisson GLM, no intercept)
#
# Uses observed-Hessian coordinate descent (MM / Cyclops-style):
#   gradient:     L1_j = Σ_i A_ij (1 - y_i / μ_i)
#   obs. Hessian: L2_j = Σ_i A_ij² y_i / μ_i²
#   update:       x_j ← max(x_j - L1/L2, 0)
#
# Reference: Suchard et al. "Massive parallelization of serial
# inference algorithms for complex generalized linear models."
# ACM TOMACS 23(10), 2013.

const POISSON_MU_FLOOR = 1f-6

"""
    updateMu!(Hs, μ, col, X1, X0)

Update predicted values μ after changing weight for column `col` from X0 to X1.
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

Uses `inv_μ = 1/μ_i` computed once per nonzero to avoid redundant divisions.
Falls back to Fisher Hessian (L2 = Σ A_ij²/μ_i) when observed L2 ≈ 0.
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
    solvePoissonMM!(Hs, μ, y, X₁, max_iter_outer,
                    relative_convergence_threshold;
                    max_inner_iter=5)

Coordinate descent Poisson MLE with MM (majorization-minimization) steps.
Takes `max_inner_iter` observed-Hessian Newton steps per coordinate per
sweep. Default inner iterations is 5.

Y-scaling is applied automatically: y and weights are divided by max(y) during
iteration to keep step sizes well-conditioned, then restored on exit.
"""
function solvePoissonMM!(Hs::SparseArray{Ti, T},
                          μ::Vector{T},
                          y::Vector{T},
                          X₁::Vector{T},
                          max_iter_outer::Int64,
                          relative_convergence_threshold::T;
                          max_inner_iter::Int64 = Int64(5)) where {Ti<:Integer, T<:AbstractFloat}

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
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            X_before = X₁[col]

            for _k in 1:max_inner_iter
                L1, L2 = getPoissonDerivativesObs!(Hs, μ, y, col)
                if L2 <= ε || isnan(L1)
                    break
                end
                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                updateMu!(Hs, μ, col, X₁[col], X0)
                if iszero(X₁[col]) || abs(X₁[col] - X0) / max(abs(X₁[col]), T(1e-10)) < T(1e-3)
                    break
                end
            end

            δx = abs(X₁[col] - X_before)
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end
            if X₁[col] > weight_floor
                rel_change = δx / max(abs(X₁[col]), T(1e-10))
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
