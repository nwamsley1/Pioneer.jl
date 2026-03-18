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

# ── Solver trait types ──
abstract type DeconvolutionSolver end
struct OLSSolver <: DeconvolutionSolver end
struct PoissonMMSolver <: DeconvolutionSolver end

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
    solvePoissonMM_fast!(Hs, μ, y, X₁, max_iter_outer,
                          relative_convergence_threshold;
                          max_inner_iter=5)

Optimized Poisson MLE coordinate descent with fused μ-update + derivative computation.
Returns `(converged::Bool, iterations::Int)` tuple matching `solveOLS!` interface.
"""
function solvePoissonMM_fast!(Hs::SparseArray{Ti, T},
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

    # ── Extract fields for direct access ──
    colptr = Hs.colptr
    rowval = Hs.rowval
    nzval  = Hs.nzval
    ncols  = Hs.n

    max_weight = T(0)
    ε = T(POISSON_MU_FLOOR)

    iter = 0
    converged = false
    while iter < max_iter_outer
        _diff = T(0)
        weight_floor = iter >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:ncols
            col_start = colptr[col]
            col_end   = colptr[col + 1] - 1
            X_before  = X₁[col]

            # ── Initial derivative computation (separate pass) ──
            L1 = zero(Float32)
            L2 = zero(Float32)
            @inbounds @fastmath for i in col_start:col_end
                row   = rowval[i]
                a_ij  = nzval[i]
                inv_μ = one(Float32) / max(μ[row], ε)
                y_i   = y[row]
                L1   += a_ij * (one(Float32) - y_i * inv_μ)
                L2   += a_ij * a_ij * y_i * inv_μ * inv_μ
            end
            # Fisher fallback when observed Hessian ≈ 0 (all y_i = 0 in support)
            if L2 < ε
                @inbounds @fastmath for i in col_start:col_end
                    a_ij  = nzval[i]
                    inv_μ = one(Float32) / max(μ[rowval[i]], ε)
                    L2   += a_ij * a_ij * inv_μ
                end
            end

            # ── Inner Newton iterations ──
            for _k in 1:max_inner_iter
                (L2 <= ε || isnan(L1)) && break

                X0 = X₁[col]
                X₁[col] = max(X₁[col] - L1 / L2, zero(T))
                delta = X₁[col] - X0

                done = iszero(X₁[col]) || abs(delta) / max(abs(X₁[col]), T(1e-10)) < T(1e-3)

                if !done && _k < max_inner_iter
                    # ── Fused: update μ + compute next derivatives in one pass ──
                    L1 = zero(Float32)
                    L2 = zero(Float32)
                    @inbounds @fastmath for i in col_start:col_end
                        row       = rowval[i]
                        a_ij      = nzval[i]
                        μ[row]   += a_ij * delta
                        inv_μ     = one(Float32) / max(μ[row], ε)
                        y_i       = y[row]
                        L1       += a_ij * (one(Float32) - y_i * inv_μ)
                        L2       += a_ij * a_ij * y_i * inv_μ * inv_μ
                    end
                    if L2 < ε
                        @inbounds @fastmath for i in col_start:col_end
                            a_ij  = nzval[i]
                            inv_μ = one(Float32) / max(μ[rowval[i]], ε)
                            L2   += a_ij * a_ij * inv_μ
                        end
                    end
                else
                    # ── Just update μ (converged or last inner iteration) ──
                    @inbounds @fastmath for i in col_start:col_end
                        μ[rowval[i]] += nzval[i] * delta
                    end
                    break
                end
            end

            # ── Convergence tracking ──
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
            converged = true
            break
        end
        iter += 1
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

    return (converged, iter)
end

# ── Dispatch function ──
function solve_deconvolution!(::OLSSolver, Hs, r, w, colnorm2, μ, y, max_iter, conv)
    initResiduals!(r, Hs, w)
    return solveOLS!(Hs, r, w, colnorm2, max_iter, conv)
end

function solve_deconvolution!(::PoissonMMSolver, Hs, r, w, colnorm2, μ, y, max_iter, conv)
    # Resize residuals for downstream getDistanceMetrics (which recomputes r from scratch)
    if length(r) < Hs.m
        append!(r, zeros(eltype(r), Hs.m - length(r)))
    end
    initObserved!(y, Hs)
    initMu!(μ, Hs, w)
    return solvePoissonMM_fast!(Hs, μ, y, w, max_iter, conv)
end
