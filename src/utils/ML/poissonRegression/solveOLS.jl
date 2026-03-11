# Self-contained OLS v1 coordinate descent solver for spectral deconvolution.
#
# Solves:  min_x  ‖Ax - y‖²   subject to x_j ≥ 0
# via non-negative least squares coordinate descent.
#
# Each coordinate update is exact (one Newton step) because the OLS Hessian
# per coordinate is constant: L2_j = Σ_i A_ij².
#
# Usage:
#   include("SparseArray.jl")   # provides SparseArray{Ti,T}
#   include("solveOLS.jl")
#
#   # Initialize
#   r = Vector{Float32}(undef, sa.m)
#   w = ones(Float32, sa.n)            # or warm-start weights
#   colnorm2 = Vector{Float32}(undef, sa.n)
#   initResiduals!(r, sa, w)
#
#   # Solve
#   solveOLS!(sa, r, w, colnorm2, 1000, 1f-4)

"""
    updateResiduals!(Hs, r, col, X1, X0)

After changing weight `col` from X0 → X1, update residuals:  r[i] += A_ij * (X1 - X0).
"""
function updateResiduals!(Hs::SparseArray{Ti, T}, r::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val = Hs.nzval[i]
        r[row_val] += nz_val * (X1 - X0)
    end
end

"""
    initResiduals!(r, sa, w)

Compute residuals r = A*w - y from scratch.
"""
function initResiduals!(r::Vector{T}, sa::SparseArray{Ti,T}, w::Vector{T}) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds for i in 1:sa.m
        r[i] = zero(T)
    end
    # Set r[i] = -y[i] (observed values stored in sa.x, one per nonzero)
    @inbounds for n in 1:sa.n_vals
        if iszero(r[sa.rowval[n]])
            r[sa.rowval[n]] = -sa.x[n]
        end
    end
    # Add A*w
    @inbounds for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1] - 1)
            r[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
end

"""
    solveOLS!(Hs, r, X₁, colnorm2, max_iter_outer, relative_convergence_threshold)

Non-negative OLS coordinate descent.

Arguments:
- `Hs`:       SparseArray encoding the design matrix A (m rows × n columns)
- `r`:        residual vector (length ≥ m), must be initialized via `initResiduals!`
- `X₁`:       weight vector (length ≥ n), modified in-place
- `colnorm2`: scratch vector (length ≥ n), used internally for Σ A_ij²
- `max_iter_outer`: maximum number of coordinate sweeps
- `relative_convergence_threshold`: stop when max relative weight change < this
"""
function solveOLS!(
    Hs::SparseArray{Ti, T},
    r::Vector{T},
    X₁::Vector{T},
    colnorm2::Vector{T},
    max_iter_outer::Int64,
    relative_convergence_threshold::T
) where {Ti<:Integer, T<:AbstractFloat}

    # Precompute column norms: colnorm2[j] = Σ_i A_ij²
    # For OLS the Hessian per coordinate is constant, so one Newton step is exact.
    @inbounds for col in 1:Hs.n
        s = zero(T)
        for i in Hs.colptr[col]:(Hs.colptr[col+1]-1)
            s += Hs.nzval[i] * Hs.nzval[i]
        end
        colnorm2[col] = s
    end

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
            L2 = colnorm2[col]
            iszero(L2) && continue

            # Gradient: L1 = Σ_i A_ij * r_i
            L1 = zero(T)
            @inbounds @fastmath for k in Hs.colptr[col]:(Hs.colptr[col+1]-1)
                L1 += Hs.nzval[k] * r[Hs.rowval[k]]
            end

            # Exact coordinate minimizer (one step)
            X0 = X₁[col]
            X₁[col] = max(X₁[col] - L1 / L2, zero(T))

            # Update residuals
            updateResiduals!(Hs, r, col, X₁[col], X0)

            # Track convergence
            δx = abs(X₁[col] - X0)
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end
            if X₁[col] > weight_floor
                rc = δx / abs(X₁[col])
                rc > _diff && (_diff = rc)
            end
        end

        _diff < relative_convergence_threshold && break
        i += 1
    end
    return nothing
end

"""
    solveOLS_yscaled!(Hs, r, X₁, colnorm2, max_iter_outer, relative_convergence_threshold)

Y-scaled non-negative OLS coordinate descent.

Scales the observed values (sa.x) and weights by 1/max(y) before solving, then
unscales the weights afterwards.  This keeps all intermediate values in O(1) range,
preventing Float32 overflow on problems with large observed intensities.

Same interface as `solveOLS!` — call after `initResiduals!`.
Modifies `Hs.x` in place (scaled then unscaled), so Hs is restored on return.
"""
function solveOLS_yscaled!(
    Hs::SparseArray{Ti, T},
    r::Vector{T},
    X₁::Vector{T},
    colnorm2::Vector{T},
    max_iter_outer::Int64,
    relative_convergence_threshold::T
) where {Ti<:Integer, T<:AbstractFloat}

    # ── Compute y_scale = max observed value ──
    y_scale = zero(T)
    @inbounds for n in 1:Hs.n_vals
        Hs.x[n] > y_scale && (y_scale = Hs.x[n])
    end

    # ── Scale if needed ──
    if y_scale > one(T)
        inv_scale = one(T) / y_scale
        # Scale observed values
        @inbounds for n in 1:Hs.n_vals
            Hs.x[n] *= inv_scale
        end
        # Scale weights
        @inbounds for j in 1:Hs.n
            X₁[j] *= inv_scale
        end
        # Recompute residuals with scaled y and scaled w
        initResiduals!(r, Hs, X₁)
    end

    # ── Solve in scaled space ──
    solveOLS!(Hs, r, X₁, colnorm2, max_iter_outer, relative_convergence_threshold)

    # ── Unscale ──
    if y_scale > one(T)
        # Restore observed values
        @inbounds for n in 1:Hs.n_vals
            Hs.x[n] *= y_scale
        end
        # Unscale weights
        @inbounds for j in 1:Hs.n
            X₁[j] *= y_scale
        end
    end

    return nothing
end
