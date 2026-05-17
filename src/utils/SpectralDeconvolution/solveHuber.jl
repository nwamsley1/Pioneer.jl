# Huber-loss coordinate descent for spectral deconvolution.
#
# Restores the robust solver used before the Poisson MM integration path,
# adapted to the current AbstractSparseDesignMatrix interface.

struct HuberSolver <: DeconvolutionSolver
    delta::Float32
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    reg_type::RegularizationType
end

function getRegL1(lambda::T, xk::T, ::NoNorm) where {T<:AbstractFloat}
    return zero(T)
end

function getRegL1(lambda::T, xk::T, ::L1Norm) where {T<:AbstractFloat}
    return lambda
end

function getRegL1(lambda::T, xk::T, ::L2Norm) where {T<:AbstractFloat}
    return T(2) * lambda * xk
end

function getRegL2(lambda::T, xk::T, ::NoNorm) where {T<:AbstractFloat}
    return zero(T)
end

function getRegL2(lambda::T, xk::T, ::L1Norm) where {T<:AbstractFloat}
    return zero(T)
end

function getRegL2(lambda::T, xk::T, ::L2Norm) where {T<:AbstractFloat}
    return T(2) * lambda
end

function updateHuberResiduals!(
    Hs::AbstractSparseDesignMatrix{Ti,T},
    r::Vector{T},
    col::Int64,
    X1::T,
    X0::T,
) where {Ti<:Integer,T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        r[row_val] += Hs.nzval[i] * (X1 - X0)
    end
    return nothing
end

@inline function huber_inv_sqrt_approx(x::Float32)
    bits = reinterpret(UInt32, x)
    y = reinterpret(Float32, UInt32(0x5f3759df) - (bits >> 1))
    return y * (1.5f0 - 0.5f0 * x * y^2)
end

function getHuberDerivatives!(
    Hs::AbstractSparseDesignMatrix{Ti,Float32},
    r::Vector{Float32},
    col::Int64,
    delta::Float32,
    lambda::Float32,
    xk::Float32,
    regularization_type::RegularizationType,
) where {Ti<:Integer}
    L1 = 0.0f0
    L2 = 0.0f0
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = 1.0f0 + (rval / delta)^2
        R = huber_inv_sqrt_approx(RS)
        hsval_r = hsval * R

        L1 += hsval_r * rval
        L2 += hsval * hsval_r * R^2
    end

    return L1 + getRegL1(lambda, xk, regularization_type),
           L2 + getRegL2(lambda, xk, regularization_type)
end

function getHuberL1(
    Hs::AbstractSparseDesignMatrix{Ti,Float32},
    r::Vector{Float32},
    col::Int64,
    delta::Float32,
    lambda::Float32,
    xk::Float32,
    regularization_type::RegularizationType,
) where {Ti<:Integer}
    L1 = 0.0f0
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = 1.0f0 + (rval / delta)^2
        R = huber_inv_sqrt_approx(RS)
        L1 += rval * hsval * R
    end
    return L1 + getRegL1(lambda, xk, regularization_type)
end

function huber_bisection!(
    Hs::AbstractSparseDesignMatrix{Ti,Float32},
    r::Vector{Float32},
    X1::Vector{Float32},
    col::Int64,
    delta::Float32,
    lambda::Float32,
    a::Float32,
    b::Float32,
    fa::Float32,
    max_iter::Int64,
    accuracy_bisection::Float32,
    regularization_type::RegularizationType,
) where {Ti<:Integer}
    n = 0
    c = (a + b) / 2.0f0
    updateHuberResiduals!(Hs, r, col, c, X1[col])
    X1[col] = c
    X_init = X1[col]
    X0 = X1[col]

    while n < max_iter
        fc = getHuberL1(Hs, r, col, delta, lambda, X1[col], regularization_type)
        if sign(fc) != sign(fa)
            b = c
        else
            a, fa = c, fc
        end

        c, X0 = (a + b) / 2.0f0, X1[col]
        X1[col] = c
        updateHuberResiduals!(Hs, r, col, X1[col], X0)

        abs(X1[col] - X0) < accuracy_bisection && break
        n += 1
    end

    return X1[col] - X_init
end

function huber_newton_bisection!(
    Hs::AbstractSparseDesignMatrix{Ti,Float32},
    r::Vector{Float32},
    X1::Vector{Float32},
    col::Int64,
    delta::Float32,
    lambda::Float32,
    max_iter_newton::Int64,
    max_iter_bisection::Int64,
    accuracy_newton::Float32,
    accuracy_bisection::Float32,
    regularization_type::RegularizationType,
    rel_tol::Float32 = 0.01f0,
) where {Ti<:Integer}
    n = 0
    X_init = X1[col]
    X0 = X1[col]
    max_l1, max_x1 = typemax(Float32), typemax(Float32)

    @inbounds begin
        while n < max_iter_newton
            L1, L2 = getHuberDerivatives!(Hs, r, col, delta, lambda, X1[col], regularization_type)
            update_rule = L1 / L2

            if isnan(update_rule)
                n = max_iter_newton
                break
            end

            if (sign(L1) == 1) & (L1 < max_l1)
                max_x1, max_l1 = X1[col], L1
            end

            X0 = X1[col]
            X1[col] = max(X1[col] - update_rule, 0.0f0)
            n += 1

            updateHuberResiduals!(Hs, r, col, X1[col], X0)

            abs_change = abs(X1[col] - X0)
            if !iszero(X0)
                abs_change / abs(X0) < rel_tol && break
            else
                abs_change < accuracy_newton && break
            end
        end

        if n == max_iter_newton
            X0 = X1[col]
            X1[col] = 0.0f0
            updateHuberResiduals!(Hs, r, col, X1[col], X0)
            L1 = getHuberL1(Hs, r, col, delta, lambda, X1[col], regularization_type)

            if sign(L1) != 1
                _ = huber_bisection!(
                    Hs, r, X1, col, delta, lambda, 0.0f0,
                    min(max(max_x1, 0.0f0), 1.0f11),
                    L1, max_iter_bisection, accuracy_bisection,
                    regularization_type,
                )
            end
        end
    end

    return X1[col] - X_init
end

function solveHuber!(
    Hs::AbstractSparseDesignMatrix{Ti,Float32},
    r::Vector{Float32},
    X1::Vector{Float32},
    delta::Float32,
    lambda::Float32,
    max_iter_newton::Int64,
    max_iter_bisection::Int64,
    max_iter_outer::Int64,
    accuracy_newton::Float32,
    accuracy_bisection::Float32,
    relative_convergence_threshold::Float32,
    regularization_type::RegularizationType,
) where {Ti<:Integer}
    newton_rel_tol = relative_convergence_threshold

    i = 0
    while i < max_iter_outer
        max_rel_change = 0.0f0
        for col in 1:Hs.n
            delta_x = abs(huber_newton_bisection!(
                Hs, r, X1, col, delta, lambda,
                max_iter_newton,
                max_iter_bisection,
                accuracy_newton,
                accuracy_bisection,
                regularization_type,
                newton_rel_tol,
            ))

            if !iszero(X1[col])
                rel_change = delta_x / abs(X1[col])
                max_rel_change = max(max_rel_change, rel_change)
            end
        end

        max_rel_change < relative_convergence_threshold && return (true, i)
        i += 1
    end

    return (false, max_iter_outer)
end

function with_huber_delta(solver::HuberSolver, delta::Float32)
    return HuberSolver(
        delta,
        solver.lambda,
        solver.max_iter_newton,
        solver.max_iter_bisection,
        solver.accuracy_newton,
        solver.accuracy_bisection,
        solver.reg_type,
    )
end

function solve_deconvolution!(solver::HuberSolver, Hs, r, w, colnorm2, mu, y, max_iter, conv)
    initResiduals!(r, Hs, w)
    return solveHuber!(
        Hs,
        r,
        w,
        solver.delta,
        solver.lambda,
        solver.max_iter_newton,
        solver.max_iter_bisection,
        max_iter,
        solver.accuracy_newton,
        solver.accuracy_bisection,
        conv,
        solver.reg_type,
    )
end
