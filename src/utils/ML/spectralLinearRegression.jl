#=
abstract type RegularizationType end
struct L1Norm <: RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

function getRegL1(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL1(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return λ
end

function getRegL1(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return λ*Float32(2)*xk
end

function getRegL2(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return Float32(2)*λ
end

function updateResiduals!(Hs::SparseArray{Ti, T}, r::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val = Hs.nzval[i]
        r[row_val] += nz_val*(X1 - X0)
    end
end

function getDerivatives!(Hs::SparseArray{Ti, T}, 
                            r::Vector{Float32}, 
                            col::Int64, 
                            δ::Float32, 
                            λ::Float32,
                            xk::Float32,
                            regularization_type::RegularizationType
                            ) where {Ti<:Integer,T<:AbstractFloat}
    L1 = zero(Float32)
    L2 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 #Magic 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        HSVAL_R = hsval * R

        L1 += HSVAL_R * rval
        L2 += (hsval) * HSVAL_R * ((R)^(2))
    end

    return Float32(L1) + getRegL1(λ, xk, regularization_type), Float32(L2) + getRegL2(λ, xk, regularization_type)
end

function getL1(Hs::SparseArray{Ti, Float32}, 
                r::Vector{Float32}, 
                col::Int64, 
                δ::Float32, 
                λ::Float32,
                xk::Float32,
                regularization_type::RegularizationType) where {Ti<:Integer}
    L1 = zero(Float32)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 #Magic 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        L1 += rval * hsval * R
    end
    return Float32(L1) + getRegL1(λ, xk, regularization_type) #λ*Float32(2)*xk)
end

function newton_bisection!(Hs::SparseArray{Ti, T}, 
                            r::Vector{T}, 
                            X₁::Vector{T}, 
                            col::Int64, 
                            δ::T, 
                            λ::T,
                            max_iter_newton::Int64, 
                            max_iter_bisection::Int64,
                            accuracy_newton::T,
                            accuracy_bisection::T,
                            regulariation_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    X_init = X₁[col] #Estimate prior to optimiztion
    X0 = X₁[col] #Keeps track of previous etimate at each iteration
    #Maximum Estimates X₁[col] and L1. Used as uper bound if newton-raphson fails to converge
    #And bisection method is neede. 
    max_l1, max_x1 = typemax(T), typemax(T)

    #Newton-Raphson Method Iterations until convergence or maximum iterations. 
    #If convergence fails in maximum iterations, switch to bisection method for guaranteed convergence
    #@inbounds @fastmath begin
    @inbounds begin
        while (n < max_iter_newton)

            #First and second derivatives 
            L1, L2 = getDerivatives!(Hs, r, col, δ, λ, X₁[col], regulariation_type)
            update_rule = (L1)/L2
            #Switch to bisection method
            if isnan(update_rule)
                n = max_iter_newton
                break
            end

            #Useful as boundaries for bisection method if Newton's method fails to converge. 
            #Want positive and negative L1 with smallest absolute values 
            if (sign(L1) == 1) & (L1 < max_l1) #L1 > max_l1
                max_x1, max_l1 = X₁[col], L1
            end

            X0 = X₁[col] 
            #Newton-Raphson update. Contrain X₁[col] to be non-negative
            X₁[col] = max(X₁[col] - update_rule, zero(T))
            n += 1

            #Update residuals given new estimate, X₁[col], and prior estimate, X0
            updateResiduals!(Hs, r, col, X₁[col], X0)

            ########
            #Stopping Criterion for single variable Newton-Raphson
            ########
            abs(X₁[col] - X0) < accuracy_newton ? break : nothing
        end
        #If newtons method fails to converge, switch to bisection method
        if n == max_iter_newton
            #println("bisection")
            #######
            #Lower bound is always X₁[col] == 0
            X0 = X₁[col]
            X₁[col] = zero(T)
            updateResiduals!(Hs, r, col, X₁[col], X0)
            L1 = getL1(Hs, r, col, δ, λ, X₁[col], regulariation_type)
            if sign(L1) != 1 #Otherwise the minimum is at X₁[col] < 0, so set X₁[col] == 0
                _ = bisection!(Hs, r, X₁, col, δ, λ, zero(T), 
                            min(max(max_x1, zero(Float32)), Float32(1e11)), #Maximum Plausible Value for X1
                            L1,  
                            max_iter_bisection, #Should never reach this. Convergence in (max_x1)/2^n
                            accuracy_bisection,
                            regulariation_type)#accuracy)
            end
            return X₁[col] - X_init
        else #Convergence reached. Return difference between current estimate and initial guess. 
            return X₁[col] - X_init
        end
    end
end

function bisection!(Hs::SparseArray{Ti, T}, 
                    r::Vector{T}, 
                    X₁::Vector{T},
                    col::Int64, 
                    δ::T, 
                    λ::T, 
                    a::T, 
                    b::T, 
                    fa::Float32,
                    max_iter::Int64, 
                    accuracy_bisection::T,
                    regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    c = (a + b)/2
    #Since first guess for X₁[col] will change to "c"
    #Need to update the residuals accordingly. 
    updateResiduals!(Hs, r, col, c, X₁[col])
    X0 = X₁[col]
    X₁[col] = c
    X_init, X0 = X₁[col],  X₁[col]
    while (n < max_iter)

        #Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ, λ, X₁[col], regularization_type)
        #Bisection Rule
        if (sign(fc) != sign(fa))
            b, fb = c, fc
        else
            a, fa = c, fc
        end

        c, X0 = (a + b)/2, X₁[col]
        X₁[col] = c
        #Update residuals given new estimate, X₁[col], and prior estimate, X0
        updateResiduals!(Hs, r, col, X₁[col], X0)

        ########
        #Stopping Criterion
        ########
        abs(X₁[col] - X0) < accuracy_bisection ? break : nothing
        n += 1
    end
    return X₁[col] - X_init
end

function solveHuber!(Hs::SparseArray{Ti, T}, 
                        r::Vector{T}, 
                        X₁::Vector{T}, 
                        δ::T,
                        λ::T,
                        max_iter_newton::Int64, 
                        max_iter_bisection::Int64, 
                        max_iter_outer::Int64,
                        accuracy_newton::T,
                        accuracy_bisection::T,
                        tol::U,
                        max_diff::T,
                        regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat,U<:Real}
    ΔX = Inf
    #λ = Float64(0.1)
    max_iter_newton = 50
    max_iter_bisection = 100 
    max_iter_outer = 1000
    i = 0
    while (i < max_iter_outer) & (ΔX > tol)
        ΔX = 0.0
        #Update each variable once 
        _diff = 0.0
        for col in range(1, Hs.n)
            #difference in X_1[col]
            δx = newton_bisection!(Hs, r, X₁,col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection,
                                        regularization_type)
            δx = abs(δx)
            if !iszero(X₁[col]) 
                if δx/X₁[col] > max_diff
                    _diff =  δx/X₁[col]
                end
            end
            
            ΔX += δx
        end
        if _diff < max_diff
            break
        end
       i += 1
    end
    return i
end
=#
abstract type RegularizationType end
struct L1Norm <: RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

function getRegL1(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL1(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return λ
end

function getRegL1(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return λ*Float32(2)*xk
end

function getRegL2(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return Float32(2)*λ
end

function updateResiduals!(Hs::SparseArray{Ti, T}, r::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val = Hs.nzval[i]
        r[row_val] += nz_val*(X1 - X0)
    end
end

function getDerivatives!(Hs::SparseArray{Ti, T}, 
                            r::Vector{Float32}, 
                            col::Int64, 
                            δ::Float32, 
                            λ::Float32,
                            xk::Float32,
                            regularization_type::RegularizationType
                            ) where {Ti<:Integer,T<:AbstractFloat}
    L1 = zero(Float32)
    L2 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 #Magic 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        HSVAL_R = hsval * R

        L1 += HSVAL_R * rval
        L2 += (hsval) * HSVAL_R * ((R)^(2))
    end

    return Float32(L1) + getRegL1(λ, xk, regularization_type), Float32(L2) + getRegL2(λ, xk, regularization_type)
end

function getL1(Hs::SparseArray{Ti, Float32}, 
                r::Vector{Float32}, 
                col::Int64, 
                δ::Float32, 
                λ::Float32,
                xk::Float32,
                regularization_type::RegularizationType) where {Ti<:Integer}
    L1 = zero(Float32)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 #Magic 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        L1 += rval * hsval * R
    end
    return Float32(L1) + getRegL1(λ, xk, regularization_type) #λ*Float32(2)*xk)
end

function newton_bisection!(Hs::SparseArray{Ti, T}, 
                            r::Vector{T}, 
                            X₁::Vector{T}, 
                            col::Int64, 
                            δ::T, 
                            λ::T,
                            max_iter_newton::Int64, 
                            max_iter_bisection::Int64,
                            accuracy_newton::T,
                            accuracy_bisection::T,
                            regulariation_type::RegularizationType,
                            rel_tol::T = T(0.01)) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    X_init = X₁[col]  # Initial estimate before optimization
    X0 = X₁[col]      # Previous estimate at each iteration
    
    # Track maximum estimates for bisection bounds if Newton fails
    max_l1, max_x1 = typemax(T), typemax(T)

    @inbounds begin
        # Newton-Raphson iterations
        while (n < max_iter_newton)
            # Compute first and second derivatives
            L1, L2 = getDerivatives!(Hs, r, col, δ, λ, X₁[col], regulariation_type)
            update_rule = (L1)/L2
            
            # Check for numerical issues
            if isnan(update_rule)
                n = max_iter_newton
                break
            end

            # Track positive L1 with smallest absolute value for bisection bounds
            if (sign(L1) == 1) & (L1 < max_l1)
                max_x1, max_l1 = X₁[col], L1
            end

            # Newton-Raphson update (constrained to be non-negative)
            X0 = X₁[col] 
            X₁[col] = max(X₁[col] - update_rule, zero(T))
            n += 1

            # Update residuals
            updateResiduals!(Hs, r, col, X₁[col], X0)

            # Check convergence
            abs_change = abs(X₁[col] - X0)
            
            if !iszero(X0)  # Can check relative change
                rel_change = abs_change / abs(X0)
                if rel_change < rel_tol
                    break
                end
            else
                # For zero weights, use absolute tolerance only
                if abs_change < accuracy_newton
                    break
                end
            end
        end
        
        # If Newton's method failed to converge, use bisection
        if n == max_iter_newton
            # Reset to zero (lower bound)
            X0 = X₁[col]
            X₁[col] = zero(T)
            updateResiduals!(Hs, r, col, X₁[col], X0)
            L1 = getL1(Hs, r, col, δ, λ, X₁[col], regulariation_type)
            
            # Only use bisection if minimum is at X₁[col] > 0
            if sign(L1) != 1
                _ = bisection!(Hs, r, X₁, col, δ, λ, zero(T), 
                            min(max(max_x1, zero(Float32)), Float32(1e11)),
                            L1,  
                            max_iter_bisection,
                            accuracy_bisection,
                            regulariation_type)
            end
            return X₁[col] - X_init
        else
            # Convergence reached
            return X₁[col] - X_init
        end
    end
end

function bisection!(Hs::SparseArray{Ti, T}, 
                    r::Vector{T}, 
                    X₁::Vector{T},
                    col::Int64, 
                    δ::T, 
                    λ::T, 
                    a::T, 
                    b::T, 
                    fa::Float32,
                    max_iter::Int64, 
                    accuracy_bisection::T,
                    regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    c = (a + b)/2
    #Since first guess for X₁[col] will change to "c"
    #Need to update the residuals accordingly. 
    updateResiduals!(Hs, r, col, c, X₁[col])
    X0 = X₁[col]
    X₁[col] = c
    X_init, X0 = X₁[col],  X₁[col]
    while (n < max_iter)

        #Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ, λ, X₁[col], regularization_type)
        #Bisection Rule
        if (sign(fc) != sign(fa))
            b, fb = c, fc
        else
            a, fa = c, fc
        end

        c, X0 = (a + b)/2, X₁[col]
        X₁[col] = c
        #Update residuals given new estimate, X₁[col], and prior estimate, X0
        updateResiduals!(Hs, r, col, X₁[col], X0)

        ########
        #Stopping Criterion
        ########
        abs(X₁[col] - X0) < accuracy_bisection ? break : nothing
        n += 1
    end
    return X₁[col] - X_init
end

function solveHuber!(Hs::SparseArray{Ti, T}, 
                        r::Vector{T}, 
                        X₁::Vector{T}, 
                        δ::T,
                        λ::T,
                        max_iter_newton::Int64, 
                        max_iter_bisection::Int64, 
                        max_iter_outer::Int64,
                        accuracy_newton::T,
                        accuracy_bisection::T,
                        tol::U,
                        max_diff::T,
                        regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat,U<:Real}
    
    # Override solver parameters for performance
    max_iter_newton = 100
    max_iter_bisection = 100 
    max_iter_outer = 1000
    
    # Initialize iteration counter and dynamic range tracking
    i = 0
    max_x = T(0)
    min_x = T(0)
    new_max_x = T(0)
    
    # Tolerance bounds for adaptive convergence
    min_rel_tol = T(1e-7)
    max_rel_tol = T(0.01)
    
    while i < max_iter_outer
        _diff = T(0)
        
        for col in range(1, Hs.n)
            # Determine relative tolerance for newton_bisection
            if i == 0
                newton_rel_tol = T(0.1)  # 10% relative tolerance for first iteration
            else
                # Adaptive tolerance based on coefficient magnitude
                rel_tol = T(10^(-2 - (log10(min(X₁[col], max_x)) - log10(min_x))))
                rel_tol = clamp(rel_tol, min_rel_tol, max_rel_tol)
                newton_rel_tol = rel_tol 
            end
            
            # Update coefficient
            δx = abs(newton_bisection!(Hs, r, X₁, col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection,
                                        regularization_type,
                                        newton_rel_tol))
            
            # Track maximum coefficient magnitude
            if abs(X₁[col]) > new_max_x
                new_max_x = abs(X₁[col])
            end

            # Track maximum relative change
            if !iszero(X₁[col]) && abs(X₁[col]) > min_x
                rel_change = δx / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end
        
        # Check convergence
        if _diff < max_diff
            break
        end
        
        # Update dynamic range for next iteration
        max_x = max(new_max_x, T(1))
        min_x = max_x / T(1e4)
        new_max_x = T(0)
        
        i += 1
    end
    
    # Debug output for slow convergence
    if i > 500
        println("[HUBER] Slow convergence: $i iters, $(Hs.n) vars")
    end
    
    return i
end