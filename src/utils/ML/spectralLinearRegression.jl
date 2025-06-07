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
                if X₁[col] > (accuracy_newton/max_diff)*2
                    if δx/X₁[col] > max_diff
                        _diff =  δx/X₁[col]
                    end
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
            abs_change = abs(X₁[col] - X0)
            
            # Check absolute tolerance
            #if abs_change < accuracy_newton
            #    break
            #end
            
            # Check convergence criteria
            if !iszero(X0)  # Can check relative change
                rel_change = abs_change / abs(X0)
                if rel_change < rel_tol  # Use the passed relative tolerance
                    break
                end
            else
                # For zero weights, use absolute tolerance only
                if abs_change < accuracy_newton
                    break
                end
            end
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
    
    # Debug: Print initial parameters (only for very small problems or rarely)
    #if Hs.n <= 5 || (Hs.n > 5 && rand() < 0.001)  # Sample 0.1% for larger problems
    #    println("\n[HUBER] Problem: $(Hs.n) vars, $(Hs.m) constraints, δ=$δ")
    #end
    
    ΔX = Inf
    #λ = Float64(0.1)
    max_iter_newton =100
    max_iter_bisection = 100 
    max_iter_outer = 1000#max(20, 2*Hs.n) #Outer loop iterations, at least 20 or 2*number of variables
    i = 0
    # Find initial max weight to establish dynamic range
    max_x, min_x = T(0), T(0)
    new_max_x = max_x
    # Base relative tolerance that tightens over iterations
    min_rel_tol, max_rel_tol = T(1e-7), T(0.01)#T(0.1) / T(2)^min(i, 15) # 10% → 5% → 2.5% → ... → 0.078%
    rel_tol = max_rel_tol 
    while (i < max_iter_outer) #& (ΔX > tol)
        ΔX = 0.0
        #Update each variable once 
        _diff = 0.0
        
        for col in range(1, Hs.n)
            # First iteration: use uniform tolerance for all coefficients
            if i == 0
                newton_rel_tol = T(0.1)  # 10% relative tolerance for first iteration
            else
                # Pre-compute linear scale constants for this iteration (only changes when max_x updates)
                rel_tol = T(10^(7 - log10(max_x) - log10(max(X₁[col], min_x))))
                #rel_tol = max_rel_tol + (x - min_weight_threshold)*((min_rel_tol - max_rel_tol)/(max_x - min_x))
                rel_tol = max(rel_tol, min_rel_tol)  # Ensure we don't go below minimum tolerance
                rel_tol = min(rel_tol, max_rel_tol)  # Ensure we don't exceed maximum tolerance
                newton_rel_tol = rel_tol 
            end
            
            #difference in X_1[col]
            δx = abs(newton_bisection!(Hs, r, X₁,col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection,
                                        regularization_type,
                                        newton_rel_tol))
            # Update max weight to keep dynamic range current
            # This ensures we have proper scaling after the first iteration
            if abs(X₁[col]) > new_max_x
                new_max_x = abs(X₁[col])
            end

            if !iszero(X₁[col]) 
                if abs(X₁[col]) > min_x
                    if δx/abs(X₁[col]) > max_diff
                        _diff =  δx/abs(X₁[col])
                    end
                end
            end
        end
        
        # Check convergence with one extra iteration for safety
        if _diff < max_diff
            break
        end
        max_x = max(new_max_x, T(1))
        min_x = max_x / T(1e4)  # Set minimum value to 1/10,000 of max_x
        new_max_x = T(0)  # Reset for next iteration
        # Ensure we have a reasonable minimum value for dynamic range calculation
        i += 1
    end
    
    # Debug: Final summary (only for problems that took many iterations)
    if i > 500 #&& rand() < 0.1  # 10% of problems that took > 100 iterations
        n_final_nonzero = sum(abs.(X₁) .> T(1e-10))
        println("[HUBER] Slow convergence: $i iters, $(Hs.n) vars, final ΔX=$ΔX")
    end
    
    return i
end


