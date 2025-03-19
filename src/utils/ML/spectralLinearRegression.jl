abstract type RegularizationType end
struct L1Norm <: RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

#=
function getDerivatives!(Hs::SparseArray{Ti, T}, 
                            r::Vector{Float32}, 
                            col::Int64, 
                            δ::Float32, 
                            λ::Float32) where {Ti<:Integer,T<:AbstractFloat}
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
    return L1 + λ, L2
end

function getL1(Hs::SparseArray{Ti, Float32}, r::Vector{Float32}, col::Int64, δ::Float32, λ::Float32) where {Ti<:Integer}
    L1 = zero(Float32)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #Huber
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
    return L1 + λ
end

function getL1(Hs::SparseArray{Ti, Float64}, r::Vector{Float64}, col::Int64, δ::Float64, λ::Float64) where {Ti<:Integer}
    L1 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #Huber
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int64 = reinterpret(UInt64, xₛ)
        int64 = 0x5fe6eb50c7b537a9 - int64 >> 1
        R = reinterpret(Float64, int64)
        R *= 1.5 - RS * 0.5 * R^2
        L1 += rval * hsval * R
    end
    return L1 + λ
end

function updateResiduals!(Hs::SparseArray{Ti, T}, r::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val = Hs.nzval[i]
        r[row_val] += nz_val*(X1 - X0)
    end
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
                            accuracy_bisection::T) where {Ti<:Integer,T<:AbstractFloat}
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
            L1, L2 = getDerivatives!(Hs, r, col, δ, λ)
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
            L1 = getL1(Hs, r, col, δ, λ)
            if sign(L1) != 1 #Otherwise the minimum is at X₁[col] < 0, so set X₁[col] == 0
                bisection!(Hs, r, X₁, col, δ, λ, zero(T), 
                            min(max_x1, Float32(1e11)), #Maximum Plausible Value for X1
                            L1,  
                            max_iter_bisection, #Should never reach this. Convergence in (max_x1)/2^n
                            accuracy_bisection)#accuracy)
            end
            X₁[col] - X_init
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
                    fa::T,
                    max_iter::Int64, 
                    accuracy_bisection::T) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    c = (a + b)/2
    #Since first guess for X₁[col] will change to "c"
    #Need to update the residuals accordingly. 
    updateResiduals!(Hs, r, col, c, X₁[col])

    X₁[col] = c
    X_init, X0 = X₁[col],  X₁[col]
    while (n < max_iter)

        #Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ, λ)
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
                        ::RegularizationType) where {Ti<:Integer,T<:AbstractFloat,U<:Real}
    ΔX = Inf
    i = 0
    while (i < max_iter_outer) & (ΔX > tol)
        ΔX = 0.0
        #Update each variable once 
        _diff = 0.0
        for col in range(1, Hs.n)
            #difference in X_1[col]
            δx = abs(newton_bisection!(Hs, r, X₁, col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection))
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
########
#For L2 Norm 

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
                            x2_sum::Float32,
                            xk::Float32
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

    return L1 + (λ)*xk/sqrt(x2_sum), L2 + ((λ)/sqrt(x2_sum))*(1 - (xk^2)/x2_sum)
end

function getL1(Hs::SparseArray{Ti, Float32}, r::Vector{Float32}, col::Int64, δ::Float32, λ::Float32,
    x2_sum::Float32,
    xk::Float32) where {Ti<:Integer}
    L1 = zero(Float32)
    #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #Huber
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
    return L1 + (λ)*xk/sqrt(x2_sum)
end

#=
function getL1(Hs::SparseArray{Ti, Float64}, r::Vector{Float64}, col::Int64, δ::Float64, λ::Float64) where {Ti<:Integer}
    L1 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #Huber
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int64 = reinterpret(UInt64, xₛ)
        int64 = 0x5fe6eb50c7b537a9 - int64 >> 1
        R = reinterpret(Float64, int64)
        R *= 1.5 - RS * 0.5 * R^2
        L1 += rval * hsval * R
    end
    return L1 + λ
end
=#
#=
function updateX2Norm!(x0::T, x1::T, x2_sum::T) where {T<:AbstractFloat}
    if (x2_sum + x1^2 - x0^2) < 0
        println("x2_sum $x2_sum x1 $x1 x0 $x0")
    end
    return max(x2_sum + x1^2 - x0^2, zero(T))
end
=#
function updateX2Norm!(x0::T, x1::T, x2_sum::T, X₁::Vector{T}, n::Int64) where {T<:AbstractFloat}
    x1_sq = x1^2
    x0_sq = x0^2
    #Can have numerical instability for small values of lambda. May need to recalculate the x2_sum. 
    if (x2_sum + x1_sq - x0_sq < 0)
        x2_sum_old = x2_sum
        x2_sum = x2_sum + x1_sq - x0_sq
        x2_sum_new = x2_sum
        x2_sum = zero(T)
        @turbo for col in 1:n
            x2_sum += X₁[col]*X₁[col]
        end
        println("x2_sum $x2_sum x2_sum_old $x2_sum_old x2_sum_new $x2_sum_new x1 $x1 x0 $x0 x1_sq $x1_sq x0_sq $x0_sq $n n")
        if n < 10
            println("X₁ ", X₁[1:n])
        end
        return x2_sum
    end
    return x2_sum + x1_sq - x0_sq
end

function newton_bisection!(Hs::SparseArray{Ti, T}, 
                            r::Vector{T}, 
                            X₁::Vector{T}, 
                            x2_sum::T,
                            col::Int64, 
                            δ::T, 
                            λ::T,
                            max_iter_newton::Int64, 
                            max_iter_bisection::Int64,
                            accuracy_newton::T,
                            accuracy_bisection::T,
                            norm_type::L2Norm) where {Ti<:Integer,T<:AbstractFloat}
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
            L1, L2 = getDerivatives!(Hs, r, col, δ, λ, x2_sum, X₁[col])
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
            x2_sum = updateX2Norm!(X0, X₁[col], x2_sum, X₁, Hs.n)
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
            x2_sum = updateX2Norm!(X0, X₁[col], x2_sum, X₁, Hs.n)
            updateResiduals!(Hs, r, col, X₁[col], X0)
            L1 = getL1(Hs, r, col, δ, λ, x2_sum, X₁[col])
            if sign(L1) != 1 #Otherwise the minimum is at X₁[col] < 0, so set X₁[col] == 0
                _, x2_sum = bisection!(Hs, r, X₁, x2_sum, col, δ, λ, zero(T), 
                            min(max(max_x1, zero(Float32)), Float32(1e11)), #Maximum Plausible Value for X1
                            L1,  
                            max_iter_bisection, #Should never reach this. Convergence in (max_x1)/2^n
                            accuracy_bisection,
                            norm_type)#accuracy)
            end
            return X₁[col] - X_init, x2_sum
        else #Convergence reached. Return difference between current estimate and initial guess. 
            return X₁[col] - X_init, x2_sum
        end
    end
end

function bisection!(Hs::SparseArray{Ti, T}, 
                    r::Vector{T}, 
                    X₁::Vector{T}, 
                    x2_sum::T,
                    col::Int64, 
                    δ::T, 
                    λ::T, 
                    a::T, 
                    b::T, 
                    fa::T,
                    max_iter::Int64, 
                    accuracy_bisection::T,
                    ::L2Norm) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    c = (a + b)/2
    #Since first guess for X₁[col] will change to "c"
    #Need to update the residuals accordingly. 
    updateResiduals!(Hs, r, col, c, X₁[col])
    X0 = X₁[col]
    X₁[col] = c
    x2_sum = updateX2Norm!(X0, X₁[col], x2_sum, X₁, Hs.n)
    X_init, X0 = X₁[col],  X₁[col]
    while (n < max_iter)

        #Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ, λ, x2_sum, X₁[col])
        #Bisection Rule
        if (sign(fc) != sign(fa))
            b, fb = c, fc
        else
            a, fa = c, fc
        end

        c, X0 = (a + b)/2, X₁[col]
        X₁[col] = c
        x2_sum = updateX2Norm!(X0, X₁[col], x2_sum, X₁, Hs.n)
        #Update residuals given new estimate, X₁[col], and prior estimate, X0
        updateResiduals!(Hs, r, col, X₁[col], X0)

        ########
        #Stopping Criterion
        ########
        abs(X₁[col] - X0) < accuracy_bisection ? break : nothing
        n += 1
    end
    return X₁[col] - X_init, x2_sum
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
                        norm_type::L2Norm) where {Ti<:Integer,T<:AbstractFloat,U<:Real}
    ΔX = Inf
    #sum of squared weight
    x2_sum = zero(T)
    @turbo for col in 1:Hs.n
        x2_sum += X₁[col]*X₁[col]
    end
    λ = Float32(1.0)
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
            δx, x2_sum = newton_bisection!(Hs, r, X₁, x2_sum, col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection,
                                        norm_type)
            δx = abs(δx)
            if !iszero(X₁[col]) 
                if δx/X₁[col] > max_diff
                    _diff =  δx/X₁[col]
                end
            end
            
            ΔX += δx
        end
        x2_sum = zero(T)
        @turbo for col in 1:Hs.n
            x2_sum += X₁[col]*X₁[col]
        end
        if _diff < max_diff
            break
        end
       i += 1
    end
    return i
end

