
function factorSpectrum(Wnew::Matrix{T}, Wold::Matrix{T}, HHt_diag::Vector{T}, WxHHt_VHt::Matrix{T}, HHt::SparseMatrixCSC{T, Int64}, λs::Vector{T}, max_iter::Int, tol::T) where {T<:AbstractFloat}
    #Fast Coordinate Descent Methods with Variable Selection for Non-negative Matrix Factorization
    #Cho-Jui Hsieh and Inderjit S. Dhillon
    a = Inf
    i = 1
    while (abs(a) > tol) & (i < max_iter)
        a = 0
        for r in 1:length(Wnew)
            
            #if iszero(Wnew[r]) #Is this legitimate? Several X speed boost from ignoring weights after they are set to zero. 
            #    continue
            #end
            Wnew[r] = max(zero(T), Wold[r] - (WxHHt_VHt[r] + (λs[r]))/HHt_diag[r])
            for i in HHt.colptr[r]:(HHt.colptr[r+1] -1)
                WxHHt_VHt[HHt.rowval[i]] += HHt.nzval[i]*(Wnew[r] - Wold[r])
            end
            a += abs(Wnew[r] - Wold[r])
            Wold[r] = Wnew[r]
        end
        i += 1
    end
end

function sparseNMF(H::SparseMatrixCSC{T, Int64}, X::Vector{T}, λ::T, γ::T, regularize::Bool = true; max_iter::Int = 1000, tol::T = 100*one(T)) where {T<:AbstractFloat}

    Wnew = 100*ones(T, (1, H.n))
    #Wnew = 100*ones(T, H.n)
    Wold = copy(Wnew)
    λs = zeros(T, H.n)
    #initW!(Wnew, Ht, X)
    #Wnew, Wold = copy(W[:]), copy(W[:])

    HHt = H'*H#H*Ht
    HHt_diag = collect(diag(HHt))
    VHt = X'*H
    WxHHt_VHt = collect(Wnew*HHt - VHt)

    ##OLS estimate with non-negative constraint since penalties are zero 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);

    regularize ? nothing : return Wnew 
    #Set adaptive weights 
    setLambdas!(λs, Float32(λ*sqrt(H.m)), γ, Wnew)
    #setLambdas!(λs, Float32(λ), γ, Wnew)

    #Addaptive LASSO estimation 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);
    return Wnew

end

#sparseNMF(H::SparseMatrixCSC{T, Int64}, X::Vector{T}) where {T<:AbstractFloat}= sparseNMF(H, X, zero(T), zero(T), regularize = false)
#=function sparseNMF(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}; λ::T = zero(T), γ::T = zero(T)/2, max_iter::Int = 1000, tol::T = 100*one(T)) where {T<:AbstractFloat}

    Wnew = 100*ones(T, (1, H.m))
    Wold = copy(Wnew)
    λs = zeros(T, H.m)
    #initW!(Wnew, Ht, X)
    #Wnew, Wold = copy(W[:]), copy(W[:])

    HHt = H*Ht
    HHt_diag = collect(diag(HHt))
    VHt = X'*Ht
    WxHHt_VHt = collect(Wnew*HHt - VHt)

    ##OLS estimate with non-negative constraint since penalties are zero 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);

    #Set adaptive weights 
    setLambdas!(λs, Float32(λ*sqrt(H.m)), γ, Wnew)
    #setLambdas!(λs, Float32(λ), γ, Wnew)

    #Addaptive LASSO estimation 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);
    return Wnew

end=#

function setLambdas!(λs::Vector{T}, λ::T, γ::T, W::Matrix{T}, adaptive::Bool = true) where {T<:AbstractFloat}
    if adaptive
        for i in 1:length(λs)
            if !isinf(λ/(W[i]))
                λs[i] = λ/(W[i]^γ)
            end
        end
    else
        for i in 1:length(λs)
                λs[i] = λ
        end
    end
end

function getDerivatives!(Hs::SparseMatrixCSC{Float64, Int64}, r::Vector{Float64}, X₁::Vector{Float64}, col::Int64, δ::Float64, X0::Float32)
    L1 = zzero(Float64)
    L2 = zero(Float64)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        #Huber
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        #R = (1 + (rval/δ)^2)^(-1/2)
        RS = (1 + (rval/δ)^2)
        #Quake's Fast Inverse Square Root Algorighm
        #Different magic for float64. 
        R = RS
        int64 = reinterpret(UInt64, xₛ)
        int64 = 0x5fe6eb50c7b537a9 - int64 >> 1
        R = reinterpret(Float64, int64)
        R *= 1.5 - RS * 0.5 * R^2
        HSVAL_R = hsval * R

        L1 += HSVAL_R * rval
        L2 += (hsval) * HSVAL_R * ((R)^(2))
        #L0 += (δ^2)*((1/R) - 1)

        #=
        #R = Hs.nzval[i]*(1 + (r[Hs.rowval[i]]/δ)^2)^(-1/2)
        #L1 += r[Hs.rowval[i]]*((R))
        #L2 += ((R)^(3))/Hs.nzval[i]
        #r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        =#
    end
    return L0, L1, L2
end

function getDerivatives!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{Float32}, col::Int64, δ::Float32)
    #L0 = zero(Float32)
    L1 = zero(Float32)
    L2 = zero(Float32)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        #Huber
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        #R = (1 + (rval/δ)^2)^(-1/2)
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
        #L0 += (δ^2)*((1/R) - 1)

        #=
        #R = Hs.nzval[i]*(1 + (r[Hs.rowval[i]]/δ)^2)^(-1/2)
        #L1 += r[Hs.rowval[i]]*((R))
        #L2 += ((R)^(3))/Hs.nzval[i]
        #r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        =#
    end
    return L1, L2
end

function getL1(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{Float32}, col::Int64, δ::Float32)
    L1 = zero(Float32)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
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
    return L1
end

function getL1(Hs::SparseMatrixCSC{Float64, Int64}, r::Vector{Float64}, col::Int64, δ::Float64)
    L1 = zero(Float32)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
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
    return L1
end

function getL0!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{Float32}, X₁::Vector{Float32}, col::Int64, δ::Float32, X0::Float32)
    L0 = zero(Float32)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        rval = r[Hs.rowval[i]]
        RS = (1 + (rval/δ)^2)
        L0 += (δ^2)*(sqrt(RS) - 1)
    end
    return L0
end

function updateResiduals!(Hs::SparseMatrixCSC{T, Int64}, r::Vector{T}, col::Int64, X1, X0) where {T<:AbstractFloat}
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        r[Hs.rowval[i]] += Hs.nzval[i]*(X1 - X0)
    end
end

function newton_bisection!(Hs::SparseMatrixCSC{T, Int64}, r::Vector{T}, X₁::Vector{T}, col::Int64, δ::T; max_iter::Int64 = 20, accuracy::T = 0.01) where {T<:AbstractFloat}
    n = 0
    X_init = X₁[col] #Estimate prior to optimiztion
    X0 = X₁[col] #Keeps track of previous etimate at each iteration
    #Maximum Estimates X₁[col] and L1. Used as uper bound if newton-raphson fails to converge
    #And bisection method is neede. 
    max_l1, max_x1 = typemax(T), typemax(T)

    #Newton-Raphson Method Iterations until convergence or maximum iterations. 
    #If convergence fails in maximum iterations, switch to bisection method for guaranteed convergence
    while (n < max_iter)

        #First and second derivatives 
        L1, L2 = getDerivatives!(Hs, r, col, δ)
        update_rule = L1/L2

        #Switch to bisection method
        if isnan(update_rule)
            n = max_iter
            break
        end

        #Useful as boundaries for bisection method if Newton's method fails to converge. 
        #Want positive and negative L1 with smallest absolute values 
        if (sign(L1) == 1) & (L1 < max_l1) #L1 > max_l1
            max_x1, max_l1 = X₁[col], L1
        end

        X0 = X₁[col] 

        #Newton-Raphson update. Contrain X₁[col] to be non-negative
        @fastmath X₁[col] = max(X₁[col] - update_rule, zero(T))

        n += 1

        #Update residuals given new estimate, X₁[col], and prior estimate, X0
        updateResiduals!(Hs, r, col, X₁[col], X0)

        ########
        #Stopping Criterion for single variable Newton-Raphson
        ########
        abs(X₁[col] - X0) < accuracy ? break : nothing
    end

    #If newtons method fails to converge, switch to bisection method
    if n == max_iter
        #######
        #Lower bound is always X₁[col] == 0
        X0 = X₁[col]
        X₁[col] = zero(T)
        updateResiduals!(Hs, r, col, X₁[col], X0)
        L1 = getL1(Hs, r, col, δ)

        if sign(L1) != 1 #Otherwise the minimum is at X₁[col] < 0, so set X₁[col] == 0
            bisection!(Hs, r, X₁, col, δ, zero(T), 
                        min(max_x1, Float32(1e11)), #Maximum Plausible Value for X1
                        L1,  
                        max_l1, 
                        max_iter_inner = 100, #Should never reach this. Convergence in (max_x1)/2^n
                        accuracy = accuracy)
        end

        return X₁[col] - X_init
    else #Convergence reached. Return difference between current estimate and initial guess. 
        return X₁[col] - X_init
    end
end

function bisection!(Hs::SparseMatrixCSC{T, Int64}, r::Vector{T}, X₁::Vector{T}, col::Int64, δ::T, a::T, b::T, fa::T, fb::T; max_iter_inner::Int64 = 20, accuracy::T = 0.01) where {T<:AbstractFloat}
    n = 0
    c = (a + b)/2

    #Since first guess for X₁[col] will change to "c"
    #Need to update the residuals accordingly. 
    updateResiduals!(Hs, r, col, c, X₁[col])

    X₁[col] = c
    X_init, X0 = X₁[col],  X₁[col]

    while (n < max_iter_inner)

        #Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ)
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
        abs(X₁[col] - X0) < accuracy ? break : nothing
        n += 1
    end
    return X₁[col] - X_init
end

function solveHuber!(Hs::SparseMatrixCSC{T, Int64}, r::Vector{T}, X₁::Vector{T}, δ::T; max_iter_outer::Int = 1000, max_iter_inner::Int = 20, tol::U = 100) where {T<:AbstractFloat,U<:Real}
    ΔX = Inf
    i = 0
    while (i < max_iter_outer) & (ΔX > tol)
        ΔX = 0.0
        #Update each variable once 
        max_diff = 0.0
        for col in range(1, Hs.n)
            #δx = abs(newtonRaphson!(Hs, r, X₁, col, δ, max_iter_inner = (min(1 + i*2, max_iter_inner)), accuracy = T(100)))
            δx = abs(newton_bisection!(Hs, r, X₁, col, δ, max_iter = max_iter_inner, accuracy = T(100)))
            if δx/X₁[col] > max_diff
                max_diff =  δx/X₁[col]
            end
            ΔX += δx
        end

       i += 1
    end
    return i
end

#mean(X.*sum(Hs_mat.>0, dims = 1)[:])*200*200
#=
function newtonRaphson!(Hs::SparseMatrixCSC{T, Int64}, r::Vector{T}, X₁::Vector{T}, col::Int64, δ::T; max_iter_inner::Int64 = 20, accuracy::T = 0.01) where {T<:AbstractFloat}
    n = 0
    a = 0
    x = X₁[col]
    X0 = X₁[col]
    while (n < max_iter_inner)
        L1, L2 = getDerivatives!(Hs, r, X₁, col, δ, X0)
        #col  ∈ (1,2) ? println("col $col, L1 $L1, L2 $L2, X₁[col] ", X₁[col], "n $n") : nothing
        #abs(L1) > δ*0.90 ? println("col $col; L1 $L1") : nothing
        X0 = X₁[col] 
        #if (abs(L1) < δ*0.9) #& (X₁[col]>0.0)#& (abs(L2) > 1e-7)
            #col  ∈ (1,2) ? println("col $col, X₁[col] - 0.5*(L1)/(sqrt(n+1)*L2) ", X₁[col] - 0.5*(L1)/(sqrt(n+1)*L2)) : nothing
            @fastmath X₁[col] = max(X₁[col] - 0.5*(L1)/(sqrt(n+1)*L2), zero(T))
            #@fastmath X₁[col] = max(X₁[col] - 0.5*one(T)*(L1)/(L2), zero(T))
        #else
        #    L0 = getL0!(Hs, r, X₁, col, δ, X0)
        #    #col ∈ (1,2) ? println("col $col, L0 $L0, L0/L1 ", L0/L1, " X₁[col] - L0/(L1) ", X₁[col] - L0/(L1)) : nothing
        #    @fastmath X₁[col] = max(X₁[col] - L0/((a + 1)*L1), zero(T))
        #    a += 1
        #end
        n += 1
        #Shouldn't need to check for this. 
        #Seems to happen when L1 or L2 blow out to Inf .
        isnan(X₁[col]) ? X₁[col] = zero(T) : nothing

        ########
        #Stopping Criterion for single variable Newton-Raphson
        ########
        abs(X₁[col] - X0) < accuracy ? break : nothing
        
    end

    #Need to update residuals after the last iteration. 
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
    end

    if n == max_iter_inner
        return bisection(Hs, r, X₁, col, δ, min_max_l1, max_iter_inner = max_iter_inner, accuracy = accuracy)
    else
        return X₁[col]-x
    end
end
function solveHuber3!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₁::Vector{T}, δ::T; max_iter_outer::Int = 1000, max_iter_inner::Int = 20, tol::U = 100, λ::Float32 = zero(Float32)) where {T<:AbstractFloat,U<:Real}
    L1 = zero(T)
    L2 = zero(T)
    Hst = sparse(transpose(Hs))
    ΔX = T(Inf)

    #Update each variable once 
    for col in range(1, Hs.n)
        newtonRaphson!(Hs, r, X₁, col, δ, max_iter_inner = max_iter_inner, accuracy = T(100), λ = λ)
    end

    #Setup greedy coordinate descnet 
    max_gradient = MutableBinaryMaxHeap([Inf for x in 1:Hs.n])
    for col in range(1, Hs.n)
        L1 = getGradient!(Hs, r, col, δ)
        update!(max_gradient, col, abs(L1))
    end
    #gradients = Float32[]
    #diffs = Float32[]
    #ns = Int64[]
    i = 0
    
    while (i < max_iter_outer) & (ΔX > tol)
        L1,L2 = zero(T), zero(T)
        val, col = top_with_handle(max_gradient) #Coefficient with the largest gradient. 
        #println("col $col; val $val; X₁[col] ", X₁[col])
        #push!(gradients, val)
        #println("val is $val and col is $col and i is $i")

        ########
        #Newton-Raphson to optimize w.r.t X₁[col]. 
        ########
        #X0 = X₁[col]
        n = newtonRaphson!(Hs, r, X₁, col, δ, max_iter_inner = max_iter_inner, accuracy = T(10), λ = λ)
        #=push!(ns, n)
        if !iszero(X₁[col])
            push!(diffs, abs(X0 - X₁[col])/(X₁[col] + 1000))
            #push!(diffs, abs(X0 - X₁[col]))
        end=#

        ########
        #Which gradients need to be updated?
        ########
        gradients_to_update = getGradientsToUpdate!(Hs, Hst, col);

        ########
        #calculate any gradients that have changed
        ########
        for column in gradients_to_update
            #println("column $column")
            L1 = getGradient!(Hs, r,  column, δ)
            #=if column == 100
                println("col is $column and L1 is $L1")
            end=#
            iszero(X₁[ column]) & (sign(L1) == 1) ? update!(max_gradient, column, 0.0) :  update!(max_gradient, column, abs(L1))
        end
        i += 1
    end
    #return gradients, diffs, ns
    return i
end

function solveHuber2!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₁::Vector{T}, δ::T; max_iter_outer::Int = 1000, max_iter_inner::Int64 = 20, tol::U = 100, λ::Float32 = zero(Float32)) where {T<:AbstractFloat,U<:Real}
    L1 = zero(T)
    L2 = zero(T)
    ΔX = T(Inf)
    #X₁ .= 1000*ones(T, length(X₁))
    #r = Hs*X₁ .- b
    i = 0
    while (i < max_iter_outer) & (ΔX > tol)
        #Sum of differences bewteen i+1 and i'th iteration
        ΔX = 0.0
        #Loop through single variable optimizations
        for col in range(1, Hs.n)
            L1,L2 = zero(T), zero(T)
            n = 0
            X0 = X₁[col]
            X₀ = X0
            ########
            #Newton-Raphson to optimize w.r.t X₁[col]. 
            ########
            while (n < max_iter_inner)
                X0 = X₁[col]
                L1 = zero(T)
                L2 = zero(T)
                @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)

                    #Least-Squares 
                    #L1 += Hs.nzval[i]*r[Hs.rowval[i]]
                    #L2 += Hs.nzval[i]^2

                    #Huber
                    R = (r[Hs.rowval[i]]/δ)^2
                    L1 += Hs.nzval[i]*r[Hs.rowval[i]]*((1 + R)^(-1/2))
                    L2 += (Hs.nzval[i]^2)/((1 + R)^(3/2))
                end
                #Apply Damping Factor
                X₁[col] = max(X₁[col] - 0.5*(L1+λ)/L2, 0.0)
                #No Damping Factor
                #X₁[col] = max(X₁[col] - (L1+λ)/L2, 0.0)


                #Update residuals 
                @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                    r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
                end

                ########
                #Stopping Criterion for single variable Newton-Raphson
                #Accuracy requirement increases each outer-loop
                #Need to find a generally acceptable parameter value
                ########
                if abs((X₁[col]-X0)/X0) < 0.1/i
                    break
                end

                if X₁[col]==0.0
                    break
                end
                n += 1
            end
            ΔX += abs(X₁[col]-X₀)
        end
        i += 1
    end
    return i
end
function newtonRaphson!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₁::Vector{T}, col::Int64, δ::T; max_iter_inner::Int64 = 20, accuracy::T = 0.01, λ::Float32 = zero(Float32)) where {T<:AbstractFloat}
    n = 0
    x = X₁[col]
    X0 = X₁[col]
    δX = T(0.0)
    while (n < max_iter_inner)
        L1 = zero(T)
        L2 = zero(T)
        
        @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            r[Hs.rowval[i]] += Hs.nzval[i]*(δX)

            #Huber
            rval = r[Hs.rowval[i]]
            hsval = Hs.nzval[i]
            #R = (1 + (rval/δ)^2)^(-1/2)
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

            #R = Hs.nzval[i]*(1 + (r[Hs.rowval[i]]/δ)^2)^(-1/2)
            #L1 += r[Hs.rowval[i]]*((R))
            #L2 += ((R)^(3))/Hs.nzval[i]
            #r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        end

        #Newton-Raphson update with damping factor
        #X₁[col] = max(X₁[col] - (0.5/(n + 1))*(L1+λ)/(L2), 0.0)
        X0 = X₁[col]
        @fastmath X₁[col] = max(X₁[col] - one(T)*(L1+λ)/(sqrt(n + 1)*L2), zero(T))
        δX = X₁[col] - X0
        #Update residuals 
        #@turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #    r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
        #end
    
        ########
        #Stopping Criterion for single variable Newton-Raphson
        #Accuracy requirement increases each outer-loop
        #Need to find a generally acceptable parameter value
        ########
        #if (abs((X₁[col]-X0)/(X0)) < accuracy)
        if abs(δX) < accuracy
            break
        end
        #=if iszero(X₁[col]) & (sign(L1) == 1)
            break
        end=#
        n += 1
    end
    return X₁[col]-x
end



function getGradient!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, col::Int64, δ::T) where {T<:AbstractFloat}
    L1 = zero(T)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        R = (r[Hs.rowval[i]]/δ)^2
        L1 += Hs.nzval[i]*r[Hs.rowval[i]]*((1 + R)^(-1/2))
    end
    return L1
end

function getGradientsToUpdate!(Hs::SparseMatrixCSC{Float32, Int64}, Hst::SparseMatrixCSC{Float32, Int64}, col::Int64)     
    gradients_to_update = Set{Int64}()  
    for row in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #Get all non-zero elements in the row
        for i in Hst.colptr[Hs.rowval[row]]:(Hst.colptr[Hs.rowval[row] + 1] - 1)
            #This gradient needs to be updated. 
            if i != col
                push!(gradients_to_update, Hst.rowval[i])
            end
        end
    end
    return gradients_to_update
end

function firstOrder!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₁::Vector{T}, col::Int64, δ::T; max_iter_inner::Int64 = 20, accuracy::T = 0.01, λ::Float32 = zero(Float32)) where {T<:AbstractFloat}

    X0 = X₁[col]
    L0 = zero(T)
    L1 = zero(T)
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        #R = Hs.nzval[i]*(1 + (r[Hs.rowval[i]]/δ)^2)^(-1/2)
        R = (1 + (r[Hs.rowval[i]]/δ)^2)^(1/2)
        L0 += (δ^2)*(R - 1)
        L1 += r[Hs.rowval[i]]*Hs.nzval[i]/R
        #r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
    end

    @fastmath X₁[col] = max(X₁[col] - L0/L1, zero(T))
    @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
    end

    return
end


function factorSpectrum(Wnew::Vector{T}, Wold::Vector{T}, HHt_diag::Vector{T}, WxHHt_VHt::Matrix{T}, HHt::SparseMatrixCSC{T, Int64}, λ::T, max_iter::Int, tol::T) where {T<:AbstractFloat}
    a = Inf
    i = 1
    while (abs(a) > tol) & (i < max_iter)
        a = 0
        for r in 1:length(Wnew)
            Wnew[r] = max(zero(T), Wold[r] - (WxHHt_VHt[r] + λ)/HHt_diag[r])
            for i in HHt.colptr[r]:(HHt.colptr[r+1] -1)
                WxHHt_VHt[HHt.rowval[i]] += HHt.nzval[i]*(Wnew[r] - Wold[r])
            end
            a += abs(Wnew[r] - Wold[r])
            Wold[r] = Wnew[r]
        end
        i += 1
    end
end

non_zeros = Int64[]
test_lambdas = Float64[]
t_λ = 1e7
for i in 1:150 #λ in range(1e4, 3e6, step = 1.5e4)
    t_λ = t_λ*1.1
    testW =  sparseNMF(Hs, Hst, X, λ=Float32(t_λ), max_iter = 1000)
    push!(non_zeros, sum(testW.>100))
    push!(test_lambdas, t_λ)
end
plot(log.(test_lambdas), non_zeros, seriestype=:scatter)

function sparseNMF(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}, λ::T = zero(T), max_iter::Int = 1000, tol::T = 10*one(T)) where {T<:AbstractFloat}
    W = [abs(x) for x in randn(T, (1, H.m))] 
    Wnew, Wold = copy(W[:]), copy(W[:])
    #Initialize
    #HHt = H*H'
    HHt = H*Ht
    HHt_diag = collect(diag(HHt))
    VHt = X'*Ht
    WxHHt_VHt = collect(W*HHt - VHt)
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λ, max_iter, tol);
    return Wnew
end

non_zeros = Int64[]
test_lambdas = Float64[]
t_λ = 1e4
for i in 1:100 #λ in range(1e4, 3e6, step = 1.5e4)
    t_λ = t_λ*1.1
    testW = sparseNMF(Hs, X, λ=Float32(t_λ));
    push!(non_zeros, sum(testW.<100))
    push!(test_lambdas, t_λ)
end
plot(log2.(test_lambdas), non_zeros, seriestype=:scatter)
=#
######
#Get seperate regulrization coefficients for each template. Termed "Adaptive LASSO"
######

#=function initW!(W::Matrix{T}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}) where {T<:AbstractFloat}
    for col in 1:(length(Ht.colptr)-1)
        n = 0
        for row in Ht.colptr[col]:(Ht.colptr[col+1] - 1)
            W[col] += X[Ht.rowval[row]]/Ht.nzval[row]
            n += 1
        end
        W[col] = W[col]/n
    end
end

sparseNMF(Hs, Hst, X, λ=one(Float32), max_iter = 1000)



norm_facs = zeros(Float32, 200)
for r in 1:(length(H_peak_pep.colptr)-1)
    n = 1
    for i in H_peak_pep.colptr[r]:(H_peak_pep.colptr[r+1] -1)
        norm_facs[r] += X[H_peak_pep.rowval[i]]
        n += 1
    end
    norm_facs[r]/n
end

Hs_mat_half = Hs_mat[1:100,:]
@benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
tol = 100, #Need a reasonable way to choose lambda?
update_H = false, #Important to keep H constant. 
regularization=:both,
l₁ratio = 1.0,
α = 1e5,
shuffle=true
), X_, W, Hs_mat).W)

W_half = W[:,1:100]
@benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
tol = 100, #Need a reasonable way to choose lambda?
update_H = false, #Important to keep H constant. 
regularization=:both,
l₁ratio = 1.0,
α = 1e5,
shuffle=true
), X_, W_half, Hs_mat_half).W[1,:])
#=@btime sparseNMF(Hs, X, λ=Float32(1e5))

@profview for i in 1:100
    sparseNMF(Hs, X, λ=Float32(1e5))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/H.jld2"  H

Hs_mat = Matrix(HS)
X_ = Matrix(X')
weights = (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
    tol = 100, #Need a reasonable way to choose lambda?
    update_H = false, #Important to keep H constant. 
    regularization=:both,
    l₁ratio = 1.0,
    α = 1e5
    ), X_, W, Hs_mat).W[1,:])

    @benchmark (NMF.solve!(NMF.GreedyCD{Float32}(maxiter=1000, verbose = false, 
           tol = 100, #Need a reasonable way to choose lambda?
           update_H = false, #Important to keep H constant. 
           lambda_w = 1e5,
           lambda_h = 1e5
           ), X_, W, Hs_mat).W[1,:])

    @benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
           tol = 100, #Need a reasonable way to choose lambda?
           update_H = false, #Important to keep H constant. 
           regularization=:both,
           l₁ratio = 1.0,
           α = 1e5,
           shuffle=true
           ), X_, W, Hs_mat).W[1,:])
           =#


=#