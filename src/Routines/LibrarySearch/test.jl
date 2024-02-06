function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::BitVector) where {T<:AbstractFloat}
    Threads.@threads for i in range(1, length(y))
        ϕ = Distributions.pdf(Normal(), η[i]) #+ 1e-6
        μ = Distributions.cdf(Normal(), η[i]) #+ 1e-6

        #ϕ = Distributions.pdf(Normal(), η)
        #μ = Distributions.cdf(Normal(), η)
        #ϕ2 = Distributions.logpdf(Normal(), η[i])
        #μ2 = Distributions.logcdf(Normal(), η[i])
        #ϕ2 = Distributions.logpdf(Normal(), η[i])
        #μ2 = Distributions.logcdf(Normal(), η[i])
        if y[i]
            #Z[i] = η[i] + exp(-ϕ2) - exp(μ2 - ϕ2) #(1 - μ)/(ϕ + 1e-6)
            Z[i] = η[i] + (1 - μ)/ϕ
        else
            Z[i] = η[i] - μ/ϕ
        end
        #if isinf(Z[i])
        #    Z[i] = sign(Z)*T(10.0)
        #elseif isnan(Z[i])
        #end
        #Z[i] = max(min(Z[i], T(10.0)), T(-10.0))
        #Z[i] = min(η[i] + exp(log(y[i] + 1e-6) - ϕ2) - exp(μ2 - ϕ2), 10.0)
        W[i] = (ϕ^2)/(μ*(1 - μ))
        
        #W[i] = exp(2*ϕ2 - ( (μ2 - 2*μ2)))
    end
end

function fillXWX!(XWX::Matrix{T}, X::Matrix{U}, W::Vector{T}) where {T,U<:AbstractFloat}
    for i in range(1, size(XWX, 1))
        Threads.@threads for j in range(1, size(XWX, 1))
            @inbounds @fastmath for n in range(1, size(X, 1))
                XWX[i, j] += T(X[n,i])*W[i]*T(X[n,j])
            end
        end
    end
end



function fillY!(Y::Vector{T}, X::Matrix{U}, W::Vector{T}, Z::Vector{T}) where {T,U<:AbstractFloat}
        Threads.@threads for col in range(1, size(X, 2))
            Y[col] = zero(T)
            for row in range(1, size(X, 1))
                Y[col] += T(X[row, col])*W[row]*Z[row]
            end
        end
end

function fillη!(η::Vector{T}, X::Matrix{U}, β::Vector{T}) where {T,U<:AbstractFloat}
    @turbo for row in range(1, size(X, 1))
        η[row] = zero(T)
    end
    for col in range(1, size(X, 2))
        for row in range(1, size(X, 1))
                η[row] += T(X[row, col])*β[col]
        end
    end
    #for row in range(1, size(X, 1))
    #    η[row] = max(min(η[row], T(3.0)), T(-3.0))
    #end
end

function ProbitRegression(β::Vector{T}, X::Matrix{U}, y::BitVector; max_iter::Int = 3) where {T,U<:AbstractFloat}
    W = zeros(T, length(y))
    η = zeros(T, length(y))
    Z = zeros(T, length(y))
    Y = zeros(T, size(X, 2))
    XWX = zeros(T, (size(X, 2), size(X, 2)))
    old_β = copy(β)
    @time for i in range(1, max_iter)#ProgressBar(range(1, max_iter))
        fillη!(η, X, β)
        fillZandW!(Z, W, η, y)
        fillXWX!(XWX, X, W)
        fillY!(Y, X, W, Z)
        println("max(W) ", maximum(W), " ", minimum(W))
        println("maximum(Z) ", maximum(Z)," ", minimum(Z))
        println("maximum(η) ", maximum(η), " ", minimum(Z))
        println("Y", maximum(Y), " ", minimum(Y))
        #println("Y $Y")
        #println(" β $β")
        β = T.(XWX\Y)
        #println(LinearAlgebra.norm2(β .- old_β)/(LinearAlgebra.norm2(β)))
        old_β = copy(β)
    end
    return β
end