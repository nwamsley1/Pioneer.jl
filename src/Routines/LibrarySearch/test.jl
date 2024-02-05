function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
    @inbounds @simd for i in range(1, length(y))
        ϕ = Distributions.pdf(Normal(), η[i])
        μ = Distributions.cdf(Normal(), η[i])

        ϕ2 = Distributions.logpdf(Normal(), η[i])
        μ2 = Distributions.logcdf(Normal(), η[i])
        Z[i] = η[i] + (y[i] - μ)/ϕ
        #W[i] = (ϕ^2)/(μ*(1 - μ))
        W[i] = exp(2*ϕ2 - ( (μ2 - 2*μ2)))
    end
end

function fillXWX!(XWX::Matrix{T}, X::Matrix{T}, W::Vector{T}) where {T<:AbstractFloat}
    for i in range(1, size(XWX, 1))
        Threads.@threads for j in range(1, size(XWX, 1))
            @inbounds @fastmath for n in range(1, size(X, 1))
                XWX[i, j] += X[n,i]*W[i]*X[n,j]
            end
        end
    end
end



function fillY!(Y::Vector{T}, X::Matrix{T}, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat}
        Threads.@threads for col in range(1, size(X, 2))
            Y[col] = zero(T)
            @inbounds @fastmath for row in range(1, size(X, 1))
                Y[col] += X[row, col]*W[row]*Z[row]
            end
        end
end

function fillη!(η::Vector{T}, X::Matrix{T}, β::Vector{T}) where {T<:AbstractFloat}
    @turbo for row in range(1, size(X, 1))
        η[row] = zero(T)
    end
    @turbo for col in range(1, size(X, 2))
        for row in range(1, size(X, 1))
                η[row] += X[row, col]*β[col]
        end
    end
end

function ProbitRegression(β::Vector{T}, X::Matrix{T}, y::Vector{T}; max_iter::Int = 3) where {T<:AbstractFloat}
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
        println("max(W) ", maximum(W))
        println("maximum(Z) ", maximum(Z))
        println("Y $Y")
        println(" β $β")
        β = T.(XWX\Y)
        #println(LinearAlgebra.norm2(β .- old_β)/(LinearAlgebra.norm2(β)))
        old_β = copy(β)
    end
    return β
end