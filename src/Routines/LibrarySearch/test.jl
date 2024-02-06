function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::BitVector) where {T<:AbstractFloat}
    @inbounds @fastmath begin
        Threads.@threads for i in range(1, length(y))
            ϕ = exp(-(η[i]^2)/2)/sqrt(2*π)
            μ = (1 + SpecialFunctions.erf(η[i]/sqrt(2)))/2
            if y[i]
                Z[i] = η[i] + (1 - μ)/ϕ
                η[i] = log(μ)
            else
                Z[i] = η[i] - μ/ϕ
                η[i] = 1- log(μ)
            end
            W[i] = (ϕ^2)/(μ*(1 - μ))

        end
    end
end

function fillXWX!(XWX::Matrix{T}, X::Matrix{U}, W::Vector{T}) where {T,U<:AbstractFloat}
    @noinline function fillCell(XWX::T, Xi::AbstractArray{U}, Xj::AbstractArray{U}, W::Vector{T}) where {T,U<:AbstractFloat}
        #N = 8
        #lane = VecRange{N}(0)
        #@inbounds for i in 1:N:length(W)
        #    XWX += Xi[lane + i]*W[lane + i]*Xj[lane + i]
        #end
        @turbo for i in 1:length(W)
            XWX += Xi[i]*W[i]*Xj[i]
        end
        return XWX
    end
    for i in range(1, size(XWX, 1))
        Threads.@threads for j in range(1, size(XWX, 1))
            XWX[i, j] = fillCell(zero(T), @view(X[:, i]),@view(X[:, j]), W)
        end
    end
end



function fillY!(Y::Vector{T}, X::Matrix{U}, W::Vector{T}, Z::Vector{T}) where {T,U<:AbstractFloat}

        @noinline function fillCell(Y::T, X::AbstractArray{U}, W::Vector{T}, Z::Vector{T}) where {T,U<:AbstractFloat}
            @turbo for i in range(1, length(W))
                Y += X[i]*W[i]*Z[i]
            end
            return Y
        end

        Threads.@threads for col in range(1, size(X, 2))
            Y[col] = zero(T)
            Y[col] = fillCell(Y[col], @view(X[:,col]), W, Z)
        end
end

function fillη!(η::Vector{T}, X::Matrix{U}, β::Vector{T}) where {T,U<:AbstractFloat}
    @turbo for row in range(1, size(X, 1))
        η[row] = zero(T)
    end
    for col in range(1, size(X, 2))
        @turbo for row in range(1, size(X, 1))
                η[row] += X[row, col]*β[col]
        end
    end
    for row in range(1, size(X, 1))
        η[row] = max(min(η[row], T(8.0)), T(-8.0))
    end
end

function ProbitRegression(β::Vector{T}, X::Matrix{U}, y::BitVector; max_iter::Int = 30) where {T,U<:AbstractFloat}
    W = zeros(T, length(y))
    η = zeros(T, length(y))
    Z = zeros(T, length(y))
    Y = zeros(T, size(X, 2))
    XWX = zeros(T, (size(X, 2), size(X, 2)))
    old_β = copy(β)
    old_loss = 0.0
    @time for i in range(1, max_iter)#ProgressBar(range(1, max_iter))
        fillη!(η, X, β)        
        fillZandW!(Z, W, η, y)
        loss = 0.0
        @turbo for i in range(1, length(η))
            loss += η[i]
        end
        fillXWX!(XWX, X, W)
        fillY!(Y, X, W, Z)
        β = T.(XWX\Y)
        if abs(1 - exp(loss - old_loss)) < 1e-2
            return β
        end
        old_loss = loss

        #println(LinearAlgebra.norm2(β .- old_β)/(LinearAlgebra.norm2(β)))
        old_β = copy(β)
    end
    return β
end