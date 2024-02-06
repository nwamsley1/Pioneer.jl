function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::Vector{Bool},data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T<:AbstractFloat}
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk#range(1, length(y))
                @inbounds @fastmath begin
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
    end
    fetch(tasks)
end

function fillXWX!(XWX::Matrix{T}, X::DataFrame, W::Vector{T}) where {T<:AbstractFloat}
    @noinline function fillCell(XWX::T, Xi::Vector{R}, Xj::Vector{V}, W::Vector{T}) where {T<:AbstractFloat, R,V<:Real}
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
            XWX[i, j] = fillCell(zero(T), X[!, i], X[!, j], W)
        end
    end
end

function fillY!(Y::Vector{T}, X::DataFrame, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat}

        @noinline function fillCell(Y::T, X::Vector{R}, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat,R<:Real}
            @turbo for i in range(1, length(W))
                Y += X[i]*W[i]*Z[i]
            end
            return Y
        end

        Threads.@threads for col in range(1, size(X, 2))
            Y[col] = zero(T)
            Y[col] = fillCell(Y[col], X[!,col], W, Z)
        end
end

function fillη!(η::Vector{T}, X::DataFrame, β::Vector{T}, bounds::Tuple{T, T}) where {T}


    function fillColumn!(η::Vector{T}, X::Vector{R}, β::T) where {T<:AbstractFloat,R<:Real}
        @turbo for row in range(1, length(η))
            η[row] += X[row]*β
        end
    end

    @turbo for row in range(1, size(X, 1))
        η[row] = zero(T)
    end

    for col in range(1, size(X, 2))
        fillColumn!(η, X[!,col], β[col])
    end

    for row in range(1, size(X, 1))
        η[row] = max(min(η[row], last(bounds)), first(bounds))
    end
end

function ProbitRegression(β::Vector{T}, X::DataFrame, y::Vector{Bool}; max_iter::Int = 30) where {T<:AbstractFloat}
    @time begin
    W = zeros(T, length(y))
    η = zeros(T, length(y))
    Z = zeros(T, length(y))
    Y = zeros(T, size(X, 2))
    XWX = zeros(T, (size(X, 2), size(X, 2)))
    old_β = copy(β)
    old_loss = 0.0

    tasks_per_thread = 10
    chunk_size = max(1, size(X, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(X, 1), chunk_size) # partition your data into chunks that
    end

    @time for i in range(1, max_iter)#ProgressBar(range(1, max_iter))
        fillη!(η, X, β, (-8.0, 8.0))        
        fillZandW!(Z, W, η, y, data_chunks)
        loss = 0.0
        @turbo for i in range(1, length(η))
            loss += η[i]
        end
        fillXWX!(XWX, X, W)
        fillY!(Y, X, W, Z)
        β = T.(XWX\Y)
        if abs(1 - exp(loss - old_loss)) < 1e-2
            break
        end
        old_loss = loss
        old_β = copy(β)
    end
    return β
end

function ModelPredict!(psms::DataFrame, β::Vector{T}) where {T<:AbstractFloat}
    scores = zeros(Float32, size(psms, 1))

    function fillColumn!(score::Vector{T}, X::Vector{R}, β::U) where {T,U<:AbstractFloat,R<:Real}
        @turbo for row in range(1, length(score))
            score[row] += X[row]*β
        end
    end

    for col in range(1, length(β))
        fillColumn!(scores, psms[!,col], β[col])
    end

    return scores
    
end