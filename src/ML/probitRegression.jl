#y = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\y.jld2")["y"]
#X = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\X.jld2")["X"]
#PSMs = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\PSMs.jld2")["PSMs"]
function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::Vector{Bool},data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T<:AbstractFloat}
    tasks = map(data_chunks) do chunk
    #@inbounds @fastmath begin
    Threads.@spawn begin
        @inbounds @fastmath for i in chunk
            #@inbounds @fastmath for i in chunk#range(1, length(y))
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
    fetch.(tasks)
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
        tasks = map(range(1, size(X, 2))) do j
            Threads.@spawn begin #for j in range(1, size(XWX, 1))
                XWX[i, j] = fillCell(zero(T), X[!, i], X[!, j], W)
            end
        end
        fetch.(tasks)
    end
end

function fillY!(Y::Vector{T}, X::DataFrame, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat}

        @noinline function fillCell(Y::T, X::Vector{R}, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat,R<:Real}
            @turbo for i in range(1, length(W))
                Y += X[i]*W[i]*Z[i]
            end
            return Y
        end

        tasks = map(range(1, size(X, 2))) do col
            Threads.@spawn begin
                Y[col] = zero(T)
                Y[col] = fillCell(Y[col], X[!,col], W, Z)
            end
        end

        fetch.(tasks)
end

function fillη!(η::Vector{T}, X::DataFrame, β::Vector{T}, bounds::Tuple{T, T}, data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T}


    function fillColumn!(η::Vector{T}, X::Vector{R}, β::T, data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T<:AbstractFloat,R<:Real}
        tasks = map(data_chunks) do row_chunk
            Threads.@spawn begin
                @turbo for row in row_chunk
                    η[row] += X[row]*β
                end
            end
        end
        fetch.(tasks)
    end

    #Initialize Z-scores to zero
    tasks = map(data_chunks) do row_chunk
        Threads.@spawn begin
            @turbo for row in row_chunk
                η[row] = zero(T)
            end
        end
    end
    fetch.(tasks)

    #Calcualte Z-scores as X*β
    for col in range(1, size(X, 2))
        fillColumn!(η, X[!,col], β[col], data_chunks)
    end

    #Enforce boundary constraint (to avoid Infs and NaNs in subsequent calcs)
    tasks = map(data_chunks) do row_chunk
        Threads.@spawn begin
            for row in row_chunk
                η[row] = max(min(η[row], last(bounds)), first(bounds))
            end
        end
    end
    fetch.(tasks)
    
end

function vecSum!(v::Vector{T}, data_chunks) where {T<:AbstractFloat}
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            vsum = zero(T)
            @turbo for i in chunk
                vsum += v[i]
            end
            return vsum
        end
    end
    return sum(fetch.(tasks))
end

function ProbitRegression(β::Vector{T}, X::DataFrame, y::Vector{Bool},
                            data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}; 
                            max_iter::Int = 30, z_score_bounds::Tuple{Float64, Float64} = (-8.0, 8.0)
                            ) where {T<:AbstractFloat}

    W, η, Z = zeros(T, length(y)), zeros(T, length(y)), zeros(T, length(y))
    Y = zeros(T, size(X, 2))
    XWX = zeros(T, (size(X, 2), size(X, 2)))
    old_loss = 0.0

    for i in range(1, max_iter)#ProgressBar(range(1, max_iter))
        fillη!(η, X, β, z_score_bounds, data_chunks)        
        fillZandW!(Z, W, η, y, data_chunks)
        loss = vecSum!(η, data_chunks)
        fillXWX!(XWX, X, W)
        fillY!(Y, X, W, Z)
        β = T.(XWX\Y)
        if abs(1 - exp(loss - old_loss)) < 1e-2
            break
        end
        old_loss = loss
    end
    return β

end

function ModelPredict!(scores::Vector{U}, 
                        psms::DataFrame,
                        β::Vector{T}, 
                        data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T,U<:AbstractFloat}
    
    function fillColumn!(scores::Vector{T}, X::Vector{R}, β::U, data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T,U<:AbstractFloat,R<:Real}
        tasks = map(data_chunks) do row_chunk
            Threads.@spawn begin
                @turbo for row in row_chunk
                    scores[row] += X[row]*β
                end
            end
        end
        fetch.(tasks)
    end

    fill!(scores, zero(U));
    for col in range(1, size(psms, 2))
        fillColumn!(scores, psms[!,col], β[col], data_chunks)
    end

end

function ModelPredictCVFold!(scores::Vector{U}, 
                        cv_folds::Vector{UInt8},
                        cv_fold::UInt8,
                        psms::DataFrame,
                        β::Vector{T}, 
                        data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T,U<:AbstractFloat}
    
    function fillColumn!(scores::Vector{T},   
                            cv_folds::Vector{UInt8},
                            cv_fold::UInt8,
                            X::Vector{R}, 
                            β::U, 
                            data_chunks::Base.Iterators.PartitionIterator{UnitRange{Int64}}) where {T,U<:AbstractFloat,R<:Real}
        tasks = map(data_chunks) do row_chunk
            Threads.@spawn begin
                for row in row_chunk
                    if cv_folds[row] != cv_fold
                        scores[row] += X[row]*β
                    end
                end
            end
        end
        fetch.(tasks)
    end

    for col in range(1, size(psms, 2))
        fillColumn!(scores, cv_folds, cv_fold, psms[!,col], β[col], data_chunks)
    end

    tasks = map(data_chunks) do row_chunk
        Threads.@spawn begin
            for row in row_chunk
                if cv_folds[row] != cv_fold
                    scores[row] = (1 + SpecialFunctions.erf(scores[row]/sqrt(2)))/2
                end
            end
        end
    end
    fetch.(tasks)

end

