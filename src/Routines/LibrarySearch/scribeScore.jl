function getScribeScore(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function scribeScore(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
       -1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
    end

    N = size(H)[1]
    scribe_score = Vector{T}(undef, N)

    for row in range(1, N)
        non_zero = (H[row,:].!=0) .&(X[1,:].!=0) 
        scribe_score[row] = scribeScore(sqrt.(H[row, non_zero]), sqrt.(X[1, non_zero]))
    end

    return scribe_score
end

#=function getScribeScore(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function scribeScore(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
       dot(a, b)
    end

    N = size(H)[1]
    scribe_score = Vector{T}(undef, N)

    for row in range(1, N)
        non_zero = (H[row,:].!=0) .&(X[1,:].!=0) 
        scribe_score[row] = scribeScore(H[row, non_zero], X[1, non_zero])
    end

    return scribe_score
end=#