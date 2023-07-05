function getScribeScore(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function scribeScore(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
       -1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
    end

    function cityBlockDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        log(sum(abs.(a./norm(a) .- b./norm(b))/length(a)))
    end

    function chebyshevDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        #log(max(abs.(norm(a) .- norm(b))))
        dot(a, b)/(norm(a)*norm(b))
    end

    N = size(H)[1]
    scribe_score = Vector{T}(undef, N)
    city_block_dist = Vector{T}(undef, N)
    chebyshev_dist = Vector{T}(undef, N)
    unmatched = Vector{T}(undef, N)

    for row in range(1, N)
        non_zero = (H[row,:].!=0) .&(X[1,:].!=0) 
        scribe_score[row] = scribeScore(sqrt.(H[row, non_zero]), sqrt.(X[1, non_zero]))
        city_block_dist[row] = cityBlockDist(H[row, non_zero], X[1, non_zero])
        chebyshev_dist[row] = chebyshevDist(H[row, non_zero], X[1, non_zero])
        unmatched[row] = sum(H[row,(X[1,:].!=0)])/sum(H[row,(X[1,:].==0)])

    end

    return scribe_score, city_block_dist, chebyshev_dist, unmatched
end

function getScribeScore(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function scribeScore(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
       -1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
    end

    function cityBlockDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        log(sum(abs.(a./norm(a) .- b./norm(b))/length(a)))
    end

    function chebyshevDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        #log(max(abs.(norm(a) .- norm(b))))
        dot(a, b)/(norm(a)*norm(b))
    end

    N = size(H)[1]
    scribe_score = Vector{T}(undef, N)
    city_block_dist = Vector{T}(undef, N)
    chebyshev_dist = Vector{T}(undef, N)
    unmatched_ratio = Vector{T}(undef, N)

    for row in range(1, N)
        matched = (H[row,:].!=0) .&(X[1,:].!=0) 
        unmatched = (H[row,:].!=0) .& (X[1,:].==0) 
        scribe_score[row] = scribeScore(sqrt.(H[row, matched]), sqrt.(X[1, matched]))
        city_block_dist[row] = cityBlockDist(H[row, matched], X[1, matched])
        chebyshev_dist[row] = chebyshevDist(H[row, matched], X[1, matched])
        unmatched_ratio[row] = sum(H[row,unmatched])/sum(H[row,unmatched])

    end

    return scribe_score, city_block_dist, chebyshev_dist, unmatched_ratio
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