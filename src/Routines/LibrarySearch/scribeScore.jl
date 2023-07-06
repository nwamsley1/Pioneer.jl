function getDistanceMetrics(H::Matrix{T}, X::Matrix{T}, UNMATCHED::Matrix{T}) where {T<:AbstractFloat}

    function scribeScore(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
       -1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
    end

    function cityBlockDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        log(sum(abs.(a./norm(a) .- b./norm(b))/length(a)))
    end

    function chebyshevDist(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        log(maximum(abs.(a./norm(a) .- b./norm(b))))
        #dot(a, b)/(norm(a)*norm(b))
    end

    function spectralContrast(a::Vector{T}, b::Vector{T}, c::Vector{T}) where {T<:AbstractFloat}
        norma_matched = 0.0
        norma_all = 0.0
        normb = 0.0
        dot = 0.0
        for i in eachindex(a)
            dot += a[i]*b[i]
            norma_matched += a[i]*a[i]
            normb += b[i]*b[i]
        end
        norma_all = norma_matched
        for i in eachindex(c)
            norma_all += c[i]*c[i]
        end
        return dot/(sqrt(norma_matched)*sqrt(normb)), dot/(sqrt(norma_all)*sqrt(normb))
    end

    N = size(H)[1]
    scribe_score = zeros(T, N) #Vector{T}(undef, N)
    city_block_dist = zeros(T, N)
    chebyshev_dist = zeros(T, N)
    matched_ratio = zeros(T, N)
    spectral_contrast_matched = zeros(T, N)
    spectral_contrast_all = zeros(T, N)
    for row in range(1, N)#range(1, N))
        non_zero = (H[row,:].!=0)
        h = H[row,non_zero]
        if isempty(h)
            continue
        end
        x = X[1, non_zero]
        #println("h \n $h")
        unmatched = UNMATCHED[row,:]
        scribe_score[row] = scribeScore(sqrt.(h), sqrt.(x))
        city_block_dist[row] = cityBlockDist(h,x)
        chebyshev_dist[row] = chebyshevDist(h, x)
        matched_ratio[row] = sum(h)/sum(unmatched)
        spectral_contrast_matched[row], spectral_contrast_all[row] = spectralContrast(h, x, unmatched)
    end

    return scribe_score, city_block_dist, chebyshev_dist, matched_ratio, spectral_contrast_matched, spectral_contrast_all
end