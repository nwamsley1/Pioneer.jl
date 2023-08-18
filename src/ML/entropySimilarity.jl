function unweightedEntropy(X::Vector{Union{Missing, T}}) where {T<:AbstractFloat}
    p = X./sum(X)
    p, StatsBase.entropy(p)
end

function weightedEntropy(X::Vector{Union{Missing, T}}) where {T<:AbstractFloat}
    p = X./sum(X)
    S = StatsBase.entropy(p)
    if S >=3
        return p, S
    else
        w = 0.25*(1 + S)
        intensity = p.^w
        return intensity, StatsBase.entropy(intensity./sum(intensity))
    end
end