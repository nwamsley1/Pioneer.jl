using Combinatorics, StatsBase, Distributions

struct KendallTau{T<:AbstractFloat}
    ecdfs::Vector{ECDF{Vector{T}, Weights{Float64, Float64, Vector{Float64}}}}
end

function KendallTau(dt::DataType)
    kt = KendallTau(Vector{ECDF{Vector{dt}, Weights{Float64, Float64, Vector{Float64}}}}(undef, 10))
    setECDFs!(kt)
    return kt
end

function setECDFs!(kt::KendallTau{T}) where {T<:AbstractFloat}
    function setECDF!(kt::KendallTau{T}, N::Int)
        τs = Vector{T}(undef, factorial(N))
        cannonical_ordering = [n for n in 1:N]
        for (i, perm) in enumerate(permutations(cannonical_ordering))
            τs[i] = corkendall(perm, cannonical_ordering)
        end
        kt.ecdfs[N] = ecdf(τs)
    end

    for N in 3:10
        setECDF!(kt, N)
    end
end 

function getPValue(kt::KendallTau{T}, x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
    N = length(x)
    τ = corkendall(x, y)
    if N <= 10 #Use exact distribution
        if N > 2
            return 1 - (kt.ecdfs[N](τ) - 1/factorial(N))
        else
            return one(T)
        end
    else #Use approximate distribution
        σ = sqrt(2*(2*N + 5)/(9*N*(N - 1)))
        return 1 - cdf(Normal(), τ/σ)
    end
end
