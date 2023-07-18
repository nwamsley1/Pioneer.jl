target_frag_list, target_frag_detailed, target_precursor_list, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_targets.csv",
                                                                                        precision = Float32,
                                                                                        isDecoys = false)

decoy_frag_list, decoy_frag_detailed, decoy_precursor_list, prec_id = readPrositLib("/Users/n.t.wamsley/Projects/PROSIT/myPrositLib_decoys.csv",
                                                                                        precision = Float32,
                                                                                        isDecoys = false,
                                                                                        first_prec_id = prec_id)
#Merge target and decoy lists
frag_list = append!(target_frag_list, decoy_frag_list)
target_frag_list = nothing
decoy_frag_list = nothing
frag_detailed = append!(target_frag_detailed, decoy_frag_detailed)
target_frag_detailed = nothing
decoy_frag_detailed = nothing
precursor_list = append!(target_precursor_list, decoy_precursor_list)
target_precursor_list = nothing
decoy_precursor_list = nothing
#Normalize
function NormalizeIntensities!(frag_list::Vector{FragmentIon{T}}, prec_list::Vector{LibraryPrecursor}) where {T<:AbstractFloat}
    for i in 1:length(frag_list)
        prec_id = frag_list[i].prec_id
        frag_list[i].prec_intensity[] = frag_list[i].prec_intensity[]/prec_list[prec_id].total_intensity[]
    end
end
#Save to Disc
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_list.jld2" frag_list
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2" frag_detailed
@save "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list

#Build Fragment Index 
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_list.jld2" frag_list
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/frag_detailed.jld2" frag_detailed
@load "/Users/n.t.wamsley/Projects/PROSIT/CombinedNormalized_071823/precursor_list.jld2" precursor_list
for i in 1:length(frag_list)

end


#h = ecdf(rand(100))
#h()
using Combinatorics, StatsBase, Distributions
struct KendallTau{T<:AbstractFloat}
    ecdfs::Vector{ECDF{Vector{T}, Weights{Float64, Float64, Vector{Float64}}}}
end
KendallTau(dt::DataType) = KendallTau(Vector{ECDF{Vector{dt}, Weights{Float64, Float64, Vector{Float64}}}}(undef, 10))

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
        return 1 - kt.ecdfs[N](τ)
    else #Use approximate distribution
        σ = sqrt(2*(2*N + 5)/(9*N*(N - 1)))
        return 1 - cdf(Normal(), τ/σ)
    end
end

KendallTauECDFs = Vector{ECDF{Vector{Float32}, Weights{Float64, Float64, Vector{Float64}}}}(undef, 10)
A = [x for x in 1:9]
τs = Vector{Float32}(undef, factorial(9))

for (i, perm) in enumerate(permutations(A))
    τs[i] = corkendall(perm, A)
    if i == 10000
        println(perm)
    end
end

kk9 = ecdf(τs)