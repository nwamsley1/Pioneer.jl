function getS(peptides::AbstractVector{String}, peptides_dict::Dict{String, Int64}, experiments::AbstractVector{UInt32}, experiments_dict::Dict{UInt32, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where T <: Real
    S = Matrix(zeros(M, N))

    for i in eachindex(peptides)
            S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = coalesce(abundance[i], 0.0)
    end
    return S
end
getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N)
function getB(S::Matrix{Float64}, N::Int)
    B = zeros(N + 1)
    for i in 1:N
        for j in (i+1):N
                r_i_j = median(-log2.(S[:,i]) + log2.(S[:,j]))
                if (S[i, j] != 0.0)
                    B[i] = B[i] - r_i_j# M[i, j]
                    B[j] = B[j] + r_i_j
                end
        end 
    end
    B.*=2
    B[end] = sum(S[S.!=0])*N/length(S) #Normalizing factor
    B
end

function getA(N::Int)
    A = ones(N+1, N+1)
    for i in 1:(N)
        for j in 1:N
            if i == j
                A[i, j] = 2*(N - 1)
            else
                A[i, j] = -1*(N-1)
            end
        end
    end
    A[end, end] = 0
    A
end

function getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Union{T, Missing}},
                          protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where T <: Real

    unique_experiments = unique(experiments)
    unique_peptides = unique(peptides)

    N = length(unique_experiments)
    M = length(unique_peptides)

    peptides_dict = Dict(zip(unique_peptides, 1:M))
    experiments_dict = Dict(zip(unique_experiments, 1:N))

    function appendResults!(protein::String, peptides::Vector{String}, log2_abundances::Vector{Float64}, protein_out::Vector{String}, peptides_out::Vector{String}, log2_abundance_out::Vector{Float64}, S::Matrix{Float64})
        
        function appendPeptides!(peptides_out::Vector{String}, peptides::Vector{String}, S::Matrix{Float64})
            for j in eachindex(eachcol(S))
                peps = String[]
                for i in eachindex(@view(S[:,j]))
                    if S[i,j] > 0.0
                        push!(peps, peptides[i])
                    end
                end
                push!(peptides_out, join(peps, ";"))
            end
        end
        append!(log2_abundance_out, log2_abundances)
        append!(experiments_out, unique_experiments)
        append!(protein_out, [protein for x in 1:N])
        appendPeptides!(peptides_out, peptides, S)
    end
    if M <=1
        appendResults!(protein, unique_peptides, coalesce.(log2.(abundance), 0.0), protein_out, peptides_out, log2_abundance_out, ones(1, N))
        return 
    end

    S = getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N)
    B = getB(S, N)
    A = getA(N)
    log2_abundances = (A\B)[1:(end - 1)]
    appendResults!(protein, unique_peptides, log2_abundances, protein_out, peptides_out, log2_abundance_out, S)

end

function getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Missing},
    protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where T <: Real
    return 
end
