using Statistics 

function getS(peptides::AbstractVector{String}, peptides_dict::Dict{String, Int64}, experiments::AbstractVector{UInt32}, experiments_dict::Dict{UInt32, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where T <: Real
    #S = allowmissing(Matrix(missings(M, N)))
    S = Array{Union{Missing,Float64}}(undef, (M, N))
    for i in eachindex(peptides)
            if !ismissing(abundance[i])#isinf(log2(coalesce(abundance[i], 0.0)))
               if abundance[i] != 0.0
                    S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = abundance[i]
               else
                    S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = missing
               end
            end
    end
    return S
end

function getB(S::Matrix{Union{Missing, Float64}}, N::Int, M::Int)
    B = zeros(N + 1)
    for i in 1:N
        for j in 1:N
            if j != i
                r_i_j = skipmissing(-log2.( @view(S[:,j]) ) .+ log2.(@view(S[:,i])))
                if !isempty(r_i_j)
                    B[i] += median(r_i_j)
                end
            end
        end
    end
    B.*=2
    norm = 0.0
    for row in eachrow(S)
        norm += sum(log2.(skipmissing(row)))*(length(row)/(length(row) - sum(ismissing.(row))))
    end
    B[end] = norm/M 
    B
end

function getA(N::Int)
    A = ones(N+1, N+1)
    for i in 1:(N)
        for j in 1:N
            if i == j
                A[i, j] = 2*(N - 1)
            else
                A[i, j] = -2
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

    function appendResults!(protein::String, peptides::Vector{String}, log2_abundances::Vector{Float64}, protein_out::Vector{String}, peptides_out::Vector{String}, log2_abundance_out::Vector{Float64}, S::Matrix{Union{Missing, Float64}})
        
        function appendPeptides!(peptides_out::Vector{String}, peptides::Vector{String}, S::Matrix{Union{Missing, Float64}})
            for j in eachindex(eachcol(S))
                peps = String[]
                for i in eachindex(@view(S[:,j]))
                    if !ismissing(S[i,j])
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
    #=if M <=1
        appendResults!(protein, unique_peptides, coalesce.(log2.(abundance) .- mean(log2.(abundance)), 0.0), protein_out, peptides_out, log2_abundance_out, allowmissing(ones(1, N)))
        return 
    end=#

    S = allowmissing(getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N))

    B = getB(S, N, M)
    A = getA(N)
    log2_abundances = (A\B)[1:(end - 1)]
    appendResults!(protein, unique_peptides, log2_abundances, protein_out, peptides_out, log2_abundance_out, S)

end

#=function getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Missing},
    protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where T <: Real
    return 
end=#
