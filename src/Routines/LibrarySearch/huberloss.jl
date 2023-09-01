###############
###############

@load "/Users/n.t.wamsley/Desktop/X_test.jld2" X
@load "/Users/n.t.wamsley/Desktop/Hs_test.jld2" Hs
@load "/Users/n.t.wamsley/Desktop/IDtoROW.jld2" IDtoROW
δ = 10000.0
δ = 10000.0
weights = sparseNMF(Hs, X, zero(Float32), zero(Float32));
b = X[:];
X₁ = 1000*ones(eltype(Hs), Hs.n);
X₀ = 1000*ones(eltype(Hs), Hs.n);
r = Hs*X₁ .- b
function solveHuber2!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₀::Vector{T}, X₁::Vector{T}, δ::T; max_iter::Int = 1000, tol::AbstractFloat = 100.0, λ::Float32 = zero(Float32)) where {T<:AbstractFloat}
    L1 = zero(T)
    L2 = zero(T)
    ΔX = Inf
    i = 0
    E = zeros(T, max_iter)
    while (ΔX > tol) & (i < max_iter)
        #ΔX = 0.0
        for col in range(1, Hs.n)
            L1,L2 = zero(T), zero(T)
            n = 0
            X0 = X₁[col]

            ########
            #Newton-Raphson to optimize w.r.t X_col. 
            ########
            while n < 20
                X0 = X₁[col]
                L1 = zero(T)
                L2 = zero(T)
                @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)

                    #Least-Squares 
                    #L1 += Hs.nzval[i]*r[Hs.rowval[i]]
                    #L2 += Hs.nzval[i]^2

                    #Huber
                    #R = 1 + (r[Hs.rowval[i]]/δ)^2
                    #L1 += Hs.nzval[i]*r[Hs.rowval[i]]*(R^(-1/2))
                    #L2 += (Hs.nzval[i]^2)*( R^(-1/2) - (r[Hs.rowval[i]]/δ)^2*(R^(-3/2)) )
                    
                    R = (r[Hs.rowval[i]]/δ)^2
                    L1 += Hs.nzval[i]*r[Hs.rowval[i]]*((1 + R)^(-1/2))
                    L2 += (Hs.nzval[i]^2)/((1 + R)^(3/2))
                end
                #X₁[col] = max(X₁[col] - (L1/L2), 0.0)
                X₁[col] = max(X₁[col] - (L1+λ)/L2, 0.0)
                @turbo for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                    r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X0)
                end

                ########
                #Stopping Criterion for single variable Newton-Raphson
                #Accuracy requirement increases each outer-loop
                ########
                if abs((X₁[col]-X0)/X0) < 0.1/i
                    break
                end
                n += 1
            end
            #println(n)
            #X₁[col] = max(X₁[col], 0.0)
            #ΔX += abs(X₁[col] - X₀[col])
        end


        ###########
        #Stopping Criterion
        ###########
        E[i + 1] = sum(abs.(X₁ .- X₀))
        if E[i + 1] < Hs.n*1
            println("stop at i = $i")
            return E
        end
        
        ###########
        #Reset X₁ and X₀
        ###########
        @turbo for col in eachindex(X₁)
            X₀[col] = X₁[col]
        end

        i += 1
    end
    println(i)
    return E
end
E = solveHuber2!(Hs, r, X₀, X₁, Float32(δ), max_iter = 100, λ = Float32(0.0));
dot(X₁, weights[:])/(norm(X₁)*norm(weights[:]))
corspearman(X₁, weights[:])
dot(X₁.==0.0, weights[:].==0.0)/(norm(X₁.==0.0)*norm(weights[:].==0.0))



plot(log10.(X₁ .+ 1), log10.(weights[:] .+ 1), seriestype=:scatter)

#plot(X₁,weights[:], seriestype=:scatter)


sum(X₁'.==0.0)

describe(X₁[(X₁.!=0.0).&(weights[:].!=0.0)])
describe(X₁[(X₁.!=0.0).&(weights[:].==0.0)])
plot([x for x in range(1, length(E[E.!=0.0]))], E[E.!=0.0])

N = 10
plot([x for x in range(N, length(E[E.!=0.0]))], E[E.!=0.0][N:end])
###############
###############
δ = 10000.0
Hs = sparse(
Float32[1 0; 
 1 0;
 1 1;
 1 0;
 0 1])
b = Float32[100000,
     10000,
     10000,
     0.0,
     0.0]
X₁ = 1000*ones(eltype(Hs), Hs.n);
X₀ = 1000*ones(eltype(Hs), Hs.n);
r = Hs*X₁ .- b


#= Least Squares solution
julia> Matrix(Hs)\b
2-element Vector{Float32}:
  32857.156
 -11428.574
=#

#= NNLS Solution
julia> sparseNMF(Hs, b, zero(Float32), zero(Float32))
1×2 Matrix{Float32}:
30000.0  0.0
=#

#= Huber loss
julia> E = solveHuber2!(Hs, r, X₀, X₁, Float32(δ), max_iter = 100);
stop at i = 2
julia> println(X₁)
Float32[11242.144, 0.0]
=#

###############
###############
δ = 10000.0
Hs = sparse(
Float32[1 0; 
 1 0;
 1 1;
 1 0;
 0 1])
b = Float32[100000,
     10000,
     10000,
     0.0,
     10000]
X₁ = 1000*ones(eltype(Hs), Hs.n);
X₀ = 1000*ones(eltype(Hs), Hs.n);
r = Hs*X₁ .- b


#= Least Squares solution
2-element Vector{Float32}:
 31428.582
 -5714.2886
=#

#= NNLS Solution
julia> sparseNMF(Hs, b, zero(Float32), zero(Float32))
1×2 Matrix{Float32}:
30000.0  0.0
=#

#= Huber loss
julia> E = solveHuber2!(Hs, r, X₀, X₁, Float32(δ), max_iter = 100, λ = Float32(0.0));
stop at i = 3

julia> println(X₁)
Float32[9089.632, 5455.184]
=#

