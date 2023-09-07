###############
###############
using JLD2, LinearAlgebra, Dictionaries, SparseArrays, LoopVectorization, StatsBase
@load "/Users/n.t.wamsley/Desktop/X_test.jld2" X
@load "/Users/n.t.wamsley/Desktop/Hs_test.jld2" Hs
@load "/Users/n.t.wamsley/Desktop/IDtoROW.jld2" IDtoROW

include("src/Routines/LibrarySearch/searchRAW.jl")
include("src/ML/sparseNNLS.jl")
Hs, weights, X, IDtoROW = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    #scan_range = (0, length(MS_TABLE[:scanNumber]))
    scan_range = (40000, 70000)
#scan_range = (101357, 102357)
);


sum(X[Hs[:,1].!=0.0])./sum(Hs[:,1][Hs[:,1].!=0.0])
sum(X[Hs[:,8].!=0.0])./sum(Hs[:,8][Hs[:,8].!=0.0])
X[Hs[:,2].!=0.0]
include("src/ML/sparseNNLS.jl")
#weights[:] = weights0
#weights0 = weights[:]
#weights = Float32.(max.(Hs\X, Float32(0.0)))
weights = sparseNMF(Hs, X, zero(Float32), zero(Float32), false, max_iter=100, tol=Float32(1000.0))[:]
weights0 = weights[:]
i_ = @time solveHuber!(Hs, Hs*weights .- X,  weights, Float32(10000), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n*100);
weights0[1:2]
weights[1:2]
corspearman(weights0, weights)
plot(log.(weights0), log.(weights), seriestype=:scatter)

include("src/ML/sparseNNLS.jl")
δ = Float32(20000.0)
#δ = 1000000.0
weights = sparseNMF(Hs, X, zero(Float32), zero(Float32));
b = X[:];
@profview begin
for i in range(1, 100)
    X₁ = 0*ones(eltype(Hs), Hs.n);
    r = Hs*X₁ .- b
    solveHuber4!(Hs, r, X₁, Float32(δ), max_iter_outer = 100, max_iter_inner = 10, tol = Hs.n*1000, λ = Float32(0.0));
end
end

println("reps $reps")
X₁ = X₁*Float32(1.1);
r = Hs*X₁ .- b;
reps = @time solveHuber4!(Hs, r, X₁, Float32(δ), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n*1000, λ = Float32(0.0));
println("reps $reps")

X₁ = 0*ones(eltype(Hs), Hs.n);
r = Hs*X₁ .- b;
reps = @time solveHuber!(Hs, r, X₁, Float32(δ), max_iter_outer = 100, max_iter_inner = 10, tol = Hs.n*1000);
println("reps $reps")
dot(X₁, weights[:])/(norm(X₁)*norm(weights[:]))
corspearman(X₁, weights[:])
dot(X₁.==0.0, weights[:].==0.0)/(norm(X₁.==0.0)*norm(weights[:].==0.0))
plot(log10.(X₁ .+ 1), log10.(weights[:] .+ 1), seriestype=:scatter, alpha = 0.2)

X₁ = X₁*Float32(1.1);
r = Hs*X₁ .- b;
@time solveHuber!(Hs, r, X₁, Float32(δ), max_iter_outer = 100, max_iter_inner = 10, tol = Hs.n*1000);
dot(X₁, weights[:])/(norm(X₁)*norm(weights[:]))
corspearman(X₁, weights[:])
dot(X₁.==0.0, weights[:].==0.0)/(norm(X₁.==0.0)*norm(weights[:].==0.0))
#Hs_test, X₁_test = @time solveHuber4!(Hs, r, X₁, Float32(δ), max_iter_outer = 100, max_iter_inner = 20, tol = Hs.n*1000, λ = Float32(0.0));

test_vals = LinRange(-1e7, 3e7, 1000)
r = Hs*weights .- X;
X₁_test = weights[:]
Hs_test = Hs
costs = Float64[]
L1s = Float64[]
L2s = Float64[]
δ = Float32(5000)
COL = 1
for val in test_vals
    X₁_test[COL] = val
    r = Hs*X₁_test .- X
    cost = 0.0
    L1 = 0.0
    L2 = 0.0
    for i in Hs_test.colptr[COL]:(Hs_test.colptr[COL + 1] - 1)
        cost += (δ^2)*(sqrt(1 + (r[Hs_test.rowval[i]]/δ)^2) - 1)
        R = Hs_test.nzval[i]*(1 + (r[Hs_test.rowval[i]]/δ)^2)^(-1/2)
        L1 += r[Hs.rowval[i]]*((R))
        L2 += ((R)^(3))/Hs.nzval[i]
    end
    push!(costs, cost)
    push!(L1s, L1)
    push!(L2s, L2)
end

test_vals[end] - L1s[end]/L2s[end]

plot(test_vals, costs)
plot!(test_vals, costs[end] .+ L1s[end].*(test_vals .- test_vals[end]) .+ (L2s[end]/2).*(test_vals .- test_vals[end]).^2)
vline!([test_vals[end] - costs[end]/L1s[end]])
vline!([weights[1]])

#plot(test_vals, costs)
#plot!(test_vals, costs[end] .+ L1s[end].*(test_vals .- test_vals[end]) .+ (L2s[end]*150/2).*(test_vals .- test_vals[end]).^2)
a = (L2s[end]/2)
b =  L1s[end]
c = costs[end]
vline!([(-b - sqrt(b^2 - 4*a*c))/(2*a)])

#vline!([test_vals[end] - L1s[end]/L2s[end]])
vline!([test_vals[end] - costs[end]/L1s[end]])
plot!(test_vals, L1s[end].*test_vals)
plot(test_vals, log2.(costs))
plot(test_vals, L1s)
hline!([9000])
plot(test_vals, log.(abs.(L1s)))
plot(test_vals, L2s)


histogram(ns)





plot(log10.(X₁ .+ 1), log10.(weights[:] .+ 1), seriestype=:scatter, alpha = 0.2)

function rsqrt(x::Float32)
    xₛ = x
    int32 = reinterpret(UInt32, xₛ)
    int32 = 0x5f3759df - int32 >> 1
    xₛ = reinterpret(Float32, int32)
    xₛ *= 1.5f0 - x * 0.5f0 * xₛ^2
    return xₛ
end
function rsqrt(x::Float64)
    xₛ = x
    int64 = reinterpret(UInt64, xₛ)
    int64 = 0x5fe6eb50c7b537a9 - int64 >> 1  # See https://stackoverflow.com/a/11644533/3260253
    xₛ = reinterpret(Float64, int64)
    xₛ *= 1.5 - x * 0.5 * xₛ^2
    return xₛ
end
rsqrt(x::Real) = rsqrt(float(x))

test_nums = 10000*rand(100000)

@btime (test_nums).^(-1/2)

@btime rsqrt.(test_nums)

histogram(ns)


X₁ = X₁ = 100*ones(eltype(Hs), Hs.n);
#X₁ = [Float32(x) for x in weights[:]]
r = Hs*X₁ .- b;
#=@time grad_, diffs_, ns_ = solveHuber3!(Hs, r, X₁, Float32(δ), max_iter_outer = 5000, max_iter_inner = 20, tol = Hs.n*10, λ = Float32(0.0));
plot([i for i in range(1, length(grad_[grad_.<Inf]))], grad_[grad_.<Inf])
plot([i for i in range(1, length(diffs_))], log.(diffs_))
plot!([i for i in range(1, length(cumsum(diffs_)))], 
log2.(cumsum(diffs_)./([i for i in range(1, length(cumsum(diffs_)))])))
plot!([i for i in range(1, length(moving_average(diffs_, 100)))], 
log2.(moving_average(diffs_, 100))
)=#
@time solveHuber3!(Hs, r, X₁, Float32(δ), max_iter_outer = 1000, max_iter_inner = 20, tol = Hs.n*10, λ = Float32(0.0));
@profview solveHuber3!(Hs, r, X₁, Float32(δ), max_iter_outer = 1000, max_iter_inner = 20, tol = Hs.n*10, λ = Float32(0.0));
#moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]


dot(X₁, weights[:])/(norm(X₁)*norm(weights[:]))
corspearman(X₁, weights[:])
dot(X₁.==0.0, weights[:].==0.0)/(norm(X₁.==0.0)*norm(weights[:].==0.0))
plot(log10.(X₁ .+ 1), log10.(weights[:] .+ 1), seriestype=:scatter, alpha = 0.2)

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

