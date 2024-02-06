

function model_predict(psms::DataFrame, model_fit::Any, column_names::Vector{Symbol})
    coefs = coef(model_fit)
    psms[!,:prob] = zeros(Float16, size(psms, 1))
    function addColumn(probs::Vector{Float16}, vals::Vector{R}, x::Float64) where {R<:Real}
        Threads.@threads for i in range(1, size(psms, 1))
            probs[i] += vals[i]*x
        end
    end
    for (i, name) in enumerate(column_names)
        x = coefs[i+1]
        addColumn(psms[!,:prob], psms[!,name], x)
    end
    intercept = coefs[1]::Float64
    Threads.@threads for i in range(1, size(psms, 1))
        psms[i,:prob] += Float16(intercept)
    end

    Threads.@threads for i in range(1, size(psms, 1))
        psms[i,:prob] = exp( psms[i,:prob] )/(1 + exp( psms[i,:prob] ))
    end

end

randomforest_surrogate = RandomForestSurrogate(Matrix(PSMs[!,column_names]),PSMs[!,:target],[0.0, 2.0num_round = 2)
LDA = MLJ.@load LDA pkg=MultivariateStats
?LDA
model = LDA()
mach = machine(model, X, y)

function fit()

X = Float16.(X)
X = hcat(X, ones(Float16, size(X, 1)))

X = Float64.(X)
y_hat = [0.5 for _ in range(1, length(y))]
y = Float64.(PSMs[!,:target])
Z = zeros(Float64, size(X, 1))
θ = X\y
∇θ = zeros(Float64, size(X, 2))
m = size(X, 1)
n = size(X, 2)

function costFunc(z::Vector{T}, y::Vector{Bool}) where {T<:AbstractFloat}
    cost = zero(T)
    for i in eachindex(z)
        if y[i]
            cost += log(logistic(z[i]) + 1e-6)/length(z)
        else
            cost += log(1 - logistic(z[i]) + 1e-6)/length(z)
        end
    end
    return cost
end

function update∇θ!(∇θ::Vector{T}, X::Matrix{Float64}, z::Vector{T}, y::Vector{Bool}) where {T<:AbstractFloat}
    fill!(∇θ, zero(T))
    m = size(X, 1)
    Threads.@threads for col in range(1, size(X, 2))
        for row in range(1, size(X, 1))
            #if y[row]
            #    ∇θ[col] -= X[row, col]/(1 + exp(z[row]))
            #else
            #    ∇θ[col] += X[row, col]*exp(z[row])/(1 + exp(z[row]))
            #end
            y_hat = logistic(z[row])#1/(1 + exp(-z[row]))
            if y[row]
                ∇θ[col] += X[row, col]*(one(T) - y_hat)/m
            else
                ∇θ[col] -= X[row, col]*y_hat/m
            end
        end
    end
end

function logistic(x::T) where {T<:AbstractFloat}
    return 1/(1 + exp(-x))
end

function updateθ!(θ::Vector{Float64}, ∇θ::Vector{Float64}, X::Matrix{Float64}, y::Vector{Bool}; α = 0.01, max_iter::Int64 = 10) #where {T<:AbstractFloat}
    costs = zeros(max_iter)
    m = size(X, 1)
    z = zeros(m)
    for i in range(1, max_iter)

        mul!(z, X, θ)
        costs[i] = costFunc(z, y)
        println("costs[i] ", costs[i])
        update∇θ!(∇θ, X, z, y)
        ∇θ = ∇θ#LinearAlgebra.norm2(∇θ)
        #println("∇θ $∇θ")
        for n in range(1, length(θ))
            θ[n] -= α*(∇θ[n])
        end
    end
    mul!(z, X, θ)
    return logistic.(z)
end


N = Distributions.Normal()

β = zeros(size(X, 2))

XWX = zeros((size(X, 2), size(X, 2)))
W = diag(W)
function fillXWX!(XWX::Matrix{T}, X::Matrix{T}, W::Vector{T}) where {T<:AbstractFloat}
    Threads.@threads for i in range(1, size(XWX, 1))
        @inbounds @fastmath for j in range(1, size(XWX, 1))
            for n in range(1, size(X, 1))
                XWX[i, j] += X[n,i]*W[i]*X[n,j]
            end
        end
    end
end

XWX = zeros((size(X, 2), size(X, 2)))
getXWX!(XWX, X, W)

η = X*β
W = zeros(length(η))
Z = zeros(length(η))
  
function fillZandW!(Z::Vector{T}, W::Vector{T}, η::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
        @inbounds @simd for i in range(1, length(y))
            ϕ = Distributions.pdf(Normal(), η[i])
            μ = Distributions.cdf(Normal(), η[i])
            Z[i] = η[i] + (y[i] - μ)/ϕ
            W[i] = (ϕ^2)/(μ*(1 - μ))
        end
end

function fillXWX!(XWX::Matrix{T}, X::Matrix{T}, W::Vector{T}) where {T<:AbstractFloat}
    for i in range(1, size(XWX, 1))
        Threads.@spawn for j in range(1, size(XWX, 1))
            for n in range(1, size(X, 1))
                XWX[i, j] += X[n,i]*W[i]*X[n,j]
            end
        end
    end
end

function fillXWX!(XWX::Matrix{T}, X::Matrix{T}, W::Vector{T}) where {T<:AbstractFloat}
    for i in range(1, size(XWX, 1))
        for j in range(1, size(XWX, 1))
            for n in range(1, size(X, 1))
                XWX[i, j] += X[n,i]*W[i]*X[n,j]
            end
        end
    end
end

function fillY!(Y::Vector{T}, X::Matrix{T}, W::Vector{T}, Z::Vector{T}) where {T<:AbstractFloat}
    col = 0
    Threads.@spawn for col in range(1, size(X, 2))
        Y[col] = zero(T)
        col += 1
        for row in range(1, size(X, 1))
            Y[col] += X[row, col]*W[row]*Z[row]
        end
    end
    return col
end

function fillη!(η::Vector{T}, X::Matrix{T}, β::Vector{T}) where {T<:AbstractFloat}
        @turbo for row in range(1, size(X, 1))
            η[row] = zero(T)
        end
        @turbo for col in range(1, size(X, 2))
            for row in range(1, size(X, 1))
                    η[row] += X[row, col]*β[col]
                end
            end
        end
end

function ProbitRegression(β::Vector{T}, X::Matrix{T}, y::Vector{T}; max_iter::Int = 3) where {T<:AbstractFloat}

    W = zeros(T, length(y))
    η = zeros(T, length(y))
    Z = zeros(T, length(y))
    Y = zeros(T, size(X, 2))
    XWX = zeros(T, (size(X, 2), size(X, 2)))
    old_β = copy(β)
    @time for i in ProgressBar(range(1, max_iter))
        fillη!(η, X, β)
        fillZandW!(Z, W, η, y)
        fillXWXMulti!(XWX, X, W)
        fillY!(Y, X, W, Z)
        β = T.(XWX\Y)
        #println(LinearAlgebra.norm2(β .- old_β)/(LinearAlgebra.norm2(β)))
        old_β = copy(β)
    end
    return β
end


X = hcat(X[:,1:11], X[:,1:11].^2)
X = hcat(X, ones(eltype(X), size(X, 1)))
y = Bool.(y)
Threads.@threads for i in range(1, 11)
    X[:,i] = X[:,i]/std(X[:,i])
end
include("src/Routines/LibrarySearch/test.jl")
β = zeros(Float64, size(X, 2));
@time β = ProbitRegression(β, X, y, max_iter = 15);
PSMs[!,:prob] = Float16.(X*β);
PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

best = ((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)) .| PSMs.decoy;
β = ProbitRegression(β, X[best,:], y[best], max_iter = 15);
PSMs[!,:prob] = Float16.(X*β);
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))


jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs.jld2"); PSMs)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "X.jld2"); X)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "y.jld2"); y)
joinpath(MS_DATA_DIR, "Search", "RESULTS", "y.jld2")


exp(-((-7.0)^2)/2)/sqrt(2*π)
(1 + SpecialFunctions.erf((-7.0)/sqrt(2)))/2


exp(-((7.0)^2)/2)/sqrt(2*π)
(1 + SpecialFunctions.erf((7.0)/sqrt(2)))/2

y = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\y.jld2")["y"]
X = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\X.jld2")["X"]
PSMs = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\RAW\\TEST_y4b3_nOf5\\Search\\RESULTS\\PSMs.jld2")["PSMs"]
PSMs[!,:prob] = Float16.(X*β)
PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

plot( LinRange(-5, 5, 1000),
[logcdf(Normal(), x) for x in LinRange(-5, 5, 1000)])
XWX = zeros(T, (size(X, 2), size(X, 2)))
@time fillXWX!(XWX, X, W)
fillY!(Y, X, W, Z)
β = T.(XWX\Y)

XWX = zeros(T, (size(X, 2), size(X, 2)))
@time fillXWXMulti!(XWX, X, W)
fillY!(Y, X, W, Z)
β = T.(XWX\Y)



best = ((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)) .| PSMs.decoy;
β = ProbitRegression(β, X[best,:], y[best], max_iter = 3);
PSMs[!,:prob] = Float16.(X*β);
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

best = ((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)) .| PSMs.decoy;
β = ProbitRegression(X[best,:], y[best], max_iter = 20);
PSMs[!,:prob] = Float16.(X*θ);
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))


μ = Distributions.cdf(N, η)
z = η .+ (y .- μ)./Distributions.pdf(N, η)
w = (Distributions.pdf(N, η).^2)./(μ.*(1 .- μ))
U = UniformScaling(1.0)
W = Diagonal(U, size(X, 1))
for i in range(1, length(w))
    W[i, i] = w[i]
end
#Z = ones(size(X, 1))
WX = W*X
COV = transpose(X)*WX
Y = (transpose(X)*W)*z
β = COV\Y


θ = zeros(size(X, 2))
∇θ = copy(θ)
y_hat = updateθ!(θ, ∇θ, X, PSMs[!,:target], α = 0.01, max_iter = 15)


θ_test = X\y


X\y



PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

best = ((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)) .| PSMs.decoy;
model_fit = glm(FORM, PSMs[best,:], 
                       Binomial(), 
                       ProbitLink(),
                       verbose = false)
column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]
model_predict(PSMs, model_fit, column_names)
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))


X[:,[1, 2, 3, 4, end]]

column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,
:Mox,:charge,:TIC,:total_ions,:err_norm]
X = Matrix(PSMs[!,column_names])
X = Float64.(X)
X = hcat(X, ones(Float64, size(X, 1)))
X = Float64.(X)
#[0.5 for _ in range(1, length(y))]
y = Float64.(PSMs[!,:target])
#θ = X\y
#y_hat = X*θ
θ = 0.01*ones(size(X, 2))
y_hat = logistic.(X*θ)
∇θ = zeros(length(θ))
θ = updateθ!(θ, ∇θ, X, y_hat, y)


θ = updateθ!(θ, X, ŷ, y)


X = Matrix(PSMs[!,column_names])
model = build_forest(y, X, 2, 10, 0.5, 6)

X = Matrix(PSMs[!,column_names])

@time model = build_forest(PSMs[!,:target],X , 4, 1000, 0.2, 3)

@time PSMs[!,:prob] = Float16.(apply_forest_proba(model, X, [true, false])[:,1])
PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));



=

getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

X = Float64.(Matrix(PSMs[!,column_names]))
lda = MultivariateStats.fit(MulticlassLDA, X, y)


lda = MultivariateStats.fit(LinearDiscriminant, X[PSMs[!,:target],:], X[PSMs[!,:decoy],:])
PSMs[!,:target]
reshape(PSMs[!,:target], size(X, 1), 1)
lda = MultivariateStats.fit(MulticlassLDA, X,  PSMs[!,:target]; outdim = 1)

mean(X, 1)
lda = MultivariateStats.fit(LinearDiscriminant, X, PSMs[!,:target])

μ = zeros((2, size(X, 2)))
Σ = zeros((1, size(X, 2)))

function LDA(X::Matrix{Float16}, targets::Vector{Bool}, decoys::Vector{Bool}, pview::AbstractArray{Int64})
    μ = zeros(Float64, (3, size(X, 2)))
    Σ = zeros(Float64, (1, size(X, 2)))
    πk = zeros(Float64, 2)
    Nk = zeros(Int64, 2)
    δk = zeros(Float64, (2, size(X, 1)))

    @time Threads.@threads for col in range(1, size(X, 2))
        @inbounds @fastmath for row in pview#range(1, size(X, 1))
            if targets[row]
                Nk[1] += 1
                μ[1, col] += X[row, col]
            elseif decoys[row]
                Nk[2] += 1
                μ[2, col] += X[row, col]
            end
        end
    end
    μ[1,:] = μ[1,:]/Nk[1]
    μ[2,:] = μ[2,:]/Nk[2]
    μ[3,:] = (μ[1,:] .+ μ[2,:])./2
    @time Threads.@threads for col in range(1, size(X, 2))
        @inbounds @fastmath for row in pview#range(1, size(X, 1))
            if (targets[row] | decoys[row])
                Σ[1, col] += (X[row, col] - μ[3, col])^2
            end
        end
    end
    Σ = Σ./(sum(Nk) + 1)
    Σ_inv = zeros((size(X, 2), size(X, 2)))
    for i in range(1, size(X, 2))
        Σ_inv[i, i] = 1/Σ[1, i]
    end
    πk[1] = Nk[1]/sum(Nk)
    πk[2] = Nk[2]/sum(Nk)
    return Float16.(log(πk[1]/πk[2]) - transpose(μ[1,:] .+ μ[2,:])*Σ_inv*(μ[1,:] .- μ[2,:])/2 .+ X*Σ_inv*diff_μ)
end
PSMs[!,:prob] = LDA(X, PSMs[!,:target], PSMs[!,:decoy], @view(p[1:end]))
PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
targets = [(true==x) for x in ((PSMs[!,:q_value].<=0.01).&(PSMs[!,:target]))];
DataFrames.shuffle!(p)
PSMs[!,:prob] = LDA(X, targets, PSMs[!,:decoy],@view(p[1:end]));
PSMs[!,:q_value] = zeros(Float16, size(PSMs, 1));
getQvalues!(PSMs[!,:prob],  PSMs[!,:decoy],PSMs[!,:q_value]);
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))

logoddsscore = log(πk[1]/πk[2]) - transpose(μ[1,:] .+ μ[2,:])*Σ_inv*(μ[1,:] .- μ[2,:])/2 .+ X*Σ_inv*diff_μ


diff_μ = μ[1,:] .- μ[2,:]
X*(Σ_inv*((μ[1,:] .- μ[2,:])))


μ = μ[1:2,:]


transpose(Σ_inv*transpose(μ))*transpose(μ) #.+ log.(πk)

μ

(μ[1,:] .- μ[2,:])
Σ_inv = inv(cov(Float32.(X)))
C1 = transpose(μ[1,:])*Σ_inv*μ[1,:] + log(πk[1])
C2 = transpose(μ[2,:])*Σ_inv*μ[2,:] + log(πk[2])
p1 = X*(Σ_inv*μ[1,:]) .- C1
p2 = X*(Σ_inv*μ[2,:]) .- C2
tranpose(μ)*transpose(Σ_inv)*μ
X = Float16.(X)

PSMs[!,:prob] = p1


μ[1,:] = μ[1,:]/sum(PSMs[!,:target])
μ[2,:] = μ[2,:]/sum(PSMs[!,:decoy])
Threads.@threads for col in range(1, size(X, 2))
    for row in range(1, size(X, 1))
        if PSMs[row,:target2] | PSMs[row,:decoy]
            Σ[1, col] += (X[row, col] - μ[1, col])^2
        end
        #else
        #    Σ[2, col] += (X[row, col] - μ[2, col])^2
        #end
    end
end
Σ[1,:] = Σ[1,:]/sum(PSMs[!,:target2])

C1 = MvNormal(μ[1,:], Σ[1,:])
C2 = MvNormal(μ[2,:], Σ[1,:])


Distributions.logpdf(C1,X[1,:])
Distributions.logpdf(C2,X[1,:])

probs = zeros(Float64, size(X, 1))
Threads.@threads for i in ProgressBar(range(1, size(X, 1)))
    probs[i] = Distributions.logpdf(C1,X[i,:]) - Distributions.logpdf(C2,X[i,:])
end

PSMs[!,:target2] = (PSMs[!,:q_value] .<= 0.1).&(PSMs[!,:target])
PSMs[!,:decoy2] = (PSMs[!,:decoy])

PSMs[!,:prob] = Float16.(probs)
using ScikitLearn

using DecisionTree

build_forest()
#=

X = Matrix(PSMs[!,:])
logistic = LogisticRegression(0.0, fit_intercept = false)

solver = MLJLinearModels.NewtonCG(
    optim_options = Optim.Options(time_limit = 20),
    newtoncg_options = (eta = 0.2,)
)
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)

solver = MLJLinearModels.LBFGS(
    optim_options = Optim.Options(time_limit = 20),
    lbfgs_options = (linesearch = Optim.LineSearches.HagerZhang(),))


@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)

PSMs[!,:intercept] = ones(Float16, size(PSMs, 1))
X = Matrix(PSMs[!,[:intercept,:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])
y = PSMs[!,:target]

logistic = LogisticRegression(fit_intercept = false)
solver = MLJLinearModels.NewtonCG(
    optim_options = Optim.Options(time_limit = 20),
    newtoncg_options = (eta =0.0,)
)
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)
PSMs[!,:prob] = X*theta
getQvalues!(PSMs, allowmissing(PSMs[!,:prob]),  allowmissing(PSMs[:,:decoy]));
println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
  
@time theta = MLJLinearModels.fit(logistic, X, y, solver = solver)



PSMs[!,:prob] = X*theta[2:end]



function generate_dataset(n_samples = 100, n_features = 10; shift = 0.0)
    X = randn(n_samples, n_features)
    w = randn(n_features)
    y = sign.(X * w)
    X .+= 0.8 * randn(n_samples, n_features) # add noise
    X .+= shift # shift the points in the feature space
    X = hcat(X, ones(n_samples, 1))
    return X, y
end

function softplus(model, t, u)
    z = @variable(model, [1:2], lower_bound = 0.0)
    @constraint(model, sum(z) <= 1.0)
    @constraint(model, [u - t, 1, z[1]] in MOI.ExponentialCone())
    @constraint(model, [-t, 1, z[2]] in MOI.ExponentialCone())
end

function build_logit_model(X, y, λ)
    n, p = size(X)
    model = Model()
    @variable(model, θ[1:p])
    @variable(model, t[1:n])
    for i in 1:n
        u = -(X[i, :]' * θ) * y[i]
        softplus(model, t[i], u)
    end
    # Add ℓ2 regularization
    @variable(model, 0.0 <= reg)
    @constraint(model, [reg; θ] in SecondOrderCone())
    # Define objective
    @objective(model, Min, sum(t) + λ * reg)
    return model
end

n, p = 200000, 10
X, y = generate_dataset(n, p; shift = 10.0);

λ = 10.0
model = build_logit_model(X, y, λ)
set_optimizer(model, SCS.Optimizer)
set_silent(model)
@time JuMP.optimize!(model)

θ♯ = JuMP.value.(model[:θ])


#=
y = Float16.(PSMs[!,:target])
PSMs[!,[:spectral_contrast,:scribe,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]]

@time X = Matrix(PSMs[!,[:spectral_contrast,:scribe,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])

theta_guess = ones(Float64, size(X, 2))

@time theta, flag, history = gradient_descent_probit(y, X, theta_guess, max_iter=1000, tol=0.00001);

hessian_log_likelihood(y, X, theta_guess)
=#
logistic = LogisticRegression(0.5)
theta = fit(logistic, X, y)



PSMs[!,:intercept] = ones(Float16, size(PSMs, 1))
X = Matrix(PSMs[!,[:intercept,:spectral_contrast,:scribe,:city_block,:entropy_score,:iRT_error,:missed_cleavage,:Mox,:charge,:TIC,:total_ions,:err_norm,:spectrum_peak_count]])
Y = PSMs[!,:target]
n, p = size(X)
beta = Variable(p)
problem = minimize(logisticloss(-Y .* (X * beta)))
solve!(problem, SCS.Optimizer; silent_solver = true)

=#