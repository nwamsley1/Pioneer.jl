using GLM
#Design matrix from an actual scan 
Hs = load("data/Hs.jld2")["Hs"]
N = Hs.n_vals

#Get Data
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end
#Design Matrix 
H = Matrix(sparse(Hs.rowval[1:N],
Hs.colval[1:N],
Hs.nzval[1:N]))
#OLS Regression 
rowvals = copy(Hs.rowval)
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end

_residuals_ = zeros(Float32, Hs.m)
_weights_ = zeros(Float32, Hs.n)
initResiduals!(_residuals_, Hs, _weights_)
solveHuber!(Hs, _residuals_, _weights_, 
            Float32(1e5), #δ
            Float32(0.0), #λ
            100,
            100,
            100, 
            10.0f0,
            10.0f0,
            10.0f0,#Hs.n/10.0,
            0.01f0
            );

Ht =  Matrix(transpose(H))
z = reshape(y, (1, Hs.m))
ŷ = zeros(Float32, (1, Hs.n))
test_alg = NMF.CoordinateDescent{Float32}(
                maxiter = 1000,
                tol = 1e-6,
                α=0.0, 
                l₁ratio=1.0,
                update_H = false,
                regularization=:none)
NMF.solve!(test_alg, z, ŷ, Matrix(transpose(H)))

sort(abs.(ŷ[:,] .- _weights_)./ŷ[:,1])

@test cor(ŷ[:,], _weights_) > 0.9999

plot(log2.(ŷ[:,]), log2.(_weights_), seriestype=:scatter)
Xs = sparse(X')

