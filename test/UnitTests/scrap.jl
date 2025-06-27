# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


#Design matrix from an actual scan 
#=
using Base.Order
using LoopVectorization
using SparseArrays
using JLD2
using NMF
using StatsBase
using Test
include("src/structs/SparseArray.jl")
include("src/utils/ML/spectralLinearRegression.jl")
=#
# include("src/structs/SparseArray.jl")
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


Hs.x[sortperm(unique(Hs.rowval[1:Hs.n_vals]))]

solveHuber!(Hs, _residuals_, _weights_, 
            Float32(1e9), #δ large so effectively squared error 
            Float32(0.0), #λ
            100,
            100,
            1000, 
            10.0f0,
            10.0f0,
            10.0f0,#Hs.n/10.0,
            0.01f0,
            L2Norm()
            );

Ht =  Matrix(transpose(H))
rowvals = copy(Hs.rowval)
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end
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

@test (cor(ŷ[:,], _weights_) - 1.0) < 1e-6
@test maximum(abs.(ŷ[:,] .- _weights_)) < 10
plot(log2.(ŷ[:,]), log2.(_weights_), seriestype=:scatter)

Ht =  Matrix(transpose(H))
rowvals = copy(Hs.rowval)
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end
z_l2 = reshape(y, (1, Hs.m))
ŷ_l2 = zeros(Float32, (1, Hs.n))

test_alg = NMF.CoordinateDescent{Float32}(
                maxiter = 1000,
                tol = 1e-6,
                α=100.0, 
                l₁ratio=0.0,
                update_H = false,
                regularization=:transformation)
                
NMF.solve!(test_alg, z_l2, ŷ_l2, Matrix(transpose(H)))
plot(log2.(ŷ[:,]), log2.(ŷ_l2[:,]), seriestype=:scatter)

_residuals_ = zeros(Float32, Hs.m)
_weights_2 = zeros(Float32, Hs.n)
initResiduals!(_residuals_, Hs, _weights_2)


Hs.x[sortperm(unique(Hs.rowval[1:Hs.n_vals]))]

solveHuber!(Hs, _residuals_, _weights_2, 
            Float32(1e9), #δ large so effectively squared error 
            Float32(100.0), #λ
            100,
            100,
            5,
            10.0f0,
            10.0f0,
            10.0f0,#Hs.n/10.0,
            0.01f0,
            L2Norm()
            );

plot(log2.(ŷ_l2[:,]), log2.(_weights_2), seriestype=:scatter)
plot(log2.(_weights_), log2.(_weights_2), seriestype=:scatter)

#=
Testing L2 norm with different lambda 
=#

_residuals_ = zeros(Float32, Hs.m)
_weights_ = zeros(Float32, Hs.n)
initResiduals!(_residuals_, Hs, _weights_)
solveHuber!(Hs, _residuals_, _weights_, 
            Float32(1e9), #δ large so effectively squared error 
            Float32(1e3), #λ
            1000,
            1000,
            1000, 
            0.01f0,
            0.01f0,
            0.001f0,#Hs.n/10.0,
            0.01f0
            );
argmax(_weights_ .- w_old)
(_weights_ .- w_old)[argmax(_weights_ .- w_old)]
Ht =  Matrix(transpose(H))
z = reshape(y, (1, Hs.m))
ŷ = zeros(Float32, (1, Hs.n))
#=
test_alg = NMF.CoordinateDescent{Float32}(
                maxiter = 1000,
                tol = 1e-6,
                α=1e3, 
                l₁ratio=1.0,
                update_H = false,
                regularization=:transformation)
NMF.solve!(test_alg, z, ŷ, Matrix(transpose(H)))


@test (cor(ŷ[:,], _weights_) - 1.0) < 1e-6
@test maximum(abs.(ŷ[:,] .- _weights_)) < 20
=#

maximum(abs.(ŷ[:,] .- _weights_)./(ŷ[:,].+1e-6))
plot(log2.(ŷ[:,]), log2.(_weights_), seriestype=:scatter)


plot(log2.(ŷ[:,]), log2.(_weights_), seriestype=:scatter)
plot!([log2.(ŷ[329])], [log2.(_weights_[329])], seriestype=:scatter)

Xs = sparse(X')

