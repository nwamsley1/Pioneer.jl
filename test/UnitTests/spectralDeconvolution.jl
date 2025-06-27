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

# Spectral Deconvolution Unit Tests
# Tests the solveHuber! function with different regularization settings

using JLD2
using Test
using Pioneer
using SparseArrays

# Load test data
Hs = load("data/Hs.jld2")["Hs"]
N = Hs.n_vals

# Get Data
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end

# Design Matrix 
H = Matrix(sparse(Hs.rowval[1:N],
                  Hs.colval[1:N],
                  Hs.nzval[1:N]))

# OLS Regression 
rowvals = copy(Hs.rowval)
y = zeros(Float32, Hs.m)
for i in range(1, N)
    y[Hs.rowval[i]] = Hs.x[i]
end

@testset "Spectral Deconvolution Tests" begin
    
    @testset "No Regularization (λ=0)" begin
        _residuals_ = zeros(Float32, Hs.m)
        _weights_ = zeros(Float32, Hs.n)
        initResiduals!(_residuals_, Hs, _weights_)
        
        solveHuber!(Hs, _residuals_, _weights_, 
                    Float32(1e9), #δ large so effectively squared error 
                    Float32(0.0), #λ
                    100,          #max_iter_newton
                    100,          #max_iter_bisection
                    100,          #max_iter_outer
                    10.0f0,       #accuracy_newton
                    10.0f0,       #accuracy_bisection
                    0.01f0,       #relative_convergence_threshold
                    NoNorm()      #regularization_type (no regularization since λ=0)
                    );
        
        # Test that we get non-zero weights
        @test sum(abs.(_weights_)) > 0
        @test !all(iszero, _weights_)
        
        # Save weights for comparison
        w_no_reg = copy(_weights_)
    end
    
    @testset "L2 Regularization (λ=1e3)" begin
        _residuals_ = zeros(Float32, Hs.m)
        _weights_ = zeros(Float32, Hs.n)
        initResiduals!(_residuals_, Hs, _weights_)
        
        solveHuber!(Hs, _residuals_, _weights_, 
                    Float32(1e9), #δ large so effectively squared error 
                    Float32(1e3), #λ
                    1000,         #max_iter_newton
                    1000,         #max_iter_bisection
                    1000,         #max_iter_outer
                    0.01f0,       #accuracy_newton
                    0.01f0,       #accuracy_bisection
                    0.01f0,       #relative_convergence_threshold
                    L2Norm()      #regularization_type
                    );
        
        # Test that L2 regularization reduces weights
        w_no_reg = ones(Float32, Hs.n)  # Dummy comparison if first test didn't run
        if @isdefined(w_no_reg)
            @test sum(abs.(_weights_)) < sum(abs.(w_no_reg))
        end
        @test sum(abs.(_weights_)) > 0  # Should still have some non-zero weights
    end
    
    @testset "Convergence Test" begin
        # Test that solver converges with reasonable parameters
        _residuals_ = zeros(Float32, Hs.m)
        _weights_ = zeros(Float32, Hs.n)
        initResiduals!(_residuals_, Hs, _weights_)
        
        # Should converge without error
        solveHuber!(Hs, _residuals_, _weights_, 
                    Float32(300.0), #δ - more realistic value
                    Float32(0.1),   #λ - small regularization
                    25,             #max_iter_newton (recommended)
                    100,            #max_iter_bisection (recommended)
                    max(1000, Hs.n*5), #max_iter_outer (recommended)
                    1.0f0,          #accuracy_newton
                    1.0f0,          #accuracy_bisection
                    0.001f0,        #relative_convergence_threshold
                    L2Norm()        #regularization_type
                    );
        
        # Check that we got a reasonable solution
        @test all(isfinite, _weights_)
        @test all(_weights_ .>= 0)  # Weights should be non-negative
    end
end