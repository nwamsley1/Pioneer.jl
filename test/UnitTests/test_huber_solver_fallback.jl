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

using Test
using Pioneer: SparseArrayFused, HuberSolver, NoNorm, solve_deconvolution!,
               huber_max_weight, HUBER_DEFAULT_MAX_WEIGHT

# One precursor, 8 observed fragments + 2 predicted-but-unobserved, observed = h * w_true.
function one_column_problem(w_true::Float32)
    h = Float32[1e-3, 8e-4, 6e-4, 4.5e-4, 3e-4, 2e-4, 1.2e-4, 8e-5, 5e-5, 5e-5]
    y = Float32[i <= 8 ? h[i] * w_true : 0f0 for i in eachindex(h)]
    n = length(h)
    H = SparseArrayFused{Int64, Float32}(n, n, 1, collect(1:n), h, y, zeros(UInt8, n), zeros(UInt8, n), [1, n + 1])
    return H, h, y
end

function solve(H; max_weight = HUBER_DEFAULT_MAX_WEIGHT)
    w = zeros(Float32, 1); r = zeros(Float32, H.m)
    solver = HuberSolver(300f0, 0f0, 50, 100, 10f0, 10f0, NoNorm())
    solve_deconvolution!(solver, H, r, w, zeros(Float32, 1), zeros(Float32, H.m), zeros(Float32, H.m),
        1000, 0.01f0; max_weight = max_weight)
    return w[1], r
end

@testset "Huber fallback: diverging Newton steps" begin
    # |r| >> delta makes the Newton step overflow; the weight used to come back 0 with a
    # residual vector that no longer equalled H*w - y.
    for w_true in (1f9, 1f10, 1f11)
        H, h, y = one_column_problem(w_true)
        w, r = solve(H)
        @test isapprox(w, w_true; rtol = 0.02)
        @test maximum(abs.(r .- (h .* w .- y))) <= 1f-3 * w_true * maximum(h)
    end
    H, _, _ = one_column_problem(1f11)
    @test solve(H; max_weight = 1f10)[1] <= 1f10
end

@testset "huber_max_weight" begin
    @test huber_max_weight(1f5) == 1f10
    for bad in (missing, nothing, NaN32, Inf32, -Inf, 0f0, -1.0, 1e38)
        @test huber_max_weight(bad) == HUBER_DEFAULT_MAX_WEIGHT
    end
end
