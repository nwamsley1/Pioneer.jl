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

using SavitzkyGolay
using Plots
@testset "integratePrecursors.jl" begin

    ##########
    #pad()
    ##########
    @test pad(Float32[1, 2, 3, 2, 1], one(Float32), zero(Float32)) == Float32[1, 1, 2, 3, 2, 1, 0]

    ##########
    #getSmoothingParams()
    ##########
    n_data_points, max_window, min_order = 14, 15, 3
    @test getSmoothingParams(n_data_points, max_window, min_order) == (5, 3)

    n_data_points, max_window, min_order = 21, 15, 3
    @test getSmoothingParams(n_data_points, max_window, min_order) == (7, 3)

    n_data_points, max_window, min_order = 50, 14, 3
    @test getSmoothingParams(n_data_points, max_window, min_order) == (15, 3)

    n_data_points, max_window, min_order = 50, 15, 3
    @test getSmoothingParams(n_data_points, max_window, min_order) == (15, 3)

    ##########
    #fillNonZero()!
    ##########
    non_zero = BitVector([false for x in 1:11])
    intensity = Float32[0,0,0,1,0,1,0,1,1,1,0]
    fillNonZero!(non_zero, intensity)
    @test non_zero == BitVector([1,1,0,1,0,1,0,1,1,1,1])
    @test intensity[non_zero] == Float32[0, 0, 1, 1, 1, 1, 1, 0]

    ##########
    #getZeroCrossings()
    ##########
    Y = Float32[1, 1, 1, -1, 10, 1, 1, -10, -1]
    X = Float32[1, 2, 3,  4,  5, 6, 7, 8, 9]
    cross_idx, cross_slope = getZeroCrossings(Y, X)
    @test cross_idx == Int64[3, 4, 7]
    @test cross_slope == Float32[-2, 11, -11]

    ##########
    #getSmoothDerivative
    ##########
    Y = Float32[0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 0, 0, 0]
    X = Float32[x for x in range(1, length(Y))]
    window, order = getSmoothingParams(length(Y), 15, 3)
    Y′ = getSmoothDerivative(X, Y, window, order)
    cross_idx, cross_slope = getZeroCrossings(Y′, X[1:end - 1])
    slope, left, right = getPeakBounds(Y, X, cross_idx, cross_slope, 1.0, 3)
    @test Tol(slope, 2/3.0)
    @test left == 2
    @test right == 19

    #####Try locating left and right bounds. 

    Y = Float32[0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 1, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 0]
    X = Float32[x for x in range(1, length(Y))]
    window, order = getSmoothingParams(length(Y), 15, 3)
    Y′ = getSmoothDerivative(X, Y, window, order)
    cross_idx, cross_slope = getZeroCrossings(Y′, X[1:end - 1])
    slope, left, right = getPeakBounds(Y, X, cross_idx, cross_slope, 1.0, 3)
    @test left == 2
    @test right == 16


    Y = Float32[1, 2, 3, 2, 3, 2, 3, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 1, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 0]
    X = Float32[x for x in range(1, length(Y))]
    window, order = getSmoothingParams(length(Y), 15, 3)
    Y′ = getSmoothDerivative(X, Y, window, order)
    cross_idx, cross_slope = getZeroCrossings(Y′, X[1:end - 1])
    slope, left, right, mid = getPeakBounds(Y, X, cross_idx, cross_slope, 1.0, 3)
    #plot(X, Y)
    #plot!(X[1:end - 1],Y′)
    #vline!([10])
    #vline!([23])
    @test left == 10
    @test right == 23

    slope, left, right, mid = getIntegrationBounds(X, Y, window, order, 1.0, 3)

    @test left == 10
    @test right == 23

    chrom = DataFrame([X, Y],:auto)
    rename!(chrom, Symbol.(["rt","weight"]))

    area, count, SN, slope, error, base_with, FWHM = integratePrecursor(chrom, isplot = true)
    println(FWHM)
    #@test Tol(FWHM, 7)
    @test Tol(area, 49)
    @test count == 14
end