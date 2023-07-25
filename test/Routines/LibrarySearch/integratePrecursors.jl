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
    @test left == 1
    @test right == 19

    #####Try locating left and right bounds. 

    Y = Float32[0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 1, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 0]
    X = Float32[x for x in range(1, length(Y))]
    window, order = getSmoothingParams(length(Y), 15, 3)
    Y′ = getSmoothDerivative(X, Y, window, order)
    cross_idx, cross_slope = getZeroCrossings(Y′, X[1:end - 1])
    slope, left, right = getPeakBounds(Y, X, cross_idx, cross_slope, 1.0, 3)
    @test left == 1
    @test right == 16


    Y = Float32[1, 2, 3, 2, 3, 2, 3, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 1, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 0]
    X = Float32[x for x in range(1, length(Y))]
    window, order = getSmoothingParams(length(Y), 15, 3)
    Y′ = getSmoothDerivative(X, Y, window, order)
    cross_idx, cross_slope = getZeroCrossings(Y′, X[1:end - 1])
    slope, left, right = getPeakBounds(Y, X, cross_idx, cross_slope, 1.0, 3)
    plot(X, Y)
    plot!(X[1:end - 1],Y′)
    vline!([10])
    vline!([23])
    @test left == 10
    @test right == 23




   

end