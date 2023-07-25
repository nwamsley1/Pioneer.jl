@testset "counter.jl" begin
    test_counter = Counter(UInt32, UInt8, Float32, 5)

    for i in 1:10
        inc!(test_counter, UInt32(4), Float32(500.0))
        inc!(test_counter, UInt32(5), Float32(1000.0))
    end

    @test test_counter.size == 3
    @test test_counter.matches == 0
    countFragMatches(test_counter, 1, Float32(1.0));
    @test test_counter.matches == 2

    @test test_counter.ids == UInt32[4, 5, 0, 0, 0]
    sort!(test_counter, 2)
    @test test_counter.ids == UInt32[5, 4, 0, 0, 0]

    test_counter = Counter(UInt32, UInt8, Float32, 5)

    for i in 1:10
        inc!(test_counter, UInt32(4), Float32(500.0))
        inc!(test_counter, UInt32(5), Float32(1000.0))
        inc!(test_counter, UInt32(1), Float32(250.0))
        inc!(test_counter, UInt32(2), Float32(750.0))
    end

    @test test_counter.ids == UInt32[4, 5, 1, 2, 0]
    sort!(test_counter, 2)
    @test test_counter.ids == UInt32[4, 5, 1, 2, 0]
    countFragMatches(test_counter, 1, Float32(1.0));
    sort!(test_counter, 2)
    @test test_counter.ids == UInt32[5, 2, 4, 1, 0]
    @test test_counter.matches == 4
    @test test_counter.size == 5
    @test test_counter.counts == [ (0x0a, 2500.0),
    (0x0a, 7500.0),
    (0x00, 0.0),
    (0x0a, 5000.0),
    (0x0a, 10000.0)]
    reset!(test_counter)
    @test test_counter.matches == 0
    @test test_counter.size == 1
    @test test_counter.counts == [ (0x00, 0.0),
                            (0x00, 0.0),
                            (0x00, 0.0),
                            (0x00, 0.0),
                            (0x00, 0.0)]

end