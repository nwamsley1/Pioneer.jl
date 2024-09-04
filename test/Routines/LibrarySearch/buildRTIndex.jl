function Tol(a, b, ppm = 2)
    abs(a-b)<=(ppm*minimum((a, b))/1000000)
end

@testset "buildFragmentIndex.jl" begin

    RTs = Float32[5, 5.1, 5.2, 5.2, 6, 6.1, 6.2, 6.2, 6.2]
    prec_mzs = Float32[110, 100, 101, 102, 100, 121, 102, 102, 102]
    prec_ids = [x for x in 1:length(RTs)]

    rt_index = buildRtIndex(RTs, prec_mzs, prec_ids, 0.5)

    @test length(rt_index.rt_bins) == 2
    @test [length(x.prec) for x in rt_index.rt_bins] == [4, 5]
    @test all([issorted(rt_bin.prec, by = x->last(x)) for rt_bin in rt_index.rt_bins])

    RTs = Float32[5, 5.1, 5.2, 5.2, 6, 6.1, 6.2, 6.2, 10.0]
    prec_mzs = Float32[100, 100, 101, 102, 100, 101, 102, 102, 102]
    prec_ids = [x for x in 1:length(RTs)]

    rt_index = buildRtIndex(RTs, prec_mzs, prec_ids, 0.5)

    @test length(rt_index.rt_bins) == 3
    @test [length(x.prec) for x in rt_index.rt_bins] == [4, 4, 1]
    @test all([issorted(rt_bin.prec, by = x->last(x)) for rt_bin in rt_index.rt_bins])

    #Try different RT bin width
    rt_index = buildRtIndex(RTs, prec_mzs, prec_ids, 1.0)
    @test length(rt_index.rt_bins) == 3
    @test [length(x.prec) for x in rt_index.rt_bins] == [5, 3, 1]
    @test all([issorted(rt_bin.prec, by = x->last(x)) for rt_bin in rt_index.rt_bins])

    #Try different RT bin width
    rt_index = buildRtIndex(RTs, prec_mzs, prec_ids, 20.0)
    @test length(rt_index.rt_bins) == 1
    @test [length(x.prec) for x in rt_index.rt_bins] == [9]
    @test all([issorted(rt_bin.prec, by = x->last(x)) for rt_bin in rt_index.rt_bins])

end