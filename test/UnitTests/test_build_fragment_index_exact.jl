using Test
using Arrow

using Pioneer: SimpleFrag, buildFragmentIndex!

@testset "buildFragmentIndex! exact synthetic layout" begin
    test_dir = mktempdir()
    try
        frag_ions = SimpleFrag{Float32}[
            SimpleFrag{Float32}(200.0f0, UInt32(1), 600.0f0, 0.0f0, UInt8(2), UInt8(10)),
            SimpleFrag{Float32}(150.0f0, UInt32(4), 520.0f0, 5.0f0, UInt8(2), UInt8(40)),
            SimpleFrag{Float32}(100.0005f0, UInt32(3), 550.0f0, 0.75f0, UInt8(3), UInt8(30)),
            SimpleFrag{Float32}(100.0f0, UInt32(2), 500.0f0, 0.5f0, UInt8(2), UInt8(20)),
        ]

        buildFragmentIndex!(
            test_dir,
            frag_ions,
            10.0f0,
            1.0f0,
            index_name = "exact_",
        )

        fragments = collect(Arrow.Table(joinpath(test_dir, "exact_f_index_fragments.arrow")).IndexFragment)
        rt_bins = collect(Arrow.Table(joinpath(test_dir, "exact_f_index_rt_bins.arrow")).FragIndexBin)
        frag_bins = collect(Arrow.Table(joinpath(test_dir, "exact_f_index_fragment_bins.arrow")).FragIndexBin)

        @test length(rt_bins) == 2
        @test [(bin.lb, bin.ub, bin.first_bin, bin.last_bin) for bin in rt_bins] == [
            (0.0f0, 0.75f0, UInt32(1), UInt32(2)),
            (5.0f0, 5.0f0, UInt32(3), UInt32(3)),
        ]

        @test length(frag_bins) == 3
        @test [(bin.lb, bin.ub, bin.first_bin, bin.last_bin) for bin in frag_bins] == [
            (100.0f0, 100.0005f0, UInt32(1), UInt32(2)),
            (200.0f0, 200.0f0, UInt32(3), UInt32(3)),
            (150.0f0, 150.0f0, UInt32(4), UInt32(4)),
        ]

        @test [(f.prec_id, f.prec_mz, f.score, f.charge) for f in fragments] == [
            (UInt32(2), 500.0f0, UInt8(20), UInt8(2)),
            (UInt32(3), 550.0f0, UInt8(30), UInt8(3)),
            (UInt32(1), 600.0f0, UInt8(10), UInt8(2)),
            (UInt32(4), 520.0f0, UInt8(40), UInt8(2)),
        ]
    finally
        safe_rmdir(test_dir)
    end
end
