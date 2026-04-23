using Test
using Pioneer
using DataStructures: OrderedDict

@testset "PioneerPredict parses dict-like Koina responses" begin
    rt_response = OrderedDict(
        "outputs" => Any[
            OrderedDict(
                "name" => "rt",
                "data" => Any[1.5, 2.5],
            ),
        ],
    )
    rt_batch = Pioneer.parse_koina_batch(Pioneer.RetentionTimeModel("chronologer"), rt_response)
    @test rt_batch.fragments.rt == Float32[1.5, 2.5]
    @test rt_batch.frags_per_precursor == 1

    ms2_response = OrderedDict(
        "outputs" => Any[
            OrderedDict(
                "name" => "intensities",
                "shape" => Any[1, 2],
                "data" => Any[10.0, 20.0],
            ),
            OrderedDict(
                "name" => "mz",
                "shape" => Any[1, 2],
                "data" => Any[100.0, 200.0],
            ),
            OrderedDict(
                "name" => "annotation",
                "shape" => Any[1, 2],
                "data" => Any["y1", "y2"],
            ),
        ],
    )
    ms2_batch = Pioneer.parse_koina_batch(Pioneer.InstrumentSpecificModel("unispec"), ms2_response)
    @test ms2_batch.frags_per_precursor == 2
    @test ms2_batch.fragments.intensities == Float32[10.0, 20.0]
    @test ms2_batch.fragments.mz == Float32[100.0, 200.0]
    @test ms2_batch.fragments.annotation == ["y1", "y2"]
end
