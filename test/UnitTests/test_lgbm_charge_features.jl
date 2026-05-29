using Test
using DataFrames
using Pioneer

@testset "LightGBM charge-state features" begin
    @test :charge in Pioneer.PRESCORE_FEATURES
    @test :charge in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(DataFrame(charge = UInt8[2, 3, 4]), [:charge])
    @test X[:, 1] == Float32[2, 3, 4]
end

@testset "LightGBM spectral contrast feature" begin
    @test :spectral_contrast in Pioneer.PRESCORE_FEATURES
    @test :spectral_contrast in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(
        DataFrame(spectral_contrast = Float32[0.0, 0.5, 1.0]),
        [:spectral_contrast],
    )
    @test X[:, 1] == Float32[0.0, 0.5, 1.0]
end

@testset "LightGBM isotope-trace collapse features" begin
    @test !(:precursor_fraction_transmitted in Pioneer.PRESCORE_FEATURES)
    @test !(:n_scans_other_traces in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_weight_corr in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_frag_sum_corr in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_apex_delta_irt in Pioneer.PRESCORE_FEATURES)
    @test :precursor_fraction_transmitted in Pioneer.ADVANCED_FEATURE_SET
    @test :n_scans_other_traces in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_weight_corr in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_frag_sum_corr in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_apex_delta_irt in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(
        DataFrame(
            precursor_fraction_transmitted = Float32[0.25, 0.75],
            n_scans_other_traces = UInt32[0, 3],
            trace_other_weight_corr = Float32[-1, 0.9],
            trace_other_frag_sum_corr = Float32[-1, 0.8],
            trace_other_apex_delta_irt = Float32[100, 0.2],
        ),
        [
            :precursor_fraction_transmitted,
            :n_scans_other_traces,
            :trace_other_weight_corr,
            :trace_other_frag_sum_corr,
            :trace_other_apex_delta_irt,
        ],
    )
    @test X == Float32[0.25 0 -1 -1 100; 0.75 3 0.9 0.8 0.2]
end
