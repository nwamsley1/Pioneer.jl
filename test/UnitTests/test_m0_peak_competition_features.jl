using Test
using DataFrames
using Pioneer

@testset "MS1 M0 peak fragment-intensity competition feature" begin
    psms = DataFrame(
        scan_idx = UInt32[10, 10, 10, 11, 11, 12, 12],
        frag1_int = Float32[100, 10, 5, 3, 0, 0, 0],
        frag2_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag3_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag4_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag5_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag6_int = Float32[0, 0, 0, 0, 0, 0, 0],
    )
    m0_peak_keys = UInt64[100, 100, 200, 100, 0, 300, 300]

    Pioneer._add_m0_peak_fragment_competition_feature!(psms, m0_peak_keys)

    @test :ms1_m0_peak_frag_intensity_fraction in propertynames(psms)
    @test psms.ms1_m0_peak_frag_intensity_fraction ≈ Float32[
        100 / 110,
        10 / 110,
        1,
        1,
        0,
        0,
        0,
    ]
    @test :ms1_m0_peak_frag_intensity_fraction in Pioneer.PRESCORE_FEATURES
    @test :ms1_m0_peak_frag_intensity_fraction in Pioneer.ADVANCED_FEATURE_SET
end
