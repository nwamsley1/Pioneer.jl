using Test
using DataFrames
using Pioneer

@testset "MS1 isolation-window signal features" begin
    mz = Float32[399, 400, 401, 402, 403]
    intens = Float32[5, 100, 50, 25, 25]

    @test Pioneer._ms1_sum_window_intensity(mz, intens, 399.5f0, 402.5f0) ≈ 175f0
    @test Pioneer._ms1_m0_m1_m2_window_fraction(100f0, 50f0, 25f0, 175f0) ≈ 1f0
    @test Pioneer._ms2_explained_fraction(log2(0.75f0)) ≈ 0.75f0
    @test Pioneer._ms1_ms2_explained_delta(Float32(150 / 175), log2(0.75f0)) ≈ Float32(150 / 175 - 0.75)
    @test Pioneer._ms1_noise_floor(Float32[0, -1, 20, 5, 100]) == 5f0
    @test Pioneer._ms1_m0_m1_m2_window_fraction_pseudocount(0f0, 0f0, 0f0, false, false, false, 0f0, 5f0) ≈ 1f0
    @test Pioneer._ms1_m0_m1_m2_window_fraction_pseudocount(0f0, 0f0, 0f0, false, false, false, 1000f0, 5f0) ≈ Float32(15 / 1015)
    @test Pioneer._ms1_m0_m1_m2_window_fraction_pseudocount(100f0, 0f0, 0f0, true, false, false, 100f0, 5f0) ≈ Float32(110 / 115)
    @test Pioneer._ms1_m0_m1_m2_window_fraction_pseudocount(100f0, 50f0, 25f0, true, true, true, 175f0, 5f0) ≈ Float32(175 / 190)

    psms = DataFrame(log2_intensity_explained = Float16[log2(0.5f0), log2(0.75f0)])
    Pioneer._initialize_ms1_isolation_window_features!(psms)
    @test psms.ms1_m0_m1_m2_window_fraction == Float32[0, 0]
    @test psms.ms1_ms2_explained_delta == Float32[0, 0]
    @test psms.ms1_m0_m1_m2_window_fraction_pc == Float32[0, 0]
    @test psms.ms1_ms2_explained_delta_pc == Float32[0, 0]

    @test :ms1_m0_m1_m2_window_fraction in Pioneer.PRESCORE_FEATURES
    @test :ms1_ms2_explained_delta in Pioneer.PRESCORE_FEATURES
    @test :ms1_m0_m1_m2_window_fraction_pc in Pioneer.PRESCORE_FEATURES
    @test :ms1_ms2_explained_delta_pc in Pioneer.PRESCORE_FEATURES
    @test :ms1_m0_m1_m2_window_fraction in Pioneer.ADVANCED_FEATURE_SET
    @test :ms1_ms2_explained_delta in Pioneer.ADVANCED_FEATURE_SET
    @test :ms1_m0_m1_m2_window_fraction_pc in Pioneer.ADVANCED_FEATURE_SET
    @test :ms1_ms2_explained_delta_pc in Pioneer.ADVANCED_FEATURE_SET
end
