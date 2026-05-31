@testset "wide-window cross-run feature values" begin
    ms1_m0 = Float32[4, 10, 20, 10, 1]
    fragments = Float32[
        3 1 0 0 0 0
        5 2 0 0 0 0
        10 4 0 0 0 0
        5 2 0 0 0 0
        2 1 0 0 0 0
    ]
    candidate_mask = Bool[false, true, true, true, false]

    features = Pioneer._wide_window_feature_values(ms1_m0, fragments, candidate_mask)

    @test features.wide_ms1_m0_candidate_fraction ≈ Float32(40 / 45)
    @test features.wide_frag_candidate_fraction ≈ Float32(28 / 35)
    @test features.wide_ms1_frag_sum_corr > 0.9f0
    @test features.wide_frag_corr_mean > 0.99f0
    @test features.wide_n_correlated_fragments == UInt8(2)
    @test features.wide_frag_corr_strength > 1.9f0
    @test features.wide_frag_corr_effective_n > 1.99f0
    @test features.wide_frag_corr_best_m0 > 0.9f0
    @test features.wide_signal_support == 1f0
end

@testset "wide-window features handle missing raw evidence" begin
    features = Pioneer._wide_window_feature_values(
        zeros(Float32, 4),
        zeros(Float32, 4, 6),
        Bool[true, false, true, false],
    )

    for value in values(features)
        @test isfinite(Float32(value))
        @test value == zero(value)
    end
end

@testset "wide-window expansion uses persisted core bounds" begin
    key = (Int32(50000), Int32(2400))
    scan_to_window_key = fill((Int32(0), Int32(0)), 120)
    scan_to_window_pos = zeros(Int32, 120)
    for (pos, scan) in enumerate(Int32[100, 105, 110, 115])
        scan_to_window_key[Int(scan)] = key
        scan_to_window_pos[Int(scan)] = Int32(pos)
    end
    scan_index = (
        window_scans = Dict(key => Int32[100, 105, 110, 115]),
        scan_to_window_key = scan_to_window_key,
        scan_to_window_pos = scan_to_window_pos,
    )
    tbl = (
        wide_core_scan_min = UInt32[100],
        wide_core_scan_max = UInt32[110],
    )

    scans = Pioneer._wide_candidate_scans_from_core(tbl, [1], scan_index, Int32[105])

    @test scans == Int32[100, 105, 110]
end

@testset "wide-window fragment lookup uses learned mass-error model" begin
    scan_mz = Union{Missing, Float32}[500.025f0]
    scan_int = Union{Missing, Float32}[123f0]
    corrected = Float32[]
    obs_low = Float32[]
    obs_high = Float32[]
    mem = Pioneer.MassErrorModel(50f0, (5f0, 5f0))

    n_peaks = Pioneer.prepare_scan_peaks!(corrected, obs_low, obs_high, mem, scan_mz, scan_int)

    @test Pioneer._wide_peak_intensity(scan_mz, scan_int, 500f0, 20f0) == 0f0
    @test Pioneer._wide_fragment_peak_intensity(
        corrected,
        obs_low,
        obs_high,
        scan_int,
        n_peaks,
        500f0,
        mem,
    ) == 123f0
end

@testset "wide-window fragment lookup uses monoisotopic peak only" begin
    frag = Pioneer.DetailedFrag{Float32}(
        UInt32(1),
        1000.0f0,
        Float16(1000.0),
        UInt16(0),
        true,
        false,
        false,
        false,
        UInt8(1),
        UInt8(7),
        UInt8(2),
        UInt8(1),
        UInt8(0),
    )

    m0 = Float32(Pioneer.getMz(frag))
    m1 = m0 + Pioneer.C13_C12_MASS_DIFF_F32 / Float32(Pioneer.getFragCharge(frag))
    scan_mz = Union{Missing, Float32}[m0, m1]
    scan_int = Union{Missing, Float32}[1500.0f0, 900.0f0]
    corrected = Float32[]
    obs_low = Float32[]
    obs_high = Float32[]
    mem = Pioneer.MassErrorModel(0.0f0, (5.0f0, 5.0f0))
    n_peaks = Pioneer.prepare_scan_peaks!(corrected, obs_low, obs_high, mem, scan_mz, scan_int)

    @test Pioneer._wide_fragment_mono_peak_intensity(
        scan_mz,
        corrected,
        obs_low,
        obs_high,
        scan_int,
        n_peaks,
        frag,
        mem,
        true,
    ) == 1500.0f0
end

@testset "wide-window learned fragment model defaults on" begin
    env_key = "PIONEER_WIDE_WINDOW_USE_LEARNED_FRAGMENT_MODEL"
    old_value = get(ENV, env_key, nothing)
    haskey(ENV, env_key) && delete!(ENV, env_key)
    try
        @test Pioneer._wide_window_use_learned_fragment_model()
        ENV[env_key] = "false"
        @test !Pioneer._wide_window_use_learned_fragment_model()
    finally
        if old_value === nothing
            haskey(ENV, env_key) && delete!(ENV, env_key)
        else
            ENV[env_key] = old_value
        end
    end
end

@testset "wide-window features are cross-run only" begin
    @test !(:wide_frag_corr_roll_peak_effective_n in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_frag_corr_roll_peak_to_mean_delta in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_frag_corr_roll_apex_delta_scan in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_frag_corr_roll_peak_effective_n in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_frag_corr_roll_peak_to_mean_delta in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_frag_corr_roll_apex_delta_scan in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_left_frag_corr_effective_n in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_right_frag_corr_effective_n in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_flank_max_frag_corr_effective_n in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_core_to_flank_effective_n_delta in Pioneer.WIDE_WINDOW_FEATURES)
    @test !(:wide_left_frag_corr_effective_n in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_right_frag_corr_effective_n in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_flank_max_frag_corr_effective_n in Pioneer.ADVANCED_FEATURE_SET)
    @test !(:wide_core_to_flank_effective_n_delta in Pioneer.ADVANCED_FEATURE_SET)
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, Pioneer.WIDE_WINDOW_FEATURES)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), Pioneer.WIDE_WINDOW_FEATURES)
end
