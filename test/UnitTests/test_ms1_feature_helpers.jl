@testset "MS1 isolation-window helper features" begin
    @test Pioneer._ms1_m0_m1_m2_window_fraction(10f0, 5f0, 0f0, 100f0) == 0.15f0
    @test Pioneer._ms1_m0_m1_m2_window_fraction(10f0, 5f0, 0f0, 0f0) == 0f0
    @test Pioneer._ms1_m0_m1_m2_window_fraction(100f0, 100f0, 100f0, 10f0) == 1f0

    frac_pc = Pioneer._ms1_m0_m1_m2_window_fraction_pseudocount(
        10f0, 0f0, 0f0,
        true, false, false,
        100f0,
        2f0,
    )
    @test frac_pc ≈ Float32(14 / 106)

    @test Pioneer._ms1_isotope_dotp_m0_m1_m2(10f0, 5f0, 0f0, 2f0, 1f0, 0f0) ≈ 1f0
    @test Pioneer._ms1_isotope_dotp_m0_m1_m2(0f0, 0f0, 0f0, 2f0, 1f0, 0f0) == 0f0
    @test Pioneer._ms2_explained_fraction(-1f0) ≈ 0.5f0
    @test Pioneer._ms2_explained_fraction(Inf32) == 0f0
    @test Pioneer._ms1_ms2_explained_delta(0.75f0, -1f0) ≈ 0.25f0
end

@testset "MS1 peak lookup advances from previous cursor" begin
    mz = Float32[100.0003, 100.5020, 101.0034, 102.0]
    intens = Float32[10, 20, 30, 40]

    m0_hit, m0_int, m0_mz, m0_idx, m0_cursor =
        Pioneer._ms1_find_peak_from_cursor(mz, intens, 100.0000f0, 10f0, 1)
    @test m0_hit
    @test m0_int == 10f0
    @test m0_mz == 100.0003f0
    @test m0_idx == 1
    @test m0_cursor == 1

    miss_hit, miss_int, miss_mz, miss_idx, miss_cursor =
        Pioneer._ms1_find_peak_from_cursor(mz, intens, 100.2500f0, 1f0, m0_cursor)
    @test !miss_hit
    @test miss_int == 0f0
    @test miss_mz == 0f0
    @test miss_idx == 0
    @test miss_cursor == 2

    m1_hit, m1_int, m1_mz, m1_idx, m1_cursor =
        Pioneer._ms1_find_peak_from_cursor(mz, intens, 100.5020f0, 10f0, miss_cursor)
    @test m1_hit
    @test m1_int == 20f0
    @test m1_mz == 100.5020f0
    @test m1_idx == 2
    @test m1_cursor == 2

    floor_cursor = Pioneer.bsearch_hybrid(mz, 100.0f0, 1, length(mz))
    ceiling_cursor = Pioneer.bsearch_hybrid(mz, nextfloat(101.0034f0), floor_cursor, length(mz)) - 1
    lower_hit, lower_int, lower_mz, lower_idx, lower_cursor =
        Pioneer._ms1_find_peak_from_bounds(mz, intens, 100.5020f0, 10f0, floor_cursor, ceiling_cursor)
    @test lower_hit
    @test lower_int == 20f0
    @test lower_mz == 100.5020f0
    @test lower_idx == 2
    @test lower_cursor == 2

    outside_hit, _, _, _, _ =
        Pioneer._ms1_find_peak_from_bounds(mz, intens, 102.0f0, 10f0, floor_cursor, ceiling_cursor)
    @test !outside_hit
end

@testset "MS1 M0 peak competition scan-run features" begin
    peak_a = (UInt64(10) << 32) | UInt64(2)
    peak_b = (UInt64(10) << 32) | UInt64(7)

    frac_out = zeros(Float32, 4)
    peak_count_out = zeros(UInt16, 4)
    same_mz_count_out = zeros(UInt16, 4)
    precursor_idx = UInt32[1, 2, 1, 3]
    prec_mzs = Float32[500, 600, 500]
    m0_peak_keys = UInt64[peak_a, peak_a, peak_a, peak_b]
    frag_cols = (
        Float32[10, 30, 10, 5],
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
    )
    scratch = Pioneer._M0PeakCompetitionScratch()

    Pioneer._add_m0_peak_fragment_competition_run!(
        frac_out,
        peak_count_out,
        same_mz_count_out,
        precursor_idx,
        prec_mzs,
        m0_peak_keys,
        frag_cols,
        1:4,
        scratch,
    )

    @test frac_out ≈ Float32[0.2, 0.6, 0.2, 1.0]
    @test peak_count_out == UInt16[2, 2, 2, 1]
    @test same_mz_count_out == UInt16[2, 1, 2, 2]
end
