using Test
using DataFrames
using Pioneer

@testset "MS1 M0 peak fragment-intensity competition feature" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 2, 3, 4, 5, 6, 6],
        scan_idx = UInt32[10, 10, 10, 11, 11, 12, 12],
        frag1_int = Float32[100, 10, 5, 3, 0, 0, 0],
        frag2_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag3_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag4_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag5_int = Float32[0, 0, 0, 0, 0, 0, 0],
        frag6_int = Float32[0, 0, 0, 0, 0, 0, 0],
    )
    m0_peak_keys = UInt64[100, 100, 200, 100, 0, 300, 300]
    precursor_mzs = Float32[500, 501, 500, 502, 502, 503]

    Pioneer._add_m0_peak_fragment_competition_feature!(psms, m0_peak_keys, precursor_mzs)

    @test :ms1_m0_peak_frag_intensity_fraction in propertynames(psms)
    @test :ms1_m0_peak_n_precursors in propertynames(psms)
    @test :scan_prec_mz_n_precursors in propertynames(psms)
    @test psms.ms1_m0_peak_frag_intensity_fraction ≈ Float32[
        100 / 110,
        10 / 110,
        1,
        1,
        0,
        0,
        0,
    ]
    @test psms.ms1_m0_peak_n_precursors == UInt16[2, 2, 1, 1, 0, 1, 1]
    @test psms.scan_prec_mz_n_precursors == UInt16[2, 1, 2, 2, 2, 1, 1]
    @test !(:ms1_m0_peak_frag_intensity_fraction in Pioneer.PRESCORE_FEATURES)
    @test !(:ms1_m0_peak_n_precursors in Pioneer.PRESCORE_FEATURES)
    @test !(:scan_prec_mz_n_precursors in Pioneer.PRESCORE_FEATURES)
    @test :ms1_m0_peak_frag_intensity_fraction in Pioneer.ADVANCED_FEATURE_SET
    @test :ms1_m0_peak_n_precursors in Pioneer.ADVANCED_FEATURE_SET
    @test :scan_prec_mz_n_precursors in Pioneer.ADVANCED_FEATURE_SET
end

@testset "MS2 fragment peak competition features" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 2, 3, 4, 5],
        scan_idx = UInt32[10, 11, 10, 11, 10],
        frag1_peak_idx = UInt32[100, 100, 100, 0, 300],
        frag2_peak_idx = UInt32[200, 200, 400, 0, 300],
        frag3_peak_idx = UInt32[300, 0, 0, 0, 500],
        frag4_peak_idx = UInt32[0, 0, 0, 0, 0],
        frag5_peak_idx = UInt32[0, 0, 0, 0, 0],
        frag6_peak_idx = UInt32[0, 0, 0, 0, 0],
    )

    Pioneer._add_fragment_peak_competition_features!(psms)

    @test :frag_competition_num_unique_fragments in propertynames(psms)
    @test :frag_competition_mean_candidates in propertynames(psms)
    @test !(:frag_competition_max_candidates in propertynames(psms))
    @test psms.frag_competition_num_unique_fragments == UInt8[3, 2, 2, 0, 2]
    @test psms.frag_competition_mean_candidates ≈ Float32[
        5 / 3,
        1,
        1.5,
        0,
        1.5,
    ]
    Pioneer.drop_fragment_peak_index_columns!(psms)
    @test !any(c -> c in propertynames(psms), Pioneer.FRAGMENT_PEAK_INDEX_COLUMNS)
    @test :frag_competition_mean_candidates in propertynames(psms)
    @test !(:frag_competition_num_unique_fragments in Pioneer.PRESCORE_FEATURES)
    @test !(:frag_competition_mean_candidates in Pioneer.PRESCORE_FEATURES)
    @test !(:frag_competition_max_candidates in Pioneer.PRESCORE_FEATURES)
    @test :frag_competition_num_unique_fragments in Pioneer.ADVANCED_FEATURE_SET
    @test :frag_competition_mean_candidates in Pioneer.ADVANCED_FEATURE_SET
    @test !(:frag_competition_max_candidates in Pioneer.ADVANCED_FEATURE_SET)
end
