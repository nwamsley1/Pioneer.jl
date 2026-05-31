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

@testset "LightGBM top-8 raw fragment intensity features" begin
    for feature in (:frag1_int, :frag2_int, :frag3_int, :frag4_int,
                    :frag5_int, :frag6_int, :frag7_int, :frag8_int)
        @test feature in Pioneer.PRESCORE_FEATURES
        @test feature in Pioneer.ADVANCED_FEATURE_SET
    end

    X = Pioneer.feature_matrix(
        DataFrame(frag7_int = Float32[0, 7], frag8_int = Float32[8, 0]),
        [:frag7_int, :frag8_int],
    )
    @test X == Float32[0 8; 7 0]
end

@testset "BitVec excess-rank companion features" begin
    tc = fill(100, 256)
    dc = fill(100, 256)
    tc[Int(0x03) + 1] = 400
    dc[Int(0x03) + 1] = 100
    tc[Int(0x01) + 1] = 250
    dc[Int(0x01) + 1] = 100

    ranks = Pioneer._bitvec_excess_rank_table(tc, dc, 1.0)
    @test ranks[Int(0x03) + 1] == UInt16(1)
    @test ranks[Int(0x01) + 1] == UInt16(2)
    @test all(==(UInt16(0)), ranks) == false

    @test :n_correlated_fragments_bitvec_rank in Pioneer.PRESCORE_FEATURES
    @test :n_correlated_fragments_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET
    @test :frag_corr_strength in Pioneer.PRESCORE_FEATURES
    @test :frag_corr_strength in Pioneer.ADVANCED_FEATURE_SET
    @test :frag_corr_effective_n in Pioneer.PRESCORE_FEATURES
    @test :frag_corr_effective_n in Pioneer.ADVANCED_FEATURE_SET
    @test !(:wide_n_correlated_fragments_bitvec_rank in Pioneer.PRESCORE_FEATURES)
    @test :wide_n_correlated_fragments_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET
    @test !(:wide_frag_corr_strength in Pioneer.PRESCORE_FEATURES)
    @test :wide_frag_corr_strength in Pioneer.ADVANCED_FEATURE_SET
    @test !(:wide_frag_corr_effective_n in Pioneer.PRESCORE_FEATURES)
    @test :wide_frag_corr_effective_n in Pioneer.ADVANCED_FEATURE_SET
    @test !(:n_frags_detected_union_bitvec_rank in Pioneer.PRESCORE_FEATURES)
    @test !(:n_frags_detected_intersection_bitvec_rank in Pioneer.PRESCORE_FEATURES)
    @test :n_frags_detected_union_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET
    @test :n_frags_detected_intersection_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(
        DataFrame(
            n_correlated_fragments_bitvec_rank = UInt16[7, 1],
            wide_n_correlated_fragments_bitvec_rank = UInt16[9, 2],
            n_frags_detected_union_bitvec_rank = UInt16[11, 3],
            n_frags_detected_intersection_bitvec_rank = UInt16[13, 4],
        ),
        [
            :n_correlated_fragments_bitvec_rank,
            :wide_n_correlated_fragments_bitvec_rank,
            :n_frags_detected_union_bitvec_rank,
            :n_frags_detected_intersection_bitvec_rank,
        ],
    )
    @test X == Float32[7 9 11 13; 1 2 3 4]
end

@testset "fragment pattern rank features use bitvec mask ranks" begin
    rank_table = fill(UInt16(256), 256)
    rank_table[Int(0x03) + 1] = UInt16(7)

    psms = DataFrame(
        precursor_idx = fill(UInt32(1), 3),
        frag1_int = Float32[1, 2, 1],
        frag2_int = Float32[2, 4, 2],
        frag3_int = Float32[1, 0, 1],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[0, 0, 0],
        frag6_int = Float32[0, 0, 0],
        frag7_int = Float32[0, 0, 0],
        frag8_int = Float32[0, 0, 0],
        weight = Float32[1, 2, 1],
        irt_obs = Float32[1, 2, 3],
    )

    Pioneer._add_fragment_chromatogram_features!(psms; bitvec_rank_table = rank_table)
    @test psms.n_correlated_fragments == UInt8[2, 2, 2]
    @test psms.n_correlated_fragments_bitvec_rank == UInt16[7, 7, 7]
    @test psms.frag_corr_strength ≈ Float32[2, 2, 2]
    @test psms.frag_corr_effective_n ≈ Float32[2, 2, 2]

    features = Pioneer._wide_window_feature_values(
        Float32[0, 0, 0],
        Float32[1 2 0 0 0 0 0 0;
                2 4 0 0 0 0 0 0;
                1 2 0 0 0 0 0 0],
        Bool[true, true, true];
        bitvec_rank_table = rank_table,
    )
    @test features.wide_n_correlated_fragments == UInt8(2)
    @test features.wide_n_correlated_fragments_bitvec_rank == UInt16(7)
    @test features.wide_frag_corr_strength ≈ 2f0
    @test features.wide_frag_corr_effective_n ≈ 2f0
end

@testset "LightGBM isotope-trace collapse features" begin
    @test !(:precursor_fraction_transmitted in Pioneer.PRESCORE_FEATURES)
    @test !(:n_scans_other_traces in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_weight_corr in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_frag_sum_corr in Pioneer.PRESCORE_FEATURES)
    @test !(:trace_other_apex_delta_irt in Pioneer.PRESCORE_FEATURES)
    @test !(:n_frags_detected_union in Pioneer.PRESCORE_FEATURES)
    @test !(:n_frags_detected_intersection in Pioneer.PRESCORE_FEATURES)
    @test !(:n_frags_detected_union_bitvec_rank in Pioneer.PRESCORE_FEATURES)
    @test !(:n_frags_detected_intersection_bitvec_rank in Pioneer.PRESCORE_FEATURES)
    @test !(:frag_observed_sum_spectral_angle in Pioneer.PRESCORE_FEATURES)
    @test !(:frag_observed_sum_hellinger in Pioneer.PRESCORE_FEATURES)
    @test :precursor_fraction_transmitted in Pioneer.ADVANCED_FEATURE_SET
    @test :n_scans_other_traces in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_weight_corr in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_frag_sum_corr in Pioneer.ADVANCED_FEATURE_SET
    @test :trace_other_apex_delta_irt in Pioneer.ADVANCED_FEATURE_SET
    @test :n_frags_detected_union in Pioneer.ADVANCED_FEATURE_SET
    @test :n_frags_detected_intersection in Pioneer.ADVANCED_FEATURE_SET
    @test :n_frags_detected_union_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET
    @test :n_frags_detected_intersection_bitvec_rank in Pioneer.ADVANCED_FEATURE_SET
    @test !(:frag_observed_sum_spectral_angle in Pioneer.ADVANCED_FEATURE_SET)
    @test :frag_observed_sum_hellinger in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(
        DataFrame(
            precursor_fraction_transmitted = Float32[0.25, 0.75],
            n_scans_other_traces = UInt32[0, 3],
            trace_other_weight_corr = Float32[-1, 0.9],
            trace_other_frag_sum_corr = Float32[-1, 0.8],
            trace_other_apex_delta_irt = Float32[100, 0.2],
            n_frags_detected_union = UInt8[2, 6],
            n_frags_detected_intersection = UInt8[1, 4],
            n_frags_detected_union_bitvec_rank = UInt16[20, 5],
            n_frags_detected_intersection_bitvec_rank = UInt16[30, 6],
            frag_observed_sum_hellinger = Float32[0.5, 0.1],
        ),
        [
            :precursor_fraction_transmitted,
            :n_scans_other_traces,
            :trace_other_weight_corr,
            :trace_other_frag_sum_corr,
            :trace_other_apex_delta_irt,
            :n_frags_detected_union,
            :n_frags_detected_intersection,
            :n_frags_detected_union_bitvec_rank,
            :n_frags_detected_intersection_bitvec_rank,
            :frag_observed_sum_hellinger,
        ],
    )
    @test X == Float32[
        0.25 0 -1 -1 100 2 1 20 30 0.5
        0.75 3 0.9 0.8 0.2 6 4 5 6 0.1
    ]
end
