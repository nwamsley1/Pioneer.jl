using Test
using DataFrames
using Pioneer

function _test_rank_weights(n::Int)
    weights = Float32[1f0 / sqrt(Float32(r)) for r in 1:n]
    return weights .* (Float32(n) / sum(weights))
end

function _test_weighted_effective_n(weights)
    strength = sum(weights)
    sumsq = sum(w -> w * w, weights)
    return Float32((strength * strength) / sumsq)
end

function _test_weighted_std(values, weights)
    sum_w = sum(weights)
    mean_v = sum(values .* weights) / sum_w
    var_v = sum(weights .* ((values .- mean_v) .^ 2)) / sum_w
    return Float32(sqrt(var_v))
end

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

@testset "LightGBM top-10 raw fragment intensity features" begin
    for feature in (:frag1_int, :frag2_int, :frag3_int, :frag4_int,
                    :frag5_int, :frag6_int, :frag7_int, :frag8_int,
                    :frag9_int, :frag10_int)
        @test feature in Pioneer.PRESCORE_FEATURES
        @test feature in Pioneer.ADVANCED_FEATURE_SET
    end

    X = Pioneer.feature_matrix(
        DataFrame(frag9_int = Float32[0, 9], frag10_int = Float32[10, 0]),
        [:frag9_int, :frag10_int],
    )
    @test X == Float32[0 10; 9 0]
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
    @test !(:n_anticorrelated_fragments in Pioneer.PRESCORE_FEATURES)
    @test !(:n_anticorrelated_fragments in Pioneer.ADVANCED_FEATURE_SET)
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
    weights = _test_rank_weights(8)
    expected_strength = weights[1] + weights[2]
    expected_effective_n = _test_weighted_effective_n(weights[1:2])
    @test psms.n_correlated_fragments == UInt8[2, 2, 2]
    @test psms.n_correlated_fragments_bitvec_rank == UInt16[7, 7, 7]
    @test psms.frag_corr_strength ≈ fill(expected_strength, 3)
    @test psms.frag_corr_effective_n ≈ fill(expected_effective_n, 3)

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

@testset "fragment peak candidate counts use scan-local unique precursors" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 2, 2, 3],
        scan_idx = UInt32[1, 1, 1, 2],
        frag1_peak_idx = UInt32[101, 101, 101, 101],
        frag2_peak_idx = UInt32[201, 0, 201, 201],
        frag3_peak_idx = zeros(UInt32, 4),
        frag4_peak_idx = zeros(UInt32, 4),
        frag5_peak_idx = zeros(UInt32, 4),
        frag6_peak_idx = zeros(UInt32, 4),
        frag7_peak_idx = zeros(UInt32, 4),
        frag8_peak_idx = zeros(UInt32, 4),
    )

    Pioneer._add_fragment_peak_candidate_counts!(psms)
    @test psms.frag1_peak_candidate_count == UInt16[2, 2, 2, 1]
    @test psms.frag2_peak_candidate_count == UInt16[2, 1, 2, 1]
    for rank in 3:8
        @test psms[!, Symbol("frag$(rank)_peak_candidate_count")] == ones(UInt16, 4)
    end
end

@testset "fragment correlation summaries weight by fragment uniqueness" begin
    psms = DataFrame(
        precursor_idx = fill(UInt32(1), 3),
        frag1_int = Float32[1, 2, 1],
        frag2_int = Float32[2, 4, 2],
        frag3_int = Float32[0, 0, 0],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[0, 0, 0],
        frag6_int = Float32[0, 0, 0],
        frag7_int = Float32[0, 0, 0],
        frag8_int = Float32[0, 0, 0],
        weight = Float32[1, 2, 1],
        irt_obs = Float32[1, 2, 3],
        frag1_peak_candidate_count = UInt16[1, 1, 1],
        frag2_peak_candidate_count = UInt16[2, 2, 2],
        frag3_peak_candidate_count = UInt16[1, 1, 1],
        frag4_peak_candidate_count = UInt16[1, 1, 1],
        frag5_peak_candidate_count = UInt16[1, 1, 1],
        frag6_peak_candidate_count = UInt16[1, 1, 1],
        frag7_peak_candidate_count = UInt16[1, 1, 1],
        frag8_peak_candidate_count = UInt16[1, 1, 1],
    )

    Pioneer._add_fragment_chromatogram_features!(psms)
    weights = _test_rank_weights(8)
    corr_weights = Float32[weights[1], weights[2] / 2f0]
    @test psms.n_correlated_fragments == UInt8[2, 2, 2]
    @test psms.frag_corr_strength ≈ fill(sum(corr_weights), 3)
    @test psms.frag_corr_effective_n ≈ fill(_test_weighted_effective_n(corr_weights), 3)
end

@testset "fragment correlation summaries include ranks 9 and 10" begin
    rank_table = fill(UInt16(1024), 1024)
    rank_table[Int(0x300) + 1] = UInt16(42)

    psms = DataFrame(
        precursor_idx = fill(UInt32(1), 3),
        frag1_int = Float32[0, 0, 0],
        frag2_int = Float32[0, 0, 0],
        frag3_int = Float32[0, 0, 0],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[0, 0, 0],
        frag6_int = Float32[0, 0, 0],
        frag7_int = Float32[0, 0, 0],
        frag8_int = Float32[0, 0, 0],
        frag9_int = Float32[5, 10, 5],
        frag10_int = Float32[3, 6, 3],
        weight = Float32[1, 2, 1],
        irt_obs = Float32[1, 2, 3],
        ms1_m0_intensity = Float32[7, 14, 7],
    )

    Pioneer._add_fragment_chromatogram_features!(psms; bitvec_rank_table = rank_table)
    weights = _test_rank_weights(10)
    expected_strength = weights[9] + weights[10]
    expected_effective_n = _test_weighted_effective_n(weights[9:10])
    @test psms.n_correlated_fragments == UInt8[2, 2, 2]
    @test psms.n_correlated_fragments_bitvec_rank == UInt16[42, 42, 42]
    @test psms.frag_corr_strength ≈ fill(expected_strength, 3)
    @test psms.frag_corr_effective_n ≈ fill(expected_effective_n, 3)
    @test psms.frag_corr_best_m0 ≈ Float32[1, 1, 1]
end

@testset "fragment apex dispersion is rank weighted" begin
    psms = DataFrame(
        precursor_idx = fill(UInt32(1), 3),
        frag1_int = Float32[9, 1, 1],
        frag2_int = Float32[0, 0, 0],
        frag3_int = Float32[0, 0, 0],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[0, 0, 0],
        frag6_int = Float32[0, 0, 0],
        frag7_int = Float32[0, 0, 0],
        frag8_int = Float32[0, 0, 0],
        frag9_int = Float32[0, 0, 0],
        frag10_int = Float32[1, 1, 9],
        weight = Float32[1, 1, 1],
        irt_obs = Float32[1, 2, 3],
    )

    Pioneer._add_fragment_chromatogram_features!(psms)
    weights = _test_rank_weights(10)
    expected = _test_weighted_std(Float32[1, 3], Float32[weights[1], weights[10]])
    @test psms.frag_apex_dispersion_irt ≈ fill(expected, 3)
    @test only(unique(psms.frag_apex_dispersion_irt)) < 1f0
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
