using Test
using Pioneer
using DataFrames

function _test_fragment_rank_weights(n_frags::Integer)
    weights = [1.0f0 / sqrt(Float32(r)) for r in 1:n_frags]
    scale = Float32(n_frags) / sum(weights)
    return Float32.(weights .* scale)
end

@testset "fragment chromatogram positive correlation summaries" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10],
        frag1_int = Float32[1, 2, 3],
        frag2_int = Float32[2, 4, 6],
        frag3_int = Float32[3, 2, 1],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[1, 1, 1],
        frag6_int = Float32[0, 2, 4],
        frag7_int = Float32[3, 2, 1],
        frag8_int = Float32[2, 3, 4],
        weight = Float32[1, 2, 3],
        irt_obs = Float32[0, 1, 2],
        ms1_m0_intensity = Float32[1, 2, 3],
    )

    Pioneer._add_fragment_chromatogram_features!(psms)

    rank_weights = _test_fragment_rank_weights(8)
    expected_strength = rank_weights[1] + rank_weights[2] + rank_weights[6] + rank_weights[8]
    expected_effective_n = expected_strength^2 /
        (rank_weights[1]^2 + rank_weights[2]^2 + rank_weights[6]^2 + rank_weights[8]^2)

    @test psms.n_correlated_fragments == UInt8[4, 4, 4]
    @test psms.frag_corr_strength ≈ fill(expected_strength, 3)
    @test psms.frag_corr_effective_n ≈ fill(expected_effective_n, 3)
end

@testset "fragment correlated mask bitvec rank" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10],
        frag1_int = Float32[1, 2, 3],
        frag2_int = Float32[2, 4, 6],
        frag3_int = Float32[3, 2, 1],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[1, 1, 1],
        frag6_int = Float32[0, 2, 4],
        frag7_int = Float32[3, 2, 1],
        frag8_int = Float32[2, 3, 4],
        weight = Float32[1, 2, 3],
        irt_obs = Float32[0, 1, 2],
        ms1_m0_intensity = Float32[1, 2, 3],
    )
    rank_table = fill(UInt16(99), 256)
    rank_table[0xa3 + 1] = UInt16(7)

    Pioneer._add_fragment_chromatogram_features!(psms; bitvec_rank_table=rank_table)

    @test psms.n_correlated_fragments == UInt8[4, 4, 4]
    @test psms.n_correlated_fragments_bitvec_rank == UInt16[7, 7, 7]
end

@testset "post-filter fragment union and intersection bitvec ranks" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 20, 30, 30],
        scan_idx = UInt32[1, 2, 3, 4, 5, 6],
        cycle_idx = UInt32[1, 2, 3, 4, 5, 6],
        lgbm_score = Float32[0.1, 0.9, 0.3, 0.4, 0.8, 0.2],
        weight = ones(Float32, 6),
        irt_obs = Float32[1, 2, 3, 4, 5, 6],
        frag1_int = Float32[1, 0, 1, 0, 5, 0],
        frag2_int = Float32[0, 2, 2, 0, 0, 0],
        frag3_int = Float32[3, 3, 0, 4, 0, 0],
        frag4_int = Float32[0, 0, 0, 0, 6, 7],
        frag5_int = zeros(Float32, 6),
        frag6_int = zeros(Float32, 6),
        frag7_int = zeros(Float32, 6),
        frag8_int = zeros(Float32, 6),
    )
    rank_table = fill(UInt16(99), 256)
    rank_table[0x07 + 1] = UInt16(11)
    rank_table[0x00 + 1] = UInt16(22)
    rank_table[0x04 + 1] = UInt16(33)
    rank_table[0x09 + 1] = UInt16(44)
    rank_table[0x08 + 1] = UInt16(55)

    best = Pioneer.select_best_per_precursor!(psms, :lgbm_score)
    deleteat!(best, 2)
    Pioneer.add_trace_and_fragment_features!(
        best,
        psms,
        falses(nrow(psms));
        bitvec_rank_table = rank_table,
    )

    @test best.precursor_idx == UInt32[10, 30]
    @test best.n_frags_detected_union == UInt8[3, 2]
    @test best.n_frags_detected_intersection == UInt8[0, 1]
    @test best.n_frags_detected_union_bitvec_rank == UInt16[11, 44]
    @test best.n_frags_detected_intersection_bitvec_rank == UInt16[22, 55]
end

@testset "trace feature pass skips slow other-trace summaries" begin
    other_trace_features = (
        :n_scans_other_traces,
        :trace_other_weight_corr,
        :trace_other_frag_sum_corr,
        :trace_other_apex_delta_irt,
    )
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[101, 102, 103, 104],
        lgbm_score = Float32[0.8, 0.9, 0.7, 0.6],
        weight = Float32[1, 2, 2, 4],
        irt_obs = Float32[1.0, 2.0, 1.2, 2.2],
        rt = Float32[10, 11, 12, 13],
        cycle_idx = UInt32[1, 2, 1, 2],
        frag1_int = Float32[1, 2, 2, 4],
        frag2_int = zeros(Float32, 4),
        frag3_int = zeros(Float32, 4),
        frag4_int = zeros(Float32, 4),
        frag5_int = zeros(Float32, 4),
        frag6_int = zeros(Float32, 4),
        frag7_int = zeros(Float32, 4),
        frag8_int = zeros(Float32, 4),
    )
    rank_table = fill(UInt16(99), 256)
    rank_table[0x01 + 1] = UInt16(8)

    best = Pioneer.select_best_per_precursor!(psms, :lgbm_score)
    Pioneer.add_trace_and_fragment_features!(
        best,
        psms,
        trues(nrow(psms));
        bitvec_rank_table = rank_table,
    )

    @test best.precursor_idx == UInt32[10]
    @test all(feature -> !(feature in Pioneer.ADVANCED_FEATURE_SET), other_trace_features)
    @test all(feature -> !hasproperty(best, feature), other_trace_features)
    @test best.n_frags_detected_union == UInt8[1]
    @test best.n_frags_detected_intersection == UInt8[1]
    @test best.n_frags_detected_union_bitvec_rank == UInt16[8]
    @test best.n_frags_detected_intersection_bitvec_rank == UInt16[8]
end
