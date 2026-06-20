using Test
using Pioneer
using DataFrames

function _test_fragment_rank_weights(n_frags::Integer)
    weights = [1.0f0 / sqrt(Float32(r)) for r in 1:n_frags]
    scale = Float32(n_frags) / sum(weights)
    return Float32.(weights .* scale)
end

function _add_shadow_peak_helper_columns_fragment_tests!(psms::DataFrame)
    n = nrow(psms)
    for col in Pioneer.FITTED_FRAGMENT_INTENSITY_COLUMNS
        psms[!, col] = zeros(Float32, n)
    end
    for col in Pioneer.SHADOW_FRAGMENT_INTENSITY_COLUMNS
        psms[!, col] = zeros(Float32, n)
    end
    return psms
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

@testset "fragment chromatogram summaries group by precursor and isolation window" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[1, 2, 3, 4],
        frag1_int = Float32[1, 2, 10, 1],
        frag2_int = Float32[2, 4, 10, 1],
        frag3_int = zeros(Float32, 4),
        frag4_int = zeros(Float32, 4),
        frag5_int = zeros(Float32, 4),
        frag6_int = zeros(Float32, 4),
        frag7_int = zeros(Float32, 4),
        frag8_int = zeros(Float32, 4),
        weight = Float32[1, 2, 1, 2],
        irt_obs = Float32[1, 2, 3, 4],
        ms1_m0_intensity = Float32[1, 2, 1, 2],
    )
    center_mzs = Union{Missing, Float32}[500, 500, 600, 600]
    isolation_widths = Union{Missing, Float32}[4, 4, 4, 4]
    groups = Pioneer._build_precursor_window_groups(
        psms.precursor_idx,
        psms.scan_idx,
        center_mzs,
        isolation_widths,
    )

    Pioneer._add_fragment_chromatogram_features!(psms; groups=groups)

    @test psms.n_scans == UInt32[2, 2, 2, 2]
    @test psms.n_correlated_fragments == UInt8[2, 2, 0, 0]
    @test psms.frag_corr_strength[1] > psms.frag_corr_strength[3]
end

@testset "selected precursor shape features use selected isolation window" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[1, 2, 3, 4],
        lgbm_score = Float32[0.8, 0.9, 0.1, 0.2],
        weight = Float32[10, 20, 100, 100],
        irt_obs = Float32[1, 2, 100, 200],
        rt = Float32[10, 20, 1000, 2000],
    )
    center_mzs = Union{Missing, Float32}[500, 500, 600, 600]
    isolation_widths = Union{Missing, Float32}[4, 4, 4, 4]

    best = Pioneer.select_best_per_precursor!(
        psms,
        :lgbm_score;
        center_mzs = center_mzs,
        isolation_widths = isolation_widths,
    )

    @test best.scan_idx == UInt32[2]
    @test best.n_above_hm == UInt16[2]
    @test best.irt_fwhm == Float32[1]
    @test best.rt_fwhm == Float32[10]
    @test best.best_rt == Float32[20]
    @test best.smoothness ≈ Float32[0.0225]
end

@testset "post-filter fragment union and intersection bitvec ranks" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 20, 30, 30],
        scan_idx = UInt32[1, 2, 3, 4, 5, 6],
        cycle_idx = UInt32[1, 2, 3, 4, 5, 6],
        lgbm_score = Float32[0.1, 0.9, 0.3, 0.4, 0.8, 0.2],
        weight = ones(Float32, 6),
        irt_obs = Float32[1, 2, 3, 4, 5, 6],
        ms1_m0_intensity = zeros(Float32, 6),
        frag1_int = Float32[1, 0, 1, 0, 5, 0],
        frag2_int = Float32[0, 2, 2, 0, 0, 0],
        frag3_int = Float32[3, 3, 0, 4, 0, 0],
        frag4_int = Float32[0, 0, 0, 0, 6, 7],
        frag5_int = zeros(Float32, 6),
        frag6_int = zeros(Float32, 6),
        frag7_int = zeros(Float32, 6),
        frag8_int = zeros(Float32, 6),
    )
    _add_shadow_peak_helper_columns_fragment_tests!(psms)
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

@testset "other-window summaries use lowest-PEP other window" begin
    other_window_features = (
        :n_scans_other_windows,
        :other_window_weight_corr,
        :other_window_apex_delta_irt,
    )
    psms = DataFrame(
        precursor_idx = fill(UInt32(10), 6),
        scan_idx = UInt32[101, 102, 201, 202, 301, 302],
        lgbm_score = Float32[0.8, 0.9, 0.7, 0.6, 0.5, 0.4],
        weight = Float32[1, 2, 2, 4, 100, 1],
        irt_obs = Float32[1.0, 2.0, 2.2, 3.2, 9.0, 10.0],
        rt = Float32[10, 11, 12, 13, 14, 15],
        cycle_idx = UInt32[1, 2, 1, 2, 1, 2],
        ms1_m0_intensity = Float32[1, 2, 2, 4, 100, 1],
        frag1_int = Float32[1, 2, 2, 4, 100, 1],
        frag2_int = zeros(Float32, 6),
        frag3_int = zeros(Float32, 6),
        frag4_int = zeros(Float32, 6),
        frag5_int = zeros(Float32, 6),
        frag6_int = zeros(Float32, 6),
        frag7_int = zeros(Float32, 6),
        frag8_int = zeros(Float32, 6),
    )
    _add_shadow_peak_helper_columns_fragment_tests!(psms)
    center_mzs = Vector{Union{Missing, Float32}}(fill(missing, 302))
    isolation_widths = Vector{Union{Missing, Float32}}(fill(missing, 302))
    center_mzs[101:102] .= 500.0f0
    center_mzs[201:202] .= 502.0f0
    center_mzs[301:302] .= 504.0f0
    isolation_widths[[101, 102, 201, 202, 301, 302]] .= 2.0f0
    peps = Float32[0.2, 0.1, 0.03, 0.04, 0.5, 0.6]
    rank_table = fill(UInt16(99), 256)
    rank_table[0x01 + 1] = UInt16(8)

    best = Pioneer.select_best_per_precursor!(psms, :lgbm_score)
    Pioneer.add_trace_and_fragment_features!(
        best,
        psms,
        trues(nrow(psms));
        bitvec_rank_table = rank_table,
        center_mzs = center_mzs,
        isolation_widths = isolation_widths,
        pep_values = peps,
    )

    @test best.precursor_idx == UInt32[10]
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, other_window_features)
    @test best.n_scans_other_windows == UInt32[4]
    @test best.other_window_weight_corr ≈ Float32[1.0]
    @test best.other_window_apex_delta_irt ≈ Float32[1.2]
    @test best.n_frags_detected_union == UInt8[1]
    @test best.n_frags_detected_intersection == UInt8[1]
    @test best.n_frags_detected_union_bitvec_rank == UInt16[8]
    @test best.n_frags_detected_intersection_bitvec_rank == UInt16[8]
end
