using Test
using Pioneer
using DataFrames

const TEST_FLANKING_WINDOW_FEATURES = [
    :flanking_ms1_m0_candidate_fraction,
    :flanking_frag_candidate_fraction,
    :flanking_ms1_frag_sum_corr,
    :flanking_frag_corr_mean,
    :flanking_n_correlated_fragments,
    :flanking_n_correlated_fragments_bitvec_rank,
    :flanking_frag_corr_strength,
    :flanking_frag_corr_effective_n,
    :flanking_frag_corr_best_m0,
    :flanking_signal_support,
    :n_contiguous_scans,
    :frag_apex_gt2x_flank_bitvec_rank,
]

function _reference_flanking_pcor(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n < 2 && return 0f0
    mean_x = 0f0
    mean_y = 0f0
    @inbounds for i in 1:n
        mean_x += Float32(x[i])
        mean_y += Float32(y[i])
    end
    inv_n = 1f0 / Float32(n)
    mean_x *= inv_n
    mean_y *= inv_n

    sum_xx = 0f0
    sum_yy = 0f0
    sum_xy = 0f0
    @inbounds for i in 1:n
        dx = Float32(x[i]) - mean_x
        dy = Float32(y[i]) - mean_y
        sum_xx += dx * dx
        sum_yy += dy * dy
        sum_xy += dx * dy
    end
    denom = sqrt(sum_xx * sum_yy)
    return denom > 0f0 ? Float32(sum_xy / denom) : 0f0
end

function _reference_flanking_feature_values(
    core_ms1_m0_signal::Float32,
    core_frag_signal::Float32,
    core_fragments::NTuple{8, Float32},
    flank_ms1_m0::AbstractVector,
    flank_fragments::AbstractMatrix;
    bitvec_rank_table = nothing,
)
    n_scans = length(flank_ms1_m0)
    n_scans == 0 && return Pioneer._wide_window_zero_feature_values(
        max(core_ms1_m0_signal, 0f0),
        max(core_frag_signal, 0f0),
    )

    n_frags = size(flank_fragments, 2)
    frag_sum = zeros(Float32, n_scans)
    ms1_flank_signal = 0f0
    frag_flank_signal = 0f0
    signal_scans = 0

    @inbounds for scan_i in 1:n_scans
        m0 = max(Float32(flank_ms1_m0[scan_i]), 0f0)
        summed_frag = 0f0
        for frag_i in 1:n_frags
            summed_frag += max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        frag_sum[scan_i] = summed_frag
        ms1_flank_signal += m0
        frag_flank_signal += summed_frag
        (m0 > 0f0 || summed_frag > 0f0) && (signal_scans += 1)
    end

    has_signal = falses(n_frags)
    @inbounds for frag_i in 1:n_frags
        for scan_i in 1:n_scans
            if Float32(flank_fragments[scan_i, frag_i]) > 0f0
                has_signal[frag_i] = true
                break
            end
        end
    end

    pair_corr = zeros(Float32, n_frags, n_frags)
    pair_sum = 0f0
    pair_count = 0
    @inbounds for frag_a in 1:n_frags
        has_signal[frag_a] || continue
        for frag_b in (frag_a + 1):n_frags
            has_signal[frag_b] || continue
            corr = _reference_flanking_pcor(view(flank_fragments, :, frag_a), view(flank_fragments, :, frag_b))
            pair_corr[frag_a, frag_b] = corr
            pair_corr[frag_b, frag_a] = corr
            pair_sum += corr
            pair_count += 1
        end
    end

    n_correlated = UInt8(0)
    correlated_mask = UInt16(0)
    apex_mask = UInt16(0)
    corr_strength = 0f0
    corr_sumsq = 0f0
    other_sum = Vector{Float32}(undef, n_scans)
    @inbounds for frag_i in 1:n_frags
        flank_avg = 0f0
        for scan_i in 1:n_scans
            flank_avg += max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        flank_avg /= Float32(n_scans)
        core_apex = max(core_fragments[frag_i], 0f0)
        if core_apex > 0f0 && core_apex >= 2f0 * flank_avg
            apex_mask |= UInt16(1) << (frag_i - 1)
        end

        has_signal[frag_i] || continue
        for scan_i in 1:n_scans
            other_sum[scan_i] = frag_sum[scan_i] - max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        corr = _reference_flanking_pcor(view(flank_fragments, :, frag_i), other_sum)
        positive_corr = min(max(corr, 0f0), 1f0)
        corr_strength += positive_corr
        corr_sumsq += positive_corr * positive_corr
        if corr > Pioneer.WIDE_WINDOW_CORR_THRESHOLD
            n_correlated += UInt8(1)
            correlated_mask |= UInt16(1) << (frag_i - 1)
        end
    end

    best_frag = 0
    best_consensus = typemin(Float32)
    @inbounds for frag_i in 1:n_frags
        has_signal[frag_i] || continue
        consensus = 0f0
        n_pairs = 0
        for other_frag_i in 1:n_frags
            (other_frag_i == frag_i || !has_signal[other_frag_i]) && continue
            consensus += pair_corr[frag_i, other_frag_i]
            n_pairs += 1
        end
        avg_consensus = n_pairs > 0 ? Float32(consensus / n_pairs) : typemin(Float32)
        if avg_consensus > best_consensus
            best_consensus = avg_consensus
            best_frag = frag_i
        end
    end

    return (
        flanking_ms1_m0_candidate_fraction = Pioneer._wide_candidate_signal_fraction(max(core_ms1_m0_signal, 0f0), ms1_flank_signal),
        flanking_frag_candidate_fraction = Pioneer._wide_candidate_signal_fraction(max(core_frag_signal, 0f0), frag_flank_signal),
        flanking_ms1_frag_sum_corr = _reference_flanking_pcor(flank_ms1_m0, frag_sum),
        flanking_frag_corr_mean = pair_count > 0 ? Float32(pair_sum / pair_count) : 0f0,
        flanking_n_correlated_fragments = n_correlated,
        flanking_n_correlated_fragments_bitvec_rank = Pioneer._bitvec_pattern_rank(bitvec_rank_table, correlated_mask),
        flanking_frag_corr_strength = corr_strength,
        flanking_frag_corr_effective_n = corr_sumsq > 0f0 ? Float32((corr_strength * corr_strength) / corr_sumsq) : 0f0,
        flanking_frag_corr_best_m0 = best_frag > 0 ? _reference_flanking_pcor(view(flank_fragments, :, best_frag), flank_ms1_m0) : 0f0,
        flanking_signal_support = Float32(signal_scans / n_scans),
        frag_apex_gt2x_flank_bitvec_rank = Pioneer._bitvec_pattern_rank(bitvec_rank_table, apex_mask),
    )
end

@testset "flanking-window features are cross-run only" begin
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, TEST_FLANKING_WINDOW_FEATURES)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), TEST_FLANKING_WINDOW_FEATURES)
    @test collect(Pioneer.FLANKING_WINDOW_FEATURES) == TEST_FLANKING_WINDOW_FEATURES
end

@testset "flanking-core bounds use contiguous passing cycles around the best score" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[101, 102, 103, 104],
        cycle_idx = UInt32[1, 2, 3, 4],
        lgbm_score = Float32[0.1, 0.9, 0.8, 0.2],
        weight = ones(Float32, 4),
        irt_obs = Float32[1, 2, 3, 4],
        ms1_m0_intensity = Float32[1, 10, 20, 40],
        frag1_int = Float32[0, 3, 5, 0],
        frag2_int = Float32[0, 1, 2, 0],
        frag3_int = zeros(Float32, 4),
        frag4_int = zeros(Float32, 4),
        frag5_int = zeros(Float32, 4),
        frag6_int = zeros(Float32, 4),
        frag7_int = zeros(Float32, 4),
        frag8_int = zeros(Float32, 4),
    )
    best = Pioneer.select_best_per_precursor!(psms, :lgbm_score)
    rank_table = fill(UInt16(99), 256)
    rank_table[0x00 + 1] = UInt16(1)

    Pioneer.add_trace_and_fragment_features!(
        best,
        psms,
        Bool[false, true, true, false];
        bitvec_rank_table = rank_table,
    )

    @test hasproperty(best, :flanking_core_scan_min)
    @test hasproperty(best, :flanking_core_scan_max)
    @test hasproperty(best, :n_contiguous_scans)
    @test !hasproperty(best, :flanking_core_n_scans)
    @test hasproperty(best, :flanking_core_ms1_m0_signal)
    @test hasproperty(best, :flanking_core_frag_signal)
    @test best.flanking_core_scan_min == UInt32[102]
    @test best.flanking_core_scan_max == UInt32[103]
    @test best.n_contiguous_scans == UInt16[2]
    @test best.flanking_core_ms1_m0_signal == Float32[30]
    @test best.flanking_core_frag_signal == Float32[11]
end

@testset "flanking-window scans exclude core scans" begin
    key = (Int32(50000), Int32(200))
    scans = Int32[100, 105, 110, 115, 120, 125]
    scan_to_window_key = fill((Int32(0), Int32(0)), 125)
    scan_to_window_pos = zeros(Int32, 125)
    for (pos, scan) in pairs(scans)
        scan_to_window_key[Int(scan)] = key
        scan_to_window_pos[Int(scan)] = Int32(pos)
    end
    scan_index = (
        window_scans = Dict(key => scans),
        scan_to_window_key = scan_to_window_key,
        scan_to_window_pos = scan_to_window_pos,
    )

    flank_scans = Int32[]
    Pioneer._wide_flank_window_scans_from_core!(
        flank_scans,
        scan_index,
        UInt32(105),
        UInt32(110),
    )
    @test flank_scans == Int32[100, 115, 120, 125]

    Pioneer._wide_group_flank_window_scans_from_core!(
        flank_scans,
        scan_index,
        UInt32[105, 115],
        UInt32[110, 120],
        UInt32[105, 115],
        1,
        2,
    )
    @test flank_scans == Int32[100, 125]
end

@testset "flanking-window scan work is grouped by MS2 scan" begin
    ms2_work = Dict{Int32, Vector{Tuple{Int, Int}}}()

    Pioneer._wide_add_group_ms2_work!(ms2_work, 1, Int32[100, 105])
    Pioneer._wide_add_group_ms2_work!(ms2_work, 2, Int32[105, 110])

    @test sort(collect(keys(ms2_work))) == Int32[100, 105, 110]
    @test ms2_work[Int32(100)] == [(1, 1)]
    @test ms2_work[Int32(105)] == [(1, 2), (2, 1)]
    @test ms2_work[Int32(110)] == [(2, 2)]
end

@testset "flanking-window MS1 work is grouped by nearest MS1 scan" begin
    scan_to_ms1 = zeros(Int32, 110)
    scan_to_ms1[100] = Int32(50)
    scan_to_ms1[105] = Int32(50)
    scan_to_ms1[110] = Int32(55)
    scan_index = (scan_to_ms1 = scan_to_ms1,)
    ms1_work = Dict{Int32, Vector{Tuple{Int, Int}}}()

    Pioneer._wide_add_group_ms1_work!(ms1_work, scan_index, 1, Int32[100, 105])
    Pioneer._wide_add_group_ms1_work!(ms1_work, scan_index, 2, Int32[105, 110])

    @test sort(collect(keys(ms1_work))) == Int32[50, 55]
    @test ms1_work[Int32(50)] == [(1, 1), (1, 2), (2, 1)]
    @test ms1_work[Int32(55)] == [(2, 2)]
end

@testset "flanking-window feature values combine core signal with flank support" begin
    flank_ms1_m0 = Float32[4, 1]
    flank_fragments = zeros(Float32, 2, 8)
    flank_fragments[:, 1] .= Float32[3, 2]
    flank_fragments[:, 2] .= Float32[1, 0.5]
    rank_table = fill(UInt16(99), 256)
    rank_table[0x03 + 1] = UInt16(7)
    core_fragments = (6f0, 2f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0)

    features = Pioneer._wide_flank_feature_values(
        40f0,
        28f0,
        core_fragments,
        flank_ms1_m0,
        flank_fragments;
        bitvec_rank_table = rank_table,
    )

    @test features.flanking_ms1_m0_candidate_fraction ≈ Float32(40 / 45)
    @test features.flanking_frag_candidate_fraction ≈ Float32(28 / 34.5)
    @test features.flanking_ms1_frag_sum_corr > 0.9f0
    @test features.flanking_frag_corr_mean > 0.99f0
    @test features.flanking_n_correlated_fragments == UInt8(2)
    @test features.flanking_n_correlated_fragments_bitvec_rank == UInt16(7)
    @test features.flanking_frag_corr_strength > 1.9f0
    @test features.flanking_frag_corr_effective_n > 1.99f0
    @test features.flanking_frag_corr_best_m0 > 0.9f0
    @test features.flanking_signal_support == 1f0
end

@testset "flanking-window fragment apex bitvec rank uses individual flanks" begin
    flank_ms1_m0 = Float32[0, 0]
    flank_fragments = zeros(Float32, 2, 8)
    flank_fragments[:, 1] .= Float32[2, 2]
    flank_fragments[:, 2] .= Float32[3, 3]
    flank_fragments[:, 3] .= Float32[0, 0]
    core_fragments = (4f0, 5f0, 1f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    rank_table = fill(UInt16(99), 256)
    rank_table[0x05 + 1] = UInt16(11)

    features = Pioneer._wide_flank_feature_values(
        0f0,
        0f0,
        core_fragments,
        flank_ms1_m0,
        flank_fragments;
        bitvec_rank_table = rank_table,
    )

    @test features.frag_apex_gt2x_flank_bitvec_rank == UInt16(11)
end

@testset "flanking-window feature values match PR418 reference" begin
    flank_ms1_m0 = Float32[1.0f7, 1.0f7 + 12f0, 1.0f7 + 24f0, 1.0f7 + 36f0]
    flank_fragments = zeros(Float32, 4, 8)
    flank_fragments[:, 1] .= Float32[2.0f7, 2.0f7 + 8f0, 2.0f7 + 16f0, 2.0f7 + 24f0]
    flank_fragments[:, 2] .= Float32[1.5f7, 1.5f7 + 4f0, 1.5f7 + 8f0, 1.5f7 + 12f0]
    flank_fragments[:, 3] .= Float32[3.0f7 + 30f0, 3.0f7 + 20f0, 3.0f7 + 10f0, 3.0f7]
    core_fragments = (2.0f7 + 48f0, 1.0f0, 8.0f7, 0f0, 0f0, 0f0, 0f0, 0f0)
    rank_table = fill(UInt16(99), 256)

    expected = _reference_flanking_feature_values(
        40f0,
        28f0,
        core_fragments,
        flank_ms1_m0,
        flank_fragments;
        bitvec_rank_table = rank_table,
    )
    observed = Pioneer._wide_flank_feature_values(
        40f0,
        28f0,
        core_fragments,
        flank_ms1_m0,
        flank_fragments;
        bitvec_rank_table = rank_table,
    )

    for name in keys(expected)
        @test getproperty(observed, name) == getproperty(expected, name)
    end
end
