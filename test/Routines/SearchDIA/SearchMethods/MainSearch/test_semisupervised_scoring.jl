using Test
using DataFrames

@testset "MainSearch semi-supervised LightGBM helpers" begin
    targets = vcat(fill(true, 6), Bool[false, true, false])
    precursor_idx = UInt32.(vcat(1:6, 100, 7, 101))
    scores = collect(range(1.0, 0.10; length=length(targets)))

    initial_mask = Pioneer._mainsearch_training_mask(targets, nothing)
    @test all(initial_mask)

    ss_mask = Pioneer._mainsearch_training_mask(targets, scores)
    @test ss_mask == Bool[true, true, true, true, true, true, true, false, true]
    @test Pioneer._mainsearch_label_counts(targets, ss_mask) == (positives = 6, negatives = 2)
    @test Pioneer._mainsearch_scale_pos_weight(targets, ss_mask) ≈ 2 / 6

    counts = Pioneer._mainsearch_target_fdr_counts(scores, targets, precursor_idx)
    @test counts[Float32(0.01)] == (target_psms = 6, unique_target_precursors = 6)
    @test counts[Float32(0.05)] == (target_psms = 6, unique_target_precursors = 6)
    @test counts[Float32(0.10)] == (target_psms = 6, unique_target_precursors = 6)
end

@testset "MainSearch prescore scan-count features" begin
    scan_cycle_idx = zeros(UInt32, 17)
    scan_idxs = [2, 5, 8, 11, 14, 17]
    center_mzs = Union{Missing, Float32}[missing for _ in 1:17]
    isolation_widths = Union{Missing, Float32}[missing for _ in 1:17]
    center_mzs[scan_idxs] = Float32[400, 500, 400, 500, 400, 500]
    isolation_widths[scan_idxs] .= 25.0f0

    Pioneer._fill_ms2_window_cycle_indices!(scan_cycle_idx, scan_idxs, center_mzs, isolation_widths)
    @test scan_cycle_idx[scan_idxs] == UInt32[1, 1, 2, 2, 3, 3]

    psms = DataFrame(
        precursor_idx = UInt32[1, 1, 1, 1, 2, 2, 2, 2, 2, 3],
        scan_idx = UInt32[100, 200, 350, 900, 10, 20, 21, 50, 90, 500],
        cycle_idx = UInt32[1, 3, 4, 5, 10, 11, 11, 12, 14, 30],
    )

    Pioneer._add_mainsearch_prescore_scan_features!(psms)

    @test psms[!, :n_scans] == UInt32[4, 4, 4, 4, 5, 5, 5, 5, 5, 1]
    @test psms[!, :n_consecutive_cycles] == UInt32[1, 3, 3, 3, 3, 3, 3, 3, 1, 1]
    @test :n_consecutive_cycles in Pioneer.PRESCORE_FEATURES
    @test !hasproperty(psms, :max_n_consecutive_scans)
end

@testset "MainSearch best scan selection restores apex summary features" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 1, 1],
        lgbm_score = Float32[0.3, 0.9, 0.2],
        weight = Float32[1, 4, 1],
        rt = Float32[0, 1, 2],
        irt_obs = Float32[0, 1, 2],
        gof = Float16[2, 5, 3],
        fitted_manhattan_distance = Float16[1, 7, 4],
        fitted_spectral_contrast = Float16[0.1, 0.8, 0.5],
        matched_ratio = Float16[0.2, 0.4, 0.9],
        scribe = Float16[1, 6, 3],
        y_count = UInt8[2, 5, 3],
    )

    selected = Pioneer.select_best_per_precursor!(psms, :lgbm_score)

    @test nrow(selected) == 1
    @test !hasproperty(selected, :num_scans)
    @test selected.max_weight == Float32[4]
    @test selected.max_gof == Float16[5]
    @test selected.max_fitted_manhattan_distance == Float16[7]
    @test selected.max_fitted_spectral_contrast == Float16[0.8]
    @test selected.max_matched_ratio == Float16[0.9]
    @test selected.max_scribe == Float16[6]
    @test selected.y_ions_sum == UInt16[10]
    @test selected.max_y_ions == UInt16[5]
    @test isfinite(selected.smoothness[1])
    @test selected.smoothness[1] > 0
end

@testset "MainSearch feature aliases avoid duplicate old names" begin
    @test :best_rank in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :best_rank_iso in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :topn in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :topn_iso in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :spectral_contrast in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :fitted_spectral_contrast in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :matched_ratio in fieldnames(Pioneer.MainSearchScoredPSM)
    @test :total_ions_iso in Pioneer.ADVANCED_FEATURE_SET
    @test !(:isotope_count in Pioneer.ADVANCED_FEATURE_SET)
    @test :n_scans in Pioneer.ADVANCED_FEATURE_SET
    @test !(:num_scans in Pioneer.ADVANCED_FEATURE_SET)
end

@testset "MainSearch MS1 feature handoff schema" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 2],
        rt = Float32[10.0, 12.0],
        irt_obs = Float32[100.0, 120.0],
    )
    ms1_psms = DataFrame(
        precursor_idx = UInt32[1],
        m0 = Bool[true],
        n_iso = UInt8[3],
        big_iso = UInt8[3],
        m0_error = Float16[0.5],
        error = Float16[1.5],
        spectral_contrast = Float16[0.9],
        fitted_spectral_contrast = Float16[0.8],
        gof = Float16[4.0],
        max_matched_residual = Float16[5.0],
        max_unmatched_residual = Float16[0.0],
        fitted_manhattan_distance = Float16[3.0],
        matched_ratio = Float16[2.0],
        weight = Float32[100.0],
        ms_file_idx = UInt32[1],
        scan_idx = UInt32[7],
        rt = Float32[10.25],
        rt_max_intensity = Float32[10.5],
        rt_diff_max_intensity = Float32[0.5],
        pair_idx = UInt32[42],
    )

    Pioneer._join_ms1_features!(psms, ms1_psms, x -> Float32(10x))

    @test psms.ms1_features_missing == Bool[false, true]
    @test psms.n_iso_ms1 == UInt8[3, 0]
    @test psms.big_iso_ms1 == UInt8[3, 0]
    @test psms.error_ms1 == Float16[1.5, -1.0]
    @test psms.rt_ms1 == Float32[10.25, -1.0]
    @test psms.ms1_ms2_rt_diff == Float32[2.5, -1.0]
    @test hasproperty(psms, :weight_ms1)
    @test all(f -> hasproperty(psms, f), filter(f -> occursin("ms1", String(f)), Pioneer.ADVANCED_FEATURE_SET))
end

@testset "MainSearch MS1 isotope scan scoring" begin
    isotopes = Pioneer.Isotope{Float32}[
        Pioneer.Isotope(500.0f0, 0.8f0, UInt8(1), UInt32(1)),
        Pioneer.Isotope(500.5017f0, 0.2f0, UInt8(2), UInt32(1)),
        Pioneer.Isotope(501.0034f0, 0.05f0, UInt8(3), UInt32(1)),
    ]
    scan_mz = Union{Missing, Float32}[500.0001f0, 500.5016f0, 700.0f0]
    scan_int = Union{Missing, Float32}[800.0f0, 200.0f0, 10.0f0]
    corrected = Float32[]
    obs_low = Float32[]
    obs_high = Float32[]

    score = Pioneer._score_ms1_candidate_scan(
        isotopes,
        scan_mz,
        scan_int,
        Pioneer.MassErrorModel(0.0f0, (20.0f0, 20.0f0)),
        corrected,
        obs_low,
        obs_high,
    )

    @test score.has_features
    @test score.m0
    @test score.n_iso == UInt8(2)
    @test score.big_iso == UInt8(2)
    @test score.weight > 900.0f0
    @test abs(score.m0_error) < Float16(1.0)
    @test isfinite(score.gof)
    @test isfinite(score.fitted_spectral_contrast)
end

@testset "MainSearch semi-supervised LightGBM smoke" begin
    n = 80
    target = Bool[i % 4 != 0 for i in 1:n]
    psms = DataFrame(
        precursor_idx = UInt32.(1:n),
        target = target,
        cv_fold = UInt8.(i % 2 for i in 1:n),
        x = Float32[target[i] ? 1 + 0.01i : 0.01i for i in 1:n],
    )

    scores, _, info = Pioneer.train_mainsearch_semisupervised_lgbm(psms; features=[:x])

    @test length(scores) == n
    @test all(isfinite, scores)
    @test info.iterations == 1
    @test info.winner == "lgbm_semisupervised"
end
