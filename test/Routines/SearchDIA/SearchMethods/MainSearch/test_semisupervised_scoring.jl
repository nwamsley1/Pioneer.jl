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
