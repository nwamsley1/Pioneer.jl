function _with_debug_chromatogram_plot_state(f::Function)
    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_targets = copy(Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[])

    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = Set(UInt32[909488])
        withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
            f()
        end
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = old_targets
    end
end

@testset "default chromatogram debug targets include requested precursors" begin
    @test UInt32(304430) in Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[]
end

@testset "mainsearch PEP intervals are captured before best scan selection" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 1, 99, 1, 2],
        scan_idx = UInt32[10, 20, 25, 30, 40],
        rt = Float32[1.0, 2.0, 2.5, 3.0, 4.0],
        target = Bool[true, true, false, true, true],
        lgbm_score = Float32[0.99, 0.98, 0.97, 0.96, 0.95],
    )

    Pioneer.add_mainsearch_integration_interval_columns!(psms, :lgbm_score)

    @test psms.mainsearch_1pct_start_scan == UInt32[10, 10, 0, 10, 40]
    @test psms.mainsearch_1pct_stop_scan == UInt32[30, 30, 0, 30, 40]
    @test isequal(psms.mainsearch_1pct_start_rt, Float32[1.0, 1.0, NaN, 1.0, 4.0])
    @test isequal(psms.mainsearch_1pct_stop_rt, Float32[3.0, 3.0, NaN, 3.0, 4.0])
end

@testset "mainsearch integration interval defaults to fifty percent PEP" begin
    psms = DataFrame(
        precursor_idx = fill(UInt32(1), 6),
        scan_idx = UInt32.(1:6),
        rt = Float32.(1:6),
        target = Bool[true, true, false, true, true, true],
        lgbm_score = Float32[4, 3, 6, 5, 2, 1],
        weight = Float32[1, 1, 1, 100, 1, 1],
    )

    Pioneer.add_mainsearch_integration_interval_columns!(psms, :lgbm_score)

    @test psms.mainsearch_1pct_start_scan == fill(UInt32(4), 6)
    @test psms.mainsearch_1pct_stop_scan == fill(UInt32(6), 6)
end

@testset "mainsearch PEP evidence is retained for debug plots" begin
    rt = Float32.(1:10)
    scan_idx = UInt32.(1:10)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[10, 20, 500, 1000, 150, 100, 102, 90, 200, 300]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        4,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        mainsearch_1pct_start_scan = UInt32(3),
        mainsearch_1pct_stop_scan = UInt32(8),
        boundary_candidate_data = candidate_data,
    )

    @test debug_plot_data[].mainsearch_evidence_range == 3:8
    candidates = DataFrame(candidate_data[])
    @test nrow(candidates) > 0
    @test :shape_fit_loss in propertynames(candidates)
    @test :mainsearch_exclusion_penalty ∉ propertynames(candidates)
    @test :mainsearch_left_bound_delta ∉ propertynames(candidates)
    @test :mainsearch_right_bound_delta ∉ propertynames(candidates)
end

@testset "single-scan evidence still produces at least three integration points" begin
    rt = Float32.(1:5)
    scan_idx = UInt32.(1:5)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 0, 10, 0, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        3,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        mainsearch_1pct_start_scan = UInt32(3),
        mainsearch_1pct_stop_scan = UInt32(3),
    )

    @test debug_plot_data[].scan_range == 2:4
    @test debug_plot_data[].points_integrated == UInt32(1)
end

@testset "wide mainsearch PEP evidence does not override fallback bounds" begin
    rt = Float32.(1:19)
    scan_idx = UInt32.(1:19)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[
        0.4, 0.35, 0.3, 0.05, 0.0, 0.6, 0.3, 0.8, 1.6,
        7.6, 3.0, 6.0, 1.4, 0.4, 0.7, 1.5, 1.9, 1.4, 0.9,
    ]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)
    fallback_plot_data = Ref{Any}(nothing)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        10,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = fallback_plot_data,
    )
    Pioneer.reset!(state)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        10,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        mainsearch_1pct_start_scan = UInt32(first(scan_idx)),
        mainsearch_1pct_stop_scan = UInt32(last(scan_idx)),
        boundary_candidate_data = candidate_data,
    )

    @test debug_plot_data[].scan_range == fallback_plot_data[].scan_range
    candidates = DataFrame(candidate_data[])
    @test !any(
        (candidates.candidate_start_idx .== UInt16(first(scan_idx))) .&
        (candidates.candidate_stop_idx .== UInt16(last(scan_idx)))
    )
end

@testset "baseline subtraction does not extrapolate above smoothed signal" begin
    rt = Float32.(1:19)
    scan_idx = UInt32.(1:19)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[
        0.6, 0.2, 0.4, 0.1, 0.1, 1.2, 0.2, 0.05, 0.3,
        1.4, 6.2, 4.3, 4.8, 2.4, 1.5, 1.6, 1.8, 1.6, 1.2,
    ]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        11,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        mainsearch_1pct_start_scan = UInt32(1),
        mainsearch_1pct_stop_scan = UInt32(15),
    )

    scan_range = debug_plot_data[].scan_range
    @test all(debug_plot_data[].baseline_subtracted[scan_range] .<= debug_plot_data[].wh_smoothed[scan_range] .+ eps(Float32))
    @test debug_plot_data[].baseline_subtracted[first(scan_range)] == 0.0f0
end

@testset "baseline subtraction uses selected boundary endpoints" begin
    rt = Float32.(1:8)
    scan_idx = UInt32.(1:8)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 0.1, 1.0, 0.05, 0.0, 100.0, 50.0, 0.0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        3,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        forced_boundary_start_scan = UInt32(1),
        forced_boundary_stop_scan = UInt32(8),
        debug_plot_data = debug_plot_data,
    )

    @test debug_plot_data[].scan_range == 1:8
    @test debug_plot_data[].baseline_subtracted[1] == 0.0f0
    @test debug_plot_data[].baseline_subtracted[8] == 0.0f0
    @test debug_plot_data[].baseline_subtracted[6] >
        0.5f0 * debug_plot_data[].wh_smoothed[6]
end

@testset "mainsearch evidence does not force whole interval as candidate" begin
    rt = Float32.(1:19)
    scan_idx = UInt32.(1:19)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[
        0.6, 0.2, 0.4, 0.1, 0.1, 1.2, 0.2, 0.05, 0.3,
        1.4, 6.2, 4.3, 4.8, 2.4, 1.5, 1.6, 1.8, 1.6, 1.2,
    ]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    debug_plot_data = Ref{Any}(nothing)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        11,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        mainsearch_1pct_start_scan = UInt32(1),
        mainsearch_1pct_stop_scan = UInt32(15),
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    @test debug_plot_data[].mainsearch_evidence_range == 1:15
    @test !any(
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(15))
    )
    @test :shape_fit_loss in propertynames(candidates)
end

@testset "boundary candidates include full chromatogram valley combinations" begin
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 4, 5.5, 6.5, 7.5, 10, 7.5, 6.5, 5.5, 4, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        6,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        mainsearch_1pct_start_scan = UInt32(first(scan_idx)),
        mainsearch_1pct_stop_scan = UInt32(last(scan_idx)),
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    candidate_ranges = Set(
        (Int(row.candidate_start_idx), Int(row.candidate_stop_idx))
        for row in eachrow(candidates)
    )

    @test candidate_ranges == Set([(1, 11), (2, 10)])
    @test nrow(candidates) == length(candidate_ranges)
    @test count(candidates.is_fallback) == 1
    @test :candidate_category in propertynames(candidates)
    @test all(candidates.candidate_category .∈ Ref(Set([
        "fallback",
        "valley_combination",
    ])))
    @test !any(startswith.(candidates.candidate_category, "fw_"))
    category_by_range = Dict(
        (Int(row.candidate_start_idx), Int(row.candidate_stop_idx)) => row.candidate_category
        for row in eachrow(candidates)
    )
    @test category_by_range[(2, 10)] == "fallback"
    @test category_by_range[(1, 11)] == "valley_combination"

    @test :mainsearch_left_bound_delta ∉ propertynames(candidates)
    @test :mainsearch_right_bound_delta ∉ propertynames(candidates)
end

@testset "recovered outside peak contributes valley combination candidates" begin
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 4, 10, 5, 1, 5, 6, 3, 1, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)
    debug_plot_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        4,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        debug_plot_data = debug_plot_data,
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    fallback = only(candidates[candidates.is_fallback .& candidates.is_primary_apex, :])
    @test debug_plot_data[].scan_range == first(fallback.candidate_start_idx):last(fallback.candidate_stop_idx)
    @test any(
        (candidates.candidate_category .== "valley_combination") .&
        (candidates.candidate_start_idx .< fallback.candidate_start_idx) .&
        (candidates.candidate_stop_idx .> fallback.candidate_stop_idx)
    )
end

@testset "boundary valley candidates use nearest local minima per side" begin
    rt = Float32.(1:19)
    scan_idx = UInt32.(1:19)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 3, 0, 3, 0, 3, 0, 3, 0, 10, 0, 3, 0, 3, 0, 3, 0, 3, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        10,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    valley_candidates = candidates[
        (candidates.candidate_category .== "valley_combination") .& candidates.is_primary_apex,
        :,
    ]
    allowed_left = Set(UInt16[5, 7, 9])
    allowed_right = Set(UInt16[11, 13, 15])

    @test nrow(valley_candidates) > 0
    @test all(in.(valley_candidates.candidate_start_idx, Ref(allowed_left)))
    @test all(in.(valley_candidates.candidate_stop_idx, Ref(allowed_right)))
    @test length(unique(valley_candidates.candidate_start_idx)) <= 3
    @test length(unique(valley_candidates.candidate_stop_idx)) <= 3
end

@testset "boundary candidates include local apexes across deconvolution span" begin
    rt = Float32.(1:21)
    scan_idx = UInt32.(1:21)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[
        0, 0, 0, 0,
        1, 5, 1,
        0, 0, 0, 0, 0,
        4, 30, 120, 30, 4,
        0, 0, 0, 0,
    ]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        2,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        mainsearch_1pct_start_scan = UInt32(13),
        mainsearch_1pct_stop_scan = UInt32(17),
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    apex_scans = sort(unique(candidates.new_best_scan))

    @test nrow(candidates) > 0
    @test all(
        (candidates.candidate_start_idx .<= UInt16(2)) .&
            (UInt16(2) .<= candidates.candidate_stop_idx),
    )
    @test apex_scans == UInt32[2, 6, 15]
    @test any(candidates.is_primary_apex)
end

@testset "mainsearch PEP intervals merge across trace rows" begin
    psms = DataFrame(
        precursor_idx = UInt32[7, 7, 8],
        ms_file_idx = UInt32[3, 3, 3],
        trace_prob = Float32[0.6, 0.7, 0.8],
        mainsearch_1pct_start_scan = UInt32[120, 100, 300],
        mainsearch_1pct_stop_scan = UInt32[130, 150, 320],
        mainsearch_1pct_start_rt = Float32[1.2, 1.0, 3.0],
        mainsearch_1pct_stop_rt = Float32[1.3, 1.5, 3.2],
    )

    Pioneer._aggregate_trace_to_precursor_probs!(psms)

    @test psms.mainsearch_1pct_start_scan == UInt32[100, 100, 300]
    @test psms.mainsearch_1pct_stop_scan == UInt32[150, 150, 320]
    @test psms.mainsearch_1pct_start_rt == Float32[1.0, 1.0, 3.0]
    @test psms.mainsearch_1pct_stop_rt == Float32[1.5, 1.5, 3.2]
end

@testset "EGH boundary shape fit penalizes disconnected internal signal" begin
    clean_irt = Float32[-2, -1, 0, 1, 2]
    clean_signal = Float32[0, 0.35, 1.0, 0.35, 0]
    disconnected_irt = Float32[-2, -1, 0, 1, 2, 3, 4]
    disconnected_signal = Float32[0, 0.35, 1.0, 0, 0.7, 0.35, 0]

    clean = Pioneer.fit_egh_boundary_shape(clean_irt, clean_signal, 3)
    disconnected = Pioneer.fit_egh_boundary_shape(disconnected_irt, disconnected_signal, 3)

    @test clean.valid
    @test disconnected.valid
    @test clean.fit_loss < disconnected.fit_loss
    @test clean.overshoot_penalty == 0.0f0
    @test disconnected.fit_loss > clean.fit_loss * 2.0f0
end

@testset "EGH boundary shape fit loss uses normalized squared error" begin
    irt = Float32[-1, 0, 1]
    signal = Float32[0, 1, 1]
    widths = Pioneer._shape_point_widths(irt)
    loss = Pioneer._egh_fit_loss(irt, signal, 2, 1.0f0, widths, 1.0f0, 0.0f0)
    predicted = Float32[
        Pioneer.egh_boundary_value(irt[i], 0.0f0, 1.0f0, 1.0f0, 0.0f0)
        for i in eachindex(irt)
    ]
    residuals = (signal .- predicted) ./ 1.0f0
    expected_loss = sum(widths .* residuals .* residuals) / sum(widths)

    @test loss.fit_loss ≈ expected_loss
    @test loss.fit_loss > 0.1f0
end

@testset "shape boundary selection uses lowest deterministic shape score" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2],
        candidate_index = UInt16[1, 2, 1, 2],
        ms_file_idx = UInt16[1, 1, 1, 1],
        is_fallback = Bool[true, false, true, false],
        shape_model_valid = Bool[true, true, false, false],
        shape_fit_loss = Float32[0.4, 0.1, Inf32, Inf32],
        shape_edge_fraction = Float32[0.0, 0.0, 0.0, 0.0],
        shape_overshoot_penalty = Float32[0.0, 0.0, 0.0, 0.0],
        shape_sigma_irt = Float32[1.0, 1.0, 1.0, 1.0],
        shape_tau_irt = Float32[0.0, 0.0, 0.0, 0.0],
        shape_apex_irt = Float32[0.0, 0.0, 0.0, 0.0],
        shape_apex_height = Float32[10.0, 10.0, 8.0, 8.0],
        shape_deconvolution_area_fraction = Float32[0.5, 0.5, 0.5, 0.5],
    )

    Pioneer.score_boundary_candidates_by_shape!(candidates)
    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test candidates.boundary_shape_score[2] < candidates.boundary_shape_score[1]
    @test selected.candidate_index == UInt16[2, 1]
    @test selected.boundary_group_id == UInt64[1, 2]
end

@testset "temporary fallback boundary override selects fallback rows" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2],
        candidate_index = UInt16[1, 2, 1, 2],
        ms_file_idx = UInt16[1, 1, 1, 1],
        is_fallback = Bool[true, false, true, false],
        quant_trace_selected = Bool[true, true, true, true],
        boundary_model_score = Float32[0.1, 0.9, 0.2, 0.8],
    )

    selected = Pioneer.select_fallback_boundary_candidate_rows(candidates)

    @test selected.candidate_index == UInt16[1, 1]
    @test all(selected.is_fallback)
end

@testset "shape score best-margin normalizes shape loss against plausible alternatives" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1],
        candidate_index = UInt16[1, 2],
        ms_file_idx = UInt16[1, 1],
        is_fallback = Bool[false, false],
        shape_model_valid = Bool[true, true],
        shape_fit_loss = Float32[0.10, 0.08],
        shape_edge_fraction = Float32[0.0, 0.0],
        shape_overshoot_penalty = Float32[0.0, 0.0],
        shape_sigma_irt = Float32[1.0, 1.0],
        shape_tau_irt = Float32[0.0, 0.0],
        shape_apex_irt = Float32[0.0, 2.0],
        shape_apex_height = Float32[10.0, 10.0],
        shape_deconvolution_area_fraction = Float32[0.9, 0.1],
    )

    Pioneer.score_boundary_candidates_by_shape!(candidates)
    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test candidates.shape_unexplained_area_penalty[1] ≈ 0.1f0
    @test candidates.shape_unexplained_area_penalty[2] ≈ 0.9f0
    @test candidates.shape_normalized_fit_prior_score[1] ≈ 1.0f0
    @test candidates.shape_normalized_fit_prior_score[2] ≈ 0.0f0
    @test Pioneer.BOUNDARY_SHAPE_FIT_PRIOR_SCORE_WEIGHT >= 1.0f0
    @test candidates.boundary_shape_score[1] ≈
        Pioneer.BOUNDARY_SHAPE_FIT_PRIOR_SCORE_WEIGHT * 1.0f0 + 0.1f0
    @test candidates.boundary_shape_score[2] ≈
        Pioneer.BOUNDARY_SHAPE_FIT_PRIOR_SCORE_WEIGHT * 0.0f0 + 0.9f0
    @test candidates.boundary_shape_score[2] < candidates.boundary_shape_score[1]
    @test selected.candidate_index == UInt16[2]
end

@testset "shape score is not compressed by awful candidates" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 1],
        candidate_index = UInt16[1, 2, 3],
        ms_file_idx = UInt16[1, 1, 1],
        is_fallback = Bool[false, false, false],
        shape_model_valid = Bool[true, true, true],
        shape_fit_loss = Float32[0.10, 0.12, 4.0],
        shape_edge_fraction = Float32[0.0, 0.0, 0.0],
        shape_overshoot_penalty = Float32[0.0, 0.0, 0.0],
        shape_sigma_irt = Float32[1.0, 1.0, 1.0],
        shape_tau_irt = Float32[0.0, 0.0, 0.0],
        shape_apex_irt = Float32[0.0, 0.0, 0.0],
        shape_apex_height = Float32[10.0, 10.0, 10.0],
        shape_deconvolution_area_fraction = Float32[0.5, 1.0, 1.0],
    )

    Pioneer.score_boundary_candidates_by_shape!(candidates)
    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test candidates.shape_normalized_fit_prior_score[1] ≈ 0.0f0
    @test candidates.shape_normalized_fit_prior_score[2] ≈ 1.0f0
    @test candidates.shape_normalized_fit_prior_score[3] ≈ 1.0f0
    @test selected.candidate_index == UInt16[1]
end

@testset "shape score has no expected apex bounds penalty" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 1],
        candidate_index = UInt16[1, 2, 3],
        ms_file_idx = UInt16[1, 1, 1],
        is_fallback = Bool[false, false, false],
        shape_model_valid = Bool[true, true, true],
        shape_fit_loss = Float32[0.10, 0.095, 0.115],
        shape_edge_fraction = Float32[0.0, 0.0, 0.0],
        shape_overshoot_penalty = Float32[0.0, 0.0, 0.0],
        shape_sigma_irt = Float32[0.05, 0.05, 0.05],
        shape_tau_irt = Float32[0.0, 0.0, 0.0],
        shape_apex_irt = Float32[0.0, 0.2, 0.0],
        shape_apex_height = Float32[10.0, 10.0, 10.0],
        shape_deconvolution_area_fraction = Float32[1.0, 1.0, 1.0],
    )

    Pioneer.score_boundary_candidates_by_shape!(candidates)
    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test :shape_expected_apex_bounds_penalty ∉ propertynames(candidates)
    @test candidates.shape_normalized_fit_prior_score[2] == 0.0f0
    @test candidates.boundary_shape_score[2] < candidates.boundary_shape_score[1]
    @test selected.candidate_index == UInt16[2]
end

@testset "shape score ignores edge and overshoot penalties" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1],
        candidate_index = UInt16[1, 2],
        ms_file_idx = UInt16[1, 1],
        is_fallback = Bool[false, false],
        shape_model_valid = Bool[true, true],
        shape_fit_loss = Float32[0.05, 0.10],
        shape_edge_fraction = Float32[10.0, 0.0],
        shape_overshoot_penalty = Float32[10.0, 0.0],
        shape_sigma_irt = Float32[1.0, 1.0],
        shape_tau_irt = Float32[0.0, 0.0],
        shape_apex_irt = Float32[0.0, 0.0],
        shape_apex_height = Float32[10.0, 10.0],
        shape_deconvolution_area_fraction = Float32[0.5, 0.5],
    )

    Pioneer.score_boundary_candidates_by_shape!(candidates)
    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test candidates.boundary_shape_score[1] < candidates.boundary_shape_score[2]
    @test selected.candidate_index == UInt16[1]
end

@testset "shape sigma prior follows apex iRT trend" begin
    rows = NamedTuple[]
    for i in 1:180
        apex_irt = Float32(24 * (i - 1) / 179)
        sigma = Float32(0.08 + 0.0012 * (apex_irt - 12)^2)
        push!(rows, (;
            boundary_group_id = UInt64(i),
            candidate_index = UInt16(1),
            ms_file_idx = UInt16(1),
            is_fallback = true,
            shape_model_valid = true,
            shape_fit_loss = 0.0f0,
            shape_edge_fraction = 0.0f0,
            shape_overshoot_penalty = 0.0f0,
            shape_sigma_irt = sigma,
            shape_tau_irt = 0.0f0,
            shape_apex_irt = apex_irt,
            shape_apex_height = 1_000.0f0,
            shape_deconvolution_area_fraction = 0.5f0,
        ))
    end
    test_sigma = 0.24f0
    push!(rows, (;
        boundary_group_id = UInt64(10_001),
        candidate_index = UInt16(1),
        ms_file_idx = UInt16(1),
        is_fallback = false,
        shape_model_valid = true,
        shape_fit_loss = 0.0f0,
        shape_edge_fraction = 0.0f0,
        shape_overshoot_penalty = 0.0f0,
        shape_sigma_irt = test_sigma,
        shape_tau_irt = 0.0f0,
        shape_apex_irt = 1.0f0,
        shape_apex_height = 1_000.0f0,
        shape_deconvolution_area_fraction = 0.5f0,
    ))
    push!(rows, (;
        boundary_group_id = UInt64(10_002),
        candidate_index = UInt16(1),
        ms_file_idx = UInt16(1),
        is_fallback = false,
        shape_model_valid = true,
        shape_fit_loss = 0.0f0,
        shape_edge_fraction = 0.0f0,
        shape_overshoot_penalty = 0.0f0,
        shape_sigma_irt = test_sigma,
        shape_tau_irt = 0.0f0,
        shape_apex_irt = 12.0f0,
        shape_apex_height = 1_000.0f0,
        shape_deconvolution_area_fraction = 0.5f0,
    ))
    candidates = DataFrame(rows)

    Pioneer.score_boundary_candidates_by_shape!(candidates)

    edge_row = findfirst(==(UInt64(10_001)), candidates.boundary_group_id)
    middle_row = findfirst(==(UInt64(10_002)), candidates.boundary_group_id)
    @test candidates.shape_prior_penalty[edge_row] < candidates.shape_prior_penalty[middle_row] / 2
end

@testset "boundary shape parameter QC plots are written per run" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 2, 3, 4, 5],
        ms_file_idx = UInt16[1, 1, 1, 2, 2],
        is_fallback = Bool[true, true, true, true, false],
        shape_model_valid = Bool[true, true, true, true, true],
        shape_sigma_irt = Float32[0.22, 0.28, 0.31, 0.45, 0.8],
        shape_tau_irt = Float32[-0.05, 0.0, 0.08, 0.12, 0.9],
        shape_apex_irt = Float32[10, 20, 30, 15, 25],
        shape_apex_height = Float32[1_000, 4_000, 16_000, 8_000, 32_000],
    )
    tmp = mktempdir()

    paths = Pioneer.write_boundary_shape_parameter_qc_plots(
        candidates,
        tmp;
        file_name_by_idx = Dict(UInt16(1) => "run A.raw", UInt16(2) => "run B.raw"),
    )

    @test length(paths) == 2
    @test all(isfile, paths)
    @test any(occursin("file_1_run_A_raw", basename(path)) for path in paths)
    @test any(occursin("file_2_run_B_raw", basename(path)) for path in paths)
    @test all(path -> filesize(path) > 0, paths)
end

@testset "boundary selected category tally logs at debug level" begin
    selected = DataFrame(
        candidate_category = [
            "fallback",
            "fallback",
            "valley_combination",
        ],
    )

    lines = Pioneer.boundary_candidate_category_tally_lines(selected)
    @test first(lines) == "Chromatogram boundary shape-model selected candidate categories (total 3):"
    @test any(occursin("fallback", line) && occursin("2", line) for line in lines)
    @test any(occursin("valley_combination", line) && occursin("1", line) for line in lines)

    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_debug_file = Pioneer.DEBUG_FILE[]
    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
        silent_path = tempname()
        open(silent_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            Pioneer.log_boundary_candidate_category_tally(selected)
        end
        @test isempty(read(silent_path, String))

        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        debug_path = tempname()
        open(debug_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            redirect_stdout(devnull) do
                Pioneer.log_boundary_candidate_category_tally(selected)
            end
        end
        @test occursin(
            "Chromatogram boundary shape-model selected candidate categories",
            read(debug_path, String),
        )
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_FILE[] = old_debug_file
    end
end

@testset "boundary final selection ignores training-only isotope groups" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2],
        candidate_index = UInt16[1, 2, 1, 2],
        candidate_category = ["fallback", "valley_combination", "fallback", "valley_combination"],
        is_fallback = Bool[true, false, true, false],
        quant_trace_selected = Bool[true, true, false, false],
        ms_file_idx = UInt16[1, 1, 1, 1],
        shape_model_valid = Bool[true, true, true, true],
        shape_fit_loss = Float32[0.3, 0.1, 0.2, 0.0],
        shape_edge_fraction = Float32[0, 0, 0, 0],
        shape_overshoot_penalty = Float32[0, 0, 0, 0],
        shape_sigma_irt = Float32[1, 1, 1, 1],
        shape_tau_irt = Float32[0, 0, 0, 0],
    )

    selected = Pioneer.select_boundary_candidate_rows_by_shape(candidates)

    @test nrow(selected) == 1
    @test selected.boundary_group_id == UInt64[1]
    @test all(selected.quant_trace_selected)
end

@testset "targeted boundary candidate debug log defaults to file one" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[10, 10, 20, 30],
        ms_file_idx = UInt16[1, 1, 2, 3],
        precursor_idx = UInt32[370714, 370714, 370714, 909488],
        isotope_key = [(0, 0), (0, 0), (0, 0), (0, 0)],
        quant_trace_selected = Bool[true, true, true, false],
        candidate_index = UInt16[1, 2, 1, 1],
        candidate_category = ["fallback", "valley_combination", "fallback", "fallback"],
        candidate_start_idx = UInt16[3, 3, 4, 5],
        candidate_stop_idx = UInt16[8, 11, 9, 12],
        candidate_start_scan = UInt32[30, 30, 40, 50],
        candidate_stop_scan = UInt32[80, 110, 90, 120],
        peak_area = Float32[100, 150, 200, 250],
        points_integrated = UInt32[5, 8, 4, 7],
        is_fallback = Bool[true, false, true, true],
        boundary_model_score = Float32[0.1, 0.9, 0.2, 0.7],
        boundary_shape_score = Float32[0.9, 0.1, 0.8, 0.3],
        shape_fit_loss = Float32[0.8, 0.05, 0.7, 0.2],
        shape_prior_penalty = Float32[0.1, 0.05, 0.1, 0.1],
        shape_sigma_irt = Float32[0.4, 0.5, 0.4, 0.6],
        shape_tau_irt = Float32[0.0, 0.1, 0.0, 0.2],
        shape_edge_fraction = Float32[0.2, 0.1, 0.2, 0.3],
        shape_overshoot_penalty = Float32[0.0, 0.0, 0.1, 0.0],
    )
    selected = candidates[[2, 3, 4], :]

    lines = Pioneer.boundary_candidate_debug_lines(
        candidates,
        selected,
    )

    @test length(lines) == 3
    @test occursin("precursor_idx=370714", first(lines))
    @test occursin("ms_file_idx=1", first(lines))
    @test occursin("rows=2", first(lines))
    @test any(occursin("candidate_index=2", line) && occursin("selected=true", line) for line in lines)
    @test !any(occursin("ms_file_idx=2", line) for line in lines)
    lines_909488 = Pioneer.boundary_candidate_debug_lines(
        candidates,
        selected;
        target_precursor_idx = UInt32(909488),
        target_ms_file_idx = UInt16(3),
    )
    @test occursin("precursor_idx=909488", first(lines_909488))
    @test occursin("ms_file_idx=3", first(lines_909488))
    @test any(occursin("selected=true", line) for line in lines_909488)

    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_debug_file = Pioneer.DEBUG_FILE[]
    old_targets = copy(Pioneer.DEBUG_BOUNDARY_CANDIDATE_TARGETS[])
    try
        Pioneer.DEBUG_BOUNDARY_CANDIDATE_TARGETS[] = Set([
            (UInt32(370714), UInt16(1)),
            (UInt32(909488), UInt16(3)),
        ])
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
        silent_path = tempname()
        open(silent_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            Pioneer.log_boundary_candidate_debug(
                candidates,
                selected,
            )
        end
        @test isempty(read(silent_path, String))

        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        debug_path = tempname()
        open(debug_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            redirect_stdout(devnull) do
                Pioneer.log_boundary_candidate_debug(
                    candidates,
                    selected,
                )
            end
        end
        debug_text = read(debug_path, String)
        @test occursin("Boundary candidate debug", debug_text)
        @test occursin("precursor_idx=909488 ms_file_idx=3", debug_text)
        @test occursin("quant_trace_selected=false", debug_text)
        @test length(collect(eachmatch(r"Boundary candidate debug", debug_text))) == 2
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_FILE[] = old_debug_file
        Pioneer.DEBUG_BOUNDARY_CANDIDATE_TARGETS[] = old_targets
    end
end

@testset "integrate_chrom exposes shape-model candidate rows" begin
    rt = Float32.(1:9)
    scan_idx = UInt32.(1:9)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 3, 9, 5, 2, 0.5, 0.1, 0]
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        4,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        boundary_candidate_data = candidate_data,
    )

    @test candidate_data[] !== nothing
    candidates = DataFrame(candidate_data[])
    @test nrow(candidates) > 1
    @test any(candidates.is_fallback)
    @test all(isfinite, candidates.peak_area)
    @test all(col in propertynames(candidates) for col in (
        :candidate_index,
        :candidate_start_idx,
        :candidate_stop_idx,
        :candidate_start_scan,
        :candidate_stop_scan,
        :candidate_category,
        :peak_area,
        :new_best_scan,
        :points_integrated,
        :is_fallback,
        :is_primary_apex,
        :shape_model_valid,
        :shape_fit_loss,
        :shape_sigma_irt,
        :shape_tau_irt,
        :shape_apex_irt,
        :shape_apex_height,
        :shape_edge_fraction,
        :shape_overshoot_penalty,
        :shape_deconvolution_area_fraction,
    ))
    @test any(candidates.shape_model_valid)
    @test :hardcoded_score ∉ propertynames(candidates)
    @test :candidate_width ∉ propertynames(candidates)
    @test :endpoint_height_fraction ∉ propertynames(candidates)
    @test :internal_dip_recovery_score ∉ propertynames(candidates)
    @test :mainsearch_left_bound_delta ∉ propertynames(candidates)
    @test :mainsearch_right_bound_delta ∉ propertynames(candidates)
    @test :fallback_start_delta ∉ propertynames(candidates)
    @test :fallback_stop_delta ∉ propertynames(candidates)
end

@testset "integrate_precursors can collect boundary candidate data" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[42, 42, 42, 42, 42, 42, 42],
        rt = Float32.(1:7),
        scan_idx = UInt32.(1:7),
        intensity = Float32[0, 1, 4, 8, 4, 1, 0],
        precursor_fraction_transmitted = fill(1.0f0, 7),
    )
    peak_area = zeros(Float32, 1)
    new_best_scan = zeros(UInt32, 1)
    points_integrated = zeros(UInt32, 1)
    boundary_candidate_data = Vector{Any}(nothing, 1)

    Pioneer.integrate_precursors(
        chromatograms,
        Pioneer.CombineTraces(0.25f0),
        0.25f0,
        UInt32[42],
        UInt32[4],
        peak_area,
        new_best_scan,
        points_integrated;
        λ = 0.0f0,
        boundary_candidate_data = boundary_candidate_data,
    )

    @test boundary_candidate_data[1] !== nothing
    @test nrow(DataFrame(boundary_candidate_data[1])) > 1
    @test peak_area[1] > 0.0f0
end

@testset "targeted debug plots use selected bounds with initial expected apex" begin
    chromatograms = DataFrame(
        precursor_idx = fill(UInt32(42), 9),
        rt = Float32.(1:9),
        scan_idx = UInt32.(10:10:90),
        intensity = Float32[0, 1, 10, 6, 1, 6, 12, 6, 0],
        precursor_fraction_transmitted = fill(1.0f0, 9),
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[42],
        scan_idx = UInt32[30],
        new_best_scan = UInt32[70],
    )
    selected = DataFrame(
        ms_file_idx = UInt16[1],
        precursor_idx = UInt32[42],
        candidate_start_scan = UInt32[20],
        candidate_stop_scan = UInt32[80],
    )
    debug_data_by_precursor = Dict{UInt32, Any}()

    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_targets = copy(Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[])
    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = Set(UInt32[42])
        tmp = mktempdir()

        Pioneer.debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            0.25f0,
            0.0f0,
            tmp,
            1,
            "run.raw";
            selected_boundary_candidates = selected,
            debug_plot_data_by_precursor = debug_data_by_precursor,
        )

        @test haskey(debug_data_by_precursor, UInt32(42))
        @test debug_data_by_precursor[UInt32(42)].scan_range == 2:8
        @test debug_data_by_precursor[UInt32(42)].apex_scan_idx == UInt32(30)
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_CHROM_TARGET_PRECURSOR_IDXS[] = old_targets
    end
end

@testset "targeted chromatogram integration debug plots" begin
    chromatograms = DataFrame(
        precursor_idx = UInt32[909488, 909488, 909488, 591443, 591443, 591443, 123],
        rt = Float32[1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0],
        scan_idx = UInt32[10, 20, 30, 10, 20, 30, 10],
        intensity = Float32[0.5, 3.0, 0.7, 0.4, 2.5, 0.6, 1.0],
        precursor_fraction_transmitted = Float32[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[909488, 591443, 123],
        scan_idx = UInt32[20, 20, 10],
    )

    _with_debug_chromatogram_plot_state() do
        tmp = mktempdir()
        Pioneer.debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            0.25f0,
            1.0f-6,
            tmp,
            2,
            "run A.raw",
        )

        plot_dir = joinpath(tmp, "qc_plots", "chromatogram_integration_debug")
        files = filter(path -> endswith(path, ".png"), readdir(plot_dir; join = true))
        @test length(files) == 1
        @test occursin("file_2", basename(only(files)))
        @test occursin("run_A_raw", basename(only(files)))
        @test occursin("precursor_909488", basename(only(files)))
        @test !occursin("row_", basename(only(files)))
        @test filesize(only(files)) > 0
    end

    _with_debug_chromatogram_plot_state() do
        tmp = mktempdir()
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
        Pioneer.debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            0.25f0,
            1.0f-6,
            tmp,
            2,
            "run A.raw",
        )

        plot_dir = joinpath(tmp, "qc_plots", "chromatogram_integration_debug")
        @test !isdir(plot_dir)
    end
end
