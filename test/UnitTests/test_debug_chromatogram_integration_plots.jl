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

@testset "mainsearch q-value intervals are captured before best scan selection" begin
    psms = DataFrame(
        precursor_idx = UInt32[1, 1, 99, 1, 2],
        scan_idx = UInt32[10, 20, 25, 30, 40],
        rt = Float32[1.0, 2.0, 2.5, 3.0, 4.0],
        target = Bool[true, true, false, true, true],
        lgbm_score = Float32[0.99, 0.98, 0.97, 0.96, 0.95],
    )

    Pioneer.add_mainsearch_integration_interval_columns!(psms, :lgbm_score)

    @test psms.mainsearch_1pct_start_scan == UInt32[10, 10, 0, 10, 0]
    @test psms.mainsearch_1pct_stop_scan == UInt32[20, 20, 0, 20, 0]
    @test isequal(psms.mainsearch_1pct_start_rt, Float32[1.0, 1.0, NaN, 1.0, NaN])
    @test isequal(psms.mainsearch_1pct_stop_rt, Float32[2.0, 2.0, NaN, 2.0, NaN])
end

@testset "mainsearch q-value evidence is retained as learned boundary features" begin
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
    @test :mainsearch_exclusion_penalty in propertynames(candidates)
    @test :mainsearch_left_bound_delta in propertynames(candidates)
    @test :mainsearch_right_bound_delta in propertynames(candidates)
    @test :mainsearch_distance_penalty ∉ propertynames(candidates)
    @test :mainsearch_center_distance_penalty ∉ propertynames(candidates)
    @test :mainsearch_width_penalty ∉ propertynames(candidates)
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

@testset "wide mainsearch q-value evidence does not override fallback bounds" begin
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

@testset "shoulder inside mainsearch evidence remains a feature, not a candidate source" begin
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
    @test :mainsearch_left_bound_delta in propertynames(candidates)
    @test :mainsearch_right_bound_delta in propertynames(candidates)
    @test all(candidates.candidate_width .>= 0)
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

    deltas_by_range = Dict(
        (Int(row.candidate_start_idx), Int(row.candidate_stop_idx)) => (
            row.mainsearch_left_bound_delta,
            row.mainsearch_right_bound_delta,
        )
        for row in eachrow(candidates)
    )
    @test deltas_by_range[(1, 11)] == (0.0f0, 0.0f0)
    @test deltas_by_range[(2, 10)] == (1.0f0, -1.0f0)
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
    fallback = only(candidates[candidates.is_fallback, :])
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
    valley_candidates = candidates[candidates.candidate_category .== "valley_combination", :]
    allowed_left = Set(UInt16[5, 7, 9])
    allowed_right = Set(UInt16[11, 13, 15])

    @test nrow(valley_candidates) > 0
    @test all(in.(valley_candidates.candidate_start_idx, Ref(allowed_left)))
    @test all(in.(valley_candidates.candidate_stop_idx, Ref(allowed_right)))
    @test length(unique(valley_candidates.candidate_start_idx)) <= 3
    @test length(unique(valley_candidates.candidate_stop_idx)) <= 3
end

@testset "mainsearch q-value intervals merge across trace rows" begin
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

@testset "boundary candidate training labels follow cross-run protein ratios" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 2, 3, 4, 4],
        precursor_idx = UInt32[10, 10, 20, 20, 20],
        protein_key = ["P1", "P1", "P1", "P1", "P1"],
        ms_file_idx = UInt16[1, 2, 1, 2, 2],
        candidate_index = UInt16[1, 1, 1, 1, 2],
        peak_area = Float32[100, 200, 50, 400, 100],
        is_fallback = Bool[true, true, true, true, false],
    )

    Pioneer.label_boundary_candidate_targets!(candidates; min_crossrun_refs = 1)

    group4 = candidates[candidates.boundary_group_id .== UInt64(4), :]
    @test group4.boundary_candidate_target == Bool[false, true]
    @test group4.boundary_consistency_loss[2] < group4.boundary_consistency_loss[1]
end

@testset "boundary candidate training labels use isotope traces from combined mode" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 2, 2],
        precursor_idx = UInt32[10, 10, 10],
        protein_key = ["P1", "P1", "P1"],
        isotope_key = Tuple{Int8, Int8}[(1, 1), (2, 2), (2, 2)],
        quant_trace_selected = Bool[false, false, false],
        ms_file_idx = UInt16[1, 1, 1],
        candidate_index = UInt16[1, 1, 2],
        peak_area = Float32[100, 400, 100],
        is_fallback = Bool[true, true, false],
    )

    Pioneer.label_boundary_candidate_targets!(
        candidates;
        min_crossrun_refs = 100,
        min_isotope_refs = 1,
    )

    group2 = candidates[candidates.boundary_group_id .== UInt64(2), :]
    @test group2.boundary_candidate_target == Bool[false, true]
    @test group2.boundary_isotope_loss[2] < group2.boundary_isotope_loss[1]
end

@testset "cross-run boundary loss uses only quant-selected traces" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 3, 4, 5],
        precursor_idx = UInt32[10, 10, 10, 20, 20, 10],
        protein_key = ["P1", "P1", "P1", "P1", "P1", "P1"],
        isotope_key = Tuple{Int8, Int8}[
            (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 1),
        ],
        quant_trace_selected = Bool[true, true, true, true, true, false],
        ms_file_idx = UInt16[1, 1, 2, 1, 2, 2],
        candidate_index = UInt16[1, 2, 1, 1, 1, 1],
        peak_area = Float32[100, 50, 100, 50, 100, 6400],
        is_fallback = Bool[true, false, true, true, true, true],
    )

    Pioneer.label_boundary_candidate_targets!(
        candidates;
        min_crossrun_refs = 1,
        min_isotope_refs = 100,
    )

    group1 = candidates[candidates.boundary_group_id .== UInt64(1), :]
    @test group1.boundary_candidate_target == Bool[false, true]
    @test group1.boundary_consistency_loss[2] ≈ 0.0f0
end

@testset "boundary labels skip groups without self-supervised objective" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 1],
        precursor_idx = UInt32[10, 10, 10],
        protein_key = ["P1", "P1", "P1"],
        ms_file_idx = UInt16[1, 1, 1],
        candidate_index = UInt16[1, 2, 3],
        peak_area = Float32[100, 200, 300],
        is_fallback = Bool[true, false, false],
    )

    Pioneer.label_boundary_candidate_targets!(
        candidates;
        min_crossrun_refs = 100,
        min_isotope_refs = 100,
    )

    @test all(isinf, candidates.boundary_training_loss)
    @test !any(candidates.boundary_candidate_target)
end

@testset "boundary candidate sampling keeps whole candidate groups" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2, 3, 3],
        candidate_index = UInt16[1, 2, 1, 2, 1, 2],
        peak_area = Float32[1, 2, 3, 4, 5, 6],
    )

    sampled = Pioneer.sample_boundary_candidate_groups(
        candidates;
        max_groups = 2,
        rng = Random.MersenneTwister(1776),
    )

    sampled_groups = unique(sampled.boundary_group_id)
    @test length(sampled_groups) == 2
    @test all(count(==(group_id), sampled.boundary_group_id) == 2 for group_id in sampled_groups)
end

@testset "boundary ranking labels preserve within-group loss ordering" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[2, 1, 1, 2, 1],
        candidate_index = UInt16[1, 1, 2, 2, 3],
        boundary_training_loss = Float32[0.5, 0.0, 4.0, 0.5, 1.0],
        candidate_width = Float32[3, 4, 5, 6, 7],
    )

    frame, relevance, group_sizes = Pioneer.boundary_ranking_training_data(
        candidates,
        [:candidate_width],
    )

    @test group_sizes == [3, 2]
    @test relevance == Int32[2, 0, 1, 1, 1]
    @test frame.candidate_width == Float32[4, 5, 7, 3, 6]
end

@testset "boundary ranking skips groups without finite objective loss" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2],
        candidate_index = UInt16[1, 2, 1, 2],
        boundary_training_loss = Float32[Inf, Inf, 0, 1],
        candidate_width = Float32[1, 2, 3, 4],
    )

    frame, relevance, group_sizes = Pioneer.boundary_ranking_training_data(
        candidates,
        [:candidate_width],
    )

    @test group_sizes == [2]
    @test relevance == Int32[1, 0]
    @test frame.candidate_width == Float32[3, 4]
end

@testset "boundary ranking labels fit LightGBM default label mapping" begin
    candidates = DataFrame(
        boundary_group_id = fill(UInt64(1), 53),
        candidate_index = UInt16.(1:53),
        boundary_training_loss = Float32.(0:52),
        candidate_width = Float32.(1:53),
    )

    _, relevance, group_sizes = Pioneer.boundary_ranking_training_data(
        candidates,
        [:candidate_width],
    )

    @test group_sizes == [53]
    @test maximum(relevance) <= 30
    @test relevance[1] == 30
    @test relevance[end] == 0
    @test all(diff(relevance) .<= 0)
end

@testset "boundary model trains a LightGBM ranker" begin
    rows = NamedTuple[]
    for group_id in UInt64(1):UInt64(30), candidate_index in UInt16(1):UInt16(3)
        push!(rows, (;
            boundary_group_id = group_id,
            candidate_index = candidate_index,
            boundary_training_loss = Float32(3 - candidate_index),
            boundary_candidate_target = candidate_index == UInt16(3),
            candidate_width = Float32(candidate_index),
        ))
    end
    candidates = DataFrame(rows)

    model = Pioneer.train_boundary_candidate_model(
        candidates;
        features = [:candidate_width],
        max_groups = 30,
        min_positive = 1,
        min_negative = 1,
        rng = Random.MersenneTwister(1776),
    )

    @test model !== nothing
    @test model.booster isa LightGBM.LGBMRanking
    @test model.features == [:candidate_width]
end

@testset "boundary model trains when candidate groups exceed default label gain length" begin
    rows = NamedTuple[]
    for group_id in UInt64(1):UInt64(30), candidate_index in UInt16(1):UInt16(53)
        push!(rows, (;
            boundary_group_id = group_id,
            candidate_index = candidate_index,
            boundary_training_loss = Float32(candidate_index - UInt16(1)),
            boundary_candidate_target = candidate_index == UInt16(1),
            candidate_width = Float32(candidate_index),
        ))
    end
    candidates = DataFrame(rows)

    model = Pioneer.train_boundary_candidate_model(
        candidates;
        features = [:candidate_width],
        max_groups = 30,
        min_positive = 1,
        min_negative = 1,
        rng = Random.MersenneTwister(1776),
    )

    @test model !== nothing
    @test model.booster isa LightGBM.LGBMRanking
end

@testset "boundary cross-fold helpers exclude held-out cv folds" begin
    candidates = DataFrame(
        boundary_group_id = UInt64[1, 1, 2, 2, 3, 3],
        candidate_index = UInt16[1, 2, 1, 2, 1, 2],
        cv_fold = UInt8[0, 0, 1, 1, 2, 2],
        boundary_training_loss = Float32[1, 0, 1, 0, 1, 0],
        boundary_candidate_target = Bool[false, true, false, true, false, true],
        candidate_width = Float32[1, 2, 1, 2, 1, 2],
    )

    fold1_training = Pioneer.boundary_training_candidates_for_cv_fold(candidates, UInt8(1))

    @test all(fold1_training.cv_fold .!= UInt8(1))
    @test sort(unique(fold1_training.boundary_group_id)) == UInt64[1, 3]
end

@testset "boundary model trains and selects candidates by held-out cv fold" begin
    rows = NamedTuple[]
    for cv_fold in UInt8[0, 1], group_idx in UInt64(1):UInt64(30), candidate_index in UInt16(1):UInt16(3)
        push!(rows, (;
            boundary_group_id = UInt64(cv_fold) << 32 | group_idx,
            cv_fold = cv_fold,
            candidate_index = candidate_index,
            is_fallback = candidate_index == UInt16(1),
            quant_trace_selected = true,
            boundary_training_loss = Float32(3 - candidate_index),
            boundary_candidate_target = candidate_index == UInt16(3),
            candidate_width = Float32(candidate_index),
        ))
    end
    candidates = DataFrame(rows)

    models = Pioneer.train_boundary_candidate_models_by_cv_fold(
        candidates;
        features = [:candidate_width],
        max_groups = 60,
        min_positive = 1,
        min_negative = 1,
        rng = Random.MersenneTwister(1776),
    )
    selected = Pioneer.select_boundary_candidate_rows_crossfold(candidates, models)

    @test sort(collect(keys(models))) == UInt8[0, 1]
    @test all(model !== nothing for model in values(models))
    @test nrow(selected) == 60
    @test all(selected.candidate_index .== UInt16(3))
end

@testset "boundary model feature importance log lines are sorted by gain" begin
    importances = [
        (:candidate_width, 3.2),
        (:peak_prominence_score, 15.8),
        (:endpoint_valley_score, 8.1),
    ]
    lines = Pioneer.boundary_model_feature_importance_lines(importances)

    @test first(lines) == "Learned chromatogram boundary LGBM feature gains (all 3):"
    @test occursin("peak_prominence_score", lines[2])
    @test occursin("endpoint_valley_score", lines[3])
    @test occursin("candidate_width", lines[4])

    old_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
    old_debug_file = Pioneer.DEBUG_FILE[]
    try
        Pioneer.DEBUG_CONSOLE_LEVEL[] = 0
        silent_path = tempname()
        open(silent_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            Pioneer.log_boundary_model_feature_importances(importances)
        end
        @test isempty(read(silent_path, String))

        Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
        debug_path = tempname()
        open(debug_path, "w") do io
            Pioneer.DEBUG_FILE[] = io
            redirect_stdout(devnull) do
                Pioneer.log_boundary_model_feature_importances(importances)
            end
        end
        @test occursin("Learned chromatogram boundary LGBM feature gains", read(debug_path, String))
    finally
        Pioneer.DEBUG_CONSOLE_LEVEL[] = old_level
        Pioneer.DEBUG_FILE[] = old_debug_file
    end
end

@testset "boundary shape features flag outside signal and recovered internal dips" begin
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 4, 10, 9, 6, 7, 8, 6, 3, 0]
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

    candidates = DataFrame(candidate_data[])
    narrow = only(candidates[
        (candidates.candidate_start_idx .== UInt16(2)) .&
        (candidates.candidate_stop_idx .== UInt16(6)),
        :,
    ])
    wide = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(11)),
        :,
    ])

    @test narrow.right_excluded_signal_fraction ≈ 0.6f0
    @test narrow.right_boundary_recovery_fraction ≈ 0.2f0
    @test narrow.right_outside_peak_fraction ≈ 0.8f0
    @test narrow.left_boundary_recovery_fraction == 0.0f0
    @test wide.right_excluded_signal_fraction == 0.0f0
    @test wide.right_boundary_recovery_fraction == 0.0f0
    @test wide.right_outside_peak_fraction == 0.0f0
    @test narrow.internal_dip_recovery_score == 0.0f0
    @test wide.internal_dip_recovery_score > 0.0f0
    @test wide.internal_dip_recovery_score > narrow.internal_dip_recovery_score
end

@testset "edge valley ratio features compare against included internal valleys" begin
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[5, 0, 4, 1, 3, 10, 3, 0, 4, 2, 5]
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
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    nearest = only(candidates[
        (candidates.candidate_start_idx .== UInt16(4)) .&
        (candidates.candidate_stop_idx .== UInt16(8)),
        :,
    ])
    lower_outer_left = only(candidates[
        (candidates.candidate_start_idx .== UInt16(2)) .&
        (candidates.candidate_stop_idx .== UInt16(8)),
        :,
    ])
    higher_outer_right = only(candidates[
        (candidates.candidate_start_idx .== UInt16(4)) .&
        (candidates.candidate_stop_idx .== UInt16(10)),
        :,
    ])

    @test nearest.left_edge_valley_log2_ratio == 0.0f0
    @test nearest.right_edge_valley_log2_ratio == 0.0f0
    @test lower_outer_left.left_edge_valley_log2_ratio < 0.0f0
    @test lower_outer_left.right_edge_valley_log2_ratio == 0.0f0
    @test higher_outer_right.left_edge_valley_log2_ratio == 0.0f0
    @test higher_outer_right.right_edge_valley_log2_ratio > 0.0f0
end

@testset "included non-apex max features describe competing signal near or far from apex" begin
    rt = Float32.(1:13)
    scan_idx = UInt32.(1:13)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[5, 3, 4, 1, 5, 10, 5, 0, 8, 20, 8, 0, 5]
    rt_to_irt = Pioneer.LinearRtConversionModel(2.0f0, 10.0f0)
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
        rt_to_irt_model = rt_to_irt,
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    expected_peak_only = only(candidates[
        (candidates.candidate_start_idx .== UInt16(4)) .&
        (candidates.candidate_stop_idx .== UInt16(8)),
        :,
    ])
    competing_peak_included = only(candidates[
        (candidates.candidate_start_idx .== UInt16(4)) .&
        (candidates.candidate_stop_idx .== UInt16(12)),
        :,
    ])

    @test expected_peak_only.included_nonapex_max_log2_ratio ≈ -1.0f0
    @test expected_peak_only.included_nonapex_max_irt_distance ≈ 2.0f0
    @test competing_peak_included.included_nonapex_max_log2_ratio ≈ 1.0f0
    @test competing_peak_included.included_nonapex_max_irt_distance ≈ 8.0f0
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
    @test first(lines) == "Learned chromatogram boundary selected candidate categories (total 3):"
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
            "Learned chromatogram boundary selected candidate categories",
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
    )

    selected = Pioneer.select_boundary_candidate_rows(candidates, nothing)

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
        boundary_training_loss = Float32[1.5, 0.2, 0.8, 0.3],
        boundary_crossrun_loss = Float32[1.0, 0.1, 0.5, 0.4],
        boundary_isotope_loss = Float32[Inf32, Inf32, Inf32, Inf32],
        candidate_width = Float32[6, 9, 6, 7],
        endpoint_valley_score = Float32[0.9, 0.4, 0.8, 0.7],
        endpoint_height_fraction = Float32[0.2, 0.3, 0.2, 0.25],
        secondary_peak_penalty = Float32[0.0, 0.3, 0.0, 0.2],
        internal_dip_recovery_score = Float32[0.0, 0.5, 0.1, 0.4],
        right_excluded_signal_fraction = Float32[0.4, 0.0, 0.0, 0.2],
        right_boundary_recovery_fraction = Float32[0.3, 0.0, 0.0, 0.1],
        right_outside_peak_fraction = Float32[0.6, 0.0, 0.0, 0.3],
        mainsearch_left_bound_delta = Float32[1, 1, 0, -1],
        mainsearch_right_bound_delta = Float32[-2, 1, 0, 2],
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

@testset "integrate_chrom exposes candidate features for boundary learner" begin
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
    @test all(candidates.candidate_width .>= 0)
    @test all(isfinite, candidates.peak_area)
    @test :candidate_category in propertynames(candidates)
    @test :hardcoded_score ∉ propertynames(candidates)
    @test :mainsearch_distance_penalty ∉ propertynames(candidates)
    @test :mainsearch_center_distance_penalty ∉ propertynames(candidates)
    @test :mainsearch_width_penalty ∉ propertynames(candidates)
    @test :fallback_start_delta ∉ propertynames(candidates)
    @test :fallback_stop_delta ∉ propertynames(candidates)
    @test all(Symbol(feature) in propertynames(candidates) for feature in Pioneer.BOUNDARY_CANDIDATE_FEATURES)
    @test :mainsearch_distance_penalty ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :mainsearch_center_distance_penalty ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :mainsearch_width_penalty ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :fallback_start_delta ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :fallback_stop_delta ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :endpoint_valley_score ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :secondary_peak_penalty ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :mainsearch_exclusion_penalty ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :mainsearch_left_bound_delta in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :mainsearch_right_bound_delta in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :smoothed_apex_weight ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :log2_smoothed_apex_weight in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :left_excluded_signal_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :right_excluded_signal_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :left_boundary_recovery_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :right_boundary_recovery_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :left_outside_peak_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :right_outside_peak_fraction in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :internal_dip_recovery_score in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :left_edge_valley_log2_ratio in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :right_edge_valley_log2_ratio in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :included_nonapex_max_log2_ratio in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :included_nonapex_max_irt_distance in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :irt_asymmetry_delta in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :baseline_disconnected_signal_fraction ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :baseline_largest_nonapex_lobe_log2_ratio ∉ Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :baseline_internal_trough_score in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :baseline_internal_trough_log2_ratio in Pioneer.BOUNDARY_CANDIDATE_FEATURES
    @test :baseline_disconnected_signal_fraction in propertynames(candidates)
    @test :baseline_largest_nonapex_lobe_log2_ratio in propertynames(candidates)
    @test :smoothed_apex_weight ∉ propertynames(candidates)
    @test all(candidates.log2_smoothed_apex_weight .≈ log2(9.0f0))
end

@testset "boundary candidate width uses iRT span" begin
    rt = Float32[0.0, 0.25, 1.0, 2.5, 5.0, 8.0]
    scan_idx = UInt32.(1:length(rt))
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 8, 20, 5, 0]
    rt_to_irt = Pioneer.LinearRtConversionModel(2.0f0, 10.0f0)
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
        rt_to_irt_model = rt_to_irt,
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    @test nrow(candidates) > 0
    for candidate in eachrow(candidates)
        expected_width = abs(
            Float32(rt_to_irt(rt[Int(candidate.candidate_stop_idx)])) -
            Float32(rt_to_irt(rt[Int(candidate.candidate_start_idx)])),
        )
        scan_count_width = Float32(
            Int(candidate.candidate_stop_idx) - Int(candidate.candidate_start_idx) + 1,
        )
        @test candidate.candidate_width ≈ expected_width
        @test candidate.candidate_width != scan_count_width
    end
end

@testset "boundary distance features use iRT geometry" begin
    rt = Float32[0.0, 0.2, 0.55, 1.4, 2.0, 3.8, 4.1, 6.5]
    scan_idx = UInt32.(10:10:(10 * length(rt)))
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 4, 10, 6, 3, 1, 0]
    rt_to_irt = Pioneer.LinearRtConversionModel(2.0f0, 10.0f0)
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)
    apex_idx = 4
    evidence_range = 2:7

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        apex_idx,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        rt_to_irt_model = rt_to_irt,
        mainsearch_1pct_start_scan = scan_idx[first(evidence_range)],
        mainsearch_1pct_stop_scan = scan_idx[last(evidence_range)],
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    evidence_start_irt = Float32(rt_to_irt(rt[first(evidence_range)]))
    evidence_stop_irt = Float32(rt_to_irt(rt[last(evidence_range)]))
    evidence_width = evidence_stop_irt - evidence_start_irt
    apex_irt = Float32(rt_to_irt(rt[apex_idx]))
    function expected_area(start_idx::Int, stop_idx::Int)
        start_idx >= stop_idx && return 0.0f0
        area = 0.0f0
        for i in start_idx:(stop_idx - 1)
            left_val = max(0.0f0, intensity[i])
            right_val = max(0.0f0, intensity[i + 1])
            span = abs(Float32(rt_to_irt(rt[i + 1])) - Float32(rt_to_irt(rt[i])))
            area += 0.5f0 * (left_val + right_val) * span
        end
        return area
    end
    for candidate in eachrow(candidates)
        start_idx = Int(candidate.candidate_start_idx)
        stop_idx = Int(candidate.candidate_stop_idx)
        start_irt = Float32(rt_to_irt(rt[start_idx]))
        stop_irt = Float32(rt_to_irt(rt[stop_idx]))

        overlap_start = max(start_idx, first(evidence_range))
        overlap_stop = min(stop_idx, last(evidence_range))
        overlap_width = overlap_start <= overlap_stop ?
            max(
                Float32(rt_to_irt(rt[overlap_stop])) -
                Float32(rt_to_irt(rt[overlap_start])),
                0.0f0,
            ) :
            0.0f0
        expected_exclusion = 1.0f0 - overlap_width / evidence_width
        expected_left_delta = start_irt - evidence_start_irt
        expected_right_delta = stop_irt - evidence_stop_irt
        expected_irt_asymmetry = (stop_irt - apex_irt) - (apex_irt - start_irt)
        left_area = expected_area(start_idx, apex_idx)
        right_area = expected_area(apex_idx, stop_idx)
        expected_asymmetry = log2(max(right_area, eps(Float32)) / max(left_area, eps(Float32)))

        @test candidate.mainsearch_exclusion_penalty ≈ expected_exclusion
        @test candidate.mainsearch_left_bound_delta ≈ expected_left_delta
        @test candidate.mainsearch_right_bound_delta ≈ expected_right_delta
        @test candidate.asymmetry_penalty ≈ expected_asymmetry
        @test candidate.irt_asymmetry_delta ≈ expected_irt_asymmetry
    end
end

@testset "boundary asymmetry preserves left/right area direction" begin
    rt = Float32.(1:11)
    scan_idx = UInt32.(1:11)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 8, 7, 10, 2, 1, 0, 0, 0, 0, 0]
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

    candidates = DataFrame(candidate_data[])
    left_heavy = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(7)),
        :,
    ])

    @test left_heavy.asymmetry_penalty ≈ log2(8.0f0 / 20.0f0)
end

@testset "boundary shape features use WH-smoothed signal" begin
    rt = Float32.(1:7)
    scan_idx = UInt32.(1:7)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[5, 6, 8, 15, 8, 6, 5]
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

    candidates = DataFrame(candidate_data[])
    wide = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(7)),
        :,
    ])

    @test wide.endpoint_height_fraction ≈ (5.0f0 + 5.0f0) / (2.0f0 * 15.0f0)
    @test wide.peak_prominence_score ≈ ((15.0f0 - 5.0f0) / 15.0f0) / sqrt(wide.candidate_width)
    @test wide.included_nonapex_max_log2_ratio ≈ log2(8.0f0 / 15.0f0)
    @test wide.included_nonapex_max_irt_distance ≈ 1.0f0
end

@testset "baseline lobe features flag disconnected non-apex residual peaks" begin
    rt = Float32.(1:10)
    scan_idx = UInt32.(1:10)
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[20, 22, 30, 50, 35, 7, 25, 35, 20, 0]
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

    candidates = DataFrame(candidate_data[])
    connected = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(6)),
        :,
    ])
    disconnected = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(10)),
        :,
    ])

    @test connected.baseline_disconnected_signal_fraction == 0.0f0
    @test connected.baseline_largest_nonapex_lobe_log2_ratio == -10.0f0
    @test connected.baseline_internal_trough_score == 0.0f0
    @test connected.baseline_internal_trough_log2_ratio == 0.0f0
    @test disconnected.baseline_disconnected_signal_fraction > 0.2f0
    @test disconnected.baseline_largest_nonapex_lobe_log2_ratio > -2.0f0
    @test disconnected.baseline_internal_trough_score > 0.2f0
    @test disconnected.baseline_internal_trough_log2_ratio < -2.0f0
end

@testset "internal dip recovery uses iRT distance" begin
    rt = Float32[0.0, 0.2, 0.6, 1.0, 1.5, 2.4, 2.55, 4.2, 4.4, 4.6, 6.0]
    scan_idx = UInt32.(1:length(rt))
    fraction = fill(1.0f0, length(rt))
    intensity = Float32[0, 1, 4, 10, 9, 6, 7, 8, 6, 3, 0]
    rt_to_irt = Pioneer.LinearRtConversionModel(2.0f0, 10.0f0)
    ws = Pioneer.WHWorkspace(length(rt))
    state = Pioneer.Chromatogram(zeros(Float32, length(rt)), zeros(Float32, length(rt)), 0)
    candidate_data = Ref{Any}(nothing)
    apex_idx = 4

    Pioneer.integrate_chrom(
        rt,
        scan_idx,
        intensity,
        fraction,
        apex_idx,
        ws,
        state,
        1.0f0,
        0.0f0;
        min_fraction_transmitted = 0.25f0,
        rt_to_irt_model = rt_to_irt,
        boundary_candidate_data = candidate_data,
    )

    candidates = DataFrame(candidate_data[])
    wide = only(candidates[
        (candidates.candidate_start_idx .== UInt16(1)) .&
        (candidates.candidate_stop_idx .== UInt16(11)),
        :,
    ])

    apex_height = intensity[apex_idx]
    expected = 0.0f0
    min_seen = apex_height
    min_idx = apex_idx
    for i in (apex_idx - 1):-1:Int(wide.candidate_start_idx)
        val = max(0.0f0, intensity[i])
        if val < min_seen
            min_seen = val
            min_idx = i
        else
            distance = abs(Float32(rt_to_irt(rt[i])) - Float32(rt_to_irt(rt[min_idx])))
            expected = max(expected, (val - min_seen) / (apex_height * sqrt(max(distance, eps(Float32)))))
        end
    end
    min_seen = apex_height
    min_idx = apex_idx
    for i in (apex_idx + 1):Int(wide.candidate_stop_idx)
        val = max(0.0f0, intensity[i])
        if val < min_seen
            min_seen = val
            min_idx = i
        else
            distance = abs(Float32(rt_to_irt(rt[i])) - Float32(rt_to_irt(rt[min_idx])))
            expected = max(expected, (val - min_seen) / (apex_height * sqrt(max(distance, eps(Float32)))))
        end
    end

    @test wide.internal_dip_recovery_score ≈ expected
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

@testset "targeted debug plots use learned selected candidate bounds" begin
    chromatograms = DataFrame(
        precursor_idx = fill(UInt32(42), 9),
        rt = Float32.(1:9),
        scan_idx = UInt32.(10:10:90),
        intensity = Float32[0, 1, 3, 6, 10, 6, 3, 1, 0],
        precursor_fraction_transmitted = fill(1.0f0, 9),
    )
    passing_psms = DataFrame(
        precursor_idx = UInt32[42],
        scan_idx = UInt32[50],
    )
    selected = DataFrame(
        ms_file_idx = UInt16[1],
        precursor_idx = UInt32[42],
        candidate_start_scan = UInt32[30],
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
        @test debug_data_by_precursor[UInt32(42)].scan_range == 3:8
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
