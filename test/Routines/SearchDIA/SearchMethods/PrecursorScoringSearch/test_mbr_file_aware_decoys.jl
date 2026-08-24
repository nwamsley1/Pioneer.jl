using DataFrames
using Test
using Pioneer

function _test_integrated_mbr_donor(
    precursor_idx::UInt32,
    probability::Float32,
    file_idx::UInt32;
    irt::Float32 = 10.0f0,
)
    spectrum = (
        1.0f0, 0.0f0, 0.0f0, 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    )
    return Pioneer._MBRDonorEntry(
        probability,
        precursor_idx,
        100.0f0,
        -1.0f0,
        0.0f0,
        irt,
        5.0f0,
        spectrum,
        UInt8(0x03),
        UInt16(2),
        file_idx,
    )
end

@testset "post-integration MBR donor selection prefers similar runs" begin
    atlas = Pioneer.build_run_similarity(Dict(
        UInt32(1) => UInt32[1, 2],
        UInt32(2) => UInt32[1],
        UInt32(3) => UInt32[2],
    ))
    donor_dict = Dict(
        UInt32(10) => Pioneer._MBRDonorEntry[
            _test_integrated_mbr_donor(UInt32(10), 0.99f0, UInt32(2)),
            _test_integrated_mbr_donor(UInt32(10), 0.90f0, UInt32(1)),
        ],
    )

    donor = Pioneer._mbr_select_donor(
        donor_dict,
        UInt32(10),
        UInt32(3),
        atlas,
    )
    @test donor !== nothing
    @test donor.ms_file_idx == UInt32(1)
end

@testset "MBR counterfactuals match charge/length and prefer the real donor run" begin
    charge = fill(UInt8(2), 5)
    sequence_length = fill(UInt8(9), 5)
    irt = Float32[10.0, 10.01, 10.02, 10.03, 10.04]
    same_file_pool = (
        pids = UInt32[1, 2, 3],
        irts = irt[1:3],
    )
    global_pool = (
        pids = UInt32[1, 2, 3, 4, 5],
        irts = irt,
    )
    pools = Pioneer._MBRPartnerPools(
        fill(UInt8(0), 5),
        fill(UInt8(1), 5),
        charge,
        sequence_length,
        irt,
        Dict{Tuple{Int, Int, Int, Int}, Pioneer._MBRIrtPool}(),
        Dict{Tuple{Int, Int, Int}, Pioneer._MBRIrtPool}(),
        Dict((2, 9) => global_pool),
        Dict((UInt32(2), 2, 9) => same_file_pool),
    )
    eligibility = Pioneer._MBRCounterfactualEligibility(
        BitSet(1:5),
        Dict(UInt32(1) => BitSet()),
    )
    donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
        UInt32(1) => [
            _test_integrated_mbr_donor(UInt32(1), 0.95f0, UInt32(2)),
        ],
        UInt32(2) => [
            _test_integrated_mbr_donor(UInt32(2), 0.90f0, UInt32(2)),
        ],
        UInt32(3) => [
            _test_integrated_mbr_donor(UInt32(3), 0.89f0, UInt32(2)),
        ],
        UInt32(4) => [
            _test_integrated_mbr_donor(UInt32(4), 0.88f0, UInt32(3)),
        ],
    )
    true_donor = only(donor_dict[UInt32(1)])
    # `_mbr_false_donors` writes into a caller-owned buffer (it `fill!`s it with `nothing` itself and
    # returns it), so the buffer must be MBR_N_COUNTERFACTUALS long. fe8a00f76 introduced this
    # parameter to stop allocating a fresh vector per row, but did not update this call.
    # Explicitly `nothing`-filled rather than `undef`: the collector counts non-nothing entries to
    # find its write position, and only works on an `undef` buffer because `_mbr_false_donors`
    # happens to `fill!` first. Not worth depending on from a test.
    donors_buf = fill!(
        Vector{Union{Nothing, Pioneer._MBRDonorEntry}}(undef, Pioneer.MBR_N_COUNTERFACTUALS),
        nothing,
    )
    false_donors = Pioneer._mbr_false_donors(
        donor_dict,
        pools,
        eligibility,
        UInt32(1),
        UInt32(1),
        10.0f0,
        true_donor,
        nothing,
        donors_buf,
    )

    @test length(false_donors) == Pioneer.MBR_N_COUNTERFACTUALS
    @test false_donors[1].precursor_idx == UInt32(2)
    @test false_donors[2].precursor_idx == UInt32(3)
    @test false_donors[1].ms_file_idx == true_donor.ms_file_idx
    @test false_donors[2].ms_file_idx == true_donor.ms_file_idx
    @test false_donors[3].precursor_idx == UInt32(4)
end

@testset "combined MBR controller spends the remaining total error budget" begin
    result = Pioneer._mbr_combined_error_recovery(
        Float32[0.9, 0.8],
        BitVector([true, false]),
        Float32[0.1, 0.1];
        base_targets = 99,
        base_decoys = 0,
        alpha = 0.01f0,
    )
    @test result.recovered == BitVector([true, true])
    @test result.mbr_targets == 1
    @test result.mbr_decoys == 1
    @test result.total_errors == 1
    @test result.total_targets == 100
    @test result.combined_error_rate == 0.01f0

    no_budget = Pioneer._mbr_combined_error_recovery(
        Float32[0.9],
        BitVector([true]),
        Float32[0.1];
        base_targets = 99,
        base_decoys = 1,
        alpha = 0.01f0,
    )
    @test !any(no_budget.recovered)
end

@testset "MBR training caps are deterministic and counterfactual-stratified" begin
    n_candidates = 12
    folds = fill(UInt8(1), n_candidates)
    positives = trues(n_candidates)
    receiver_decoys = falses(n_candidates)
    present = falses(n_candidates, Pioneer.MBR_N_COUNTERFACTUALS)
    present[:, 1] .= true
    present[1:8, 2] .= true
    present[1:4, 3] .= true

    sampled = Pioneer._mbr_training_rows(
        folds,
        positives,
        receiver_decoys,
        present,
        UInt8(1);
        max_positives = 5,
        max_negatives = 7,
    )
    sampled_again = Pioneer._mbr_training_rows(
        folds,
        positives,
        receiver_decoys,
        present,
        UInt8(1);
        max_positives = 5,
        max_negatives = 7,
    )

    @test sampled.available_positives == 12
    @test sampled.used_positives == 5
    @test sampled.available_negatives == 24
    @test sampled.used_negatives == 7
    @test sampled.used_negatives_by_counterfactual == Int[4, 2, 1]
    @test sampled.rows == sampled_again.rows
    @test count(sampled.labels) == 5
    @test count(!, sampled.labels) == 7
    @test length(unique(sampled.rows)) == length(sampled.rows)

    uncapped = Pioneer._mbr_training_rows(
        folds,
        positives,
        receiver_decoys,
        present,
        UInt8(1);
        max_positives = 100,
        max_negatives = 100,
    )
    expected_rows = Int[]
    expected_labels = Bool[]
    for candidate_idx in 1:n_candidates
        push!(expected_rows, candidate_idx)
        push!(expected_labels, true)
        for counterfactual_idx in 1:Pioneer.MBR_N_COUNTERFACTUALS
            present[candidate_idx, counterfactual_idx] || continue
            push!(
                expected_rows,
                counterfactual_idx * n_candidates + candidate_idx,
            )
            push!(expected_labels, false)
        end
    end
    @test uncapped.rows == expected_rows
    @test uncapped.labels == expected_labels
end

@testset "MBR receiver decoys remain explicit training negatives" begin
    folds = fill(UInt8(1), 6)
    positives = BitVector([true, true, false, false, false, false])
    receiver_decoys =
        BitVector([false, false, true, true, true, true])
    present = falses(6, Pioneer.MBR_N_COUNTERFACTUALS)
    present[1, 1] = true
    present[3, 1] = true

    training = Pioneer._mbr_training_rows(
        folds,
        positives,
        receiver_decoys,
        present,
        UInt8(1);
        max_positives = 100,
        max_negatives = 100,
    )

    @test training.available_positives == 2
    @test training.available_receiver_decoys == 4
    @test training.available_negatives_by_counterfactual == [2, 0, 0]
    @test training.rows == [1, 7, 2, 3, 9, 4, 5, 6]
    @test training.labels ==
        Bool[true, false, true, false, false, false, false, false]
end

@testset "MBR candidates do not require a counterfactual donor" begin
    frame = DataFrame(
        global_qval = Float32[0.005, 0.005, 0.02, 0.005],
        qval = Float32[0.2, 0.005, 0.2, 0.2],
        MBR_best_is_missing_true = Bool[false, false, false, true],
    )

    @test Pioneer._mbr_candidate_mask(frame, 0.01f0) ==
        BitVector([true, false, false, false])
end

@testset "MBR evaluation retains no-counterfactual candidates" begin
    scores = Float32[
        0.9, 0.8, 0.7,
        0.1, 0.6, 0.2,
        0.3, 0.4, 0.5,
        0.2, 0.1, 0.4,
    ]
    present = BitMatrix([
        false false false
        true  true  false
        true  false true
    ])
    top_labels = BitVector([true, false, false])
    top_mask = BitVector([true, true, false])

    eval_mask, labels = Pioneer._mbr_eval_mask_and_labels(
        scores,
        present,
        top_labels,
        top_mask,
    )

    @test eval_mask[1]
    @test labels[1]
    @test eval_mask[2]
    @test !labels[2]
    @test !eval_mask[3]
    @test count(eval_mask) == 4
end

@testset "MBR Hellinger contrasts rank real and counterfactual donors" begin
    base = "MBR_best_temporal_frag_hellinger"
    frame = DataFrame(
        Pioneer._mbr_true_feature(base) => Float32[0.10, -1.0],
        Pioneer._mbr_false_feature(base, 1) => Float32[0.20, 0.40],
        Pioneer._mbr_false_feature(base, 2) => Float32[0.05, 0.30],
        Pioneer._mbr_false_feature(base, 3) => Float32[-1.0, 0.20],
    )

    Pioneer._mbr_add_hellinger_contrasts!(frame)

    @test frame.MBR_best_temporal_frag_hellinger_rank_true ==
        Float32[2, -1]
    @test frame.MBR_best_temporal_frag_hellinger_rank_false2 ==
        Float32[1, 2]
    @test frame.MBR_best_temporal_frag_hellinger_rank_false3 ==
        Float32[-1, 1]
    @test frame.MBR_best_temporal_frag_hellinger_margin_true[1] ≈ -0.05f0
    @test frame.MBR_best_temporal_frag_hellinger_margin_false2[1] ≈
        0.05f0
end

@testset "MBR receiver run-cluster support excludes the receiver run" begin
    clusters = Pioneer._MBRReceiverRunClusters(
        Dict(UInt32(1) => UInt32(1), UInt32(2) => UInt32(1),
             UInt32(3) => UInt32(1)),
        UInt32[3],
        Dict(UInt32(10) => Dict(UInt32(1) => UInt32(2))),
        Dict(UInt32(1) => BitSet([10])),
    )

    support = Pioneer._mbr_receiver_cluster_features(
        clusters,
        UInt32(10),
        UInt32(1),
    )
    @test support.support_count == 1.0f0
    @test support.peer_count == 2.0f0
    @test support.support_fraction == 0.5f0
end

@testset "MBR donor feature tuple stays aligned with its sidecar schema" begin
    donor = _test_integrated_mbr_donor(
        UInt32(10),
        0.95f0,
        UInt32(2),
    )
    donor_dict = Dict(UInt32(10) => [donor])
    spectrum = (
        1.0f0, 0.0f0, 0.0f0, 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    )
    temporal_trace = Float32[
        1.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
    ]

    values = Pioneer._mbr_feature_values(
        UInt32(10),
        100.0f0,
        -1.0f0,
        10.0f0,
        10.0f0,
        5.0f0,
        spectrum,
        temporal_trace,
        UInt8(0x03),
        donor,
        UInt32(3),
        nothing,
        Pioneer._MBRReceiverRunClusters(),
        nothing,
        donor_dict,
    )

    @test length(values) == length(Pioneer.MBR_PAIRED_FEATURE_STEMS)
    @test values[19] == 1.0f0
    @test values[20] == donor.trace_prob
end

@testset "paired post-integration MBR model produces OOF transfer scores" begin
    n_baseline = 100
    n_candidates = 40
    n = n_baseline + n_candidates
    frame = DataFrame(
        precursor_idx = UInt32.(1:n),
        scan_idx = UInt32.(1001:(1000 + n)),
        ms_file_idx = UInt32.(1 .+ ((0:(n - 1)) .% 4)),
        cv_fold = UInt8.((0:(n - 1)) .% 2),
        target = vcat(trues(n_baseline), trues(n_candidates)),
        qval = vcat(
            fill(0.005f0, n_baseline),
            fill(0.20f0, n_candidates),
        ),
        global_qval = fill(0.005f0, n),
        trace_prob_prepass = vcat(
            fill(0.95f0, n_baseline),
            fill(0.60f0, n_candidates),
        ),
        trace_prob_infold = vcat(
            fill(0.95f0, n_baseline),
            Float32.(range(0.55, 0.75; length = n_candidates)),
        ),
    )
    for feature in Pioneer.MBR_RECEIVER_FEATURES
        feature === :trace_prob_infold && continue
        frame[!, feature] = fill(0.5f0, n)
    end
    frame[!, :MBR_best_is_missing_true] =
        vcat(trues(n_baseline), falses(n_candidates))
    for stem in Pioneer.MBR_PAIRED_FEATURE_STEMS
        frame[!, Pioneer._mbr_true_feature(stem)] =
            vcat(fill(-1.0f0, n_baseline), fill(0.9f0, n_candidates))
    end
    for counterfactual_idx in 1:Pioneer.MBR_N_COUNTERFACTUALS
        frame[
            !,
            Pioneer._mbr_missing_feature(counterfactual_idx),
        ] = vcat(trues(n_baseline), falses(n_candidates))
        for stem in Pioneer.MBR_PAIRED_FEATURE_STEMS
            frame[
                !,
                Pioneer._mbr_false_feature(stem, counterfactual_idx),
            ] = vcat(
                fill(-1.0f0, n_baseline),
                fill(0.1f0 + 0.01f0 * counterfactual_idx, n_candidates),
            )
        end
    end

    summary = Pioneer.apply_postintegration_mbr_rescoring!(
        frame;
        alpha = 0.01f0,
        q_value_threshold = 0.01f0,
    )
    candidate_rows = (n_baseline + 1):n
    @test summary.n_candidates == n_candidates
    @test all(frame.MBR_transfer_candidate[candidate_rows])
    @test all(isfinite, frame.ftr_qval_true[candidate_rows])
    @test all(isfinite, frame.ftr_pep_true[candidate_rows])
end

@testset "hardest MBR counterfactual control retains score and block index" begin
    # Scores are block-stacked: real, false-1, false-2, false-3.
    scores = Float32[
        0.9, 0.8,
        0.2, 0.7,
        0.6, 0.4,
        0.5, 0.95,
    ]
    present = BitMatrix(Bool[
        true true true
        true false true
    ])
    controls = Pioneer._mbr_top_counterfactual_controls(scores, present)

    @test controls.scores == Float32[0.6, 0.95]
    @test controls.indices == UInt8[2, 3]
    @test Pioneer._mbr_top_counterfactual_scores(scores, present) ==
        controls.scores
end

@testset "MBR training caps bound work, not just output" begin
    # The capped branch used to push every eligible row index into per-class Vector{Int}s and then
    # sample those down, so hitting the cap allocated 8 bytes per *available* row to keep `limit` of
    # them. Selecting 1.875M negatives out of 60M eligible would have cost ~480 MB of index vectors.
    # Allocation must now track the cap, not the eligible population.
    function _cap_alloc(n_eligible::Int)
        ncf = Pioneer.MBR_N_COUNTERFACTUALS
        folds = fill(UInt8(1), n_eligible)
        positives = trues(n_eligible)
        receivers = falses(n_eligible)
        present = trues(n_eligible, ncf)
        args = (folds, positives, receivers, present, UInt8(1))
        Pioneer._mbr_training_rows(args...; max_positives = 50, max_negatives = 50)
        GC.gc()
        return @allocated Pioneer._mbr_training_rows(
            args...; max_positives = 50, max_negatives = 50)
    end
    small = _cap_alloc(20_000)
    large = _cap_alloc(200_000)
    # 10x the eligible rows at the same cap must not cost materially more memory.
    @test large < 2 * small

    # And the cap is genuinely enforced at that scale.
    ncf = Pioneer.MBR_N_COUNTERFACTUALS
    n = 100_000
    training = Pioneer._mbr_training_rows(
        fill(UInt8(1), n), trues(n), falses(n), trues(n, ncf), UInt8(1);
        max_positives = 37, max_negatives = 41,
    )
    @test training.used_positives == 37
    @test training.used_negatives == 41
    @test count(training.labels) == 37
    @test count(!, training.labels) == 41
    @test length(unique(training.rows)) == length(training.rows)
    @test all(1 .<= training.rows .<= (1 + ncf) * n)
    @test sum(training.used_negatives_by_counterfactual) == 41
end

@testset "MBR row gather matches the block-stacked matrix" begin
    # `_mbr_feature_matrix` used to build the whole (1 + NCF) * n_candidates x n_features expansion,
    # but it was only ever consumed as x[train_rows, :] / x[test_rows, :]. `_mbr_gather_feature_rows`
    # produces those subsets directly, so it must agree with the full matrix exactly -- including the
    # missing -> 0.0f0 convention and the global row numbering.
    ncf = Pioneer.MBR_N_COUNTERFACTUALS
    n_candidates = 37
    n_features = 4
    true_features = [Symbol("gt_$j") for j in 1:n_features]
    false_features = [[Symbol("gf_$(k)_$j") for j in 1:n_features] for k in 1:ncf]
    df = DataFrames.DataFrame()
    for j in 1:n_features
        df[!, true_features[j]] = Float32.(collect(1:n_candidates) .+ 100j)
        for k in 1:ncf
            col = Vector{Union{Missing, Float32}}(Float32.(collect(1:n_candidates) .+ 1000k .+ 10j))
            col[k] = missing                      # exercise the missing -> 0.0f0 path
            df[!, false_features[k][j]] = col
        end
    end
    x = Pioneer._mbr_feature_matrix(df, true_features, false_features)
    n_rows = (1 + ncf) * n_candidates
    @test size(x) == (n_rows, n_features)
    @test Pioneer._mbr_gather_feature_rows(
        df, true_features, false_features, collect(1:n_rows), n_candidates) == x
    scattered = [n_rows, 1, n_candidates + 2, 7, 7, n_rows - 1]
    @test Pioneer._mbr_gather_feature_rows(
        df, true_features, false_features, scattered, n_candidates) == x[scattered, :]
    view_rows = @view collect(1:n_rows)[3:11]
    @test Pioneer._mbr_gather_feature_rows(
        df, true_features, false_features, view_rows, n_candidates) == x[collect(view_rows), :]
    @test_throws Exception Pioneer._mbr_gather_feature_rows(
        df, true_features, false_features, [n_rows + 1], n_candidates)
end
