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
        charge,
        sequence_length,
        irt,
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
    false_donors = Pioneer._mbr_false_donors(
        donor_dict,
        pools,
        eligibility,
        UInt32(1),
        UInt32(1),
        10.0f0,
        true_donor,
        nothing,
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
    present = falses(n_candidates, Pioneer.MBR_N_COUNTERFACTUALS)
    present[:, 1] .= true
    present[1:8, 2] .= true
    present[1:4, 3] .= true

    sampled = Pioneer._mbr_training_rows(
        folds,
        positives,
        present,
        UInt8(1);
        max_positives = 5,
        max_negatives = 7,
    )
    sampled_again = Pioneer._mbr_training_rows(
        folds,
        positives,
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
