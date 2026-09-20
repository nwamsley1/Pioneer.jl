using Arrow, DataFrames, Pioneer, Test

function _write_mbr_candidate_fixture(directory, name, qval, global_qval, donor_missing, target)
    n = length(qval)
    path = joinpath(directory, name * ".arrow")
    ids = (precursor_idx = UInt32.(1:n), scan_idx = UInt32.(101:(100 + n)))
    main = DataFrame(;
        ids..., ms_file_idx = fill(UInt32(1), n), cv_fold = UInt8.((1:n) .% 2),
        target, qval, global_qval, irt_error = Float32.(1:n),
    )
    open(Arrow.Writer, path; file=true) do writer
        Arrow.write(writer, main[1:1, :])
        Arrow.write(writer, main[2:end, :])
    end
    Arrow.write(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(;
        ids..., trace_prob_prepass = fill(0.7f0, n), trace_prob_infold = fill(0.8f0, n),
    ))
    features = DataFrame(;
        ids..., MBR_best_is_missing_true = donor_missing,
        MBR_best_is_missing_false = falses(n),
        MBR_best_pair_prob_true = Float32.(1:n) ./ 10,
        MBR_log2_weight_lod_ratio = Float32.(1:n),
    )
    Arrow.write(path * Pioneer.MBR_SIDECAR_SUFFIX, features)
    return path, features
end

@testset "candidate-only MBR feature sidecars" begin
    mktempdir() do directory
        empty_path, empty_features = _write_mbr_candidate_fixture(
            directory, "empty", fill(0.005f0, 3), fill(0.005f0, 3),
            falses(3), Bool[true, false, true],
        )
        mixed_path, mixed_features = _write_mbr_candidate_fixture(
            directory, "mixed", Float32[0.005, 0.2, 0.2, 0.3, 0.005, NaN, 0.2],
            Float32[0.005, 0.005, 0.2, 0.005, 0.005, 0.005, NaN],
            Bool[false, false, false, true, false, false, false],
            Bool[true, true, true, true, false, false, true],
        )
        all_path, all_features = _write_mbr_candidate_fixture(
            directory, "all", fill(0.2f0, 3), fill(0.005f0, 3),
            falses(3), Bool[true, false, true],
        )
        paths = [empty_path, mixed_path, all_path]
        rows_by_file = [Int[], [2, 6], [1, 2, 3]]
        dense = Pioneer.load_postintegration_mbr_candidates(paths, 0.01f0)
        sparse_features = DataFrame[]
        for (path, features, rows) in zip(paths, [empty_features, mixed_features, all_features], rows_by_file)
            sparse = features[rows, :]
            insertcols!(sparse, 1, :row_idx => Int64.(rows))
            Arrow.write(path * Pioneer.MBR_SIDECAR_SUFFIX, sparse)
            push!(sparse_features, sparse)
        end

        sparse = Pioneer.load_postintegration_mbr_candidates(paths, 0.01f0)
        @test isequal(sparse.candidates, dense.candidates)
        @test sparse.masks == dense.masks
        @test findall.(sparse.masks) == rows_by_file
        @test sparse.n_rows == [3, 7, 3]
        @test (sparse.base_targets, sparse.base_decoys) == (3, 2)
        @test !hasproperty(sparse.candidates, :row_idx)

        full = Pioneer.load_postintegration_mbr_frame(paths)
        full_mask = vcat(sparse.masks...)
        @test nrow(full) == sum(sparse.n_rows)
        @test isequal(full[full_mask, :], sparse.candidates)
        @test all(full.MBR_best_is_missing_true[.!full_mask])
        @test all(full.MBR_best_is_missing_false[.!full_mask])
        @test all(==(-1.0f0), full.MBR_best_pair_prob_true[.!full_mask])
        @test all(==(-1.0f0), full.MBR_log2_weight_lod_ratio[.!full_mask])

        candidates = copy(sparse.candidates)
        n = nrow(candidates)
        candidates[!, :mbr_recovered] = trues(n)
        candidates[!, :MBR_transfer_candidate] = trues(n)
        for column in (:mbr_target_decoy_prob, :ftr_qval_true, :ftr_pep_true,
                       :mbr_total_error_qval_true, :mbr_total_error_rate_true,
                       :mbr_counterfactual_decoy_prob)
            candidates[!, column] = Float32.(1:n) ./ 10
        end
        candidates[!, :mbr_counterfactual_decoy_index] = ones(UInt8, n)
        @test Pioneer._write_mbr_recovery_sidecars_from_candidates!(
            candidates, sparse.masks, sparse.n_rows, paths,
        ) == length(paths)
        cursor = 0
        for (path, rows) in zip(paths, rows_by_file)
            main = Arrow.Table(path)
            recovered = DataFrame(Arrow.Table(path * Pioneer.RECOVERY_SIDECAR_SUFFIX))
            @test recovered.precursor_idx == main.precursor_idx
            @test recovered.scan_idx == main.scan_idx
            @test findall(recovered.MBR_transfer_candidate) == rows
            @test findall(recovered.mbr_recovered) == rows
            @test recovered.mbr_target_decoy_prob[rows] ==
                candidates.mbr_target_decoy_prob[(cursor + 1):(cursor + length(rows))]
            @test all(isnan, recovered.mbr_target_decoy_prob[.!recovered.MBR_transfer_candidate])
            cursor += length(rows)
        end

        empty = Pioneer.load_postintegration_mbr_candidates([empty_path], 0.01f0)
        @test nrow(empty.candidates) == 0
        Pioneer.apply_postintegration_mbr_rescoring!(empty.candidates;
            alpha = 0.01f0, q_value_threshold = 0.01f0,
            baseline_counts = (empty.base_targets, empty.base_decoys), frame_is_candidates = true,
        )
        @test Pioneer._write_mbr_recovery_sidecars_from_candidates!(
            empty.candidates, empty.masks, empty.n_rows, [empty_path],
        ) == 1

        @testset "invalid row mappings are rejected" begin
            sidecar_path = mixed_path * Pioneer.MBR_SIDECAR_SUFFIX
            for rows in ([0, 6], [2, 8], [2, 2], [6, 2])
                invalid = copy(sparse_features[2])
                invalid.row_idx = Int64.(rows)
                Arrow.write(sidecar_path, invalid)
                @test_throws ErrorException Pioneer.load_postintegration_mbr_candidates([mixed_path], 0.01f0)
                @test_throws ErrorException Pioneer.load_postintegration_mbr_frame([mixed_path])
            end
            invalid = copy(sparse_features[2])
            invalid.scan_idx[1] += UInt32(1)
            Arrow.write(sidecar_path, invalid)
            @test_throws ErrorException Pioneer.load_postintegration_mbr_candidates([mixed_path], 0.01f0)
            @test_throws ErrorException Pioneer.load_postintegration_mbr_frame([mixed_path])
            for row in (1, 3, 4, 7)
                invalid = mixed_features[[row], :]
                insertcols!(invalid, 1, :row_idx => Int64[row])
                Arrow.write(sidecar_path, invalid)
                @test_throws ErrorException Pioneer.load_postintegration_mbr_candidates([mixed_path], 0.01f0)
            end
        end

        @testset "feature producer writes candidate rows only" begin
            main = DataFrame(Arrow.Table(mixed_path); copycols=true)
            n = nrow(main)
            spectrum = (1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0)
            trace = Float32[1, 1, 0, 0, 0, 0, 0, 0, 0]
            for (column, value) in zip(Pioneer.MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS, spectrum)
                main[!, column] = fill(value, n)
            end
            main[!, Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN] = [copy(trace) for _ in 1:n]
            for (column, value) in (
                Pioneer.MBR_INTEGRATED_WEIGHT_COLUMN => 100.0f0,
                Pioneer.MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN => -1.0f0,
                Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN => 10.0f0,
                Pioneer.MBR_INTEGRATED_N_SCANS_COLUMN => 5.0f0,
                Pioneer.MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN => UInt8(0x03),
                Pioneer.MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN => UInt16(2),
            )
                main[!, column] = fill(value, n)
            end
            Arrow.write(mixed_path, main)
            donors = Dict(pid => [Pioneer._MBRDonorEntry(
                0.9f0, pid, 100.0f0, -1.0f0, 0.0f0, 10.0f0, 5.0f0,
                spectrum, UInt8(0x03), UInt16(2), UInt32(2),
            )] for pid in UInt32[1, 2, 3, 6, 7])
            pools = Pioneer._MBRPartnerPools(
                zeros(UInt8, n), ones(UInt8, n), fill(UInt8(2), n), fill(UInt8(9), n),
                fill(10.0f0, n),
                Dict{NTuple{4, Int}, Pioneer._MBRIrtPool}(),
                Dict{NTuple{3, Int}, Pioneer._MBRIrtPool}(),
                Dict{NTuple{2, Int}, Pioneer._MBRIrtPool}(),
                Dict{Tuple{UInt32, Int, Int}, Pioneer._MBRIrtPool}(),
            )
            eligibility = Pioneer._MBRCounterfactualEligibility(BitSet(1:n), Dict(UInt32(1) => BitSet()))
            clusters = Pioneer._MBRReceiverRunClusters()
            kwargs = (run_similarity_atlas = nothing, receiver_run_clusters = clusters,
                      lod_log2_weight_by_file = Dict{UInt32, Float32}(),
                      lod_log2_weight_global = 0.0f0, q_value_threshold = 0.01f0)
            stats = Ref{Any}()
            Pioneer.compute_postintegration_mbr_features!(
                mixed_path, Pioneer._MBRDonorIndex(donors), pools, eligibility;
                kwargs..., stats=stats,
            )
            @test stats[].rows == n
            @test stats[].selection_seconds >= 0
            @test stats[].feature_seconds >= 0
            @test stats[].write_seconds >= 0
            produced = DataFrame(Arrow.Table(mixed_path * Pioneer.MBR_SIDECAR_SUFFIX))
            @test stats[].candidates == 2
            @test produced.row_idx == Int64[2, 6]
            @test produced.precursor_idx == main.precursor_idx[[2, 6]]
            @test produced.scan_idx == main.scan_idx[[2, 6]]
            @test !any(produced.MBR_best_is_missing_true)
            @test all(produced.MBR_best_is_missing_false)
            expected = Pioneer._mbr_feature_values(
                UInt32(2), 100.0f0, -1.0f0, 10.0f0, 10.0f0, 5.0f0,
                spectrum, trace, UInt8(0x03), only(donors[UInt32(2)]),
                UInt32(1), nothing, clusters, nothing, donors,
            )
            true_columns = Pioneer.MBR_PAIRED_COLUMN_NAMES[1:length(expected)]
            @test isequal(Tuple(produced[1, true_columns]), expected)

            main.qval .= 0.005f0
            Arrow.write(mixed_path, main)
            stats = Ref{Any}()
            Pioneer.compute_postintegration_mbr_features!(
                mixed_path, Pioneer._MBRDonorIndex(donors), pools, eligibility;
                kwargs..., stats=stats,
            )
            @test stats[].rows == n
            @test stats[].selection_seconds >= 0
            @test stats[].feature_seconds >= 0
            @test stats[].write_seconds >= 0
            produced_empty = DataFrame(Arrow.Table(mixed_path * Pioneer.MBR_SIDECAR_SUFFIX))
            @test stats[].candidates == 0
            @test nrow(produced_empty) == 0
            @test names(produced_empty) == names(produced)
            @test nrow(Pioneer.load_postintegration_mbr_candidates([mixed_path], 0.01f0).candidates) == 0
        end
    end
end
