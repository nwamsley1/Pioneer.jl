using Test
using DataFrames
using Arrow
using Interpolations
using Pioneer

function _apply_pipeline_to_df(df::DataFrame, pipeline)
    result = copy(df)
    for (_, op) in pipeline.operations
        result = op(result)
        nrow(result) == 0 && break
    end
    return result
end

@testset "MBR rescue candidate helpers" begin
    @testset "main-search PEP partition bounds MBR rescue candidates" begin
        @test Pioneer.MAIN_MBR_RESCUE_PEP_MAX == 0.99f0

        peps = Float32[0.0, 0.9, 0.9001, 0.98, 0.99, 1.0]

        part = Pioneer._mainsearch_pep_partition(peps, 0.9f0, 0.99f0)

        @test part.keep_mask == Bool[true, true, false, false, false, false]
        @test part.rescue_mask == Bool[false, false, true, true, true, false]
        @test part.discard_mask == Bool[false, false, false, false, false, true]
    end

    @testset "raw RT MBR features are no longer emitted or modeled" begin
        @test :MBR_best_rt_diff_true ∉ Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_rt_diff_false ∉ Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_best_rt_diff_true ∉ Pioneer._MBR_SIDECAR_OUT_COLS
        @test :MBR_best_rt_diff_false ∉ Pioneer._MBR_SIDECAR_OUT_COLS
    end

    @testset "FTR includes per-run and cross-run probability features" begin
        @test :main_search_prob in Pioneer.FTR_FEATURES_F_TRUE
        @test :trace_prob_infold in Pioneer.FTR_FEATURES_F_TRUE
        @test :main_search_prob in Pioneer.FTR_FEATURES_F_FALSE
        @test :trace_prob_infold in Pioneer.FTR_FEATURES_F_FALSE
    end

    @testset "MBR semi-supervised training keeps only FTR/FDR target positives" begin
        positive_top = Bool[true, false, false, false]
        target_top = Bool[true, true, false, false]

        @test Pioneer.MBR_SEMISUPERVISED_FTR_THRESHOLD == 0.03f0
        @test Pioneer.MBR_SEMISUPERVISED_FDR_THRESHOLD == 0.01f0
        @test Pioneer._mbr_transfer_training_labels(positive_top) ==
              Bool[true, false, false, false, false, false, false, false]
        @test Pioneer._mbr_transfer_training_mask(positive_top, target_top) ==
              Bool[true, false, true, true, true, true, true, true]

        scores_double = Float32[1.0; collect(range(0.99f0, 0.5f0; length=40)); zeros(Float32, 41)]
        eval_labels = Bool[false; trues(40); falses(41)]
        target_top = Bool[false; trues(40)]
        metrics = Pioneer._mbr_transfer_iteration_metrics(
            scores_double,
            eval_labels,
            target_top,
            41;
            ftr_threshold = 1.0f0,
        )
        @test !metrics.positive_top[41]
    end

    @testset "FTR includes MBR scan-count comparison features" begin
        @test :MBR_abs_n_scans_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_log2_n_scans_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_abs_n_scans_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_log2_n_scans_ratio_false in Pioneer.FTR_FEATURES_F_FALSE
    end

    @testset "FTR includes observed iRT comparison features" begin
        @test :MBR_best_irt_obs_diff_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_best_irt_obs_diff_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_best_irt_obs_diff_true in Pioneer._MBR_SIDECAR_OUT_COLS
        @test :MBR_best_irt_obs_diff_false in Pioneer._MBR_SIDECAR_OUT_COLS
    end

    @testset "FTR features include MBR and row-level evidence for transfer ranking" begin
        @test :main_search_prob in Pioneer.FTR_FEATURES_F_TRUE
        @test :trace_prob_infold in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_max_pair_prob_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_log2_weight_ratio_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_log2_weight_lod_ratio in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_smoothed_frag_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_donor_library_hellinger_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_donor_library_hellinger_worst_true in Pioneer.FTR_FEATURES_F_TRUE
        @test :MBR_max_pair_prob_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_donor_library_hellinger_false in Pioneer.FTR_FEATURES_F_FALSE
        @test :MBR_donor_library_hellinger_worst_false in Pioneer.FTR_FEATURES_F_FALSE
    end

    @testset "counterfactual iRT window uses robust observed residuals" begin
        @test Pioneer._mbr_counterfactual_irt_window!(Float32[]) == 1.0f0
        @test Pioneer._mbr_counterfactual_irt_window!(Float32[0.01, 0.02, 0.03]) == 0.25f0
        @test Pioneer._mbr_counterfactual_irt_window!(Float32[0.2, 0.4, 0.6, 8.0]) == 1.0f0
        @test Pioneer._mbr_counterfactual_irt_window!(Float32[3.0, 4.0, 5.0]) == 2.0f0
    end

    @testset "rescue path discovery mirrors main_search_psms layout" begin
        mktempdir() do dir
            normal_dir = joinpath(dir, "main_search_psms")
            rescue_dir = joinpath(dir, "mbr_rescue_psms")
            mkpath(normal_dir)
            mkpath(rescue_dir)
            normal_path = joinpath(normal_dir, "run1_fold0.arrow")
            rescue_path = joinpath(rescue_dir, "run1_fold0.arrow")
            Pioneer.writeArrow(normal_path, DataFrame(precursor_idx = UInt32[1]))
            Pioneer.writeArrow(rescue_path, DataFrame(precursor_idx = UInt32[2]))

            @test Pioneer.mbr_rescue_fold_path(normal_path) == rescue_path
            @test Pioneer.get_existing_mbr_rescue_fold_paths([normal_path]) == [rescue_path]
        end
    end

    @testset "rescue candidate slim rows use per-run 1-PEP confidence without dense sidecars" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "run1_fold0.arrow")
            smoothed_cols = Pair{Symbol, Vector{Float32}}[
                col => Float32[9, 4] for col in Pioneer.SMOOTHED_FRAGMENT_INTENSITY_COLUMNS
            ]
            Pioneer.writeArrow(rescue_path, DataFrame(
                :precursor_idx => UInt32[1, 3],
                :scan_idx => UInt32[10, 30],
                :ms_file_idx => UInt32[1, 1],
                :cv_fold => UInt8[0, 0],
                :target => Bool[true, true],
                :main_pep => Float32[0.93, 0.94],
                :weight => Float32[8, 4],
                :log2_intensity_explained => Float16[3, 2],
                :irt_pred => Float32[10, 30],
                :irt_obs => Float32[9, 29],
                :log_by_ratio_m0 => Float16[0.5, 0.2],
                :n_scans => Float32[6, 2],
                :smoothed_2d_shadow_hellinger => Float32[0.12, 0.22],
                smoothed_cols...,
            ))
            donor_sqrt = ntuple(_ -> sqrt(1.0f0 / 8.0f0), 8)
            donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
                UInt32(1) => [Pioneer._MBRDonorEntry(
                    0.95f0, UInt32(1), 10.0f0, 2.0f0, 0.25f0, 8.8f0,
                    0.1f0, 4.0f0, donor_sqrt, 0.12f0, UInt32(2), false,
                )],
                UInt32(2) => [Pioneer._MBRDonorEntry(
                    0.94f0, UInt32(2), 12.0f0, 1.0f0, 0.5f0, 8.5f0,
                    0.3f0, 3.0f0, donor_sqrt, 0.22f0, UInt32(2), false,
                )],
            )
            pool = (pids = UInt32[2], irts = Float32[10])
            partner_pools = Pioneer._CounterfactualPartnerPools(
                Bool[true, false, true],
                UInt8[0, 0, 0],
                UInt8[1, 1, 1],
                UInt8[2, 2, 2],
                UInt8[10, 10, 10],
                Float32[10, 10, 30],
                Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}((0, 1, 2, 10) => pool),
                Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 10) => pool),
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 10) => pool),
            )
            fragment_keys = Pioneer._MBRFragmentAnnotationKeys(zeros(UInt16, 8 * 3))

            slim = Pioneer.load_mbr_rescue_candidate_slim_dataframe(
                [rescue_path],
                donor_dict,
                partner_pools,
                fragment_keys,
                0.90f0;
                lod_log2_weight_by_file = Dict(UInt32(1) => Float32(log2(10.0f0))),
                lod_log2_weight_global = Float32(log2(10.0f0)),
            )

            @test !isfile(rescue_path * Pioneer.PASS1_SIDECAR_SUFFIX)
            @test slim.precursor_idx == UInt32[1]
            @test slim.scan_idx == UInt32[10]
            @test slim.trace_prob_prepass == Float32[0.07]
            @test slim.main_search_prob == Float32[0.07]
            @test slim.trace_prob_infold == Float32[0.0]
            @test slim.mbr_rescue_candidate == Bool[true]
            @test slim._mbr_rescue_file_idx == UInt32[1]
            @test slim._mbr_rescue_row_idx == UInt32[1]
            @test slim.MBR_is_missing_true == Bool[false]
            @test slim.MBR_is_missing_false == Bool[false]
            @test slim.MBR_abs_n_scans_diff_true == Float32[2]
            @test slim.MBR_abs_n_scans_diff_false == Float32[3]
            @test isapprox(
                slim.MBR_log2_n_scans_ratio_true[1],
                Float32(log2((6.0f0 + 1.0f0) / (4.0f0 + 1.0f0)));
                atol = 1.0f-6,
            )
            @test isapprox(
                slim.MBR_log2_n_scans_ratio_false[1],
                Float32(log2((6.0f0 + 1.0f0) / (3.0f0 + 1.0f0)));
                atol = 1.0f-6,
            )
        end
    end

    @testset "rescue rows get the smoothed fragment columns required by MBR" begin
        psms = DataFrame(
            precursor_idx = UInt32[1, 1, 1],
            scan_idx = UInt32[10, 20, 30],
            cycle_idx = UInt32[1, 2, 3],
            frag1_int = Float32[2, 10, 6],
            frag2_int = Float32[0, 4, 0],
            frag3_int = Float32[1, 1, 1],
            frag4_int = Float32[0, 0, 0],
            frag5_int = Float32[0, 0, 0],
            frag6_int = Float32[0, 0, 0],
            frag7_int = Float32[0, 0, 0],
            frag8_int = Float32[0, 0, 0],
        )
        rescue_psms = DataFrame(
            precursor_idx = UInt32[1],
            scan_idx = UInt32[20],
        )
        center_mzs = fill(500.0f0, 30)
        isolation_widths = fill(8.0f0, 30)

        Pioneer.add_smoothed_fragment_intensities!(
            rescue_psms,
            psms,
            trues(nrow(psms));
            center_mzs = center_mzs,
            isolation_widths = isolation_widths,
        )

        @test rescue_psms.frag1_smoothed_intensity == Float32[7]
        @test rescue_psms.frag2_smoothed_intensity == Float32[2]
        @test rescue_psms.frag3_smoothed_intensity == Float32[1]
        for col in Pioneer.SMOOTHED_FRAGMENT_INTENSITY_COLUMNS
            @test hasproperty(rescue_psms, col)
        end
    end

    @testset "rescue fold files get deferred precursor metadata" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "run1_fold0.arrow")
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[2],
                irt_obs = Float32[12],
                irt_pred = Float32[10],
                rt = Float32[5],
            ))

            result = Pioneer._enrich_mainsearch_phase2_fold_file!(
                rescue_path,
                Float32[0, 10],
                Float32[0, 450],
                UInt32[0, 123],
                UInt8[0, 7],
            )
            tbl = DataFrame(Arrow.Table(rescue_path))

            @test result.n_rows == 1
            @test tbl.irt_diff == Float32[2]
            @test tbl.prec_mz == Float32[450]
            @test tbl.pair_id == UInt32[123]
            @test tbl.entrapment_group_id == UInt8[7]
            @test tbl.trace_prob == Float32[0]
            @test tbl.q_value == Float32[0]
        end
    end

    @testset "candidate gate uses normal rows for donor threshold and admits rescue rows" begin
        psms = DataFrame(
            trace_prob_prepass = Float32[0.99, 0.20, 0.99],
            target = Bool[true, false, true],
            mbr_rescue_candidate = Bool[false, false, true],
            MBR_max_pair_prob_true = Float32[0.0, 0.0, 0.99],
            MBR_is_missing_true = Bool[true, true, false],
            MBR_is_missing_false = Bool[true, true, false],
        )

        gate = Pioneer._mbr_candidate_gate(psms; q_thresh = 0.01f0)

        @test gate.rescue_mask == Bool[false, false, true]
        @test gate.prob_thresh == 0.99f0
        @test gate.pre_qvals[3] == 1.0f0
        @test gate.candidate_mask == Bool[false, false, true]
    end

    @testset "normal MBR candidates are sparse and do not require dense MBR sidecars" begin
        mktempdir() do dir
            normal_path = joinpath(dir, "run1_fold0.arrow")
            smoothed_cols = Pair{Symbol, Vector{Float32}}[
                col => Float32[9, 4, 7] for col in Pioneer.SMOOTHED_FRAGMENT_INTENSITY_COLUMNS
            ]
            Pioneer.writeArrow(normal_path, DataFrame(
                :precursor_idx => UInt32[1, 2, 3],
                :scan_idx => UInt32[10, 20, 30],
                :ms_file_idx => UInt32[1, 1, 1],
                :cv_fold => UInt8[0, 0, 0],
                :target => Bool[true, false, true],
                :main_pep => Float32[0.42, 0.10, 0.01],
                :weight => Float32[8, 4, 9],
                :log2_intensity_explained => Float16[3, 2, 4],
                :irt_pred => Float32[10, 20, 30],
                :irt_obs => Float32[9, 19, 29],
                :log_by_ratio_m0 => Float16[0.5, 0.2, 0.7],
                :n_scans => Float32[7, 3, 9],
                :smoothed_2d_shadow_hellinger => Float32[0.05, 0.20, 0.01],
                smoothed_cols...,
            ))
            Pioneer.writeArrow(normal_path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = UInt32[1, 2, 3],
                scan_idx = UInt32[10, 20, 30],
                trace_prob_prepass = Float32[0.80, 0.95, 0.99],
                trace_prob_infold = Float32[0.81, 0.96, 0.995],
            ))
            computed_prepass = Pioneer._normal_mbr_prepass_qvals_and_threshold([normal_path])
            @test computed_prepass.lod_log2_weight_by_file[UInt32(1)] == Float32(log2(9.0f0))
            @test computed_prepass.lod_log2_weight_global == Float32(log2(9.0f0))
            @test computed_prepass.counterfactual_irt_window == 2.0f0
            @test computed_prepass.n_counterfactual_irt_window_samples == 1

            donor_sqrt = ntuple(_ -> sqrt(1.0f0 / 8.0f0), 8)
            worst_donor_sqrt = (1.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0)
            donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
                UInt32(1) => [
                    Pioneer._MBRDonorEntry(
                        0.99f0, UInt32(1), 10.0f0, 2.0f0, 0.25f0, 8.8f0,
                        0.1f0, 4.0f0, donor_sqrt, 0.11f0, UInt32(2), false,
                    ),
                    Pioneer._MBRDonorEntry(
                        0.98f0, UInt32(1), 4.0f0, 1.0f0, 0.4f0, 8.2f0,
                        0.2f0, 2.0f0, worst_donor_sqrt, 0.33f0, UInt32(3), false,
                    ),
                ],
                UInt32(2) => [
                    Pioneer._MBRDonorEntry(
                        0.99f0, UInt32(2), 12.0f0, 1.0f0, 0.5f0, 8.5f0,
                        0.3f0, 3.0f0, donor_sqrt, 0.22f0, UInt32(2), false,
                    ),
                    Pioneer._MBRDonorEntry(
                        0.98f0, UInt32(2), 2.0f0, 0.5f0, 0.6f0, 8.1f0,
                        0.4f0, 1.0f0, worst_donor_sqrt, 0.44f0, UInt32(3), false,
                    ),
                ],
            )
            pool = (pids = UInt32[2], irts = Float32[10])
            partner_pools = Pioneer._CounterfactualPartnerPools(
                Bool[true, false, true],
                UInt8[0, 0, 0],
                UInt8[1, 1, 1],
                UInt8[2, 2, 2],
                UInt8[10, 10, 10],
                Float32[10, 10, 30],
                Dict{Tuple{Int, Int, Int, Int}, Pioneer._IrtPool}((0, 1, 2, 10) => pool),
                Dict{Tuple{Int, Int, Int}, Pioneer._IrtPool}((0, 2, 10) => pool),
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}((2, 10) => pool),
            )
            fragment_key_values = zeros(UInt16, 8 * 3)
            for pid in 1:3, rank in 1:8
                fragment_key_values[8 * (pid - 1) + rank] = UInt16(rank)
            end
            fragment_keys = Pioneer._MBRFragmentAnnotationKeys(fragment_key_values)
            prepass = (
                scores = Float32[0.80, 0.95, 0.99],
                qvals = Float32[0.02, 0.0, 0.0],
                row_counts = Int[3],
                prob_thresh = 0.95f0,
                counterfactual_irt_window = 2.0f0,
                n_counterfactual_irt_window_samples = 1,
                lod_log2_weight_by_file = Dict(UInt32(1) => Float32(log2(9.0f0))),
                lod_log2_weight_global = Float32(log2(9.0f0)),
            )

            result = Pioneer.load_normal_mbr_candidate_slim_dataframe(
                [normal_path],
                donor_dict,
                partner_pools,
                fragment_keys;
                q_thresh = 0.01f0,
                prepass = prepass,
            )

            @test !isfile(normal_path * ".mbr_sidecar.arrow")
            @test result.prob_thresh == 0.95f0
            @test result.n_rows == 3
            @test result.candidates.precursor_idx == UInt32[1]
            @test result.candidates.scan_idx == UInt32[10]
            @test result.candidates.trace_prob_prepass == Float32[0.80]
            @test result.candidates.main_search_prob == Float32[0.58]
            @test result.candidates.trace_prob_infold == Float32[0.81]
            @test result.candidates.mbr_rescue_candidate == Bool[false]
            @test result.candidates._mbr_normal_file_idx == UInt32[1]
            @test result.candidates._mbr_normal_row_idx == UInt32[1]
            @test result.candidates.MBR_is_missing_true == Bool[false]
            @test result.candidates.MBR_is_missing_false == Bool[false]
            @test result.candidates.MBR_abs_n_scans_diff_true == Float32[3]
            @test result.candidates.MBR_abs_n_scans_diff_false == Float32[4]
            @test result.candidates.MBR_abs_n_scans_diff_worst_true[1] == Float32(5)
            @test result.candidates.MBR_abs_n_scans_diff_worst_false[1] == Float32(6)
            @test isapprox(result.candidates.MBR_best_irt_obs_diff_true[1], 0.2f0; atol = 1.0f-6)
            @test isapprox(result.candidates.MBR_best_irt_obs_diff_false[1], 0.5f0; atol = 1.0f-6)
            @test result.candidates.MBR_log2_weight_lod_ratio[1] == Float32(log2(8.0f0) - log2(9.0f0))
            @test result.candidates.MBR_log2_weight_ratio_worst_true[1] == Float32(log2(8.0f0 / 4.0f0))
            @test result.candidates.MBR_log2_weight_ratio_worst_false[1] == Float32(log2(8.0f0 / 2.0f0))
            @test result.candidates.MBR_log2_explained_ratio_worst_true[1] == Float32(3.0f0 - 1.0f0)
            @test result.candidates.MBR_log2_explained_ratio_worst_false[1] == Float32(3.0f0 - 0.5f0)
            @test result.candidates.MBR_smoothed_frag_hellinger_worst_true[1] >
                  result.candidates.MBR_smoothed_frag_hellinger_true[1]
            @test result.candidates.MBR_smoothed_frag_hellinger_worst_false[1] >
                  result.candidates.MBR_smoothed_frag_hellinger_false[1]
            @test result.candidates.MBR_donor_library_hellinger_true == Float32[0.11]
            @test result.candidates.MBR_donor_library_hellinger_false == Float32[0.22]
            @test result.candidates.MBR_donor_library_hellinger_worst_true == Float32[0.33]
            @test result.candidates.MBR_donor_library_hellinger_worst_false == Float32[0.44]
            @test isapprox(
                result.candidates.MBR_log2_n_scans_ratio_true[1],
                Float32(log2((7.0f0 + 1.0f0) / (4.0f0 + 1.0f0)));
                atol = 1.0f-6,
            )
            @test isapprox(
                result.candidates.MBR_log2_n_scans_ratio_false[1],
                Float32(log2((7.0f0 + 1.0f0) / (3.0f0 + 1.0f0)));
                atol = 1.0f-6,
            )
        end
    end

    @testset "sparse normal FTR results expand to dense recovery sidecars" begin
        mktempdir() do dir
            normal_path = joinpath(dir, "normal_fold0.arrow")
            Pioneer.writeArrow(normal_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
            ))
            candidates = DataFrame(
                precursor_idx = UInt32[2],
                scan_idx = UInt32[20],
                _mbr_normal_file_idx = UInt32[1],
                _mbr_normal_row_idx = UInt32[2],
                mbr_recovered = Bool[true],
                MBR_transfer_candidate = Bool[true],
                mbr_target_decoy_prob = Float32[0.88],
                ftr_qval_true = Float32[0.005],
                ftr_pep_true = Float32[0.006],
            )

            result = Pioneer.write_sparse_normal_recovery_sidecars!(candidates, [normal_path])
            recovery = Arrow.Table(normal_path * Pioneer.RECOVERY_SIDECAR_SUFFIX)

            @test result.n_files == 1
            @test result.n_rows == 2
            @test collect(recovery.precursor_idx) == UInt32[1, 2]
            @test collect(recovery.scan_idx) == UInt32[10, 20]
            @test collect(recovery.mbr_recovered) == Bool[false, true]
            @test collect(recovery.MBR_transfer_candidate) == Bool[false, true]
            @test collect(recovery.mbr_target_decoy_prob)[1] |> isnan
            @test collect(recovery.mbr_target_decoy_prob)[2] == 0.88f0
            @test !hasproperty(recovery, :ftr_score_true)
            @test collect(recovery.ftr_qval_true)[1] |> isnan
            @test collect(recovery.ftr_qval_true)[2] == 0.005f0
            @test collect(recovery.ftr_pep_true)[1] |> isnan
            @test collect(recovery.ftr_pep_true)[2] == 0.006f0
        end
    end

    @testset "recovered rescue files contain only transferred rescue rows" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "rescue_fold0.arrow")
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                ms_file_idx = UInt32[1, 1],
                cv_fold = UInt8[0, 0],
                target = Bool[true, true],
                rt = Float32[1, 2],
            ))
            slim_df = DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                _mbr_rescue_file_idx = UInt32[1, 1],
                _mbr_rescue_row_idx = UInt32[1, 2],
                trace_prob_prepass = Float32[0.07, 0.08],
                trace_prob_infold = Float32[0.07, 0.08],
                mbr_recovered = Bool[false, true],
                MBR_transfer_candidate = Bool[true, true],
                mbr_target_decoy_prob = Float32[0.10, 0.88],
                ftr_qval_true = Float32[0.20, 0.005],
                ftr_pep_true = Float32[0.20, 0.006],
            )

            result = Pioneer.write_mbr_rescue_recovered_files!(slim_df, [rescue_path])
            recovered_path = Pioneer.mbr_rescue_recovered_path(rescue_path)
            recovered = DataFrame(Arrow.Table(recovered_path))

            @test result.n_files == 1
            @test result.n_rows == 1
            @test recovered.precursor_idx == UInt32[2]
            @test recovered.scan_idx == UInt32[20]
            @test recovered.trace_prob_prepass == Float32[0.08]
            @test recovered.trace_prob == Float32[0.08]
            @test recovered.mbr_recovered == Bool[true]
            @test recovered.MBR_transfer_candidate == Bool[true]
            @test recovered.mbr_target_decoy_prob == Float32[0.88]
            @test !hasproperty(recovered, :ftr_score_true)
            @test recovered.ftr_qval_true == Float32[0.005]
            @test recovered.ftr_pep_true == Float32[0.006]

            empty_result = Pioneer.write_mbr_rescue_recovered_files!(slim_df[1:0, :], [rescue_path])
            @test empty_result.n_files == 0
            @test empty_result.n_rows == 0
            @test !isfile(recovered_path)
        end
    end

    @testset "paired FTR emits target-decoy probabilities for recovered rows" begin
        psms = DataFrame(
            target = Bool[true, false, true, false, true, false, true, false],
            cv_fold = UInt8[0, 0, 0, 0, 1, 1, 1, 1],
            mbr_rescue_candidate = Bool[false, true, false, false, true, false, false, false],
            main_search_prob = Float32[0.95, 0.20, 0.90, 0.30, 0.94, 0.25, 0.91, 0.35],
            trace_prob_infold = Float32[0.96, 0.15, 0.88, 0.25, 0.93, 0.22, 0.87, 0.32],
            MBR_max_pair_prob_true = Float32[0.98, 0.30, 0.90, 0.20, 0.97, 0.35, 0.91, 0.25],
            MBR_max_pair_prob_false = Float32[0.40, 0.80, 0.35, 0.75, 0.45, 0.82, 0.36, 0.78],
            MBR_log2_weight_ratio_true = Float32[0.1, -2.0, 0.2, -1.5, 0.0, -1.8, 0.3, -1.2],
            MBR_log2_weight_ratio_false = Float32[-1.0, 0.1, -1.2, 0.2, -1.1, 0.0, -1.3, 0.3],
            MBR_log2_explained_ratio_true = Float32[0.2, -1.0, 0.1, -0.8, 0.2, -0.9, 0.1, -0.7],
            MBR_log2_explained_ratio_false = Float32[-0.8, 0.2, -0.7, 0.1, -0.9, 0.3, -0.6, 0.2],
            MBR_abs_n_scans_diff_true = Float32[0.0, 5.0, 1.0, 4.0, 0.0, 5.0, 1.0, 4.0],
            MBR_abs_n_scans_diff_false = Float32[4.0, 0.0, 5.0, 1.0, 4.0, 0.0, 5.0, 1.0],
            MBR_log2_n_scans_ratio_true = Float32[0.0, -2.0, 0.1, -1.5, 0.0, -1.8, 0.2, -1.2],
            MBR_log2_n_scans_ratio_false = Float32[-1.5, 0.0, -2.0, 0.1, -1.6, 0.0, -1.9, 0.2],
            MBR_best_irt_diff_true = Float32[0.1, 3.0, 0.2, 2.5, 0.1, 3.2, 0.2, 2.6],
            MBR_best_irt_diff_false = Float32[2.5, 0.1, 3.0, 0.2, 2.6, 0.1, 3.1, 0.2],
            MBR_log_by_diff_true = Float32[0.0, 1.5, 0.1, 1.2, 0.0, 1.6, 0.1, 1.3],
            MBR_log_by_diff_false = Float32[1.2, 0.0, 1.5, 0.1, 1.3, 0.0, 1.6, 0.1],
            MBR_smoothed_frag_hellinger_true = Float32[0.05, 0.6, 0.08, 0.5, 0.04, 0.55, 0.09, 0.45],
            MBR_smoothed_frag_hellinger_false = Float32[0.5, 0.05, 0.6, 0.08, 0.55, 0.04, 0.45, 0.09],
        )

        Pioneer.apply_mbr_filter_paired!(
            psms;
            alpha = 0.5f0,
            pregated = true,
            pregated_prob_thresh = 0.0f0,
        )

        @test hasproperty(psms, :mbr_target_decoy_prob)
        @test any((psms.ftr_qval_true .<= 0.5f0) .& (psms.ftr_pep_true .> 0.5f0))
        @test psms.mbr_recovered == (psms.ftr_qval_true .<= 0.5f0)
        @test all(!isnan, psms.mbr_target_decoy_prob[psms.mbr_recovered])
        @test all(isnan, psms.mbr_target_decoy_prob[.!psms.mbr_recovered])
        @test all(!isnan, psms.mbr_target_decoy_prob[psms.mbr_rescue_candidate .& psms.mbr_recovered])
    end

    @testset "MBR recovered rows bypass only the initial row q-value filter" begin
        initial_qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pep_interp = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer._precursor_scoring_qvalue_filter_pipeline(
            Dict(UInt32(1) => 0.95f0, UInt32(2) => 0.05f0),
            Dict(UInt32(1) => 0.01f0, UInt32(2) => 0.01f0),
            Dict(UInt32(1) => 0.01f0, UInt32(2) => 0.01f0),
            initial_qval_spline,
            pep_interp,
            0.1f0,
        )
        psms = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prec_prob = Float32[0.95, 0.05, 0.05],
            target = Bool[true, true, false],
            mbr_recovered = Bool[false, true, false],
            mbr_target_decoy_prob = Float32[NaN32, 0.90, NaN32],
        )

        initial_filtered = _apply_pipeline_to_df(psms, pipeline)

        @test initial_filtered.precursor_idx == UInt32[1, 2]
        @test initial_filtered.prec_prob[2] == 0.05f0
        @test initial_filtered.qval[2] > 0.1f0
    end

    @testset "MBR recovered rows are ranked by target-decoy score below direct row-level passers" begin
        qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer.TransformPipeline() |>
            Pioneer.add_interpolated_column(:qval, :prec_prob, qval_spline) |>
            Pioneer.remap_mbr_recovered_prec_probs(qval_spline, 0.1f0)
        psms = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],
            prec_prob = Float32[0.95, 0.05, 0.05, 0.04],
            target = Bool[true, true, true, false],
            mbr_recovered = Bool[false, true, true, false],
            mbr_target_decoy_prob = Float32[NaN32, 0.90, 0.30, NaN32],
            ftr_pep_true = Float32[NaN32, 0.009, 0.001, NaN32],
        )

        remapped = _apply_pipeline_to_df(psms, pipeline)
        direct_floor = Pioneer._score_floor_for_qvalue(qval_spline, 0.1f0)

        @test remapped.prec_prob[1] == 0.95f0
        @test remapped.prec_prob[4] == 0.04f0
        @test remapped.prec_prob[2] < direct_floor
        @test remapped.prec_prob[3] < direct_floor
        @test remapped.prec_prob[2] > remapped.prec_prob[3]
    end

    @testset "MBR recovered rows require non-MBR global q-value support" begin
        qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer.TransformPipeline() |>
            Pioneer.gate_mbr_recovered_by_global_qvalue(
                Dict(UInt32(2) => 0.01f0, UInt32(3) => 0.20f0),
                0.1f0,
            ) |>
            Pioneer.remap_mbr_recovered_prec_probs(qval_spline, 0.1f0)
        psms = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            prec_prob = Float32[0.95, 0.05, 0.06],
            target = Bool[true, true, true],
            mbr_recovered = Bool[false, true, true],
            mbr_target_decoy_prob = Float32[NaN32, 0.90, 0.90],
        )

        remapped = _apply_pipeline_to_df(psms, pipeline)

        @test remapped.mbr_recovered == Bool[false, true, false]
        @test remapped.prec_prob[2] > 0.05f0
        @test remapped.prec_prob[3] == 0.06f0
        @test isnan(remapped.mbr_target_decoy_prob[3])
    end

    @testset "non-MBR global scores ignore recovered rescue-only rows" begin
        mktempdir() do dir
            path = joinpath(dir, "psms.arrow")
            Pioneer.writeArrow(path, DataFrame(
                precursor_idx = UInt32[1, 2, 3],
                target = Bool[true, true, true],
                prec_prob = Float32[0.90, 0.70, 0.80],
                mbr_recovered = Bool[false, true, true],
                mbr_rescue_candidate = Bool[false, false, true],
            ))

            global_prob, target_dict = Pioneer.build_precursor_global_prob_dicts(
                [Pioneer.PSMFileReference(path)],
                1,
                3;
                exclude_mbr_rescue_recovered = true,
            )

            @test global_prob[UInt32(1)] == 0.90f0
            @test global_prob[UInt32(2)] == 0.70f0
            @test !haskey(global_prob, UInt32(3))
            @test !haskey(target_dict, UInt32(3))
        end
    end

    @testset "MBR recovered rows are filtered by recalculated row q-values" begin
        recalc_qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        recalc_pep_interp = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer._precursor_scoring_recalculated_qvalue_filter_pipeline(
            recalc_qval_spline,
            recalc_pep_interp,
            0.1f0,
        )
        psms = DataFrame(
            precursor_idx = UInt32[1, 2],
            prec_prob = Float32[0.95, 0.05],
            target = Bool[true, true],
            mbr_recovered = Bool[false, true],
        )

        filtered = _apply_pipeline_to_df(psms, pipeline)

        @test filtered.precursor_idx == UInt32[1]
        @test filtered.qval ≈ Float32[0.01]
    end

    @testset "recalculated row-level filter recomputes emitted PEP" begin
        recalc_qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        recalc_pep_interp = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.02];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer._precursor_scoring_recalculated_qvalue_filter_pipeline(
            recalc_qval_spline,
            recalc_pep_interp,
            0.1f0,
        )
        psms = DataFrame(
            precursor_idx = UInt32[1, 2],
            prec_prob = Float32[0.95, 0.05],
            pep = Float32[0.9999, 0.9999],
            target = Bool[true, true],
            mbr_recovered = Bool[true, true],
        )

        filtered = _apply_pipeline_to_df(psms, pipeline)

        @test filtered.precursor_idx == UInt32[1]
        @test filtered.pep[1] < 0.1f0
        @test filtered.pep[1] ≈ Float32(recalc_pep_interp(0.95f0))
    end
end
