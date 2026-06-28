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
        @test Pioneer.MAIN_MBR_RESCUE_PEP_MAX == 0.98f0

        peps = Float32[0.0, 0.9, 0.9001, 0.98, 0.9801, 1.0]

        part = Pioneer._mainsearch_pep_partition(peps, 0.9f0, 0.98f0)

        @test part.keep_mask == Bool[true, true, false, false, false, false]
        @test part.rescue_mask == Bool[false, false, true, true, false, false]
        @test part.discard_mask == Bool[false, false, false, false, true, true]
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
                smoothed_cols...,
            ))
            donor_sqrt = ntuple(_ -> sqrt(1.0f0 / 8.0f0), 8)
            donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
                UInt32(1) => [Pioneer._MBRDonorEntry(
                    0.95f0, UInt32(1), 10.0f0, 2.0f0, 0.25f0, 8.8f0,
                    0.1f0, donor_sqrt, UInt32(2), false,
                )],
                UInt32(2) => [Pioneer._MBRDonorEntry(
                    0.94f0, UInt32(2), 12.0f0, 1.0f0, 0.5f0, 8.5f0,
                    0.3f0, donor_sqrt, UInt32(2), false,
                )],
            )
            pool = (pids = UInt32[2], irts = Float32[10])
            partner_pools = Pioneer._CounterfactualPartnerPools(
                Bool[true, false, true],
                UInt8[0, 0, 0],
                UInt8[1, 1, 1],
                Float32[10, 10, 30],
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}(),
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => pool),
                Dict{Int, Pioneer._IrtPool}(0 => Pioneer._empty_pool(), 1 => Pioneer._empty_pool()),
                Dict{Int, Pioneer._IrtPool}(0 => pool, 1 => Pioneer._empty_pool()),
                Pioneer._empty_pool(),
                pool,
            )
            fragment_keys = Pioneer._MBRFragmentAnnotationKeys(zeros(UInt16, 8 * 3))

            slim = Pioneer.load_mbr_rescue_candidate_slim_dataframe(
                [rescue_path],
                donor_dict,
                partner_pools,
                fragment_keys,
                0.90f0,
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
                smoothed_cols...,
            ))
            Pioneer.writeArrow(normal_path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = UInt32[1, 2, 3],
                scan_idx = UInt32[10, 20, 30],
                trace_prob_prepass = Float32[0.80, 0.95, 0.99],
                trace_prob_infold = Float32[0.81, 0.96, 0.995],
            ))
            donor_sqrt = ntuple(_ -> sqrt(1.0f0 / 8.0f0), 8)
            donor_dict = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}(
                UInt32(1) => [Pioneer._MBRDonorEntry(
                    0.99f0, UInt32(1), 10.0f0, 2.0f0, 0.25f0, 8.8f0,
                    0.1f0, donor_sqrt, UInt32(2), false,
                )],
                UInt32(2) => [Pioneer._MBRDonorEntry(
                    0.99f0, UInt32(2), 12.0f0, 1.0f0, 0.5f0, 8.5f0,
                    0.3f0, donor_sqrt, UInt32(2), false,
                )],
            )
            pool = (pids = UInt32[2], irts = Float32[10])
            partner_pools = Pioneer._CounterfactualPartnerPools(
                Bool[true, false, true],
                UInt8[0, 0, 0],
                UInt8[1, 1, 1],
                Float32[10, 10, 30],
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}(),
                Dict{Tuple{Int, Int}, Pioneer._IrtPool}((0, 1) => pool),
                Dict{Int, Pioneer._IrtPool}(0 => Pioneer._empty_pool(), 1 => Pioneer._empty_pool()),
                Dict{Int, Pioneer._IrtPool}(0 => pool, 1 => Pioneer._empty_pool()),
                Pioneer._empty_pool(),
                pool,
            )
            fragment_keys = Pioneer._MBRFragmentAnnotationKeys(zeros(UInt16, 8 * 3))

            result = Pioneer.load_normal_mbr_candidate_slim_dataframe(
                [normal_path],
                donor_dict,
                partner_pools,
                fragment_keys;
                q_thresh = 0.01f0,
            )

            @test !isfile(normal_path * ".mbr_sidecar.arrow")
            @test result.prob_thresh == 0.99f0
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
            @test recovered.ftr_qval_true == Float32[0.005]
            @test recovered.ftr_pep_true == Float32[0.006]

            empty_result = Pioneer.write_mbr_rescue_recovered_files!(slim_df[1:0, :], [rescue_path])
            @test empty_result.n_files == 0
            @test empty_result.n_rows == 0
            @test !isfile(recovered_path)
        end
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
            mbr_recovered = Bool[false, true, false],
        )

        initial_filtered = _apply_pipeline_to_df(psms, pipeline)

        @test initial_filtered.precursor_idx == UInt32[1, 2]
        @test initial_filtered.qval[2] > 0.1f0
    end

    @testset "MBR recovered rows are filtered by recalculated row q-values" begin
        recalc_qval_spline = linear_interpolation(
            Float32[0.0, 1.0],
            Float32[0.2, 0.0];
            extrapolation_bc = Interpolations.Flat(),
        )
        pipeline = Pioneer._precursor_scoring_recalculated_qvalue_filter_pipeline(
            recalc_qval_spline,
            0.1f0,
        )
        psms = DataFrame(
            precursor_idx = UInt32[1, 2],
            prec_prob = Float32[0.95, 0.05],
            mbr_recovered = Bool[false, true],
        )

        filtered = _apply_pipeline_to_df(psms, pipeline)

        @test filtered.precursor_idx == UInt32[1]
        @test filtered.qval ≈ Float32[0.01]
    end
end
