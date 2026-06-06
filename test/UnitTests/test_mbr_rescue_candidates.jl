using Test
using DataFrames
using Arrow
using Pioneer

@testset "MBR rescue candidate helpers" begin
    @testset "main-search PEP partition bounds MBR rescue candidates" begin
        peps = Float32[0.0, 0.9, 0.9001, 0.95, 0.9501, 1.0]

        part = Pioneer._mainsearch_pep_partition(peps, 0.9f0, 0.95f0)

        @test part.keep_mask == Bool[true, true, false, false, false, false]
        @test part.rescue_mask == Bool[false, false, true, true, false, false]
        @test part.discard_mask == Bool[false, false, false, false, true, true]
    end

    @testset "main-search PEP thresholds are configurable" begin
        user = (main_pep_filter_threshold = 0.8,
                mbr_rescue_pep_threshold = 0.98)

        thresholds = Pioneer._resolve_mainsearch_pep_thresholds(user, NamedTuple())

        @test thresholds.pep_filter_threshold == 0.8f0
        @test thresholds.mbr_rescue_pep_max == 0.98f0
    end

    @testset "rescue pass1 sidecars use per-run 1-PEP confidence" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "run1_fold0.arrow")
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[1, 2, 3],
                scan_idx = UInt32[10, 20, 30],
                main_pep = Float32[0.91, 0.95, 1.20],
            ))

            result = Pioneer.write_mbr_rescue_main_pep_pass1_sidecars!([rescue_path])
            side = Arrow.Table(rescue_path * Pioneer.PASS1_SIDECAR_SUFFIX)

            @test result.n_rows == 3
            @test collect(side.precursor_idx) == UInt32[1, 2, 3]
            @test collect(side.scan_idx) == UInt32[10, 20, 30]
            @test collect(side.trace_prob_prepass) == Float32[0.09, 0.05, 0.0]
            @test collect(side.trace_prob_infold) == Float32[0.09, 0.05, 0.0]
        end
    end

    @testset "FTR slim loader exposes per-run 1-PEP confidence" begin
        @test :mbr_recipient_main_confidence in Pioneer.FTR_FEATURES_F_TRUE
        @test :trace_prob_infold ∉ Pioneer.FTR_FEATURES_F_TRUE

        mktempdir() do dir
            path = joinpath(dir, "run1_fold0.arrow")
            Pioneer.writeArrow(path, DataFrame(
                precursor_idx = UInt32[1],
                scan_idx = UInt32[10],
                ms_file_idx = UInt32[1],
                cv_fold = UInt8[0],
                target = Bool[true],
                main_pep = Float32[0.93],
                mbr_rescue_candidate = Bool[true],
            ))
            Pioneer.writeArrow(path * Pioneer.PASS1_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = UInt32[1],
                scan_idx = UInt32[10],
                trace_prob_prepass = Float32[0.07],
                trace_prob_infold = Float32[0.07],
            ))
            Pioneer.writeArrow(path * Pioneer.MBR_SIDECAR_SUFFIX, DataFrame(
                precursor_idx = UInt32[1],
                scan_idx = UInt32[10],
                MBR_max_pair_prob_true = Float32[0.99],
                MBR_max_pair_prob_false = Float32[0.98],
                MBR_log2_weight_ratio_true = Float32[0.1],
                MBR_log2_weight_ratio_false = Float32[0.2],
                MBR_log2_explained_ratio_true = Float32[0.3],
                MBR_log2_explained_ratio_false = Float32[0.4],
                MBR_best_irt_diff_true = Float32[0.5],
                MBR_best_irt_diff_false = Float32[0.6],
                MBR_is_missing_true = Bool[false],
                MBR_is_missing_false = Bool[false],
                MBR_log_by_diff_true = Float32[0.7],
                MBR_log_by_diff_false = Float32[0.8],
                MBR_best_rt_diff_true = Float32[0.9],
                MBR_best_rt_diff_false = Float32[1.0],
            ))

            slim = Pioneer.load_ftr_slim_dataframe([path])

            @test slim.mbr_recipient_main_confidence == Float32[0.07]
            @test slim.mbr_rescue_candidate == Bool[true]
        end
    end

    @testset "precursor scoring merge tolerates slim recovered rescue rows" begin
        normal = DataFrame(
            precursor_idx = UInt32[1],
            scan_idx = UInt32[10],
            wide_ms1_m0_candidate_fraction = Float32[0.5],
        )
        rescue = DataFrame(
            precursor_idx = UInt32[2],
            scan_idx = UInt32[20],
            mbr_recovered = Bool[true],
        )

        combined = Pioneer._combine_precursor_scoring_fold_dfs([normal, rescue])

        @test nrow(combined) == 2
        @test collect(combined.precursor_idx) == UInt32[1, 2]
        @test ismissing(combined.wide_ms1_m0_candidate_fraction[2])
        @test ismissing(combined.mbr_recovered[1])
        @test combined.mbr_recovered[2] == true
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

    @testset "candidate gate uses normal rows for donor threshold and admits rescue rows" begin
        psms = DataFrame(
            trace_prob_prepass = Float32[0.99, 0.20, 0.05],
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

    @testset "recovery sidecars preserve duplicate normal and rescue fold keys" begin
        mktempdir() do dir
            normal_path = joinpath(dir, "normal_fold0.arrow")
            rescue_path = joinpath(dir, "rescue_fold0.arrow")
            Pioneer.writeArrow(normal_path, DataFrame(
                precursor_idx = UInt32[1],
                scan_idx = UInt32[10],
                ms_file_idx = UInt32[1],
                cv_fold = UInt8[0],
            ))
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[2],
                scan_idx = UInt32[20],
                ms_file_idx = UInt32[1],
                cv_fold = UInt8[0],
            ))

            slim_df = DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                ms_file_idx = UInt32[1, 1],
                cv_fold = UInt8[0, 0],
                mbr_recovered = Bool[false, true],
                MBR_transfer_candidate = Bool[false, true],
                ftr_qval_true = Float32[NaN32, 0.005],
                ftr_pep_true = Float32[NaN32, 0.006],
            )

            Pioneer.write_recovery_sidecars(slim_df, [normal_path, rescue_path])

            normal_side = Arrow.Table(normal_path * Pioneer.RECOVERY_SIDECAR_SUFFIX)
            rescue_side = Arrow.Table(rescue_path * Pioneer.RECOVERY_SIDECAR_SUFFIX)
            @test collect(normal_side.precursor_idx) == UInt32[1]
            @test collect(normal_side.mbr_recovered) == Bool[false]
            @test collect(rescue_side.precursor_idx) == UInt32[2]
            @test collect(rescue_side.mbr_recovered) == Bool[true]
        end
    end
end
