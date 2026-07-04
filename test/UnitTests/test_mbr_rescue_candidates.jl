using Test
using DataFrames
using Arrow
using Pioneer

@testset "MBR rescue candidate helpers" begin
    @testset "main-search PEP partition bounds MBR rescue candidates" begin
        @test Pioneer.MAIN_MBR_RESCUE_PEP_MAX == 0.99f0

        peps = Float32[0.0, 0.9, 0.9001, 0.98, 0.99, 1.0]
        part = Pioneer._mainsearch_pep_partition(peps, 0.9f0, 0.99f0)

        @test part.keep_mask == Bool[true, true, false, false, false, false]
        @test part.rescue_mask == Bool[false, false, true, true, true, false]
        @test part.discard_mask == Bool[false, false, false, false, false, true]
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
            @test Pioneer.mbr_rescue_recovered_path(rescue_path) ==
                  joinpath(rescue_dir, "run1_fold0.recovered.arrow")
        end
    end

    @testset "recovered rescue files contain only accepted rescue rows" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "run1_fold0.arrow")
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                ms_file_idx = UInt32[1, 1],
                cv_fold = UInt8[0, 0],
                target = Bool[true, true],
                rt = Float32[1, 2],
                main_pep = Float32[0.93, 0.94],
            ))
            slim_df = DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                _mbr_rescue_file_idx = UInt32[1, 1],
                _mbr_rescue_row_idx = UInt32[1, 2],
                trace_prob_prepass = Float32[0.07, 0.08],
                trace_prob_infold = Float32[0.0, 0.0],
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
            @test :prec_prob ∉ propertynames(recovered)
            @test recovered.mbr_recovered == Bool[true]
            @test recovered.MBR_transfer_candidate == Bool[true]
            @test recovered.mbr_target_decoy_prob == Float32[0.88]
            @test recovered.ftr_qval_true == Float32[0.005]
            @test recovered.ftr_pep_true == Float32[0.006]

            empty_result = Pioneer.write_mbr_rescue_recovered_files!(slim_df[1:0, :], [rescue_path])
            @test empty_result.n_files == 0
            @test empty_result.n_rows == 0
            @test !isfile(recovered_path)
        end
    end

    @testset "normal/rescue fold unions normalize optional MBR columns" begin
        df = vcat(
            DataFrame(
                precursor_idx = UInt32[1],
                mbr_rescue_candidate = Bool[false],
            ),
            DataFrame(
                precursor_idx = UInt32[2],
                mbr_recovered = Bool[true],
                MBR_transfer_candidate = Bool[true],
                mbr_target_decoy_prob = Float32[0.8],
                ftr_qval_true = Float32[0.005],
                ftr_pep_true = Float32[0.006],
            );
            cols = :union,
        )

        Pioneer._normalize_optional_mbr_columns!(df)

        @test df.mbr_rescue_candidate == Bool[false, false]
        @test df.mbr_recovered == Bool[false, true]
        @test df.MBR_transfer_candidate == Bool[false, true]
        @test isnan(df.mbr_target_decoy_prob[1])
        @test df.mbr_target_decoy_prob[2] == Float32(0.8)
        @test isnan(df.ftr_qval_true[1])
        @test df.ftr_qval_true[2] == Float32(0.005)
        @test isnan(df.ftr_pep_true[1])
        @test df.ftr_pep_true[2] == Float32(0.006)
    end

    @testset "rescue slim rows use synthetic Pass-1 confidence from main-search PEP" begin
        mktempdir() do dir
            rescue_path = joinpath(dir, "run1_fold0.arrow")
            Pioneer.writeArrow(rescue_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                ms_file_idx = UInt32[1, 1],
                cv_fold = UInt8[0, 0],
                target = Bool[true, false],
                main_pep = Float32[0.93, 0.99],
            ))

            Pioneer._write_rescue_pass1_sidecars_from_main!([rescue_path])

            sidecar_columns = Dict{Symbol, Any}(
                :precursor_idx => UInt32[1, 2],
                :scan_idx => UInt32[10, 20],
            )
            for col in Pioneer._MBR_SIDECAR_OUT_COLS
                col in (:precursor_idx, :scan_idx) && continue
                sidecar_columns[col] = occursin("is_missing", String(col)) ?
                                       Bool[false, true] :
                                       Float32[0.5, -1.0]
            end
            Pioneer.writeArrow(rescue_path * Pioneer.MBR_SIDECAR_SUFFIX, DataFrame(sidecar_columns))

            slim = Pioneer.load_ftr_slim_dataframe(
                [rescue_path];
                mbr_rescue_candidate = true,
            )

            @test slim.trace_prob_prepass == Float32[0.07, 0.01]
            @test slim.trace_prob_infold == Float32[0.0, 0.0]
            @test slim.main_search_prob == Float32[0.07, 0.01]
            @test slim.mbr_rescue_candidate == Bool[true, true]
            @test slim._mbr_rescue_file_idx == UInt32[1, 1]
            @test slim._mbr_rescue_row_idx == UInt32[1, 2]
            @test slim.MBR_best_pair_prob_true == Float32[0.5, -1.0]
            @test slim.MBR_best_is_missing_false == Bool[false, true]
        end
    end
end
