using Arrow
using DataFrames
using Test
using Pioneer

function _mbr_staging_table(
    file_idx::UInt32,
    precursor_idx::Vector{UInt32},
    qval::Vector{Float32},
    global_qval::Vector{Float32},
    scores::Vector{Float32},
)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = precursor_idx,
        scan_idx = UInt32.(1000 .+ (1:n)),
        ms_file_idx = fill(file_idx, n),
        target = trues(n),
        cv_fold = UInt8.(precursor_idx .% UInt32(2)),
        qval = qval,
        global_qval = global_qval,
        trace_prob_prepass = scores,
        trace_prob_infold = scores,
        prec_prob = scores,
        irt_pred = Float32.(precursor_idx),
        rt = Float32.(precursor_idx),
    )
end

@testset "post-integration MBR never drops baseline IDs for missing peaks" begin
    psms = DataFrame(
        qval = Float32[0.005, 0.20, 0.20],
        global_qval = Float32[0.005, 0.005, 0.005],
        peak_area = Float32[0.0, 0.0, 10.0],
    )
    psms[!, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] =
        Float32[NaN, NaN, 12.0]

    @test Pioneer._select_postintegration_mbr_rows(psms, 0.01f0) ==
        Int[1, 3]
end

@testset "post-integration MBR removes internal evidence columns" begin
    mktempdir() do directory
        path = joinpath(directory, "integrated.arrow")
        table = DataFrame(
            precursor_idx = UInt32[1],
            scan_idx = UInt32[2],
            mbr_recovered = Bool[true],
            integration_start_scan = UInt32[10],
            integration_stop_scan = UInt32[12],
        )
        table[!, Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN] =
            [Float32[1, 2, 3]]
        Arrow.write(path, table)

        cleaned = DataFrame(Arrow.Table(path))
        Pioneer._drop_internal_mbr_columns!(cleaned)
        @test hasproperty(cleaned, :mbr_recovered)
        @test !hasproperty(
            cleaned,
            Pioneer.MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN,
        )
        @test !hasproperty(cleaned, :integration_start_scan)
        @test !hasproperty(cleaned, :integration_stop_scan)
    end
end

@testset "MBR recovery scores use the frozen pre-MBR q-value boundary" begin
    mktempdir() do directory
        path = joinpath(directory, "recoveries.arrow")
        Arrow.write(
            path,
            DataFrame(
                precursor_idx = UInt32[1, 2],
                scan_idx = UInt32[10, 20],
                prec_prob = Float32[0.90, 0.10],
                mbr_recovered = Bool[false, true],
                mbr_target_decoy_prob = Float32[NaN, 0.50],
                mbr_counterfactual_decoy_prob = Float32[NaN, 0.25],
            ),
        )
        frozen_qvalue = score ->
            Float32(score) >= 0.80f0 ? 0.005f0 : 0.10f0

        Pioneer._remap_mbr_scores!(
            Pioneer.PSMFileReference[
                Pioneer.PSMFileReference(path),
            ],
            joinpath(directory, "unused.arrow");
            q_value_threshold = 0.01f0,
            min_pep_points_per_bin = 10,
            fdr_scale_factor = 1.0f0,
            pre_mbr_qval_spline = frozen_qvalue,
        )
        remapped = DataFrame(Arrow.Table(path))
        @test remapped.prec_prob[1] == 0.90f0
        @test remapped.prec_prob[2] < 0.80f0
        @test remapped.prec_prob[2] ≈ 0.40f0 atol = 1.0f-4
        @test remapped.mbr_counterfactual_decoy_prec_prob[1] |> isnan
        @test remapped.mbr_counterfactual_decoy_prec_prob[2] ≈
            0.20f0 atol = 1.0f-4
    end
end

@testset "post-integration MBR staging preserves baseline and eligible candidates" begin
    mktempdir() do directory
        annotated_directory = joinpath(directory, "annotated")
        donor_directory = joinpath(directory, "donors")
        output_directory = joinpath(directory, "passing")
        mkpath(annotated_directory)
        mkpath(donor_directory)

        annotated_run1 = joinpath(annotated_directory, "run1.arrow")
        annotated_run2 = joinpath(annotated_directory, "run2.arrow")
        donor_run1 = joinpath(donor_directory, "run1.arrow")
        donor_run2 = joinpath(donor_directory, "run2.arrow")

        Arrow.write(
            annotated_run1,
            _mbr_staging_table(
                UInt32(1),
                UInt32[1, 2, 3],
                Float32[0.005, 0.20, 0.20],
                Float32[0.005, 0.005, 0.20],
                Float32[0.95, 0.60, 0.50],
            ),
        )
        Arrow.write(
            annotated_run2,
            _mbr_staging_table(
                UInt32(2),
                UInt32[2],
                Float32[0.005],
                Float32[0.005],
                Float32[0.96],
            ),
        )
        Arrow.write(
            donor_run1,
            _mbr_staging_table(
                UInt32(1),
                UInt32[1],
                Float32[0.005],
                Float32[0.005],
                Float32[0.95],
            ),
        )
        Arrow.write(
            donor_run2,
            _mbr_staging_table(
                UInt32(2),
                UInt32[2],
                Float32[0.005],
                Float32[0.005],
                Float32[0.96],
            ),
        )

        result = Pioneer.prepare_postintegration_mbr!(
            Pioneer.PSMFileReference[
                Pioneer.PSMFileReference(annotated_run1),
                Pioneer.PSMFileReference(annotated_run2),
            ],
            Pioneer.PSMFileReference[
                Pioneer.PSMFileReference(donor_run1),
                Pioneer.PSMFileReference(donor_run2),
            ],
            output_directory;
            q_value_threshold = 0.01f0,
        )

        @test result.n_candidates == 1
        @test result.n_rows == 3
        @test length(result.integration_refs) == 2

        staged_run1 = DataFrame(Arrow.Table(
            joinpath(output_directory, "run1.arrow"),
        ))
        @test staged_run1.precursor_idx == UInt32[1, 2]
        staged_pass1 = DataFrame(Arrow.Table(
            joinpath(output_directory, "run1.arrow") *
            Pioneer.PASS1_SIDECAR_SUFFIX,
        ))
        @test staged_pass1.precursor_idx == staged_run1.precursor_idx
        @test staged_pass1.scan_idx == staged_run1.scan_idx
    end
end
