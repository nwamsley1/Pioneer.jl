using Test
using DataFrames
using Pioneer

using Pioneer: IrtLinearRefinement, _precursor_token_counts, apply_iteration_postprocess!
using Pioneer: DataFramePSMContainer, LibraryPrecursors, default_scoring_config, setup_scoring_workspace

struct MockPrecursors <: LibraryPrecursors
    sequence::Vector{String}
    structural_mods::Vector{Union{Missing, String}}
end

Pioneer.getSequence(precursors::MockPrecursors) = precursors.sequence
Pioneer.getStructuralMods(precursors::MockPrecursors) = precursors.structural_mods

@testset "Precursor iRT Refinement" begin
    @testset "Modified Residues Get Unique Tokens" begin
        counts = _precursor_token_counts(
            "ACDM",
            "(1,n,Acetyl)(3,D,Phospho)(4,M,Unimod:35)(4,c,Amidated)"
        )

        @test counts["A"] == 1.0f0
        @test counts["C"] == 1.0f0
        @test counts["D|Phospho"] == 1.0f0
        @test counts["M|Unimod:35"] == 1.0f0
        @test counts["n|Acetyl"] == 1.0f0
        @test counts["c|Amidated"] == 1.0f0
    end

    @testset "Iteration One Updates Held-Out Folds Per Run" begin
        precursors = MockPrecursors(
            [
                "AAAA", "CCCC", "AAAA", "CCCC", "DDDD", "DDDD",
                "AAAA", "CCCC", "AAAA", "CCCC", "DDDD", "DDDD"
            ],
            fill(missing, 12)
        )

        psms = DataFrame(
            target = Bool[
                true, true, false, true, true, false,
                true, true, false, true, true, false
            ],
            cv_fold = UInt8[
                0, 0, 0, 1, 1, 1,
                0, 0, 0, 1, 1, 1
            ],
            precursor_idx = UInt32[
                1, 2, 5, 3, 4, 6,
                7, 8, 11, 9, 10, 12
            ],
            ms_file_idx = UInt32[
                1, 1, 1, 1, 1, 1,
                2, 2, 2, 2, 2, 2
            ],
            irt_pred = Float32[
                10, 10, 10, 10, 10, 10,
                20, 20, 20, 20, 20, 20
            ],
            isotopes_captured = fill((Int8(0), Int8(0)), 12),
            trace_prob = Float32[
                0.99, 0.98, 0.01, 0.99, 0.98, 0.01,
                0.99, 0.98, 0.01, 0.99, 0.98, 0.01
            ],
            irt_obs = Float32[
                4, 8, 6, 4, 8, 6,
                18, 14, 16, 18, 14, 16
            ],
            irt_error = fill(0.0f0, 12)
        )
        psms[!, :irt_error] = abs.(psms.irt_obs .- psms.irt_pred)

        workspace = setup_scoring_workspace(
            DataFramePSMContainer(psms, Val(:unsafe)),
            default_scoring_config()
        )
        strategy = IrtLinearRefinement(precursors)

        apply_iteration_postprocess!(strategy, workspace, 1, 4)

        @test psms.irt_pred[1] ≈ 4.0f0 atol = 1e-4
        @test psms.irt_pred[2] ≈ 8.0f0 atol = 1e-4
        @test psms.irt_pred[4] ≈ 4.0f0 atol = 1e-4
        @test psms.irt_pred[5] ≈ 8.0f0 atol = 1e-4
        @test psms.irt_pred[7] ≈ 18.0f0 atol = 1e-4
        @test psms.irt_pred[8] ≈ 14.0f0 atol = 1e-4
        @test psms.irt_pred[10] ≈ 18.0f0 atol = 1e-4
        @test psms.irt_pred[11] ≈ 14.0f0 atol = 1e-4
        @test psms.irt_error[1] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[2] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[7] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[8] ≈ 0.0f0 atol = 1e-4
    end

    @testset "Correction Model Uses Current Predicted iRT Per Run And Writes QC Plots" begin
        mktempdir() do qc_dir
            precursors = MockPrecursors(
                fill("AAAA", 12),
                fill(missing, 12)
            )

            psms = DataFrame(
                target = Bool[
                    true, true, false, true, true, false,
                    true, true, false, true, true, false
                ],
                cv_fold = UInt8[
                    0, 0, 0, 1, 1, 1,
                    0, 0, 0, 1, 1, 1
                ],
                precursor_idx = UInt32[1, 2, 5, 3, 4, 6, 7, 8, 11, 9, 10, 12],
                ms_file_idx = UInt32[
                    1, 1, 1, 1, 1, 1,
                    2, 2, 2, 2, 2, 2
                ],
                irt_pred = Float32[
                    12, 16, 14, 10, 14, 12,
                    22, 26, 24, 20, 24, 22
                ],
                isotopes_captured = fill((Int8(0), Int8(0)), 12),
                trace_prob = Float32[
                    0.99, 0.98, 0.01, 0.99, 0.98, 0.01,
                    0.99, 0.98, 0.01, 0.99, 0.98, 0.01
                ],
                irt_obs = Float32[
                    8, 8, 8, 8, 8, 8,
                    18, 18, 18, 18, 18, 18
                ],
                irt_error = fill(0.0f0, 12)
            )
            psms[!, :irt_error] = abs.(psms.irt_obs .- psms.irt_pred)

            workspace = setup_scoring_workspace(
                DataFramePSMContainer(psms, Val(:unsafe)),
                default_scoring_config()
            )
            strategy = IrtLinearRefinement(precursors; qc_plot_dir = qc_dir, qc_run_id = UInt32(1))

            apply_iteration_postprocess!(strategy, workspace, 1, 4)

            @test psms.irt_pred[1] ≈ 8.0f0 atol = 1e-4
            @test psms.irt_pred[2] ≈ 8.0f0 atol = 1e-4
            @test psms.irt_pred[4] ≈ 8.0f0 atol = 1e-4
            @test psms.irt_pred[5] ≈ 8.0f0 atol = 1e-4
            @test psms.irt_pred[7] ≈ 18.0f0 atol = 1e-4
            @test psms.irt_pred[8] ≈ 18.0f0 atol = 1e-4
            @test psms.irt_pred[10] ≈ 18.0f0 atol = 1e-4
            @test psms.irt_pred[11] ≈ 18.0f0 atol = 1e-4
            @test psms.irt_error[4] ≈ 0.0f0 atol = 1e-4
            @test psms.irt_error[5] ≈ 0.0f0 atol = 1e-4
            @test psms.irt_error[10] ≈ 0.0f0 atol = 1e-4
            @test psms.irt_error[11] ≈ 0.0f0 atol = 1e-4
            @test isfile(joinpath(qc_dir, "irt_error_refinement_run_1.png"))
            @test !isfile(joinpath(qc_dir, "irt_error_refinement_run_2.png"))
        end
    end
end
