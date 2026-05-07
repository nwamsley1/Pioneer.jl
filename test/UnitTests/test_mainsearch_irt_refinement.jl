using DataFrames
using Test

struct MockMainSearchIrtPrecursors <: Pioneer.LibraryPrecursors
    sequence::Vector{String}
    structural_mods::Vector{Union{Missing, String}}
end

Pioneer.getSequence(precursors::MockMainSearchIrtPrecursors) = precursors.sequence
Pioneer.getStructuralMods(precursors::MockMainSearchIrtPrecursors) = precursors.structural_mods

@testset "MainSearch per-run iRT refinement" begin
    @test isdefined(Pioneer, :_mainsearch_precursor_token_counts)
    @test isdefined(Pioneer, :_fit_mainsearch_irt_refinement_model)
    @test isdefined(Pioneer, :_predict_mainsearch_irt_refinement)
    @test isdefined(Pioneer, :_refine_mainsearch_irt_after_first_iteration!)
    @test isdefined(Pioneer, :MAINSEARCH_ENABLE_IRT_REFINEMENT)
    @test !Pioneer.MAINSEARCH_ENABLE_IRT_REFINEMENT

    @testset "modified residues get distinct iRT tokens" begin
        counts = Pioneer._mainsearch_precursor_token_counts(
            "ACDM",
            "(1,n,Acetyl)(3,D,Phospho)(4,M,Unimod:35)(4,c,Amidated)",
        )

        @test counts["A"] == 1.0f0
        @test counts["C"] == 1.0f0
        @test counts["D|Phospho"] == 1.0f0
        @test counts["M|Unimod:35"] == 1.0f0
        @test counts["n|Acetyl"] == 1.0f0
        @test counts["c|Amidated"] == 1.0f0
        @test counts["NTERM|A"] == 1.0f0
        @test counts["CTERM|M|Unimod:35"] == 1.0f0
    end

    @testset "quadratic iRT basis fits nonlinear correction" begin
        precursors = MockMainSearchIrtPrecursors(
            fill("AAAA", 4),
            fill(missing, 4),
        )
        cache = Dict{UInt32, Dict{String, Float32}}()
        model = Pioneer._fit_mainsearch_irt_refinement_model(
            precursors,
            UInt32[1, 2, 3],
            Float32[1, 2, 3],
            Float32[1, 4, 9],
            cache,
        )

        @test !isnothing(model)
        correction = Pioneer._predict_mainsearch_irt_refinement(
            precursors,
            model,
            UInt32(4),
            4.0f0,
            cache,
        )
        @test correction ≈ 16.0f0 atol = 1e-3
    end

    @testset "training precursors use first predicted iRT and average observed iRT" begin
        precursor_ids, irt_pred_inputs, irt_corrections =
            Pioneer._mainsearch_passing_irt_training_precursors(
                UInt32[1, 1, 1, 2],
                Bool[true, true, true, true],
                Float32[0.99, 0.98, 0.97, 0.96],
                Float32[10, 14, 20, 20],
                Float32[12, 19, 25, 25],
            )

        @test precursor_ids == UInt32[1, 2]
        @test irt_pred_inputs[1] == 10f0
        @test irt_pred_inputs[2] == 20f0
        @test irt_corrections[1] ≈ Float32(((12 + 19 + 25) / 3) - 10)
        @test irt_corrections[2] == 5f0
    end

    @testset "first-iteration scores refine held-out folds only within a run" begin
        precursors = MockMainSearchIrtPrecursors(
            ["AAAA", "CCCC", "AAAA", "CCCC", "DDDD", "DDDD"],
            fill(missing, 6),
        )
        psms = DataFrame(
            target = Bool[true, true, false, true, true, false],
            cv_fold = UInt8[0, 0, 0, 1, 1, 1],
            precursor_idx = UInt32[1, 2, 5, 3, 4, 6],
            irt_pred = Float32[10, 10, 10, 10, 10, 10],
            irt_obs = Float32[4, 8, 6, 4, 8, 6],
            irt_error = Float32[6, 2, 4, 6, 2, 4],
        )
        scores = Float64[0.99, 0.98, 0.01, 0.99, 0.98, 0.01]

        summary = Pioneer._refine_mainsearch_irt_after_first_iteration!(
            psms,
            precursors,
            scores;
            label = "unit",
        )

        @test summary.refined
        @test summary.models_fit == 2
        @test summary.training_precursors == 4
        @test psms.irt_pred[1] ≈ 4.0f0 atol = 1e-4
        @test psms.irt_pred[2] ≈ 8.0f0 atol = 1e-4
        @test psms.irt_pred[4] ≈ 4.0f0 atol = 1e-4
        @test psms.irt_pred[5] ≈ 8.0f0 atol = 1e-4
        @test psms.irt_error[1] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[2] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[4] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[5] ≈ 0.0f0 atol = 1e-4
        @test psms.irt_error[3] < 4.0f0
        @test psms.irt_error[6] < 4.0f0
    end
end
