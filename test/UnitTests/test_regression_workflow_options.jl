using Pioneer: SCORING_SEMISUPERVISED_MAX_ITERATIONS

@testset "Regression workflow options" begin
    @testset "cross-run scoring semi-supervised iteration budget" begin
        @test SCORING_SEMISUPERVISED_MAX_ITERATIONS == 1
    end

    @testset "APMS Astral split regression sets are selectable" begin
        project_root = normpath(joinpath(@__DIR__, "..", ".."))
        workflow_paths = [
            joinpath(project_root, ".github", "workflows", "regression.yml"),
            joinpath(project_root, ".github", "workflows", "regression_slurm.yml"),
        ]
        expected_sets = [
            "APMS_Astral_1",
            "APMS_Astral_2",
            "APMS_Astral_5",
            "APMS_Astral_10",
            "APMS_Astral_20",
        ]

        for workflow_path in workflow_paths
            workflow_lines = strip.(split(read(workflow_path, String), '\n'))
            for regression_set in expected_sets
                @test "- $regression_set" in workflow_lines
            end
        end
    end
end
