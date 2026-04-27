using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "Required QC plot failures are fatal" begin
    parameter_tuning = read(
        joinpath(
            REPO_ROOT,
            "src",
            "Routines",
            "SearchDIA",
            "SearchMethods",
            "ParameterTuningSearch",
            "ParameterTuningSearch.jl",
        ),
        String,
    )
    nce_tuning = read(
        joinpath(
            REPO_ROOT,
            "src",
            "Routines",
            "SearchDIA",
            "SearchMethods",
            "NceTuningSearch",
            "NceTuningSearch.jl",
        ),
        String,
    )

    @test !occursin("Failed to generate plots for file", parameter_tuning)
    @test !occursin("Failed to merge QC plots", parameter_tuning)
    @test !occursin("catch\n        nothing\n    end\n    # Could add NCE model statistics or plots here", nce_tuning)
    @test occursin("save_multipage_pdf(results.rt_plot_specs, rt_combined_path)", parameter_tuning)
    @test occursin("save_multipage_pdf(results.mass_plot_specs, mass_combined_path)", parameter_tuning)
    @test occursin("save_multipage_pdf(results.nce_plot_specs, output_path)", nce_tuning)
end
