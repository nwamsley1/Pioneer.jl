using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "PioneerPlots package boundary" begin
    plots_project = joinpath(REPO_ROOT, "packages", "PioneerPlots", "Project.toml")
    plots_module = joinpath(REPO_ROOT, "packages", "PioneerPlots", "src", "PioneerPlots.jl")
    plot_specs = joinpath(REPO_ROOT, "src", "shared", "plot_specs.jl")
    common_module = joinpath(REPO_ROOT, "packages", "PioneerCommon", "src", "PioneerCommon.jl")
    search_project = joinpath(REPO_ROOT, "packages", "PioneerSearch", "Project.toml")
    search_module = joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl")
    search_bootstrap = joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "bootstrap.jl")
    parameter_tuning_types = joinpath(
        REPO_ROOT,
        "src",
        "Routines",
        "SearchDIA",
        "SearchMethods",
        "ParameterTuningSearch",
        "types.jl",
    )
    first_pass_search = joinpath(
        REPO_ROOT,
        "src",
        "Routines",
        "SearchDIA",
        "SearchMethods",
        "FirstPassSearch",
        "FirstPassSearch.jl",
    )
    nce_tuning_search = joinpath(
        REPO_ROOT,
        "src",
        "Routines",
        "SearchDIA",
        "SearchMethods",
        "NceTuningSearch",
        "NceTuningSearch.jl",
    )
    quad_tuning_search = joinpath(
        REPO_ROOT,
        "src",
        "Routines",
        "SearchDIA",
        "SearchMethods",
        "QuadTuningSearch",
        "QuadTuningSearch.jl",
    )
    predict_project = joinpath(REPO_ROOT, "packages", "PioneerPredict", "Project.toml")
    predict_module = joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "PioneerPredict.jl")
    predict_bootstrap = joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "bootstrap.jl")

    @test isfile(plots_project)
    @test isfile(plots_module)
    @test isfile(plot_specs)

    @test occursin("module PioneerPlots", read(plots_module, String))
    @test occursin("struct RtAlignmentPlotSpec", read(plot_specs, String))
    @test occursin("struct MassErrorPlotSpec", read(plot_specs, String))
    @test occursin("export RtAlignmentPlotSpec", read(common_module, String))
    @test occursin("export MassErrorPlotSpec", read(common_module, String))

    @test occursin("PioneerPlots", read(search_project, String))
    search_content = read(search_module, String)
    @test !occursin("\nusing PioneerPlots\n", search_content)
    @test occursin("Base.require(PIONEER_PLOTS_PKGID)", search_content)
    @test occursin("get_pioneer_plots_module()", search_content)
    @test !occursin("\nusing Measures\n", search_content)
    @test !occursin("\nusing Plots", search_content)
    @test !occursin("\nusing StatsPlots", search_content)
    @test !occursin("\nusing LaTeXStrings\n", search_content)
    @test occursin("const Plots = LazyPlotsProxy", search_content)

    bootstrap_content = read(search_bootstrap, String)
    @test !occursin("pdfUtils.jl", bootstrap_content)
    @test !occursin("qcPlots.jl", bootstrap_content)
    @test !occursin("plotRTAlignment.jl", bootstrap_content)

    parameter_tuning_types_content = read(parameter_tuning_types, String)
    @test !occursin("rt_plots::Vector{Plots.Plot}", parameter_tuning_types_content)
    @test !occursin("mass_plots::Vector{Plots.Plot}", parameter_tuning_types_content)
    @test occursin("rt_plot_specs::Vector{RtAlignmentPlotSpec}", parameter_tuning_types_content)
    @test occursin("mass_plot_specs::Vector{MassErrorPlotSpec}", parameter_tuning_types_content)

    first_pass_content = read(first_pass_search, String)
    @test !occursin("ms1_mass_plots::Vector{Plots.Plot}", first_pass_content)
    @test occursin("ms1_mass_plot_specs::Vector{MassErrorPlotSpec}", first_pass_content)
    @test occursin("save_multipage_pdf(results.ms1_mass_plot_specs, output_path)", first_pass_content)

    nce_content = read(nce_tuning_search, String)
    @test !occursin("nce_plots::Vector{Plots.Plot}", nce_content)
    @test occursin("nce_plot_specs::Vector{MultiSeriesPlotSpec}", nce_content)
    @test occursin("save_multipage_pdf(results.nce_plot_specs, output_path)", nce_content)

    quad_content = read(quad_tuning_search, String)
    @test !occursin("quad_model_plots::Vector{Plots.Plot}", quad_content)
    @test !occursin("quad_data_plots::Vector{Plots.Plot}", quad_content)
    @test occursin("quad_model_plot_specs::Vector{MultiSeriesPlotSpec}", quad_content)
    @test occursin("quad_data_plot_specs::Vector{MultiSeriesPlotSpec}", quad_content)
    @test occursin("save_multipage_pdf(results.quad_model_plot_specs, models_path)", quad_content)
    @test occursin("save_multipage_pdf(results.quad_data_plot_specs, data_path)", quad_content)

    @test !occursin("\nPlots =", read(predict_project, String))
    @test !occursin("\nStatsPlots =", read(predict_project, String))
    @test !occursin("\nMeasures =", read(predict_project, String))
    @test !occursin("\nLaTeXStrings =", read(predict_project, String))
    @test !occursin("using Plots", read(predict_module, String))
    @test !occursin("using StatsPlots", read(predict_module, String))
    @test !occursin("using Measures", read(predict_module, String))
    @test !occursin("using LaTeXStrings", read(predict_module, String))
    @test !occursin("pdfUtils.jl", read(predict_bootstrap, String))
end
