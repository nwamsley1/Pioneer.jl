using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const PACKAGE_FILES = Dict(
    "PioneerCommon" => (
        joinpath(REPO_ROOT, "packages", "PioneerCommon", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerCommon", "src", "PioneerCommon.jl"),
    ),
    "PioneerParams" => (
        joinpath(REPO_ROOT, "packages", "PioneerParams", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerParams", "src", "PioneerParams.jl"),
    ),
    "PioneerPlots" => (
        joinpath(REPO_ROOT, "packages", "PioneerPlots", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerPlots", "src", "PioneerPlots.jl"),
    ),
    "PioneerPredict" => (
        joinpath(REPO_ROOT, "packages", "PioneerPredict", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "PioneerPredict.jl"),
    ),
    "PioneerSearch" => (
        joinpath(REPO_ROOT, "packages", "PioneerSearch", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl"),
    ),
)

@testset "Split packages own their runtime without using Pioneer" begin
    common_project = read(PACKAGE_FILES["PioneerCommon"][1], String)
    common_module = read(PACKAGE_FILES["PioneerCommon"][2], String)
    params_project = read(PACKAGE_FILES["PioneerParams"][1], String)
    params_module = read(PACKAGE_FILES["PioneerParams"][2], String)
    params_generate = read(joinpath(REPO_ROOT, "packages", "PioneerParams", "src", "routines", "GenerateParams.jl"), String)
    plots_project = read(PACKAGE_FILES["PioneerPlots"][1], String)
    plots_module = read(PACKAGE_FILES["PioneerPlots"][2], String)
    predict_project = read(PACKAGE_FILES["PioneerPredict"][1], String)
    predict_module = read(PACKAGE_FILES["PioneerPredict"][2], String)
    search_project = read(PACKAGE_FILES["PioneerSearch"][1], String)
    search_module = read(PACKAGE_FILES["PioneerSearch"][2], String)
    search_load_spectral_library = read(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "ParseInputs", "loadSpectralLibrary.jl"), String)

    @test occursin("Dates", common_project)
    @test occursin("Interpolations", common_project)
    @test occursin("asset_path", common_module)
    @test occursin("export get_pioneer_version", common_module)
    @test occursin("export asset_path", common_module)
    @test occursin("export @user_info", common_module)
    @test occursin("export InterpolationTypeAlias", common_module)

    @test occursin("PioneerCommon", params_project)
    @test occursin("using PioneerCommon", params_module)
    @test occursin("include(joinpath(@__DIR__, \"routines\", \"GenerateParams.jl\"))", params_module)
    @test occursin("isa(regex_codes, AbstractDict)", params_generate)

    @test occursin("Plots", plots_project)
    @test occursin("module PioneerPlots", plots_module)
    @test occursin("export qcPlots", plots_module)
    @test occursin("export save_multipage_pdf", plots_module)
    @test occursin("export plotRTAlign", plots_module)

    @test !occursin("Pioneer =", predict_project)
    @test !occursin("\nusing Pioneer\n", predict_module)
    @test occursin("using PioneerCommon", predict_module)
    @test occursin("using PioneerParams", predict_module)
    @test occursin("const Pioneer = PioneerCommon", predict_module)
    @test !occursin("packages\", \"PioneerSearch\", \"src\", \"routines\", \"GenerateParams.jl", predict_module)
    @test occursin("GetBuildLibParams(args...; kwargs...) = PioneerParams.GetBuildLibParams(args...; kwargs...)", predict_module)

    @test !occursin("Pioneer =", search_project)
    @test !occursin("\nPlots =", search_project)
    @test !occursin("\nStatsPlots =", search_project)
    @test !occursin("\nMeasures =", search_project)
    @test !occursin("\nLaTeXStrings =", search_project)
    @test !occursin("\nGR =", search_project)
    @test !occursin("\nusing Pioneer\n", search_module)
    @test occursin("using PioneerCommon", search_module)
    @test occursin("using PioneerParams", search_module)
    @test occursin("\nusing PioneerPlots\n", search_module)
    @test !occursin("Base.require(PIONEER_PLOTS_PKGID)", search_module)
    @test occursin("const Plots = PioneerPlots.Plots", search_module)
    @test !occursin("save_multipage_pdf(args...; kwargs...)", search_module)
    @test !occursin("qcPlots(args...; kwargs...)", search_module)
    @test occursin("const Runtime = PioneerCommon", search_module)
    @test occursin("const Pioneer = @__MODULE__", search_module)
    @test occursin("include(joinpath(@__DIR__, \"bootstrap.jl\"))", search_module)
    @test !occursin("include(joinpath(@__DIR__, \"routines\", \"GenerateParams.jl\"))", search_module)
    @test occursin("GetSearchParams(args...; kwargs...) = PioneerParams.GetSearchParams(args...; kwargs...)", search_module)
    @test occursin("GetParseSpecLibParams(args...; kwargs...) = PioneerParams.GetParseSpecLibParams(args...; kwargs...)", search_module)
    @test occursin("include(joinpath(@__DIR__, \"routines\", \"SearchDIA.jl\"))", search_module)
    @test occursin("JSON.parsefile(config_path, dicttype=Dict{String,Any})", search_load_spectral_library)
end
