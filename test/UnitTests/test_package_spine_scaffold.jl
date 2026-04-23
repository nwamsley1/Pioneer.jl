using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const EXPECTED_PACKAGE_FILES = Dict(
    "PioneerCommon" => (
        joinpath(REPO_ROOT, "packages", "PioneerCommon", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerCommon", "src", "PioneerCommon.jl")
    ),
    "PioneerParams" => (
        joinpath(REPO_ROOT, "packages", "PioneerParams", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerParams", "src", "PioneerParams.jl")
    ),
    "PioneerPlots" => (
        joinpath(REPO_ROOT, "packages", "PioneerPlots", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerPlots", "src", "PioneerPlots.jl")
    ),
    "PioneerPredict" => (
        joinpath(REPO_ROOT, "packages", "PioneerPredict", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "PioneerPredict.jl")
    ),
    "PioneerSearch" => (
        joinpath(REPO_ROOT, "packages", "PioneerSearch", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl")
    ),
    "PioneerSIMD" => (
        joinpath(REPO_ROOT, "packages", "PioneerSIMD", "Project.toml"),
        joinpath(REPO_ROOT, "packages", "PioneerSIMD", "src", "PioneerSIMD.jl")
    ),
    "PioneerParamsApp" => (
        joinpath(REPO_ROOT, "apps", "PioneerParamsApp", "Project.toml"),
        joinpath(REPO_ROOT, "apps", "PioneerParamsApp", "src", "PioneerParamsApp.jl")
    ),
    "PioneerPredictApp" => (
        joinpath(REPO_ROOT, "apps", "PioneerPredictApp", "Project.toml"),
        joinpath(REPO_ROOT, "apps", "PioneerPredictApp", "src", "PioneerPredictApp.jl")
    ),
    "PioneerSearchApp" => (
        joinpath(REPO_ROOT, "apps", "PioneerSearchApp", "Project.toml"),
        joinpath(REPO_ROOT, "apps", "PioneerSearchApp", "src", "PioneerSearchApp.jl")
    )
)

@testset "Package spine scaffold" begin
    for (module_name, (project_file, module_file)) in EXPECTED_PACKAGE_FILES
        @test isfile(project_file)
        @test isfile(module_file)
        @test occursin("module $module_name", read(module_file, String))
    end
end
