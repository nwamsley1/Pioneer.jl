using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "macOS combined app PackageCompiler experiment" begin
    cli_project = joinpath(REPO_ROOT, "packages", "PioneerCLI", "Project.toml")
    cli_module = joinpath(REPO_ROOT, "packages", "PioneerCLI", "src", "PioneerCLI.jl")
    macos_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)
    unix_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer"), String)

    @test isfile(cli_project)
    @test isfile(cli_module)

    module_text = isfile(cli_module) ? read(cli_module, String) : ""
    @test occursin("module PioneerCLI", module_text)
    @test !occursin("__precompile__(false)", module_text)
    @test !occursin("import PioneerSearch", module_text)
    @test !occursin("import PioneerPredict", module_text)
    @test occursin("_pioneer_search()", module_text)
    @test occursin("_pioneer_predict()", module_text)
    @test occursin("_pioneer_params()", module_text)
    @test occursin("_pioneer_convert()", module_text)
    @test occursin("Base.invokelatest", module_text)

    @test occursin("\"packages/PioneerCLI\"", macos_workflow)
    @test occursin("joinpath(apps_root, \"runtime\")", macos_workflow)
    @test occursin("\"SearchDIA\"=>\"main_SearchDIA\"", macos_workflow)
    @test !occursin("create_app(\n              \"packages/PioneerParams\"", macos_workflow)
    @test !occursin("create_app(\n              \"packages/PioneerSearch\"", macos_workflow)
    @test occursin("apps/runtime/bin/\$executable_name", unix_launcher)
    @test occursin("apps/runtime/share/julia", unix_launcher)
end
