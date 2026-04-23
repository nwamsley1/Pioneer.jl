using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const PACKAGE_FILES = Dict(
    "PioneerCommon" => (
        project = joinpath(REPO_ROOT, "packages", "PioneerCommon", "Project.toml"),
        module = joinpath(REPO_ROOT, "packages", "PioneerCommon", "src", "PioneerCommon.jl"),
    ),
    "PioneerPredict" => (
        project = joinpath(REPO_ROOT, "packages", "PioneerPredict", "Project.toml"),
        module = joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "PioneerPredict.jl"),
    ),
    "PioneerSearch" => (
        project = joinpath(REPO_ROOT, "packages", "PioneerSearch", "Project.toml"),
        module = joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl"),
    ),
)

@testset "Split packages own their runtime without using Pioneer" begin
    common_project = read(PACKAGE_FILES["PioneerCommon"].project, String)
    common_module = read(PACKAGE_FILES["PioneerCommon"].module, String)
    predict_project = read(PACKAGE_FILES["PioneerPredict"].project, String)
    predict_module = read(PACKAGE_FILES["PioneerPredict"].module, String)
    search_project = read(PACKAGE_FILES["PioneerSearch"].project, String)
    search_module = read(PACKAGE_FILES["PioneerSearch"].module, String)

    @test occursin("Dates", common_project)
    @test occursin("Interpolations", common_project)
    @test occursin("asset_path", common_module)
    @test occursin("export get_pioneer_version", common_module)
    @test occursin("export asset_path", common_module)
    @test occursin("export @user_info", common_module)
    @test occursin("export InterpolationTypeAlias", common_module)

    @test !occursin("Pioneer =", predict_project)
    @test !occursin("using Pioneer", predict_module)
    @test occursin("using PioneerCommon", predict_module)
    @test occursin("const Pioneer = PioneerCommon", predict_module)
    @test occursin("include(joinpath(@__DIR__, \"bootstrap.jl\"))", predict_module)
    @test occursin("include(joinpath(@__DIR__, \"owned\", \"BuildSpecLib.jl\"))", predict_module)

    @test !occursin("Pioneer =", search_project)
    @test !occursin("using Pioneer", search_module)
    @test occursin("using PioneerCommon", search_module)
    @test occursin("const Pioneer = PioneerCommon", search_module)
    @test occursin("include(joinpath(@__DIR__, \"bootstrap.jl\"))", search_module)
    @test occursin("include(joinpath(@__DIR__, \"owned\", \"GenerateParams.jl\"))", search_module)
    @test occursin("include(joinpath(@__DIR__, \"owned\", \"SearchDIA.jl\"))", search_module)
end
