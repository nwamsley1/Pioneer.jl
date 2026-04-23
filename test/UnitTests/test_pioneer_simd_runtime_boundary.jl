using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "PioneerSIMD owns runtime-only SIMD acceleration" begin
    simd_project = joinpath(REPO_ROOT, "packages", "PioneerSIMD", "Project.toml")
    simd_module = joinpath(REPO_ROOT, "packages", "PioneerSIMD", "src", "PioneerSIMD.jl")
    common_module = read(joinpath(REPO_ROOT, "packages", "PioneerCommon", "src", "PioneerCommon.jl"), String)
    predict_project = read(joinpath(REPO_ROOT, "packages", "PioneerPredict", "Project.toml"), String)
    predict_module = read(joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "PioneerPredict.jl"), String)
    search_project = read(joinpath(REPO_ROOT, "packages", "PioneerSearch", "Project.toml"), String)
    search_module = read(joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl"), String)
    search_entrypoint = read(joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "owned", "SearchDIA.jl"), String)

    @test isfile(simd_project)
    @test isfile(simd_module)
    @test occursin("include(joinpath(@__DIR__, \"..\", \"..\", \"..\", \"src\", \"shared\", \"simd_kernels.jl\"))", common_module)
    @test occursin("export load_pioneer_simd!", common_module)

    @test occursin("PioneerSIMD", predict_project)
    @test !occursin("LoopVectorization", predict_project)
    @test !occursin("LoopVectorization", predict_module)
    @test occursin("function __init__()", predict_module)
    @test occursin("load_pioneer_simd!()", predict_module)

    @test occursin("PioneerSIMD", search_project)
    @test !occursin("LoopVectorization", search_project)
    @test !occursin("LoopVectorization", search_module)
    @test occursin("function __init__()", search_module)
    @test occursin("load_pioneer_simd!()", search_module)
    @test !occursin("function asset_path(", search_entrypoint)

    @test !isfile(joinpath(REPO_ROOT, "apps", "PioneerPredictApp", "LocalPreferences.toml"))
    @test !isfile(joinpath(REPO_ROOT, "apps", "PioneerSearchApp", "LocalPreferences.toml"))
end
