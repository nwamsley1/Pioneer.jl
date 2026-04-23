using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "PioneerPredict bootstrap loads build defaults before validation" begin
    bootstrap = read(joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "bootstrap.jl"), String)
    build_defaults = read(joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib", "utils", "buildParamDefaults.jl"), String)
    check_params = read(joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib", "utils", "check_params.jl"), String)

    @test occursin("buildParamDefaults.jl", bootstrap)
    @test findfirst("buildParamDefaults.jl", bootstrap) < findfirst("check_params.jl", bootstrap)
    @test occursin("function get_build_default_parameters()", build_defaults)
    @test occursin("function merge_with_build_defaults(user_params::AbstractDict, defaults::AbstractDict)", build_defaults)
    @test occursin("JSON.parse(json_string, dicttype=Dict{String,Any})", check_params)
    @test occursin("JSON.parsefile(json_path, dicttype=Dict{String,Any})", check_params)
end
