using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "Split package bootstraps fail fast on required source failures" begin
    bootstraps = [
        joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "bootstrap.jl"),
        joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "bootstrap.jl"),
    ]

    for bootstrap_path in bootstraps
        bootstrap = read(bootstrap_path, String)

        @test !occursin("@user_warn \"File not found", bootstrap)
        @test !occursin("@user_warn \"Failed to include", bootstrap)
        @test !occursin("catch e", bootstrap)
        @test occursin("throw(ArgumentError(\"Required source file not found:", bootstrap)
        @test occursin("include(file_path)", bootstrap)
        @test occursin("push!(files_loaded, file_path)", bootstrap)
    end
end
