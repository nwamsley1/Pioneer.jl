using Test
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "Portable app environments pin HostCPUFeatures to generic features" begin
    for app_name in ("PioneerPredictApp", "PioneerSearchApp")
        prefs_path = joinpath(REPO_ROOT, "apps", app_name, "LocalPreferences.toml")
        project_path = joinpath(REPO_ROOT, "apps", app_name, "Project.toml")
        @test isfile(prefs_path)

        prefs = TOML.parsefile(prefs_path)
        project = TOML.parsefile(project_path)
        @test haskey(prefs, "HostCPUFeatures")
        host_cpu_prefs = prefs["HostCPUFeatures"]
        @test get(get(project, "extras", Dict{String,Any}()), "HostCPUFeatures", nothing) == "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"

        @test get(host_cpu_prefs, "cpu_target", nothing) == "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"
        @test get(host_cpu_prefs, "freeze_cpu_target", nothing) === true
    end
end
