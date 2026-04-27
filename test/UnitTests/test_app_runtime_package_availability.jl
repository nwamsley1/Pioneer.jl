using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

include(joinpath(REPO_ROOT, "src", "build", "runtime_bundle_utils.jl"))

function run_app_package_check(project_name::String, script::String)
    project_path = joinpath(REPO_ROOT, "apps", project_name)
    cmd = Cmd([JULIA_BIN, "--startup-file=no", "--project=$project_path", "-e", script])
    run(cmd)
    return nothing
end

@testset "App projects expose runtime-only packages for lazy loading" begin
    search_script = """
    using PioneerSearchApp
    Base.require(Base.PkgId(Base.UUID("9adec8df-eb60-4bac-9b0a-4137287d3bbd"), "PioneerPlots"))
    Base.require(Base.PkgId(Base.UUID("05f169b2-1f58-4f3e-a5d2-7bfdd2b4f8e3"), "PioneerSIMD"))
    @assert Base.identify_package("PioneerPlots") !== nothing
    @assert Base.identify_package("PioneerSIMD") !== nothing
    """
    predict_script = """
    using PioneerPredictApp
    Base.require(Base.PkgId(Base.UUID("05f169b2-1f58-4f3e-a5d2-7bfdd2b4f8e3"), "PioneerSIMD"))
    @assert Base.identify_package("PioneerSIMD") !== nothing
    """

    run_app_package_check("PioneerSearchApp", search_script)
    run_app_package_check("PioneerPredictApp", predict_script)
    @test true
end

@testset "Bundled app environments vendor runtime-only packages" begin
    mktempdir() do tmp_dir
        search_bundle_dir = joinpath(tmp_dir, "search")
        predict_bundle_dir = joinpath(tmp_dir, "predict")

        search_share_dir = bundle_runtime_project!(
            joinpath(REPO_ROOT, "apps", "PioneerSearchApp"),
            search_bundle_dir,
        )
        predict_share_dir = bundle_runtime_project!(
            joinpath(REPO_ROOT, "apps", "PioneerPredictApp"),
            predict_bundle_dir,
        )

        search_project = TOML.parsefile(joinpath(search_share_dir, "Project.toml"))
        predict_project = TOML.parsefile(joinpath(predict_share_dir, "Project.toml"))

        @test !haskey(search_project, "name")
        @test !haskey(search_project, "uuid")
        @test search_project["sources"]["PioneerPlots"]["path"] == joinpath("dev", "PioneerPlots")
        @test search_project["sources"]["PioneerSIMD"]["path"] == joinpath("dev", "PioneerSIMD")
        @test isfile(joinpath(search_share_dir, "src", "shared", "version_utils.jl"))
        @test isfile(joinpath(search_share_dir, "dev", "PioneerPlots", "Project.toml"))
        @test isfile(joinpath(search_share_dir, "dev", "PioneerSIMD", "Project.toml"))

        @test !haskey(predict_project, "name")
        @test !haskey(predict_project, "uuid")
        @test predict_project["sources"]["PioneerPredict"]["path"] == joinpath("dev", "PioneerPredict")
        @test predict_project["sources"]["PioneerSIMD"]["path"] == joinpath("dev", "PioneerSIMD")
        @test isfile(joinpath(predict_share_dir, "src", "shared", "version_utils.jl"))
        @test isfile(joinpath(predict_share_dir, "dev", "PioneerPredict", "Project.toml"))
        @test isfile(joinpath(predict_share_dir, "dev", "PioneerSIMD", "Project.toml"))
    end

    @test true
end
