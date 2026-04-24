using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

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
