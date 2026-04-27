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

@testset "Runtime-heavy apps disable automatic package precompile" begin
    unix_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer"), String)
    windows_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat"), String)
    linux_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_linux.yml"), String)
    macos_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)
    windows_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_windows.yml"), String)

    @test occursin("export JULIA_PKG_PRECOMPILE_AUTO=0", unix_launcher)
    @test occursin("set JULIA_PKG_PRECOMPILE_AUTO=0", windows_launcher)
    @test occursin("JULIA_PKG_PRECOMPILE_AUTO=0", linux_workflow)
    @test occursin("JULIA_PKG_PRECOMPILE_AUTO=0", macos_workflow)
    @test occursin("JULIA_PKG_PRECOMPILE_AUTO=0", windows_workflow)
end

@testset "Runtime-heavy apps use a writable runtime depot ahead of bundled sources" begin
    unix_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer"), String)
    windows_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat"), String)
    linux_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_linux.yml"), String)
    macos_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)
    windows_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_windows.yml"), String)

    @test occursin("runtime_depot", unix_launcher)
    @test occursin("mkdir -p \"\$runtime_depot\"", unix_launcher)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_depot:\$bundle_share", unix_launcher)

    @test occursin("RUNTIME_DEPOT", windows_launcher)
    @test occursin("LOCALAPPDATA", windows_launcher)
    @test occursin("JULIA_DEPOT_PATH=%RUNTIME_DEPOT%;%BUNDLE_SHARE%", windows_launcher)

    @test occursin("runtime_cache", linux_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$search_share", linux_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$predict_share", linux_workflow)

    @test occursin("runtime_cache", macos_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$search_share", macos_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$predict_share", macos_workflow)

    @test occursin("runtime_cache", windows_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache;\$search_share", windows_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache;\$predict_share", windows_workflow)
end

@testset "Installers do not make bundled Julia depots globally writable" begin
    linux_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_linux.yml"), String)
    macos_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)

    @test !occursin("chmod -R a+rwx", linux_workflow)
    @test !occursin("chmod -R a+w", macos_workflow)
end

@testset "Packaging smoke tests the packaged launcher runtime" begin
    linux_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_linux.yml"), String)
    macos_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)
    windows_workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_windows.yml"), String)

    for workflow in (linux_workflow, macos_workflow, windows_workflow)
        @test occursin("Smoke test packaged wrapper runtime", workflow)
    end

    @test occursin("\"\$pioneer_root/pioneer\" search --help", linux_workflow)
    @test occursin("\"\$pioneer_root/pioneer\" predict --help", linux_workflow)
    @test occursin("\"\$pioneer_root/pioneer\" search --help", macos_workflow)
    @test occursin("\"\$pioneer_root/pioneer\" predict --help", macos_workflow)
    @test occursin("pioneer.bat\" search --help", windows_workflow)
    @test occursin("pioneer.bat\" predict --help", windows_workflow)
end

@testset "Search launcher precompiles plotting runtime portably" begin
    unix_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer"), String)
    windows_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat"), String)

    @test occursin("precompile_search_plot_runtime", unix_launcher)
    @test occursin("-C generic", unix_launcher)
    @test occursin("--pkgimages=no", unix_launcher)
    @test occursin("PioneerPlots", unix_launcher)
    @test occursin("save_multipage_pdf", unix_launcher)

    @test occursin("precompile_search_plot_runtime", windows_launcher)
    @test occursin("-C generic", windows_launcher)
    @test occursin("--pkgimages=no", windows_launcher)
    @test occursin("PioneerPlots", windows_launcher)
    @test occursin("save_multipage_pdf", windows_launcher)
end
