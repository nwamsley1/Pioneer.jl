using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

include(joinpath(REPO_ROOT, "src", "build", "runtime_bundle_utils.jl"))

function run_project_check(project_path::String, script::String)
    cmd = Cmd([JULIA_BIN, "--startup-file=no", "--project=$project_path", "-e", script])
    run(cmd)
    return nothing
end

@testset "App projects expose compiled and runtime package dependencies" begin
    search_script = """
    using PioneerSearch, PioneerPlots
    @assert isdefined(PioneerSearch, :PioneerPlots)
    @assert isdefined(PioneerSearch, :Plots)
    @assert PioneerSearch.qcPlots === PioneerPlots.qcPlots
    @assert PioneerSearch.save_multipage_pdf === PioneerPlots.save_multipage_pdf
    @assert PioneerSearch.Plots === PioneerPlots.Plots
    @assert Base.identify_package("PioneerPlots") !== nothing
    """
    predict_script = """
    using PioneerPredict
    @assert Base.identify_package("PioneerSIMD") !== nothing
    """

    run_project_check(joinpath(REPO_ROOT, "packages", "PioneerSearch"), search_script)
    run_project_check(joinpath(REPO_ROOT, "packages", "PioneerPredict"), predict_script)
    @test true
end

@testset "Bundled app environments vendor app package dependencies" begin
    mktempdir() do tmp_dir
        search_bundle_dir = joinpath(tmp_dir, "search")
        predict_bundle_dir = joinpath(tmp_dir, "predict")
        cli_bundle_dir = joinpath(tmp_dir, "runtime")

        search_share_dir = bundle_runtime_project!(
            joinpath(REPO_ROOT, "packages", "PioneerSearch"),
            search_bundle_dir,
        )
        predict_share_dir = bundle_runtime_project!(
            joinpath(REPO_ROOT, "packages", "PioneerPredict"),
            predict_bundle_dir,
        )
        cli_share_dir = bundle_runtime_project!(
            joinpath(REPO_ROOT, "packages", "PioneerCLI"),
            cli_bundle_dir,
        )

        search_project = TOML.parsefile(joinpath(search_share_dir, "Project.toml"))
        predict_project = TOML.parsefile(joinpath(predict_share_dir, "Project.toml"))
        cli_project = TOML.parsefile(joinpath(cli_share_dir, "Project.toml"))

        @test !haskey(search_project, "name")
        @test !haskey(search_project, "uuid")
        @test search_project["sources"]["PioneerSearch"]["path"] == joinpath("dev", "PioneerSearch")
        @test search_project["sources"]["PioneerPlots"]["path"] == joinpath("dev", "PioneerPlots")
        @test search_project["sources"]["PioneerSIMD"]["path"] == joinpath("dev", "PioneerSIMD")
        @test isfile(joinpath(search_share_dir, "src", "shared", "version_utils.jl"))
        @test isfile(joinpath(search_share_dir, "stdlib", "v$(VERSION.major).$(VERSION.minor)", "Pkg", "src", "Pkg.jl"))
        @test isfile(joinpath(search_share_dir, "test", "testhelpers", "FakePTYs.jl"))
        @test isfile(joinpath(search_share_dir, "dev", "PioneerPlots", "Project.toml"))
        @test isfile(joinpath(search_share_dir, "dev", "PioneerSIMD", "Project.toml"))

        @test !haskey(predict_project, "name")
        @test !haskey(predict_project, "uuid")
        @test predict_project["sources"]["PioneerPredict"]["path"] == joinpath("dev", "PioneerPredict")
        @test predict_project["sources"]["PioneerSIMD"]["path"] == joinpath("dev", "PioneerSIMD")
        @test predict_project["sources"]["PioneerCommon"]["path"] == joinpath("dev", "PioneerCommon")
        @test predict_project["sources"]["PioneerParams"]["path"] == joinpath("dev", "PioneerParams")
        @test isfile(joinpath(predict_share_dir, "src", "shared", "version_utils.jl"))
        @test isfile(joinpath(predict_share_dir, "stdlib", "v$(VERSION.major).$(VERSION.minor)", "Pkg", "src", "Pkg.jl"))
        @test isfile(joinpath(predict_share_dir, "test", "testhelpers", "FakePTYs.jl"))
        @test isfile(joinpath(predict_share_dir, "dev", "PioneerPredict", "Project.toml"))
        @test isfile(joinpath(predict_share_dir, "dev", "PioneerSIMD", "Project.toml"))

        @test !haskey(cli_project, "name")
        @test !haskey(cli_project, "uuid")
        @test cli_project["sources"]["PioneerCLI"]["path"] == joinpath("dev", "PioneerCLI")
        @test cli_project["sources"]["PioneerSearch"]["path"] == joinpath("dev", "PioneerSearch")
        @test cli_project["sources"]["PioneerPredict"]["path"] == joinpath("dev", "PioneerPredict")
        @test cli_project["sources"]["PioneerSIMD"]["path"] == joinpath("dev", "PioneerSIMD")
        @test isfile(joinpath(cli_share_dir, "stdlib", "v$(VERSION.major).$(VERSION.minor)", "Pkg", "src", "Pkg.jl"))
        @test isfile(joinpath(cli_share_dir, "test", "testhelpers", "FakePTYs.jl"))
        @test isfile(joinpath(cli_share_dir, "dev", "PioneerCLI", "Project.toml"))
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
    @test occursin("runtime_share=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/runtime/share/julia\"", macos_workflow)
    @test occursin("search_share=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/search/share/julia\"", macos_workflow)
    @test occursin("predict_share=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/predict/share/julia\"", macos_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$runtime_share", macos_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$search_share", macos_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache:\$predict_share", macos_workflow)

    @test occursin("runtime_cache", windows_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache;\$search_share", windows_workflow)
    @test occursin("JULIA_DEPOT_PATH=\"\$runtime_cache;\$predict_share", windows_workflow)
end

@testset "Bundled asset lookup is safe in REPL-style sessions" begin
    script = """
    using Pioneer
    @assert isempty(PROGRAM_FILE)
    @assert isdefined(Pioneer, :Runtime)
    @assert Pioneer.Runtime === Pioneer
    missing_asset = Pioneer.asset_path("__definitely_missing_pioneer_asset__")
    @assert endswith(missing_asset, joinpath("assets", "__definitely_missing_pioneer_asset__"))
    println(missing_asset)
    """
    cmd = Cmd([JULIA_BIN, "--startup-file=no", "--project=$REPO_ROOT", "-e", script])
    @test success(cmd)
end

@testset "Windows launcher preserves subcommand arguments with spaces" begin
    windows_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat"), String)

    @test occursin(":collect_subcommand_args", windows_launcher)
    @test occursin("set SUBCOMMAND_ARGS=\"%~1\"", windows_launcher)
    @test occursin("set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% \"%~1\"", windows_launcher)
    @test !occursin("set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% %~1", windows_launcher)
    @test !occursin("set SUBCOMMAND_ARGS=%~1", windows_launcher)
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

    @test occursin("search_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/search/bin/julia\"", linux_workflow)
    @test occursin("predict_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/predict/bin/julia\"", linux_workflow)
    @test occursin("\"\$search_julia\" --startup-file=no --pkgimages=no --project=\"\$search_share\"", linux_workflow)
    @test occursin("\"\$predict_julia\" --startup-file=no --pkgimages=no --project=\"\$predict_share\"", linux_workflow)

    @test occursin("runtime_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/runtime/bin/julia\"", macos_workflow)
    @test occursin("search_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/search/bin/julia\"", macos_workflow)
    @test occursin("predict_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/predict/bin/julia\"", macos_workflow)
    @test occursin("\"\$runtime_julia\" --startup-file=no --pkgimages=no --project=\"\$runtime_share\"", macos_workflow)
    @test occursin("\"\$search_julia\" --startup-file=no --pkgimages=no --project=\"\$search_share\"", macos_workflow)
    @test occursin("\"\$predict_julia\" --startup-file=no --pkgimages=no --project=\"\$predict_share\"", macos_workflow)

    @test occursin("search_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/search/bin/julia.exe\"", windows_workflow)
    @test occursin("predict_julia=\"build/Pioneer_\${{ matrix.identifier }}/Applications/Pioneer/apps/predict/bin/julia.exe\"", windows_workflow)
    @test occursin("\"\$search_julia\" --startup-file=no --pkgimages=no --project=\"\$search_share\"", windows_workflow)
    @test occursin("\"\$predict_julia\" --startup-file=no --pkgimages=no --project=\"\$predict_share\"", windows_workflow)
end

@testset "Search package compiles mandatory plotting into the package boundary" begin
    search_module = read(joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "PioneerSearch.jl"), String)
    unix_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer"), String)
    windows_launcher = read(joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat"), String)

    @test occursin("\nusing PioneerPlots\n", search_module)
    @test occursin("const Plots = PioneerPlots.Plots", search_module)
    @test !occursin("qcPlots(args...; kwargs...)", search_module)
    @test !isfile(joinpath(REPO_ROOT, "apps", "PioneerSearchApp", "Project.toml"))
    @test !isfile(joinpath(REPO_ROOT, "apps", "PioneerSearchApp", "src", "PioneerSearchApp.jl"))

    @test !occursin("precompile_search_plot_runtime", unix_launcher)
    @test !occursin("pioneer-plots-runtime-smoke", unix_launcher)
    @test !occursin("precompile_search_plot_runtime", windows_launcher)
    @test !occursin("pioneer-plots-runtime-smoke", windows_launcher)
end
