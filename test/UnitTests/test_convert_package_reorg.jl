using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONVERT_PACKAGE_DIR = joinpath(REPO_ROOT, "packages", "PioneerConvert")
const CONVERT_APP_DIR = joinpath(REPO_ROOT, "apps", "PioneerConvertApp")
const DEV_JULIA = "/Users/dennisgoldfarb/.julia/juliaup/julia-1.12.6+0.x64.apple.darwin14/Julia-1.12.app/Contents/Resources/julia/bin/julia"
const TEST_DEPOT = mktempdir()

function run_julia_code(code::String; workdir::String=REPO_ROOT, project::Union{Nothing, String}=REPO_ROOT)
    args = String[DEV_JULIA, "--startup-file=no"]
    if project !== nothing
        push!(args, "--project=$(project)")
    end
    append!(args, ["-e", code])

    cmd = Cmd(Cmd(args); dir=workdir)
    env = copy(ENV)
    env["JULIA_DEPOT_PATH"] = "$(TEST_DEPOT):$(joinpath(homedir(), ".julia"))"
    return read(setenv(cmd, env), String)
end

@testset "PioneerConvert package extraction" begin
    @test isfile(joinpath(CONVERT_PACKAGE_DIR, "Project.toml"))
    @test isfile(joinpath(CONVERT_PACKAGE_DIR, "src", "PioneerConvert.jl"))
end

@testset "PioneerConvertApp package extraction" begin
    @test isfile(joinpath(CONVERT_APP_DIR, "Project.toml"))
    @test isfile(joinpath(CONVERT_APP_DIR, "src", "PioneerConvertApp.jl"))
end

subprocess_output = run_julia_code(
    """
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.develop(path=raw\"$CONVERT_PACKAGE_DIR\")
    Pkg.develop(path=raw\"$CONVERT_APP_DIR\")
    Pkg.instantiate()
    using PioneerConvert
    using PioneerConvertApp
    println("convert=" * PioneerConvert.CONVERT_MZML_APP_NAME)
    println("app=" * string(isdefined(PioneerConvertApp, :main_convertMzML)))
    """;
    project=nothing
)

@testset "Extracted package loadability" begin
    @test occursin("convert=convertMzML", subprocess_output)
    @test occursin("app=true", subprocess_output)
end
