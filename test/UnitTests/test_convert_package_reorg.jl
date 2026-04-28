using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CONVERT_PACKAGE_DIR = joinpath(REPO_ROOT, "packages", "PioneerConvert")
const DEV_JULIA = joinpath(Sys.BINDIR, Base.julia_exename())
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

subprocess_output = run_julia_code(
    """
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.develop(path=raw\"$CONVERT_PACKAGE_DIR\")
    Pkg.instantiate()
    using PioneerConvert
    println("convert=" * PioneerConvert.CONVERT_MZML_APP_NAME)
    println("main=" * string(isdefined(PioneerConvert, :main_convertMzML)))
    """;
    project=nothing
)

@testset "Extracted package loadability" begin
    @test occursin("convert=convertMzML", subprocess_output)
    @test occursin("main=true", subprocess_output)
end
