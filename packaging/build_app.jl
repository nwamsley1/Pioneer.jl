# Reproducible PackageCompiler entrypoint shared by all native build workflows.
#
# Run from the repository root:
#   julia --threads auto --project=packaging/compiler packaging/build_app.jl OUTPUT_DIR
#
# PIONEER_APP_CPU_TARGET may override the reviewed architecture default for a
# diagnostic build. Release workflows should leave it unset.

using PackageCompiler

const PACKAGING_DIR = @__DIR__
const REPOSITORY_ROOT = normpath(joinpath(PACKAGING_DIR, ".."))
const APP_PROJECT = REPOSITORY_ROOT
const COMPILER_PROJECT = joinpath(PACKAGING_DIR, "compiler", "Project.toml")
const PRECOMPILE_EXECUTION_FILE = joinpath(REPOSITORY_ROOT, "src", "build", "snoop.jl")
const PRECOMPILE_STATEMENTS_FILE = joinpath(
    REPOSITORY_ROOT,
    "src",
    "build",
    "precompile_statements.jl",
)

include(joinpath(PACKAGING_DIR, "app_config.jl"))

function app_cpu_target()
    override = strip(get(ENV, "PIONEER_APP_CPU_TARGET", ""))
    isempty(override) || return override

    haskey(PIONEER_APP_CPU_TARGETS, Sys.ARCH) || error(
        "no reviewed Pioneer application CPU target for $(Sys.ARCH); " *
        "set PIONEER_APP_CPU_TARGET explicitly for a diagnostic build",
    )
    return PIONEER_APP_CPU_TARGETS[Sys.ARCH]
end

function check_build_environment()
    active_project = Base.active_project()
    isnothing(active_project) && error(
        "no active compiler project; run with --project=packaging/compiler",
    )
    normpath(abspath(active_project)) == normpath(abspath(COMPILER_PROJECT)) || error(
        "wrong compiler project $(active_project); expected $(COMPILER_PROJECT)",
    )
    pkgversion(PackageCompiler) == v"2.4.0" || error(
        "PackageCompiler $(pkgversion(PackageCompiler)) is active; expected locked version 2.4.0",
    )

    for required_file in (
        joinpath(REPOSITORY_ROOT, "Manifest.toml"),
        joinpath(dirname(COMPILER_PROJECT), "Manifest.toml"),
        PRECOMPILE_EXECUTION_FILE,
        PRECOMPILE_STATEMENTS_FILE,
    )
        isfile(required_file) || error("missing application build input: $required_file")
    end
    return nothing
end

function build_app(output_dir::AbstractString)
    check_build_environment()
    output_path = abspath(output_dir)
    cpu_target = app_cpu_target()

    println("Pioneer application project: ", APP_PROJECT)
    println("Pioneer application output:  ", output_path)
    println("Pioneer application target:  ", cpu_target)

    cd(REPOSITORY_ROOT) do
        PackageCompiler.create_app(
            APP_PROJECT,
            output_path;
            incremental=false,
            force=true,
            executables=PIONEER_APP_EXECUTABLES,
            precompile_execution_file=PRECOMPILE_EXECUTION_FILE,
            precompile_statements_file=PRECOMPILE_STATEMENTS_FILE,
            cpu_target=cpu_target,
        )
    end
    return nothing
end

function main(args=ARGS)
    if args == ["--print-cpu-target"]
        println(app_cpu_target())
        return 0
    end
    length(args) == 1 || error(
        "usage: julia --project=packaging/compiler packaging/build_app.jl OUTPUT_DIR",
    )
    build_app(only(args))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
