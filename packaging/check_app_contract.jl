# Fast check for drift between the canonical Pioneer module and the executable
# list consumed by PackageCompiler.
#
# Run:
#   julia --project=. packaging/check_app_contract.jl

using Pioneer

const PACKAGING_DIR = @__DIR__
const REPOSITORY_ROOT = normpath(joinpath(PACKAGING_DIR, ".."))

include(joinpath(PACKAGING_DIR, "app_config.jl"))

realpath(pkgdir(Pioneer)) == realpath(REPOSITORY_ROOT) || error(
    "Pioneer loaded from $(pkgdir(Pioneer)); expected $REPOSITORY_ROOT",
)

for (_, entrypoint_name) in PIONEER_APP_EXECUTABLES
    entrypoint = Symbol(entrypoint_name)
    isdefined(Pioneer, entrypoint) || error(
        "Pioneer application entrypoint is missing: $entrypoint_name",
    )
    applicable(getfield(Pioneer, entrypoint), String[]) || error(
        "Pioneer.$entrypoint_name must accept an argument vector",
    )
end

println(
    "Pioneer application contract: PASS (",
    length(PIONEER_APP_EXECUTABLES),
    " entrypoints)",
)
