# Fast PackageCompiler guard for portable application builds.
#
# This incremental sysimage compiles a dependency-free kernel package that
# mirrors Pioneer's Base @simd reduction and fixed-width LLVM vector primitive,
# then launches the generated image and repeats the workload. The separate
# check_cpu_dependencies.jl guard proves Pioneer does not load the packages
# that previously leaked build-host features. The full non-incremental app
# build remains the release gate.
#
# Run:
#   julia --project=packaging/compiler src/build/aot_canary.jl
# Defaults to PackageCompiler's application target. Override when the full app
# build uses a different target:
#   PIONEER_AOT_CPU_TARGET='generic;haswell,-rdrnd' \
#       julia --project=packaging/compiler src/build/aot_canary.jl

using Libdl
using PackageCompiler

const CANARY_PROJECT = joinpath(@__DIR__, "aot_canary")
const WORKLOAD = joinpath(CANARY_PROJECT, "workload.jl")
const DEFAULT_CPU_TARGET = PackageCompiler.default_app_cpu_target()
const CPU_TARGET = get(ENV, "PIONEER_AOT_CPU_TARGET", DEFAULT_CPU_TARGET)

println("Pioneer AOT canary: building target ", repr(CPU_TARGET))
mktempdir(prefix="pioneer-aot-canary-") do output_dir
    sysimage_path = joinpath(output_dir, "PioneerAOTCanary.$(Libdl.dlext)")
    PackageCompiler.create_sysimage(
        ["PioneerAOTCanary"];
        project=CANARY_PROJECT,
        sysimage_path=sysimage_path,
        incremental=true,
        cpu_target=CPU_TARGET,
        precompile_execution_file=WORKLOAD,
    )

    command = `$(Base.julia_cmd()) --startup-file=no --project=$CANARY_PROJECT --sysimage=$sysimage_path $WORKLOAD`
    run(command)
end
println("Pioneer AOT canary: PASS")
