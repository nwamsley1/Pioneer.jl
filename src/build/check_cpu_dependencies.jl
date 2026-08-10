# Fail fast if a portable application build reintroduces a package that makes
# AOT code generation depend on features of the build host.
#
# Run from a freshly resolved environment:
#   julia --project=. src/build/check_cpu_dependencies.jl

using Pioneer

const FORBIDDEN_CPU_PACKAGES = Set([
    "HostCPUFeatures",
    "LinearSolve",
    "LoopVectorization",
    "RecursiveFactorization",
    "TriangularSolve",
    "VectorizationBase",
])

loaded_packages = Set(pkgid.name for pkgid in keys(Base.loaded_modules))
unexpected_packages = intersect(loaded_packages, FORBIDDEN_CPU_PACKAGES)
isempty(unexpected_packages) ||
    error("portable AOT boundary loaded forbidden packages: $(sort!(collect(unexpected_packages)))")

project_text = read(joinpath(@__DIR__, "..", "..", "Project.toml"), String)
for package in ("HostCPUFeatures", "LinearSolve", "LoopVectorization")
    isnothing(match(Regex("(?m)^" * package * raw"\s*="), project_text)) ||
        error("$package must not be a direct Pioneer dependency")
end

println("Pioneer CPU dependency boundary: PASS")
