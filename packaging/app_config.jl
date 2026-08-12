# Shared declarative configuration for Pioneer's application contract and
# PackageCompiler. Keep executable names here instead of duplicating them in
# each operating-system workflow.

const PIONEER_APP_EXECUTABLES = Pair{String,String}[
    "GetSearchParams" => "main_GetSearchParams",
    "GetBuildLibParams" => "main_GetBuildLibParams",
    "GetParseSpecLibParams" => "main_GetParseSpecLibParams",
    "BuildSpecLib" => "main_BuildSpecLib",
    "SearchDIA" => "main_SearchDIA",
    "convertMzML" => "main_convertMzML",
]

# HostCPUFeatures-free code can safely participate in Julia's target cloning.
# Avoid PackageCompiler's Sandy Bridge clone_all tier for now: it substantially
# increases codegen work and exercised the LLVM legalization failure on Windows.
const PIONEER_APP_CPU_TARGETS = Dict{Symbol,String}(
    :x86_64 => "generic;haswell,-rdrnd",
    :aarch64 => "generic",
)
