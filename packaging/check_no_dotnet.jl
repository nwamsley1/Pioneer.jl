# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Fails when the application directory ARGS[1] holds .NET runtime or Thermo Fisher library files.
#   julia packaging/check_no_dotnet.jl <app_dir>

const FORBIDDEN = [
    r"^thermofisher\."i,            # Thermo Fisher libraries
    r"^pioneerconverter"i,          # the former .NET converter
    r"\.deps\.json$"i,              # .NET dependency manifests
    r"\.runtimeconfig\.json$"i,     # .NET runtime configuration
    r"^(lib)?coreclr\."i,           # .NET runtime
    r"^(lib)?hostfxr\."i,
    r"^(lib)?hostpolicy\."i,
    r"^(lib)?clrjit\."i,
    r"^system\.private\.corelib\.dll$"i,
    r"^mscorlib\.dll$"i,
    r"^createdump(\.exe)?$"i,
]

function main(dir)
    isdir(dir) || error("not a directory: $dir")
    found = String[]
    for (root, _, files) in walkdir(dir), f in files
        any(p -> occursin(p, f), FORBIDDEN) && push!(found, joinpath(root, f))
    end
    if !isempty(found)
        println(stderr, "Found $(length(found)) .NET / Thermo Fisher file(s) in $dir:")
        foreach(f -> println(stderr, "  ", f), first(found, 50))
        exit(1)
    end
    println("No .NET or Thermo Fisher files in $dir")
end

main(ARGS[1])
