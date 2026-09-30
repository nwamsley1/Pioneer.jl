# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# SCIEX `.wiff` + `.wiff.scan` SWATH runs -> `.scxs` runs (SciexWiff, src/vendor/SciexWiff), which SearchDIA reads directly.


const CONVERT_SCIEX_APP_NAME = "convertSciex"

"""
    convertSciex(path; output_dir = "", zt_scan = nothing) -> Vector{String}

Convert SCIEX `.wiff` + `.wiff.scan` runs to `.scxs` runs for SearchDIA. `path` is one `.wiff` file or a folder
containing them; each `X.wiff` needs its `X.wiff.scan` (or `X.scan`) beside it. Each run becomes
`<output_dir>/<name>.scxs`; `output_dir` defaults to `scxs_out` next to the input. `.wiff2`-only runs are not
supported (the `.wiff2` format is encrypted); convert those with msConvert and `convertMzML`. Uses SciexWiff's
default centroiding. Returns the `.scxs` paths.

`zt_scan`: whether the batch is ZT Scan DIA. The `.wiff` does not record the scan mode, so it is always asked:
`nothing` prompts on an interactive terminal (default No) and is an error otherwise. ZT runs are written as
`<name>.zt.scxs` with `acquisition_type = zt_scan_dia` in their metadata.
"""
function convertSciex(path::AbstractString; output_dir::AbstractString = "",
                      zt_scan::Union{Nothing, Bool} = nothing)
    src = rstrip(expanduser(String(path)), ['/', '\\'])
    is_wiff(p) = isfile(p) && endswith(lowercase(p), ".wiff")
    runs = if is_wiff(src)
        [src]
    elseif isdir(src)
        sort!(filter(is_wiff, readdir(src; join = true)))
    else
        throw(ArgumentError("$path is not a SCIEX .wiff file or a folder containing them"))
    end
    isempty(runs) && throw(ArgumentError("No SCIEX .wiff files in $path"))
    zt = zt_scan === nothing ? ask_zt_scan() : zt_scan
    out = isempty(output_dir) ? joinpath(is_wiff(src) ? dirname(src) : src, "scxs_out") : expanduser(String(output_dir))
    mkpath(out)

    println("$(CONVERT_SCIEX_APP_NAME) $(get_pioneer_version())")
    println("Input : $(length(runs)) .wiff file$(length(runs) == 1 ? "" : "s") from $src")
    println("Output: $out  (threads: $(Threads.nthreads()))")
    println("Mode  : ", zt ? "ZT Scan DIA (runs written as <name>.zt.scxs)" : "SWATH / stepped DIA")
    params = SciexWiff.ConvertParams(format = :scxs, zt_scan = zt)
    written = String[]
    for (i, w) in enumerate(runs)
        name = splitext(basename(w))[1]
        println("\n[$i/$(length(runs))] $name")
        t = @elapsed paths = SciexWiff.convert_run(w, out; params = params, name = name)
        push!(written, paths.scxs)
        println("[$i/$(length(runs))] $name done in $(round(t; digits = 1)) s")
    end
    return written
end

"""
    ask_zt_scan(; input = stdin, output = stdout, interactive = input isa Base.TTY) -> Bool

Ask whether the batch being converted is ZT Scan DIA (default No). The `.wiff` does not record the scan mode,
so SCIEX conversion always asks; without an interactive terminal pass `--zt` / `--no-zt` (CLI) or `zt_scan`.
"""
function ask_zt_scan(; input::IO = stdin, output::IO = stdout, interactive::Bool = input isa Base.TTY)
    interactive || throw(ArgumentError(
        "Is this batch ZT Scan DIA data? No terminal to ask: pass --zt or --no-zt " *
        "(or zt_scan = true / false to convertSciex)"))
    while true
        print(output, "Is this batch ZT Scan DIA data? [y/N] "); flush(output)
        a = lowercase(strip(readline(input)))
        a in ("", "n", "no") && return false
        a in ("y", "yes") && return true
        println(output, "Please answer y or n.")
    end
end

function show_convert_sciex_help(io::IO = stdout)
    println(io, "$(CONVERT_SCIEX_APP_NAME) $(get_pioneer_version())")
    println(io)
    println(io, "Usage: $(CONVERT_SCIEX_APP_NAME) WIFF_PATH [options]")
    println(io)
    println(io, "Arguments:")
    println(io, "  WIFF_PATH                  SCIEX .wiff file (with its .wiff.scan), or a folder containing them")
    println(io)
    println(io, "Options:")
    println(io, "  -o, --output-dir <path>   Output directory for .scxs runs (default: <input_dir>/scxs_out)")
    println(io, "      --zt                  The runs are ZT Scan DIA (written as <name>.zt.scxs)")
    println(io, "      --no-zt               The runs are not ZT Scan DIA (SWATH / stepped DIA)")
    println(io, "                            Without either, you are asked (default No); required when not on a terminal")
    println(io, "      --version             Show version information")
    println(io, "  -h, --help                Show help information")
end

# Entry point for PackageCompiler
function main_convertSciex(argv = ARGS)::Cint
    args = String[a for a in argv]
    try
        if isempty(args) || any(in(("-h", "--help")), args)
            show_convert_sciex_help()
            return isempty(args) ? 1 : 0
        end
        if "--version" in args
            println("$(CONVERT_SCIEX_APP_NAME) $(get_pioneer_version())")
            return 0
        end
        path = ""; output_dir = ""; zt_scan = nothing
        i = 1
        while i <= length(args)
            a = args[i]
            if a == "-o" || a == "--output-dir"
                i += 1
                i > length(args) && throw(ArgumentError("Missing value for $a"))
                output_dir = args[i]
            elseif a == "--zt" || a == "--no-zt"
                zt_scan === nothing || zt_scan == (a == "--zt") ||
                    throw(ArgumentError("--zt and --no-zt cannot both be given"))
                zt_scan = a == "--zt"
            elseif startswith(a, "-")
                throw(ArgumentError("Unknown option: $a"))
            elseif isempty(path)
                path = a
            else
                throw(ArgumentError("Unexpected argument: $a"))
            end
            i += 1
        end
        convertSciex(path; output_dir = output_dir, zt_scan = zt_scan)
    catch e
        println(sprint(showerror, e))
        e isa ArgumentError && show_convert_sciex_help()
        return 1
    end
    return 0
end
