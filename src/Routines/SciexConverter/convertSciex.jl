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

# SCIEX `.wiff` + `.wiff.scan` SWATH runs -> `.scxs` runs (SciexWiff.jl), which SearchDIA reads directly.

import SciexWiff

const CONVERT_SCIEX_APP_NAME = "convertSciex"

"""
    convertSciex(path; output_dir = "") -> Vector{String}

Convert SCIEX `.wiff` + `.wiff.scan` SWATH runs to `.scxs` runs for SearchDIA. `path` is one `.wiff` file or a
folder containing them; each `X.wiff` needs its `X.wiff.scan` (or `X.scan`) beside it. Each run becomes
`<output_dir>/<name>.scxs`; `output_dir` defaults to `scxs_out` next to the input. `.wiff2`-only runs are not
supported (the `.wiff2` format is encrypted); convert those with msConvert and `convertMzML`. Uses SciexWiff's
default centroiding. Returns the `.scxs` paths.
"""
function convertSciex(path::AbstractString; output_dir::AbstractString = "")
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
    out = isempty(output_dir) ? joinpath(is_wiff(src) ? dirname(src) : src, "scxs_out") : expanduser(String(output_dir))
    mkpath(out)

    println("$(CONVERT_SCIEX_APP_NAME) $(get_pioneer_version())")
    println("Input : $(length(runs)) .wiff file$(length(runs) == 1 ? "" : "s") from $src")
    println("Output: $out  (threads: $(Threads.nthreads()))")
    params = SciexWiff.ConvertParams(format = :scxs)
    written = String[]
    for (i, w) in enumerate(runs)
        name = splitext(basename(w))[1]
        println("\n[$i/$(length(runs))] $name")
        t = @elapsed paths = SciexWiff.convert(w, out; params = params, name = name)
        push!(written, paths.scxs)
        println("[$i/$(length(runs))] $name done in $(round(t; digits = 1)) s")
    end
    return written
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
        path = ""; output_dir = ""
        i = 1
        while i <= length(args)
            a = args[i]
            if a == "-o" || a == "--output-dir"
                i += 1
                i > length(args) && throw(ArgumentError("Missing value for $a"))
                output_dir = args[i]
            elseif startswith(a, "-")
                throw(ArgumentError("Unknown option: $a"))
            elseif isempty(path)
                path = a
            else
                throw(ArgumentError("Unexpected argument: $a"))
            end
            i += 1
        end
        convertSciex(path; output_dir = output_dir)
    catch e
        println(sprint(showerror, e))
        e isa ArgumentError && show_convert_sciex_help()
        return 1
    end
    return 0
end
