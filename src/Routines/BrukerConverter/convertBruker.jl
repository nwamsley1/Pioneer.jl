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

# Bruker timsTOF diaPASEF `.d` bundles -> `.tdfs` runs (TimsSlices.jl), the format SearchDIA reads for timsTOF data.

const CONVERT_BRUKER_APP_NAME = "convertBruker"

"""
    convertBruker(path; output_dir = "") -> Vector{String}

Convert Bruker timsTOF diaPASEF `.d` bundles to `.tdfs` runs for SearchDIA. `path` is one `.d` folder or a folder
containing them. Each bundle becomes `<output_dir>/<name>.tdfs`; `output_dir` defaults to `tdfs_out` next to the
bundles. Uses TimsSlices' default slicing (8 IM scans per slice, IM sigma 5 scans, m/z sigma 3 bins, summed
intensities, no cull), the settings Pioneer's timsTOF search was validated with. Returns the `.tdfs` paths.
"""
function convertBruker(path::AbstractString; output_dir::AbstractString = "")
    src = rstrip(expanduser(String(path)), ['/', '\\'])
    is_d(p) = isdir(p) && endswith(lowercase(p), ".d")
    bundles = if is_d(src)
        [src]
    elseif isdir(src)
        sort!(filter(is_d, readdir(src; join = true)))
    else
        throw(ArgumentError("$path is not a Bruker .d folder or a folder containing them"))
    end
    isempty(bundles) && throw(ArgumentError("No Bruker .d folders in $path"))
    out = isempty(output_dir) ? joinpath(is_d(src) ? dirname(src) : src, "tdfs_out") : expanduser(String(output_dir))
    mkpath(out)

    println("$(CONVERT_BRUKER_APP_NAME) $(get_pioneer_version())")
    println("Input : $(length(bundles)) .d bundle$(length(bundles) == 1 ? "" : "s") from $src")
    println("Output: $out  (threads: $(Threads.nthreads()))")
    written = String[]
    for (i, d) in enumerate(bundles)
        name = replace(basename(d), r"\.d$"i => "")
        println("\n[$i/$(length(bundles))] $name")
        t = @elapsed paths = TimsSlices.convert(d, out; name = name)
        push!(written, paths.tdfs)
        println("[$i/$(length(bundles))] $name done in $(round(t; digits = 1)) s")
    end
    return written
end

function show_convert_bruker_help(io::IO = stdout)
    println(io, "$(CONVERT_BRUKER_APP_NAME) $(get_pioneer_version())")
    println(io)
    println(io, "Usage: $(CONVERT_BRUKER_APP_NAME) D_PATH [options]")
    println(io)
    println(io, "Arguments:")
    println(io, "  D_PATH                     Bruker timsTOF .d folder, or a folder containing .d folders")
    println(io)
    println(io, "Options:")
    println(io, "  -o, --output-dir <path>   Output directory for .tdfs runs (default: <input_dir>/tdfs_out)")
    println(io, "      --version             Show version information")
    println(io, "  -h, --help                Show help information")
end

# Entry point for PackageCompiler
function main_convertBruker(argv = ARGS)::Cint
    args = String[a for a in argv]
    try
        if isempty(args) || any(in(("-h", "--help")), args)
            show_convert_bruker_help()
            return isempty(args) ? 1 : 0
        end
        if "--version" in args
            println("$(CONVERT_BRUKER_APP_NAME) $(get_pioneer_version())")
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
        convertBruker(path; output_dir = output_dir)
    catch e
        println(sprint(showerror, e))
        e isa ArgumentError && show_convert_bruker_help()
        return 1
    end
    return 0
end
