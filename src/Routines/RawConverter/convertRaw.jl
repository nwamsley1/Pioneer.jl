# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Orbitrap-family `.raw` acquisition files -> Pioneer `.arrow` runs, read by the vendored SpectrumDecoder module.

const CONVERT_RAW_APP_NAME = "convertRaw"

"""
    convertRaw(path; output_dir = "", skip_existing = false, threads_per_file = 0) -> Vector{String}

Convert `.raw` acquisition files to Pioneer `.arrow` runs for SearchDIA. `path` is one `.raw` file or a folder
containing them (not searched recursively). Each file becomes `<output_dir>/<name>.arrow`; `output_dir` defaults to
`arrow_out` inside the folder, or beside the file when `path` is a file. Files are converted one at a time, each by
`threads_per_file` tasks (0: `min(4, Threads.nthreads())`), and written to a `.partial` file that is renamed when
complete. With `skip_existing`, files whose `.arrow` already exists are not converted.

A file that cannot be read is reported and the remaining files are still converted; an `ErrorException` listing the
failed files is thrown at the end. Returns the `.arrow` paths written (or skipped as existing).
"""
function convertRaw(path::AbstractString; output_dir::AbstractString = "", skip_existing::Bool = false,
                    threads_per_file::Integer = 0)
    src = rstrip(expanduser(String(path)), ['/', '\\'])
    is_raw(p) = isfile(p) && endswith(lowercase(p), ".raw")
    files = if is_raw(src)
        [src]
    elseif isdir(src)
        sort!(filter(is_raw, readdir(src; join = true)))
    else
        throw(ArgumentError("$path is not a .raw file or a folder containing them"))
    end
    isempty(files) && throw(ArgumentError("No .raw files in $path"))
    out = isempty(output_dir) ? joinpath(is_raw(src) ? dirname(src) : src, "arrow_out") : expanduser(String(output_dir))
    mkpath(out)
    threads_per_file < 0 && throw(ArgumentError("threads_per_file must be at least 0"))
    nt = threads_per_file == 0 ? SpectrumDecoder.default_threads() : min(Int(threads_per_file), Threads.nthreads())

    println("$(CONVERT_RAW_APP_NAME) $(get_pioneer_version())")
    println("Input : $(length(files)) .raw file$(length(files) == 1 ? "" : "s") from $src")
    println("Output: $out  (threads per file: $nt)")
    written = String[]
    failed = String[]
    for (i, raw) in enumerate(files)
        name = replace(basename(raw), r"\.raw$"i => "")
        dest = joinpath(out, name * ".arrow")
        tag = "[$i/$(length(files))] $name"
        if skip_existing && isfile(dest)
            println("$tag: $dest exists, skipped")
            push!(written, dest)
            continue
        end
        partial = dest * ".partial"
        try
            t = @elapsed begin
                f = SpectrumDecoder.open_acquisition(raw)
                SpectrumDecoder.write_pioneer_arrow(f, partial; threads = nt)
                mv(partial, dest; force = true)
            end
            n = length(SpectrumDecoder.scans(f))
            println("$tag: $n scans in $(round(t; digits = 1)) s")
            push!(written, dest)
        catch e
            rm(partial; force = true)
            println("$tag: FAILED: $(first(split(sprint(showerror, e), '\n')))")
            push!(failed, raw)
        end
    end
    isempty(failed) || error("$(length(failed)) of $(length(files)) files could not be converted:\n  " *
                             join(failed, "\n  "))
    return written
end

function show_convert_raw_help(io::IO = stdout)
    println(io, "$(CONVERT_RAW_APP_NAME) $(get_pioneer_version())")
    println(io)
    println(io, "Usage: $(CONVERT_RAW_APP_NAME) RAW_PATH [options]")
    println(io)
    println(io, "Arguments:")
    println(io, "  RAW_PATH                   .raw file, or a folder containing .raw files")
    println(io)
    println(io, "Options:")
    println(io, "  -o, --output-dir <path>    Output directory for .arrow runs (default: <input_dir>/arrow_out)")
    println(io, "      --skip-existing        Skip files whose .arrow output already exists")
    println(io, "  -t, --threads-per-file <n> Tasks used to convert each file (default: min(4, Julia threads))")
    println(io, "      --version              Show version information")
    println(io, "  -h, --help                 Show help information")
end

# Entry point for PackageCompiler
function main_convertRaw(argv = ARGS)::Cint
    args = String[a for a in argv]
    try
        if isempty(args) || any(in(("-h", "--help")), args)
            show_convert_raw_help()
            return isempty(args) ? 1 : 0
        end
        if "--version" in args
            println("$(CONVERT_RAW_APP_NAME) $(get_pioneer_version())")
            return 0
        end
        path = ""; output_dir = ""; skip_existing = false; threads_per_file = 0
        i = 1
        while i <= length(args)
            a = args[i]
            if a in ("-o", "--output-dir", "-t", "--threads-per-file")
                i += 1
                i > length(args) && throw(ArgumentError("Missing value for $a"))
                if a in ("-o", "--output-dir")
                    output_dir = args[i]
                else
                    n = tryparse(Int, args[i])
                    (n === nothing || n < 1) && throw(ArgumentError("$a needs a positive integer, got $(args[i])"))
                    threads_per_file = n
                end
            elseif a == "--skip-existing"
                skip_existing = true
            elseif startswith(a, "-")
                throw(ArgumentError("Unknown option: $a"))
            elseif isempty(path)
                path = a
            else
                throw(ArgumentError("Unexpected argument: $a"))
            end
            i += 1
        end
        convertRaw(path; output_dir = output_dir, skip_existing = skip_existing, threads_per_file = threads_per_file)
    catch e
        println(sprint(showerror, e))
        e isa ArgumentError && show_convert_raw_help()
        return 1
    end
    return 0
end
