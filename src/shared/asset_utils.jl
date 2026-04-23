"""
    asset_path(parts...)

Return the path to a bundled asset. The function first checks the compile-time
repository `assets/` directory and falls back to a path relative to the
installed executable.
"""
function asset_path(parts...)
    compile_dir = joinpath(@__DIR__, "..", "..", "assets", parts...)
    if ispath(compile_dir)
        return normpath(compile_dir)
    end

    exe = PROGRAM_FILE
    if !isabspath(exe)
        exe_full = Sys.which(exe)
        exe = exe_full !== nothing ? exe_full : exe
    end

    exe_dir = abspath(dirname(realpath(exe)))
    return normpath(joinpath(exe_dir, "..", "data", parts...))
end
