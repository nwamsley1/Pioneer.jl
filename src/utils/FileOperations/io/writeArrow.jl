# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Audit harness: when PIONEER_AUDIT_WRITES is set to a writable path,
# every writeArrow call appends a TSV row recording {path, caller, nrows,
# ncols, in-mem bytes, semicolon-joined column names}. Lets us answer
# (1) is this write necessary? (2) are all columns necessary or expired?
const _AUDIT_WRITE_LOCK = ReentrantLock()

function _audit_log_write(fpath::AbstractString, df::AbstractDataFrame)
    out = get(ENV, "PIONEER_AUDIT_WRITES", "")
    isempty(out) && return
    n = nrow(df); m = ncol(df)
    cb = 0
    @inbounds for c in 1:m
        col = df[!, c]
        et = eltype(col)
        cb += isbitstype(et) ? n * sizeof(et) : Base.summarysize(col)
    end
    bt = stacktrace()
    caller = "?"
    for i in 2:min(length(bt), 8)
        s = string(bt[i].file)
        occursin("writeArrow.jl", s) && continue
        caller = "$(basename(s)):$(bt[i].line) ($(bt[i].func))"
        break
    end
    cols_str = join(string.(names(df)), ";")
    lock(_AUDIT_WRITE_LOCK) do
        open(out, "a") do io
            println(io, fpath, "\t", caller, "\t", n, "\t", m, "\t", cb, "\t", cols_str)
        end
    end
    return
end

function writeArrow(fpath::String, df::AbstractDataFrame; temp_dir::AbstractString=tempdir())
    _audit_log_write(fpath, df)
    fpath = normpath(fpath)
    if Sys.iswindows()
        # Create a unique temporary file
        tpath = tempname(temp_dir) * ".arrow"
        # Write to the temporary file
        Arrow.write(tpath, df)
        # Route replacement through the same normalized, retrying deletion
        # path as every other Arrow cleanup operation.
        safeRm(fpath; force=true)

        # Move the temporary file to the final location
        try
            mv(tpath, fpath, force=true)
        catch e
            # If move fails, try copy and delete
            try
                cp(tpath, fpath, force=true)
                rm(tpath, force=true)
            catch
                error("Unable to write to file: $fpath")
            end
        end
    else
        # For Linux/MacOS, use temp file approach for safety
        # This avoids Bus errors when writing to a file that may still be memory-mapped
        tpath = tempname(temp_dir) * ".arrow"
        Arrow.write(tpath, df)
        mv(tpath, fpath, force=true)
    end
    return nothing
end
