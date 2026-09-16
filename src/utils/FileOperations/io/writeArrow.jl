# Copyright (C) 2024 Nathan Wamsley
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

function writeArrow(fpath::String, df::AbstractDataFrame)
    _audit_log_write(fpath, df)
    fpath = normpath(fpath)
    if Sys.iswindows()
        # Create a unique temporary file
        tpath = tempname() * ".arrow"
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
        tpath = tempname() * ".arrow"
        Arrow.write(tpath, df)
        mv(tpath, fpath, force=true)
    end
    return nothing
end
