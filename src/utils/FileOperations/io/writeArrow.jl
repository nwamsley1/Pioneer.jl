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

# Thread-safe lock for garbage collection calls to prevent race conditions
# when multiple threads call GC.gc() concurrently on Windows
const GC_LOCK = ReentrantLock()

# Audit harness: when PIONEER_AUDIT_WRITES is set to a writable path,
# every writeArrow call appends a TSV row recording {path, caller, nrows,
# ncols, in-mem bytes, semicolon-joined column names}. Lets us answer
# (1) is this write necessary? (2) are all columns necessary or expired?
const _AUDIT_WRITE_LOCK = ReentrantLock()

"""
    OUTPUT_DICT_ENCODED_COLUMNS

String columns in the final long-format outputs that are worth Arrow dictionary encoding. All are
massively repeated: on a 6-file KEAP1 run `file_name` stores 6 distinct strings across 729,008 rows
(21.6 MB for six values), `species` 3, `structural_mods` 1,367, `inferred_protein_group` 9,319,
`accession_numbers` 11,392, `sequence` 107,461 -- 46.3 MB of string payload in total.

Encoding is applied only at the FINAL output write, not to intermediate files: a `DictEncoded` column
reads back with a different Julia type, and confining it to the last write means no downstream code sees
a type change. Verified on the 141-column precursors_long: 366.3 -> 313.9 MB (-14.3%), every column
round-trips identically, and `DataFrame(Tables.columntable(...))` -- the access pattern qcPlots.jl uses
-- still yields `eltype == String`.
"""
const OUTPUT_DICT_ENCODED_COLUMNS = (:file_name, :species, :accession_numbers,
                                     :sequence, :structural_mods, :inferred_protein_group)

"""
    dict_encode_output_columns(tbl, encode = OUTPUT_DICT_ENCODED_COLUMNS)

Return `tbl` with the named string columns wrapped in `Arrow.DictEncode`. Returns `tbl` untouched when
none are present, so it is safe to apply to any output table. Column order and every other column are
preserved exactly.
"""
function dict_encode_output_columns(tbl, encode = OUTPUT_DICT_ENCODED_COLUMNS)
    all_names = Symbol.(Tables.columnnames(tbl))
    any(n -> n in encode, all_names) || return tbl
    cols = Any[]
    for n in all_names
        c = Tables.getcolumn(tbl, n)
        push!(cols, n in encode ? Arrow.DictEncode(c) : c)
    end
    return NamedTuple{Tuple(all_names)}(Tuple(cols))
end

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

#=
function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname()
        Arrow.write(tpath, df)
        if isfile(fpath)
            run(`cmd /c del /f "$fpath"`)
        end
        mv(tpath, fpath, force = true)
    else     #If Linux/MacOS easy
        Arrow.write(fpath, df)
    end
    return nothing
end
=#
function writeArrow(fpath::String, df::AbstractDataFrame)
    _audit_log_write(fpath, df)
    fpath = normpath(fpath)
    if Sys.iswindows()
        # Create a unique temporary file
        tpath = tempname() * ".arrow"
        # Write to the temporary file
        Arrow.write(tpath, df)
        # Try to delete the existing file with retries
        if isfile(fpath)
            try
                run(`cmd /c del /f /q "$fpath"`)#rm(fpath, force=true)
            catch e 
                #@user_info "Initial deletion failed for $fpath, attempting retries. \n"
                @debug_l1 "Windows specifiec deletion failed for $fpath. \n"
                # If that also fails, rename the old file instead of deleting
                backup_path = fpath * ".backup_" * string(time_ns())
                try
                    mv(fpath, backup_path, force=true)
                catch
                    @debug_l1 "Renaming original failed for $fpath, attempting GC. \n"
                    try
                        lock(GC_LOCK) do
                            GC.gc()
                        end
                        # Try to delete using Julia's rm with force flag
                        rm(fpath, force=true)
                        #break
                    catch
                        error("Unable to remove or rename existing file: $fpath")                            
                    end
                end
                #max_retries = 5
                #=
                for i in 1:max_retries
                    try
                        # Force garbage collection to release any file handles
                        # Use lock to prevent concurrent GC calls from multiple threads
                        lock(GC_LOCK) do
                            GC.gc()
                        end

                        # Try to delete using Julia's rm with force flag
                        rm(fpath, force=true)
                        break
                    catch e
                        if i == max_retries
                            # If all retries failed, try Windows-specific deletion
                            try
                                run(`cmd /c del /f /q "$fpath"`)
                            catch
                                # If that also fails, rename the old file instead of deleting
                                backup_path = fpath * ".backup_" * string(time_ns())
                                try
                                    mv(fpath, backup_path, force=true)
                                catch
                                    error("Unable to remove or rename existing file: $fpath")
                                end
                            end
                        else
                            # Wait a bit before retrying
                            sleep(0.1 * i)
                        end
                    end
                    @debug_l1 "Retry $i to delete $fpath failed"
                end
                =#
            end
        end
        
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