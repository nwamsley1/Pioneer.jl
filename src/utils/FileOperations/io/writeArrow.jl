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

"""
    OUTPUT_DICT_ENCODED_COLUMNS

Low-cardinality string columns in the final long-format outputs that are safe and worthwhile to Arrow
dictionary-encode. Their distinct-value count is bounded by run structure regardless of library size:
`file_name` (one per MS file), `species` (1-3), `structural_mods` (~hundreds of mod combinations) --
all far below 32,767. On a 6-file KEAP1 run `file_name` stores 6 distinct strings across 729,008 rows
(21.6 MB for six values), `species` 3, `structural_mods` 1,367; these massively-repeated columns give
the bulk of the size win.

⚠️ High-cardinality identity columns (`sequence`, `accession_numbers`, `inferred_protein_group`) are
DELIBERATELY NOT encoded. The precursors_long output is written to a persistent `Arrow.Writer` one merge
chunk at a time, and each chunk reads back as a `ChainedVector` -> multiple record batches. Arrow.jl locks
a column's dictionary index integer type from the FIRST batch's cardinality
(`PooledArray(...; compress=true)`); if a later batch pushes the cumulative unique count past that type's
max it throws "fatal error writing arrow data" (e.g. Int16 -> >32,767). Because the chunks are
protein-group-sorted, the first batch sees only a slice of peptides, so `sequence` (100k-300k+ uniques)
locks `Int16` and overflows on the next batch. Keeping these three as plain strings is robust across yeast
and human/isoform libraries; they compress the worst under dict-encoding anyway, so little size is lost.

Encoding is applied only at the FINAL output write, not to intermediate files: a `DictEncoded` column
reads back with a different Julia type, and confining it to the last write means no downstream code sees
a type change. `DataFrame(Tables.columntable(...))` -- the access pattern qcPlots.jl uses -- still yields
`eltype == String` for both encoded and plain columns.
"""
const OUTPUT_DICT_ENCODED_COLUMNS = (:file_name, :species, :structural_mods)

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
