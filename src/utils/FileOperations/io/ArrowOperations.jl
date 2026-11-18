# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
Arrow-specific file operations and utilities.

Provides specialized operations for Arrow files including
sorting, writing, and file management with Windows compatibility.
"""

using Arrow, DataFrames, Tables
using Mmap: Mmap

#==========================================================
Arrow File Operations
==========================================================#

"""
    with_arrow_table(f, path)

Memory-map `path`, call `f(Arrow.Table(...))`, and unmap the file before returning.

`Arrow.Table(path)` releases its mapping only when the garbage collector finalizes it,
and that can be much later than the table going out of scope: Arrow reads record batches
on spawned tasks, and a worker thread keeps the last task it ran (with the batch, and so
the mapped bytes) reachable until it runs another one. A dead temporary can also stay
rooted in the calling function's frame until that function returns. Until then the file
counts toward resident memory, keeps its disk space on Linux, and cannot be deleted or
replaced on a network share. `f` must copy what it needs: its result must not reference
the table's columns, which become invalid when the file is unmapped.
"""
function with_arrow_table(f, path::AbstractString)
    bytes = Mmap.mmap(path)
    try
        return f(Arrow.Table(bytes))
    finally
        # Mmap attaches its unmap finalizer to the array's backing Memory.
        Base.finalize(bytes.ref.mem)
    end
end

"""
    with_arrow_stream(f, path)

Like `with_arrow_table`, for record-batch iteration: call `f(Arrow.Stream(...))` over a
memory map of `path` and unmap the file before returning. Nothing that references a
batch's columns may be used afterwards.
"""
function with_arrow_stream(f, path::AbstractString)
    bytes = Mmap.mmap(path)
    try
        return f(Arrow.Stream(bytes))
    finally
        Base.finalize(bytes.ref.mem)
    end
end

"""
    with_arrow_tables(f)

Call `f(open_table)`, where `open_table(path)` memory-maps `path` and returns its
`Arrow.Table`, then unmap every file opened that way before returning. Use it when several
tables must stay mapped together (for example while they are concatenated and written);
as with `with_arrow_table`, nothing that references their columns may be used afterwards.
"""
function with_arrow_tables(f)
    mapped = Vector{Vector{UInt8}}()
    open_table = function (path::AbstractString)
        bytes = Mmap.mmap(path)
        push!(mapped, bytes)
        return Arrow.Table(bytes)
    end
    try
        return f(open_table)
    finally
        foreach(bytes -> Base.finalize(bytes.ref.mem), mapped)
    end
end

"""
    load_arrow_dataframe(path; cols = nothing) -> DataFrame

Read `path` into a DataFrame that owns all of its memory, then unmap the file. Column types are
those of `DataFrame(Tables.columntable(Arrow.Table(path)))`, or with `cols` those of
`DataFrame([c => collect(Tables.getcolumn(tbl, c)) for c in cols])` (columns missing from the
file are skipped), except that list columns, whose cells are views into the file, get their
cells copied. Every column is checked before the file is unmapped; a column type that might
still point into the file is an error rather than a read of released memory.
"""
function load_arrow_dataframe(path::AbstractString; cols = nothing)
    return with_arrow_table(path) do tbl
        df = cols === nothing ? DataFrame(Tables.columntable(tbl)) :
            DataFrame([c => collect(Tables.getcolumn(tbl, c)) for c in cols if hasproperty(tbl, c)];
                      copycols = false)
        for name in names(df)
            col = _owned_column(df[!, name])
            _owns_memory(col) || error("load_arrow_dataframe: column $name of $path has type " *
                                       "$(typeof(col)), which may still reference the mapped file")
            df[!, name] = col
        end
        df
    end
end

# List cells come back as views into the file; copy each into a Vector of the cell's element
# type. A nested list keeps views in its copied cells, so `_owns_memory` rejects it.
function _owned_column(col::AbstractVector)
    S = nonmissingtype(eltype(col))
    S <: AbstractVector || return col
    C = Vector{eltype(S)}
    out = Vector{S === eltype(col) ? C : Union{Missing, C}}(undef, length(col))
    @inbounds for i in eachindex(col, out)
        x = col[i]
        out[i] = x === missing ? missing : C(x)
    end
    return out
end

_owned_eltype(T::Type) =
    isbitstype(T) || T === String || T === Symbol ||
    (T isa Union && all(_owned_eltype, Base.uniontypes(T))) ||
    (T <: Vector && T !== Vector && _owned_eltype(eltype(T)))
_owns_memory(col::Vector) = _owned_eltype(eltype(col))
# PooledVector (dictionary-encoded columns): its own refs and pool.
_owns_memory(col::AbstractVector) =
    nameof(typeof(col)) === :PooledArray && col.refs isa Vector && _owns_memory(col.pool)

"""
    sort_file_by_keys!(ref::FileReference, sort_keys::Symbol...; 
                      reverse::Union{Bool, Vector{Bool}}=false)

Sort a file in-place by the specified keys.
Updates the reference's sorted_by metadata.

# Arguments
- `ref`: FileReference to sort
- `sort_keys`: Column names to sort by
- `reverse`: Either a single Bool (applies to all keys) or a Vector{Bool} (one per key)

# Returns
- The updated FileReference

# Examples
```julia
# Sort by single key descending
sort_file_by_keys!(ref, :score; reverse=true)

# Sort by multiple keys with different directions
sort_file_by_keys!(ref, :score, :target; reverse=[true, true])
```
"""
function sort_file_by_keys!(ref::FileReference, sort_keys::Symbol...;
                           reverse::Union{Bool, Vector{Bool}}=false,
                           show_progress::Bool=true)
    validate_exists(ref)

    # Validate that all sort keys exist in schema (or in any sidecar, for
    # PSMFileReference). This is the natural consolidation point: the sort
    # produces a flat file that includes previously-sidecarred columns.
    for key in sort_keys
        if ref isa PSMFileReference
            has_column_anywhere(ref, key) ||
                error("Sort key $key not found in file schema or any sidecar")
        else
            has_column(schema(ref), key) ||
                error("Sort key $key not found in file schema")
        end
    end

    # Load file (and any registered sidecars) into memory.
    df = load_with_sidecars(ref)

    # Build sort order vector
    rev_vec = if reverse isa Bool
        fill(reverse, length(sort_keys))
    else
        # Validate length matches
        if length(reverse) != length(sort_keys)
            error("Length of reverse vector ($(length(reverse))) must match number of sort keys ($(length(sort_keys)))")
        end
        reverse
    end

    # Sort dataframe
    fast_df_sort!(df, collect(sort_keys), rev=rev_vec)

    # Write back to same file using writeArrow for Windows compatibility.
    # This is now the consolidation point — main + sidecar columns are
    # merged into one file.
    writeArrow(file_path(ref), df)

    # Sidecars are now stale (their data is in main); clear them and
    # unlink the on-disk sidecar files.
    if ref isa PSMFileReference && !isempty(ref.sidecars)
        clear_sidecars!(ref; delete_files=true)
    end

    # Update reference metadata
    mark_sorted!(ref, sort_keys...)

    return ref
end

"""
    sort_file_by_keys!(refs::Vector{<:FileReference}, keys::Symbol...; 
                      reverse::Union{Bool, Vector{Bool}}=false, parallel::Bool=true)

Sort multiple files by the specified keys, optionally in parallel.

# Arguments
- `refs`: Vector of FileReferences to sort
- `keys`: Column names to sort by
- `reverse`: Either a single Bool (applies to all keys) or a Vector{Bool} (one per key)
- `parallel`: If true, sort files in parallel using threads

# Examples
```julia
# Sort all files by score descending
sort_file_by_keys!(refs, :score; reverse=true)

# Sort all files by score descending, then target descending
sort_file_by_keys!(refs, :score, :target; reverse=[true, true])
```
"""
function sort_file_by_keys!(refs::Vector{<:FileReference}, keys::Symbol...; 
                           reverse::Union{Bool, Vector{Bool}}=false, parallel::Bool=true,
                           show_progress::Bool=true)
    # Validate reverse vector length if it's a vector
    if reverse isa Vector{Bool} && length(reverse) != length(keys)
        error("Length of reverse vector ($(length(reverse))) must match number of sort keys ($(length(keys)))")
    end
    
    if parallel && length(refs) > 1
        if show_progress
            Threads.@threads for ref in refs
                if exists(ref)
                    sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                end
            end
        else
            Threads.@threads for ref in refs  # No ProgressBar
                if exists(ref)
                    sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                end
            end
        end
    else
        # Sequential processing
        if show_progress
            for ref in refs
                if exists(ref)
                    sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                end
            end
        else
            for ref in refs  # No ProgressBar
                if exists(ref)
                    sort_file_by_keys!(ref, keys...; reverse=reverse, show_progress=false)
                end
            end
        end
    end
    return refs
end


"""
    write_arrow_file(ref::FileReference, df::DataFrame) -> FileReference
    
Write a DataFrame to the file referenced by ref, updating all metadata.
Uses writeArrow from utils/writeArrow.jl to handle Windows file locking issues.
"""
function write_arrow_file(ref::FileReference, df::DataFrame)
    # Use writeArrow which handles Windows-specific issues
    writeArrow(file_path(ref), df)
    
    # Update reference metadata
    new_ref = create_reference(file_path(ref), typeof(ref))
    ref.schema = schema(new_ref)
    ref.row_count = row_count(new_ref)
    ref.sorted_by = ()  # Reset sort state as we don't know if df is sorted
    
    return ref
end

"""
    transform_and_write!(transform_fn::Function, ref::FileReference) -> FileReference
    
Load entire file, apply transformation, and write back.
For operations that need full dataset access (like sorting).
"""
function transform_and_write!(transform_fn::Function, ref::FileReference)
    validate_exists(ref)

    # Load main + any registered sidecars. The in-place rewrite is the
    # natural consolidation point for the sidecar architecture: after this
    # function returns, all sidecar columns have been baked into the main
    # file and the stale sidecar files on disk are deleted.
    had_sidecars = ref isa PSMFileReference && !isempty(ref.sidecars)
    df = load_with_sidecars(ref)

    # Apply transformation
    transformed_df = transform_fn(df)

    # Write back using writeArrow for Windows compatibility
    write_arrow_file(ref, transformed_df)

    # Sidecar data is now in main; remove the orphaned sidecar files.
    had_sidecars && clear_sidecars!(ref; delete_files=true)

    return ref
end

"""
    transform_and_write!(transform_fn::Function, ref::FileReference, output_path::String) -> FileReference
    
Load entire file, apply transformation, and write to a new location.
Does not modify the original file.
"""
function transform_and_write!(transform_fn::Function, ref::FileReference, output_path::String)
    validate_exists(ref)

    # Ensure output directory exists
    output_dir = dirname(output_path)
    !isdir(output_dir) && mkpath(output_dir)

    # Load main file + any registered sidecars (PSMFileReference specialization)
    df = load_with_sidecars(ref)
    
    # Apply transformation
    transformed_df = transform_fn(df)
    
    # Write to output path using writeArrow for Windows compatibility
    writeArrow(output_path, transformed_df)
    
    # Create reference for output of same type as input
    return create_reference(output_path, typeof(ref))
end

#==========================================================
Helper Functions for Arrow Files
==========================================================#

"""
    load_dataframe(ref::FileReference) -> DataFrame

Load entire Arrow file as DataFrame. Hides direct Arrow.Table access.
"""
function load_dataframe(ref::FileReference)
    validate_exists(ref)
    # Sidecar-aware read: for PSMFileReference with registered sidecars, this
    # returns a DataFrame containing main + all sidecar columns. For plain
    # FileReferences, identical to the original DataFrame(Tables.columntable(...))
    # behavior.
    return load_with_sidecars(ref)
end

"""
    column_names(ref::FileReference) -> Vector{Symbol}

Get column names without loading the full dataset.
"""
function column_names(ref::FileReference)
    validate_exists(ref)
    return with_arrow_table(table -> Symbol.(Tables.columnnames(table)), file_path(ref))
end

"""
    has_columns(ref::FileReference, cols::Symbol...) -> Bool

Check if all specified columns exist in the file.
"""
function has_columns(ref::FileReference, cols::Symbol...)
    available = column_names(ref)
    return all(col ∈ available for col in cols)
end

