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

"""
Arrow-specific file operations and utilities.

Provides specialized operations for Arrow files including
sorting, writing, and file management with Windows compatibility.
"""

using Arrow, DataFrames, Tables

#==========================================================
Arrow File Operations
==========================================================#

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
    compute_sortperm(ref::FileReference, sort_keys::Symbol...;
                     reverse::Union{Bool, Vector{Bool}}=false) -> Vector{Int32}

Compute a permutation that would sort `ref` by `sort_keys` without rewriting
the file. Only the sort-key columns are read from disk. The returned vector
can be applied to one or more files (e.g. main + sidecars) via
`apply_sortperm!` so they end up in lockstep order without each being sorted
independently.

`reverse` accepts a single Bool (applied to all keys) or a Vector{Bool} (one
per key). Returns a `Vector{Int32}` permutation (sufficient for ≤ 2^31 rows).
"""
function compute_sortperm(ref::FileReference, sort_keys::Symbol...;
                          reverse::Union{Bool, Vector{Bool}}=false)
    validate_exists(ref)
    for k in sort_keys
        has_column(schema(ref), k) || error("Sort key $k not found in file schema")
    end
    rev_vec = reverse isa Bool ? fill(reverse, length(sort_keys)) : reverse
    length(rev_vec) == length(sort_keys) ||
        error("Length of reverse vector ($(length(rev_vec))) must match number of sort keys ($(length(sort_keys)))")

    tbl = Arrow.Table(file_path(ref))
    key_cols = Any[collect(Tables.getcolumn(tbl, k)) for k in sort_keys]
    n = isempty(key_cols) ? 0 : length(key_cols[1])
    perm = collect(Int32(1):Int32(n))
    nkeys = length(key_cols)

    @inline function key_lt(i::Int32, j::Int32)
        @inbounds for c in 1:nkeys
            a = key_cols[c][i]; b = key_cols[c][j]
            if !isequal(a, b)
                return rev_vec[c] ? isless(b, a) : isless(a, b)
            end
        end
        return false
    end
    sort!(perm; lt=key_lt)
    return perm
end

"""
    apply_sortperm!(ref::FileReference, perm::AbstractVector{<:Integer};
                    mark_sort_keys::Tuple{Vararg{Symbol}}=()) -> FileReference

Apply `perm` to every column of the file referenced by `ref`. The file is
read, permuted, and written back via the Windows-safe `writeArrow` path.
If `mark_sort_keys` is non-empty, the reference's `sorted_by` metadata is
updated to that tuple.

The intended use is to compute one permutation via `compute_sortperm` and
then apply it to the main file plus any registered sidecars so they all share
the same row order.
"""
function apply_sortperm!(ref::FileReference, perm::AbstractVector{<:Integer};
                         mark_sort_keys::Tuple{Vararg{Symbol}}=())
    validate_exists(ref)
    _permute_arrow_file_inplace!(file_path(ref), perm)
    if !isempty(mark_sort_keys)
        mark_sorted!(ref, mark_sort_keys...)
    else
        ref.sorted_by = ()
    end
    return ref
end

"""
    apply_sortperm!(ref::PSMFileReference, perm; mark_sort_keys=()) -> PSMFileReference

PSM-specialization: permutes the main file AND every registered sidecar in
lockstep. The row-aligned invariant between main and sidecars is preserved.
"""
function apply_sortperm!(ref::PSMFileReference, perm::AbstractVector{<:Integer};
                         mark_sort_keys::Tuple{Vararg{Symbol}}=())
    validate_exists(ref)
    _permute_arrow_file_inplace!(file_path(ref), perm)
    for s in ref.sidecars
        _permute_arrow_file_inplace!(s.path, perm)
    end
    if !isempty(mark_sort_keys)
        mark_sorted!(ref, mark_sort_keys...)
    else
        ref.sorted_by = ()
    end
    return ref
end

# Internal: permute every column of an Arrow file in place. Used by the
# FileReference and PSMFileReference apply_sortperm! methods.
function _permute_arrow_file_inplace!(path::String, perm::AbstractVector{<:Integer})
    df = DataFrame(Tables.columntable(Arrow.Table(path)))
    n = nrow(df)
    length(perm) == n ||
        error("perm length $(length(perm)) ≠ file row count $n for $path")
    df_sorted = df[perm, :]
    writeArrow(path, df_sorted)
    return path
end

"""
    apply_sortperm!(refs::Vector{<:FileReference}, perm::AbstractVector{<:Integer};
                    mark_sort_keys::Tuple{Vararg{Symbol}}=(), parallel::Bool=true)

Apply the same permutation to multiple files (e.g. main + sidecars). All
files must share the same row count and original row order. Optionally
parallelized across files.
"""
function apply_sortperm!(refs::Vector{<:FileReference}, perm::AbstractVector{<:Integer};
                         mark_sort_keys::Tuple{Vararg{Symbol}}=(), parallel::Bool=true)
    if parallel && length(refs) > 1
        Threads.@threads for ref in refs
            exists(ref) && apply_sortperm!(ref, perm; mark_sort_keys=mark_sort_keys)
        end
    else
        for ref in refs
            exists(ref) && apply_sortperm!(ref, perm; mark_sort_keys=mark_sort_keys)
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
    table = Arrow.Table(file_path(ref))
    return Symbol.(Tables.columnnames(table))
end

"""
    has_columns(ref::FileReference, cols::Symbol...) -> Bool

Check if all specified columns exist in the file.
"""
function has_columns(ref::FileReference, cols::Symbol...)
    available = column_names(ref)
    return all(col ∈ available for col in cols)
end

