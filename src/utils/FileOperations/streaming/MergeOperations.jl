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
High-performance merge operations with heap-based sorting.

Provides type-stable, memory-efficient merging of multiple
sorted files with support for arbitrary numbers of sort keys
and mixed sort directions.
"""

using Arrow, DataFrames, Tables
using DataStructures: BinaryMinHeap, BinaryMaxHeap, BinaryHeap
using Base.Order: Ordering, Forward, Reverse, By, Lt

#==========================================================
Helper Functions for Heap-based Merge
==========================================================#

# Create empty DataFrame with same schema as Arrow table
function _create_empty_dataframe(table::Arrow.Table, n_rows::Int)
    df = DataFrame()
    col_names = Tables.columnnames(table)
    schema_types = Tables.schema(table).types
    
    for (col_name, col_type) in zip(col_names, schema_types)
        # Handle Union types properly
        if col_type isa Union && Missing <: col_type
            # Extract the non-missing type
            non_missing_type = Base.uniontypes(col_type)
            actual_type = first(t for t in non_missing_type if t !== Missing)
            df[!, col_name] = Vector{Union{Missing, actual_type}}(undef, n_rows)
        else
            df[!, col_name] = Vector{col_type}(undef, n_rows)
        end
    end
    return df
end

# Fill batch columns from sorted (table, row) locations.
#
# Each column is looked up once per batch for every input table, then gathered in a function
# specialised on the concrete column types. Looking the column up by name inside the row loop
# (one Symbol lookup, an untyped read and a boxed value per cell) made the merge ~3x slower and
# allocate ~3x more.
function _fill_batch_columns_nkey!(
    batch_df::DataFrame,
    tables::Vector{Arrow.Table},
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    n_rows::Int
)
    for col_name in names(batch_df)
        col_symbol = Symbol(col_name)
        # Concretely typed when every input stores the column with the same Arrow array type.
        sources = map(table -> Tables.getcolumn(table, col_symbol), tables)
        _gather_column!(batch_df[!, col_symbol], sources, sorted_tuples, n_rows)
    end
    return batch_df
end

# Function barrier: one specialisation per (batch column, source column) type pair.
function _gather_column!(
    dest::AbstractVector,
    sources::AbstractVector,
    sorted_tuples::Vector{Tuple{Int64, Int64}},
    n_rows::Int
)
    _check_gather_types(dest, sources)
    @inbounds for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        dest[i] = sources[table_idx][row_idx]
    end
    return dest
end

# Keep a mismatched input an error, as the per-cell type assertions did, rather than letting
# assignment convert it (e.g. Float64 into a Float32 batch column). A nullable source may feed a
# non-nullable batch column: an actual `missing` value still fails on assignment.
function _check_gather_types(dest::AbstractVector, sources::AbstractVector)
    T = eltype(dest)
    for source in sources
        S = eltype(source)
        # List columns: Arrow yields views of the batch's array type.
        (nonmissingtype(S) <: nonmissingtype(T) || (T <: AbstractArray && S <: AbstractArray)) ||
            throw(ArgumentError("Merge input column type $S does not match the batch column type $T"))
    end
    return nothing
end

#==========================================================
Type-Stable DataFrame Creation
==========================================================#

"""
    unwrap_array_eltype(col_type)

Unwrap nested array types to get the scalar element type.
E.g., SubArray{Float32, ...} -> Float32, Vector{Int32} -> Int32
Strings are not unwrapped as they are valid scalar types.
"""
function unwrap_array_eltype(col_type)
    # Don't unwrap strings - they are valid scalar types
    if col_type <: AbstractString
        return col_type
    end

    # Unwrap array types to get the inner element type
    while col_type isa DataType && col_type <: AbstractArray
        col_type = eltype(col_type)
    end

    return col_type
end

"""
Create a type-stable empty DataFrame with pre-allocated vectors.
Types are determined from the schema of the first table.

Note: Array element types (e.g., SubArray{Float32, ...}) are preserved as-is,
since some Arrow columns store array data per row (List columns).
"""
function create_typed_dataframe(reference_table::Arrow.Table, batch_size::Int)
    df = DataFrame()
    col_names = Tables.columnnames(reference_table)

    for col_name in col_names
        col_type = eltype(Tables.getcolumn(reference_table, col_name))

        # Create appropriately typed vector (preserve array types as-is)
        if col_type isa Union && Missing <: col_type
            # Handle Union{Missing, T} types
            non_missing_types = Base.uniontypes(col_type)
            actual_type = first(t for t in non_missing_types if t !== Missing)
            df[!, col_name] = Vector{Union{Missing, actual_type}}(undef, batch_size)
        else
            df[!, col_name] = Vector{col_type}(undef, batch_size)
        end
    end

    return df
end

#==========================================================
N-Key Implementation Support Functions
==========================================================#

"""
Normalize reverse specification to vector format.
"""
function _normalize_reverse_spec(reverse::Union{Bool,Vector{Bool}}, n_keys::Int)
    if reverse isa Bool
        return fill(reverse, n_keys)
    elseif length(reverse) == 1
        return fill(reverse[1], n_keys)
    elseif length(reverse) == n_keys
        return reverse
    else
        error("reverse must be Bool, single-element vector, or have same length as sort keys ($n_keys)")
    end
end

"""
Compare two tuples with mixed reverse directions for each key.
"""
function _mixed_reverse_compare(a, b, reverse_vec::Vector{Bool})
    # Compare each sort key according to its reverse setting
    for i in 1:(length(a)-1)  # -1 to skip table_idx (last element)
        val_a, val_b = a[i], b[i]
        isequal(val_a, val_b) && continue

        # `isless` provides a total order for `missing`, unlike `<`/`>`, which
        # return `missing`. This keeps mixed-direction merges consistent with
        # Julia/DataFrames sorting: missing values are last for ascending keys
        # and first for descending keys.
        return reverse_vec[i] ? isless(val_b, val_a) : isless(val_a, val_b)
    end
    # All sort keys are equal, use table_idx as tiebreaker (always ascending)
    return a[end] < b[end]
end

"""
Generate heap type dynamically based on sort types and reverse specification.
Supports mixed reverse directions for different keys.
"""
function _create_typed_heap(::Type{SortTypes}, reverse_vec::Vector{Bool}) where SortTypes
    if length(reverse_vec) == 1
        # Single key - use simple min/max heap
        if reverse_vec[1]
            return BinaryMaxHeap{SortTypes}()
        else
            return BinaryMinHeap{SortTypes}()
        end
    elseif all(reverse_vec)
        # All keys reversed - use max heap
        return BinaryMaxHeap{SortTypes}()
    elseif !any(reverse_vec)
        # No keys reversed - use min heap
        return BinaryMinHeap{SortTypes}()
    else
        # Mixed reverse directions - use custom comparison with BinaryHeap
        # Create comparison function that returns true if a < b in our desired order
        comp_func = (a, b) -> _mixed_reverse_compare(a, b, reverse_vec)
        ordering = Lt(comp_func)
        return BinaryHeap(ordering, SortTypes[])
    end
end

"""
Add entry to N-key heap with variadic sort keys.
"""
function _add_to_nkey_heap!(
    heap::H,
    table::Arrow.Table,
    table_idx::Int,
    row_idx::Int,
    sort_keys::NTuple{N,Symbol}
) where {N, H<:Union{BinaryMinHeap, BinaryMaxHeap, BinaryHeap}}
    # Build tuple: (val1, val2, ..., valN, table_idx)
    values = tuple((Tables.getcolumn(table, key)[row_idx] for key in sort_keys)..., table_idx)
    push!(heap, values)
end

# A merge owns exactly one active output. IO ownership stays here rather than
# inside Arrow so even a failed writer finalization cannot leave the file open.
mutable struct _MergeOutput
    path::String
    temp_path::String
    io::IOStream
    writer::Arrow.Writer{IOStream}
    closed::Bool
    complete::Bool
    publish_attempted::Bool
    published::Bool
end

function _open_merge_output(path::String)
    temp_path, io = mktemp(dirname(abspath(path)); cleanup=false)
    try
        writer = open(Arrow.Writer, io; file=false, ntasks=0, closeio=false)
        # Arrow 2.x leaves its consumer task unbound. An IO failure must close
        # the channel and wake blocked producers, rather than hang the merge.
        Base.bind(writer.msgs, writer.task)
        return _MergeOutput(path, temp_path, io, writer, false, false, false, false)
    catch
        try
            close(io)
            safeRm(temp_path; force=true)
        catch cleanup_error
            @warn "Failed to clean up merge writer construction" path=temp_path exception=(cleanup_error, catch_backtrace())
        end
        rethrow()
    end
end

function _close_merge_output!(output::_MergeOutput)
    output.closed && return nothing
    # Do not retry Arrow finalization if it failed after closing its channel.
    output.closed = true
    io_was_open = isopen(output.io)
    try
        close(output.writer)
    finally
        close(output.io)
    end
    # Arrow 2.x's close returns normally when its consumer task has already failed, and
    # ignores the end-of-stream write to a closed IO, so neither failure reaches us on its
    # own. A stream missing batches must never be marked complete and published.
    io_was_open || error("Merge output IO was closed before finalization: $(output.temp_path)")
    istaskfailed(output.writer.task) && throw(TaskFailedException(output.writer.task))
    output.complete = true
    return nothing
end

function _publish_merge_output!(output::_MergeOutput)
    _close_merge_output!(output)
    output.complete || error("Cannot publish an incomplete merge output: $(output.temp_path)")
    output.publish_attempted = true
    # mv(...; force=true) can recursively remove an existing directory on Unix.
    # A merge may replace a file, never a directory or a link to a directory.
    isdir(output.path) && error("Merge destination is a directory: $(output.path)")
    # Windows cannot replace a file while we still own its writer handle.
    # safeRm handles old, unreachable mappings using its serialized GC fallback.
    if Sys.iswindows()
        safeRm(output.path; force=true)
        ispath(output.path) && error("Merge destination is still present: $(output.path)")
    end
    mv(output.temp_path, output.path; force=!Sys.iswindows())
    output.published = true
    return nothing
end

function _cleanup_merge_output!(output::_MergeOutput)
    output.published && return nothing
    try
        _close_merge_output!(output)
    catch error
        @warn "Failed to finalize merge output during cleanup" path=output.temp_path exception=(error, catch_backtrace())
    end
    if output.publish_attempted && output.complete
        # Keep a completed stream recoverable if replacement/move failed.
        @warn "Merge publication failed; completed temporary output retained" path=output.temp_path destination=output.path
    else
        try
            safeRm(output.temp_path; force=true)
        catch error
            @warn "Failed to remove temporary merge output" path=output.temp_path exception=(error, catch_backtrace())
        end
    end
    return nothing
end

function _with_merge_output(f, path::String)
    output = _open_merge_output(path)
    try
        f(output)
        _publish_merge_output!(output)
    finally
        # Preserve the original merge/publication exception if cleanup fails.
        _cleanup_merge_output!(output)
    end
    return nothing
end

"""
Submit filled rows to a persistent stream writer using owned column vectors.

Arrow's unbuffered channel bounds pending messages, but the consumer can still
be writing when Arrow.write returns. The merge may immediately refill batch_df,
so even a full batch needs a snapshot. Nested Arrow values refer to read-only
source tables; the merge only replaces outer vector entries, never their contents.

Arrow 2.8 retains per-record-batch block metadata even in stream mode. This
removes quadratic append scans, but does not claim constant writer metadata.
"""
function _write_batch_typed(output::_MergeOutput, batch_df::DataFrame, n_rows::Int)
    Arrow.write(output.writer, batch_df[1:n_rows, :])
    return nothing
end

function _validate_merge_destination(path::String, refs::Vector{<:FileReference})
    destination = abspath(normpath(path))
    for ref in refs
        source = abspath(normpath(file_path(ref)))
        same_path = Sys.iswindows() ? lowercase(destination) == lowercase(source) : destination == source
        same_file = ispath(destination) && ispath(source) && Base.Filesystem.samefile(destination, source)
        (same_path || same_file) && throw(ArgumentError("Merge output aliases an input file: $path"))
    end
    return nothing
end

#==========================================================
Hierarchical Merge Support (FD-safe staging)
==========================================================#

# Each concurrent staging merge holds its own output batch (up to `batch_size` rows), so concurrency is capped.
const MAX_CONCURRENT_STAGE_MERGES = 8

"""
    _stage_merge(refs, sort_keys...; max_fanin, reverse, batch_size)

Merge `refs` in groups of `max_fanin` into temporary Arrow files.
Returns a vector of FileReferences pointing to the staged temp files, in group order.
Used internally to keep the number of inputs per merge bounded.

The groups are independent, so up to `MAX_CONCURRENT_STAGE_MERGES` of them are merged at once. The result is
identical to merging them one after another: each staged file keeps its group's position, and ties within a merge
are broken by input order.
"""
function _stage_merge(
    refs::Vector{<:FileReference},
    sort_keys::Symbol...;
    max_fanin::Int,
    reverse::Union{Bool,Vector{Bool}},
    batch_size::Int
)
    started = time()
    groups = collect(Iterators.partition(refs, max_fanin))
    n_groups = length(groups)
    n_workers = min(n_groups, Threads.nthreads(), MAX_CONCURRENT_STAGE_MERGES)
    @debug_l1 "Staged merge starting: keys=$(join(sort_keys, ',')) files=$(length(refs)) groups=$n_groups max_fanin=$max_fanin workers=$n_workers"
    temp_dir = mktempdir()
    staged_refs = similar(refs, n_groups)
    # Concurrent merges share the output-batch budget, so buffered rows stay what one serial merge would hold.
    worker_batch_size = cld(batch_size, n_workers)
    next_group = Threads.Atomic{Int}(1)
    @sync for _ in 1:n_workers
        Threads.@spawn while true
            i = Threads.atomic_add!(next_group, 1)
            i > n_groups && break
            staged_refs[i] = stream_sorted_merge(
                collect(groups[i]), joinpath(temp_dir, "stage_$i.arrow"), sort_keys...;
                reverse, batch_size = worker_batch_size, max_fanin
            )
        end
    end
    rows_processed = sum(row_count, staged_refs; init = 0)
    @debug_l1 "Staged merge complete: keys=$(join(sort_keys, ',')) files=$(length(refs)) groups=$n_groups rows=$rows_processed elapsed=$(round(time() - started, digits=2))s"
    return staged_refs
end

#==========================================================
Main Stream Sorted Merge Function
==========================================================#

"""
    stream_sorted_merge(refs::Vector{<:FileReference}, 
                       output_path::String,
                       sort_keys::Symbol...;
                       reverse::Union{Bool,Vector{Bool}}=false,
                       batch_size::Int=1_000_000)

High-performance type-stable merge function that handles arbitrary numbers of sort keys.
Determines types at compile time and uses specialized operations for maximum performance.

Provides 4-20x speedup over the original implementation through:
- Type-stable heap operations with compile-time optimizations
- Specialized column filling methods for different data types  
- Memory-efficient batch processing
- Support for arbitrary number of sort keys (2, 3, 4, 5+)

# Arguments
- `refs`: Vector of file references to merge
- `output_path`: Path for merged output file
- `sort_keys...`: Variable number of sort key symbols
- `reverse`: Boolean or vector specifying sort direction(s)
- `batch_size`: Number of rows to process in each batch

# Examples
```julia
# 2 keys (protein groups)
stream_sorted_merge(refs, path, :pg_score, :target; reverse=[true, true])

# 4 keys (MaxLFQ scenario)  
stream_sorted_merge(refs, path, :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx; reverse=true)

# Mixed reverse directions
stream_sorted_merge(refs, path, :protein, :score, :name; reverse=[false, true, false])
```
"""
function stream_sorted_merge(
    refs::Vector{<:FileReference},
    output_path::String,
    sort_keys::Symbol...;
    reverse::Union{Bool,Vector{Bool}}=false,
    batch_size::Int=1_000_000,
    max_fanin::Int=64
)
    # Validate inputs
    isempty(refs) && error("No files to merge")
    isempty(sort_keys) && error("At least one sort key must be specified")
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    max_fanin >= 2 || throw(ArgumentError("max_fanin must be at least 2"))
    _validate_merge_destination(output_path, refs)

    # Hierarchical merge: stage in batches to avoid FD exhaustion
    if length(refs) > max_fanin
        staged_refs = _stage_merge(refs, sort_keys...; max_fanin, reverse, batch_size)
        return stream_sorted_merge(
            staged_refs, output_path, sort_keys...;
            reverse, batch_size, max_fanin
        )
    end

    # Convert sort_keys to tuple for type stability
    sort_keys_tuple = tuple(sort_keys...)

    # Normalize reverse specification
    reverse_vec = _normalize_reverse_spec(reverse, length(sort_keys))

    # Determine types from first file
    sort_types = with_arrow_table(file_path(first(refs))) do first_table
        tuple((eltype(Tables.getcolumn(first_table, key)) for key in sort_keys_tuple)...)
    end

    # Dispatch to N-key implementation. The sources are unmapped once the output is
    # written, so callers can delete or replace them straight away.
    return with_arrow_tables() do open_table
        _stream_sorted_merge_nkey_impl(
            refs, output_path, sort_keys_tuple, sort_types, reverse_vec, batch_size, open_table
        )
    end
end

"""
Internal N-key implementation with dynamic type dispatch.
"""
function _stream_sorted_merge_nkey_impl(
    refs::Vector{<:FileReference},
    output_path::String,
    sort_keys::NTuple{N,Symbol},
    sort_types::NTuple{M,Type},
    reverse_vec::Vector{Bool},
    batch_size::Int,
    open_table = Arrow.Table
) where {N, M}
    started = last_progress = time()
    # Validate all files exist and have compatible schemas
    for ref in refs
        validate_exists(ref)
    end

    # Validate that all files are sorted by the requested keys
    for ref in refs
        if !is_sorted_by(ref, sort_keys...)
            error("File $(file_path(ref)) is not sorted by the required keys: $(sort_keys). Use sort_file_by_keys! or mark_sorted! first.")
        end
    end

    # Load all tables
    tables = [open_table(file_path(ref)) for ref in refs]
    
    # Validate that all tables have the required sort columns
    for (i, table) in enumerate(tables)
        available_columns = Set(Tables.columnnames(table))
        for key in sort_keys
            if key ∉ available_columns
                throw(BoundsError("Column $key not found in file $(file_path(refs[i])). Available columns: $(collect(available_columns))"))
            end
        end
    end
    
    table_sizes = [length(Tables.getcolumn(table, 1)) for table in tables]
    total_rows = sum(table_sizes)
    table_indices = ones(Int64, length(tables))
    
    # Create type-stable batch DataFrame
    batch_df = create_typed_dataframe(first(tables), batch_size)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, batch_size)
    
    # Create heap with full type information
    heap_tuple_type = Tuple{sort_types..., Int64}
    heap = _create_typed_heap(heap_tuple_type, reverse_vec)
    
    # Initialize heap with first row from each table
    for (i, table) in enumerate(tables)
        if table_sizes[i] > 0
            _add_to_nkey_heap!(heap, table, i, 1, sort_keys)
        end
    end
    
    _with_merge_output(output_path) do output
        row_idx = 1
        n_writes = 0

        while !isempty(heap)
            heap_entry = pop!(heap)
            table_idx = heap_entry[end]
            current_row_idx = table_indices[table_idx]
            sorted_tuples[row_idx] = (table_idx, current_row_idx)

            table_indices[table_idx] += 1
            next_row_idx = table_indices[table_idx]
            if next_row_idx <= table_sizes[table_idx]
                _add_to_nkey_heap!(heap, tables[table_idx], table_idx, next_row_idx, sort_keys)
            end

            row_idx += 1
            if row_idx > batch_size
                _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, batch_size)
                _write_batch_typed(output, batch_df, batch_size)
                n_writes += 1
                row_idx = 1
                if time() - last_progress >= 60
                    @debug_l1 "Sorted merge ($(basename(output_path))): keys=$(join(sort_keys, ',')) files=$(length(refs)) rows=$(n_writes * batch_size)/$total_rows elapsed=$(round(time() - started, digits=2))s"
                    last_progress = time()
                end
            end
        end

        if row_idx > 1
            final_rows = row_idx - 1
            _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, final_rows)
            _write_batch_typed(output, batch_df, final_rows)
        elseif n_writes == 0
            # An all-empty merge still publishes a valid stream with its schema.
            _write_batch_typed(output, batch_df, 0)
        end
    end
    
    # Create output reference
    output_ref = create_reference(output_path, typeof(first(refs)))

    # Mark as sorted - heap-based merge maintains sort order for all cases
    mark_sorted!(output_ref, sort_keys...)

    if time() - started >= 60
        @debug_l1 "Sorted merge ($(basename(output_path))) complete: keys=$(join(sort_keys, ',')) files=$(length(refs)) rows=$total_rows elapsed=$(round(time() - started, digits=2))s"
    end
    return output_ref
end

"""
    stream_sorted_merge_chunked(refs, output_dir, group_key, sort_keys...;
                                reverse, batch_size, max_chunk_bytes)

Like `stream_sorted_merge` but splits the merged output into multiple chunk
files using `max_chunk_bytes` as an estimated byte target. The estimate uses
source bytes per row, avoiding probes of an asynchronously written output.
Chunk boundaries are placed only at `group_key` transitions so every chunk
contains only complete groups (when the group key is a leading sort key).
Each chunk always gets at least one complete group even if it exceeds the limit.

Returns `Vector{<:FileReference}` — one per chunk, each marked as sorted.
"""
function stream_sorted_merge_chunked(
    refs::Vector{<:FileReference},
    output_dir::String,
    group_key::Symbol,
    sort_keys::Symbol...;
    reverse::Union{Bool,Vector{Bool}}=false,
    batch_size::Int=1_000_000,
    max_chunk_bytes::Int=1_000_000_000,
    max_fanin::Int=64
)
    isempty(refs) && error("No files to merge")
    isempty(sort_keys) && error("At least one sort key must be specified")
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    max_fanin >= 2 || throw(ArgumentError("max_fanin must be at least 2"))
    max_chunk_bytes > 0 || throw(ArgumentError("max_chunk_bytes must be positive"))

    # Validate prospective chunk destinations before staging loses the original
    # input references. Later chunks are validated again before opening them.
    for ref in refs
        source_path = file_path(ref)
        name = basename(ispath(source_path) ? realpath(source_path) : source_path)
        candidate_name = Sys.iswindows() ? lowercase(name) : name
        if occursin(r"^chunk_[0-9]{4,}\.arrow$", candidate_name)
            _validate_merge_destination(joinpath(output_dir, name), [ref])
        end
    end

    # Hierarchical merge: stage in batches to avoid FD exhaustion
    if length(refs) > max_fanin
        staged_refs = _stage_merge(refs, sort_keys...; max_fanin, reverse, batch_size)
        return stream_sorted_merge_chunked(
            staged_refs, output_dir, group_key, sort_keys...;
            reverse, batch_size, max_chunk_bytes, max_fanin
        )
    end

    sort_keys_tuple = tuple(sort_keys...)
    reverse_vec = _normalize_reverse_spec(reverse, length(sort_keys))

    sort_types = with_arrow_table(file_path(first(refs))) do first_table
        tuple((eltype(Tables.getcolumn(first_table, key)) for key in sort_keys_tuple)...)
    end

    # The sources are unmapped once every chunk is written.
    return with_arrow_tables() do open_table
        _stream_sorted_merge_chunked_impl(
            refs, output_dir, group_key, sort_keys_tuple, sort_types,
            reverse_vec, batch_size, max_chunk_bytes, open_table
        )
    end
end

function _stream_sorted_merge_chunked_impl(
    refs::Vector{<:FileReference},
    output_dir::String,
    group_key::Symbol,
    sort_keys::NTuple{N,Symbol},
    sort_types::NTuple{M,Type},
    reverse_vec::Vector{Bool},
    batch_size::Int,
    max_chunk_bytes::Int,
    open_table = Arrow.Table
) where {N, M}
    # Validate
    for ref in refs
        validate_exists(ref)
        if !is_sorted_by(ref, sort_keys...)
            error("File $(file_path(ref)) is not sorted by the required keys: $(sort_keys)")
        end
    end

    tables = [open_table(file_path(ref)) for ref in refs]
    for (i, table) in enumerate(tables)
        available_columns = Set(Tables.columnnames(table))
        for key in sort_keys
            key ∉ available_columns && throw(BoundsError("Column $key not found in file $(file_path(refs[i]))"))
        end
        group_key ∉ available_columns && throw(BoundsError("Group key $group_key not found in file $(file_path(refs[i]))"))
    end

    table_sizes = [length(Tables.getcolumn(table, 1)) for table in tables]
    total_source_bytes = sum(filesize(file_path(ref)) for ref in refs)
    total_source_rows = sum(table_sizes)
    estimated_bytes_per_row = total_source_rows > 0 ? total_source_bytes / total_source_rows : 0.0
    table_indices = ones(Int64, length(tables))
    # Looked up once per merge, not by name for every row.
    group_columns = map(table -> Tables.getcolumn(table, group_key), tables)

    batch_df = create_typed_dataframe(first(tables), batch_size)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, batch_size)

    heap_tuple_type = Tuple{sort_types..., Int64}
    heap = _create_typed_heap(heap_tuple_type, reverse_vec)

    for (i, table) in enumerate(tables)
        if table_sizes[i] > 0
            _add_to_nkey_heap!(heap, table, i, 1, sort_keys)
        end
    end

    # Chunking state
    mkpath(output_dir)
    chunk_paths = String[]
    chunk_idx = 1
    chunk_path() = joinpath(output_dir, @sprintf("chunk_%04d.arrow", chunk_idx))

    row_idx = 1
    n_writes_in_chunk = 0
    current_chunk_bytes = 0.0
    prev_group = nothing      # group_key value of the most recently queued row
    rows_processed = 0
    rows_with_missing_group = 0
    output = nothing

    try
        while !isempty(heap)
            heap_entry = pop!(heap)
            src_table_idx = heap_entry[end]
            src_row_idx = table_indices[src_table_idx]

            row_group = group_columns[src_table_idx][src_row_idx]
            if ismissing(row_group)
                rows_with_missing_group += 1
            end

            # The byte target is soft: a complete group is never split. Keep
            # source-based estimates rather than inspecting a file whose Arrow
            # consumer may still be writing the most recently submitted batch.
            if prev_group !== nothing && !isequal(row_group, prev_group) &&
               current_chunk_bytes >= max_chunk_bytes && (n_writes_in_chunk > 0 || row_idx > 1)
                if output === nothing
                    _validate_merge_destination(chunk_path(), refs)
                    output = _open_merge_output(chunk_path())
                end
                if row_idx > 1
                    pending = row_idx - 1
                    _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, pending)
                    _write_batch_typed(output, batch_df, pending)
                    row_idx = 1
                end
                _publish_merge_output!(output)
                push!(chunk_paths, chunk_path())
                output = nothing
                chunk_idx += 1
                n_writes_in_chunk = 0
                current_chunk_bytes = 0.0
            end

            sorted_tuples[row_idx] = (src_table_idx, src_row_idx)
            prev_group = row_group
            rows_processed += 1
            current_chunk_bytes += estimated_bytes_per_row

            table_indices[src_table_idx] += 1
            next_src_row = table_indices[src_table_idx]
            if next_src_row <= table_sizes[src_table_idx]
                _add_to_nkey_heap!(heap, tables[src_table_idx], src_table_idx, next_src_row, sort_keys)
            end

            row_idx += 1
            if row_idx > batch_size
                _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, batch_size)
                if output === nothing
                    _validate_merge_destination(chunk_path(), refs)
                    output = _open_merge_output(chunk_path())
                end
                _write_batch_typed(output, batch_df, batch_size)
                n_writes_in_chunk += 1
                row_idx = 1
            end
        end

        if output === nothing
            _validate_merge_destination(chunk_path(), refs)
            output = _open_merge_output(chunk_path())
        end
        if row_idx > 1
            final_rows = row_idx - 1
            _fill_batch_columns_nkey!(batch_df, tables, sorted_tuples, final_rows)
            _write_batch_typed(output, batch_df, final_rows)
        elseif n_writes_in_chunk == 0
            _write_batch_typed(output, batch_df, 0)
        end
        _publish_merge_output!(output)
        push!(chunk_paths, chunk_path())
        output = nothing
    catch
        if !isempty(chunk_paths)
            @warn "Chunked merge failed after publishing completed chunks" directory=abspath(output_dir) completed_chunks=length(chunk_paths)
        end
        rethrow()
    finally
        output === nothing || _cleanup_merge_output!(output)
    end

    if rows_with_missing_group > 0
        # Expected, not a fault: the message itself says these rows are filtered downstream. It fired on
        # every run, so it belongs in the debug log rather than the warnings sidecar.
        @debug_l1 "Chunked merge: $rows_with_missing_group / $rows_processed rows had missing $group_key (shared peptides — filtered downstream by LFQ)"
    end

    # Create output references
    ref_type = typeof(first(refs))
    chunk_refs = map(chunk_paths) do path
        ref = create_reference(path, ref_type)
        mark_sorted!(ref, sort_keys...)
        ref
    end

    @debug_l1 "Chunked merge: $(length(chunk_refs)) chunk(s), $rows_processed total rows"
    return chunk_refs
end
