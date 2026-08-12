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
File reference types for managing Arrow files with metadata tracking.

These types provide lightweight references to files on disk with:
- Schema tracking (immutable column information)
- Sort state tracking (which keys the file is sorted by)
- Validation methods to ensure data integrity
"""

using Arrow, DataFrames, Tables

#==========================================================
Abstract Types
==========================================================#

"""
    FileReference

Abstract base type for all file references.
Provides common interface for file metadata and operations.
"""
abstract type FileReference end

# Common interface that all FileReferences must implement
file_path(ref::FileReference) = ref.file_path
schema(ref::FileReference) = ref.schema
sorted_by(ref::FileReference) = ref.sorted_by
row_count(ref::FileReference) = ref.row_count
exists(ref::FileReference) = ref.file_exists

#==========================================================
File Reference Types
==========================================================#

"""
    Sidecar

A narrow auxiliary file holding a few columns aligned by row position with
its parent main file. Row positions in main and every sidecar must match
exactly at all times. Use `register_sidecar!` to attach one to a
`PSMFileReference`; permutation-aware ops like `apply_sortperm!` walk all
sidecars together so the alignment is preserved.
"""
struct Sidecar
    path::String
    cols::Vector{Symbol}
end

"""
    PSMFileReference <: FileReference

Reference to a PSM (Peptide-Spectrum Match) Arrow file with metadata.
Tracks schema, sort state, basic statistics, and any registered sidecars.
"""
mutable struct PSMFileReference <: FileReference
    file_path::String
    schema::FileSchema            # Immutable schema for the main file
    sorted_by::Tuple{Vararg{Symbol}}  # Which keys main+sidecars are sorted by
    row_count::Int64
    file_exists::Bool
    sidecars::Vector{Sidecar}     # Row-aligned auxiliary files (added Phase 2)

    # Constructor that validates file and extracts metadata. Any on-disk
    # sidecar files matching `{file_path}.<tag>.sidecar.arrow` and whose
    # row count matches the main file are auto-discovered and registered.
    # This ensures fresh PSMFileReferences built from a path inherit any
    # sidecars produced upstream by `add_columns_via_sidecar!`.
    function PSMFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false, Sidecar[])
        end

        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))

        ref = new(file_path, schema, (), row_count, true, Sidecar[])
        _discover_sidecars_from_disk!(ref)
        return ref
    end
end

# Scan the directory containing `ref.file_path` for files named
# "{basename}.<tag>.sidecar.arrow" and register any that have the matching
# row count + no schema collisions. Used by the PSMFileReference constructor.
function _discover_sidecars_from_disk!(ref::PSMFileReference)
    dir = dirname(ref.file_path)
    base = basename(ref.file_path)
    isdir(dir) || return ref
    pattern = base * "."
    for entry in readdir(dir)
        startswith(entry, pattern) || continue
        endswith(entry, ".sidecar.arrow") || continue
        side_path = joinpath(dir, entry)
        # Skip self-collisions: if the main path itself has a .sidecar.arrow
        # tail (shouldn't, but be defensive), don't register it as its own sidecar.
        side_path == ref.file_path && continue
        try
            tbl = Arrow.Table(side_path)
            side_cols = Symbol.(Tables.columnnames(tbl))
            n = isempty(side_cols) ? 0 : length(Tables.getcolumn(tbl, 1))
            n == ref.row_count || continue
            # Filter out cols that already exist in main or another sidecar.
            # Discovery is opportunistic — colliding cols are just skipped
            # (the file may legitimately exist with extra cols we already
            # have). Only register the non-colliding ones.
            keep_cols = Symbol[c for c in side_cols if !has_column(ref.schema, c) &&
                               !any(s -> c in s.cols, ref.sidecars)]
            isempty(keep_cols) && continue
            push!(ref.sidecars, Sidecar(side_path, keep_cols))
        catch
            # Skip unreadable / malformed sidecar candidates
            continue
        end
    end
    return ref
end

"""
    register_sidecar!(ref::PSMFileReference, path::String, cols::Vector{Symbol})

Attach a row-aligned sidecar file holding `cols` to `ref`. The file at `path`
must have the same row count as the main file and the same row order at the
time of registration (sidecars stay aligned because `apply_sortperm!` walks
them in lockstep with the main file).

Errors if `path` does not exist, if its row count disagrees with the main
file, or if any of `cols` overlaps with the main file's schema or an
already-registered sidecar.
"""
function register_sidecar!(ref::PSMFileReference, path::String, cols::Vector{Symbol})
    validate_exists(ref)
    isfile(path) || error("Sidecar file does not exist: $path")
    tbl = Arrow.Table(path)
    side_cols = Symbol.(Tables.columnnames(tbl))
    n_side = isempty(side_cols) ? 0 : length(Tables.getcolumn(tbl, 1))
    n_side == ref.row_count ||
        error("Sidecar row count $n_side ≠ main row count $(ref.row_count) for $path")
    for c in cols
        c in side_cols ||
            error("Sidecar at $path is missing declared column $c")
        has_column(ref.schema, c) &&
            error("Sidecar column $c already exists in main file schema")
        for s in ref.sidecars
            c in s.cols && error("Sidecar column $c already registered on another sidecar ($(s.path))")
        end
    end
    push!(ref.sidecars, Sidecar(path, copy(cols)))
    return ref
end

"""
    clear_sidecars!(ref::PSMFileReference; delete_files::Bool=false)

Detach all registered sidecars from `ref`. If `delete_files=true` (and the
files are safe to delete — typically after consolidation), the sidecar files
on disk are removed via `safeRm`.
"""
function clear_sidecars!(ref::PSMFileReference; delete_files::Bool=false)
    if delete_files
        for s in ref.sidecars
            isfile(s.path) && safeRm(s.path)
        end
    end
    empty!(ref.sidecars)
    return ref
end

"""
    sidecar_paths(ref::PSMFileReference) -> Vector{String}

Return all registered sidecar file paths.
"""
sidecar_paths(ref::PSMFileReference) = String[s.path for s in ref.sidecars]

"""
    has_column_anywhere(ref::PSMFileReference, col::Symbol) -> Bool

True iff `col` is in the main file's schema OR in any registered sidecar.
"""
function has_column_anywhere(ref::PSMFileReference, col::Symbol)
    has_column(ref.schema, col) && return true
    for s in ref.sidecars
        col in s.cols && return true
    end
    return false
end

"""
    locate_column(ref::PSMFileReference, col::Symbol) -> String

Return the path (main or sidecar) that contains `col`. Errors if not found.
"""
function locate_column(ref::PSMFileReference, col::Symbol)
    has_column(ref.schema, col) && return ref.file_path
    for s in ref.sidecars
        col in s.cols && return s.path
    end
    error("Column $col not found in main file or any sidecar of $(ref.file_path)")
end

"""
    load_with_sidecars(ref::FileReference) -> DataFrame

Read every column from the main file (and from any sidecars registered on a
`PSMFileReference`) into a single DataFrame, joined by row position. For
plain `FileReference` (no sidecars), behaves identically to
`DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))`.
"""
function load_with_sidecars(ref::FileReference)
    # DataFrame's default copycols=true gives us mutable Vector{T} columns
    # (downstream pipeline ops like filter/empty! require this).
    return DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
end

function load_with_sidecars(ref::PSMFileReference)
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    for s in ref.sidecars
        side = Arrow.Table(s.path)
        for c in s.cols
            # collect ensures the column is a mutable Vector{T}, matching the
            # mutability invariant used by the rest of the pipeline (filter!,
            # empty!, transform!).
            df[!, c] = collect(Tables.getcolumn(side, c))
        end
    end
    return df
end

"""
    add_columns_via_sidecar!(ref::PSMFileReference, name_to_col::Pair{Symbol, <:AbstractVector}...;
                             tag::AbstractString="cols") -> PSMFileReference

Write the supplied column(s) to a row-aligned sidecar Arrow file and
register them on `ref`. Avoids rewriting the main file when the only
operation is to append new columns. All columns must have length equal to
`ref.row_count`. `tag` is used to disambiguate the sidecar filename.
"""
function add_columns_via_sidecar!(ref::PSMFileReference,
                                  name_to_col::Pair{Symbol, <:AbstractVector}...;
                                  tag::AbstractString="cols")
    validate_exists(ref)
    isempty(name_to_col) && return ref
    n = ref.row_count
    cols = Symbol[]
    pairs = Pair[]
    for (k, v) in name_to_col
        length(v) == n ||
            error("Column $k has length $(length(v)); expected $n (main row count)")
        push!(cols, k)
        push!(pairs, k => v)
    end
    side_path = file_path(ref) * "." * tag * ".sidecar.arrow"
    side_df = DataFrame(pairs...)
    writeArrow(side_path, side_df)
    register_sidecar!(ref, side_path, cols)
    return ref
end

"""
    materialize_columns(ref::PSMFileReference, cols::Vector{Symbol}) -> DataFrame

Read just `cols` from wherever they live (main or sidecars) and return a
DataFrame. Used by readers that want to consume specific columns without
caring about the sidecar layout. Sidecars are aligned by row position with
main, so no join key is needed.
"""
function materialize_columns(ref::PSMFileReference, cols::Vector{Symbol})
    validate_exists(ref)
    # Group requested cols by source path so we only read each file once
    cols_per_path = Dict{String, Vector{Symbol}}()
    for c in cols
        push!(get!(cols_per_path, locate_column(ref, c), Symbol[]), c)
    end
    df = DataFrame()
    for (path, ps) in cols_per_path
        tbl = Arrow.Table(path)
        for c in ps
            df[!, c] = collect(Tables.getcolumn(tbl, c))
        end
    end
    return df[!, cols]  # preserve caller's column order
end

"""
    ProteinGroupFileReference <: FileReference

Reference to a protein group Arrow file with metadata.
Tracks schema, sort state, and basic statistics.
"""
mutable struct ProteinGroupFileReference <: FileReference
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
    
    # Constructor that validates file and extracts metadata
    function ProteinGroupFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false)
        end
        
        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))
        
        new(file_path, schema, (), row_count, true)
    end
end

"""
    ProteinQuantFileReference <: FileReference

Reference to a protein quantification Arrow file with metadata.
Specialized for MaxLFQ quantification results with additional metadata.
"""
mutable struct ProteinQuantFileReference <: FileReference
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
    n_protein_groups::Int64
    n_experiments::Int64
    
    # Constructor that validates file and extracts metadata
    function ProteinQuantFileReference(file_path::String)
        if !isfile(file_path)
            return new(file_path, FileSchema(Symbol[]), (), 0, false, 0, 0)
        end
        
        # Read schema from file
        tbl = Arrow.Table(file_path)
        schema = FileSchema(collect(Symbol.(Tables.columnnames(tbl))))
        # Get row count from first column if columns exist
        col_names = Tables.columnnames(tbl)
        row_count = isempty(col_names) ? 0 : length(Tables.getcolumn(tbl, 1))
        
        # Estimate protein groups and experiments for MaxLFQ context
        n_protein_groups = 0
        n_experiments = 0
        
        # Try to extract meaningful counts if the file has appropriate columns
        if has_column(schema, :protein)
            protein_col = Tables.getcolumn(tbl, :protein)
            n_protein_groups = length(unique(skipmissing(protein_col)))
        end
        
        if has_column(schema, :experiments) || has_column(schema, :file_name)
            exp_col = has_column(schema, :experiments) ? 
                     Tables.getcolumn(tbl, :experiments) : 
                     Tables.getcolumn(tbl, :file_name)
            n_experiments = length(unique(skipmissing(exp_col)))
        end
        
        new(file_path, schema, (), row_count, true, n_protein_groups, n_experiments)
    end
end

# Specialized accessors for ProteinQuantFileReference
n_protein_groups(ref::ProteinQuantFileReference) = ref.n_protein_groups
n_experiments(ref::ProteinQuantFileReference) = ref.n_experiments

#==========================================================
Paired File References
==========================================================#

"""
    PairedSearchFiles

Enforces 1:1 pairing between PSM and protein group files.
Ensures consistency between related files.
"""
struct PairedSearchFiles
    psm_ref::PSMFileReference
    protein_ref::ProteinGroupFileReference
    ms_file_idx::Int64
    
    # Constructor validates pairing
    function PairedSearchFiles(psm_path::String, protein_path::String, idx::Int64)
        psm_ref = PSMFileReference(psm_path)
        protein_ref = ProteinGroupFileReference(protein_path)
        
        # Validate both files exist or both don't
        if psm_ref.file_exists != protein_ref.file_exists
            error("Inconsistent file existence: PSM exists=$(psm_ref.file_exists), " *
                  "Protein exists=$(protein_ref.file_exists)")
        end
        
        new(psm_ref, protein_ref, idx)
    end
    
    # Alternative constructor from existing references
    function PairedSearchFiles(psm_ref::PSMFileReference, 
                              protein_ref::ProteinGroupFileReference, 
                              idx::Int64)
        # Validate consistency
        if psm_ref.file_exists != protein_ref.file_exists
            error("Inconsistent file existence: PSM exists=$(psm_ref.file_exists), " *
                  "Protein exists=$(protein_ref.file_exists)")
        end
        
        new(psm_ref, protein_ref, idx)
    end
end

#==========================================================
File Operations
==========================================================#

"""
    validate_exists(ref::FileReference)

Validate that a referenced file exists and update the file_exists flag.
"""
function validate_exists(ref::FileReference)
    ref.file_exists = isfile(file_path(ref))
    exists(ref) || error("File does not exist: $(file_path(ref))")
    return true
end

"""
    create_psm_reference(file_path::String) -> PSMFileReference

Create a PSM file reference, extracting metadata from the file.
"""
create_psm_reference(file_path::String) = PSMFileReference(file_path)

"""
    create_protein_reference(file_path::String) -> ProteinGroupFileReference

Create a protein group file reference, extracting metadata from the file.
"""
create_protein_reference(file_path::String) = ProteinGroupFileReference(file_path)

"""
    create_protein_quant_reference(file_path::String) -> ProteinQuantFileReference

Create a protein quantification file reference, extracting metadata from the file.
"""
create_protein_quant_reference(file_path::String) = ProteinQuantFileReference(file_path)

"""
    create_reference(file_path::String, ::Type{T}) where T <: FileReference

Create a file reference of the specified type.
"""
function create_reference(file_path::String, ::Type{PSMFileReference})
    return PSMFileReference(file_path)
end

function create_reference(file_path::String, ::Type{ProteinGroupFileReference})
    return ProteinGroupFileReference(file_path)
end

function create_reference(file_path::String, ::Type{ProteinQuantFileReference})
    return ProteinQuantFileReference(file_path)
end

"""
    describe_reference(ref::FileReference)

Print a human-readable description of a file reference.
"""
function describe_reference(ref::FileReference)
    println("File Reference ($(typeof(ref))):")
    println("  Path: $(file_path(ref))")
    println("  Exists: $(exists(ref))")
    println("  Sorted by: $(isempty(sorted_by(ref)) ? "none" : join(sorted_by(ref), ", "))")
    println("  Row count: $(row_count(ref))")
    println("  Columns: $(join(sort(collect(schema(ref).columns)), ", "))")
end

"""
    validate_schema(ref::FileReference, required_cols::Set{Symbol})

Validate that a file reference has all required columns.
"""
function validate_schema(ref::FileReference, required_cols::Set{Symbol})
    validate_required_columns(schema(ref), required_cols)
end
