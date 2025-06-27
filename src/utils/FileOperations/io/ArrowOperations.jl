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
                           reverse::Union{Bool, Vector{Bool}}=false)
    validate_exists(ref)
    
    # Validate that all sort keys exist in schema
    for key in sort_keys
        if !has_column(schema(ref), key)
            error("Sort key $key not found in file schema")
        end
    end
    
    # Load file into memory (not memory-mapped)
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    
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
    sort!(df, collect(sort_keys), rev=rev_vec)
    
    # Write back to same file using writeArrow for Windows compatibility
    writeArrow(file_path(ref), df)
    
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
                           reverse::Union{Bool, Vector{Bool}}=false, parallel::Bool=true)
    # Validate reverse vector length if it's a vector
    if reverse isa Vector{Bool} && length(reverse) != length(keys)
        error("Length of reverse vector ($(length(reverse))) must match number of sort keys ($(length(keys)))")
    end
    
    if parallel && length(refs) > 1
        Threads.@threads for ref in ProgressBar(refs)
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse)
            end
        end
    else
        # Sequential with traditional progress bar
        for ref in ProgressBar(refs)
            if exists(ref)
                sort_file_by_keys!(ref, keys...; reverse=reverse)
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
    
    # Load entire file into memory
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    
    # Apply transformation
    transformed_df = transform_fn(df)
    
    # Write back using writeArrow for Windows compatibility
    return write_arrow_file(ref, transformed_df)
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
    
    # Load entire file into memory
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    
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
    return DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
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

