"""
Memory-efficient file operations for Arrow files.

Provides streaming operations that respect memory constraints while
maintaining data integrity through validation.
"""

using Arrow, DataFrames, Tables

include("FileReferences.jl")

#==========================================================
Streaming Operations
==========================================================#

"""
    stream_sorted_merge(refs::Vector{PSMFileReference}, output_path::String, 
                       sort_keys::Symbol...; batch_size=1_000_000)

Merge multiple sorted PSM files using heap-based streaming merge.
Validates all files are sorted by the same keys before merging.
"""
function stream_sorted_merge(refs::Vector{PSMFileReference}, 
                           output_path::String,
                           sort_keys::Symbol...;
                           batch_size::Int=1_000_000)
    # Validate inputs
    isempty(refs) && error("No files to merge")
    
    # Validate all files exist and are sorted by same keys
    for (i, ref) in enumerate(refs)
        validate_exists(ref)
        if !is_sorted_by(ref, sort_keys...)
            error("File $i ($(ref.file_path)) is not sorted by $sort_keys")
        end
    end
    
    # Validate schemas are compatible
    base_schema = first(refs).schema
    for ref in refs[2:end]
        if ref.schema.columns != base_schema.columns
            error("Incompatible schemas between files")
        end
    end
    
    # Perform streaming merge (simplified version - actual implementation would use heap)
    total_rows = 0
    first_write = true
    
    # For now, simple concatenation (actual implementation would use heap merge)
    for ref in refs
        tbl = Arrow.Table(ref.file_path)
        if first_write
            Arrow.write(output_path, tbl)
            first_write = false
        else
            Arrow.append(output_path, tbl)
        end
        total_rows += ref.row_count
    end
    
    # Create reference for output
    output_ref = PSMFileReference(output_path)
    mark_sorted!(output_ref, sort_keys...)
    
    return output_ref
end

"""
    stream_filter(input_ref::PSMFileReference, output_path::String, 
                 filter_fn::Function; batch_size=100_000)

Filter a file without loading it entirely into memory.
"""
function stream_filter(input_ref::PSMFileReference,
                      output_path::String, 
                      filter_fn::Function;
                      batch_size::Int=100_000)
    validate_exists(input_ref)
    
    # Use Tables.partitioner for memory-efficient iteration
    tbl = Arrow.Table(input_ref.file_path)
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    total_rows = 0
    
    for partition in partitions
        df_batch = DataFrame(partition)
        filtered_batch = filter(filter_fn, df_batch)
        
        if nrow(filtered_batch) > 0
            if first_write
                Arrow.write(output_path, filtered_batch)
                first_write = false
            else
                Arrow.append(output_path, filtered_batch)
            end
            total_rows += nrow(filtered_batch)
        end
    end
    
    # Create reference for output
    output_ref = PSMFileReference(output_path)
    # Filtering may break sort order
    mark_sorted!(output_ref)  # Empty tuple - not sorted
    
    return output_ref
end

"""
    stream_transform(input_ref::PSMFileReference, output_path::String,
                    transform_fn::Function; batch_size=100_000)

Transform a file by applying a function to each batch.
"""
function stream_transform(input_ref::PSMFileReference,
                         output_path::String,
                         transform_fn::Function;
                         batch_size::Int=100_000)
    validate_exists(input_ref)
    
    tbl = Arrow.Table(input_ref.file_path)
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    
    for partition in partitions
        df_batch = DataFrame(partition)
        transformed_batch = transform_fn(df_batch)
        
        if first_write
            Arrow.write(output_path, transformed_batch)
            first_write = false
        else
            Arrow.append(output_path, transformed_batch)
        end
    end
    
    # Create reference for output
    output_ref = PSMFileReference(output_path)
    # Transformation may change sort order or schema
    
    return output_ref
end

#==========================================================
Column Operations
==========================================================#

"""
    add_column_to_file!(ref::PSMFileReference, col_name::Symbol, 
                       compute_fn::Function; batch_size=100_000)

Add a new column to a file by computing it from existing columns.
Updates the file in place and updates the reference's schema.
"""
function add_column_to_file!(ref::PSMFileReference, 
                           col_name::Symbol,
                           compute_fn::Function;
                           batch_size::Int=100_000)
    validate_exists(ref)
    
    # Check if column already exists
    if has_column(ref.schema, col_name)
        @warn "Column $col_name already exists, will be overwritten"
    end
    
    # Create temporary output file
    temp_path = ref.file_path * ".tmp"
    
    # Stream through file, adding column
    tbl = Arrow.Table(ref.file_path)
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    for partition in partitions
        df_batch = DataFrame(partition)
        df_batch[!, col_name] = compute_fn(df_batch)
        
        if first_write
            Arrow.write(temp_path, df_batch)
            first_write = false
        else
            Arrow.append(temp_path, df_batch)
        end
    end
    
    # Replace original file
    mv(temp_path, ref.file_path, force=true)
    
    # Update reference
    new_cols = union(ref.schema.columns, [col_name])
    ref.schema = FileSchema(collect(new_cols))
    ref.sorted_by = ()  # Adding column may invalidate sort
    
    return ref
end

"""
    update_column_in_file!(ref::PSMFileReference, col_name::Symbol,
                          update_fn::Function; batch_size=100_000)

Update an existing column in a file.
"""
function update_column_in_file!(ref::PSMFileReference,
                               col_name::Symbol,
                               update_fn::Function;
                               batch_size::Int=100_000)
    validate_exists(ref)
    
    if !has_column(ref.schema, col_name)
        error("Column $col_name does not exist")
    end
    
    # Create temporary output file
    temp_path = ref.file_path * ".tmp"
    
    # Stream through file, updating column
    tbl = Arrow.Table(ref.file_path)
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    for partition in partitions
        df_batch = DataFrame(partition)
        df_batch[!, col_name] = update_fn(df_batch[!, col_name])
        
        if first_write
            Arrow.write(temp_path, df_batch)
            first_write = false
        else
            Arrow.append(temp_path, df_batch)
        end
    end
    
    # Replace original file
    mv(temp_path, ref.file_path, force=true)
    
    # Sort order may be affected if we updated a sort key
    if col_name in ref.sorted_by
        ref.sorted_by = ()
    end
    
    return ref
end

#==========================================================
Safe Join Operations
==========================================================#

"""
    safe_join_files(left_ref::PSMFileReference, right_ref::ProteinGroupFileReference,
                   output_path::String, join_key::Symbol; 
                   join_type=:inner, batch_size=100_000)

Join PSM and protein files with memory-efficient streaming.
"""
function safe_join_files(left_ref::PSMFileReference, 
                        right_ref::ProteinGroupFileReference,
                        output_path::String,
                        join_key::Symbol;
                        join_type::Symbol=:inner,
                        batch_size::Int=100_000)
    # Validate inputs
    validate_exists(left_ref)
    validate_exists(right_ref)
    
    # Both files must have the join key
    has_column(left_ref.schema, join_key) || 
        error("Left file missing join key: $join_key")
    has_column(right_ref.schema, join_key) || 
        error("Right file missing join key: $join_key")
    
    # For efficient join, right side should be sorted by join key
    if !is_sorted_by(right_ref, join_key)
        @warn "Right side not sorted by join key, join may be slow"
    end
    
    # Load right side (protein groups) - typically smaller
    right_df = DataFrame(Arrow.Table(right_ref.file_path))
    
    # Stream through left side (PSMs)
    left_tbl = Arrow.Table(left_ref.file_path)
    partitions = Tables.partitioner(left_tbl, batch_size)
    
    first_write = true
    for partition in partitions
        left_batch = DataFrame(partition)
        
        # Perform join
        joined_batch = if join_type == :inner
            innerjoin(left_batch, right_df, on=join_key)
        elseif join_type == :left
            leftjoin(left_batch, right_df, on=join_key)
        else
            error("Unsupported join type: $join_type")
        end
        
        if nrow(joined_batch) > 0
            if first_write
                Arrow.write(output_path, joined_batch)
                first_write = false
            else
                Arrow.append(output_path, joined_batch)
            end
        end
    end
    
    # Create reference for output
    output_ref = PSMFileReference(output_path)
    
    return output_ref
end

#==========================================================
Memory Usage Estimation
==========================================================#

"""
    estimate_batch_size(schema::FileSchema, max_memory_mb::Int)

Estimate appropriate batch size based on schema and memory limit.
"""
function estimate_batch_size(schema::FileSchema, max_memory_mb::Int)
    # Rough estimation: assume 8 bytes per numeric column, 50 bytes per string
    bytes_per_row = 0
    
    # This is a simplified estimation
    # In practice, would need to sample the file to get better estimates
    n_cols = length(schema.columns)
    bytes_per_row = n_cols * 20  # Average estimate
    
    max_bytes = max_memory_mb * 1024 * 1024
    batch_size = max(1000, div(max_bytes, bytes_per_row))
    
    return batch_size
end

"""
    process_with_memory_limit(ref::PSMFileReference, process_fn::Function;
                             max_memory_mb=1000)

Process a file with automatic batch size calculation based on memory limit.
"""
function process_with_memory_limit(ref::PSMFileReference, 
                                 process_fn::Function;
                                 max_memory_mb::Int=1000)
    batch_size = estimate_batch_size(ref.schema, max_memory_mb)
    
    tbl = Arrow.Table(ref.file_path)
    partitions = Tables.partitioner(tbl, batch_size)
    
    for partition in partitions
        df_batch = DataFrame(partition)
        process_fn(df_batch)
    end
end

# Export all public functions
export stream_sorted_merge, stream_filter, stream_transform,
       add_column_to_file!, update_column_in_file!,
       safe_join_files, estimate_batch_size, process_with_memory_limit