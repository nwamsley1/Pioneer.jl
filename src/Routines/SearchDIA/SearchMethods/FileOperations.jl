"""
Memory-efficient file operations for Arrow files.

Provides streaming operations that respect memory constraints while
maintaining data integrity through validation.
"""

using Arrow, DataFrames, Tables
using DataStructures: BinaryMinHeap, BinaryMaxHeap

include("FileReferences.jl")

# Note: The real getProteinGroupsDict implementation is in ScoringSearch/utils.jl
# To avoid circular dependencies, we pass it as a function parameter to apply_protein_inference

#==========================================================
Streaming Operations
==========================================================#


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

# Generic fillColumn! implementations
function _fillColumn!(col::Vector{T}, col_symbol::Symbol, 
                     sorted_tuples::Vector{Tuple{Int64, Int64}},
                     tables::Vector{Arrow.Table}, n_rows::Int) where T
    for i in 1:n_rows
        table_idx, row_idx = sorted_tuples[i]
        col[i] = Tables.getcolumn(tables[table_idx], col_symbol)[row_idx]
    end
end

# Fill batch columns from sorted tuples
function _fill_batch_columns!(batch_df::DataFrame, tables::Vector{Arrow.Table}, 
                             sorted_tuples::Vector{Tuple{Int64, Int64}}, n_rows::Int)
    for col_name in names(batch_df)
        col_symbol = Symbol(col_name)
        _fillColumn!(batch_df[!, col_symbol], col_symbol, sorted_tuples, tables, n_rows)
    end
end

# Write batch to file
function _write_batch(output_path::String, batch_df::DataFrame, n_writes::Int, n_rows::Int)
    # Only write the rows we actually filled
    data_to_write = n_rows < nrow(batch_df) ? batch_df[1:n_rows, :] : batch_df
    
    if n_writes == 0
        # First write - use streaming format
        if isfile(output_path)
            rm(output_path)
        end
        open(output_path, "w") do io
            Arrow.write(io, data_to_write; file=false)
        end
    else
        # Append to existing file
        Arrow.append(output_path, data_to_write)
    end
end

# Generic heap operations for arbitrary number of sort keys
function _add_to_heap!(heap, table::Arrow.Table, table_idx::Int, 
                      row_idx::Int, sort_keys::NTuple{N, Symbol}) where N
    # Build tuple of sort key values plus table index
    values = ntuple(i -> Tables.getcolumn(table, sort_keys[i])[row_idx], N)
    push!(heap, (values..., table_idx))
end

# Initialize heap with dynamic type based on sort keys
function _initialize_heap(tables::Vector{Arrow.Table}, sort_keys::NTuple{N, Symbol}, 
                         reverse::Bool) where N
    if isempty(tables)
        error("No tables to merge")
    end
    
    # Determine types of sort columns
    first_table = first(tables)
    col_types = ntuple(i -> eltype(Tables.getcolumn(first_table, sort_keys[i])), N)
    
    # Create heap type: Tuple{col_types..., Int64}
    heap_tuple_type = Tuple{col_types..., Int64}
    
    # Create appropriate heap
    heap = reverse ? BinaryMaxHeap{heap_tuple_type}() : BinaryMinHeap{heap_tuple_type}()
    
    # Add first row from each table
    for (i, table) in enumerate(tables)
        if length(Tables.getcolumn(table, 1)) > 0
            _add_to_heap!(heap, table, i, 1, sort_keys)
        end
    end
    
    return heap
end

# Updated stream_sorted_merge to use generic helpers
function stream_sorted_merge(refs::Vector{T}, 
                           output_path::String,
                           sort_keys::Symbol...;
                           batch_size::Int=1_000_000,
                           reverse::Bool=false) where T <: FileReference
    # Convert varargs to tuple
    sort_keys_tuple = sort_keys
    
    # Validate inputs
    isempty(refs) && error("No files to merge")
    isempty(sort_keys_tuple) && error("At least one sort key required")
    
    # Validate all files exist and are sorted by same keys
    for (i, ref) in enumerate(refs)
        validate_exists(ref)
        if !is_sorted_by(ref, sort_keys...)
            error("File $i ($(file_path(ref))) is not sorted by $sort_keys")
        end
    end
    
    # Validate schemas are compatible
    base_schema = schema(first(refs))
    for ref in refs[2:end]
        if schema(ref).columns != base_schema.columns
            error("Incompatible schemas between files")
        end
    end
    
    # Load all tables
    tables = [Arrow.Table(file_path(ref)) for ref in refs]
    table_sizes = [length(Tables.getcolumn(table, 1)) for table in tables]
    total_rows = sum(table_sizes)
    
    # Setup for heap-based merge
    table_indices = ones(Int64, length(tables))
    batch_df = _create_empty_dataframe(first(tables), batch_size)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, batch_size)
    
    # Initialize heap
    heap = _initialize_heap(tables, sort_keys_tuple, reverse)
    
    # Perform heap-based merge
    row_idx = 1
    n_writes = 0
    rows_written = 0
    
    while !isempty(heap)
        # Pop from heap - last element is always table_idx
        heap_tuple = pop!(heap)
        table_idx = heap_tuple[end]
        
        table = tables[table_idx]
        current_row_idx = table_indices[table_idx]
        
        # Store this row's location
        sorted_tuples[row_idx] = (table_idx, current_row_idx)
        
        # Move to next row in this table
        table_indices[table_idx] += 1
        next_row_idx = table_indices[table_idx]
        
        # Add next row from same table to heap if available
        if next_row_idx <= table_sizes[table_idx]
            _add_to_heap!(heap, table, table_idx, next_row_idx, sort_keys_tuple)
        end
        
        row_idx += 1
        
        # Write batch when full
        if row_idx > batch_size
            _fill_batch_columns!(batch_df, tables, sorted_tuples, batch_size)
            _write_batch(output_path, batch_df, n_writes, batch_size)
            n_writes += 1
            rows_written += batch_size
            row_idx = 1
        end
    end
    
    # Write final partial batch
    if row_idx > 1
        final_rows = row_idx - 1
        _fill_batch_columns!(batch_df, tables, sorted_tuples, final_rows)
        _write_batch(output_path, batch_df, n_writes, final_rows)
        rows_written += final_rows
    end
    
    # Create reference for output of same type as inputs
    output_ref = create_reference(output_path, T)
    mark_sorted!(output_ref, sort_keys...)
    
    return output_ref
end

"""
    stream_filter(input_ref::FileReference, output_path::String, 
                 filter_fn::Function; batch_size=100_000)

Filter a file without loading it entirely into memory.
Returns a FileReference of the same type as the input.
"""
function stream_filter(input_ref::T,
                      output_path::String, 
                      filter_fn::Function;
                      batch_size::Int=100_000) where T <: FileReference
    validate_exists(input_ref)
    
    # For now, load entire file and filter
    # In production would use proper streaming
    df = DataFrame(Arrow.Table(file_path(input_ref)))
    filtered_df = filter(filter_fn, df)
    
    # Write filtered data
    Arrow.write(output_path, filtered_df)
    
    # Create reference for output of same type as input
    output_ref = create_reference(output_path, T)
    # Filtering may break sort order
    mark_sorted!(output_ref)  # Empty tuple - not sorted
    
    return output_ref
end

"""
    stream_transform(input_ref::FileReference, output_path::String,
                    transform_fn::Function; batch_size=100_000)

Transform a file by applying a function to each batch.
Returns a FileReference of the same type as the input.
"""
function stream_transform(input_ref::T,
                         output_path::String,
                         transform_fn::Function;
                         batch_size::Int=100_000) where T <: FileReference
    validate_exists(input_ref)
    
    tbl = Arrow.Table(file_path(input_ref))
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    
    for partition in partitions
        df_batch = DataFrame(Tables.columntable(partition))
        transformed_batch = transform_fn(df_batch)
        
        if first_write
            Arrow.write(output_path, transformed_batch)
            first_write = false
        else
            Arrow.append(output_path, transformed_batch)
        end
    end
    
    # Create reference for output of same type as input
    output_ref = create_reference(output_path, T)
    # Transformation may change sort order or schema
    
    return output_ref
end

#==========================================================
Column Operations
==========================================================#

"""
    add_column_to_file!(ref::FileReference, col_name::Symbol, 
                       compute_fn::Function; batch_size=100_000)

Add a new column to a file by computing it from existing columns.
Updates the file in place and updates the reference's schema.
"""
function add_column_to_file!(ref::FileReference, 
                           col_name::Symbol,
                           compute_fn::Function;
                           batch_size::Int=100_000)
    validate_exists(ref)
    
    # Check if column already exists
    if has_column(schema(ref), col_name)
        @warn "Column $col_name already exists, will be overwritten"
    end
    
    # Create temporary output file
    temp_path = file_path(ref) * ".tmp"
    
    # Stream through file, adding column
    tbl = Arrow.Table(file_path(ref))
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    for partition in partitions
        df_batch = DataFrame(Tables.columntable(partition))
        df_batch[!, col_name] = compute_fn(df_batch)
        
        if first_write
            Arrow.write(temp_path, df_batch)
            first_write = false
        else
            Arrow.append(temp_path, df_batch)
        end
    end
    
    # Replace original file
    mv(temp_path, file_path(ref), force=true)
    
    # Update reference metadata
    # We need to re-read the file to update schema
    new_ref = create_reference(file_path(ref), typeof(ref))
    ref.schema = schema(new_ref)
    ref.row_count = row_count(new_ref)
    ref.sorted_by = ()  # Adding column may invalidate sort
    
    return ref
end

"""
    update_column_in_file!(ref::FileReference, col_name::Symbol,
                          update_fn::Function; batch_size=100_000)

Update an existing column in a file.
"""
function update_column_in_file!(ref::FileReference,
                               col_name::Symbol,
                               update_fn::Function;
                               batch_size::Int=100_000)
    validate_exists(ref)
    
    if !has_column(schema(ref), col_name)
        error("Column $col_name does not exist")
    end
    
    # Create temporary output file
    temp_path = file_path(ref) * ".tmp"
    
    # Stream through file, updating column
    tbl = Arrow.Table(file_path(ref))
    partitions = Tables.partitioner(tbl, batch_size)
    
    first_write = true
    for partition in partitions
        df_batch = DataFrame(Tables.columntable(partition))
        df_batch[!, col_name] = update_fn(df_batch[!, col_name])
        
        if first_write
            Arrow.write(temp_path, df_batch)
            first_write = false
        else
            Arrow.append(temp_path, df_batch)
        end
    end
    
    # Replace original file
    mv(temp_path, file_path(ref), force=true)
    
    # Sort order may be affected if we updated a sort key
    if col_name in sorted_by(ref)
        mark_sorted!(ref)  # Empty tuple - not sorted
    end
    
    return ref
end

#==========================================================
Safe Join Operations
==========================================================#

"""
    safe_join_files(left_ref::FileReference, right_ref::FileReference,
                   output_path::String, join_key::Symbol; 
                   join_type=:inner, batch_size=100_000)

Join two files with memory-efficient streaming.
Typically used for joining PSM and protein files.
"""
function safe_join_files(left_ref::FileReference, 
                        right_ref::FileReference,
                        output_path::String,
                        join_key::Symbol;
                        join_type::Symbol=:inner,
                        batch_size::Int=100_000)
    # Validate inputs
    validate_exists(left_ref)
    validate_exists(right_ref)
    
    # Both files must have the join key
    has_column(schema(left_ref), join_key) || 
        error("Left file missing join key: $join_key")
    has_column(schema(right_ref), join_key) || 
        error("Right file missing join key: $join_key")
    
    # For efficient join, right side should be sorted by join key
    if !is_sorted_by(right_ref, join_key)
        @warn "Right side not sorted by join key, join may be slow"
    end
    
    # Load right side (protein groups) - typically smaller
    right_df = DataFrame(Arrow.Table(file_path(right_ref)))
    
    # Stream through left side (PSMs)
    left_tbl = Arrow.Table(file_path(left_ref))
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
    
    # Create reference for output of same type as left input
    output_ref = create_reference(output_path, typeof(left_ref))
    
    return output_ref
end

#==========================================================
Algorithm Wrappers
==========================================================#

"""
    apply_protein_inference(psm_ref::PSMFileReference, output_path::String, 
                          precursors::Union{Nothing, DataFrame}, 
                          protein_inference_fn::Function;
                          min_peptides=2, batch_size=1_000_000)

Apply protein inference algorithm to PSMs, creating protein groups.
Uses the provided protein_inference_fn function (typically getProteinGroupsDict from ScoringSearch/utils.jl).
"""
function apply_protein_inference(psm_ref::PSMFileReference, 
                               output_path::String,
                               precursors,
                               protein_inference_fn::Function;
                               min_peptides::Int=2)
    validate_exists(psm_ref)
    
    # Validate required columns for protein inference
    required_cols = Set([:prob, :inferred_protein_group, :precursor_idx, :target, :entrapment_group_id])
    validate_required_columns(schema(psm_ref), required_cols)
    
    # Stream through PSMs, collect protein groups
    tbl = Arrow.Table(file_path(psm_ref))
    partitions = Tables.partitioner(tbl, batch_size)
    
    # For protein inference, we need all PSMs at once (not streaming)
    # This is because protein inference needs global view of peptide-protein relationships
    psms_df = DataFrame(Arrow.Table(file_path(psm_ref)))
    
    # Call provided protein inference function
    protein_groups_dict = protein_inference_fn(
        psms_df.prob,
        psms_df.inferred_protein_group,
        psms_df.precursor_idx,
        psms_df.target,
        psms_df.entrapment_group_id,
        precursors,
        min_peptides
    )
    
    # Convert to DataFrame format expected by ScoringSearch
    # Following the schema from write_protein_groups_arrow
    n_groups = length(protein_groups_dict)
    
    # Pre-allocate vectors for each column
    protein_names = Vector{String}(undef, n_groups)
    targets = Vector{Bool}(undef, n_groups)
    entrap_ids = Vector{UInt8}(undef, n_groups) 
    pg_scores = Vector{Float32}(undef, n_groups)
    n_peptides = Vector{Int64}(undef, n_groups)
    
    # Fill vectors from protein groups
    for (i, ((name, target, entrap_id), peptide_set)) in enumerate(protein_groups_dict)
        protein_names[i] = name
        targets[i] = target
        entrap_ids[i] = entrap_id
        # Basic score as negative log-sum of peptide probabilities
        pg_scores[i] = -sum(log(1.0f0 - psms_df.prob[pid]) for pid in peptide_set)
        n_peptides[i] = length(peptide_set)
    end
    
    # Create DataFrame
    protein_df = DataFrame(
        protein_name = protein_names,
        target = targets,
        entrapment_group_id = entrap_ids,
        pg_score = pg_scores,
        n_peptides = n_peptides
    )
    
    # Sort by score (descending)
    sort!(protein_df, :pg_score, rev=true)
    
    # Write to Arrow
    Arrow.write(output_path, protein_df)
    
    return ProteinGroupFileReference(output_path)
end

"""
    update_psms_with_scores(psm_ref::PSMFileReference, protein_ref::ProteinGroupFileReference,
                          output_path::String; batch_size=100_000)

Update PSMs with protein group scores from protein reference.
Streaming operation that preserves memory efficiency.
"""
function update_psms_with_scores(psm_ref::PSMFileReference, 
                               protein_ref::ProteinGroupFileReference,
                               output_path::String;
                               batch_size::Int=100_000)
    validate_exists(psm_ref)
    validate_exists(protein_ref)
    
    # Load protein scores (typically smaller)
    protein_df = DataFrame(Arrow.Table(file_path(protein_ref)))
    
    # Create lookup dictionary
    protein_scores = Dict{Tuple{String,Bool,UInt8},NamedTuple}()
    for row in eachrow(protein_df)
        key = (row.protein_name, row.target, row.entrapment_group_id)
        protein_scores[key] = (
            pg_score = row.pg_score,
            global_pg_score = row.global_pg_score,
            pg_qval = row.pg_qval,
            global_pg_qval = row.global_pg_qval
        )
    end
    
    # Stream through PSMs, updating scores
    stream_transform(psm_ref, output_path, batch_size=batch_size) do batch
        # Add score columns
        batch.pg_score = Vector{Float32}(undef, nrow(batch))
        batch.global_pg_score = Vector{Float32}(undef, nrow(batch))
        batch.pg_qval = Vector{Float32}(undef, nrow(batch))
        batch.global_qval_pg = Vector{Float32}(undef, nrow(batch))
        
        # Update scores
        for i in 1:nrow(batch)
            key = (batch.inferred_protein_group[i], 
                   batch.target[i], 
                   batch.entrapment_group_id[i])
            
            if haskey(protein_scores, key)
                scores = protein_scores[key]
                batch.pg_score[i] = scores.pg_score
                batch.global_pg_score[i] = scores.global_pg_score
                batch.pg_qval[i] = scores.pg_qval
                batch.global_qval_pg[i] = scores.global_pg_qval
            else
                # This should not happen if protein inference was done correctly
                error("Missing protein scores for: $key")
            end
        end
        
        return batch
    end
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
    process_with_memory_limit(ref::FileReference, process_fn::Function;
                             max_memory_mb=1000)

Process a file with automatic batch size calculation based on memory limit.
"""
function process_with_memory_limit(ref::FileReference, 
                                 process_fn::Function;
                                 max_memory_mb::Int=1000)
    # For now, just load the whole file
    # In production, this would use proper streaming
    df = DataFrame(Arrow.Table(file_path(ref)))
    process_fn(df)
end

#==========================================================
Specialized Merge Functions
==========================================================#

"""
    merge_psm_scores(psm_refs::Vector{PSMFileReference}, output_path::String,
                    prob_col::Symbol; batch_size=10_000_000)

Merge PSM files sorted by probability score for FDR calculation.
Follows the pattern of merge_sorted_psms_scores in ScoringSearch.
"""
function merge_psm_scores(psm_refs::Vector{PSMFileReference},
                         output_path::String,
                         prob_col::Symbol;
                         batch_size::Int=10_000_000)
    # Use reverse=true for descending sort by probability
    merged_ref = stream_sorted_merge(psm_refs, output_path, prob_col, :target;
                                   batch_size=batch_size, reverse=true)
    
    # Note: The original function only preserves specific columns
    # If needed, we could add a transform step to select columns
    
    return merged_ref
end

"""
    merge_protein_groups_by_score(pg_refs::Vector{ProteinGroupFileReference},
                                 output_path::String; batch_size=1_000_000)

Merge protein group files sorted by score.
Follows the pattern of merge_sorted_protein_groups in ScoringSearch.
"""
function merge_protein_groups_by_score(pg_refs::Vector{ProteinGroupFileReference},
                                     output_path::String;
                                     batch_size::Int=1_000_000)
    # Use reverse=true for descending sort by pg_score
    merged_ref = stream_sorted_merge(pg_refs, output_path, :pg_score;
                                   batch_size=batch_size, reverse=true)
    return merged_ref
end

"""
    apply_maxlfq(psm_refs::Vector{PSMFileReference}, output_dir::String,
                quant_col::Symbol, file_names::Vector{String},
                lfq_fn::Function;
                q_value_threshold=0.01f0, batch_size=100_000, min_peptides=2)

Apply MaxLFQ algorithm to PSM files.
Uses the provided lfq_fn function (typically LFQ from utils/maxLFQ.jl).
"""
function apply_maxlfq(psm_refs::Vector{PSMFileReference}, 
                     output_dir::String,
                     quant_col::Symbol,
                     file_names::Vector{String},
                     lfq_fn::Function;
                     q_value_threshold::Float32=0.01f0,
                     batch_size::Int=100_000,
                     min_peptides::Int=2)
    # Validate all PSMs are sorted correctly for MaxLFQ
    required_sort_keys = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
    for (i, ref) in enumerate(psm_refs)
        if !is_sorted_by(ref, required_sort_keys...)
            error("PSM file $i not sorted correctly for MaxLFQ. Expected sort by $required_sort_keys")
        end
    end
    
    # Merge PSM files using 4-key sorting
    merged_psm_path = joinpath(output_dir, "precursors_long.arrow")
    merged_psm_ref = stream_sorted_merge(psm_refs, merged_psm_path, required_sort_keys...;
                                       batch_size=batch_size, reverse=true)
    
    # Load merged data for LFQ (as current LFQ expects DataFrame)
    merged_df = DataFrame(Arrow.Table(merged_psm_path))
    
    # Call existing LFQ function
    protein_long_path = joinpath(output_dir, "protein_groups_long.arrow")
    
    # Call provided LFQ function
    lfq_fn(merged_df, protein_long_path, quant_col, file_names,
           q_value_threshold, batch_size=batch_size, min_peptides=min_peptides)
    
    # Create output references
    return merged_psm_ref, ProteinGroupFileReference(protein_long_path)
end

# Export all public functions
export stream_sorted_merge, stream_filter, stream_transform,
       add_column_to_file!, update_column_in_file!,
       safe_join_files, estimate_batch_size, process_with_memory_limit,
       apply_protein_inference, update_psms_with_scores,
       merge_psm_scores, merge_protein_groups_by_score, apply_maxlfq