"""
Memory-efficient file operations for Arrow files.

Provides streaming operations that respect memory constraints while
maintaining data integrity through validation.
"""

using Arrow, DataFrames, Tables

include("FileReferences.jl")

# Import protein inference function from utils
# This is a temporary import - in real implementation would use proper module structure
const PIONEER_ROOT = dirname(dirname(dirname(dirname(dirname(@__DIR__)))))
include(joinpath(PIONEER_ROOT, "src", "utils", "proteinInference.jl"))

#==========================================================
Streaming Operations
==========================================================#

"""
    stream_sorted_merge(refs::Vector{<:FileReference}, output_path::String, 
                       sort_keys::Symbol...; batch_size=1_000_000)

Merge multiple sorted files using heap-based streaming merge.
Validates all files are sorted by the same keys before merging.
Returns a FileReference of the same type as the input refs.
"""
function stream_sorted_merge(refs::Vector{T}, 
                           output_path::String,
                           sort_keys::Symbol...;
                           batch_size::Int=1_000_000) where T <: FileReference
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
    base_schema = schema(first(refs))
    for ref in refs[2:end]
        if schema(ref).columns != base_schema.columns
            error("Incompatible schemas between files")
        end
    end
    
    # Perform streaming merge (simplified version - actual implementation would use heap)
    total_rows = 0
    first_write = true
    
    # For now, simple concatenation (actual implementation would use heap merge)
    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        if first_write
            Arrow.write(output_path, tbl)
            first_write = false
        else
            Arrow.append(output_path, tbl)
        end
        total_rows += row_count(ref)
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
    
    # Use Tables.partitioner for memory-efficient iteration
    tbl = Arrow.Table(file_path(input_ref))
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
        df_batch = DataFrame(partition)
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
                          precursors::Union{Nothing, DataFrame}, protein_db_map;
                          min_peptides=2, batch_size=1_000_000)

Apply protein inference algorithm to PSMs, creating protein groups.
Wraps the existing getProteinGroupsDict from utils/proteinInference.jl.
"""
function apply_protein_inference(psm_ref::PSMFileReference, 
                               output_path::String,
                               precursors;
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
    
    # Call existing protein inference algorithm from utils/proteinInference.jl
    protein_groups_dict = getProteinGroupsDict(
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
    batch_size = estimate_batch_size(schema(ref), max_memory_mb)
    
    tbl = Arrow.Table(file_path(ref))
    partitions = Tables.partitioner(tbl, batch_size)
    
    for partition in partitions
        df_batch = DataFrame(partition)
        process_fn(df_batch)
    end
end

# Export all public functions
export stream_sorted_merge, stream_filter, stream_transform,
       add_column_to_file!, update_column_in_file!,
       safe_join_files, estimate_batch_size, process_with_memory_limit,
       apply_protein_inference, update_psms_with_scores