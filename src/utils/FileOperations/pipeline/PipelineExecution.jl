# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
Pipeline execution and batch processing.

Provides utilities for executing transformation pipelines
on single files or batches of files with parallel support.
"""

using DataFrames

#==========================================================
Pipeline Execution
==========================================================#

"""
    apply_pipeline!(ref::FileReference, pipeline::TransformPipeline) -> FileReference

Apply all operations in the pipeline to the file in a single pass.
Automatically runs post-actions after transformation completes.
"""
function apply_pipeline!(ref::FileReference, pipeline::TransformPipeline)
    if isempty(pipeline.operations)
        return ref
    end

    # Check if file exists and has data before processing
    if !exists(ref)
        @debug_l1 "Skipping non-existent file: $(file_path(ref))"
        return ref
    end

    # Apply all operations in single pass with empty data detection
    transform_and_write!(ref) do df
        # Check for empty input data
        if nrow(df) == 0
            @debug_l1 "Skipping empty file: $(file_path(ref))"
            return df  # Return empty DataFrame unchanged
        end

        for (desc, op) in pipeline.operations
            try
                df = op(df)
                # Check if operation resulted in empty data
                if nrow(df) == 0
                    @debug_l1 "File became empty after operation '$desc': $(file_path(ref))"
                    break  # No point in continuing pipeline
                end
            catch e
                @user_warn "Pipeline operation '$desc' failed for file $(file_path(ref)): $e"
                rethrow(e)  # Let caller handle the error appropriately
            end
        end
        return df
    end

    # Run post-actions (like mark_sorted!) only if file still has data
    if exists(ref) && row_count(ref) > 0
        for action in pipeline.post_actions
            action(ref)
        end
    end

    return ref
end

"""
    apply_pipeline!(refs::Vector{<:FileReference}, pipeline::TransformPipeline)

Apply pipeline to multiple files, optionally in parallel.
"""
function apply_pipeline!(refs::Vector{<:FileReference}, pipeline::TransformPipeline;
                        parallel::Bool=true)
    if parallel && length(refs) > 1
        Threads.@threads for ref in refs
            if exists(ref)
                try
                    apply_pipeline!(ref, pipeline)
                catch e
                    @user_warn "Pipeline failed for file $(file_path(ref)): $e"
                    # Continue processing other files
                end
            end
        end
    else
        for ref in refs
            if exists(ref)
                try
                    apply_pipeline!(ref, pipeline)
                catch e
                    @user_warn "Pipeline failed for file $(file_path(ref)): $e"
                    # Continue processing other files
                end
            end
        end
    end
    return refs
end

"""
    transform_and_write!(ref::FileReference, output_path::String, pipeline::TransformPipeline)

Apply a pipeline to a file and write the result to a new location.
Does not modify the original file.

Example:
```julia
transform_and_write!(input_ref, "output.arrow", my_pipeline)
```
"""
function transform_and_write!(ref::FileReference, output_path::String, pipeline::TransformPipeline)
    # Validate input
    validate_exists(ref)
    
    # Ensure output directory exists
    output_dir = dirname(output_path)
    !isdir(output_dir) && mkpath(output_dir)
    
    # Create transform function that applies all operations
    transform_fn = function(df)
        for (desc, op) in pipeline.operations
            df = op(df)
        end
        return df
    end
    
    # Use existing transform_and_write! with output path
    transform_and_write!(transform_fn, ref, output_path)
    
    return nothing
end

"""
    apply_pipeline_batch(refs::Vector{<:FileReference}, pipeline::TransformPipeline,
                        output_folder::String; preserve_basename::Bool = true)

Apply a pipeline to multiple files and write results to a new folder.
Returns vector of new FileReferences.

Example:
```julia
passing_refs = apply_pipeline_batch(
    filtered_refs,
    qvalue_pipeline,
    "output/passing_psms"
)
```
"""
function apply_pipeline_batch(refs::Vector{<:FileReference}, 
                            pipeline::TransformPipeline,
                            output_folder::String;
                            preserve_basename::Bool = true)
    # Ensure output folder exists
    !isdir(output_folder) && mkpath(output_folder)
    
    new_refs = similar(refs, 0)  # Empty vector of same type
    
    for ref in refs
        if exists(ref)
            # Determine output path
            output_name = preserve_basename ? basename(file_path(ref)) : 
                         "processed_$(length(new_refs) + 1).arrow"
            output_path = joinpath(output_folder, output_name)
            
            # Apply pipeline and write to new location in one pass
            # This avoids modifying the original file
            transform_and_write!(ref, output_path, pipeline)
            
            # Create new reference
            new_ref = if ref isa PSMFileReference
                PSMFileReference(output_path)
            elseif ref isa ProteinGroupFileReference
                ProteinGroupFileReference(output_path)
            else
                FileReference(output_path)
            end
            
            # Copy metadata if applicable
            if !isempty(sorted_by(ref))
                mark_sorted!(new_ref, sorted_by(ref)...)
            end
            
            push!(new_refs, new_ref)
        end
    end
    
    return new_refs
end


