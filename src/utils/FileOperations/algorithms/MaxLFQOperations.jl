"""
MaxLFQ algorithm operations and validation.

Provides specialized operations for MaxLFQ quantification
including validation, preparation, and algorithm application.
"""

using Arrow, DataFrames

#==========================================================
MaxLFQ Operations
==========================================================#

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

"""
    prepare_maxlfq_input(psm_refs::Vector{PSMFileReference}, output_dir::String;
                        q_value_threshold::Float32=0.01f0)

Prepare PSM files for MaxLFQ by applying necessary preprocessing and validation.
Returns preprocessed file references.
"""
function prepare_maxlfq_input(psm_refs::Vector{PSMFileReference}, output_dir::String;
                             q_value_threshold::Float32=0.01f0)
    !isdir(output_dir) && mkpath(output_dir)
    
    # Create preprocessing pipeline using existing operations
    preprocessing_pipeline = TransformPipeline() |>
        filter_by_multiple_thresholds([
            (:pg_qval, q_value_threshold),
            (:global_qval_pg, q_value_threshold)
        ]) |>
        filter_rows(row -> row.use_for_protein_quant; desc="filter_for_protein_quant")
    
    # Apply preprocessing to all files
    @info "Applying MaxLFQ preprocessing to $(length(psm_refs)) files..."
    processed_refs = apply_pipeline_batch(psm_refs, preprocessing_pipeline, output_dir)
    
    # Validate and sort each processed file
    required_sort_keys = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
    
    @info "Ensuring correct sort order for MaxLFQ..."
    for ref in processed_refs
        if exists(ref)
            # Validate required columns exist
            validate_maxlfq_input(ref)
            
            # Ensure proper sorting
            if !is_sorted_by(ref, required_sort_keys...)
                @info "Sorting $(basename(file_path(ref))) for MaxLFQ..."
                sort_file_by_keys!(ref, required_sort_keys...; reverse=true)
            end
        end
    end
    
    @info "MaxLFQ preprocessing complete. $(length(processed_refs)) files ready."
    return processed_refs
end

# Export MaxLFQ functions
export apply_maxlfq, prepare_maxlfq_input