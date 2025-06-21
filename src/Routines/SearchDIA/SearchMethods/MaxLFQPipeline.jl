"""
    MaxLFQPipeline.jl

Simple pipeline operations for MaxLFQ processing integration.
Provides basic utilities for FileReference and TransformPipeline integration
without over-optimizing the working MaxLFQ algorithms.
"""

using DataFrames

# Include required modules (these are included by importScripts.jl in SearchMethods)
# include("FileReferences.jl")
# include("FileOperations.jl")

#==========================================================
Simple MaxLFQ Pipeline Integration
==========================================================#

"""
    create_maxlfq_preprocessing_pipeline(q_value_threshold::Float32)

Create a standardized preprocessing pipeline for MaxLFQ.
Uses simple composition of existing operations.
Note: This is now handled directly in the LFQ function with FileReference support.
"""
function create_maxlfq_preprocessing_pipeline(q_value_threshold::Float32)
    return TransformPipeline() |>
        filter_by_multiple_thresholds([
            (:pg_qval, q_value_threshold),
            (:global_qval_pg, q_value_threshold)
        ]) |>
        filter_rows(row -> row.use_for_protein_quant; desc="filter_for_protein_quant")
end

"""
    create_maxlfq_results_pipeline()

Create a post-processing pipeline for MaxLFQ results.
Adds computed columns and performs basic validation.
"""
function create_maxlfq_results_pipeline()
    return TransformPipeline() |>
        add_column(:abundance, row -> ismissing(row.log2_abundance) ? missing : exp2(row.log2_abundance)) |>
        filter_rows(row -> !ismissing(row.n_peptides); desc="filter_valid_peptide_counts")
end

#==========================================================
Exports (Simplified)
==========================================================#

export create_maxlfq_preprocessing_pipeline, create_maxlfq_results_pipeline