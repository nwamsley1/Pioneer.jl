# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
Validation utilities for file operations and algorithms.

Provides validation functions for algorithm parameters,
memory estimation, and input validation for various
file operations.
"""

#==========================================================
MaxLFQ Validation Functions
==========================================================#

"""
    validate_maxlfq_input(ref::PSMFileReference)

Validate that a PSM file reference has all required columns for MaxLFQ processing
and is sorted correctly.
"""
function validate_maxlfq_input(ref::PSMFileReference)
    validate_exists(ref)
    
    # Check required columns for MaxLFQ
    required_columns = Set([
        :inferred_protein_group, :target, :entrapment_group_id, :precursor_idx,
        :pg_qval, :global_pg_qval, :use_for_protein_quant, :peak_area
    ])
    
    validate_required_columns(schema(ref), required_columns)
    
    # Validate sort order (correct sort is critical for MaxLFQ batching)
    expected_sort = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
    if !is_sorted_by(ref, expected_sort...)
        error("Input must be sorted by $expected_sort for MaxLFQ processing. " *
              "Current sort: $(sorted_by(ref))")
    end
    
    return true
end

"""
    validate_maxlfq_parameters(params::Dict)

Validate MaxLFQ parameters for safety and correctness.
"""
function validate_maxlfq_parameters(params::Dict)
    # Check for required parameters
    required_params = [:q_value_threshold, :batch_size, :min_peptides]
    for param in required_params
        if !haskey(params, param)
            error("Missing required parameter: $param")
        end
    end
    
    # Validate q-value threshold
    if haskey(params, :q_value_threshold)
        q_val = params[:q_value_threshold]
        if !(0.0 < q_val <= 1.0)
            error("q_value_threshold must be between 0 and 1, got: $q_val")
        end
    end
    
    # Validate batch size
    if haskey(params, :batch_size)
        batch_size = params[:batch_size]
        if batch_size <= 0
            error("batch_size must be positive, got: $batch_size")
        end
        if batch_size < 1000
            @user_warn "Very small batch_size ($batch_size) may cause performance issues"
        end
    end
    
    # Validate min_peptides
    if haskey(params, :min_peptides)
        min_pep = params[:min_peptides]
        if min_pep < 1
            error("min_peptides must be at least 1, got: $min_pep")
        end
    end
    
    return true
end


#==========================================================
General Validation Functions
==========================================================#


