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
        :pg_qval, :global_qval_pg, :use_for_protein_quant, :peak_area
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
            @warn "Very small batch_size ($batch_size) may cause performance issues"
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

"""
    check_maxlfq_memory_requirements(ref::PSMFileReference, params::Dict)

Estimate memory requirements for MaxLFQ processing and warn if potentially problematic.
Returns estimated memory usage in MB.
"""
function check_maxlfq_memory_requirements(ref::PSMFileReference, params::Dict)
    validate_exists(ref)
    
    # Get basic file statistics
    n_rows = row_count(ref)
    
    # Estimate number of protein groups and experiments
    # This is rough - in practice would need to sample the file
    estimated_protein_groups = max(1000, div(n_rows, 10))  # Conservative estimate
    estimated_experiments = get(params, :n_experiments, 100)  # Default assumption
    
    # Memory estimation based on maxLFQ.jl pre-allocation pattern
    # Lines 354-375: nrows = nfiles*ngroups, multiple vectors of this size
    batch_memory_mb = (estimated_experiments * estimated_protein_groups * 8 * 10) / (1024 * 1024)  # ~10 columns, 8 bytes each
    
    # Add memory for input data (DataFrame materialization)
    input_memory_mb = (n_rows * 20 * 8) / (1024 * 1024)  # ~20 columns average
    
    total_memory_mb = batch_memory_mb + input_memory_mb
    
    # Issue warnings for potential memory issues
    if total_memory_mb > 8000  # 8GB threshold
        @warn "MaxLFQ processing may require significant memory" estimated_memory_mb=total_memory_mb n_protein_groups=estimated_protein_groups n_experiments=estimated_experiments
        
        if total_memory_mb > 32000  # 32GB threshold
            @warn "Very high memory usage expected - consider reducing batch_size or using streaming approach"
        end
    end
    
    return total_memory_mb
end

#==========================================================
General Validation Functions
==========================================================#

"""
    validate_join_compatibility(left_ref::FileReference, right_ref::FileReference, join_key::Symbol)

Validate that two file references are compatible for joining on the specified key.
"""
function validate_join_compatibility(left_ref::FileReference, right_ref::FileReference, join_key::Symbol)
    # Both files must exist
    validate_exists(left_ref)
    validate_exists(right_ref)
    
    # Both files must have the join key
    has_column(schema(left_ref), join_key) || 
        error("Left file missing join key: $join_key")
    has_column(schema(right_ref), join_key) || 
        error("Right file missing join key: $join_key")
    
    return true
end

"""
    validate_sort_compatibility(refs::Vector{<:FileReference}, sort_keys::Symbol...)

Validate that all file references have the required sort keys for merging.
"""
function validate_sort_compatibility(refs::Vector{<:FileReference}, sort_keys::Symbol...)
    for (i, ref) in enumerate(refs)
        validate_exists(ref)
        
        # Check that all sort keys exist in the schema
        for key in sort_keys
            if !has_column(schema(ref), key)
                error("File $i missing sort key: $key")
            end
        end
    end
    
    return true
end

# Export validation functions
export validate_maxlfq_input, validate_maxlfq_parameters, check_maxlfq_memory_requirements,
       validate_join_compatibility, validate_sort_compatibility