"""
Search result reference types for managing method outputs.

These types provide structured references to search method outputs:
- Enforce relationships between PSM and protein files
- Track processing state across methods
- Enable type-safe method communication
"""

using Arrow, DataFrames

include("FileReferences.jl")

#==========================================================
Search Result References
==========================================================#

"""
    ScoringSearchResultRefs

References to all outputs from ScoringSearch method.
Maintains paired PSM-protein files and merged scores.
"""
struct ScoringSearchResultRefs
    paired_files::Vector{PairedSearchFiles}  # Enforces 1:1 pairing
    merged_scores_ref::Union{Nothing,PSMFileReference}
    
    # Constructor with validation
    function ScoringSearchResultRefs(paired_files::Vector{PairedSearchFiles}, 
                                   merged_scores_ref::Union{Nothing,PSMFileReference}=nothing)
        # Validate paired files are consistent
        for (i, paired) in enumerate(paired_files)
            if paired.ms_file_idx != i
                @warn "Paired file index mismatch" expected=i actual=paired.ms_file_idx
            end
        end
        new(paired_files, merged_scores_ref)
    end
end

# Helper accessors
get_psm_refs(r::ScoringSearchResultRefs) = [p.psm_ref for p in r.paired_files]
get_protein_refs(r::ScoringSearchResultRefs) = [p.protein_ref for p in r.paired_files]
get_paired_file(r::ScoringSearchResultRefs, idx::Int) = r.paired_files[idx]
num_files(r::ScoringSearchResultRefs) = length(r.paired_files)

"""
    MaxLFQSearchResultRefs

References to all outputs from MaxLFQSearch method.
Tracks input files and quantification outputs.
"""
struct MaxLFQSearchResultRefs
    input_paired_files::Vector{PairedSearchFiles}
    protein_quant_ref::Union{Nothing,ProteinGroupFileReference}
    precursors_long_ref::Union{Nothing,PSMFileReference}
    precursors_wide_ref::Union{Nothing,PSMFileReference}
    proteins_long_ref::Union{Nothing,ProteinGroupFileReference} 
    proteins_wide_ref::Union{Nothing,ProteinGroupFileReference}
    
    # Constructor with defaults
    function MaxLFQSearchResultRefs(input_paired_files::Vector{PairedSearchFiles};
                                   protein_quant_ref=nothing,
                                   precursors_long_ref=nothing,
                                   precursors_wide_ref=nothing,
                                   proteins_long_ref=nothing,
                                   proteins_wide_ref=nothing)
        new(input_paired_files, protein_quant_ref, precursors_long_ref, 
            precursors_wide_ref, proteins_long_ref, proteins_wide_ref)
    end
end

# Helper accessors
get_input_psm_refs(r::MaxLFQSearchResultRefs) = [p.psm_ref for p in r.input_paired_files]
get_input_protein_refs(r::MaxLFQSearchResultRefs) = [p.protein_ref for p in r.input_paired_files]

#==========================================================
Validation Functions
==========================================================#

"""
    validate_scoring_output(refs::ScoringSearchResultRefs)

Validate that ScoringSearch outputs have expected structure.
"""
function validate_scoring_output(refs::ScoringSearchResultRefs)
    # Check each paired file
    for paired in refs.paired_files
        # Validate pairing
        validate_pairing(paired)
        
        # Validate PSM schema has required columns
        validate_scoring_psm_schema(paired.psm_ref)
        
        # Validate protein schema has required columns  
        validate_scoring_protein_schema(paired.protein_ref)
    end
    
    # Validate merged scores if present
    if !isnothing(refs.merged_scores_ref)
        validate_exists(refs.merged_scores_ref)
        required_cols = Set([:precursor_idx, :global_prob, :prec_prob, :target])
        validate_required_columns(refs.merged_scores_ref.schema, required_cols)
    end
    
    return true
end

"""
    validate_pairing(paired::PairedSearchFiles)

Validate that paired PSM and protein files are consistent.
"""
function validate_pairing(paired::PairedSearchFiles)
    psm_exists = paired.psm_ref.file_exists
    pg_exists = paired.protein_ref.file_exists
    
    if psm_exists != pg_exists
        error("Inconsistent pairing for file $(paired.ms_file_idx): " *
              "PSM exists=$psm_exists, PG exists=$pg_exists")
    end
    
    # If files exist, could validate cross-references
    # (e.g., protein groups in PSMs match those in protein file)
    
    return true
end

"""
    validate_scoring_psm_schema(psm_ref::PSMFileReference)

Validate PSM file has all columns required after scoring.
"""
function validate_scoring_psm_schema(psm_ref::PSMFileReference)
    required = Set([
        :precursor_idx, :prob, :prec_prob, :global_prob,
        :target, :entrapment_group_id, 
        :inferred_protein_group, :pg_score, :global_pg_score,
        :pg_qval, :global_qval_pg, :use_for_protein_quant
    ])
    validate_required_columns(psm_ref.schema, required)
end

"""
    validate_scoring_protein_schema(protein_ref::ProteinGroupFileReference)

Validate protein file has all columns required after scoring.
"""
function validate_scoring_protein_schema(protein_ref::ProteinGroupFileReference)
    required = Set([
        :protein_name, :target, :entrapment_group_id,
        :pg_score, :global_pg_score, :n_peptides,
        :pg_qval, :global_pg_qval
    ])
    validate_required_columns(protein_ref.schema, required)
end

"""
    validate_maxlfq_input(refs::Vector{PSMFileReference})

Validate PSM files are ready for MaxLFQ processing.
"""
function validate_maxlfq_input(refs::Vector{PSMFileReference})
    for ref in refs
        # Check file exists
        validate_exists(ref)
        
        # Check required columns
        required = Set([
            :inferred_protein_group, :target, :entrapment_group_id,
            :precursor_idx, :peak_area, :use_for_protein_quant,
            :pg_qval, :global_qval_pg
        ])
        validate_required_columns(ref.schema, required)
        
        # Check sorted correctly
        required_sort = (:inferred_protein_group, :target, :entrapment_group_id, :precursor_idx)
        if !is_sorted_by(ref, required_sort...)
            error("PSM file $(ref.file_path) must be sorted by $required_sort for MaxLFQ")
        end
    end
    
    return true
end

#==========================================================
Utility Functions
==========================================================#

"""
    create_result_refs_from_paths(psm_paths::Vector{String}, 
                                 protein_paths::Vector{String}) -> ScoringSearchResultRefs

Create ScoringSearchResultRefs from file paths.
"""
function create_result_refs_from_paths(psm_paths::Vector{String}, 
                                     protein_paths::Vector{String})
    length(psm_paths) == length(protein_paths) || 
        error("PSM and protein paths must have same length")
    
    paired_files = PairedSearchFiles[]
    for (idx, (psm_path, protein_path)) in enumerate(zip(psm_paths, protein_paths))
        push!(paired_files, PairedSearchFiles(psm_path, protein_path, idx))
    end
    
    return ScoringSearchResultRefs(paired_files)
end

"""
    describe_results(refs::ScoringSearchResultRefs)

Print summary of scoring search results.
"""
function describe_results(refs::ScoringSearchResultRefs)
    println("ScoringSearchResultRefs:")
    println("  Number of paired files: $(num_files(refs))")
    
    for (i, paired) in enumerate(refs.paired_files)
        println("  File $i:")
        println("    PSM: $(paired.psm_ref.file_path) ($(paired.psm_ref.row_count) rows)")
        println("    Protein: $(paired.protein_ref.file_path) ($(paired.protein_ref.row_count) rows)")
    end
    
    if !isnothing(refs.merged_scores_ref)
        println("  Merged scores: $(refs.merged_scores_ref.file_path)")
    end
end

"""
    describe_results(refs::MaxLFQSearchResultRefs)

Print summary of MaxLFQ search results.
"""
function describe_results(refs::MaxLFQSearchResultRefs)
    println("MaxLFQSearchResultRefs:")
    println("  Number of input files: $(length(refs.input_paired_files))")
    
    if !isnothing(refs.precursors_long_ref)
        println("  Precursors long: $(refs.precursors_long_ref.file_path)")
    end
    
    if !isnothing(refs.proteins_long_ref)
        println("  Proteins long: $(refs.proteins_long_ref.file_path)")
    end
end

# Export all public types and functions
export ScoringSearchResultRefs, MaxLFQSearchResultRefs,
       get_psm_refs, get_protein_refs, get_paired_file, num_files,
       get_input_psm_refs, get_input_protein_refs,
       validate_scoring_output, validate_pairing, validate_maxlfq_input,
       validate_scoring_psm_schema, validate_scoring_protein_schema,
       create_result_refs_from_paths, describe_results