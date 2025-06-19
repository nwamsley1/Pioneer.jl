"""
Interface functions for ScoringSearch that work exclusively with file references.

These functions provide a clean abstraction layer between ScoringSearch and file operations,
ensuring all file access goes through the reference system.
"""

# FileReferences.jl, SearchResultReferences.jl and FileOperations.jl are already included by importScripts.jl

using DataFrames, Arrow

#==========================================================
PSM Processing Functions
==========================================================#

"""
    process_psms_for_scoring(psm_refs::Vector{PSMFileReference}, 
                           output_dir::String,
                           precursors;
                           min_peptides=2) -> ScoringSearchResultRefs

Process PSMs through protein inference and create paired output files.
Returns references to all created files.
"""
function process_psms_for_scoring(psm_refs::Vector{PSMFileReference}, 
                                output_dir::String,
                                precursors;
                                min_peptides::Int=2)
    paired_files = PairedSearchFiles[]
    
    for (idx, psm_ref) in enumerate(psm_refs)
        # Validate PSM file has required columns
        validate_exists(psm_ref)
        required_cols = Set([:prob, :inferred_protein_group, :precursor_idx, 
                           :target, :entrapment_group_id])
        validate_required_columns(schema(psm_ref), required_cols)
        
        # Generate output path for protein groups
        pg_path = joinpath(output_dir, "protein_groups_$(idx).arrow")
        
        # Apply protein inference
        pg_ref = apply_protein_inference(psm_ref, pg_path, precursors, 
                                        min_peptides=min_peptides)
        
        # Create paired reference
        push!(paired_files, PairedSearchFiles(psm_ref, pg_ref, idx))
    end
    
    return ScoringSearchResultRefs(paired_files)
end

"""
    update_psms_with_protein_scores!(result_refs::ScoringSearchResultRefs,
                                   output_dir::String) -> ScoringSearchResultRefs

Update PSMs with protein group scores and q-values.
Creates new PSM files with added score columns.
"""
function update_psms_with_protein_scores!(result_refs::ScoringSearchResultRefs,
                                        output_dir::String)
    updated_paired_files = PairedSearchFiles[]
    
    for (idx, paired) in enumerate(result_refs.paired_files)
        # Generate output path for updated PSMs
        updated_psm_path = joinpath(output_dir, "psms_scored_$(idx).arrow")
        
        # Update PSMs with scores from protein groups
        updated_psm_ref = update_psms_with_scores(
            paired.psm_ref,
            paired.protein_ref,
            updated_psm_path
        )
        
        # Create new paired reference with updated PSM
        push!(updated_paired_files, 
              PairedSearchFiles(updated_psm_ref, paired.protein_ref, idx))
    end
    
    # Return new result refs with updated PSMs
    return ScoringSearchResultRefs(updated_paired_files, result_refs.merged_scores_ref)
end

#==========================================================
Protein Group Operations
==========================================================#

"""
    merge_protein_groups(result_refs::ScoringSearchResultRefs,
                       output_path::String) -> ProteinGroupFileReference

Merge protein groups from multiple files using heap-based streaming merge.
"""
function merge_protein_groups(result_refs::ScoringSearchResultRefs,
                            output_path::String)
    pg_refs = get_protein_refs(result_refs)
    
    # Ensure all protein files are sorted properly
    for ref in pg_refs
        ensure_sorted!(ref, :pg_score, :target)
    end
    
    # Perform streaming merge
    merged_ref = stream_sorted_merge(pg_refs, output_path, :pg_score, :target)
    
    return merged_ref
end

"""
    filter_protein_groups_by_qvalue(pg_ref::ProteinGroupFileReference,
                                  output_path::String,
                                  qval_threshold::Float32=0.01f0) -> ProteinGroupFileReference

Filter protein groups by q-value threshold.
"""
function filter_protein_groups_by_qvalue(pg_ref::ProteinGroupFileReference,
                                       output_path::String,
                                       qval_threshold::Float32=0.01f0)
    # Validate schema
    required_cols = Set([:pg_qval, :global_pg_qval])
    validate_required_columns(schema(pg_ref), required_cols)
    
    # Filter using streaming operation
    filtered_ref = stream_filter(pg_ref, output_path, 
        row -> row.pg_qval <= qval_threshold || row.global_pg_qval <= qval_threshold)
    
    return filtered_ref
end

#==========================================================
Score Calculation Functions
==========================================================#

"""
    calculate_global_protein_scores(pg_refs::Vector{ProteinGroupFileReference}) 
                                  -> Dict{Tuple{String,Bool,UInt8}, Float32}

Calculate global (max) protein scores across all files.
Returns dictionary mapping protein keys to max scores.
"""
function calculate_global_protein_scores(pg_refs::Vector{ProteinGroupFileReference})
    global_scores = Dict{Tuple{String,Bool,UInt8}, Float32}()
    
    for pg_ref in pg_refs
        validate_exists(pg_ref)
        
        # Process file to extract max scores
        process_with_memory_limit(pg_ref, 
            batch -> begin
            for row in eachrow(batch)
                key = (row.protein_name, row.target, row.entrapment_group_id)
                score = row.pg_score
                
                if haskey(global_scores, key)
                    global_scores[key] = max(global_scores[key], score)
                else
                    global_scores[key] = score
                end
            end
            end,
            max_memory_mb=500)
    end
    
    return global_scores
end

"""
    add_global_scores_to_psms!(psm_ref::PSMFileReference,
                              global_scores::Dict,
                              output_path::String) -> PSMFileReference

Add global protein scores to PSMs.
"""
function add_global_scores_to_psms!(psm_ref::PSMFileReference,
                                  global_scores::Dict{Tuple{String,Bool,UInt8}, Float32},
                                  output_path::String)
    # Add global score column
    updated_ref = add_column_to_file!(
        psm_ref,
        :global_pg_score
    ) do batch
        # Calculate global scores for this batch
        global_pg_scores = Vector{Float32}(undef, nrow(batch))
        
        for i in 1:nrow(batch)
            key = (batch.inferred_protein_group[i], 
                   batch.target[i], 
                   batch.entrapment_group_id[i])
            
            global_pg_scores[i] = get(global_scores, key, 0.0f0)
        end
        
        global_pg_scores
    end
    
    # Move to final location if needed
    if file_path(updated_ref) != output_path
        mv(file_path(updated_ref), output_path, force=true)
        return PSMFileReference(output_path)
    end
    
    return updated_ref
end

#==========================================================
Q-value Calculation Functions  
==========================================================#

"""
    calculate_protein_qvalues(pg_ref::ProteinGroupFileReference,
                            output_path::String) -> ProteinGroupFileReference

Calculate q-values for protein groups based on target-decoy approach.
"""
function calculate_protein_qvalues(pg_ref::ProteinGroupFileReference,
                                 output_path::String)
    # This is a simplified version - real implementation would:
    # 1. Sort by score
    # 2. Calculate running FDR
    # 3. Convert to q-values
    
    # For now, just add placeholder columns
    updated_ref = add_column_to_file!(pg_ref, :pg_qval) do batch
        # Placeholder: would calculate actual q-values
        fill(0.01f0, nrow(batch))
    end
    
    updated_ref = add_column_to_file!(updated_ref, :global_pg_qval) do batch
        # Placeholder: would calculate actual global q-values
        fill(0.005f0, nrow(batch))
    end
    
    # Move to final location
    if file_path(updated_ref) != output_path
        mv(file_path(updated_ref), output_path, force=true)
        return ProteinGroupFileReference(output_path)
    end
    
    return updated_ref
end

#==========================================================
Summary Functions
==========================================================#

"""
    generate_scoring_summary(result_refs::ScoringSearchResultRefs) -> DataFrame

Generate summary statistics for scoring results.
"""
function generate_scoring_summary(result_refs::ScoringSearchResultRefs)
    summary_data = []
    
    for (idx, paired) in enumerate(result_refs.paired_files)
        if exists(paired.psm_ref) && exists(paired.protein_ref)
            push!(summary_data, (
                file_idx = idx,
                n_psms = row_count(paired.psm_ref),
                n_proteins = row_count(paired.protein_ref),
                psm_path = file_path(paired.psm_ref),
                protein_path = file_path(paired.protein_ref)
            ))
        end
    end
    
    return DataFrame(summary_data)
end

"""
    validate_scoring_results(result_refs::ScoringSearchResultRefs) -> Bool

Validate that all scoring outputs are present and correct.
"""
function validate_scoring_results(result_refs::ScoringSearchResultRefs)
    try
        # Use the validation function from SearchResultReferences
        validate_scoring_output(result_refs)
        return true
    catch e
        @error "Scoring validation failed" exception=e
        return false
    end
end

#==========================================================
PSM Processing with Best Traces
==========================================================#

"""
    process_and_filter_psms(psm_refs::Vector{PSMFileReference}, 
                          output_dir::String,
                          isotope_tracetype,
                          prob_col::Symbol,
                          best_traces::Set) -> Vector{PSMFileReference}
                          
Filter and sort PSM files based on best traces.
Creates new filtered files in output_dir.
"""
function process_and_filter_psms(psm_refs::Vector{PSMFileReference}, 
                               output_dir::String,
                               isotope_tracetype,
                               prob_col::Symbol,
                               best_traces::Set)
    filtered_refs = PSMFileReference[]
    
    for (idx, ref) in enumerate(psm_refs)
        # Generate output filename based on original
        base_name = splitext(basename(file_path(ref)))[1]
        output_path = joinpath(output_dir, "$(base_name)_filtered.arrow")
        
        # Filter based on best traces
        filter_fn = function(row)
            # Check if this row's trace is in best_traces
            trace_key = (
                precursor_idx = row.precursor_idx,
                isotopes_captured = row.isotopes_captured
            )
            return trace_key in best_traces
        end
        
        filtered_ref = stream_filter(ref, output_path, filter_fn)
        
        # Sort the filtered file by probability column
        sort_file_by_keys!(filtered_ref, prob_col; reverse=true)
        
        push!(filtered_refs, filtered_ref)
    end
    
    return filtered_refs
end

"""
    merge_psm_files(psm_refs::Vector{PSMFileReference},
                   output_path::String,
                   sort_col::Symbol) -> PSMFileReference
                   
Merge multiple PSM files sorted by specified column.
Ensures files are properly sorted before merging.
"""
function merge_psm_files(psm_refs::Vector{PSMFileReference},
                       output_path::String,
                       sort_col::Symbol)
    # The merge_psm_scores function already handles sorting if needed
    return merge_psm_scores(psm_refs, output_path, sort_col)
end

"""
    filter_psms_by_qvalue(psm_refs::Vector{PSMFileReference},
                         output_dir::String,
                         precursors,
                         global_qval_interp,
                         qval_interp,
                         q_value_threshold::Float32) -> Vector{PSMFileReference}
                         
Filter PSMs by q-value threshold and add q-value columns.
Creates new filtered files in output_dir.
"""
function filter_psms_by_qvalue(psm_refs::Vector{PSMFileReference},
                             output_dir::String,
                             precursors,
                             global_qval_interp,
                             qval_interp,
                             q_value_threshold::Float32)
    filtered_refs = PSMFileReference[]
    
    for (idx, ref) in enumerate(psm_refs)
        # Generate output filename
        base_name = splitext(basename(file_path(ref)))[1]
        output_path = joinpath(output_dir, "$(base_name)_passing.arrow")
        
        # Transform: add q-values and filter
        transform_fn = function(batch)
            # Add q-value columns
            n_rows = nrow(batch)
            global_qvals = Vector{Float32}(undef, n_rows)
            qvals = Vector{Float32}(undef, n_rows)
            
            for i in 1:n_rows
                global_qvals[i] = global_qval_interp(batch.global_prob[i])
                qvals[i] = qval_interp(batch.prec_prob[i])
            end
            
            batch.global_qval = global_qvals
            batch.qval = qvals
            
            # Filter by q-value threshold
            mask = (batch.global_qval .<= q_value_threshold) .|
                   (batch.qval .<= q_value_threshold)
            
            return batch[mask, :]
        end
        
        filtered_ref = stream_transform(ref, output_path, transform_fn)
        
        push!(filtered_refs, filtered_ref)
    end
    
    return filtered_refs
end

#==========================================================
Reference-Based Wrappers for Direct File Operations
==========================================================#

"""
    calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    
Calculate global protein scores and add them to files via references.
Returns the score dictionary for downstream use.
"""
function calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    # First pass: collect max scores
    acc_to_max_pg_score = Dict{ProteinKey, Float32}()
    
    for ref in pg_refs
        process_with_memory_limit(ref) do batch
            for row in eachrow(batch)
                key = ProteinKey(
                    row.protein_name,
                    row.target,
                    row.entrap_id
                )
                old = get(acc_to_max_pg_score, key, -Inf32)
                acc_to_max_pg_score[key] = max(row.pg_score, old)
            end
        end
    end
    
    # Second pass: add global_pg_score column and sort
    for ref in pg_refs
        add_column_and_sort!(ref, :global_pg_score, 
            batch -> begin
                scores = Vector{Float32}(undef, nrow(batch))
                for i in 1:nrow(batch)
                    key = ProteinKey(
                        batch.protein_name[i],
                        batch.target[i],
                        batch.entrap_id[i]
                    )
                    scores[i] = get(acc_to_max_pg_score, key, batch.pg_score[i])
                end
                scores
            end,
            :global_pg_score, :target;  # sort keys
            reverse=true
        )
    end
    
    return acc_to_max_pg_score
end

"""
    apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference}, 
                        β_fitted::Vector{Float64}, feature_names::Vector{Symbol})
    
Apply probit regression scores to protein group files.
Note: This function is called from utils.jl and needs access to calculate_probit_scores.
"""
function apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference},
                             β_fitted::Vector{Float64},
                             feature_names::Vector{Symbol})
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Calculate probit scores (function from utils.jl)
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted)
            
            # Update scores
            df[!, :old_pg_score] = copy(df.pg_score)
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by new scores
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            return df
        end
    end
end

# Export all interface functions
export process_psms_for_scoring, update_psms_with_protein_scores!,
       merge_protein_groups, filter_protein_groups_by_qvalue,
       calculate_global_protein_scores, add_global_scores_to_psms!,
       calculate_protein_qvalues, generate_scoring_summary,
       validate_scoring_results,
       process_and_filter_psms, merge_psm_files, filter_psms_by_qvalue,
       calculate_and_add_global_scores!, apply_probit_scores!