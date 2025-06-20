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



#==========================================================
Protein Group Operations
==========================================================#



#==========================================================
Score Calculation Functions
==========================================================#



#==========================================================
Q-value Calculation Functions  
==========================================================#


#==========================================================
Summary Functions
==========================================================#



#==========================================================
PSM Processing with Best Traces
==========================================================#

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
        process_with_memory_limit(ref, 
            batch -> begin
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
        )
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

#==========================================================
Reference-based Sort and Filter Functions
==========================================================#



"""
    update_psms_with_probit_scores_refs(paired_refs::Vector{PairedSearchFiles},
                                       acc_to_max_pg_score::Dict{ProteinKey,Float32},
                                       pg_score_to_qval::Interpolations.Extrapolation,
                                       global_pg_score_to_qval::Interpolations.Extrapolation)

Update PSMs with probit-scored pg_score values and q-values using references.

# Arguments
- `paired_refs`: Paired PSM/protein group file references
- `acc_to_max_pg_score`: Dictionary mapping protein keys to global scores
- `pg_score_to_qval`: Interpolation function for pg_score to q-value
- `global_pg_score_to_qval`: Interpolation function for global_pg_score to q-value
"""
function update_psms_with_probit_scores_refs(
    paired_refs::Vector{PairedSearchFiles},
    acc_to_max_pg_score::Dict{ProteinKey,Float32},
    pg_score_to_qval::Interpolations.Extrapolation,
    global_pg_score_to_qval::Interpolations.Extrapolation
)
    start_time = time()
    @info "[PERF] update_psms_with_probit_scores_refs: Starting"
    
    total_psms_updated = 0
    files_processed = 0
    
    for paired_ref in paired_refs
        psm_ref = paired_ref.psm_ref
        pg_ref = paired_ref.protein_ref
        
        # Verify both files exist
        if !exists(psm_ref) || !exists(pg_ref)
            @warn "Missing file" psm_path=file_path(paired_ref.psm_ref) pg_path=file_path(paired_ref.protein_ref)
            continue
        end
        
        # Load protein groups to get probit-scored pg_score values
        pg_table = Arrow.Table(file_path(pg_ref))
        
        # Create lookup dictionary: ProteinKey -> probit pg_score
        pg_score_lookup = Dict{ProteinKey, Float32}()
        n_pg_rows = length(pg_table[:protein_name])
        
        for i in 1:n_pg_rows
            key = ProteinKey(
                pg_table[:protein_name][i],
                pg_table[:target][i],
                pg_table[:entrap_id][i]
            )
            pg_score_lookup[key] = pg_table[:pg_score][i]  # This is now the probit score
        end
        
        # Transform PSM file
        transform_and_write!(psm_ref) do psms_df
            n_psms = nrow(psms_df)
            
            # Update pg_score column with probit scores
            probit_pg_scores = Vector{Float32}(undef, n_psms)
            global_pg_scores = Vector{Float32}(undef, n_psms)
            pg_qvals = Vector{Float32}(undef, n_psms)
            global_pg_qvals = Vector{Float32}(undef, n_psms)
            
            for i in 1:n_psms
                # Skip if missing inferred protein group
                if ismissing(psms_df[i, :inferred_protein_group])
                    throw("Missing Inferred Protein Group!!!")
                end
                
                # Skip if peptide didn't match to a distinct protein group
                if psms_df[i,:use_for_protein_quant] == false
                    probit_pg_scores[i] = 0.0f0
                    global_pg_scores[i] = 0.0f0
                    pg_qvals[i] = 1.0f0
                    global_pg_qvals[i] = 1.0f0
                    continue
                end
                
                # Create key for lookup
                key = ProteinKey(
                    psms_df[i, :inferred_protein_group],
                    psms_df[i, :target],
                    psms_df[i, :entrapment_group_id]
                )
                
                # Get scores
                if !haskey(pg_score_lookup, key)
                    throw("Missing pg score lookup key!!!")
                end
                probit_pg_scores[i] = pg_score_lookup[key]
                
                if !haskey(acc_to_max_pg_score, key)
                    throw("Missing global pg score lookup key!!!")
                end
                global_pg_scores[i] = acc_to_max_pg_score[key]
                
                # Calculate q-values
                pg_qvals[i] = pg_score_to_qval(probit_pg_scores[i])
                global_pg_qvals[i] = global_pg_score_to_qval(global_pg_scores[i])
            end
            
            # Update columns
            psms_df[!, :pg_score] = probit_pg_scores
            psms_df[!, :global_pg_score] = global_pg_scores
            psms_df[!, :pg_qval] = pg_qvals
            psms_df[!, :global_qval_pg] = global_pg_qvals
            
            total_psms_updated += n_psms
            
            return psms_df
        end
        
        files_processed += 1
    end
    
    elapsed = time() - start_time
    @info "[PERF] update_psms_with_probit_scores_refs: Completed" elapsed=round(elapsed, digits=3) files_processed=files_processed total_psms_updated=total_psms_updated
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

#==========================================================
Scoring-Specific Pipeline Operations
==========================================================#

"""
    add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)

Add best trace indicator based on isotope trace type.
"""
function add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)
    desc = "add_best_trace_indicator"
    op = function(df)  # df is passed by transform_and_write!
        if seperateTraces(isotope_type)
            # Efficient vectorized operation for separate traces
            df[!,:best_trace] = [
                (precursor_idx=row.precursor_idx, 
                 isotopes_captured=row.isotopes_captured) ∈ best_traces
                for row in eachrow(df)
            ]
        else
            # Group-based operation for combined traces
            transform!(groupby(df, :precursor_idx),
                      :prob => (p -> p .== maximum(p)) => :best_trace)
        end
        return df
    end
    return desc => op
end

"""
    get_quant_necessary_columns() -> Vector{Symbol}

Get the standard columns needed for quantification analysis.
"""
function get_quant_necessary_columns()
    return [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :global_qval,
        :run_specific_qval,
        :prec_mz,
        #:pep,  # Commented out as in original
        :weight,
        :target,
        :rt,
        :irt_obs,
        :missed_cleavage,
        :isotopes_captured,
        :scan_idx,
        :entrapment_group_id,
        :ms_file_idx
    ]
end

# Export all interface functions
export merge_psm_files,
       calculate_and_add_global_scores!, apply_probit_scores!,
       update_psms_with_probit_scores_refs,
       # Pipeline operations
       add_best_trace_indicator, get_quant_necessary_columns