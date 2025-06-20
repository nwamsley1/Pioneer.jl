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
       # Pipeline operations
       add_best_trace_indicator, get_quant_necessary_columns