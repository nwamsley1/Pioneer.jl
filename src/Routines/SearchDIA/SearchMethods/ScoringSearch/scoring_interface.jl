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
        # Generate output filename - keep original basename
        base_name = basename(file_path(ref))
        output_path = joinpath(output_dir, base_name)
        
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
            
            # Filter by q-value threshold - use AND logic (both conditions must be met)
            mask = (batch.global_qval .<= q_value_threshold) .&
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
    sort_and_filter_quant_tables_refs(psm_refs::Vector{PSMFileReference},
                                     isotope_trace_type::IsotopeTraceType,
                                     prob_col::Symbol,
                                     best_traces::Set{@NamedTuple})
                                     
Sort and filter PSM tables using file references.
Modifies files in place to optimize storage.
"""
function sort_and_filter_quant_tables_refs(
    psm_refs::Vector{PSMFileReference},
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}
)
    @warn """
    sort_and_filter_quant_tables_refs is deprecated!
    
    Use the Pipeline API directly for better clarity:
    
    pipeline = TransformPipeline() |>
        add_best_trace_indicator(isotope_trace_type, best_traces) |>
        rename_column(:prob, :trace_prob) |>
        select_columns(necessary_cols) |>
        filter_rows(row -> row.best_trace) |>
        remove_columns(:best_trace) |>
        sort_by([prob_col, :target], rev=[true, true])
    
    apply_pipeline!(psm_refs, pipeline)
    """
    
    @info "[PERF] sort_and_filter_quant_tables_refs: Starting" n_files=length(psm_refs)
    start_time = time()
    total_rows_processed = 0
    total_rows_kept = 0
    files_processed = 0
    
    for ref in psm_refs
        if !exists(ref)
            continue
        end
        
        # Load file
        initial_rows = row_count(ref)
        total_rows_processed += initial_rows
        files_processed += 1
        
        # Apply transformation
        # Use a local variable to track rows kept for this file
        local rows_kept_this_file = 0
        
        transform_and_write!(ref) do psms_table
            if seperateTraces(isotope_trace_type)
                # Add indicator variable for best trace
                psms_table[!,:best_trace] = zeros(Bool, nrow(psms_table))
                for i in 1:nrow(psms_table)
                    key = (precursor_idx = psms_table[i, :precursor_idx], 
                          isotopes_captured = psms_table[i, :isotopes_captured])
                    if key ∈ best_traces
                        psms_table[i,:best_trace] = true
                    end
                end
            else
                transform!(
                    groupby(psms_table, :precursor_idx),
                    :prob => (p -> p .== maximum(p)) => :best_trace
                )
            end
            rename!(psms_table, :prob => :trace_prob)
                        # Select only necessary columns to prevent column proliferation
            necessary_cols = [
                :precursor_idx,
                :global_prob,
                :prec_prob,
                :trace_prob,
                :global_qval,
                :run_specific_qval,
                :prec_mz,
                #:pep,
                :weight,
                :target,
                :rt,
                :irt_obs,
                :missed_cleavage,
                :isotopes_captured,
                :scan_idx,
                :entrapment_group_id,
                :ms_file_idx,
                :best_trace  # Temporarily keep for filtering
            ]
            
            # Only keep columns that exist in the table
            available_cols = intersect(necessary_cols, Symbol.(names(psms_table)))
            select!(psms_table, available_cols)
            
            # Filter out unused traces
            filter!(x->x.best_trace, psms_table)
            rows_kept_this_file = nrow(psms_table)
            
            # Remove best_trace column after filtering
            select!(psms_table, Not(:best_trace))
            
            # Sort in descending order of probability
            sort!(psms_table, [prob_col, :target], rev = [true, true], alg=QuickSort)
            
            return psms_table
        end
        
        # Update the total after transformation completes
        total_rows_kept += rows_kept_this_file
        
        # Update reference sort state
        mark_sorted!(ref, prob_col, :target)
    end
    
    elapsed = time() - start_time
    @info "[PERF] sort_and_filter_quant_tables_refs: Completed" elapsed=round(elapsed, digits=3) files_processed=files_processed total_rows=total_rows_processed kept_rows=total_rows_kept retention_rate=round(total_rows_kept/max(total_rows_processed,1), digits=3)
    
    return psm_refs
end

"""
    add_global_scores_to_psms_refs(psm_refs::Vector{PSMFileReference},
                                   acc_to_max_pg_score::Dict{ProteinKey,Float32})

Add global_pg_score to PSM files using references.

# Arguments
- `psm_refs`: PSM file references
- `acc_to_max_pg_score`: Dictionary mapping protein keys to global scores
"""
function add_global_scores_to_psms_refs(
    psm_refs::Vector{PSMFileReference},
    acc_to_max_pg_score::Dict{ProteinKey,Float32}
)
    fifth_pass_start = time()
    @info "[PERF] add_global_scores_to_psms_refs: Starting"
    
    for ref in psm_refs
        if !exists(ref)
            continue
        end
        
        # Add global_pg_score column
        add_column_to_file!(ref, :global_pg_score) do batch
            global_pg_scores = Vector{Float32}(undef, nrow(batch))
            
            for i in 1:nrow(batch)
                key = ProteinKey(
                    batch[i, :inferred_protein_group],
                    batch[i, :target],
                    batch[i, :entrapment_group_id]
                )
                global_pg_scores[i] = get(acc_to_max_pg_score, key, 0.0f0)
            end
            
            global_pg_scores
        end
    end
    
    fifth_pass_time = time() - fifth_pass_start
    @info "[PERF] add_global_scores_to_psms_refs: Completed" elapsed=round(fifth_pass_time, digits=3)
end

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
    get_proteins_passing_qval_refs(pg_refs::Vector{ProteinGroupFileReference},
                                   global_qval_interp::Interpolations.Extrapolation,
                                   experiment_wide_qval_interp::Interpolations.Extrapolation,
                                   global_prob_col::Symbol,
                                   experiment_wide_prob_col::Symbol,
                                   global_qval_col::Symbol,
                                   experiment_wide_qval_col::Symbol,
                                   q_val_threshold::Float32)

Add q-values and passing status to protein group files using references.

# Arguments
- `pg_refs`: Protein group file references
- `global_qval_interp`: Interpolation function for global q-values
- `experiment_wide_qval_interp`: Interpolation function for experiment-wide q-values
- `global_prob_col`: Column name for global probability scores
- `experiment_wide_prob_col`: Column name for experiment-wide probability scores
- `global_qval_col`: Column name for global q-values
- `experiment_wide_qval_col`: Column name for experiment-wide q-values
- `q_val_threshold`: Q-value threshold for filtering
"""
function get_proteins_passing_qval_refs(
    pg_refs::Vector{ProteinGroupFileReference},
    global_qval_interp::Interpolations.Extrapolation,
    experiment_wide_qval_interp::Interpolations.Extrapolation,
    global_prob_col::Symbol,
    experiment_wide_prob_col::Symbol,
    global_qval_col::Symbol,
    experiment_wide_qval_col::Symbol,
    q_val_threshold::Float32
)
    for ref in pg_refs
        if !exists(ref)
            continue
        end
        
        # Transform the file to add q-values and passing status
        transform_and_write!(ref) do passing_proteins
            # Add q-value columns
            passing_proteins[!,global_qval_col] = global_qval_interp.(passing_proteins[!,global_prob_col])
            passing_proteins[!,experiment_wide_qval_col] = experiment_wide_qval_interp.(passing_proteins[!,experiment_wide_prob_col])
            
            # Add a column to mark which protein groups pass the q-value threshold
            # Instead of filtering, we preserve all groups for lookups
            passing_proteins[!, :passes_qval] = (passing_proteins[!,global_qval_col] .<= q_val_threshold) .& 
                                                (passing_proteins[!,experiment_wide_qval_col] .<= q_val_threshold)
            
            return passing_proteins
        end
    end
    
    return nothing
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
export process_psms_for_scoring, update_psms_with_protein_scores!,
       merge_protein_groups, filter_protein_groups_by_qvalue,
       calculate_global_protein_scores, add_global_scores_to_psms!,
       calculate_protein_qvalues, generate_scoring_summary,
       validate_scoring_results,
       process_and_filter_psms, merge_psm_files, filter_psms_by_qvalue,
       calculate_and_add_global_scores!, apply_probit_scores!,
       # Pipeline operations
       add_best_trace_indicator, get_quant_necessary_columns