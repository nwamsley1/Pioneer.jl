"""
Modified protein group scoring utilities that integrate ML-based scoring

Note: This now uses the standalone getProteinGroupsDict and writeProteinGroups 
functions from utils.jl instead of duplicating them.

Uses memory-efficient out-of-memory (OOM) approach for ML protein scoring
to handle large datasets without loading all protein groups into memory.
"""

"""
    get_protein_groups_with_ml(
        passing_psms_paths::Vector{String},
        passing_pg_paths::Vector{String},
        protein_groups_folder::String,
        temp_folder::String,
        precursors::LibraryPrecursors;
        min_peptides = 2,
        protein_q_val_threshold::Float32 = 0.01f0,
        use_ml_scoring::Bool = true,
        n_top_precursors::Int = 5,
        num_parallel_tree::Int = 100
    )

Enhanced version of get_protein_groups that optionally uses ML-based protein scoring.
"""
function get_protein_groups_with_ml(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    protein_groups_folder::String,
    temp_folder::String,
    precursors::LibraryPrecursors;
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0,
    use_ml_scoring::Bool = true,
    n_top_precursors::Int = 5,
    num_parallel_tree::Int = 100
)
    # First, run the standard protein inference and initial scoring
    # This is mostly copied from the original get_protein_groups function
    
    # Load all passing PSMs
    passing_psms = Arrow.Table(passing_psms_paths)
    
    # Build protein_peptide_rows using data from PSMs
    protein_peptide_rows = Set{NamedTuple{(:sequence, :protein_name, :decoy, :entrap_id), Tuple{String, String, Bool, UInt8}}}()
    
    # Get data from PSMs
    passing_precursor_idx = passing_psms[:precursor_idx]
    
    # Get other data from precursors
    accession_numbers = getAccessionNumbers(precursors)
    decoys = getIsDecoy(precursors)
    entrap_ids = getEntrapmentGroupId(precursors)
    sequences = getSequence(precursors)
    
    for pid in passing_precursor_idx
        push!(
            protein_peptide_rows, 
            (
                sequence = sequences[pid],
                protein_name = accession_numbers[pid],
                decoy = decoys[pid],
                entrap_id = entrap_ids[pid]
            )
        )
    end
    
    protein_peptide_rows = collect(protein_peptide_rows)
    peptides = [row.sequence for row in protein_peptide_rows]
    proteins = [(protein_name = row.protein_name, decoy = row.decoy, entrap_id = row.entrap_id) for row in protein_peptide_rows]
    
    # Perform protein inference
    protein_inference_dict = infer_proteins(proteins, peptides)
    
    # Initialize storage
    acc_to_max_pg_score = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Float32}()
    run_to_protein_groups = Dict{UInt64, Dictionary}()
    
    # First pass to compute initial protein group scores
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        _, extension = splitext(file_path)
        if extension != ".arrow"
            continue
        end
        
        psms_table = Arrow.Table(file_path)
        
        # Get initial protein groups using standard scoring
        pg_score, inferred_protein_group_names, protein_groups = getProteinGroupsDict(
            protein_inference_dict,
            psms_table[:precursor_idx],
            psms_table[:prob],
            psms_table[:target],
            psms_table[:entrapment_group_id],
            precursors;
            min_peptides = min_peptides
        )
        
        # Convert to DataFrame and add columns
        psms_df = DataFrame(Tables.columntable(psms_table))
        psms_df[!, :pg_score] = pg_score
        psms_df[!, :inferred_protein_group] = inferred_protein_group_names
        psms_df[!, :ms_file_idx] .= ms_file_idx
        # Add CV fold information from precursors
        psms_df[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psms_df.precursor_idx]
        
        # Add sequence and charge for grouping
        psms_df[!, :sequence] = [sequences[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :charge] = [getCharge(precursors)[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :proteins] = [accession_numbers[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :entrapment_group_id] = [entrap_ids[pid] for pid in psms_df.precursor_idx]
        
        # Store initial results
        writeArrow(file_path, psms_df)
        run_to_protein_groups[ms_file_idx] = protein_groups
        
        # Update max scores
        for (k, v) in pairs(protein_groups)
            old = get(acc_to_max_pg_score, k, -Inf32)
            acc_to_max_pg_score[k] = max(v.pg_score, old)
        end
    end
    
    # Write initial protein group files BEFORE ML scoring
    @info "Writing initial protein group files..."
    pg_count = 0
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        _, extension = splitext(file_path)
        if extension != ".arrow"
            continue
        end
        
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        passing_pg_paths[ms_file_idx] = protein_groups_path
        protein_groups = run_to_protein_groups[ms_file_idx]
        
        # Write protein groups for this file
        pg_count += writeProteinGroups(
            acc_to_max_pg_score,
            protein_groups,
            protein_groups_path
        )
    end
    @info "Wrote $pg_count protein groups to $(length(passing_pg_paths)) files"
    
    # Apply ML-based protein scoring if requested
    if use_ml_scoring
        @info "Applying memory-efficient ML-based protein group scoring..."
        
        # Use OOM approach to avoid loading all protein groups into memory
        apply_ml_protein_scoring_oom!(
            protein_groups_folder,
            passing_psms_paths,
            passing_pg_paths,
            precursors;
            max_proteins_for_training = 50000,
            n_top_precursors = n_top_precursors,
            num_parallel_tree = num_parallel_tree,
            n_rounds = 1,
            subsample = 0.5,
            colsample_bynode = 0.3
        )
        
        # Recompute max scores after ML scoring
        @info "Recomputing global protein group scores after ML scoring..."
        acc_to_max_pg_score = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Float32}()
        
        # Read updated protein group files to get ML scores
        for pg_path in passing_pg_paths
            if isfile(pg_path)
                pg_table = Arrow.Table(pg_path)
                for i in 1:length(pg_table[1])
                    key = (
                        protein_name = pg_table.protein_name[i],
                        target = pg_table.target[i],
                        entrap_id = pg_table.entrap_id[i]
                    )
                    # Use global_pg_score if available, otherwise pg_score
                    score = hasproperty(pg_table, :global_pg_score) ? 
                           pg_table.global_pg_score[i] : pg_table.pg_score[i]
                    
                    old_score = get(acc_to_max_pg_score, key, -Inf32)
                    acc_to_max_pg_score[key] = max(score, old_score)
                end
            end
        end
        
        @info "Updated $(length(acc_to_max_pg_score)) protein groups with ML scores"
    end
    
    # Second pass to write results with updated scores
    pg_count = 0
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        _, extension = splitext(file_path)
        if extension != ".arrow"
            continue
        end
        
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        passing_pg_paths[ms_file_idx] = protein_groups_path
        protein_groups = run_to_protein_groups[ms_file_idx]
        
        # Update PSMs with global scores
        psms_table = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        psms_table[!, :global_pg_score] = [
            get(acc_to_max_pg_score, (protein_name = prot, target = tgt, entrap_id = entrap_id), 0.0f0)
            for (prot, tgt, entrap_id) in zip(
                psms_table.inferred_protein_group, 
                psms_table.target, 
                psms_table.entrapment_group_id
            )
        ]
        
        writeArrow(file_path, psms_table)
        
        # Write protein groups
        pg_count += writeProteinGroups(
            acc_to_max_pg_score,
            protein_groups,
            protein_groups_path
        )
    end
    
    return protein_inference_dict
end