"""
Modified protein group scoring utilities that integrate ML-based scoring
"""

include("../../../utils/ML/proteinGroupScoring.jl")

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
    
    # If ML scoring is requested, we need to collect all PSMs with CV fold information
    all_psms_df = DataFrame()
    
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
        psms_df[!, :ms_file_idx] = ms_file_idx
        
        # Add CV fold information from precursors
        psms_df[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psms_df.precursor_idx]
        
        # Add sequence and charge for grouping
        psms_df[!, :sequence] = [sequences[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :charge] = [getCharge(precursors)[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :proteins] = [accession_numbers[pid] for pid in psms_df.precursor_idx]
        psms_df[!, :entrapment_group_id] = [entrap_ids[pid] for pid in psms_df.precursor_idx]
        
        # Collect for ML scoring if enabled
        if use_ml_scoring
            append!(all_psms_df, psms_df)
        end
        
        # Store initial results
        writeArrow(file_path, psms_df)
        run_to_protein_groups[ms_file_idx] = protein_groups
        
        # Update max scores
        for (k, v) in pairs(protein_groups)
            old = get(acc_to_max_pg_score, k, -Inf32)
            acc_to_max_pg_score[k] = max(v.pg_score, old)
        end
    end
    
    # Apply ML-based protein scoring if requested
    if use_ml_scoring && nrow(all_psms_df) > 0
        @info "Applying ML-based protein group scoring..."
        
        # Collect all protein groups across runs
        all_protein_groups = Dictionary{
            @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
            @NamedTuple{pg_score::Float32, peptides::Set{String}}
        }()
        
        for (_, protein_groups) in run_to_protein_groups
            for (k, v) in pairs(protein_groups)
                if haskey(all_protein_groups, k)
                    # Merge peptides and use max score
                    existing = all_protein_groups[k]
                    merged_peptides = union(existing.peptides, v.peptides)
                    max_score = max(existing.pg_score, v.pg_score)
                    all_protein_groups[k] = (pg_score = max_score, peptides = merged_peptides)
                else
                    all_protein_groups[k] = v
                end
            end
        end
        
        # Apply ML scoring
        protein_features, models = integrate_protein_scoring!(
            all_psms_df,
            all_protein_groups,
            precursors;
            n_top_precursors = n_top_precursors,
            num_parallel_tree = num_parallel_tree
        )
        
        # Export protein features for debugging/analysis
        protein_features_path = joinpath(temp_folder, "protein_features.arrow")
        export_protein_features(protein_features, protein_features_path)
        
        # Update protein groups in each run with ML scores
        for (ms_file_idx, protein_groups) in run_to_protein_groups
            for (k, v) in pairs(protein_groups)
                if haskey(all_protein_groups, k)
                    ml_data = all_protein_groups[k]
                    # Update with ML score
                    protein_groups[k] = (
                        pg_score = ml_data.pg_score,
                        peptides = v.peptides  # Keep run-specific peptides
                    )
                end
            end
        end
        
        # Update max scores with ML scores
        for (k, v) in pairs(all_protein_groups)
            acc_to_max_pg_score[k] = v.pg_score
        end
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