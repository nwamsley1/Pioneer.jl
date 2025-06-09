# Efficient OOM implementation for protein group scoring with integrated ML

"""
    get_protein_groups_efficient(passing_psms_paths::Vector{String}, passing_pg_paths::Vector{String},
                                protein_groups_folder::String, temp_folder::String,
                                precursors::LibraryPrecursors; kwargs...)

Efficient OOM implementation that integrates ML scoring into protein group creation.
Reduces file I/O by combining operations into fewer passes.

# Process
1. Train ML model on sampled data (if enabled)
2. Create protein groups and apply ML scoring in single pass per file
3. Update global scores and filter in final single pass
"""
function get_protein_groups_efficient(
    passing_psms_paths::Vector{String},
    passing_pg_paths::Vector{String},
    protein_groups_folder::String,
    temp_folder::String,
    precursors::LibraryPrecursors;
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0,
    max_psms_in_memory::Int64 = 10000000
)
    @info "Starting efficient protein group analysis with integrated ML scoring..."
    start_time = time()
    
    # Count all possible peptides for each protein in the library
    peptide_count_start = time()
    protein_to_possible_peptides = count_possible_peptides(precursors)
    @info "Peptide counting completed" elapsed=round(time()-peptide_count_start, digits=3) n_proteins=length(protein_to_possible_peptides)
    
    # Perform protein inference
    passing_psms = Arrow.Table(passing_psms_paths)
    protein_peptide_rows = build_protein_peptide_rows(passing_psms, precursors)
    peptides = [row.sequence for row in protein_peptide_rows]
    proteins = [(protein_name = row.protein_name, decoy = row.decoy, entrap_id = row.entrap_id) for row in protein_peptide_rows]
    
    inference_start = time()
    protein_inference_dict = infer_proteins(proteins, peptides)
    @info "Protein inference completed" elapsed=round(time()-inference_start, digits=3)
    
    # Check if ML scoring should be enabled
    total_protein_groups = estimate_total_protein_groups(passing_psms_paths, protein_inference_dict, precursors, min_peptides)
    max_protein_groups_in_memory_limit = 5 * max_psms_in_memory
    use_ml_scoring = total_protein_groups > 10  # Need sufficient data for ML
    
    if use_ml_scoring
        @info "ML protein scoring enabled ($(total_protein_groups) protein groups)"
        
        # Train ML model on sampled data
        ml_model_start = time()
        β_fitted, X_mean, X_std, feature_names = train_ml_model_streaming(
            passing_psms_paths, 
            protein_inference_dict,
            precursors,
            protein_to_possible_peptides,
            min_peptides,
            max_psms_in_memory
        )
        @info "ML model training completed" elapsed=round(time()-ml_model_start, digits=3)
    else
        @info "ML protein scoring disabled (insufficient data)"
        β_fitted = X_mean = X_std = feature_names = nothing
    end
    
    # Process files with integrated scoring
    acc_to_max_pg_score = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Float32}()
    pg_count = 0
    
    # First pass: Create protein groups with ML scoring
    first_pass_start = time()
    @info "Creating protein groups with integrated ML scoring..."
    
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        if !endswith(file_path, ".arrow")
            continue
        end
        
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        passing_pg_paths[ms_file_idx] = protein_groups_path
        
        # Process single file
        pg_count_file = process_file_with_ml_scoring!(
            file_path,
            protein_groups_path,
            ms_file_idx,
            protein_inference_dict,
            precursors,
            protein_to_possible_peptides,
            acc_to_max_pg_score,
            min_peptides,
            use_ml_scoring,
            β_fitted,
            X_mean,
            X_std,
            feature_names
        )
        
        pg_count += pg_count_file
    end
    
    @info "First pass completed" elapsed=round(time()-first_pass_start, digits=3) total_protein_groups=pg_count
    
    # Second pass: Update global scores and write final PSMs
    second_pass_start = time()
    @info "Updating global scores in PSM files..."
    
    for (ms_file_idx, file_path) in enumerate(passing_psms_paths)
        if !endswith(file_path, ".arrow")
            continue
        end
        
        # Read PSMs
        psms_table = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        
        # Update global pg scores
        psms_table[!,:global_pg_score] = [
            get(acc_to_max_pg_score, 
                (protein_name = prot, target = tgt, entrap_id = entrap_id), 
                0.0f0)
            for (prot, tgt, entrap_id) in zip(
                psms_table.inferred_protein_group, 
                psms_table.target, 
                psms_table.entrapment_group_id
            )
        ]
        
        # Write back
        writeArrow(file_path, psms_table)
    end
    
    @info "Second pass completed" elapsed=round(time()-second_pass_start, digits=3)
    
    # Final step: Apply filtering (this modifies files in place)
    if pg_count > 0
        filter_protein_groups_by_qval!(
            passing_pg_paths,
            acc_to_max_pg_score,
            protein_q_val_threshold
        )
    end
    
    total_elapsed = time() - start_time
    @info "Efficient protein group analysis completed" total_elapsed=round(total_elapsed, digits=3) total_protein_groups=pg_count
    
    if use_ml_scoring
        @info "ML-based protein scoring applied. The pg_score and pg_qval columns contain ML-enhanced scores."
    end
    
    return protein_inference_dict
end

"""
    process_file_with_ml_scoring!(...)

Process a single PSM file: create protein groups, apply ML scoring if enabled, and write results.
Returns the number of protein groups created.
"""
function process_file_with_ml_scoring!(
    psm_file_path::String,
    protein_groups_path::String,
    ms_file_idx::Int64,
    protein_inference_dict::Dictionary,
    precursors::LibraryPrecursors,
    protein_to_possible_peptides::Dict,
    acc_to_max_pg_score::Dict,
    min_peptides::Int64,
    use_ml_scoring::Bool,
    β_fitted,
    X_mean,
    X_std,
    feature_names
)
    # Read PSMs
    psms_table = Arrow.Table(psm_file_path)
    
    # Create protein groups
    pg_score, inferred_protein_group_names, protein_groups = getProteinGroupsDict(
        protein_inference_dict,
        psms_table[:precursor_idx],
        psms_table[:prob],
        psms_table[:target],
        psms_table[:entrapment_group_id],
        precursors;
        min_peptides = min_peptides
    )
    
    # Update PSMs with protein group info
    psms_df = DataFrame(Tables.columntable(psms_table))
    psms_df[!,:pg_score] = pg_score
    psms_df[!,:inferred_protein_group] = inferred_protein_group_names
    writeArrow(psm_file_path, psms_df)
    
    # Create protein groups DataFrame
    keys_array = keys(protein_groups)
    values_array = values(protein_groups)
    
    pg_df = DataFrame((
        protein_name = [k[:protein_name] for k in keys_array],
        target = [k[:target] for k in keys_array],
        entrap_id = [k[:entrap_id] for k in keys_array],
        pg_score = [v[:pg_score] for v in values_array],
        n_peptides = [length(unique(v[:peptides])) for v in values_array],
        total_peptide_length = [sum(length(pep) for pep in v[:peptides]) for v in values_array]
    ))
    
    # Calculate additional features
    add_protein_group_features!(pg_df, protein_to_possible_peptides)
    
    # Apply ML scoring if enabled
    if use_ml_scoring && nrow(pg_df) > 0 && !isnothing(β_fitted)
        # Add ML features
        add_feature_columns!(pg_df)
        
        # Calculate ML scores
        X = Matrix{Float64}(pg_df[:, feature_names])
        ml_scores = calculate_probit_scores(X, β_fitted, X_mean, X_std)
        ml_qvalues = calculate_qvalues_from_scores(ml_scores, pg_df.target)
        
        # Overwrite pg_score and pg_qval with ML values
        pg_df[!, :pg_score] = Float32.(ml_scores)
        pg_df[!, :pg_qval] = Float32.(ml_qvalues)
    else
        # Calculate traditional q-values
        pg_qvalues = calculate_qvalues_from_scores(pg_df.pg_score, pg_df.target)
        pg_df[!, :pg_qval] = Float32.(pg_qvalues)
    end
    
    # Update max scores dictionary
    for (i, key) in enumerate(keys_array)
        old_score = get(acc_to_max_pg_score, key, -Inf32)
        acc_to_max_pg_score[key] = max(pg_df[i, :pg_score], old_score)
    end
    
    # Add placeholders for global scores (will be updated in second pass)
    pg_df[!, :global_pg_score] = pg_df.pg_score  # Temporary
    pg_df[!, :global_pg_qval] = pg_df.pg_qval    # Temporary
    
    # Sort by score and write
    sort!(pg_df, :pg_score, rev = true)
    writeArrow(protein_groups_path, pg_df)
    
    return nrow(pg_df)
end

"""
    train_ml_model_streaming(...)

Train ML model by streaming through files and sampling protein groups.
"""
function train_ml_model_streaming(
    passing_psms_paths::Vector{String},
    protein_inference_dict::Dictionary,
    precursors::LibraryPrecursors,
    protein_to_possible_peptides::Dict,
    min_peptides::Int64,
    max_samples::Int64
)
    # Calculate samples per file
    n_files = count(endswith(path, ".arrow") for path in passing_psms_paths)
    samples_per_file = max(100, max_samples ÷ n_files)
    
    sampled_protein_groups = DataFrame()
    
    for file_path in passing_psms_paths
        if !endswith(file_path, ".arrow")
            continue
        end
        
        # Read PSMs
        psms_table = Arrow.Table(file_path)
        
        # Create protein groups for this file
        _, _, protein_groups = getProteinGroupsDict(
            protein_inference_dict,
            psms_table[:precursor_idx],
            psms_table[:prob],
            psms_table[:target],
            psms_table[:entrapment_group_id],
            precursors;
            min_peptides = min_peptides
        )
        
        if length(protein_groups) == 0
            continue
        end
        
        # Convert to DataFrame
        keys_array = keys(protein_groups)
        values_array = values(protein_groups)
        
        file_pg_df = DataFrame((
            protein_name = [k[:protein_name] for k in keys_array],
            target = [k[:target] for k in keys_array],
            entrap_id = [k[:entrap_id] for k in keys_array],
            pg_score = [v[:pg_score] for v in values_array],
            n_peptides = [length(unique(v[:peptides])) for v in values_array],
            total_peptide_length = [sum(length(pep) for pep in v[:peptides]) for v in values_array]
        ))
        
        # Add features
        add_protein_group_features!(file_pg_df, protein_to_possible_peptides)
        
        # Sample if needed
        if nrow(file_pg_df) > samples_per_file
            sample_indices = sample(1:nrow(file_pg_df), samples_per_file, replace=false)
            file_pg_df = file_pg_df[sample_indices, :]
        end
        
        append!(sampled_protein_groups, file_pg_df)
    end
    
    # Add ML features
    add_feature_columns!(sampled_protein_groups)
    
    # Define features
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_binom_coeff]
    X = Matrix{Float64}(sampled_protein_groups[:, feature_names])
    y = sampled_protein_groups.target
    
    # Fit model
    β_fitted, X_mean, X_std = fit_probit_model(X, y)
    
    @info "ML model trained on $(nrow(sampled_protein_groups)) sampled protein groups"
    
    return β_fitted, X_mean, X_std, feature_names
end

"""
    filter_protein_groups_by_qval!(pg_paths, max_scores_dict, q_val_threshold)

Final pass to update global scores and apply FDR filtering.
Modifies files in place.
"""
function filter_protein_groups_by_qval!(
    pg_paths::Vector{String},
    max_scores_dict::Dict,
    q_val_threshold::Float32
)
    @info "Applying FDR filtering to protein groups..."
    
    # First, calculate global q-values from max scores
    all_max_scores = Float32[]
    all_targets = Bool[]
    
    for (key, score) in max_scores_dict
        push!(all_max_scores, score)
        push!(all_targets, key.target)
    end
    
    global_qvalues = calculate_qvalues_from_scores(all_max_scores, all_targets)
    
    # Create lookup dictionary
    score_to_qval = Dict{Float32, Float32}()
    for (i, score) in enumerate(all_max_scores)
        score_to_qval[score] = global_qvalues[i]
    end
    
    # Update files with global scores and filter
    total_kept = 0
    total_processed = 0
    
    for pg_path in pg_paths
        if !isfile(pg_path) || !endswith(pg_path, ".arrow")
            continue
        end
        
        # Read
        df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
        total_processed += nrow(df)
        
        # Update global scores
        for i in 1:nrow(df)
            key = (protein_name = df[i, :protein_name], 
                   target = df[i, :target], 
                   entrap_id = df[i, :entrap_id])
            max_score = max_scores_dict[key]
            df[i, :global_pg_score] = max_score
            df[i, :global_pg_qval] = get(score_to_qval, max_score, 1.0f0)
        end
        
        # Filter by q-value
        filter!(x -> x.pg_qval <= q_val_threshold && x.global_pg_qval <= q_val_threshold, df)
        total_kept += nrow(df)
        
        # Write back
        writeArrow(pg_path, df)
    end
    
    @info "FDR filtering complete" total_processed=total_processed kept=total_kept retention_rate=round(total_kept/max(total_processed,1), digits=3)
end

# Helper functions

function count_possible_peptides(precursors::LibraryPrecursors)
    protein_to_possible_peptides = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Set{String}}()
    
    all_accession_numbers = getAccessionNumbers(precursors)
    all_sequences = getSequence(precursors)
    all_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    
    for i in 1:length(all_accession_numbers)
        protein_names = split(all_accession_numbers[i], ';')
        is_decoy = all_decoys[i]
        entrap_id = all_entrap_ids[i]
        
        for protein_name in protein_names
            key = (protein_name = String(protein_name), target = !is_decoy, entrap_id = entrap_id)
            if !haskey(protein_to_possible_peptides, key)
                protein_to_possible_peptides[key] = Set{String}()
            end
            push!(protein_to_possible_peptides[key], all_sequences[i])
        end
    end
    
    return protein_to_possible_peptides
end

function build_protein_peptide_rows(passing_psms::Arrow.Table, precursors::LibraryPrecursors)
    protein_peptide_rows = Set{NamedTuple{(:sequence, :protein_name, :decoy, :entrap_id), Tuple{String, String, Bool, UInt8}}}()
    
    passing_precursor_idx = passing_psms[:precursor_idx]
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
    
    return collect(protein_peptide_rows)
end

function estimate_total_protein_groups(
    psm_paths::Vector{String},
    protein_inference_dict::Dictionary,
    precursors::LibraryPrecursors,
    min_peptides::Int64
)
    # Quick estimation by sampling first few files
    estimated_total = 0
    files_sampled = 0
    
    for path in psm_paths[1:min(3, length(psm_paths))]
        if endswith(path, ".arrow")
            psms_table = Arrow.Table(path)
            _, _, protein_groups = getProteinGroupsDict(
                protein_inference_dict,
                psms_table[:precursor_idx],
                psms_table[:prob],
                psms_table[:target],
                psms_table[:entrapment_group_id],
                precursors;
                min_peptides = min_peptides
            )
            estimated_total += length(protein_groups)
            files_sampled += 1
        end
    end
    
    if files_sampled > 0
        avg_per_file = estimated_total / files_sampled
        total_files = count(endswith(path, ".arrow") for path in psm_paths)
        return Int(round(avg_per_file * total_files))
    end
    
    return 0
end

function add_protein_group_features!(df::DataFrame, protein_to_possible_peptides::Dict)
    # Calculate n_possible_peptides and peptide_coverage
    n_possible_peptides = zeros(Int64, nrow(df))
    
    for i in 1:nrow(df)
        # Handle protein groups with multiple proteins
        protein_names_in_group = split(df[i, :protein_name], ';')
        
        # Union of all peptide sets from proteins in the group
        all_possible_peptides = Set{String}()
        for individual_protein in protein_names_in_group
            individual_key = (
                protein_name = String(individual_protein), 
                target = df[i, :target], 
                entrap_id = df[i, :entrap_id]
            )
            if haskey(protein_to_possible_peptides, individual_key)
                union!(all_possible_peptides, protein_to_possible_peptides[individual_key])
            end
        end
        
        n_possible_peptides[i] = max(length(all_possible_peptides), 1)
    end
    
    df[!, :n_possible_peptides] = n_possible_peptides
    df[!, :peptide_coverage] = df.n_peptides ./ df.n_possible_peptides
end