"""
Helper functions for protein inference and group analysis.

This module contains refactored functions that break down the complex
protein inference logic into manageable, testable components.
"""

using Arrow
using DataFrames
using Tables
using Dictionaries

"""
    build_protein_peptide_catalog(precursors::LibraryPrecursors) -> Dictionary{ProteinKey, Set{String}}

Build a catalog of all possible peptides for each protein from the spectral library.

# Arguments
- `precursors`: Library precursor information

# Returns
- Dictionary mapping each protein to its set of possible peptides
"""
function build_protein_peptide_catalog(precursors::LibraryPrecursors)
    catalog = Dictionary{ProteinKey, Set{String}}()
    
    all_accessions = getAccessionNumbers(precursors)
    all_sequences = getSequence(precursors)
    all_is_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    
    @info "[PERF] Building protein-peptide catalog" n_precursors=length(all_accessions)
    start_time = time()
    
    for i in eachindex(all_accessions)
        # Handle proteins with shared peptides (semicolon-delimited)
        protein_names = split(all_accessions[i], ';')
        sequence = all_sequences[i]
        is_target = !all_is_decoys[i]
        entrap_id = all_entrap_ids[i]
        
        for protein_name in protein_names
            key = ProteinKey(String(protein_name), is_target, entrap_id)
            
            if !haskey(catalog, key)
                insert!(catalog, key, Set{String}())
            end
            push!(catalog[key], sequence)
        end
    end
    
    elapsed = time() - start_time
    @info "[PERF] Protein-peptide catalog built" elapsed=round(elapsed, digits=3) n_proteins=length(catalog)
    
    return catalog
end

"""
    extract_peptide_protein_pairs(psms_table, precursors) -> Vector{Tuple{PeptideKey, ProteinKey}}

Extract unique peptide-protein pairs from PSMs for inference.

# Arguments
- `psms_table`: Arrow table of PSMs
- `precursors`: Library precursor information

# Returns
- Vector of (peptide, protein) key pairs
"""
function extract_peptide_protein_pairs(psms_table, precursors)
    pairs = Set{Tuple{PeptideKey, ProteinKey}}()
    
    precursor_indices = psms_table[:precursor_idx]
    accession_numbers = getAccessionNumbers(precursors)
    sequences = getSequence(precursors)
    is_decoys = getIsDecoy(precursors)
    entrap_ids = getEntrapmentGroupId(precursors)
    
    for pid in precursor_indices
        sequence = sequences[pid]
        accessions = accession_numbers[pid]
        is_target = !is_decoys[pid]
        entrap_id = entrap_ids[pid]
        
        peptide_key = PeptideKey(sequence, is_target, entrap_id)
        protein_key = ProteinKey(accessions, is_target, entrap_id)
        
        push!(pairs, (peptide_key, protein_key))
    end
    
    return collect(pairs)
end

"""
    perform_file_protein_inference(psms_table, precursors) -> InferenceResult

Perform protein inference for peptides in a single file.

# Arguments
- `psms_table`: Arrow table of PSMs
- `precursors`: Library precursor information

# Returns
- `InferenceResult` with peptide-to-protein mappings and quantification flags
"""
function perform_file_protein_inference(psms_table, precursors)::InferenceResult
    # Extract peptide-protein pairs
    pairs = extract_peptide_protein_pairs(psms_table, precursors)
    
    # Convert to format expected by infer_proteins
    proteins_vec = Vector{NamedTuple{(:protein_name, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}}()
    peptides_vec = Vector{String}()
    
    for (pep_key, prot_key) in pairs
        push!(proteins_vec, (
            protein_name = prot_key.name,
            decoy = !prot_key.is_target,
            entrap_id = prot_key.entrap_id
        ))
        push!(peptides_vec, pep_key.sequence)
    end
    
    # Run core inference algorithm
    raw_inference = infer_proteins(proteins_vec, peptides_vec)
    
    # Convert results to our clean structure
    peptide_to_protein = Dictionary{PeptideKey, ProteinKey}()
    use_for_quant = Dictionary{PeptideKey, Bool}()
    
    for (pep_nt, prot_result) in pairs(raw_inference)
        pep_key = to_peptide_key(pep_nt)
        prot_key = ProteinKey(
            prot_result.protein_name,
            !prot_result.decoy,
            prot_result.entrap_id
        )
        
        insert!(peptide_to_protein, pep_key, prot_key)
        insert!(use_for_quant, pep_key, prot_result.retain)
    end
    
    return InferenceResult(peptide_to_protein, use_for_quant)
end

"""
    update_psms_with_inference(psms_df::DataFrame, inference::InferenceResult, precursors) -> DataFrame

Update PSM DataFrame with protein inference results.

# Arguments
- `psms_df`: DataFrame of PSMs
- `inference`: Protein inference results
- `precursors`: Library precursor information

# Returns
- Updated DataFrame with inference columns added
"""
function update_psms_with_inference(psms_df::DataFrame, inference::InferenceResult, precursors)
    n_psms = nrow(psms_df)
    use_for_protein_quant = Vector{Bool}(undef, n_psms)
    inferred_protein_names = Vector{String}(undef, n_psms)
    
    sequences = getSequence(precursors)
    is_decoys = getIsDecoy(precursors)
    entrap_ids = getEntrapmentGroupId(precursors)
    
    for i in 1:n_psms
        pid = psms_df[i, :precursor_idx]
        pep_key = PeptideKey(
            sequences[pid],
            !is_decoys[pid],
            entrap_ids[pid]
        )
        
        if haskey(inference.peptide_to_protein, pep_key)
            protein_key = inference.peptide_to_protein[pep_key]
            use_for_protein_quant[i] = get(inference.use_for_quant, pep_key, false)
            inferred_protein_names[i] = protein_key.name
        else
            # This shouldn't happen, but handle gracefully
            use_for_protein_quant[i] = false
            inferred_protein_names[i] = "UNKNOWN"
            @warn "Peptide not found in inference results" peptide=pep_key
        end
    end
    
    psms_df[!, :use_for_protein_quant] = use_for_protein_quant
    psms_df[!, :inferred_protein_group] = inferred_protein_names
    
    return psms_df
end

"""
    create_protein_groups_from_psms(psms_df::DataFrame, inference::InferenceResult, 
                                   precursors, min_peptides::Int) -> Dictionary{ProteinKey, ProteinGroupBuilder}

Create protein groups from PSMs using inference results.

# Arguments
- `psms_df`: DataFrame with PSMs and inference columns
- `inference`: Protein inference results
- `precursors`: Library precursor information
- `min_peptides`: Minimum peptides required for a group

# Returns
- Dictionary of protein group builders
"""
function create_protein_groups_from_psms(psms_df::DataFrame, inference::InferenceResult, 
                                       precursors, min_peptides::Int)
    groups = Dictionary{ProteinKey, ProteinGroupBuilder}()
    
    sequences = getSequence(precursors)
    is_decoys = getIsDecoy(precursors)
    entrap_ids = getEntrapmentGroupId(precursors)
    
    for i in 1:nrow(psms_df)
        # Skip if not used for quantification
        if !psms_df[i, :use_for_protein_quant]
            continue
        end
        
        pid = psms_df[i, :precursor_idx]
        score = psms_df[i, :prob]
        
        # Get protein assignment
        pep_key = PeptideKey(
            sequences[pid],
            !is_decoys[pid],
            entrap_ids[pid]
        )
        
        if !haskey(inference.peptide_to_protein, pep_key)
            continue
        end
        
        protein_key = inference.peptide_to_protein[pep_key]
        
        # Get or create builder
        if !haskey(groups, protein_key)
            insert!(groups, protein_key, ProteinGroupBuilder(protein_key))
        end
        
        # Add peptide to group
        add_peptide!(groups[protein_key], sequences[pid], score)
    end
    
    # Filter by minimum peptides
    filter!(x -> length(x.peptides) >= min_peptides, groups)
    
    return groups
end

"""
    calculate_protein_features(builder::ProteinGroupBuilder, catalog::Dictionary{ProteinKey, Set{String}}) -> ProteinFeatures

Calculate statistical features for a protein group.

# Arguments
- `builder`: Protein group builder
- `catalog`: Catalog of all possible peptides per protein

# Returns
- `ProteinFeatures` with calculated statistics
"""
function calculate_protein_features(builder::ProteinGroupBuilder, catalog::Dictionary{ProteinKey, Set{String}})
    # Get possible peptides for this protein group
    possible_peptides = Set{String}()
    
    # Handle protein groups (semicolon-delimited names)
    protein_names = split(builder.key.name, ';')
    for protein_name in protein_names
        individual_key = ProteinKey(String(protein_name), builder.key.is_target, builder.key.entrap_id)
        if haskey(catalog, individual_key)
            union!(possible_peptides, catalog[individual_key])
        end
    end
    
    n_observed = UInt16(length(builder.peptides))
    n_possible = UInt16(max(length(possible_peptides), 1))
    coverage = Float32(n_observed / n_possible)
    total_length = UInt16(sum(length(pep) for pep in builder.peptides))
    
    # Calculate derived features for regression
    log_n_possible = log(n_possible + 1)
    log_binom = if n_observed <= n_possible && n_observed >= 0
        lgamma(n_possible + 1) - lgamma(n_observed + 1) - lgamma(n_possible - n_observed + 1)
    else
        0.0
    end
    
    return ProteinFeatures(
        n_observed,
        n_possible,
        coverage,
        total_length,
        Float32(log_n_possible),
        Float32(log_binom)
    )
end

"""
    process_single_file_protein_inference(psm_path::String, ms_file_idx::Int, 
                                        passing_pg_paths::Vector{String},
                                        catalog::Dictionary{ProteinKey, Set{String}},
                                        precursors::LibraryPrecursors,
                                        file_mappings::FileMapping,
                                        min_peptides::Int,
                                        protein_groups_folder::String)

Process protein inference for a single file.

# Arguments
- `psm_path`: Path to PSM Arrow file
- `ms_file_idx`: Index in the file list
- `passing_pg_paths`: Vector to store output paths
- `catalog`: Protein-peptide catalog
- `precursors`: Library precursor information
- `file_mappings`: File mapping tracker
- `min_peptides`: Minimum peptides for a group
- `protein_groups_folder`: Output folder

# Returns
- Number of protein groups created
"""
function process_single_file_protein_inference(psm_path::String, ms_file_idx::Int,
                                             passing_pg_paths::Vector{String},
                                             catalog::Dictionary{ProteinKey, Set{String}},
                                             precursors::LibraryPrecursors,
                                             file_mappings::FileMapping,
                                             min_peptides::Int,
                                             protein_groups_folder::String)
    # Skip non-arrow files
    if !endswith(psm_path, ".arrow")
        return 0
    end
    
    # Load PSMs
    psms_table = Arrow.Table(psm_path)
    @info "Processing file for protein inference" file=basename(psm_path) n_psms=length(psms_table[:precursor_idx])
    
    # Perform inference
    inference = perform_file_protein_inference(psms_table, precursors)
    
    # Convert to DataFrame and update with inference results
    psms_df = DataFrame(Tables.columntable(psms_table))
    psms_df = update_psms_with_inference(psms_df, inference, precursors)
    
    # Create protein groups
    group_builders = create_protein_groups_from_psms(psms_df, inference, precursors, min_peptides)
    
    # Calculate features and finalize groups
    protein_groups = Dictionary{ProteinKey, ProteinGroup}()
    for (key, builder) in pairs(group_builders)
        features = calculate_protein_features(builder, catalog)
        insert!(protein_groups, key, finalize(builder, features))
    end
    
    # Write updated PSMs
    writeArrow(psm_path, psms_df)
    
    # Write protein groups
    pg_path = joinpath(protein_groups_folder, basename(psm_path))
    write_protein_groups_arrow(protein_groups, pg_path)
    
    # Update mappings
    passing_pg_paths[ms_file_idx] = pg_path
    add_mapping!(file_mappings, psm_path, pg_path)
    
    return length(protein_groups)
end

"""
    write_protein_groups_arrow(protein_groups::Dictionary{ProteinKey, ProteinGroup}, output_path::String)

Write protein groups to Arrow file format.

# Arguments
- `protein_groups`: Dictionary of protein groups
- `output_path`: Output file path
"""
function write_protein_groups_arrow(protein_groups::Dictionary{ProteinKey, ProteinGroup}, output_path::String)
    # Extract data for Arrow table
    n_groups = length(protein_groups)
    
    protein_names = Vector{String}(undef, n_groups)
    targets = Vector{Bool}(undef, n_groups)
    entrap_ids = Vector{UInt8}(undef, n_groups)
    pg_scores = Vector{Float32}(undef, n_groups)
    n_peptides = Vector{UInt16}(undef, n_groups)
    total_peptide_lengths = Vector{UInt16}(undef, n_groups)
    n_possible_peptides = Vector{UInt16}(undef, n_groups)
    peptide_coverages = Vector{Float32}(undef, n_groups)
    log_n_possible_peptides = Vector{Float32}(undef, n_groups)
    log_binom_coeffs = Vector{Float32}(undef, n_groups)
    
    for (i, (key, group)) in enumerate(pairs(protein_groups))
        protein_names[i] = group.key.name
        targets[i] = group.key.is_target
        entrap_ids[i] = group.key.entrap_id
        pg_scores[i] = group.score
        n_peptides[i] = group.features.n_peptides
        total_peptide_lengths[i] = group.features.total_peptide_length
        n_possible_peptides[i] = group.features.n_possible_peptides
        peptide_coverages[i] = group.features.peptide_coverage
        log_n_possible_peptides[i] = group.features.log_n_possible_peptides
        log_binom_coeffs[i] = group.features.log_binom_coeff
    end
    
    # Create DataFrame and sort
    df = DataFrame(
        protein_name = protein_names,
        target = targets,
        entrap_id = entrap_ids,
        pg_score = pg_scores,
        n_peptides = n_peptides,
        total_peptide_length = total_peptide_lengths,
        n_possible_peptides = n_possible_peptides,
        peptide_coverage = peptide_coverages,
        log_n_possible_peptides = log_n_possible_peptides,
        log_binom_coeff = log_binom_coeffs
    )
    
    sort!(df, [:pg_score, :target], rev = [true, true])
    
    # Write to Arrow
    Arrow.write(output_path, df)
end