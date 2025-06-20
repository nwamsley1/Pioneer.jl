"""
    protein_inference_pipeline.jl

Composable pipeline operations for protein inference workflow.
This module provides a clear, testable interface for protein inference
while preserving the core `infer_proteins` algorithm.
"""

using DataFrames
using Arrow
using Tables
using Dictionaries

# Import the core inference algorithm
using Pioneer: infer_proteins

#==========================================================
Pipeline Operations for Protein Inference
==========================================================#

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add peptide sequence and protein information from precursor library.
Enriches PSM data with metadata needed for inference.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"
    
    op = function(df)
        # Get library data
        all_sequences = getSequence(precursors)
        all_accessions = getAccessionNumbers(precursors)
        all_is_decoys = getIsDecoy(precursors)
        all_entrap_ids = getEntrapmentGroupId(precursors)
        
        # Add columns if they don't exist
        n_rows = nrow(df)
        if !hasproperty(df, :sequence)
            sequences = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                pid = df.precursor_idx[i]
                sequences[i] = all_sequences[pid]
            end
            df.sequence = sequences
        end
        
        if !hasproperty(df, :accession_numbers)
            accessions = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                pid = df.precursor_idx[i]
                accessions[i] = all_accessions[pid]
            end
            df.accession_numbers = accessions
        end
        
        # Add target/decoy status if needed
        if !hasproperty(df, :is_decoy)
            is_decoy_vec = Vector{Bool}(undef, n_rows)
            for i in 1:n_rows
                pid = df.precursor_idx[i]
                is_decoy_vec[i] = all_is_decoys[pid]
            end
            df.is_decoy = is_decoy_vec
        end
        
        # Add entrapment group if needed
        if !hasproperty(df, :entrapment_group_id)
            entrap_vec = Vector{UInt8}(undef, n_rows)
            for i in 1:n_rows
                pid = df.precursor_idx[i]
                entrap_vec[i] = all_entrap_ids[pid]
            end
            df.entrapment_group_id = entrap_vec
        end
        
        return df
    end
    
    return desc => op
end

"""
    deduplicate_peptides()

Remove duplicate peptide-protein pairs to prepare for inference.
Keeps the highest scoring PSM for each unique peptide.
"""
function deduplicate_peptides()
    desc = "deduplicate_peptides"
    
    op = function(df)
        # Group by peptide sequence and protein, keep best scoring
        if nrow(df) == 0
            return df
        end
        
        # Create composite key for grouping
        gdf = groupby(df, [:sequence, :accession_numbers, :entrapment_group_id])
        
        # Keep row with highest probability for each group
        deduped = combine(gdf) do sdf
            if nrow(sdf) == 1
                return sdf
            else
                # Find row with max probability
                max_idx = argmax(sdf.prec_prob)
                return sdf[max_idx:max_idx, :]
            end
        end
        
        return deduped
    end
    
    return desc => op
end

"""
    validate_peptide_data()

Ensure all required columns exist for protein inference.
"""
function validate_peptide_data()
    desc = "validate_peptide_data"
    
    op = function(df)
        required_cols = [:sequence, :accession_numbers, :is_decoy, 
                        :entrapment_group_id, :precursor_idx]
        
        missing_cols = setdiff(required_cols, Symbol.(names(df)))
        if !isempty(missing_cols)
            error("Missing required columns for protein inference: $missing_cols")
        end
        
        return df
    end
    
    return desc => op
end

#==========================================================
Inference Wrapper Functions
==========================================================#

"""
    InferenceResult

Container for protein inference results.
"""
struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
    use_for_quant::Dictionary{PeptideKey, Bool}
end

"""
    apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)

Apply the core infer_proteins algorithm to a prepared DataFrame.
Returns structured inference results.
"""
function apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)
    if nrow(df) == 0
        # Return empty result for empty input
        return InferenceResult(
            Dictionary{PeptideKey, ProteinKey}(),
            Dictionary{PeptideKey, Bool}()
        )
    end
    
    # Extract unique peptide-protein pairs
    unique_pairs = unique(df, [:sequence, :accession_numbers, :is_decoy, :entrapment_group_id])
    
    # Convert to format expected by infer_proteins
    proteins_vec = Vector{NamedTuple{(:protein_name, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}}()
    peptides_vec = Vector{String}()
    
    for row in eachrow(unique_pairs)
        push!(proteins_vec, (
            protein_name = row.accession_numbers,
            decoy = row.is_decoy,
            entrap_id = row.entrapment_group_id
        ))
        push!(peptides_vec, row.sequence)
    end
    
    # Call core inference algorithm
    raw_result = infer_proteins(proteins_vec, peptides_vec)
    
    # Convert to structured result
    peptide_to_protein = Dictionary{PeptideKey, ProteinKey}()
    use_for_quant = Dictionary{PeptideKey, Bool}()
    
    for (pep_nt, prot_result) in pairs(raw_result)
        pep_key = PeptideKey(pep_nt.peptide, !pep_nt.decoy, pep_nt.entrap_id)
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

#==========================================================
Post-Inference Pipeline Operations
==========================================================#

"""
    add_inferred_protein_column(inference_result::InferenceResult)

Add inferred protein group assignments to PSMs based on inference results.
"""
function add_inferred_protein_column(inference_result::InferenceResult)
    desc = "add_inferred_protein_column"
    
    op = function(df)
        n_rows = nrow(df)
        inferred_proteins = Vector{String}(undef, n_rows)
        
        for i in 1:n_rows
            pep_key = PeptideKey(
                df.sequence[i],
                !df.is_decoy[i],
                df.entrapment_group_id[i]
            )
            
            if haskey(inference_result.peptide_to_protein, pep_key)
                protein_key = inference_result.peptide_to_protein[pep_key]
                inferred_proteins[i] = protein_key.name
            else
                # Fallback to original protein
                inferred_proteins[i] = df.accession_numbers[i]
            end
        end
        
        df.inferred_protein_group = inferred_proteins
        return df
    end
    
    return desc => op
end

"""
    add_quantification_flag(inference_result::InferenceResult)

Add flag indicating whether peptide should be used for quantification.
"""
function add_quantification_flag(inference_result::InferenceResult)
    desc = "add_quantification_flag"
    
    op = function(df)
        n_rows = nrow(df)
        use_for_quant = Vector{Bool}(undef, n_rows)
        
        for i in 1:n_rows
            pep_key = PeptideKey(
                df.sequence[i],
                !df.is_decoy[i],
                df.entrapment_group_id[i]
            )
            
            if haskey(inference_result.use_for_quant, pep_key)
                use_for_quant[i] = inference_result.use_for_quant[pep_key]
            else
                # Default to true if not in inference results
                use_for_quant[i] = true
            end
        end
        
        df.use_for_protein_quant = use_for_quant
        return df
    end
    
    return desc => op
end

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
function group_psms_by_protein(df::DataFrame)
    if nrow(df) == 0
        # Return empty protein groups DataFrame with expected schema
        return DataFrame(
            protein_name = String[],
            target = Bool[],
            entrapment_group_id = UInt8[],
            n_peptides = Int64[],
            peptide_list = Vector{String}[],
            pg_score = Float32[]
        )
    end
    
    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrapment_group_id])
    
    # Aggregate to protein groups
    protein_groups = combine(grouped) do gdf
        # Get unique peptides that are used for quantification
        quant_peptides = unique(gdf[gdf.use_for_protein_quant .== true, :sequence])
        n_peptides = length(quant_peptides)
        
        # Calculate initial protein score (log-sum)
        peptide_probs = gdf[gdf.use_for_protein_quant .== true, :prec_prob]
        if isempty(peptide_probs)
            pg_score = 0.0f0
        else
            # Use best probability per peptide
            unique_pep_probs = Float32[]
            for pep in quant_peptides
                pep_mask = (gdf.sequence .== pep) .& (gdf.use_for_protein_quant .== true)
                if any(pep_mask)
                    push!(unique_pep_probs, maximum(gdf[pep_mask, :prec_prob]))
                end
            end
            pg_score = -sum(log.(1.0f0 .- unique_pep_probs))
        end
        
        DataFrame(
            n_peptides = n_peptides,
            peptide_list = [quant_peptides],
            pg_score = pg_score
        )
    end
    
    # Rename the grouping column
    rename!(protein_groups, :inferred_protein_group => :protein_name)
    
    return protein_groups
end

"""
    filter_by_min_peptides(min_peptides::Int)

Filter protein groups that don't meet minimum peptide requirement.
"""
function filter_by_min_peptides(min_peptides::Int)
    desc = "filter_by_min_peptides(min=$min_peptides)"
    
    op = function(df)
        if hasproperty(df, :n_peptides)
            # For protein groups
            filter!(row -> row.n_peptides >= min_peptides, df)
        end
        return df
    end
    
    return desc => op
end

"""
    add_protein_features(protein_catalog::Dict)

Add protein-level features like peptide coverage.
"""
function add_protein_features(protein_catalog::Dict)
    desc = "add_protein_features"
    
    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end
        
        n_rows = nrow(df)
        n_possible = Vector{Int64}(undef, n_rows)
        peptide_coverage = Vector{Float32}(undef, n_rows)
        
        for i in 1:n_rows
            key = (
                protein_name = df.protein_name[i],
                target = df.target[i],
                entrap_id = df.entrapment_group_id[i]
            )
            
            if haskey(protein_catalog, key)
                possible_peptides = protein_catalog[key]
                n_possible[i] = length(possible_peptides)
                
                # Calculate coverage
                if n_possible[i] > 0
                    peptide_coverage[i] = Float32(df.n_peptides[i]) / Float32(n_possible[i])
                else
                    peptide_coverage[i] = 0.0f0
                end
            else
                n_possible[i] = 0
                peptide_coverage[i] = 0.0f0
            end
        end
        
        df.n_possible_peptides = n_possible
        df.peptide_coverage = peptide_coverage
        
        return df
    end
    
    return desc => op
end

#==========================================================
High-Level Interface
==========================================================#

"""
    perform_protein_inference_pipeline(psm_refs, output_folder, precursors, protein_catalog; kwargs...)

Perform protein inference using a composable pipeline approach.

# Arguments
- `psm_refs`: Vector of PSM file references
- `output_folder`: Directory for protein group output
- `precursors`: Library precursor information
- `protein_catalog`: Pre-computed protein-to-peptide mappings
- `min_peptides`: Minimum peptides per protein group (default: 2)

# Returns
- `pg_refs`: Vector of protein group file references
- `psm_to_pg_mapping`: Dictionary mapping PSM paths to protein group paths
"""
function perform_protein_inference_pipeline(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    precursors::LibraryPrecursors,
    protein_catalog::Dict;
    min_peptides::Int = 2
)
    # Ensure output folder exists
    !isdir(output_folder) && mkpath(output_folder)
    
    # Build pre-inference pipeline
    pre_inference_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors) |>
        deduplicate_peptides() |>
        validate_peptide_data()
    
    # Process each file
    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    
    @info "Performing protein inference on $(length(psm_refs)) files"
    
    for (idx, psm_ref) in enumerate(psm_refs)
        if !exists(psm_ref)
            continue
        end
        
        # Step 1: Apply pre-inference pipeline to add necessary columns
        apply_pipeline!(psm_ref, pre_inference_pipeline)
        
        # Step 2: Run inference on the prepared PSMs
        prepared_df = load_dataframe(psm_ref)
        inference_result = apply_inference_to_dataframe(prepared_df, precursors)
        
        # Step 3: Update PSMs with inference results
        psm_update_pipeline = TransformPipeline() |>
            add_inferred_protein_column(inference_result) |>
            add_quantification_flag(inference_result)
        
        # Apply updates to the same PSM file (which now has necessary columns)
        apply_pipeline!(psm_ref, psm_update_pipeline)
        
        # Step 4: Create protein groups
        # Reload updated PSMs
        updated_psms = load_dataframe(psm_ref)
        
        # Group by protein
        protein_groups_df = group_psms_by_protein(updated_psms)
        
        # Build post-inference pipeline
        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_protein_features(protein_catalog)
        
        # Apply post-processing
        for (desc, op) in post_inference_pipeline.operations
            protein_groups_df = op(protein_groups_df)
        end
        
        # Write protein groups
        pg_filename = "protein_groups_$(lpad(idx, 3, '0')).arrow"
        pg_path = joinpath(output_folder, pg_filename)
        
        if nrow(protein_groups_df) > 0
            writeArrow(pg_path, protein_groups_df)
            pg_ref = ProteinGroupFileReference(pg_path)
            push!(pg_refs, pg_ref)
            psm_to_pg_mapping[file_path(psm_ref)] = pg_path
        end
    end
    
    @info "Created $(length(pg_refs)) protein group files"
    
    return pg_refs, psm_to_pg_mapping
end

# Export functions
export perform_protein_inference_pipeline,
       add_peptide_metadata,
       deduplicate_peptides,
       validate_peptide_data,
       add_inferred_protein_column,
       add_quantification_flag,
       filter_by_min_peptides,
       add_protein_features,
       apply_inference_to_dataframe,
       group_psms_by_protein,
       InferenceResult