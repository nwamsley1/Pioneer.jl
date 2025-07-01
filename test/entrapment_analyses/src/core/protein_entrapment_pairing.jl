"""
Protein-level entrapment pairing functions for EFDR analysis.

Unlike precursor-level pairing which considers charge states and modifications,
protein-level pairing only considers the protein name itself.
"""

"""
    assign_protein_entrapment_pairs!(df::DataFrame)

Assign unique pair IDs to protein entrapment groups.
For each unique protein name:
- Pairs rows where entrap_id==0 (originals) with rows from other entrap_id groups
- Uses round-robin pairing to distribute original proteins across entrapment versions
- Adds 'entrap_pair_id' column to the dataframe

# Arguments
- `df::DataFrame`: Protein results dataframe with columns:
  - `protein`: Protein name/identifier
  - `entrap_id`: Entrapment group ID (0 = original, >0 = entrapment)
  - `file_name`: File identifier

# Notes
- Unlike precursor pairing, proteins pair by name alone (no charge/mod considerations)
- Each original protein (entrap_id=0) gets paired with one protein from each entrapment group
"""
function assign_protein_entrapment_pairs!(df::DataFrame)
    # Check required columns
    required_cols = [:protein, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(df, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    
    # Initialize pair_id column if it doesn't exist
    if !hasproperty(df, :entrap_pair_id)
        df[!, :entrap_pair_id] = Vector{Union{Missing, UInt32}}(missing, nrow(df))
    end
    
    # Counter for unique pair IDs
    pair_counter = UInt32(1)
    
    # Group by protein name only
    grouped = groupby(df, :protein)
    
    for (_, group_df) in pairs(grouped)
        # Separate by entrapment group
        group_0_indices = findall(group_df.entrap_id .== 0)
        other_indices = findall(group_df.entrap_id .!= 0)
        
        if isempty(group_0_indices) || isempty(other_indices)
            # Skip if no pairing possible
            continue
        end
        
        # Get unique entrapment groups (excluding 0)
        entrap_groups_dict = Dict{UInt8, Vector{Int}}()
        for idx in other_indices
            entrap_id = group_df.entrap_id[idx]
            if !haskey(entrap_groups_dict, entrap_id)
                entrap_groups_dict[entrap_id] = Int[]
            end
            push!(entrap_groups_dict[entrap_id], idx)
        end
        
        unique_entrap_groups = sort(collect(keys(entrap_groups_dict)))
        
        # For each row in group 0 (original proteins)
        for (idx_0_pos, idx_0) in enumerate(group_0_indices)
            # Assign pair_id to group 0 member
            group_df[idx_0, :entrap_pair_id] = pair_counter
            
            # Get one row from each non-zero entrapment group
            for entrap_id in unique_entrap_groups
                group_member_indices = entrap_groups_dict[entrap_id]
                
                # Use round-robin to distribute group 0 proteins across other groups
                member_pos = ((idx_0_pos - 1) % length(group_member_indices)) + 1
                member_idx = group_member_indices[member_pos]
                
                group_df[member_idx, :entrap_pair_id] = pair_counter
            end
            
            pair_counter += UInt32(1)
        end
    end
    
    return nothing
end

"""
    add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame)

Add entrapment pair IDs to protein results based on protein names.

This function is provided for API consistency but may not be necessary if
protein results already contain entrap_id information. It maps pair IDs
from a protein library to results based on protein names.

# Arguments
- `protein_results::DataFrame`: Protein results with at least a `protein` column
- `protein_library::DataFrame`: Protein library with `protein` and `entrap_pair_id` columns

# Notes
- This function may be unnecessary if protein results already have entrap_id
- Consider whether the protein library is needed for your use case
"""
function add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame)
    # Check for required columns
    if !hasproperty(protein_results, :protein)
        error("protein_results must have :protein column")
    end
    if !hasproperty(protein_library, :entrap_pair_id)
        error("protein_library must have :entrap_pair_id column. Run assign_protein_entrapment_pairs! first.")
    end
    if !hasproperty(protein_library, :protein)
        error("protein_library must have :protein column")
    end
    
    # Create a mapping from protein name to entrap_pair_id
    # Note: If a protein appears multiple times with different pair_ids,
    # we'll use the first occurrence
    protein_to_pair = Dict{String, Union{Missing, UInt32}}()
    for row in eachrow(protein_library)
        if !haskey(protein_to_pair, row.protein)
            protein_to_pair[row.protein] = row.entrap_pair_id
        end
    end
    
    # Map entrap_pair_ids to results
    protein_results[!, :entrap_pair_id] = [
        get(protein_to_pair, protein, missing) 
        for protein in protein_results.protein
    ]
    
    return nothing
end