"""
    getModKey(mod_string::AbstractString)

Extract modification names from a modification string, sort them, and concatenate with semicolons.

The input string contains modifications in the format (position,amino_acid,mod_name).
This function extracts all mod_names, sorts them, and returns them joined by semicolons.

# Examples
```julia
getModKey("(5,M,x)(5,M,Unimod:4)(5,M,Unimod:1)") # => "Unimod:1;Unimod:4;x"
getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") # => "Unimod:4;Unimod:35"
getModKey("(5,M,x)") # => "x"
```
"""
function getModKey(mod_string::AbstractString)
    # Regular expression to match the pattern (number,letter,mod_name)
    # Captures the mod_name part (everything after the second comma until the closing parenthesis)
    mod_pattern = r"\(\d+,[A-Z],([^)]+)\)"
    
    # Extract all modification names
    mod_names = String[]
    for match in eachmatch(mod_pattern, mod_string)
        push!(mod_names, match.captures[1])
    end
    
    # Sort the modification names
    sort!(mod_names)
    
    # Join with semicolons
    return join(mod_names, ";")
end

"""
    assign_entrapment_pairs!(df::DataFrame)

Assign unique pair IDs to entrapment groups within each base peptide group.
For each group (defined by base_pep_id, prec_charge, is_decoy, mod_key):
- Pairs sequences where entrapment_group_id==0 with sequences from other groups
- If multiple entrapment groups exist (>0), each group 0 sequence is paired with one from each other group
- Returns the dataframe with added 'entrap_pair_id' column
"""
function assign_entrapment_pairs!(df::DataFrame)
    # Initialize pair_id column if it doesn't exist
    if !hasproperty(df, :entrap_pair_id)
        df[!, :entrap_pair_id] = Vector{Union{Missing, UInt32}}(missing, nrow(df))
    end
    
    # Counter for unique pair IDs
    pair_counter = UInt32(1)
    
    # Group by base peptide characteristics
    grouped = groupby(df, [:base_pep_id, :prec_charge, :is_decoy, :mod_key])
    
    for (key, group_df) in pairs(grouped)
        # Separate by entrapment group
        group_0_indices = findall(group_df.entrapment_group_id .== 0)
        other_indices = findall(group_df.entrapment_group_id .!= 0)
        
        if isempty(group_0_indices) || isempty(other_indices)
            # Skip if no pairing possible
            continue
        end
        
        # Get unique entrapment groups (excluding 0)
        entrap_groups_dict = Dict{Int, Vector{Int}}()
        for idx in other_indices
            entrap_id = group_df.entrapment_group_id[idx]
            if !haskey(entrap_groups_dict, entrap_id)
                entrap_groups_dict[entrap_id] = Int[]
            end
            push!(entrap_groups_dict[entrap_id], idx)
        end
        
        unique_entrap_groups = sort(collect(keys(entrap_groups_dict)))
        
        # For each sequence in group 0
        for (idx_0_pos, idx_0) in enumerate(group_0_indices)
            # Assign pair_id to group 0 member
            group_df[idx_0, :entrap_pair_id] = pair_counter
            
            # Get one sequence from each non-zero entrapment group
            for entrap_id in unique_entrap_groups
                group_member_indices = entrap_groups_dict[entrap_id]
                
                # Use round-robin to distribute group 0 sequences across other groups
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
    add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)

Add entrapment pair IDs to precursor results based on precursor indices.

Maps the entrap_pair_id from library_precursors to prec_results using the precursor_idx column.
"""
function add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)
    # Check for required columns
    if !hasproperty(prec_results, :precursor_idx)
        error("prec_results must have :precursor_idx column")
    end
    if !hasproperty(library_precursors, :entrap_pair_id)
        error("library_precursors must have :entrap_pair_id column")
    end
    
    # Map entrap_pair_ids from library to results
    prec_results[!, :entrap_pair_id] = [library_precursors[!, :entrap_pair_id][pid] for pid in prec_results[!, :precursor_idx]]
    
    return nothing
end