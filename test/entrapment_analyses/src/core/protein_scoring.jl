"""
Protein-level scoring functions for EFDR analysis.

This version uses protein names directly for pairing instead of entrap_pair_id.
The key insight is that proteins with the same name are inherently paired.
"""

using DataFrames
using Dictionaries

"""
    add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)

Adds an original target score column for protein-level EFDR analysis.
For target proteins (entrap_id==0): gets their own score
For entrapment proteins (entrap_id>0): gets the score of the target with the same protein name

# Arguments
- `protein_results::DataFrame`: Protein results with columns:
  - `protein`: Protein name/identifier
  - `entrap_id`: Entrapment group ID (0 = original)
  - `file_name` or `ms_file_idx`: File identifier
  - Score column specified by `score_col`
- `score_col::Symbol`: Name of the score column to process (default: :pg_score)

# Notes
- Creates new column named `{score_col}_original_target`
- Sets score to -1.0 if no target (entrap_id==0) exists for that protein in the same file
- Handles both per-file and global analyses based on file identifier
"""
function add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)
    # Check for required columns
    required_cols = [:protein, :entrap_id]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    
    if !hasproperty(protein_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    
    # Determine file identifier column
    file_col = if hasproperty(protein_results, :ms_file_idx)
        :ms_file_idx
    elseif hasproperty(protein_results, :file_name)
        :file_name
    else
        error("DataFrame must have either :ms_file_idx or :file_name column")
    end
    
    # Create the original_target column name
    original_target_col = Symbol(String(score_col) * "_original_target")
    
    # Build a nested dictionary mapping (file_identifier, protein_name) -> target_score
    # This handles both per-file and global analyses
    if file_col == :ms_file_idx
        protein_to_target = Dictionary{Tuple{Int, String, String}, Float32}()
    else
        protein_to_target = Dictionary{Tuple{String, String, String}, Float32}()
    end
    
    # First pass: collect all target scores (entrap_id == 0)
    for row in eachrow(protein_results)
        if row.entrap_id == 0 && !ismissing(row[score_col])
            key = (row[file_col], row.species, row.protein)
            if haskey(protein_to_target, key)
                error("Duplicate target protein found: protein '$(row.protein)' appears multiple times " *
                      "with entrap_id=0 in file '$(row[file_col])'. Each protein should appear " *
                      "at most once as a target (entrap_id=0) per file.")
            end
            insert!(protein_to_target, key, Float32(row[score_col]))
        end
    end
    
    # Second pass: assign original target scores
    original_target_scores = Float32[]
    
    for row in eachrow(protein_results)
        if !ismissing(row[score_col])
            if row.entrap_id == 0
                # Target gets its own score
                push!(original_target_scores, Float32(row[score_col]))
            else
                # Entrapment gets the target's score from same protein and file
                key = (row[file_col], row.species, row.protein)
                if haskey(protein_to_target, key)
                    @info "Found target protein for entrapment protein '$(row.protein)' " *
                        "in file '$(row[file_col])'. Using score: $(protein_to_target[key])"
                    push!(original_target_scores, protein_to_target[key])
                else
                    # No target found for this protein in this file
                    @info "No target protein found for entrapment protein '$(row.protein)' " *
                        "in file '$(row[file_col])'. Setting original target score to -1.0."
                    push!(original_target_scores, -1.0f0)
                end
            end
        else
            # Missing score
            push!(original_target_scores, -1.0f0)
        end
    end
    
    protein_results[!, original_target_col] = original_target_scores
    return nothing
end

"""
    add_original_target_protein_scores!(protein_results::DataFrame, score_cols::Vector{Symbol})

Convenience function to add original target scores for multiple score columns at once.
"""
function add_original_target_protein_scores!(protein_results::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_protein_scores!(protein_results; score_col=score_col)
    end
    return nothing
end

"""
    create_global_protein_results_df(protein_results::DataFrame; score_col::Symbol=:global_pg_score)

Create a dataframe for global protein score analysis by selecting the best scoring row for each protein.

# Arguments
- `protein_results::DataFrame`: Original protein results with multiple files
- `score_col::Symbol`: Column name to use for selecting best score (default: :global_pg_score)

# Returns
- DataFrame with one row per protein (best scoring across all files)
- All rows have file identifier set to indicate global analysis
  - If ms_file_idx exists: set to 0
  - If file_name exists: set to "global"

# Example
```julia
global_df = create_global_protein_results_df(protein_results; score_col=:global_pg_score)
```
"""
function create_global_protein_results_df(protein_results::DataFrame; score_col::Symbol=:global_pg_score)
    # Check required columns
    required_cols = [:species,:protein,:file_name,score_col]
    missing_cols = [col for col in required_cols if !hasproperty(protein_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    
    # Determine file identifier column
    file_col = if hasproperty(protein_results, :ms_file_idx)
        :ms_file_idx
    elseif hasproperty(protein_results, :file_name)
        :file_name
    else
        error("DataFrame must have either :ms_file_idx or :file_name column")
    end
    
    # Create a deep copy to avoid aliasing issues
    protein_results_copy = copy(protein_results)
    
    # Group by protein name
    grouped = groupby(protein_results_copy, [:species,:protein])
    
    # For each group, keep only the row with maximum score
    global_df = combine(grouped) do group
        # Handle missing scores by filtering them out first
        valid_rows = group[.!ismissing.(group[!, score_col]), :]
        
        if nrow(valid_rows) == 0
            # If all scores are missing, return empty dataframe with same schema
            return similar(group, 0)
        end
        
        # Find row with maximum score
        best_idx = argmax(valid_rows[!, score_col])
        return valid_rows[best_idx:best_idx, :]
    end
    
    # Set file identifier to indicate global analysis
    if nrow(global_df) > 0
        if file_col == :ms_file_idx
            global_df[!, file_col] .= 0
        else
            global_df[!, file_col] .= "global"
        end
    end
    
    return global_df
end