"""
    get_complement_score(prec_results::DataFrame, row_idx::Int; score_col=:score)

For a specific row, returns the score of its complement (same entrap_pair_id, different row).
Returns 0.0 if no complement found. If multiple complements exist, returns the maximum score.
"""
function get_complement_score(prec_results::DataFrame, row_idx::Int; score_col=:score)
    # Check for required columns
    if !hasproperty(prec_results, :entrap_pair_id)
        error("DataFrame must have :entrap_pair_id column. Run add_entrap_pair_ids! first.")
    end
    if !hasproperty(prec_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    
    # Check row index
    if row_idx < 1 || row_idx > nrow(prec_results)
        error("row_idx $row_idx is out of bounds. DataFrame has $(nrow(prec_results)) rows.")
    end
    
    # Get the entrap_pair_id for this row
    pair_id = prec_results[row_idx, :entrap_pair_id]
    if ismissing(pair_id)
        return 0.0
    end
    
    # Find complement rows (same pair_id, different row)
    complement_mask = (prec_results.entrap_pair_id .== pair_id) .& (1:nrow(prec_results) .!= row_idx)
    complement_scores = prec_results[complement_mask, score_col]
    
    return isempty(complement_scores) ? 0.0 : maximum(skipmissing(complement_scores), init=0.0)
end

"""
    add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame; score_col=:score)

Adds an original target score column named as score_col with "_original_target" appended.
For target sequences (entrapment_group_id==0): gets their own score
For entrapment sequences (entrapment_group_id>0): gets the score of the target in their pair
"""
function add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame; score_col=:score)
    # Check for required columns
    print("test!")
    if !hasproperty(prec_results, :entrap_pair_id)
        error("DataFrame must have :entrap_pair_id column. Run add_entrap_pair_ids! first.")
    end
    if !hasproperty(prec_results, score_col)
        error("DataFrame must have :$score_col column.")
    end
    if !hasproperty(prec_results, :precursor_idx)
        error("DataFrame must have :precursor_idx column.")
    end
    if !hasproperty(library_precursors, :entrapment_group_id)
        error("library_precursors must have :entrapment_group_id column.")
    end
    
    # Create the original_target column name
    original_target_col = Symbol(String(score_col) * "_original_target")
    
    # Build a nested dictionary mapping ms_file_idx -> pair_id -> target info
    # This handles both per-file (ms_file_idx > 0) and global (ms_file_idx = 0) analyses
    pair_to_target = Dictionary{Int, Dictionary{UInt32, @NamedTuple{target_row::UInt32,target_score::Float32}}}()
    #@info "ms_file_idx values in prec_results: $(unique(prec_results.ms_file_idx))"
    #@info "prec_results pair id 942181 " prec_results[prec_results[!,:entrap_pair_id].==942181,:]
    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrap_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrap_pair_id
            precursor_idx = row.precursor_idx
            ms_file_idx = row.ms_file_idx 
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            
            # If this is an original target (group 0), store its score for the pair
            if entrap_group == 0
                # Initialize inner dictionary if it doesn't exist
                if !haskey(pair_to_target, ms_file_idx)
                    insert!(pair_to_target, ms_file_idx, Dictionary{UInt32, @NamedTuple{target_row::UInt32,target_score::Float32}}())
                end
                insert!(
                    pair_to_target[ms_file_idx], pair_id, 
                    (target_row = UInt32(idx), target_score = row[score_col])
                )
            end
        end
    end
    
    # Assign original target scores
    original_target_scores = zeros(Float32, nrow(prec_results))

    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrap_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrap_pair_id
            precursor_idx = row.precursor_idx
            ms_file_idx = row.ms_file_idx 
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            
            if entrap_group == 0
                # Target gets its own score
                original_target_scores[idx] = row[score_col]
            else
                # Entrapment gets the target's score from its pair
                if haskey(pair_to_target, ms_file_idx) && haskey(pair_to_target[ms_file_idx], pair_id)
                    # Original target complement found in the data
                    _, target_score = pair_to_target[ms_file_idx][pair_id]
                    original_target_scores[idx] = target_score
                else
                    # Complement target wasn't observed in the data so its score is set to -1.0
                    original_target_scores[idx] = -1.0
                end
            end
        end
    end
    prec_results[!, original_target_col] = original_target_scores
    return nothing
end

"""
    add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame, score_cols::Vector{Symbol})

Convenience function to add original target scores for multiple score columns at once.
"""
function add_original_target_scores!(prec_results::DataFrame, library_precursors::DataFrame, score_cols::Vector{Symbol})
    for score_col in score_cols
        add_original_target_scores!(prec_results, library_precursors; score_col=score_col)
    end
    return nothing
end

"""
    create_global_results_df(prec_results::DataFrame; score_col::Symbol=:global_prob)

Create a dataframe for global score analysis by selecting the best scoring row for each precursor.

# Arguments
- `prec_results::DataFrame`: Original precursor results with multiple files
- `score_col::Symbol`: Column name to use for selecting best score (default: :global_prob)

# Returns
- DataFrame with one row per precursor (best scoring across all files)
- All rows have ms_file_idx set to 0 to indicate global analysis

# Example
```julia
global_df = create_global_results_df(prec_results; score_col=:global_prob)
```
"""
function create_global_results_df(prec_results::DataFrame; score_col::Symbol=:global_prob)
    # Check required columns
    required_cols = [:precursor_idx, :ms_file_idx, score_col]
    missing_cols = [col for col in required_cols if !hasproperty(prec_results, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    
    # Create a deep copy to avoid aliasing issues
    prec_results_copy = copy(prec_results)
    
    # Group by precursor_idx
    grouped = groupby(prec_results_copy, :precursor_idx)
    
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
    
    # Set all ms_file_idx to 0 to indicate global analysis
    if nrow(global_df) > 0
        global_df[!, :ms_file_idx] .= 0
    end
    
    return global_df
end