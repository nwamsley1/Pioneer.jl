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
    
    # Build a dictionary mapping pair_id to target row indices and scores
    pair_to_target = Dict{UInt32, Tuple{Int,Float64}}()
    
    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrap_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrap_pair_id
            precursor_idx = row.precursor_idx
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            
            # If this is a target (group 0), store its score for the pair
            if entrap_group == 0
                pair_to_target[pair_id] = (idx, Float64(row[score_col]))
            end
        end
    end
    
    # Assign original target scores
    original_target_scores = zeros(Float64, nrow(prec_results))
    
    for (idx, row) in enumerate(eachrow(prec_results))
        if !ismissing(row.entrap_pair_id) && !ismissing(row[score_col])
            pair_id = row.entrap_pair_id
            precursor_idx = row.precursor_idx
            entrap_group = library_precursors.entrapment_group_id[precursor_idx]
            
            if entrap_group == 0
                # Target gets its own score
                original_target_scores[idx] = Float64(row[score_col])
            else
                # Entrapment gets the target's score from its pair
                if haskey(pair_to_target, pair_id)
                    _, target_score = pair_to_target[pair_id]
                    original_target_scores[idx] = target_score
                else
                    # No target found in pair, default to 0.0
                    original_target_scores[idx] = 0.0
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