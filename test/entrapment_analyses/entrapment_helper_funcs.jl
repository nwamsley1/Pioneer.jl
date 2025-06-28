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

# Test the function
@test getModKey("(5,M,x)(5,M,Unimod:4)(5,M,Unimod:1)") == "Unimod:1;Unimod:4;x"
@test getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") == "Unimod:35;Unimod:4"
@test getModKey("(5,M,Unimod:35)(5,M,Unimod:4)") == "Unimod:35;Unimod:4"
@test getModKey("(5,M,x)") == "x"

"""
    assign_entrapment_pairs!(df::DataFrame)

Assign unique pair IDs to entrapment groups within each base peptide group.
For each group (defined by base_pep_id, prec_charge, is_decoy, mod_key):
- Pairs sequences where entrapment_group_id==0 with sequences from other groups
- If multiple entrapment groups exist (>0), each group 0 sequence is paired with one from each other group
- Returns the dataframe with added 'pair_id' column
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
            # Get the actual row index in the original dataframe
            #row_idx_0 = rownumber(group_df[idx_0, :])
            
            # Assign pair_id to group 0 member
            group_df[idx_0, :entrap_pair_id] = pair_counter
            
            # Get one sequence from each non-zero entrapment group
            for entrap_id in unique_entrap_groups
                group_member_indices = entrap_groups_dict[entrap_id]
                
                # Use round-robin to distribute group 0 sequences across other groups
                member_pos = ((idx_0_pos - 1) % length(group_member_indices)) + 1
                member_idx = group_member_indices[member_pos]
                
                # Get the actual row index in the original dataframe
                #row_idx = rownumber(group[member_idx, :])
                group_df[member_idx, :entrap_pair_id] = pair_counter
            end
            
            pair_counter += UInt32(1)
        end
    end
    
    return nothing
end

# Tests for assign_entrapment_pairs!
@testset "assign_entrapment_pairs! tests" begin
    # Test 1: Basic pairing with one target and one entrapment group
    # Expected: 1 pair with 2 members (group 0 + one from group 1)
    test_df1 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [0, 1, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df1)
    
    @test !any(ismissing.(test_df1.entrap_pair_id[1:2]))  # First two should be paired
    @test test_df1.entrap_pair_id[1] == test_df1.entrap_pair_id[2]  # Same pair ID
    @test ismissing(test_df1.entrap_pair_id[3])  # Extra entrapment seq not considered (this shouldn't happen in real data)
    
    # Test 2: Multiple entrapment groups
    # Expected: 2 pairs (one for each group 0), each with 3 members (group 0 + one from group 1 + one from group 2)
    test_df2 = DataFrame(
        base_pep_id = fill(1, 6),
        prec_charge = fill(2, 6),
        is_decoy = fill(false, 6),
        mod_key = fill("", 6),
        entrapment_group_id = [0, 0, 1, 1, 2, 2],
        sequence = fill("PEPTIDE", 6)
    )
    assign_entrapment_pairs!(test_df2)
    
    # Each group 0 should be paired with one from group 1 and one from group 2
    @test length(unique(test_df2.entrap_pair_id[.!ismissing.(test_df2.entrap_pair_id)])) == 2  # Two unique pairs
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[1]) == 3  # First pair has 3 members
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[2]) == 3  # Second pair has 3 members
    
    # Test 3: No valid pairs (no group 0)
    # Expected: No pairs created
    test_df3 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [1, 1, 2],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df3)
    
    @test all(ismissing.(test_df3.entrap_pair_id))  # No pairs should be created
    
    # Test 4: Different base peptides should not be paired
    # Expected: 2 pairs (one for each group), each with 2 members
    test_df4 = DataFrame(
        base_pep_id = [1, 1, 2, 2],
        prec_charge = [2, 2, 2, 2],
        is_decoy = [false, false, false, false],
        mod_key = ["", "", "", ""],
        entrapment_group_id = [0, 1, 0, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PROTEIN", "PROTEIN"]
    )
    assign_entrapment_pairs!(test_df4)
    
    @test length(unique(skipmissing(test_df4.entrap_pair_id))) == 2  # Two separate pairs
    @test !ismissing(test_df4.entrap_pair_id[1]) && !ismissing(test_df4.entrap_pair_id[2])  # First group paired
    @test !ismissing(test_df4.entrap_pair_id[3]) && !ismissing(test_df4.entrap_pair_id[4])  # Second group paired
    @test test_df4.entrap_pair_id[1] == test_df4.entrap_pair_id[2]  # First peptide pair
    @test test_df4.entrap_pair_id[3] == test_df4.entrap_pair_id[4]  # Second peptide pair
    @test test_df4.entrap_pair_id[1] != test_df4.entrap_pair_id[3]  # Different pairs
    
    # Test 5: Round-robin distribution
    # Expected: 3 pairs (one for each group 0), each with 2 members
    # The group 1 members are distributed round-robin style
    test_df5 = DataFrame(
        base_pep_id = fill(1, 5),
        prec_charge = fill(2, 5),
        is_decoy = fill(false, 5),
        mod_key = fill("", 5),
        entrapment_group_id = [0, 0, 0, 1, 1],
        sequence = fill("PEPTIDE", 5)
    )
    assign_entrapment_pairs!(test_df5)
    
    # Three group 0 entries, two group 1 entries
    # Each group 0 gets paired with one group 1 (round-robin)
    # First group 0 pairs with first group 1
    # Second group 0 pairs with second group 1  
    # Third group 0 pairs with first group 1 again (round-robin)
    @test length(unique(skipmissing(test_df5.entrap_pair_id))) == 3  # Three unique pairs
    @test all(test_df5.entrap_pair_id .== UInt32[1, 2, 3, 3, 2])
    
    println("All assign_entrapment_pairs! tests passed!")
end

"""
    add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)

Add entrap_pair_id column to prec_results by mapping precursor_idx to the library.
"""
function add_entrap_pair_ids!(prec_results::DataFrame, library_precursors::DataFrame)
    # Check for required columns
    if !hasproperty(prec_results, :precursor_idx)
        error("prec_results must have :precursor_idx column.")
    end
    if !hasproperty(library_precursors, :entrap_pair_id)
        error("library_precursors must have :entrap_pair_id column.")
    end
    
    # Check that all precursor_idx values are valid
    max_idx = maximum(prec_results.precursor_idx)
    if max_idx > nrow(library_precursors)
        error("precursor_idx values exceed library_precursors rows. Max idx: $max_idx, library rows: $(nrow(library_precursors))")
    end
    
    prec_results[!,:entrap_pair_id] = [library_precursors.entrap_pair_id[pid] for pid in prec_results[!,:precursor_idx]]
    return nothing
end

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

