# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#=
function reverseSequence(sequence::AbstractString, structural_mods::String)
    # Early return if no modifications
    if isempty(structural_mods)
        return reverse(sequence[1:end-1])*sequence[end], ""
    end
    
    # Pre-count number of modifications for array pre-allocation
    n_mods = count("(", structural_mods)
    
    # Pre-allocate arrays
    mod_positions = Vector{Int}(undef, n_mods)
    mod_aas = Vector{Char}(undef, n_mods)
    mod_names = Vector{String}(undef, n_mods)
    
    # Parse modifications into parallel arrays
    mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    i = 1
    for m in eachmatch(mod_regex, structural_mods)
        mod_positions[i] = parse(Int, m.captures[1])
        mod_aas[i] = first(m.captures[2])
        mod_names[i] = m.captures[3]
        i += 1
    end
    
    # Create reversed sequence (keeping last AA in place)
    seq_length = length(sequence)
    rev_seq = reverse(sequence[1:end-1])*sequence[end]
    
    # Calculate new positions and create tuples for sorting
    n_final_mods = n_mods
    final_mods = Vector{Tuple{Int, String}}(undef, n_final_mods)
    j = 1
    for i in 1:n_mods
        pos = mod_positions[i]
        aa = mod_aas[i]
        name = mod_names[i]
        
        if aa == 'n'
            final_mods[j] = (1, "(1,n,$name)")
        elseif pos == seq_length
            final_mods[j] = (seq_length, "($pos,$aa,$name)")
        else
            new_pos = seq_length - pos
            new_aa = rev_seq[new_pos]
            final_mods[j] = (new_pos, "($new_pos,$new_aa,$name)")
        end
        j += 1
    end
    
    # Sort by new positions and join
    sort!(final_mods, by=first)
    return rev_seq, join(last.(final_mods))
end

function reverseSequence(sequence::AbstractString, structural_mods::Missing)
    return reverse(sequence[1:end-1])*sequence[end], missing
end

function parse_mods_to_set(structural_mods::String)
    # Return an empty set if no modifications
    if isempty(structural_mods)
        return Set{Tuple{Int, Char, String}}()
    end
    
    # Parse modifications into a set of tuples
    mod_set = Set{Tuple{Int, Char, String}}()
    mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    
    for m in eachmatch(mod_regex, structural_mods)
        position = parse(Int, m.captures[1])
        aa = first(m.captures[2])
        mod_name = m.captures[3]
        push!(mod_set, (position, aa, mod_name))
    end
    
    return mod_set
end

function add_precursor_partner_columns!(df::DataFrame)
    """
    Adds both 'partner_precursor_idx' and 'pair_id' columns in a single pass.
    'partner_precursor_idx' links each row to its target/decoy partner.
    'pair_id' provides a unique UInt32 identifier for each target/decoy pair.
    """
    # Make a copy to avoid modifying the original
    result_df = copy(df)

    # Create a lookup dictionary for each sequence+attributes combination
    lookup = Dict()
    for (idx, row) in enumerate(eachrow(df))
        if !ismissing(row.structural_mods)
            mods = parse_mods_to_set(row.structural_mods)
        else
            mods = missing
        end
        key = (row.proteome_identifiers, row.sequence, mods, row.prec_charge)
        lookup[key] = UInt32(idx)  # Store the row index as precursor_idx
    end
    
    # Create vectors for partner indices and pair IDs
    n = nrow(df)
    partner_indices = Vector{Union{UInt32, Missing}}(missing, n)
    pair_ids = Vector{Union{UInt32, Missing}}(missing, n)
    
    # Track which pairs have been processed
    pair_to_pair_id = Dict{Tuple{UInt32, UInt32}, UInt32}()
    next_pair_id = UInt32(0)
    missed_rows = 1
    # Process each row to find partners
    for idx in 1:n
        # Get current row data
        row = df[idx, :]
        if row.is_decoy == true
            continue
        end
        # Calculate what the partner sequence would be
        partner_seq , partner_mods = reverseSequence(row.sequence, row.structural_mods)
        if !ismissing(partner_mods)
            partner_mods = parse_mods_to_set(partner_mods)
        end
        # Look up potential partner with opposite is_decoy status
        partner_key = (row.proteome_identifiers, partner_seq, partner_mods, row.prec_charge)
        if haskey(lookup, partner_key)
            partner_idx = lookup[partner_key]
            # Only process each pair once
            pair_tuple = minmax(idx, partner_idx)
            if !haskey(pair_to_pair_id, pair_tuple)
                next_pair_id += one(UInt32)
                pair_to_pair_id[pair_tuple] = next_pair_id
            end
            pair_id = pair_to_pair_id[pair_tuple]

            # Link each row to its partner
            partner_indices[idx] = partner_idx
            partner_indices[partner_idx] = idx
            
            # Assign the same pair ID to both
            pair_ids[idx] = pair_id
            pair_ids[partner_idx] = pair_id
        else
            pair_ids[idx] = next_pair_id
            next_pair_id += one(UInt32)
            missed_rows += 1
        end
    end

    for i in 1:n 
        if ismissing(pair_ids[i])
            pair_ids[i] = next_pair_id
            next_pair_id += one(UInt32)
        end
    end
    
    # Add both columns to the result
    result_df.partner_precursor_idx = partner_indices
    result_df.pair_id = pair_ids
    
    return result_df, lookup
end
=#

#=
BuildSpecLib("/Users/nathanwamsley/Desktop/test_params.json")
ptable = DataFrame(Tables.columntable(Arrow.Table("/Users/nathanwamsley/Desktop/keap1_library.poin/precursors.arrow")));
#sort!(ptable,:pair_id)
#sort!(ptable,:base_prec_id)
ptable[1,[:proteome_identifiers,:accession_number,:sequence,:mod_key,:precursor_charge,:decoy,:entrapment_group_id, :pair_id,:base_prec_id,:base_pep_id,:precursor_pair_idx]]
ptable[ptable[1,:precursor_pair_idx],[:proteome_identifiers,:accession_number,:sequence,:mod_key,:precursor_charge,:decoy,:entrapment_group_id, :pair_id,:base_prec_id,:base_pep_id,:precursor_pair_idx]]
value_counts(df, col) = combine(groupby(df, col), nrow)
value_counts(ptable, :pair_id)
pair_id_counts = unique(value_counts(ptable, :pair_id)[!,:nrow])
@test length(pair_id_counts)==1
@test pair_id_counts[1]==2
=#

function add_charge_specific_partner_columns!(df::DataFrame)
    """
    Creates target-decoy pairs using (base_entrap_id, precursor_charge) combination.
    Each base_entrap_id/precursor_charge combo should have exactly one target and exactly one decoy.
    Uses combination of base_entrap_id and precursor_charge to assign unique pair_ids.
    """
    result_df = copy(df)
    n = nrow(df)
    
    # Diagnostic: count targets and decoys
    targets = sum(.!df.decoy)
    decoys = sum(df.decoy)
    @user_info "Starting pairing by (base_entrap_id, precursor_charge): $targets targets, $decoys decoys"
    
    # Create lookup dictionary: (base_entrap_id, precursor_charge, is_target) -> row_index
    lookup = Dict{Tuple{UInt32, UInt8, Bool}, Int}()
    for (idx, row) in enumerate(eachrow(df))
        key = (row.base_entrap_id, row.precursor_charge, !row.decoy)  # !decoy = is_target
        lookup[key] = idx
    end
    
    @user_info "Created lookup with $(length(lookup)) entries"
    
    # Initialize pair_id vector (will be completely reassigned)
    new_pair_ids = Vector{Union{UInt32, Missing}}(missing, n)
    next_pair_id = UInt32(1)
    
    # Track processed pairs to avoid duplicates
    processed_pairs = Set{Tuple{Int, Int}}()
    paired_decoys = Set{Int}()  # Track which decoys got paired
    
    for idx in 1:n
        row = df[idx, :]
        
        # Skip if already processed or if it's a decoy (process targets first)
        if row.decoy
            continue
        end
        
        # Look for decoy partner with same base_entrap_id and precursor_charge
        decoy_key = (row.base_entrap_id, row.precursor_charge, false)   # is_target = false
        
        if haskey(lookup, decoy_key)
            decoy_idx = lookup[decoy_key]
            
            # Avoid double processing
            pair_tuple = (min(idx, decoy_idx), max(idx, decoy_idx))
            if pair_tuple ∉ processed_pairs
                # Assign same pair_id to both
                new_pair_ids[idx] = next_pair_id
                new_pair_ids[decoy_idx] = next_pair_id
                next_pair_id += 1
                
                push!(processed_pairs, pair_tuple)
                push!(paired_decoys, decoy_idx)
            end
        else
            @debug_l2 "No decoy found for target: base_entrap_id=$(row.base_entrap_id), precursor_charge=$(row.precursor_charge)"
            # Unpaired target gets unique pair_id
            new_pair_ids[idx] = next_pair_id
            next_pair_id += 1
        end
    end
    
    # Find unpaired decoys and assign them unique pair_ids
    unpaired_decoys = 0
    
    for idx in 1:n
        row = df[idx, :]
        if row.decoy && idx ∉ paired_decoys
            # This decoy has no target partner - assign unique pair_id
            new_pair_ids[idx] = next_pair_id
            next_pair_id += 1
            unpaired_decoys += 1
        end
    end
    
    # Handle any remaining missing entries
    for idx in 1:n
        if ismissing(new_pair_ids[idx])
            new_pair_ids[idx] = next_pair_id
            next_pair_id += 1
        end
    end
    
    # Diagnostic reporting
    if unpaired_decoys > 0
        @user_warn "Found $unpaired_decoys decoy precursors without matching target partners"
    end
    
    successful_pairs = length(processed_pairs)
    @user_info "Three-tier pairing complete: $successful_pairs target-decoy pairs matched using (base_pep_id, base_prec_id)"
    
    # Update the DataFrame with new pair_ids
    result_df.pair_id = new_pair_ids
    
    return result_df
end

function add_entrapment_partner_columns!(df::DataFrame)
    """
    Creates entrapment pairing system - Stage 1: entrapment_pair_id only.
    - entrapment_pair_id: Groups original target with all its entrapment variants
    
    Uses (base_pep_id, precursor_charge) as the key for entrapment pairing.
    Decoys get 'missing' for entrapment_pair_id.
    
    Note: entrapment_target_idx is created later in Stage 2 by add_entrapment_indices!
    This mirrors the two-stage process used for pair_id -> partner_precursor_idx.
    
    Logic:
    - All entrapment variants of same target get same entrapment_pair_id
    - entrapment_target_idx will be added later with correct row indices
    """
    n = nrow(df)
    
    # Initialize entrapment_pair_id column only
    entrapment_pair_ids = Vector{Union{UInt32, Missing}}(missing, n)
    
    # Create lookup for targets only (exclude decoys)
    target_lookup = Dict{Tuple{UInt32, UInt8}, UInt32}()
    
    # First pass: identify all original targets (entrapment_group_id == 0, not decoy)
    for idx in 1:n
        row = df[idx, :]
        if !row.decoy && row.entrapment_group_id == 0
            key = (row.base_pep_id, row.precursor_charge)
            target_lookup[key] = UInt32(idx)
        end
    end
    
    @user_info "Found $(length(target_lookup)) original targets for entrapment pairing"
    
    # Track next pair ID
    next_pair_id = UInt32(1)
    pair_id_map = Dict{UInt32, UInt32}()  # target_idx -> pair_id
    
    # Second pass: assign pair IDs only
    for idx in 1:n
        row = df[idx, :]
        
        # Skip decoys - leave as missing
        if row.decoy
            continue
        end
        
        key = (row.base_pep_id, row.precursor_charge)
        
        if row.entrapment_group_id == 0
            # Original target
            if !haskey(pair_id_map, idx)
                pair_id_map[idx] = next_pair_id
                next_pair_id += 1
            end
            entrapment_pair_ids[idx] = pair_id_map[idx]
        else
            # Entrapment sequence
            if haskey(target_lookup, key)
                target_idx = target_lookup[key]
                
                # Get or assign pair_id for this target
                if !haskey(pair_id_map, target_idx)
                    pair_id_map[target_idx] = next_pair_id
                    next_pair_id += 1
                end
                
                entrapment_pair_ids[idx] = pair_id_map[target_idx]
            else
                @debug_l2 "No target found for entrapment: base_pep_id=$(row.base_pep_id), charge=$(row.precursor_charge)"
            end
        end
    end
    
    # Add entrapment_pair_id column to DataFrame (entrapment_target_idx added later)
    df.entrapment_pair_id = entrapment_pair_ids
    
    # Log statistics
    n_targets = sum(df.entrapment_group_id .== 0 .&& .!df.decoy)
    n_entrapments = sum(df.entrapment_group_id .> 0)
    n_paired_entrapments = sum(.!ismissing.(entrapment_pair_ids[df.entrapment_group_id .> 0]))
    n_decoys = sum(df.decoy)
    
    @user_info "Entrapment pairing Stage 1 complete (entrapment_pair_id created):"
    @user_info "  Original targets: $n_targets"
    @user_info "  Entrapment sequences: $n_entrapments" 
    @user_info "  Successfully paired entrapments: $n_paired_entrapments"
    @user_info "  Decoys (set to missing): $n_decoys"
    
    return df
end

function add_entrapment_indices!(df)
    """
    Add entrapment_target_idx column based on entrapment_pair_id values - Stage 2.
    Called AFTER loading from Arrow to ensure correct row indices.
    Mirrors the approach used by add_pair_indices! for partner_precursor_idx.
    
    Logic:
    - Find all original targets (entrapment_group_id == 0, not decoy if column exists)
    - Map entrapment_pair_id to target row index
    - Apply mapping to all rows with same entrapment_pair_id
    - Decoys remain as missing (if decoy column exists)
    
    Example:
    Row | entrapment_pair_id | entrapment_group_id | decoy | -> entrapment_target_idx
    ----|-------------------|-------------------|-------|---------------------
    1   | 1                 | 0                 | false | 1 (target, points to self)
    2   | 1                 | 1                 | false | 1 (entrap -> target)  
    3   | 1                 | 2                 | false | 1 (entrap -> target)
    4   | missing           | 0                 | true  | missing (decoy)
    """
    n = nrow(df)
    
    # Check if decoy column exists (be robust about column names)
    has_decoy_col = hasproperty(df, :decoy)
    if !has_decoy_col
        @user_warn "Decoy column not found - proceeding without decoy filtering"
        @user_warn "Available columns: $(names(df))"
    end
    
    # Create mapping: entrapment_pair_id -> target row index
    pair_to_target = Dict{UInt32, UInt32}()
    
    # First pass: find all original targets
    for i in 1:n
        # Check if this is a valid original target
        is_target = !ismissing(df.entrapment_pair_id[i]) && 
                   df.entrapment_group_id[i] == 0
        
        # Additional decoy filtering if column exists
        if has_decoy_col && is_target
            is_target = is_target && !df.decoy[i]
        end
        
        if is_target
            # This is an original target - map its pair_id to its row
            pair_to_target[df.entrapment_pair_id[i]] = UInt32(i)
        end
    end
    
    @user_info "Found $(length(pair_to_target)) entrapment targets for index mapping"
    
    # Create the entrapment_target_idx column
    entrapment_target_idx = Vector{Union{UInt32, Missing}}(missing, n)
    
    # Second pass: assign target indices based on pair_id
    n_mapped = 0
    for i in 1:n
        # Check if this row should get an entrapment_target_idx
        should_map = !ismissing(df.entrapment_pair_id[i])
        
        # Additional decoy filtering if column exists
        if has_decoy_col && should_map
            should_map = should_map && !df.decoy[i]
        end
        
        if should_map
            # Look up the target row for this entrapment_pair_id
            if haskey(pair_to_target, df.entrapment_pair_id[i])
                entrapment_target_idx[i] = pair_to_target[df.entrapment_pair_id[i]]
                n_mapped += 1
            end
        end
        # Decoys and missing entrapment_pair_id remain as missing
    end
    
    # Add the column to the DataFrame
    df[!, :entrapment_target_idx] = entrapment_target_idx
    
    @user_info "Entrapment pairing Stage 2 complete (entrapment_target_idx created):"
    @user_info "  Mapped $n_mapped entries to target indices"
    @user_info "  Max target index: $(isempty(pair_to_target) ? 0 : maximum(values(pair_to_target)))"
    @user_info "  Table size: $n rows"
    @user_info "  Decoy column used: $(has_decoy_col ? "✅ YES" : "❌ NO")"
    
    return nothing
end