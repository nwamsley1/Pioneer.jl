#=
using DataFrames
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
    println("missed_rows $missed_rows")
    
    # Add both columns to the result
    result_df.partner_precursor_idx = partner_indices
    result_df.pair_id = pair_ids
    
    return result_df, lookup
end
#=
function add_precursor_partner_columns!(df::DataFrame)
    # Add better debugging and error handling
    result_df = copy(df)
    
    # Create a lookup dictionary with more robust handling
    lookup = Dict()
    for (idx, row) in enumerate(eachrow(df))
        # Convert any special types to standard types to avoid comparison issues
        key = (
            String(row.proteome_identifiers), 
            String(row.sequence), 
            row.structural_mods, 
            row.prec_charge
        )
        lookup[key] = UInt32(idx)
    end
    
    # Create vectors for partner indices and pair IDs
    n = nrow(df)
    partner_indices = Vector{Union{UInt32, Missing}}(missing, n)
    pair_ids = Vector{Union{UInt32, Missing}}(missing, n)
    
    # Track which pairs have been processed
    pair_to_pair_id = Dict{Tuple{UInt32, UInt32}, UInt32}()
    next_pair_id = UInt32(0)
    missed_rows = 0
    missed_details = Dict{String, Int}()
    
    # Process each row to find partners
    for idx in 1:n
        # Get current row data
        row = df[idx, :]
        if row.is_decoy == true
            continue
        end
        
        # Calculate what the partner sequence would be
        try
            partner_seq, partner_mods = reverseSequence(row.sequence, row.structural_mods)
            
            # Look up potential partner with opposite is_decoy status
            partner_key = (
                String(row.proteome_identifiers), 
                partner_seq, 
                partner_mods, 
                row.prec_charge
            )
            
            if haskey(lookup, partner_key)
                partner_idx = lookup[partner_key]
                partner_row = df[partner_idx, :]
                
                # Verify the partner is actually a decoy
                if !partner_row.is_decoy
                    # Skip this pair - the partner isn't a decoy
                    reason = "Partner is not a decoy"
                    missed_details[reason] = get(missed_details, reason, 0) + 1
                    missed_rows += 1
                    continue
                end
                
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
                # No partner found - increment counter and track reason
                reason = "No matching partner key"
                missed_details[reason] = get(missed_details, reason, 0) + 1
                missed_rows += 1
                
                # Debug a sample of missing matches
                if rand() < 0.01  # Sample 1% of missed rows for debugging
                    println("Missing partner for row $idx:")
                    println("  Sequence: $(row.sequence)")
                    println("  Expected partner: $partner_seq")
                    println("  Mods: $(row.structural_mods)")
                    println("  Partner mods: $partner_mods")
                end
            end
        catch e
            # Track errors in reverseSequence
            reason = "Error: $(typeof(e))"
            missed_details[reason] = get(missed_details, reason, 0) + 1
            missed_rows += 1
            println("Error processing row $idx: $e")
        end
    end
    
    # Print detailed diagnostic information
    println("Total missed rows: $missed_rows")
    println("Reasons for missed rows:")
    for (reason, count) in missed_details
        println("  $reason: $count")
    end
    
    # Check for modification patterns in unmatched rows
    unmatched_target_indices = findall(row -> 
        (ismissing(pair_ids[row]) && !df[row, :is_decoy]), 1:n)
    
    if !isempty(unmatched_target_indices)
        println("Examining modification patterns in $(min(5, length(unmatched_target_indices))) unmatched target rows:")
        for i in unmatched_target_indices[1:min(5, length(unmatched_target_indices))]
            println("Row $i: $(df[i, :sequence]), Mods: $(df[i, :structural_mods])")
        end
    end
    
    # Add both columns to the result
    result_df.partner_precursor_idx = partner_indices
    result_df.pair_id = pair_ids
    
    return result_df, lookup
end
=#
ptable_path = "/Users/n.t.wamsley/Documents/PIONEER_PAPER/SPEC_LIBS/altimeter_3proteome_olsenAstral_paired_020625.poin/precursors_table.arrow"
precursors_table = DataFrame(Tables.columntable(Arrow.Table(
    ptable_path
    )))
    
precursors_table, lookup_table = add_precursor_partner_columns!(precursors_table)
temp_dir = "/Users/n.t.wamsley/Documents/PIONEER_PAPER/SPEC_LIBS/altimeter_3proteome_olsenAstral_paired_020625.poin/precursors_table_temp.arrow"
Arrow.write(temp_dir, precursors_table)
mv(temp_dir, ptable_path, force = true)

ptable_path = "/Users/n.t.wamsley/Documents/PIONEER_PAPER/SPEC_LIBS/altimeter_yeast_olsenAstral_paired_020725.poin/precursors_table.arrow"
precursors_table = DataFrame(Tables.columntable(Arrow.Table(
    ptable_path
    )))
precursors_table, lookup_table = add_precursor_partner_columns!(precursors_table)
temp_dir = "/Users/n.t.wamsley/Documents/PIONEER_PAPER/SPEC_LIBS/altimeter_yeast_olsenAstral_paired_020725.poin/precursors_table_temp.arrow"
Arrow.write(temp_dir, precursors_table)
mv(temp_dir, ptable_path, force = true)



ptable_path = "/Users/nathanwamsley/Projects/test_pioneer/Pioneer.jl/data/ecoli_test/Prosit_ECOLI_500_600mz_062625.poin/precursors_table.arrow"
precursors_table = DataFrame(Tables.columntable(Arrow.Table(
    ptable_path
    )))
precursors_table, lookup_table = add_precursor_partner_columns!(precursors_table)
temp_dir = "/Users/nathanwamsley/Projects/test_pioneer/Pioneer.jl/data/ecoli_test/Prosit_ECOLI_500_600mz_062625.poin/precursors_table_temp.arrow"
Arrow.write(temp_dir, precursors_table)
mv(temp_dir, ptable_path, force = true)



#Testing 
#=

=#
=#