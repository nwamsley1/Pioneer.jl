using Arrow, Tables, DataFrames

df = DataFrame(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow"))

# Check for sequences that share the same base_prec_id but should be different
println("Looking for sequences that share base_prec_id but have different sequences or mods...")
grouped = groupby(df, :base_prec_id)
problematic_count = 0

for group in grouped
    if nrow(group) > 2  # Should only be target-decoy pairs (max 2)
        global problematic_count += 1
        if problematic_count <= 3
            println("\n❌ base_prec_id $(group[1,:base_prec_id]) has $(nrow(group)) entries:")
            for row in eachrow(group)
                println("  Seq: $(row.sequence) | Mods: $(row.structural_mods) | Decoy: $(row.is_decoy)")
            end
        end
    elseif nrow(group) == 2
        row1, row2 = eachrow(group)[1], eachrow(group)[2]
        # Check if they have different sequences or different modifications
        seq_different = row1.sequence != row2.sequence
        mod_different = row1.structural_mods != row2.structural_mods
        
        if seq_different || mod_different
            global problematic_count += 1
            if problematic_count <= 3
                println("\n❌ base_prec_id $(group[1,:base_prec_id]) has mismatched pairs:")
                println("  Target: $(row1.sequence) | Mods: $(row1.structural_mods) | Decoy: $(row1.is_decoy)")
                println("  Decoy:  $(row2.sequence) | Mods: $(row2.structural_mods) | Decoy: $(row2.is_decoy)")
            end
        end
    end
end

println("\nTotal problematic base_prec_id groups: $problematic_count")