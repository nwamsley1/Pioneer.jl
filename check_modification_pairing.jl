#!/usr/bin/env julia

using Arrow, Tables, DataFrames

# Load the library data
df = DataFrame(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow"))

println("🔍 COMPREHENSIVE MODIFICATION PAIRING CHECK:")

# Group by pair_id and check for mismatched modifications
paired_groups = groupby(df, :pair_id)
mismatched_pairs = 0
total_pairs_checked = 0

for group in paired_groups
    if nrow(group) == 2  # Only check actual pairs
        global total_pairs_checked += 1
        row1, row2 = eachrow(group)[1], eachrow(group)[2]
        
        # Check if modifications match
        mod1 = ismissing(row1.structural_mods) ? "" : row1.structural_mods
        mod2 = ismissing(row2.structural_mods) ? "" : row2.structural_mods
        
        if mod1 != mod2
            global mismatched_pairs += 1
            if mismatched_pairs <= 3  # Show first 3 examples
                println("❌ PAIR_ID $(row1.pair_id):")
                println("   Row 1: '$(row1.sequence)' mods: '$mod1' decoy:$(row1.is_decoy)")
                println("   Row 2: '$(row2.sequence)' mods: '$mod2' decoy:$(row2.is_decoy)")
                println()
            end
        end
    end
end

println("FINAL RESULTS:")
println("   Total target-decoy pairs checked: $total_pairs_checked")  
println("   Pairs with mismatched modifications: $mismatched_pairs")

if mismatched_pairs == 0
    println("   ✅ SUCCESS: All target-decoy pairs have matching modifications!")
else
    println("   ❌ ISSUE: $mismatched_pairs pairs still have mismatched modifications")
    percentage = round(mismatched_pairs / total_pairs_checked * 100, digits=2)
    println("   Percentage of problematic pairs: $percentage%")
end