# Validation script for partner_precursor_idx fix
# Run this in Julia REPL after implementing the fix

using Pioneer, DataFrames, Arrow, Tables

println("=== KEAP1 Library Fix Validation ===\n")

# 1. Rebuild library with fix
#println("1. Building library with fix...")
#BuildSpecLib("test_keap1_params.json")
#println("✓ Library build completed\n")

# 2. Load and analyze results
println("2. Loading library for analysis...")
ptable = DataFrame(Tables.columntable(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow")))
println("✓ Library loaded: $(nrow(ptable)) precursors\n")

# 3. Check partner indices are valid
println("3. Validating partner_precursor_idx values...")
max_partner_idx = Int64(maximum(skipmissing(ptable.partner_precursor_idx)))
table_size = nrow(ptable)
valid_indices = max_partner_idx <= table_size

println("   Max partner index: $max_partner_idx")
println("   Table size: $table_size")
println("   Valid indices: $(valid_indices ? "✓ PASS" : "✗ FAIL")")

if !valid_indices
    println("   ❌ ERROR: Partner indices exceed table size!")
    invalid_count = sum(skipmissing(ptable.partner_precursor_idx) .> table_size)
    println("   Invalid indices count: $invalid_count")
end
println()

# 4. Verify pairing still works correctly
println("4. Verifying target-decoy pairing...")
gptable = groupby(ptable, :pair_id)
pair_sizes = [nrow(subdf) for (key, subdf) in pairs(gptable)]
unique_pair_sizes = unique(pair_sizes)

println("   Total pairs: $(length(gptable))")
println("   Pair sizes: $unique_pair_sizes")
println("   Pairing correct: $(unique_pair_sizes == [2] ? "✓ PASS" : "✗ FAIL")")

if unique_pair_sizes != [2]
    println("   ❌ ERROR: Some pairs don't have exactly 2 members!")
    size_counts = [(size, count(==(size), pair_sizes)) for size in unique_pair_sizes]
    for (size, count) in size_counts
        println("   Size $size: $count pairs")
    end
end
println()

# 5. Check target/decoy balance
println("5. Checking target/decoy balance...")
targets = sum(.!ptable.is_decoy)
decoys = sum(ptable.is_decoy)
balanced = targets == decoys

println("   Targets: $targets")
println("   Decoys: $decoys")
println("   Balanced: $(balanced ? "✓ PASS" : "✗ FAIL")")
println()

# 6. Verify specific partner relationships
println("6. Testing partner relationship integrity...")
partner_test_passed = true
test_sample = rand(1:nrow(ptable), min(100, nrow(ptable)))

for row_idx in test_sample
    partner_idx = ptable.partner_precursor_idx[row_idx]
    
    if !ismissing(partner_idx)
        if partner_idx > nrow(ptable)
            println("   ❌ ERROR: Row $row_idx partner index $partner_idx exceeds table size!")
            partner_test_passed = false
            break
        end
        
        # Check if partner points back
        partner_partner_idx = ptable.partner_precursor_idx[partner_idx]
        if !ismissing(partner_partner_idx) && partner_partner_idx != row_idx
            println("   ⚠️  WARNING: Asymmetric partnership: Row $row_idx → $partner_idx → $partner_partner_idx")
        end
        
        # Check if they have same pair_id
        if ptable.pair_id[row_idx] != ptable.pair_id[partner_idx]
            println("   ❌ ERROR: Partners have different pair_ids: $(ptable.pair_id[row_idx]) vs $(ptable.pair_id[partner_idx])")
            partner_test_passed = false
            break
        end
        
        # Check target/decoy relationship
        if ptable.is_decoy[row_idx] == ptable.is_decoy[partner_idx]
            println("   ❌ ERROR: Partners are both $(ptable.is_decoy[row_idx] ? "decoys" : "targets")!")
            partner_test_passed = false
            break
        end
    end
end

println("   Partner relationships: $(partner_test_passed ? "✓ PASS" : "✗ FAIL")")
println()

# 7. Summary
println("=== VALIDATION SUMMARY ===")
all_tests_passed = valid_indices && (unique_pair_sizes == [2]) && balanced && partner_test_passed

if all_tests_passed
    println("🎉 ALL TESTS PASSED!")
    println("   The library is ready for SearchDIA testing")
    println()
    println("Next step: Test with SearchDIA")
    println("   SearchDIA(\"path/to/search_params.json\")")
else
    println("❌ SOME TESTS FAILED!")
    println("   Review the errors above before proceeding")
end

println()
println("=== DETAILED STATISTICS ===")
println("Total precursors: $(nrow(ptable))")
println("Unique base_pep_ids: $(length(unique(ptable.base_pep_id)))")
println("Max base_pep_id: $(maximum(ptable.base_pep_id))")
println("Total pairs: $(length(gptable))")
println("Missing partners: $(sum(ismissing.(ptable.partner_precursor_idx)))")
println("Library size: $(round(Base.summarysize(ptable) / 1024^2, digits=2)) MB")