#!/usr/bin/env julia

# Test script for target-decoy pairing functionality
using Pioneer
using DataFrames
using Random

# Create test PSM data
function create_test_psms()
    # Create realistic test data with multiple precursors
    psms = DataFrame(
        precursor_idx = UInt32[1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10],
        target = [true, true, true, true, false, false, false, false, false, false, false, false, false],
        irt = [10.5, 10.5, 15.2, 15.2, 12.1, 12.1, 16.8, 25.4, 30.1, 35.5, 40.2, 45.8, 50.3],
        isotopes_captured = [(0,1), (1,2), (0,1), (1,2), (0,1), (1,2), (0,1), (0,1), (1,2), (0,1), (1,2), (0,1), (1,2)],
        pair_id = [1, 1, 2, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8]  # Original library-based pairing
    )
    return psms
end

function main()
    println("Testing Target-Decoy Pairing Implementation")
    println("=" ^ 50)
    
    # Create test data
    psms = create_test_psms()
    
    println("\\nOriginal data:")
    println(psms)
    
    println("\\nSummary before pairing:")
    println("- Total PSMs: ", nrow(psms))
    println("- Unique targets: ", length(unique(psms.precursor_idx[psms.target])))
    println("- Unique decoys: ", length(unique(psms.precursor_idx[.!psms.target])))
    
    # Apply pairing
    println("\\nApplying target-decoy pairing...")
    try
        Pioneer.assign_random_target_decoy_pairs!(psms)
        
        println("\\nPairing completed successfully!")
        println("\\nUpdated data:")
        println(psms)
        
        # Analyze results
        paired_mask = .!ismissing.(psms.pair_id)
        println("\\nResults analysis:")
        println("- PSMs with pair_id: ", sum(paired_mask))
        println("- PSMs without pair_id: ", sum(.!paired_mask))
        println("- Unique pair_ids: ", length(unique(skipmissing(psms.pair_id))))
        
        # Check target/decoy balance in pairs
        if any(paired_mask)
            paired_psms = psms[paired_mask, :]
            pair_groups = groupby(paired_psms, :pair_id)
            
            println("\\nPair validation:")
            for group in pair_groups
                pair_id = group.pair_id[1]
                n_targets = sum(group.target)
                n_decoys = sum(.!group.target)
                target_precursors = unique(group.precursor_idx[group.target])
                decoy_precursors = unique(group.precursor_idx[.!group.target])
                
                println("  Pair $pair_id: $n_targets targets ($(length(target_precursors)) unique), $n_decoys decoys ($(length(decoy_precursors)) unique)")
                
                if length(target_precursors) != 1 || length(decoy_precursors) != 1
                    println("    ⚠️  WARNING: Expected exactly 1 target and 1 decoy precursor")
                else
                    println("    ✅ Valid 1:1 precursor pairing")
                end
            end
        end
        
        println("\\n🎉 Test completed successfully!")
        
    catch e
        println("❌ Error during pairing: ", e)
        println("Stack trace:")
        for (i, frame) in enumerate(stacktrace(catch_backtrace()))
            println("  $i: $frame")
        end
        return false
    end
    
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    success = main()
    exit(success ? 0 : 1)
end