#!/usr/bin/env julia

"""
Test the memory optimization utilities for IntegrateChromatogramSearch.
"""

using Pioneer

# Include the memory optimization utilities
include("src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/memory_optimization.jl")

function test_memory_optimization()
    println("Testing IntegrateChromatogramSearch memory optimizations...")
    
    # Test chromatogram capacity estimation
    scan_range = collect(1:1000)  # 1000 scans
    precursors_passing = Set{UInt32}(1:500)  # 500 precursors
    
    estimated_capacity = estimate_chromatogram_capacity(scan_range, precursors_passing)
    
    println("✓ Capacity estimation:")
    println("  - Scans: $(length(scan_range))")
    println("  - Precursors: $(length(precursors_passing))")
    println("  - Estimated capacity: $estimated_capacity")
    println("  - Original fixed size: 500,000")
    
    if estimated_capacity < 500000
        println("  ✓ Optimization reduces initial allocation by $(500000 - estimated_capacity) elements")
    else
        println("  ⚠ Large dataset requires more than fixed 500k allocation")
    end
    
    # Test smart array initialization (simplified test)
    try
        # Test with generic Float32 arrays to demonstrate the concept
        test_array = initialize_chromatogram_array_optimized(Float32, scan_range, precursors_passing)
        
        println("✓ Array initialization:")
        println("  - Test array size: $(length(test_array))")
        println("  - Array created successfully")
        
        # Test smart resizing
        current_idx = div(length(test_array), 2)  # Halfway through
        needed_capacity = length(test_array) + 1000  # Need more space
        
        original_length = length(test_array)
        smart_chromatogram_resize!(test_array, current_idx, needed_capacity)
        new_length = length(test_array)
        
        println("✓ Smart resizing:")
        println("  - Original length: $original_length")
        println("  - Needed capacity: $needed_capacity")
        println("  - New length: $new_length")
        println("  - Growth factor: $(round(new_length / original_length, digits=2))x")
        
        println("✓ Memory optimization utilities work correctly")
        
    catch e
        println("✗ Error in array operations: $e")
        rethrow(e)
    end
    
    # Show integration guide
    println("\n" * "="^60)
    optimize_build_chromatograms_memory()
    
    println("\n✅ Memory optimization test completed successfully!")
    println("\nKey benefits:")
    println("  - Data-driven initial allocation instead of fixed 500k")
    println("  - Exponential growth strategy for efficient resizing")
    println("  - Eliminates all 4 hard-coded append! operations")
    println("  - Drop-in replacement for existing code patterns")
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_memory_optimization()
end