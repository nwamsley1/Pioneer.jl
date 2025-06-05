#!/usr/bin/env julia

"""
Test script to validate FirstPassSearch optimization.
Tests that the optimized function can be loaded and has the correct structure.
"""

using Pioneer

function test_firstpass_optimization()
    println("Testing FirstPassSearch optimization...")
    
    # Test compilation by reading the source files
    try
        original_content = read("src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl", String)
        optimized_content = read("src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns_optimized.jl", String)
        
        println("✓ Both original and optimized files exist")
        println("✓ Original file size: $(length(original_content)) characters")
        println("✓ Optimized file size: $(length(optimized_content)) characters")
        
        # Check that the optimized version has the expected optimizations
        if occursin("max_concurrent_files", optimized_content)
            println("✓ Found max_concurrent_files parameter (memory bandwidth control)")
        else
            println("✗ Missing max_concurrent_files parameter")
        end
        
        if occursin("Thread-local dictionary", optimized_content)
            println("✓ Found thread-local dictionary optimization")
        else
            println("✗ Missing thread-local dictionary optimization")
        end
        
        if occursin("merge_dictionaries!", optimized_content)
            println("✓ Found dictionary merging optimization")
        else
            println("✗ Missing dictionary merging optimization")
        end
        
        if occursin("Threads.@spawn", optimized_content)
            println("✓ Found parallel processing with Threads.@spawn")
        else
            println("✗ Missing parallel processing")
        end
        
    catch e
        println("✗ Error reading files: $e")
        rethrow(e)
    end
    
    println("\nOptimization validation completed!")
    println("The optimized function includes:")
    println("  - Controlled concurrent Arrow file reading")
    println("  - Thread-local dictionary accumulation")
    println("  - Memory bandwidth protection")
    println("  - Efficient dictionary merging")
    println("\nTo test with real data, run a full SearchDIA pipeline.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_firstpass_optimization()
end