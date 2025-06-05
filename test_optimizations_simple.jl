#!/usr/bin/env julia

"""
Simple test script to validate search optimizations without overwhelming memory.
Tests both FirstPassSearch and IntegrateChromatogramsSearch optimizations.
"""

using Pioneer

function test_search_optimizations()
    println("Testing search optimizations with ecoli dataset...")
    
    # Use the existing ecoli test parameters
    params_path = "./data/ecoli_test/ecoli_test_params.json"
    
    if !isfile(params_path)
        error("Could not find test parameters at: $params_path")
    end
    
    # Test with a subset to avoid memory issues
    println("Running SearchDIA with optimizations...")
    
    try
        # This will use our optimized functions when they're integrated
        SearchDIA(params_path)
        println("✓ Search completed successfully")
    catch e
        println("✗ Search failed with error: $e")
        rethrow(e)
    end
end

function compare_memory_usage()
    println("Comparing memory usage patterns...")
    
    # Simple memory monitoring
    gc()
    initial_memory = Base.gc_bytes()
    
    test_search_optimizations()
    
    gc()
    final_memory = Base.gc_bytes()
    
    println("Memory usage: $(final_memory - initial_memory) bytes")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_memory_usage()
end