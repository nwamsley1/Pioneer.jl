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

"""
Performance testing framework for file merging operations.

This module provides comprehensive testing and benchmarking capabilities for
comparing different merge implementations:
1. Current FileOperations.jl implementation
2. Type-stable implementation  
3. Recursive merge strategy

Usage examples:
```julia
include("test/performance/test_merge_performance.jl")
using .MergePerformanceTests

# Test with the protein groups data
test_dir = "/Users/nathanwamsley/Data/test_io/passing_proteins"

# Quick accuracy test
test_accuracy(test_dir, n_files=10)

# Performance comparison
benchmark_merge_strategies(test_dir, file_counts=[10, 50, 100])

# Full test suite
run_all_tests(test_dir)
```
"""

module MergePerformanceTests

using Arrow, DataFrames, Tables
using Statistics
using Printf
using Dates

# Simple timing function to replace BenchmarkTools
macro simple_time(expr)
    quote
        local start_time = time()
        local result = $(esc(expr))
        local elapsed = time() - start_time
        elapsed
    end
end

# Include the modules we want to test
include("../../src/utils/FileOperations/FileOperations.jl")
include("../../src/utils/writeArrow.jl")

export test_accuracy, benchmark_merge_strategies, run_all_tests, 
       profile_memory_usage, test_file_handle_limits

#==========================================================
Test Data Utilities
==========================================================#

"""
Get list of .arrow files from test directory, limited to n_files.
"""
function get_test_files(test_dir::String, n_files::Int=10)
    all_files = [joinpath(test_dir, f) for f in readdir(test_dir) if endswith(f, ".arrow")]
    if length(all_files) < n_files
        @warn "Only $(length(all_files)) files available, using all of them"
        return all_files
    end
    return all_files[1:n_files]
end

"""
Create FileReferences from file paths.
"""
function create_file_references(file_paths::Vector{String})
    return [ProteinGroupFileReference(path) for path in file_paths]
end

"""
Analyze a single Arrow file to understand its schema and content.
"""
function analyze_file(file_path::String)
    table = Arrow.Table(file_path)
    
    println("File: $(basename(file_path))")
    println("  Rows: $(length(Tables.getcolumn(table, 1)))")
    println("  Columns: $(length(Tables.columnnames(table)))")
    
    # Show column info
    for col_name in Tables.columnnames(table)
        col = Tables.getcolumn(table, col_name)
        col_type = eltype(col)
        println("    $col_name: $col_type")
    end
    println()
end

#==========================================================
Accuracy Testing
==========================================================#

"""
Test that different implementations produce identical results.
"""
function test_accuracy(test_dir::String; n_files::Int=10, verbose::Bool=true)
    println("=== ACCURACY TEST ===")
    println("Testing with $n_files files from $test_dir")
    
    file_paths = get_test_files(test_dir, n_files)
    refs = create_file_references(file_paths)
    
    if verbose
        println("Analyzing first file to understand schema:")
        analyze_file(file_paths[1])
    end
    
    # Create temporary output files
    temp_dir = mktempdir(prefix="merge_accuracy_test_")
    
    try
        # Test different sort key combinations commonly used in protein groups
        test_cases = [
            # 2-key tests (backward compatibility)
            ((:pg_score, :target), [true, true], "pg_score_target_desc"),
            ((:global_pg_score, :target), [true, true], "global_pg_score_target_desc"),
            ((:protein_name, :target), [false, false], "protein_name_target_asc"),
            
            # 3-key tests  
            ((:pg_score, :target, :protein_name), [true, true, false], "3key_mixed"),
            
            # Test single boolean reverse (should broadcast)
            ((:pg_score, :target), true, "2key_all_desc_single_bool"),
            ((:protein_name, :target), false, "2key_all_asc_single_bool")
        ]
        
        for (sort_keys, reverse, test_name) in test_cases
            println("\n--- Testing sort keys: $sort_keys, reverse: $reverse ---")
            
            # Skip if columns don't exist
            first_table = Arrow.Table(file_paths[1])
            if !all(haskey(Tables.columntable(first_table), key) for key in sort_keys)
                println("  Skipping - columns not found in test data")
                continue
            end
            
            # Test current implementation
            current_output = joinpath(temp_dir, "current_$(test_name).arrow")
            println("  Testing current implementation...")
            
            current_success = false
            current_rows = 0
            try
                @time current_ref = stream_sorted_merge(
                    refs, current_output, sort_keys...;
                    reverse=reverse, batch_size=100_000
                )
                current_success = true
                current_rows = row_count(current_ref)
            catch e
                println("    Current implementation failed: $e")
                current_success = false
                current_rows = 0
            end
            
            # Test type-stable implementation
            typed_output = joinpath(temp_dir, "typed_$(test_name).arrow")
            println("  Testing type-stable implementation...")
            
            typed_success = false
            typed_rows = 0
            try
                @time typed_ref = stream_sorted_merge(
                    refs, typed_output, sort_keys...;
                    reverse=reverse, batch_size=100_000
                )
                typed_success = true
                typed_rows = row_count(typed_ref)
            catch e
                println("    Type-stable implementation failed: $e")
                typed_success = false
                typed_rows = 0
            end
            
            # Test recursive implementation
            recursive_output = joinpath(temp_dir, "recursive_$(test_name).arrow")
            println("  Testing recursive implementation...")
            
            recursive_success = false
            recursive_rows = 0
            try
                @time recursive_ref = recursive_merge(
                    refs, recursive_output, sort_keys...;
                    reverse=reverse, group_size=3, batch_size=100_000
                )
                recursive_success = true
                recursive_rows = row_count(recursive_ref)
            catch e
                println("    Recursive implementation failed: $e")
                recursive_success = false
                recursive_rows = 0
            end
            
            # Compare results
            if current_success && typed_success
                println("  Row counts - Current: $current_rows, Typed: $typed_rows")
                if current_rows == typed_rows
                    println("  ✓ Row counts match")
                else
                    println("  ✗ Row counts differ!")
                end
                
                # Compare first few rows for detailed verification
                if current_rows > 0 && typed_rows > 0
                    compare_file_contents(current_output, typed_output, sort_keys, 5)
                end
            end
            
            if recursive_success
                println("  Recursive rows: $recursive_rows")
                if current_success && current_rows == recursive_rows
                    println("  ✓ Recursive row count matches")
                end
            end
        end
        
    finally
        # Clean up
        rm(temp_dir, recursive=true)
        println("\nAccuracy test completed.")
    end
end

"""
Compare contents of two Arrow files in detail.
"""
function compare_file_contents(file1::String, file2::String, sort_keys::Tuple, n_rows::Int=5)
    table1 = Arrow.Table(file1)
    table2 = Arrow.Table(file2)
    
    df1 = DataFrame(table1)
    df2 = DataFrame(table2)
    
    println("    Comparing first $n_rows rows:")
    
    for i in 1:min(n_rows, nrow(df1), nrow(df2))
        row1 = df1[i, :]
        row2 = df2[i, :]
        
        # Check sort key values
        match = true
        for key in sort_keys
            if haskey(row1, key) && haskey(row2, key)
                val1, val2 = row1[key], row2[key]
                if val1 != val2
                    println("      Row $i: $key differs ($val1 vs $val2)")
                    match = false
                end
            end
        end
        
        if match && i <= 3  # Show first few matching rows
            println("      Row $i: $(sort_keys[1])=$(row1[sort_keys[1]]), $(sort_keys[2])=$(row1[sort_keys[2]]) ✓")
        end
    end
end

#==========================================================
Performance Benchmarking
==========================================================#

"""
Benchmark different merge strategies with varying numbers of files.
"""
function benchmark_merge_strategies(test_dir::String; file_counts::Vector{Int}=[10, 25, 50], 
                                   batch_size::Int=100_000, group_size::Int=10)
    println("=== PERFORMANCE BENCHMARK ===")
    
    all_files = [joinpath(test_dir, f) for f in readdir(test_dir) if endswith(f, ".arrow")]
    println("Total files available: $(length(all_files))")
    
    # Determine sort keys by examining first file
    first_table = Arrow.Table(all_files[1])
    available_cols = Set(Tables.columnnames(first_table))
    
    # Choose appropriate sort keys based on available columns
    sort_keys = if :pg_score in available_cols && :target in available_cols
        (:pg_score, :target)
    elseif :protein_name in available_cols && :target in available_cols
        (:protein_name, :target)
    else
        # Use first two columns as fallback
        cols = collect(Tables.columnnames(first_table))
        (cols[1], cols[2])
    end
    
    reverse = [true, true]  # Typically want descending sorts for scores
    
    println("Using sort keys: $sort_keys with reverse=$reverse")
    println("Batch size: $batch_size, Group size: $group_size")
    
    results = []
    
    for n_files in file_counts
        if n_files > length(all_files)
            println("\nSkipping $n_files files (only $(length(all_files)) available)")
            continue
        end
        
        println("\n--- Benchmarking with $n_files files ---")
        
        file_paths = all_files[1:n_files]
        refs = create_file_references(file_paths)
        temp_dir = mktempdir(prefix="merge_benchmark_")
        
        try
            # Benchmark current implementation
            current_time = benchmark_current_implementation(
                refs, temp_dir, sort_keys, reverse, batch_size
            )
            
            # Benchmark type-stable implementation
            typed_time = benchmark_typed_implementation(
                refs, temp_dir, sort_keys, reverse, batch_size
            )
            
            # Benchmark recursive implementation
            recursive_time = benchmark_recursive_implementation(
                refs, temp_dir, sort_keys, reverse, batch_size, group_size
            )
            
            # Store results
            push!(results, (
                n_files = n_files,
                current_time = current_time,
                typed_time = typed_time,
                recursive_time = recursive_time
            ))
            
            # Print summary
            println("  Results:")
            if current_time > 0
                println(@sprintf("    Current:   %.3f seconds", current_time))
            else
                println("    Current:   FAILED")
            end
            
            if typed_time > 0
                println(@sprintf("    Typed:     %.3f seconds", typed_time))
                if current_time > 0
                    speedup = current_time / typed_time
                    println(@sprintf("    Speedup:   %.2fx", speedup))
                end
            else
                println("    Typed:     FAILED")
            end
            
            if recursive_time > 0
                println(@sprintf("    Recursive: %.3f seconds", recursive_time))
                if current_time > 0
                    speedup = current_time / recursive_time
                    println(@sprintf("    vs Current: %.2fx", speedup))
                end
            else
                println("    Recursive: FAILED")
            end
            
        finally
            #@info "temp_dir: $temp_dir"
            rm(temp_dir, recursive=true)
        end
    end
    
    # Print overall summary
    println("\n=== BENCHMARK SUMMARY ===")
    println(@sprintf("%-10s %-12s %-12s %-12s %-10s %-10s", 
                    "Files", "Current(s)", "Typed(s)", "Recursive(s)", "T-Speedup", "R-Speedup"))
    println("-" ^ 80)
    
    for result in results
        current_str = result.current_time > 0 ? @sprintf("%.3f", result.current_time) : "FAILED"
        typed_str = result.typed_time > 0 ? @sprintf("%.3f", result.typed_time) : "FAILED"
        recursive_str = result.recursive_time > 0 ? @sprintf("%.3f", result.recursive_time) : "FAILED"
        
        typed_speedup = (result.current_time > 0 && result.typed_time > 0) ? 
                       @sprintf("%.2fx", result.current_time / result.typed_time) : "-"
        recursive_speedup = (result.current_time > 0 && result.recursive_time > 0) ?
                           @sprintf("%.2fx", result.current_time / result.recursive_time) : "-"
        
        println(@sprintf("%-10d %-12s %-12s %-12s %-10s %-10s", 
                        result.n_files, current_str, typed_str, recursive_str, 
                        typed_speedup, recursive_speedup))
    end
    
    return results
end

function benchmark_current_implementation(refs, temp_dir, sort_keys, reverse, batch_size)
    output_path = joinpath(temp_dir, "current.arrow")
    try
        elapsed = @simple_time begin
            stream_sorted_merge(
                refs, output_path, sort_keys...;
                reverse=reverse, batch_size=batch_size
            )
        end
        return elapsed
    catch e
        println("    Current implementation error: $e")
        return -1.0
    end
end

function benchmark_typed_implementation(refs, temp_dir, sort_keys, reverse, batch_size)
    output_path = joinpath(temp_dir, "typed.arrow")
    try
        elapsed = @simple_time begin
            stream_sorted_merge(
                refs, output_path, sort_keys...;
                reverse=reverse, batch_size=batch_size
            )
        end
        return elapsed
    catch e
        println("    Typed implementation error: $e")
        return -1.0
    end
end

function benchmark_recursive_implementation(refs, temp_dir, sort_keys, reverse, batch_size, group_size)
    output_path = joinpath(temp_dir, "recursive.arrow")
    try
        elapsed = @simple_time begin
            recursive_merge(
                refs, output_path, sort_keys...;
                reverse=reverse, batch_size=batch_size, group_size=group_size
            )
        end
        return elapsed
    catch e
        println("    Recursive implementation error: $e")
        return -1.0
    end
end

#==========================================================
Memory Usage Testing
==========================================================#

"""
Profile memory usage during merge operations.
"""
function profile_memory_usage(test_dir::String; n_files::Int=50)
    println("=== MEMORY USAGE PROFILING ===")
    println("Note: This provides rough estimates. Use external tools for precise measurements.")
    
    file_paths = get_test_files(test_dir, n_files)
    refs = create_file_references(file_paths)
    
    # Get baseline memory usage
    GC.gc()
    baseline_memory = Sys.total_memory() - Sys.free_memory()
    
    temp_dir = mktempdir(prefix="memory_test_")
    
    try
        # Test current implementation
        println("\nTesting current implementation memory usage...")
        GC.gc()
        before_memory = Sys.total_memory() - Sys.free_memory()
        
        output_path = joinpath(temp_dir, "current_memory.arrow")
        @time stream_sorted_merge(refs, output_path, :pg_score, :target; 
                                 reverse=[true, true], batch_size=50_000)
        
        GC.gc()
        after_memory = Sys.total_memory() - Sys.free_memory()
        current_usage = (after_memory - before_memory) / (1024^3)  # GB
        
        # Test type-stable implementation
        println("\nTesting type-stable implementation memory usage...")
        GC.gc()
        before_memory = Sys.total_memory() - Sys.free_memory()
        
        output_path = joinpath(temp_dir, "typed_memory.arrow")
        @time stream_sorted_merge(refs, output_path, :pg_score, :target; 
                                 reverse=[true, true], batch_size=50_000)
        
        GC.gc()
        after_memory = Sys.total_memory() - Sys.free_memory()
        typed_usage = (after_memory - before_memory) / (1024^3)  # GB
        
        # Test recursive implementation
        println("\nTesting recursive implementation memory usage...")
        GC.gc()
        before_memory = Sys.total_memory() - Sys.free_memory()
        
        output_path = joinpath(temp_dir, "recursive_memory.arrow")
        @time recursive_merge(refs, output_path, :pg_score, :target; 
                             reverse=[true, true], batch_size=50_000, group_size=5)
        
        GC.gc()
        after_memory = Sys.total_memory() - Sys.free_memory()
        recursive_usage = (after_memory - before_memory) / (1024^3)  # GB
        
        println("\nMemory usage estimates (GB):")
        println(@sprintf("  Current:   %.3f", current_usage))
        println(@sprintf("  Typed:     %.3f", typed_usage))
        println(@sprintf("  Recursive: %.3f", recursive_usage))
        
    finally
        rm(temp_dir, recursive=true)
    end
end

#==========================================================
File Handle Testing
==========================================================#

"""
Test behavior with increasing numbers of files to identify limits.
"""
function test_file_handle_limits(test_dir::String; max_files::Int=200)
    println("=== FILE HANDLE LIMITS TEST ===")
    
    all_files = [joinpath(test_dir, f) for f in readdir(test_dir) if endswith(f, ".arrow")]
    max_available = min(max_files, length(all_files))
    
    println("Testing up to $max_available files")
    
    # Test with increasing numbers of files
    test_counts = [10, 25, 50, 100, 150, 200, max_available]
    test_counts = filter(x -> x <= max_available, test_counts)
    
    for n_files in test_counts
        println("\n--- Testing $n_files files ---")
        
        file_paths = all_files[1:n_files]
        refs = create_file_references(file_paths)
        temp_dir = mktempdir(prefix="handle_test_")
        
        try
            # Test current implementation (likely to fail with many files)
            output_path = joinpath(temp_dir, "current.arrow")
            print("  Current implementation: ")
            
            try
                @time stream_sorted_merge(refs, output_path, :pg_score, :target; 
                                         reverse=[true, true], batch_size=10_000)
                println("SUCCESS")
            catch e
                println("FAILED - $e")
            end
            
            # Test recursive implementation (should handle more files)
            output_path = joinpath(temp_dir, "recursive.arrow")
            print("  Recursive implementation: ")
            
            try
                @time recursive_merge(refs, output_path, :pg_score, :target; 
                                     reverse=[true, true], batch_size=10_000, group_size=8)
                println("SUCCESS")
            catch e
                println("FAILED - $e")
            end
            
        finally
            rm(temp_dir, recursive=true)
        end
    end
end

#==========================================================
Main Test Suite
==========================================================#

"""
Run comprehensive test suite.
"""
function run_all_tests(test_dir::String)
    println("🚀 COMPREHENSIVE MERGE PERFORMANCE TEST SUITE")
    println("=" ^ 60)
    println("Test directory: $test_dir")
    println("Started at: $(now())")
    println()
    
    # Check that test directory exists and has files
    if !isdir(test_dir)
        error("Test directory does not exist: $test_dir")
    end
    
    arrow_files = [f for f in readdir(test_dir) if endswith(f, ".arrow")]
    if isempty(arrow_files)
        error("No .arrow files found in $test_dir")
    end
    
    println("Found $(length(arrow_files)) .arrow files")
    println()
    
    try
        # 1. Accuracy testing
        test_accuracy(test_dir, n_files=15, verbose=false)
        println()
        
        # 2. Performance benchmarking
        file_counts = if length(arrow_files) >= 100
            [10, 25, 50, 100]
        elseif length(arrow_files) >= 50
            [10, 25, 50]
        else
            [10, min(25, length(arrow_files))]
        end
        
        benchmark_merge_strategies(test_dir, file_counts=file_counts)
        println()
        
        # 3. Memory profiling (if enough files)
        if length(arrow_files) >= 30
            profile_memory_usage(test_dir, n_files=30)
            println()
        end
        
        # 4. File handle limits
        if length(arrow_files) >= 50
            test_file_handle_limits(test_dir, max_files=min(150, length(arrow_files)))
            println()
        end
        
        println("✅ All tests completed successfully!")
        
    catch e
        println("❌ Test suite failed with error: $e")
        rethrow(e)
    end
end

end # module MergePerformanceTests