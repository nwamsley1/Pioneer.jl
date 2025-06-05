# ScoringSearch Performance Testing Instructions

## Overview
This document contains instructions for testing performance optimizations in Pioneer.jl's ScoringSearch module, specifically for handling large numbers of Arrow files.

## Setup

### 1. Create the Optimized Utils File

Create file: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_optimized.jl`

```julia
# Optimized versions of critical functions in utils.jl

using DataFrames, Arrow, DataStructures

"""
    load_protein_groups_optimized(passing_pg_paths)

Optimized loading using single Arrow.Table call.
"""
function load_protein_groups_optimized(passing_pg_paths::Vector{String})
    # Filter valid paths
    valid_paths = filter(passing_pg_paths) do pg_path
        isfile(pg_path) && endswith(pg_path, ".arrow")
    end
    
    if isempty(valid_paths)
        return DataFrame()
    end
    
    # Load all at once - Arrow handles multiple files efficiently
    return DataFrame(Arrow.Table(valid_paths))
end

"""
    count_total_rows_optimized(file_paths)

Count rows without loading data into memory.
"""
function count_total_rows_optimized(file_paths::Vector{String})
    total_rows = 0
    
    for path in file_paths
        if isfile(path) && endswith(path, ".arrow")
            # Use Arrow metadata to get row count without loading data
            table = Arrow.Table(path)
            if !isempty(propertynames(table))
                col = first(propertynames(table))
                total_rows += length(table[col])
            end
        end
    end
    
    return total_rows
end

"""
    process_protein_groups_in_batches(pg_paths, process_func; batch_size=100000)

Process protein groups in memory-efficient batches.
"""
function process_protein_groups_in_batches(
    pg_paths::Vector{String},
    process_func::Function;
    batch_size::Int = 100000
)
    for path in pg_paths
        if !isfile(path) || !endswith(path, ".arrow")
            continue
        end
        
        # Process file in batches
        table = Arrow.Table(path)
        n_rows = length(table[first(propertynames(table))])
        
        for start_idx in 1:batch_size:n_rows
            end_idx = min(start_idx + batch_size - 1, n_rows)
            
            # Create batch DataFrame
            batch_df = DataFrame()
            for col in propertynames(table)
                batch_df[!, col] = table[col][start_idx:end_idx]
            end
            
            # Process batch
            process_func(batch_df)
        end
    end
end

"""
    sort_and_filter_parallel(file_paths, sort_col; filter_func=nothing)

Sort and optionally filter files in parallel.
"""
function sort_and_filter_parallel(
    file_paths::Vector{String},
    sort_col::Symbol;
    filter_func=nothing
)
    Threads.@threads for fpath in file_paths
        try
            df = DataFrame(Arrow.Table(fpath))
            
            # Apply filter if provided
            if !isnothing(filter_func)
                filter!(filter_func, df)
            end
            
            # Sort
            sort!(df, sort_col, rev=true)
            
            # Write back
            Arrow.write(fpath, df)
        catch e
            @error "Failed to process file $fpath" exception=e
        end
    end
end
```

### 2. Create the Performance Test Script

Create file: `test_scoring_performance.jl`

```julia
using Pkg
Pkg.activate(".")

# Load necessary packages
using DataFrames, Arrow, BenchmarkTools, Statistics
using Base.Threads: @threads

# Include the optimized functions
include("src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_optimized.jl")

# Configuration
const TEST_FOLDER = "/Users/nathanwamsley/Data/test_io/passing_proteins"  # Update this path
const MAX_FILES_TO_TEST = 100  # Limit for initial testing

"""
Get test files from the specified folder
"""
function get_test_files(folder::String, max_files::Int=MAX_FILES_TO_TEST)
    all_files = [f for f in readdir(folder, join=true) if endswith(f, ".arrow")]
    return all_files[1:min(length(all_files), max_files)]
end

"""
Phase 1: Test basic optimizations
"""
function test_phase1_optimizations(test_files::Vector{String})
    println("\n" * "="^60)
    println("PHASE 1: BASIC OPTIMIZATIONS TEST")
    println("="^60)
    println("Testing with $(length(test_files)) files")
    
    results = Dict{String, Any}()
    
    # Test 1: Loading optimization
    println("\n1. TESTING DATAFRAME LOADING")
    println("-" * 40)
    
    # Original method
    println("Original method (append! in loop):")
    GC.gc()
    mem_before = Base.gc_live_bytes()
    original_time = @elapsed begin
        df_original = DataFrame()
        for file in test_files
            append!(df_original, DataFrame(Arrow.Table(file)))
        end
    end
    mem_after_original = Base.gc_live_bytes()
    mem_used_original = (mem_after_original - mem_before) / 1024^2
    
    println("  Time: $(round(original_time, digits=3))s")
    println("  Memory: $(round(mem_used_original, digits=2)) MB")
    println("  Rows: $(nrow(df_original))")
    
    df_original = nothing
    GC.gc()
    
    # Optimized method
    println("\nOptimized method (single Arrow.Table):")
    mem_before = Base.gc_live_bytes()
    optimized_time = @elapsed begin
        df_optimized = load_protein_groups_optimized(test_files)
    end
    mem_after_optimized = Base.gc_live_bytes()
    mem_used_optimized = (mem_after_optimized - mem_before) / 1024^2
    
    println("  Time: $(round(optimized_time, digits=3))s")
    println("  Memory: $(round(mem_used_optimized, digits=2)) MB")
    println("  Rows: $(nrow(df_optimized))")
    
    speedup = round(original_time / optimized_time, digits=2)
    mem_reduction = round((1 - mem_used_optimized/mem_used_original) * 100, digits=1)
    
    println("\nImprovement:")
    println("  Speedup: $(speedup)x")
    println("  Memory reduction: $(mem_reduction)%")
    
    results["loading"] = Dict(
        "original_time" => original_time,
        "optimized_time" => optimized_time,
        "speedup" => speedup,
        "memory_reduction" => mem_reduction
    )
    
    # Test 2: Row counting optimization
    println("\n2. TESTING ROW COUNTING")
    println("-" * 40)
    
    # Original method
    println("Original method (load table):")
    original_count_time = @elapsed begin
        total_rows_original = 0
        for file in test_files
            table = Arrow.Table(file)
            total_rows_original += length(table[first(propertynames(table))])
        end
    end
    println("  Time: $(round(original_count_time, digits=3))s")
    println("  Total rows: $total_rows_original")
    
    # Optimized method
    println("\nOptimized method:")
    optimized_count_time = @elapsed begin
        total_rows_optimized = count_total_rows_optimized(test_files)
    end
    println("  Time: $(round(optimized_count_time, digits=3))s")
    println("  Total rows: $total_rows_optimized")
    
    count_speedup = round(original_count_time / optimized_count_time, digits=2)
    println("\nSpeedup: $(count_speedup)x")
    
    results["counting"] = Dict(
        "original_time" => original_count_time,
        "optimized_time" => optimized_count_time,
        "speedup" => count_speedup
    )
    
    # Test 3: Parallel processing
    if Threads.nthreads() > 1
        println("\n3. TESTING PARALLEL PROCESSING")
        println("-" * 40)
        println("Using $(Threads.nthreads()) threads")
        
        # Create temporary copies for testing
        test_copies = String[]
        for (i, file) in enumerate(test_files[1:min(10, length(test_files))])
            copy_path = tempname() * ".arrow"
            cp(file, copy_path)
            push!(test_copies, copy_path)
        end
        
        # Sequential processing
        println("\nSequential processing:")
        seq_time = @elapsed begin
            for file in test_copies
                df = DataFrame(Arrow.Table(file))
                sort!(df, names(df)[1], rev=true)
                Arrow.write(file, df)
            end
        end
        println("  Time: $(round(seq_time, digits=3))s")
        
        # Parallel processing
        println("\nParallel processing:")
        par_time = @elapsed begin
            sort_and_filter_parallel(test_copies, Symbol(names(DataFrame(Arrow.Table(test_copies[1])))[1]))
        end
        println("  Time: $(round(par_time, digits=3))s")
        
        parallel_speedup = round(seq_time / par_time, digits=2)
        println("\nSpeedup: $(parallel_speedup)x")
        
        # Cleanup
        for file in test_copies
            rm(file, force=true)
        end
        
        results["parallel"] = Dict(
            "sequential_time" => seq_time,
            "parallel_time" => par_time,
            "speedup" => parallel_speedup
        )
    else
        println("\n3. PARALLEL PROCESSING")
        println("  Skipped (only 1 thread available)")
    end
    
    # Test 4: Batch processing
    println("\n4. TESTING BATCH PROCESSING")
    println("-" * 40)
    
    total_rows_processed = 0
    batch_time = @elapsed begin
        process_protein_groups_in_batches(test_files[1:min(5, length(test_files))], batch_size=10000) do batch
            total_rows_processed += nrow(batch)
        end
    end
    
    println("  Processed $total_rows_processed rows in $(round(batch_time, digits=3))s")
    println("  Rate: $(round(total_rows_processed / batch_time, digits=0)) rows/second")
    
    results["batch_processing"] = Dict(
        "time" => batch_time,
        "rows" => total_rows_processed,
        "rate" => total_rows_processed / batch_time
    )
    
    return results
end

"""
Generate summary report
"""
function generate_report(results::Dict)
    println("\n" * "="^60)
    println("PERFORMANCE TEST SUMMARY")
    println("="^60)
    
    if haskey(results, "loading")
        println("\nDataFrame Loading:")
        println("  Original: $(round(results["loading"]["original_time"], digits=3))s")
        println("  Optimized: $(round(results["loading"]["optimized_time"], digits=3))s")
        println("  Speedup: $(results["loading"]["speedup"])x")
        println("  Memory reduction: $(results["loading"]["memory_reduction"])%")
    end
    
    if haskey(results, "counting")
        println("\nRow Counting:")
        println("  Original: $(round(results["counting"]["original_time"], digits=3))s")
        println("  Optimized: $(round(results["counting"]["optimized_time"], digits=3))s")
        println("  Speedup: $(results["counting"]["speedup"])x")
    end
    
    if haskey(results, "parallel")
        println("\nParallel Processing:")
        println("  Sequential: $(round(results["parallel"]["sequential_time"], digits=3))s")
        println("  Parallel: $(round(results["parallel"]["parallel_time"], digits=3))s")
        println("  Speedup: $(results["parallel"]["speedup"])x")
    end
    
    if haskey(results, "batch_processing")
        println("\nBatch Processing:")
        println("  Time: $(round(results["batch_processing"]["time"], digits=3))s")
        println("  Rate: $(round(results["batch_processing"]["rate"], digits=0)) rows/second")
    end
    
    println("\n" * "="^60)
    println("RECOMMENDATIONS")
    println("="^60)
    
    if haskey(results, "loading") && results["loading"]["speedup"] > 1.5
        println("✓ Use optimized loading ($(results["loading"]["speedup"])x speedup)")
    end
    
    if haskey(results, "parallel") && results["parallel"]["speedup"] > 1.5
        println("✓ Enable parallel processing ($(results["parallel"]["speedup"])x speedup)")
    end
    
    if haskey(results, "loading") && results["loading"]["memory_reduction"] > 20
        println("✓ Significant memory savings ($(results["loading"]["memory_reduction"])% reduction)")
    end
end

"""
Main test runner
"""
function run_performance_tests()
    println("ScoringSearch Performance Testing")
    println("Test folder: $TEST_FOLDER")
    
    # Check folder exists
    if !isdir(TEST_FOLDER)
        println("\nError: Test folder not found!")
        println("Please update TEST_FOLDER in the script.")
        return
    end
    
    # Get test files
    test_files = get_test_files(TEST_FOLDER)
    println("\nFound $(length(readdir(TEST_FOLDER, join=true))) total Arrow files")
    println("Testing with $(length(test_files)) files")
    
    if isempty(test_files)
        println("\nError: No Arrow files found in test folder!")
        return
    end
    
    # Run Phase 1 tests
    results = test_phase1_optimizations(test_files)
    
    # Generate report
    generate_report(results)
    
    return results
end

# Run the tests
run_performance_tests()
```

## Running the Tests

1. **On the target computer**, save both files in the appropriate locations
2. Update the `TEST_FOLDER` path in the test script
3. Run in Julia:
   ```julia
   include("test_scoring_performance.jl")
   ```

## Expected Output

The script will test:
1. **Loading optimization** - Compare append! vs single Arrow.Table load
2. **Row counting** - Compare loading full data vs metadata-only counting
3. **Parallel processing** - If multiple threads available
4. **Batch processing** - Memory-efficient processing of large files

You'll get a summary report with speedups and recommendations.

## Interpreting Results

- **Speedup > 2x**: Significant improvement, definitely implement
- **Memory reduction > 30%**: Important for large datasets
- **Parallel speedup**: Should be close to number of threads

## Next Steps

Based on results:
- If loading speedup is significant → Implement in production
- If parallel processing helps → Add to all file operations
- If memory reduction is large → Critical for your large datasets