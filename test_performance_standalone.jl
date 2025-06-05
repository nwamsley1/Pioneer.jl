using Pkg
Pkg.activate(".")

# Load necessary packages
using DataFrames, Arrow, BenchmarkTools, Statistics
using Base.Threads: @threads

# Configuration for local test
const TEST_FOLDER = "/Users/nathanwamsley/Data/test_io/passing_proteins"
const MAX_FILES_TO_TEST = 100  # Start with 100 files for testing

# Standalone optimized functions (no Pioneer dependencies)

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
    println("Julia threads available: $(Threads.nthreads())")
    
    results = Dict{String, Any}()
    
    # Test 1: Loading optimization
    println("\n1. TESTING DATAFRAME LOADING")
    println("-"^40)
    
    # Warm up
    println("Warming up...")
    df_warmup = DataFrame(Arrow.Table(test_files[1]))
    df_warmup = nothing
    
    # Original method
    println("\nOriginal method (append! in loop):")
    GC.gc()
    
    original_bench = @benchmark begin
        df_original = DataFrame()
        for file in $test_files
            append!(df_original, DataFrame(Arrow.Table(file)))
        end
        df_original
    end samples=3 seconds=30
    
    # Get one sample for row count
    df_sample = DataFrame()
    for file in test_files
        append!(df_sample, DataFrame(Arrow.Table(file)))
    end
    n_rows = nrow(df_sample)
    df_sample = nothing
    GC.gc()
    
    println("  Time: $(round(median(original_bench.times)/1e9, digits=3))s (median of $(length(original_bench.times)) samples)")
    println("  Memory: $(round(original_bench.memory/1024^2, digits=2)) MB")
    println("  Allocations: $(original_bench.allocs)")
    println("  Rows: $n_rows")
    
    # Optimized method
    println("\nOptimized method (single Arrow.Table):")
    GC.gc()
    
    optimized_bench = @benchmark begin
        load_protein_groups_optimized($test_files)
    end samples=3 seconds=30
    
    println("  Time: $(round(median(optimized_bench.times)/1e9, digits=3))s (median of $(length(optimized_bench.times)) samples)")
    println("  Memory: $(round(optimized_bench.memory/1024^2, digits=2)) MB")
    println("  Allocations: $(optimized_bench.allocs)")
    
    speedup = round(median(original_bench.times) / median(optimized_bench.times), digits=2)
    mem_reduction = round((1 - optimized_bench.memory/original_bench.memory) * 100, digits=1)
    
    println("\nImprovement:")
    println("  Speedup: $(speedup)x")
    println("  Memory reduction: $(mem_reduction)%")
    println("  Allocation reduction: $(round((1 - optimized_bench.allocs/original_bench.allocs) * 100, digits=1))%")
    
    results["loading"] = Dict(
        "original_time" => median(original_bench.times)/1e9,
        "optimized_time" => median(optimized_bench.times)/1e9,
        "speedup" => speedup,
        "memory_reduction" => mem_reduction
    )
    
    # Test 2: Row counting optimization
    println("\n2. TESTING ROW COUNTING")
    println("-"^40)
    
    # Original method
    println("Original method (load table):")
    original_count_bench = @benchmark begin
        total_rows = 0
        for file in $test_files
            table = Arrow.Table(file)
            total_rows += length(table[first(propertynames(table))])
        end
        total_rows
    end samples=3
    
    # Optimized method
    println("\nOptimized method:")
    optimized_count_bench = @benchmark begin
        count_total_rows_optimized($test_files)
    end samples=3
    
    println("Original:")
    println("  Time: $(round(median(original_count_bench.times)/1e9, digits=3))s")
    println("  Memory: $(round(original_count_bench.memory/1024^2, digits=2)) MB")
    
    println("\nOptimized:")
    println("  Time: $(round(median(optimized_count_bench.times)/1e9, digits=3))s")
    println("  Memory: $(round(optimized_count_bench.memory/1024^2, digits=2)) MB")
    
    count_speedup = round(median(original_count_bench.times) / median(optimized_count_bench.times), digits=2)
    println("\nSpeedup: $(count_speedup)x")
    
    results["counting"] = Dict(
        "original_time" => median(original_count_bench.times)/1e9,
        "optimized_time" => median(optimized_count_bench.times)/1e9,
        "speedup" => count_speedup
    )
    
    # Test 3: Parallel processing
    if Threads.nthreads() > 1
        println("\n3. TESTING PARALLEL PROCESSING")
        println("-"^40)
        println("Using $(Threads.nthreads()) threads")
        
        # Create temporary copies for testing (just 10 files)
        test_copies = String[]
        for (i, file) in enumerate(test_files[1:min(10, length(test_files))])
            copy_path = tempname() * ".arrow"
            cp(file, copy_path)
            push!(test_copies, copy_path)
        end
        
        # Get first column name for sorting
        first_col = Symbol(names(DataFrame(Arrow.Table(test_copies[1])))[1])
        
        # Sequential processing
        println("\nSequential processing:")
        seq_bench = @benchmark begin
            for file in $test_copies
                df = DataFrame(Arrow.Table(file))
                sort!(df, $first_col, rev=true)
                Arrow.write(file, df)
            end
        end samples=3
        
        # Parallel processing
        println("\nParallel processing:")
        par_bench = @benchmark begin
            @threads for file in $test_copies
                df = DataFrame(Arrow.Table(file))
                sort!(df, $first_col, rev=true)
                Arrow.write(file, df)
            end
        end samples=3
        
        println("Sequential: $(round(median(seq_bench.times)/1e9, digits=3))s")
        println("Parallel: $(round(median(par_bench.times)/1e9, digits=3))s")
        
        parallel_speedup = round(median(seq_bench.times) / median(par_bench.times), digits=2)
        println("\nSpeedup: $(parallel_speedup)x")
        
        # Cleanup
        for file in test_copies
            rm(file, force=true)
        end
        
        results["parallel"] = Dict(
            "sequential_time" => median(seq_bench.times)/1e9,
            "parallel_time" => median(par_bench.times)/1e9,
            "speedup" => parallel_speedup
        )
    else
        println("\n3. PARALLEL PROCESSING")
        println("  Skipped (only 1 thread available)")
    end
    
    # Test 4: Batch processing
    println("\n4. TESTING BATCH PROCESSING")
    println("-"^40)
    
    test_files_batch = test_files[1:min(5, length(test_files))]
    total_rows_processed = 0
    
    batch_bench = @benchmark begin
        rows = 0
        process_protein_groups_in_batches($test_files_batch, batch_size=10000) do batch
            rows += nrow(batch)
        end
        rows
    end samples=3
    
    # Get actual row count
    process_protein_groups_in_batches(test_files_batch, batch_size=10000) do batch
        total_rows_processed += nrow(batch)
    end
    
    batch_time = median(batch_bench.times)/1e9
    println("  Processed $total_rows_processed rows in $(round(batch_time, digits=3))s")
    println("  Rate: $(round(total_rows_processed / batch_time, digits=0)) rows/second")
    println("  Memory: $(round(batch_bench.memory/1024^2, digits=2)) MB")
    
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
        println("  Efficiency: $(round(results["parallel"]["speedup"] / Threads.nthreads() * 100, digits=1))%")
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
        println("✓ Implement optimized loading ($(results["loading"]["speedup"])x speedup)")
    end
    
    if haskey(results, "parallel") && results["parallel"]["speedup"] > 1.5
        println("✓ Enable parallel processing ($(results["parallel"]["speedup"])x speedup)")
    end
    
    if haskey(results, "loading") && results["loading"]["memory_reduction"] > 20
        println("✓ Significant memory savings ($(results["loading"]["memory_reduction"])% reduction)")
    end
    
    # Estimate impact for large datasets
    if haskey(results, "loading")
        println("\nProjected impact for 1000 files:")
        orig_time_1000 = results["loading"]["original_time"] * 10
        opt_time_1000 = results["loading"]["optimized_time"] * 10
        println("  Original method: ~$(round(orig_time_1000/60, digits=1)) minutes")
        println("  Optimized method: ~$(round(opt_time_1000/60, digits=1)) minutes")
        println("  Time saved: ~$(round((orig_time_1000 - opt_time_1000)/60, digits=1)) minutes")
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
        return
    end
    
    # Get test files
    all_arrow_files = [f for f in readdir(TEST_FOLDER, join=true) if endswith(f, ".arrow")]
    test_files = get_test_files(TEST_FOLDER)
    
    println("\nFound $(length(all_arrow_files)) total Arrow files in folder")
    println("Testing with $(length(test_files)) files")
    
    if isempty(test_files)
        println("\nError: No Arrow files found in test folder!")
        return
    end
    
    # Check file sizes
    total_size = sum(filesize(f) for f in test_files) / 1024^2
    println("Total size of test files: $(round(total_size, digits=1)) MB")
    
    # Run Phase 1 tests
    results = test_phase1_optimizations(test_files)
    
    # Generate report
    generate_report(results)
    
    return results
end

# Run the tests
println("Starting performance tests...")
println("This may take a few minutes depending on file sizes...")
results = run_performance_tests()