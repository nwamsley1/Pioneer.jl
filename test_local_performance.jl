using Pkg
Pkg.activate(".")

# Load necessary packages
using DataFrames, Arrow, BenchmarkTools, Statistics
using Base.Threads: @threads

# Include the optimized functions from the existing file
include("src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_optimized.jl")

# Configuration for local test
const TEST_FOLDER = "/Users/nathanwamsley/Data/test_io/passing_proteins"
const MAX_FILES_TO_TEST = 100  # Start with 100 files for testing

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
        seq_time = @elapsed begin
            for file in test_copies
                df = DataFrame(Arrow.Table(file))
                sort!(df, first_col, rev=true)
                Arrow.write(file, df)
            end
        end
        println("  Time: $(round(seq_time, digits=3))s")
        
        # Parallel processing
        println("\nParallel processing:")
        par_time = @elapsed begin
            @threads for file in test_copies
                df = DataFrame(Arrow.Table(file))
                sort!(df, first_col, rev=true)
                Arrow.write(file, df)
            end
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
        return
    end
    
    # Get test files
    all_arrow_files = [f for f in readdir(TEST_FOLDER, join=true) if endswith(f, ".arrow")]
    test_files = get_test_files(TEST_FOLDER)
    
    println("\nFound $(length(all_arrow_files)) total Arrow files")
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
results = run_performance_tests()