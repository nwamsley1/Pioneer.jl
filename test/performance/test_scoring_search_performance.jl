using Pkg
Pkg.activate(".")

include("../../src/Pioneer.jl")
using .Pioneer
using BenchmarkTools
using Profile
using Arrow
using DataFrames
using Random

"""
Performance testing suite for ScoringSearch optimizations
"""

# Configuration for test dataset
const TEST_CONFIG_PATH = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/MtacYeastAlternating3M_1p.json"

# Load configuration to find output paths
function get_test_paths()
    config = JSON.parsefile(TEST_CONFIG_PATH)
    output_dir = config["global_settings"]["output"]["results_path"]
    temp_dir = joinpath(output_dir, "temp_data")
    
    return Dict(
        "output_dir" => output_dir,
        "temp_dir" => temp_dir,
        "passing_psms_folder" => joinpath(temp_dir, "passing_psms"),
        "passing_proteins_folder" => joinpath(temp_dir, "passing_proteins"),
        "second_pass_folder" => joinpath(temp_dir, "second_pass_psms")
    )
end

# Benchmark 1: Loading Arrow files
function benchmark_arrow_loading(folder_path::String, max_files::Int=10)
    files = [f for f in readdir(folder_path, join=true) if endswith(f, ".arrow")]
    files = files[1:min(length(files), max_files)]
    
    println("\n=== Arrow Loading Benchmark ($(length(files)) files) ===")
    
    # Method 1: Load all into memory
    print("Method 1 - Load all into memory: ")
    t1 = @benchmark begin
        dfs = []
        for file in $files
            push!(dfs, DataFrame(Arrow.Table(file)))
        end
    end samples=3
    display(t1)
    
    # Method 2: Streaming approach (to be implemented)
    print("\nMethod 2 - Arrow.Stream (lazy loading): ")
    t2 = @benchmark begin
        streams = []
        for file in $files
            open(file) do io
                push!(streams, Arrow.Stream(io))
            end
        end
    end samples=3
    display(t2)
    
    return t1, t2
end

# Benchmark 2: DataFrame append operations
function benchmark_dataframe_append(folder_path::String, max_files::Int=10)
    files = [f for f in readdir(folder_path, join=true) if endswith(f, ".arrow")]
    files = files[1:min(length(files), max_files)]
    
    println("\n=== DataFrame Append Benchmark ($(length(files)) files) ===")
    
    # Method 1: Repeated append!
    print("Method 1 - Repeated append!: ")
    t1 = @benchmark begin
        df = DataFrame()
        for file in $files
            append!(df, DataFrame(Arrow.Table(file)))
        end
    end samples=3
    display(t1)
    
    # Method 2: Pre-allocate and vcat
    print("\nMethod 2 - Pre-allocate and vcat: ")
    t2 = @benchmark begin
        dfs = [DataFrame(Arrow.Table(file)) for file in $files]
        df = vcat(dfs...)
    end samples=3
    display(t2)
    
    # Method 3: Single Arrow.Table with multiple files
    print("\nMethod 3 - Single Arrow.Table load: ")
    t3 = @benchmark begin
        df = DataFrame(Arrow.Table($files))
    end samples=3
    display(t3)
    
    return t1, t2, t3
end

# Benchmark 3: Sorting operations
function benchmark_sorting(df::DataFrame)
    println("\n=== Sorting Benchmark ($(nrow(df)) rows) ===")
    
    # Assuming there's a :prob column
    if :prob in names(df)
        print("QuickSort: ")
        t1 = @benchmark sort($df, :prob, rev=true, alg=QuickSort) samples=3
        display(t1)
        
        print("\nMergeSort: ")
        t2 = @benchmark sort($df, :prob, rev=true, alg=MergeSort) samples=3
        display(t2)
        
        print("\nDefault sort: ")
        t3 = @benchmark sort($df, :prob, rev=true) samples=3
        display(t3)
        
        return t1, t2, t3
    else
        println("No :prob column found for sorting benchmark")
        return nothing
    end
end

# Benchmark 4: File I/O patterns
function benchmark_file_io(file_path::String)
    println("\n=== File I/O Benchmark ===")
    
    # Method 1: Read-modify-write pattern
    print("Method 1 - Read, modify in memory, write: ")
    t1 = @benchmark begin
        df = DataFrame(Arrow.Table($file_path))
        df[!, :test_col] = rand(nrow(df))
        Pioneer.writeArrow(tempname() * ".arrow", df)
    end samples=3
    display(t1)
    
    # Method 2: Streaming transformation (to be implemented)
    # This would require implementing a streaming transformation approach
    
    return t1
end

# Benchmark 5: Merge operations
function benchmark_merge_operations(files::Vector{String})
    println("\n=== Merge Operations Benchmark ($(length(files)) files) ===")
    
    # Get sample of files
    sample_files = files[1:min(length(files), 5)]
    
    # Method 1: Current implementation (simplified)
    print("Method 1 - Load all tables in memory: ")
    t1 = @benchmark begin
        tables = [Arrow.Table(path) for path in $sample_files]
        # Simulate merge operation
    end samples=3
    display(t1)
    
    return t1
end

# Profile memory usage
function profile_memory_usage(func::Function, args...)
    println("\n=== Memory Profile ===")
    
    # Get baseline memory
    GC.gc()
    baseline = Base.gc_live_bytes()
    
    # Run function
    result = func(args...)
    
    # Get peak memory
    peak = Base.gc_live_bytes()
    
    println("Memory used: $(round((peak - baseline) / 1024^2, digits=2)) MB")
    
    return result
end

# Main benchmark runner
function run_benchmarks()
    paths = get_test_paths()
    
    println("Test paths:")
    for (k, v) in paths
        println("  $k: $v")
    end
    
    # Check if directories exist
    if !isdir(paths["passing_psms_folder"])
        println("\nError: passing_psms folder not found. Please run SearchDIA first.")
        return
    end
    
    # Run benchmarks
    println("\n" * "="^60)
    println("SCORING SEARCH PERFORMANCE BENCHMARKS")
    println("="^60)
    
    # 1. Arrow loading
    arrow_results = benchmark_arrow_loading(paths["passing_psms_folder"])
    
    # 2. DataFrame appends
    append_results = benchmark_dataframe_append(paths["passing_psms_folder"])
    
    # 3. Sorting (load a sample DataFrame)
    sample_files = [f for f in readdir(paths["passing_psms_folder"], join=true) if endswith(f, ".arrow")]
    if !isempty(sample_files)
        sample_df = DataFrame(Arrow.Table(sample_files[1]))
        sorting_results = benchmark_sorting(sample_df)
    end
    
    # 4. File I/O
    if !isempty(sample_files)
        io_results = benchmark_file_io(sample_files[1])
    end
    
    # 5. Merge operations
    if length(sample_files) > 1
        merge_results = benchmark_merge_operations(sample_files[1:min(5, length(sample_files))])
    end
    
    # Memory profiling example
    println("\n" * "="^60)
    println("MEMORY PROFILING")
    println("="^60)
    
    profile_memory_usage(benchmark_dataframe_append, paths["passing_psms_folder"], 5)
    
    println("\n" * "="^60)
    println("BENCHMARK COMPLETE")
    println("="^60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmarks()
end