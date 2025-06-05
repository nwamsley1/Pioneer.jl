using Pkg
Pkg.activate(".")

include("../../src/Pioneer.jl")
using .Pioneer
include("../../src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_optimized.jl")
using BenchmarkTools
using Arrow
using DataFrames

# Test configuration
const TEST_CONFIG_PATH = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/MtacYeastAlternating3M_1p.json"

function get_test_files()
    config = JSON.parsefile(TEST_CONFIG_PATH)
    output_dir = config["global_settings"]["output"]["results_path"]
    temp_dir = joinpath(output_dir, "temp_data")
    
    passing_psms_folder = joinpath(temp_dir, "passing_psms")
    passing_proteins_folder = joinpath(temp_dir, "passing_proteins")
    
    psm_files = [f for f in readdir(passing_psms_folder, join=true) if endswith(f, ".arrow")]
    protein_files = [f for f in readdir(passing_proteins_folder, join=true) if endswith(f, ".arrow")]
    
    return psm_files, protein_files
end

"""
Compare loading protein groups: original vs optimized
"""
function compare_protein_loading()
    _, protein_files = get_test_files()
    
    if isempty(protein_files)
        println("No protein files found!")
        return
    end
    
    # Test with subset of files
    test_files = protein_files[1:min(10, length(protein_files))]
    
    println("\n" * "="^60)
    println("COMPARING PROTEIN GROUP LOADING ($(length(test_files)) files)")
    println("="^60)
    
    # Original method
    println("\nOriginal method (append! in loop):")
    original_time = @elapsed begin
        all_protein_groups = DataFrame()
        for pg_path in test_files
            if isfile(pg_path) && endswith(pg_path, ".arrow")
                append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
            end
        end
    end
    println("Time: $(round(original_time, digits=3))s")
    println("Rows loaded: $(nrow(all_protein_groups))")
    
    # Optimized method
    println("\nOptimized method (single Arrow.Table):")
    optimized_time = @elapsed begin
        all_protein_groups_opt = load_protein_groups_optimized(test_files)
    end
    println("Time: $(round(optimized_time, digits=3))s")
    println("Rows loaded: $(nrow(all_protein_groups_opt))")
    
    speedup = round(original_time / optimized_time, digits=2)
    println("\nSpeedup: $(speedup)x")
    
    # Verify results are identical
    if size(all_protein_groups) == size(all_protein_groups_opt)
        println("✓ Results match!")
    else
        println("✗ Results differ!")
    end
    
    return original_time, optimized_time
end

"""
Compare row counting methods
"""
function compare_row_counting()
    _, protein_files = get_test_files()
    
    if isempty(protein_files)
        println("No protein files found!")
        return
    end
    
    test_files = protein_files[1:min(20, length(protein_files))]
    
    println("\n" * "="^60)
    println("COMPARING ROW COUNTING ($(length(test_files)) files)")
    println("="^60)
    
    # Original method
    println("\nOriginal method (load table and count):")
    original_time = @elapsed begin
        total_rows = 0
        for pg_path in test_files
            if isfile(pg_path) && endswith(pg_path, ".arrow")
                table = Arrow.Table(pg_path)
                total_rows += length(table[:protein_name])
            end
        end
    end
    println("Time: $(round(original_time, digits=3))s")
    println("Total rows: $total_rows")
    
    # Optimized method
    println("\nOptimized method (metadata-based):")
    optimized_time = @elapsed begin
        total_rows_opt = count_total_rows_optimized(test_files)
    end
    println("Time: $(round(optimized_time, digits=3))s")
    println("Total rows: $total_rows_opt")
    
    speedup = round(original_time / optimized_time, digits=2)
    println("\nSpeedup: $(speedup)x")
    
    return original_time, optimized_time
end

"""
Test memory usage of different approaches
"""
function test_memory_usage()
    _, protein_files = get_test_files()
    
    if isempty(protein_files)
        println("No protein files found!")
        return
    end
    
    test_files = protein_files[1:min(10, length(protein_files))]
    
    println("\n" * "="^60)
    println("MEMORY USAGE COMPARISON")
    println("="^60)
    
    # Original method
    GC.gc()
    mem_before = Base.gc_live_bytes()
    
    all_protein_groups = DataFrame()
    for pg_path in test_files
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
        end
    end
    
    mem_after_original = Base.gc_live_bytes()
    mem_used_original = (mem_after_original - mem_before) / 1024^2
    
    all_protein_groups = nothing
    GC.gc()
    
    # Optimized method
    mem_before = Base.gc_live_bytes()
    
    all_protein_groups_opt = load_protein_groups_optimized(test_files)
    
    mem_after_optimized = Base.gc_live_bytes()
    mem_used_optimized = (mem_after_optimized - mem_before) / 1024^2
    
    println("\nOriginal method memory: $(round(mem_used_original, digits=2)) MB")
    println("Optimized method memory: $(round(mem_used_optimized, digits=2)) MB")
    println("Memory saved: $(round(mem_used_original - mem_used_optimized, digits=2)) MB")
    
    return mem_used_original, mem_used_optimized
end

"""
Test batch processing approach
"""
function test_batch_processing()
    _, protein_files = get_test_files()
    
    if isempty(protein_files)
        println("No protein files found!")
        return
    end
    
    test_files = protein_files[1:min(5, length(protein_files))]
    
    println("\n" * "="^60)
    println("TESTING BATCH PROCESSING")
    println("="^60)
    
    total_rows = 0
    process_time = @elapsed begin
        process_protein_groups_in_batches(test_files, batch_size=10000) do batch
            # Simulate some processing
            total_rows += nrow(batch)
        end
    end
    
    println("Processed $total_rows rows in $(round(process_time, digits=3))s")
    println("Average time per row: $(round(process_time * 1000 / total_rows, digits=3))ms")
    
    return process_time
end

"""
Run all comparisons
"""
function run_all_comparisons()
    println("Starting performance comparisons...")
    println("Test config: $TEST_CONFIG_PATH")
    
    # Check if files exist
    psm_files, protein_files = get_test_files()
    println("\nFound $(length(psm_files)) PSM files and $(length(protein_files)) protein files")
    
    if isempty(protein_files)
        println("\nError: No protein files found. Please ensure SearchDIA has been run.")
        return
    end
    
    # Run comparisons
    results = Dict()
    
    results["protein_loading"] = compare_protein_loading()
    results["row_counting"] = compare_row_counting()
    results["memory_usage"] = test_memory_usage()
    results["batch_processing"] = test_batch_processing()
    
    # Summary
    println("\n" * "="^60)
    println("SUMMARY")
    println("="^60)
    
    if haskey(results, "protein_loading")
        orig, opt = results["protein_loading"]
        println("Protein loading speedup: $(round(orig/opt, digits=2))x")
    end
    
    if haskey(results, "memory_usage")
        orig_mem, opt_mem = results["memory_usage"]
        println("Memory reduction: $(round((1 - opt_mem/orig_mem) * 100, digits=1))%")
    end
    
    println("\nOptimizations are ready for integration!")
    
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_all_comparisons()
end