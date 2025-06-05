# ScoringSearch Performance Optimization Plan

## Overview
This document outlines the performance optimization strategy for the ScoringSearch module when processing large datasets (hundreds of files).

## Identified Bottlenecks

### 1. Memory-Intensive Operations
- Loading entire Arrow files into DataFrames
- DataFrame append operations in loops
- Loading all files simultaneously for merging

### 2. I/O Inefficiencies
- Repeated reading and writing of the same files
- No streaming/chunked processing
- Synchronous file operations

### 3. Algorithmic Inefficiencies
- In-memory sorting of massive datasets
- No parallel processing where applicable
- Inefficient data structure conversions

## Optimization Strategies

### Strategy 1: Streaming Processing

Replace full-file loading with streaming approaches:

```julia
# Current approach (problematic)
df = DataFrame(Arrow.Table(file_path))

# Optimized approach
function process_arrow_stream(file_path::String, process_func::Function)
    open(file_path) do io
        stream = Arrow.Stream(io)
        for batch in stream
            process_func(batch)
        end
    end
end
```

### Strategy 2: Batch Processing for Large Operations

Implement batch processing for merge operations:

```julia
function merge_sorted_psms_scores_optimized(
    input_paths::Vector{String}, 
    output_path::String,
    prob_col::Symbol;
    batch_size = 1000000,
    max_open_files = 100
)
    # Process files in batches to limit memory usage
    # Use external merge sort algorithm
    # Write results incrementally
end
```

### Strategy 3: Parallel File Processing

Use multi-threading for independent file operations:

```julia
function sort_and_filter_quant_tables_parallel(
    psms_paths::Vector{String},
    prob_col::Symbol,
    best_traces::Set
)
    Threads.@threads for fpath in psms_paths
        process_single_file(fpath, prob_col, best_traces)
    end
end
```

### Strategy 4: Memory-Mapped Files

For very large files, use memory mapping:

```julia
function process_large_arrow_file(file_path::String)
    # Use memory mapping for large files
    mmap_data = Mmap.mmap(file_path)
    # Process data in chunks
end
```

### Strategy 5: Lazy Evaluation

Implement lazy evaluation for operations that don't need full data:

```julia
# Count rows without loading full table
function count_arrow_rows(file_path::String)
    metadata = Arrow.getmetadata(file_path)
    return metadata.num_rows
end
```

### Strategy 6: External Sorting

For massive datasets, implement external merge sort:

```julia
function external_sort_arrow_files(
    input_paths::Vector{String},
    output_path::String,
    sort_col::Symbol;
    chunk_size = 1000000
)
    # Sort chunks and write to temp files
    # Merge sorted chunks
end
```

## Implementation Priority

### High Priority (Quick Wins)
1. **Replace DataFrame appends with vcat or pre-allocation**
   - Impact: High
   - Effort: Low
   - Files: utils.jl (lines 1044-1048)

2. **Implement parallel file processing**
   - Impact: High
   - Effort: Medium
   - Files: utils.jl (sort_and_filter_quant_tables)

3. **Use Arrow.Table with multiple files**
   - Impact: Medium
   - Effort: Low
   - Files: score_psms.jl (load_psms_for_xgboost)

### Medium Priority
4. **Implement streaming for large operations**
   - Impact: High
   - Effort: High
   - Files: utils.jl (merge functions)

5. **Optimize sorting with external sort**
   - Impact: High for very large data
   - Effort: High
   - Files: utils.jl (sorting operations)

### Low Priority
6. **Memory-mapped file processing**
   - Impact: Medium
   - Effort: Medium
   - For edge cases with extremely large files

## Specific Optimizations

### 1. Optimize `get_protein_groups` (utils.jl)

```julia
# Replace lines 1044-1048
# Current:
all_protein_groups = DataFrame()
for pg_path in passing_pg_paths
    if isfile(pg_path) && endswith(pg_path, ".arrow")
        append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
    end
end

# Optimized:
valid_paths = filter(pg_path -> isfile(pg_path) && endswith(pg_path, ".arrow"), passing_pg_paths)
all_protein_groups = DataFrame(Arrow.Table(valid_paths))
```

### 2. Optimize `sort_and_filter_quant_tables` (utils.jl)

```julia
# Add parallel processing
Threads.@threads for fpath in second_pass_psms_paths
    # Process each file independently
end
```

### 3. Optimize `merge_sorted_psms_scores` (utils.jl)

Implement streaming merge:
- Don't load all tables at once
- Use iterators over Arrow batches
- Write output incrementally

### 4. Optimize `load_psms_for_xgboost` (score_psms.jl)

```julia
# Current (line 81):
return DataFrame(Tables.columntable(Arrow.Table(file_paths)))

# Optimized:
return DataFrame(Arrow.Table(file_paths))  # Arrow handles multiple files efficiently
```

## Testing Strategy

1. **Benchmark current implementation** using test suite
2. **Implement optimizations incrementally**
3. **Measure improvement after each change**
4. **Validate results match original implementation**

## Expected Performance Improvements

- **Memory usage**: 50-80% reduction for large datasets
- **Processing time**: 2-5x speedup depending on dataset size
- **Scalability**: Linear scaling with number of files (currently super-linear)

## Implementation Timeline

1. Week 1: Implement high-priority optimizations
2. Week 2: Test and validate changes
3. Week 3: Implement medium-priority optimizations
4. Week 4: Performance testing and documentation