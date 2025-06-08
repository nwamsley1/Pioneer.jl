# ScoringSearch Performance Diagnostics

This document tracks all performance diagnostics added to `utils.jl` for analyzing bottlenecks in large-scale searches. These changes can be reverted once optimization is complete.

## Overview

Added comprehensive timing, memory tracking, and throughput monitoring to identify performance bottlenecks when processing hundreds of files.

## Functions Modified

### 1. `get_best_traces` (lines 16-72)

**Added diagnostics:**
- Line 20-22: Start timing and initial memory tracking
- Line 31-35: Per-file timing and row count logging
- Line 66-69: Final elapsed time, memory usage, and trace count

**Original code to restore:**
```julia
function get_best_traces(
    second_pass_psms_paths::Vector{String},
    min_prob::Float32 = 0.75f0
)
    psms_trace_scores = Dictionary{
            @NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()

    for file_path in second_pass_psms_paths
        if splitext(file_path)[end] != ".arrow"
            continue
        end
        row_score = zero(Float32)
        psms_table = Arrow.Table(file_path)
        for i in range(1, length(psms_table[1]))
            # ... rest of loop
        end
    end
    # ... rest of function
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
end
```

### 2. `sort_and_filter_quant_tables` (lines 83-145)

**Added diagnostics:**
- Lines 90-93: Start timing and row counters
- Lines 101-110: File read timing and row tracking
- Lines 121-125: Filter timing and row retention
- Lines 128-130: Sort timing
- Lines 133-138: Write timing and per-file summary
- Lines 141-142: Final summary with retention rate

**Original code to restore:**
```julia
function sort_and_filter_quant_tables(
    second_pass_psms_paths::Vector{String},
    merged_quant_path::String,
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}
)
    #Remove if present 
    if isfile(merged_quant_path)
        rm(merged_quant_path)
    end
    #file_paths = [fpath for fpath in readdir(quant_psms_folder,join=true) if endswith(fpath,".arrow")]
    #Sort and filter each psm table 
    for fpath in second_pass_psms_paths
        psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        
        # ... filtering logic ...
        
        #Filter out unused traces 
        filter!(x->x.best_trace,psms_table)
        #Sort in descending order of probability
        sort!(psms_table, [prob_col, :target], rev = [true, true], alg=QuickSort)
        #write back
        writeArrow(fpath, psms_table)
    end
    return nothing
end
```

### 3. `merge_sorted_psms_scores` (lines 188-361)

**Added diagnostics:**
- Lines 194-196: Start timing and memory tracking
- Lines 241-246: Table loading timing and size reporting
- Line 270: Row counter initialization
- Line 277: Merge start timing
- Lines 319-330: Batch write timing, GC monitoring
- Lines 355-358: Final summary with throughput

**Original code to restore:**
```julia
function merge_sorted_psms_scores(
                    input_paths::Vector{String}, 
                    output_path::String,
                    prob_col::Symbol
                    ; N = 10000000
)
    # ... function body without @info statements ...
    tables = [Arrow.Table(path) for path in input_paths]
    # ... rest of original function ...
    #println("output_path $output_path")
    return nothing
end
```

### 4. `get_protein_groups` (lines 791-1234)

**Added diagnostics:**
- Lines 801-803: Start timing and memory tracking
- Lines 1013-1022: Peptide counting timing
- Lines 1038-1039: Peptide count completion
- Lines 1047-1051: PSM loading timing
- Lines 1079-1085: Protein inference timing
- Lines 1095-1096: First pass start
- Lines 1105-1132: Per-file first pass timing (dict creation, write)
- Lines 1141-1142: First pass completion
- Lines 1145-1146: Second pass start
- Lines 1162-1176: Per-file second pass timing and completion
- Lines 1229-1231: Final summary

**Original code to restore:**
Remove all @info statements and timing variables, keeping only the core logic.

### 5. `merge_sorted_protein_groups` (lines 1819-1977)

**Added diagnostics:**
- Lines 1825-1827: Start timing and memory tracking
- Lines 1869-1877: File loading and table size reporting  
- Lines 1898-1899: Row counter and merge timing
- Lines 1942-1945: Batch write progress
- Lines 1971-1974: Final summary

**Original code to restore:**
```julia
function merge_sorted_protein_groups(
    input_dir::String, 
    output_path::String,
    sort_key::Symbol;
    N = 1000000 
) #N -> batch size 
    # ... function body without timing/logging ...
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    tables = [Arrow.Table(path) for path in input_paths]
    # ... rest of original function ...
    return nothing
end
```

## Removing Diagnostics

To remove all diagnostics:

1. Delete all lines containing `@info "[PERF]`
2. Remove all timing variables:
   - `start_time = time()`
   - `initial_memory = Base.gc_live_bytes() / 1024^2`
   - `*_start = time()`
   - `*_time = time() - *_start`
   - `elapsed = time() - start_time`
   - `final_memory = Base.gc_live_bytes() / 1024^2`
   - `memory_used = final_memory - initial_memory`

3. Remove all tracking variables:
   - `total_rows_processed`
   - `total_rows_kept`
   - `rows_written`
   - `n_rows`, `initial_rows`, `filtered_rows`
   - `file_start`, `read_start`, `filter_start`, `sort_start`, `write_start`
   - `table_sizes`, `total_rows`

4. Remove GC monitoring block in `merge_sorted_psms_scores` (lines 323-330)

5. Change enumerated loops back to regular loops where only used for diagnostics:
   - `for (file_idx, file_path) in enumerate(...)` → `for file_path in ...`

## Summary of Performance Metrics Collected

- **Timing**: Operation elapsed time at multiple granularities
- **Memory**: Initial/final usage, GC triggers
- **I/O**: Read/write operation timing
- **Data Volume**: Row counts, file counts, retention rates
- **Throughput**: Rows/second processing rates

These diagnostics will help identify:
- File I/O bottlenecks
- Memory allocation hotspots
- GC pressure points
- Inefficient sorting/filtering
- Processing rate degradation with scale