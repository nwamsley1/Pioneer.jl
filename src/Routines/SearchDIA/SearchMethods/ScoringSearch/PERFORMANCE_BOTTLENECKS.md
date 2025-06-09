# ScoringSearch Performance Bottlenecks Analysis

## Overview

The ScoringSearch module processes large volumes of PSM (Peptide-Spectrum Match) and protein group data, performing multiple file I/O operations, sorting, and merging. With hundreds of files in large-scale experiments, several performance bottlenecks can emerge.

## Major Bottlenecks

### 1. File I/O Operations

#### Arrow File Reading
- **Issue**: Each file is read entirely into memory as a DataFrame
- **Location**: `sort_and_filter_quant_tables`, `merge_sorted_psms_scores`
- **Impact**: Memory spikes, GC pressure
- **Current pattern**:
  ```julia
  psms_table = DataFrame(Tables.columntable(Arrow.Table(fpath)))
  ```

**Optimization Opportunities**:
- Use Arrow.Stream for iterative processing
- Memory-map large files
- Process in chunks without full materialization

#### Arrow File Writing
- **Issue**: Frequent file rewrites in place
- **Location**: Throughout utils.jl
- **Impact**: Disk I/O bottleneck, temporary file creation
- **Current pattern**:
  ```julia
  writeArrow(fpath, psms_table)  # Overwrites original
  ```

**Optimization Opportunities**:
- Write to new files, then atomic rename
- Batch writes when possible
- Use Arrow IPC streaming format for append operations

### 2. Memory Management

#### DataFrame Materialization
- **Issue**: Converting Arrow tables to DataFrames forces full memory allocation
- **Location**: All sorting and filtering operations
- **Impact**: 2-3x memory usage vs streaming
- **Example**:
  ```julia
  DataFrame(Tables.columntable(Arrow.Table(fpath)))
  ```

**Optimization Opportunities**:
- Use lazy evaluation with Arrow tables directly
- Implement streaming transformations
- Process subsets of columns when full row not needed

#### Heap-based Merging Memory
- **Issue**: Pre-allocated batch sizes may be suboptimal
- **Location**: `merge_sorted_psms_scores`, `merge_sorted_protein_groups`
- **Current**: Fixed batch size N = 10,000,000
- **Impact**: Memory allocation/deallocation cycles

**Optimization Opportunities**:
- Dynamic batch sizing based on available memory
- Reuse buffers between batches
- Implement zero-copy merging where possible

### 3. Sorting Operations

#### In-Memory Sorting
- **Issue**: Full table sorts require all data in memory
- **Location**: `sort_quant_tables`, `sort_protein_tables`
- **Pattern**:
  ```julia
  sort!(psms_table, [prob_col, :target], rev = [true, true], alg=QuickSort)
  ```

**Optimization Opportunities**:
- External sorting for files > memory threshold
- Parallel sorting algorithms
- Maintain sorted order during initial write

#### Repeated Sorting
- **Issue**: Files may be sorted multiple times across pipeline stages
- **Impact**: Redundant CPU cycles

**Optimization Opportunities**:
- Track sort status in metadata
- Preserve sort order through pipeline
- Use sorted merge algorithms

### 4. Protein Group Processing

#### Dictionary Operations
- **Issue**: Large dictionaries for protein inference
- **Location**: `getProteinGroupsDict`
- **Impact**: Hash table resizing, collision handling

**Optimization Opportunities**:
- Pre-size dictionaries based on expected entries
- Use more efficient key types (integers vs strings)
- Consider sorted arrays for read-heavy operations

#### Two-Pass Processing
- **Issue**: Files read twice for global scoring
- **Location**: `get_protein_groups`
- **Impact**: Doubles I/O time

**Optimization Opportunities**:
- Cache intermediate results
- Combine passes where possible
- Use memory-mapped files for second pass

### 5. Thread Contention

#### Shared Data Structures
- **Issue**: Global dictionaries accessed by multiple threads
- **Impact**: Lock contention, false sharing

**Optimization Opportunities**:
- Thread-local accumulation with final merge
- Lock-free data structures
- Partition work to minimize sharing

### 6. GC Pressure Points

#### Temporary Allocations
- **Issue**: Many temporary DataFrames and arrays
- **Location**: Throughout filtering and transformation
- **Impact**: Frequent GC pauses

**Optimization Opportunities**:
- Reuse buffers
- In-place operations where possible
- Manual memory management for hot paths

## Recommended Optimizations

### High Priority

1. **Streaming Arrow Operations**
   ```julia
   # Instead of:
   df = DataFrame(Arrow.Table(path))
   
   # Use:
   Arrow.Stream(path) do stream
       for batch in stream
           # Process batch
       end
   end
   ```

2. **Lazy Column Selection**
   ```julia
   # Only load needed columns
   table = Arrow.Table(path, columns=[:prob, :target, :precursor_idx])
   ```

3. **Memory-Mapped Files**
   ```julia
   # For large read-only files
   table = Arrow.Table(path; mmap=true)
   ```

### Medium Priority

1. **Batch Processing Pipeline**
   - Process files in batches to control memory
   - Implement backpressure mechanisms
   - Use producer-consumer patterns

2. **Optimize Data Types**
   - Use categorical arrays for repeated strings
   - Compress integer columns
   - Use bit arrays for boolean flags

### Low Priority

1. **Parallel I/O**
   - Async file operations
   - Parallel file reading with thread pool
   - I/O scheduling optimization

2. **Custom Serialization**
   - Binary format for intermediate results
   - Compression for temporary files
   - Delta encoding for sorted data

## Benchmarking Strategy

1. **Profile Key Operations**
   ```julia
   @time @allocated merge_sorted_psms_scores(...)
   ```

2. **Monitor Memory Usage**
   ```julia
   GC.gc()
   memory_before = Base.gc_live_bytes()
   # Operation
   memory_after = Base.gc_live_bytes()
   ```

3. **Track I/O Patterns**
   - File open/close frequency
   - Read/write sizes
   - Seek patterns

## Platform-Specific Considerations

### Windows
- File locking more restrictive (addressed with safeRm)
- Different memory allocation patterns
- Antivirus impact on I/O

### Linux
- Better memory-mapped file support
- More efficient file caching
- Different GC behavior

### Memory Constraints
- Adjust batch sizes based on available RAM
- Implement out-of-core algorithms for large datasets
- Consider memory pressure indicators

## Future Improvements

1. **Arrow.jl 3.0 Features**
   - Better streaming support
   - Lazy evaluation improvements
   - Native array views

2. **Julia 1.12+ Optimizations**
   - Better GC for large allocations
   - Improved thread scheduling
   - Native memory pools

3. **Algorithmic Improvements**
   - Approximate algorithms for large scales
   - Sampling strategies
   - Incremental processing