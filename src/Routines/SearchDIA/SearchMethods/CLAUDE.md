# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with SearchMethods in the SearchDIA pipeline.

## SearchMethods Overview

SearchMethods implements a 9-stage sequential pipeline where each method performs a specific aspect of DIA analysis. All methods follow a common interface pattern and share state through SearchContext.

## Interface Pattern

Every search method implements this exact interface:

```julia
abstract type SearchMethod end
abstract type SearchParameters end  
abstract type SearchResults end

# Required interface methods (exactly 6)
getSearchMethodName(method::SearchMethod) -> String
getParameters(method::SearchMethod) -> SearchParameters
getResults(method::SearchMethod) -> SearchResults
performSearch!(method::SearchMethod) -> Nothing
summarizeResults!(method::SearchMethod) -> Nothing
processChunk!(method::SearchMethod, batch_id::Int, thread_id::Int) -> Nothing
```

## Pipeline Architecture

### Stage 1: Parameter Optimization
**ParameterTuningSearch** - Fragment mass tolerance optimization
- **Purpose**: Establishes optimal fragment mass tolerance using grid search
- **Algorithm**: Tests multiple tolerance values, selects based on PSM count/quality
- **Key Output**: Calibrated fragment tolerance for downstream methods
- **Performance**: Fast, minimal memory usage

**NceTuningSearch** - Collision energy calibration  
- **Purpose**: Builds instrument-specific collision energy models
- **Algorithm**: Linear regression on CE vs observed intensity patterns
- **Key Output**: CE prediction models stored in SearchContext
- **Performance**: Moderate CPU, stores spline models

**QuadTuningSearch** - Quadrupole transmission modeling
- **Purpose**: Models quadrupole isolation efficiency across m/z range
- **Algorithm**: Isotope ratio analysis with spline fitting
- **Key Output**: Quadrupole transmission correction functions
- **Performance**: CPU-intensive spline fitting

**HuberTuningSearch** - Robust loss parameter optimization
- **Purpose**: Optimizes Huber loss parameters for robust quantification
- **Algorithm**: Cross-validation on PSM scoring performance
- **Key Output**: Huber delta parameters for downstream scoring
- **Performance**: Iterative optimization, moderate memory

### Stage 2: PSM Identification
**FirstPassSearch** - Initial PSM discovery
- **Purpose**: First-pass peptide identification with RT calibration
- **Algorithm**: Fragment matching + RT model building
- **Key Output**: Initial PSM set + RT conversion models
- **Memory**: Loads full fragment index, high memory usage
- **Threading**: Parallel chunk processing by scan ranges

**SecondPassSearch** - Comprehensive PSM identification
- **Purpose**: Exhaustive search with calibrated parameters
- **Algorithm**: Full library search using optimized tolerances/models
- **Key Output**: Complete PSM catalog for downstream analysis
- **Memory**: Highest memory usage - full library + data
- **Threading**: Most CPU-intensive stage

### Stage 3: Scoring & FDR Control
**ScoringSearch** - ML training and FDR control
- **Purpose**: EvoTrees/XGBoost training, PSM rescoring, protein grouping
- **Algorithm**: Cross-validation training + FDR filtering + protein inference
- **Key Output**: High-confidence PSMs + protein groups
- **Memory**: Stores all features for ML training
- **Threading**: Parallel model training per CV fold

### Stage 4: Quantification
**IntegrateChromatogramSearch** - Peak integration
- **Purpose**: Chromatographic peak integration and quantification
- **Algorithm**: Gaussian/spline fitting on extracted ion chromatograms
- **Key Output**: Integrated intensities for quantification
- **Memory**: Chromatogram storage, moderate usage
- **Performance**: I/O intensive for chromatogram extraction

**MaxLFQSearch** - Label-free quantification
- **Purpose**: Cross-run normalization and final quantification
- **Algorithm**: MaxLFQ algorithm with missing value handling
- **Key Output**: Final protein/peptide quantification tables
- **Memory**: Full quantification matrix in memory
- **Performance**: Linear algebra operations, CPU-intensive

## Common Development Patterns

### Method Implementation Template
```julia
# 1. Define types
struct MySearchMethod <: SearchMethod
    search_context::SearchContext
end

struct MySearchParameters <: SearchParameters
    # method-specific parameters
    param1::Float64
    param2::Int
    # Always include these common fields:
    output_folder::String
    temp_folder::String
end

struct MySearchResults <: SearchResults
    # method-specific results
    my_result::Vector{Float64}
end

# 2. Implement interface
getSearchMethodName(method::MySearchMethod) = "MySearch"

function getParameters(method::MySearchMethod)
    # Extract from PioneerParameters
    params = getParams(getSearchContext(method))
    return MySearchParameters(
        params.my_param1,
        params.my_param2,
        params.output_folder,
        params.temp_folder
    )
end

function performSearch!(method::MySearchMethod)
    # Main search logic - coordinate parallel processing
    thread_tasks = partitionThreadTasks(...)
    @threads for (batch_id, thread_id) in thread_tasks
        processChunk!(method, batch_id, thread_id)
    end
end

function processChunk!(method::MySearchMethod, batch_id::Int, thread_id::Int)
    # Thread-safe chunk processing
    # Access data via: getThreadSpecData(method, thread_id)
    # Store results in thread-local storage
end

function summarizeResults!(method::MySearchMethod)
    # Combine thread results + write outputs
    # Generate QC plots
    # Update SearchContext with results
end
```

### Parameter Extraction Pattern
```julia
# Standard parameter extraction from PioneerParameters
function getParameters(method::MySearchMethod)
    params = getParams(getSearchContext(method))
    return MySearchParameters(
        Float64(params.optimization.my_method.param1),
        Int64(params.optimization.my_method.param2),
        String(params.output_folder),
        String(params.temp_folder)
    )
end
```

### Threading Pattern
```julia
function processChunk!(method::MySearchMethod, batch_id::Int, thread_id::Int)
    # Get thread-specific data access
    spec_data = getThreadSpecData(method, thread_id)
    results = Vector{MyResult}()
    
    # Process assigned scans
    for scan_idx in getAssignedScans(batch_id)
        # Perform computation
        result = processOneScan(spec_data, scan_idx)
        push!(results, result)
    end
    
    # Store in thread-local storage (thread-safe)
    setThreadResults!(method, thread_id, results)
end
```

## Performance Optimization Guidelines

### Memory Management
1. **Avoid loading full datasets** - Use iterators and chunking
2. **Thread-local storage** - Minimize shared memory access
3. **Arrow file streaming** - Don't materialize large tables
4. **Sparse operations** - Use sparse matrices for quantification

### Threading Best Practices
1. **Partition by scan density** - Balance workload across threads
2. **Minimize synchronization** - Use thread-local storage
3. **Avoid nested threading** - Don't spawn threads from threads
4. **Pre-allocate** - Allocate buffers once per thread

### I/O Optimization
1. **Batch file operations** - Minimize file open/close
2. **Async I/O** - Use Julia's async features for large files
3. **Compression** - Use Arrow compression for temporary files
4. **Memory mapping** - For read-only large files

## Testing SearchMethods

### Unit Testing Individual Methods
```julia
# Test method in isolation
function test_my_search_method()
    # Create minimal SearchContext
    context = createTestSearchContext()
    method = MySearchMethod(context)
    
    # Run method
    performSearch!(method)
    summarizeResults!(method)
    
    # Validate results
    results = getResults(method)
    @test length(results.my_result) > 0
end
```

### Integration Testing
```julia
# Test method within pipeline
function test_method_integration()
    # Run pipeline up to your method
    context = SearchDIA("test_params.json", stop_at="MySearch")
    
    # Validate state
    @test hasResults(context, "FirstPass")
    @test hasModels(context, "RT")
end
```

### Debugging Tips
1. **Enable verbose logging** in parameters
2. **Check intermediate outputs** in temp_folder
3. **Use @show macros** in processChunk! for thread debugging
4. **Validate SearchContext state** between methods
5. **Plot QC metrics** for visual debugging

## Common Pitfalls

### Data Flow Issues
- **Missing dependencies**: Ensure previous methods completed successfully
- **State mutations**: Don't modify shared SearchContext in threads
- **Parameter extraction**: Always validate parameter types/ranges

### Threading Issues  
- **Race conditions**: Use thread-local storage exclusively
- **Memory contention**: Minimize shared memory access
- **Load balancing**: Test with different thread counts

### Memory Issues
- **Memory leaks**: Clear large temporary arrays
- **OOM errors**: Use iterators instead of materializing
- **GC pressure**: Pre-allocate frequently used objects

## Advanced Patterns

### Cross-Method Communication
```julia
# Access results from previous methods
function getFirstPassResults(method::MySearchMethod)
    context = getSearchContext(method)
    return getSearchResults(context, "FirstPass")
end

# Store results for downstream methods
function summarizeResults!(method::MySearchMethod)
    results = getResults(method)
    context = getSearchContext(method)
    setSearchResults!(context, "MySearch", results)
end
```

### Model Integration
```julia
# Use calibrated models from previous methods
function processChunk!(method::MySearchMethod, batch_id::Int, thread_id::Int)
    context = getSearchContext(method)
    rt_model = getRTModel(context)
    mass_error_model = getMassErrorModel(context)
    
    # Apply models in processing
    corrected_rt = applyRTModel(rt_model, observed_rt)
    expected_mass_error = predictMassError(mass_error_model, mz)
end
```

### Quality Control Integration
```julia
function summarizeResults!(method::MySearchMethod)
    # Always generate QC plots
    generateQCPlots(method)
    
    # Export diagnostic data
    exportDiagnostics(method)
    
    # Validate result quality
    validateResults(method)
end
```

## File Organization

- `SearchMethods.jl` - Main module with exports
- `SearchTypes.jl` - Abstract type definitions
- `{MethodName}/` - Each method in its own folder
  - `{MethodName}.jl` - Main implementation
  - `utils.jl` - Method-specific utilities
- Common utilities in parent `CommonSearchUtils/`