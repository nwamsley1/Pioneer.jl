# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with CommonSearchUtils - the algorithmic core of Pioneer.jl's SearchDIA pipeline.

## Overview

CommonSearchUtils implements the performance-critical algorithms that power DIA analysis across all search methods. These utilities provide shared functionality for spectral matching, quantification matrix construction, retention time indexing, and efficient data management.

## Core Utilities

### Peak Matching Engine (`matchPeaks.jl`)

**Purpose**: High-performance algorithm for matching theoretical fragments to empirical peaks

**Key Algorithm**: Single-pass O(T+P) algorithm where T=transitions, P=peaks
- Requires pre-sorted inputs (transitions by m/z, peaks by m/z)
- Mass-error-aware matching with intensity-dependent tolerances
- Each fragment gets at most one peak, but peaks can match multiple fragments

**Core Functions**:
```julia
matchPeaks!(fragment_matches, precursor_matches, transitions, masses, intensities, ...)
setNearest!(match, transition, peak_idx, masses, intensities)  # Best peak selection
setMatch!(match, transition, intensity, mz, peak_idx)          # Create match object
```

**Performance Features**:
- `@inbounds @fastmath` optimizations
- Pre-allocated vectors with block growth
- Branchless operations for CPU pipeline efficiency

**Usage Pattern**:
```julia
# Pre-sort inputs
sort!(transitions, by=getMZ)
sort_indices = sortperm(masses)

# Perform matching
n_fragment_matches, n_precursor_matches = matchPeaks!(
    fragment_matches, precursor_matches, transitions,
    masses[sort_indices], intensities[sort_indices], ...)
```

### Quantification Matrix Construction (`buildDesignMatrix.jl`)

**Purpose**: Constructs sparse design matrices for protein quantification

**Algorithm**: 
- Rows = unique spectrum peaks, Columns = precursors
- Uses custom `SparseArray` type for memory efficiency
- Handles matched fragments (empirical intensities) and missed fragments (predicted)

**Key Functions**:
```julia
buildDesignMatrix!(Hs::SparseArray, fragment_matches, ion_templates, ...)
addToSparseDesignMatrix!(Hs, match, spline_data, weight)
```

**Memory Management**:
- Pre-allocated block growth (default 10k elements)
- Deduplication of peak matches per precursor
- In-place sparse matrix construction with final column sorting

**Usage in Search Methods**:
```julia
# Reset and build design matrix
reset!(Hs)
buildDesignMatrix!(Hs, fragment_matches, ion_templates, ...)

# Use for quantification
weights = solve(Hs, observed_intensities)
```

### Fragment Index Queries (`queryFragmentIndex.jl`)

**Purpose**: Efficiently searches hierarchical fragment index for precursor identification

**Key Algorithms**:
- **Exponential Search**: `exponentialFragmentBinSearch()` for rapid bound expansion
- **Branchless Binary Search**: `findFirstFragmentBin()` for optimal performance
- **Hierarchical Indexing**: RT bins → Fragment bins → Individual fragments

**Search Strategy**:
```julia
function searchFragmentBin!(counter, spec_lib, rt_bin_idx, peak_mz, ...)
    # 1. Find fragment bins covering m/z tolerance
    first_bin, last_bin = exponentialFragmentBinSearch(bins, peak_mz, tolerance)
    
    # 2. Binary search within bins for precursor matches
    for bin_idx in first_bin:last_bin
        for frag_idx in bin_range
            if precursor_matches_criteria(frag, peak_mz, tolerance)
                counter[getPrecID(frag)] += 1
            end
        end
    end
end
```

**Performance Optimizations**:
- Counter-based scoring with pre-allocated structures
- Cached bounds for subsequent searches
- Branchless operations for CPU efficiency

### Retention Time Indexing (`buildRTIndex.jl`)

**Purpose**: Creates efficient RT-based indexing for chromatographic window searches

**Index Structure**:
```julia
struct RetentionTimeIndex
    rt_bins::Vector{Float32}              # Bin boundaries (0.1 min default)
    precursor_ranges::Vector{UnitRange}   # Precursor ranges per bin
    sorted_precursors::Vector{UInt32}     # Precursors sorted by m/z within bins
end
```

**Key Functions**:
```julia
buildRtIndex(rts, mzs, prec_ids; bin_width=0.1f0)
makeRTIndices(passing_psms_paths, precursors, rt_conversion_models)
```

**Imputation Strategy**:
- Uses library iRT values for low-confidence PSMs
- RT-to-iRT conversion models applied per MS file
- Fallback to library values when empirical RT unavailable

### Work Distribution (`partitionThreadTasks.jl`)

**Purpose**: Intelligent task partitioning for multithreaded processing

**Strategies**:
```julia
partitionScansToThreads(scan_ranges, n_threads)  # RT-aware distribution
partitionTasks(total_work, n_threads)            # Generic work splitting
```

**Load Balancing**:
- Considers scan density across RT range
- Maintains data locality for cache efficiency
- Separate handling for MS1 vs MS2 scans

### Data Merging (`mergeSortedArrowTables.jl`)

**Purpose**: Memory-efficient merging of large sorted Arrow tables

**Algorithm**: Min-heap-based k-way merge
- Configurable batch sizes to control memory usage
- Supports single and dual sort keys
- Streaming write to avoid memory overflow

**Usage Pattern**:
```julia
# Merge PSM tables from parallel processing
merge_sorted_arrow_tables(
    input_paths,
    output_path,
    sort_columns=["scan_idx", "score"],
    batch_size=100000
)
```

### Quantification Normalization (`normalizeQuant.jl`)

**Purpose**: RT-dependent quantification normalization using spline-based median correction

**Process**:
1. Fits splines to median quantification values across RT
2. Calculates file-specific offsets from global median  
3. Applies corrections to normalize between runs

**Implementation**:
```julia
function normalizeQuantification!(quant_matrix, rt_values, ...)
    # Fit splines to median profiles
    median_splines = fitQuantSplines(quant_matrix, rt_values)
    
    # Calculate normalization factors
    norm_factors = calculateNormFactors(median_splines)
    
    # Apply corrections
    applyNormalization!(quant_matrix, norm_factors)
end
```

## Data Flow Integration

**Typical Search Method Flow**:
1. **RT Indexing** → Fast precursor lookup by retention time
2. **Transition Selection** → Choose fragments for current window
3. **Peak Matching** → Match fragments to spectrum peaks  
4. **Design Matrix** → Construct quantification matrix from matches
5. **Fragment Index Query** → Score precursor candidates
6. **Data Merging** → Combine results across batches
7. **Normalization** → Apply quantification corrections

## Performance Patterns

### Memory Management
- **Pre-allocation**: Working arrays sized once with block growth
- **Sparse Structures**: Custom `SparseArray` for memory efficiency
- **In-place Operations**: Minimize allocations in hot paths

### Threading Optimizations  
- **Data Locality**: RT-aware task partitioning for cache efficiency
- **Thread-local Storage**: Avoid contention with per-thread data structures
- **Batch Processing**: Optimal cache utilization patterns

### Algorithmic Optimizations
- **Branchless Binary Search**: CPU pipeline efficiency
- **Single-pass Algorithms**: O(T+P) vs O(T*P) complexity
- **Exponential Search**: Rapid bound finding for sparse data
- **Cached Bounds**: Avoid redundant computations

## Key Data Structures

```julia
# Custom sparse array with metadata
struct SparseArray{Ti<:Integer, T<:AbstractFloat}
    colptr::Vector{Ti}      # Column pointers
    rowval::Vector{Ti}      # Row indices
    nzval::Vector{T}        # Non-zero values
    isotope::Vector{UInt8}  # Isotope metadata
    matched::Vector{Bool}   # Match indicators
end

# Optimized fragment match representation
struct FragmentMatch{T<:AbstractFloat}
    predicted_intensity::T
    intensity::T
    theoretical_mz::T
    match_mz::T
    # ... compact encoding fields
end

# High-performance scoring accumulator
struct Counter{K,V}
    keys::Vector{K}
    values::Vector{V}
    n::Int64
end
```

## Development Patterns

### Extending Core Utilities

**Adding New Matching Algorithm**:
```julia
function matchPeaksCustom!(fragment_matches, transitions, masses, ...)
    # Custom matching logic
    # Maintain pre-sorted input requirement
    # Use block-based growth for output arrays
    # Return (n_fragment_matches, n_precursor_matches)
end
```

**Custom Scoring Integration**:
```julia
function scorePrecursorsCustom!(counter, fragment_matches, ...)
    # Reset counter
    reset!(counter)
    
    # Accumulate scores
    for match in fragment_matches
        prec_id = getPrecID(match)
        score = calculateCustomScore(match)
        addCount!(counter, prec_id, score)
    end
end
```

### Performance Optimization Guidelines

**Memory Optimization**:
- Use `@inbounds @fastmath` in hot loops after validation
- Pre-allocate working arrays with appropriate block sizes
- Consider sparse vs dense representations based on data characteristics

**Threading Considerations**:
- Design algorithms for data locality
- Use thread-local storage for mutable state
- Partition work by RT or m/z ranges for cache efficiency

**Algorithm Design**:
- Leverage sorted data structures for binary search
- Implement branchless operations where possible
- Cache frequently computed values

## Testing and Debugging

### Unit Testing Approach
```julia
# Test peak matching with synthetic data
function test_peak_matching()
    transitions = create_test_transitions()
    masses = [100.0, 200.0, 300.0]
    intensities = [1000.0, 2000.0, 1500.0]
    
    fragment_matches = Vector{FragmentMatch{Float32}}(undef, 1000)
    n_matches = matchPeaks!(fragment_matches, transitions, masses, ...)
    
    @test n_matches > 0
    @test all(match -> getMZ(match) > 0, fragment_matches[1:n_matches])
end
```

### Debug Entry Points
- `matchPeaks.jl:306-429` - Main spectral matching loop
- `buildDesignMatrix.jl:42-91` - Matrix construction algorithm  
- `queryFragmentIndex.jl:271-327` - Fragment index search
- Remove `@inbounds @fastmath` for bounds checking during debug

### Performance Profiling
```julia
using Profile, PProf

# Profile peak matching
@profile matchPeaks!(fragment_matches, transitions, masses, ...)
pprof()

# Monitor memory allocations
@time @allocated buildDesignMatrix!(Hs, fragment_matches, ...)
```

### Common Debugging Issues
1. **Unsorted Inputs**: Algorithms assume sorted data - validate sort order
2. **Array Bounds**: Check pre-allocated array sizes vs actual data
3. **Memory Pressure**: Monitor sparse array growth and GC pressure
4. **Thread Safety**: Ensure thread-local storage for mutable state

## Integration with Search Methods

### Search Method Requirements
Each search method typically needs:
1. **Thread Partitioning**: Call `partitionThreadTasks` for work distribution
2. **Peak Matching**: Use `matchPeaks!` for spectral matching
3. **Scoring**: Apply `queryFragmentIndex` or custom scoring
4. **Quantification**: Build design matrix with `buildDesignMatrix!`

### Example Integration Pattern
```julia
function processChunk!(method::MySearchMethod, batch_id::Int, thread_id::Int)
    search_data = getSearchData(method)[thread_id]
    
    # Get spectrum data for this chunk
    spectrum = getSpectrumData(batch_id)
    
    # Select transitions for this RT window
    ion_idx = selectTransitions!(search_data.ion_templates, ...)
    
    # Match peaks to transitions
    n_matches = matchPeaks!(
        search_data.fragment_matches,
        search_data.ion_templates[1:ion_idx],
        spectrum.masses, spectrum.intensities, ...
    )
    
    # Build quantification matrix
    buildDesignMatrix!(
        search_data.Hs,
        search_data.fragment_matches[1:n_matches], ...
    )
    
    # Score precursors (if needed)
    queryFragmentIndex!(search_data.counter, spec_lib, ...)
end
```

## File Organization

- **buildDesignMatrix.jl** - Sparse matrix construction for quantification
- **buildRTIndex.jl** - Retention time indexing infrastructure  
- **entrapmentAnalysis.jl** - FDR estimation with entrapment strategies
- **matchPeaks.jl** - Core spectral matching algorithms
- **mergeSortedArrowTables.jl** - Memory-efficient data merging
- **normalizeQuant.jl** - Cross-run quantification normalization
- **partitionThreadTasks.jl** - Intelligent work distribution
- **queryFragmentIndex.jl** - Fragment index search and scoring
- **selectTransitions/** - Modular transition selection framework

## Module-Specific Documentation

- Transition Selection: See `selectTransitions/CLAUDE.md` for detailed guidance on the transition selection framework