"""
Memory optimization utilities for IntegrateChromatogramSearch.

This file provides drop-in replacement optimizations for the key memory bottlenecks
in build_chromatograms without requiring a complete rewrite of the function.
"""

"""
    estimate_chromatogram_capacity(scan_range::Vector{Int64}, 
                                  precursors_passing::Set{UInt32},
                                  avg_precursors_per_scan::Int = 50)

Estimate the total capacity needed for chromatogram arrays based on scan count 
and precursor characteristics.

# Arguments
- `scan_range`: Range of scans to process
- `precursors_passing`: Set of precursors that pass filters
- `avg_precursors_per_scan`: Average precursors per scan (default: 50)

# Returns
- Estimated total capacity needed for chromatogram array
"""
function estimate_chromatogram_capacity(scan_range::Vector{Int64}, 
                                       precursors_passing::Set{UInt32},
                                       avg_precursors_per_scan::Int = 50)
    n_scans = length(scan_range)
    n_precursors = length(precursors_passing)
    
    # Conservative estimate: use the smaller of total precursors or avg per scan
    effective_precursors_per_scan = min(avg_precursors_per_scan, n_precursors)
    
    base_estimate = n_scans * effective_precursors_per_scan
    safety_factor = 1.3  # 30% buffer for variability
    
    return ceil(Int, base_estimate * safety_factor)
end

"""
    initialize_chromatogram_array_optimized(::Type{T}, scan_range::Vector{Int64},
                                           precursors_passing::Set{UInt32}) where T

Create an optimally pre-allocated chromatogram array instead of starting with 500k elements.

# Arguments
- `T`: ChromObject type (MS1ChromObject or MS2ChromObject)  
- `scan_range`: Range of scans to process
- `precursors_passing`: Set of precursors that pass filters

# Returns
- Pre-allocated vector with estimated optimal size
"""
function initialize_chromatogram_array_optimized(::Type{T}, scan_range::Vector{Int64},
                                                 precursors_passing::Set{UInt32}) where T
    estimated_capacity = estimate_chromatogram_capacity(scan_range, precursors_passing)
    
    # Use estimated capacity instead of fixed 500k
    return Vector{T}(undef, estimated_capacity)
end

"""
    smart_chromatogram_resize!(chromatograms::Vector{T}, current_idx::Int, 
                              needed_capacity::Int) where T

Efficiently resize chromatogram array when more space is needed.
Uses exponential growth strategy instead of fixed 500k chunks.

# Arguments
- `chromatograms`: Current chromatogram array
- `current_idx`: Current index in the array
- `needed_capacity`: How much total capacity is needed

# Returns
- Nothing (modifies chromatograms in place)
"""
function smart_chromatogram_resize!(chromatograms::Vector{T}, current_idx::Int, 
                                   needed_capacity::Int) where T
    current_capacity = length(chromatograms)
    
    if needed_capacity > current_capacity
        # Calculate optimal growth amount
        # Use exponential growth: max(needed, 1.5x current)
        growth_factor = 1.5
        new_capacity = max(needed_capacity, ceil(Int, current_capacity * growth_factor))
        additional_needed = new_capacity - current_capacity
        
        # Grow the array efficiently
        append!(chromatograms, Vector{T}(undef, additional_needed))
    end
end

"""
    optimize_build_chromatograms_memory()

Instructions for integrating memory optimizations into build_chromatograms:

Replace this pattern:
```julia
chromatograms = Vector{MS1ChromObject}(undef, 500000)  # Initial size
```

With:
```julia
chromatograms = initialize_chromatogram_array_optimized(MS1ChromObject, scan_range, precursors_passing)
```

Replace this pattern:
```julia
if rt_idx + 1 > length(chromatograms)
    append!(chromatograms, Vector{MS2ChromObject}(undef, 500000))
end
```

With:
```julia
smart_chromatogram_resize!(chromatograms, rt_idx, rt_idx + prec_temp_size)
```

These changes will:
1. Start with optimal initial size instead of fixed 500k
2. Grow efficiently when needed instead of fixed 500k chunks
3. Reduce memory allocation overhead significantly
4. Maintain identical functionality
"""
function optimize_build_chromatograms_memory()
    println("""
    Memory Optimization Integration Guide:
    
    1. Replace initial allocation:
       OLD: chromatograms = Vector{MS1ChromObject}(undef, 500000)
       NEW: chromatograms = initialize_chromatogram_array_optimized(MS1ChromObject, scan_range, precursors_passing)
    
    2. Replace resize pattern:
       OLD: if rt_idx + 1 > length(chromatograms)
                append!(chromatograms, Vector{MS2ChromObject}(undef, 500000))
            end
       NEW: smart_chromatogram_resize!(chromatograms, rt_idx, rt_idx + prec_temp_size)
    
    3. Benefits:
       - Eliminates 4 hard-coded 500k append! operations
       - Uses data-driven size estimation
       - Exponential growth strategy for efficiency
       - Significantly reduced memory allocation overhead
    """)
end