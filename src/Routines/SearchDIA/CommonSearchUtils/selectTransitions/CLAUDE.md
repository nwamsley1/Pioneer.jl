# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the selectTransitions modular framework.

## Overview

The selectTransitions framework implements a sophisticated strategy pattern for transition selection in DIA proteomics searches. It provides four distinct strategies optimized for different search phases, with a unified interface supporting both partial and full precursor capture modes.

## Strategy Pattern Architecture

### Core Type Hierarchy
```julia
abstract type TransitionSelectionStrategy end
abstract type PrecEstimation end

# Strategy implementations
struct StandardTransitionSelection <: TransitionSelectionStrategy end
struct MassErrEstimationStrategy <: TransitionSelectionStrategy end  
struct RTIndexedTransitionSelection <: TransitionSelectionStrategy end
struct QuadEstimationTransitionSelection <: TransitionSelectionStrategy end

# Precursor capture modes
struct PartialPrecCapture <: PrecEstimation end  # Detailed isotope calculation
struct FullPrecCapture <: PrecEstimation end     # Simplified calculation
```

### Unified Interface
```julia
function selectTransitions!(
    transitions::Vector{DetailedFrag{Float32}},
    strategy::TransitionSelectionStrategy,
    prec_estimation_type::PrecEstimation,
    common_args...;
    kwargs...)
    
    # Strategy-specific implementation via multiple dispatch
    transition_idx, n_precursors = _select_transitions_impl!(
        transitions, strategy, prec_estimation_type, common_args...; kwargs...
    )
    
    # Common post-processing: sort by m/z
    sort!(@view(transitions[1:transition_idx]), by = x->getMZ(x))
    
    return transition_idx, n_precursors
end
```

## Strategy Implementations

### 1. StandardTransitionSelection

**Purpose**: Basic transition selection for standard library searches

**Algorithm**:
- Filters precursors by retention time tolerance around target iRT
- Applies isotope error bounds checking within quadrupole window
- Processes all fragments for qualifying precursors up to max rank
- Uses spline-based intensity prediction when available

**Performance**: O(n) where n = precursors in RT window

**Key Decision Logic**:
```julia
# RT filtering
if abs(prec_irts[prec_idx] - iRT) > iRT_tol
    continue
end

# Isotope error bounds
mz_low = getPrecMinBound(quad_func) - NEUTRON*first(isotope_err_bounds)/prec_charge
mz_high = getPrecMaxBound(quad_func) + NEUTRON*last(isotope_err_bounds)/prec_charge

if (prec_mz < mz_low) || (prec_mz > mz_high)
    continue
end
```

**Usage**: Primary library search method (LibrarySearch.jl)

### 2. RTIndexedTransitionSelection

**Purpose**: Highly efficient transition selection using pre-built RT index

**Algorithm**:
- Uses retention time index for O(log n) precursor lookup
- Supports optional precursor filtering via Set{UInt32}
- Processes RT bins within specified start/stop range
- Binary searches within RT bins for quadrupole window bounds

**Performance**: O(log n + k) where k = precursors in RT window

**Key Optimization**:
```julia
# Efficient RT bin lookup
rt_bin_range = getRTBinRange(rt_index, rt_start, rt_stop)

for rt_bin_idx in rt_bin_range
    precs = getPrecursorsInBin(rt_index, rt_bin_idx)
    
    # Binary search within RT bin
    start_idx, stop_idx = getPrecursorWindowRange(
        precs, min_prec_mz, max_prec_mz, isotope_err_bounds
    )
    
    # Process only relevant precursors
    for prec_idx in start_idx:stop_idx
        # ... process precursor
    end
end
```

**Usage**: Multi-pass searches (SecondPassSearch, IntegrateChromatograms, HuberTuning)

### 3. MassErrEstimationStrategy

**Purpose**: Lightweight fragment selection for mass error calibration

**Algorithm**:
- Minimal processing - only rank filtering (max_rank = 5)
- No isotope calculations or spline evaluations  
- Direct fragment-to-DetailedFrag conversion
- Optimized for parameter tuning phase

**Performance**: O(n*f) where f = fragments per precursor (very fast)

**Simplification**:
```julia
# Minimal processing for speed
for frag in precursor_fragments
    if getRank(frag) <= max_rank
        transition_idx += 1
        transitions[transition_idx] = DetailedFrag(
            Float32(frag.mz),
            Float32(1.0),  # Unit intensity
            frag.frag_charge,
            false,         # Not isotope
            frag.prec_id
        )
    end
end
```

**Usage**: Parameter tuning searches for mass error modeling

### 4. QuadEstimationTransitionSelection

**Purpose**: Specialized selection for quadrupole transmission modeling

**Algorithm**:
- Processes multiple isotope indices (typically 0:2) per precursor
- Uses PartialPrecCapture mode exclusively
- Generates artificial precursor IDs encoding isotope information
- Focuses on high-quality fragments (rank ≤ 7) across isotope pattern

**Performance**: O(n*i*f) where i = isotope indices (3x slower than standard)

**Isotope Processing**:
```julia
for prec_iso_idx in isotope_range  # Usually 0:2
    # Set up transmission for single isotope
    fill!(precursor_transmission, zero(Float32))
    precursor_transmission[prec_iso_idx + 1] = one(Float32)
    
    # Process fragments with isotope-specific precursor ID
    for frag in precursor_fragments
        if getRank(frag) <= 7  # High-quality fragments only
            # Encode isotope info in precursor ID
            modified_prec_id = UInt32((getPID(frag) - 1) * 3 + (prec_iso_idx + 1))
            
            fillTransitionList!(transitions, transition_idx, frag, 
                              modified_prec_id, precursor_transmission, ...)
        end
    end
end
```

**Usage**: Quadrupole tuning searches

## Core Framework Components

### Transition Population (`fillTransitionList.jl`)

**Core Algorithm**: Converts library fragments to transitions with proper isotope patterns

```julia
function fillTransitionList!(
    transitions::Vector{DetailedFrag{Float32}},
    transition_idx::Int64,
    frag::LibraryFragmentIon,
    prec_id::UInt32,
    precursor_transmission::Vector{Float32},
    prec_estimation_type::PrecEstimation,
    isotopes::Vector{Float32},
    iso_splines::IsotopeSplineModel,
    nce_data::Union{Missing, Float32},
    quad_transmission_func::QuadTransmissionModel,
    n_frag_isotopes::Int64,
    block_size::Int64
) -> Int64
```

**Key Steps**:
1. **Precursor Isotope Setup**: Calculate transmission through quadrupole
2. **Fragment Processing**: For each fragment in precursor
3. **Isotope Calculation**: Generate isotope pattern based on capture mode
4. **Transition Creation**: Convert to DetailedFrag with m/z shifts

**Isotope Pattern Generation**:
```julia
# Calculate isotopes based on capture mode
if prec_estimation_type isa PartialPrecCapture
    # Detailed abundance calculation for partial capture
    getFragAbundance!(isotopes, precursor_transmission, sulfur_count, 
                      base_mass, n_frag_isotopes, iso_splines)
else  # FullPrecCapture
    # Direct spline evaluation for full capture
    for iso_idx in 0:(n_frag_isotopes-1)
        isotopes[iso_idx + 1] = iso_splines(sulfur_count, iso_idx, base_mass)
    end
end

# Apply NCE adjustment if available
if !ismissing(nce_data)
    adjustForNCE!(isotopes, nce_data, frag)
end
```

### Array Management
```julia
function ensureTransitionCapacity!(transitions, required_size, block_size)
    if required_size > length(transitions)
        new_size = required_size + block_size
        resize!(transitions, new_size)
    end
end
```

## Precursor Capture Strategies

### PartialPrecCapture
- **Use Case**: Isolation window captures only part of precursor isotope envelope
- **Calculation**: Complex abundance calculation accounting for partial capture
- **Performance**: Slower but more accurate
- **Implementation**: Uses detailed isotope modeling with transmission functions

### FullPrecCapture  
- **Use Case**: Full precursor isotope pattern is captured
- **Calculation**: Direct isotope spline evaluation
- **Performance**: Faster approximation
- **Implementation**: Direct lookup: `iso_splines(sulfur_count, iso_idx, mass)`

## Performance Characteristics

| Strategy | Time Complexity | Memory Usage | Accuracy | Use Case |
|----------|----------------|--------------|----------|----------|
| Standard | O(n) | Moderate | High | Basic library search |
| RTIndexed | O(log n + k) | Low | High | Multi-pass efficient search |
| MassErr | O(n*f) | Minimal | Medium | Parameter calibration |
| QuadEst | O(n*i*f) | High | Very High | Quadrupole modeling |

### Performance Optimization Tips
1. **Use RTIndexed** for repeated searches with same RT index
2. **Use MassErr** for calibration phases requiring speed over accuracy
3. **Batch process** multiple spectra to amortize setup costs
4. **Pre-allocate** transition vectors with appropriate capacity

## Integration with Search Methods

### Strategy Selection by Search Method
```julia
# Typical usage patterns
search_methods_to_strategies = Dict(
    LibrarySearch => StandardTransitionSelection(),
    SecondPassSearch => RTIndexedTransitionSelection(),
    IntegrateChromatograms => RTIndexedTransitionSelection(),
    ParameterTuning => MassErrEstimationStrategy(),
    QuadTuning => QuadEstimationTransitionSelection(),
    HuberTuning => RTIndexedTransitionSelection()
)
```

### Usage Pattern in Search Methods
```julia
function processSpectrum!(search_method, spectrum_data, thread_id)
    search_data = getSearchData(search_method)[thread_id]
    
    # Select appropriate strategy
    strategy = getStrategy(search_method)
    prec_estimation = getParams(search_method).prec_estimation
    
    # Select transitions for current window
    ion_idx, n_precursors = selectTransitions!(
        getIonTemplates(search_data),
        strategy,
        prec_estimation,
        spectrum_data.rt, spectrum_data.mz_range,
        getSpecLib(search_method), ...
    )
    
    # Use selected transitions for matching
    matchPeaks!(fragment_matches, 
               view(getIonTemplates(search_data), 1:ion_idx),
               spectrum_data.masses, spectrum_data.intensities, ...)
end
```

## Configuration and Parameters

### Key Parameters from Configuration Files
```json
{
    "isotope_settings": {
        "err_bounds_first_pass": [1, 0],
        "err_bounds_quant_search": [2, 0],   
        "partial_capture": false,
        "min_fraction_transmitted": 0.25
    },
    "fragment_settings": {
        "max_rank": 25,
        "n_isotopes": 1,
        "min_count": 4
    },
    "rt_settings": {
        "irt_tol": 5.0,
        "rt_bin_tol": 0.1
    }
}
```

### Runtime Configuration
- `isotope_err_bounds`: Controls isotope error tolerance
- `max_frag_rank`: Maximum fragment rank to include
- `iRT_tol`: Retention time tolerance for precursor filtering
- `block_size`: Array growth increment (default 10000)

## Extension Patterns

### Adding New Selection Strategy

**1. Define Strategy Type**:
```julia
struct NewTransitionSelection <: TransitionSelectionStrategy end
```

**2. Implement Selection Logic**:
```julia
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::NewTransitionSelection,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    # ... custom arguments
    spec_lib::SpectralLibrary,
    target_rt::Float32,
    mz_bounds::Tuple{Float32, Float32},
    kwargs...
) 
    n_precursors = 0
    
    # Custom precursor selection logic
    for prec_idx in relevant_precursors
        # Apply custom filtering criteria
        if custom_criteria(prec_idx, target_rt, mz_bounds)
            # Process precursor fragments
            transition_idx = fillTransitionList!(
                transitions, transition_idx, 
                precursor_fragments[prec_idx], prec_idx,
                precursor_transmission, prec_estimation_type, ...)
            n_precursors += 1
        end
    end
    
    return transition_idx, n_precursors
end
```

**3. Update Search Method Integration**:
```julia
# Add to search method's strategy selection
function getStrategy(::MyNewSearchMethod)
    return NewTransitionSelection()
end
```

### Extension Guidelines
- Maintain return signature: `(transition_idx, n_precursors)`
- Use `ensureTransitionCapacity!` for array growth
- Leverage existing `fillTransitionList!` for transition creation
- Consider performance implications for your use case
- Test with different precursor estimation modes

## Testing and Validation

### Unit Testing Framework
```julia
function test_transition_selection_strategy()
    # Setup test data
    transitions = Vector{DetailedFrag{Float32}}(undef, 10000)
    spec_lib = create_test_library()
    
    # Test strategy
    ion_idx, n_precs = selectTransitions!(
        transitions,
        StandardTransitionSelection(),
        PartialPrecCapture(),
        0,  # transition_idx
        100.0f0,  # target_rt
        (400.0f0, 800.0f0),  # mz_bounds
        spec_lib,
        # ... other parameters
    )
    
    # Validate outputs
    @test ion_idx > 0
    @test n_precs > 0
    @test all(getMZ(transitions[i]) <= getMZ(transitions[i+1]) for i in 1:ion_idx-1)  # Sorted
    @test all(transitions[i].intensity > 0 for i in 1:ion_idx)  # Valid intensities
    @test all(getPrecID(transitions[i]) > 0 for i in 1:ion_idx)  # Valid precursor IDs
end
```

### Performance Benchmarking
```julia
using BenchmarkTools

# Compare strategy performance
@benchmark selectTransitions!(transitions, StandardTransitionSelection(), ...)
@benchmark selectTransitions!(transitions, RTIndexedTransitionSelection(), ...)
@benchmark selectTransitions!(transitions, MassErrEstimationStrategy(), ...)

# Memory allocation analysis
@time @allocated selectTransitions!(transitions, strategy, ...)
```

### Validation Checks
```julia
function validate_transitions(transitions, n_transitions)
    # Check sorting
    @assert all(getMZ(transitions[i]) <= getMZ(transitions[i+1]) for i in 1:n_transitions-1)
    
    # Check valid m/z values
    @assert all(getMZ(transitions[i]) > 0 for i in 1:n_transitions)
    
    # Check intensity values
    @assert all(getIntensity(transitions[i]) >= 0 for i in 1:n_transitions)
    
    # Check precursor IDs
    @assert all(getPrecID(transitions[i]) > 0 for i in 1:n_transitions)
end
```

## Common Development Patterns

### Debugging Transition Selection
```julia
# Enable debug output
function debug_selectTransitions!(transitions, strategy, ...)
    println("Strategy: $(typeof(strategy))")
    println("Input parameters: rt=$target_rt, mz_bounds=$mz_bounds")
    
    ion_idx, n_precs = selectTransitions!(transitions, strategy, ...)
    
    println("Selected transitions: $ion_idx")
    println("Precursors processed: $n_precs")
    println("m/z range: $(getMZ(transitions[1])) - $(getMZ(transitions[ion_idx]))")
    
    return ion_idx, n_precs
end
```

### Memory Management
```julia
# Optimize for your workload
function tune_transition_selection(expected_transitions_per_spectrum)
    # Adjust block size based on expected data
    optimal_block_size = max(1000, expected_transitions_per_spectrum * 2)
    
    # Pre-allocate with optimal size
    transitions = Vector{DetailedFrag{Float32}}(undef, optimal_block_size)
    
    return transitions
end
```

### Error Handling
```julia
function safe_selectTransitions!(transitions, strategy, ...)
    try
        return selectTransitions!(transitions, strategy, ...)
    catch e
        @error "Transition selection failed" strategy=typeof(strategy) exception=e
        return 0, 0  # Return empty selection
    end
end
```

## File Organization

- **selectTransitions.jl** - Main interface and strategy coordination
- **fillTransitionList.jl** - Core transition list population algorithm  
- **standardTransitionSelection.jl** - Basic library search strategy
- **rtIndexTransitionSelection.jl** - RT index-based efficient strategy
- **massErrEstimationStrategy.jl** - Fast strategy for parameter tuning
- **quadEstimationSelection.jl** - Specialized strategy for quadrupole modeling

The selectTransitions framework demonstrates excellent modular design, enabling efficient and flexible transition selection across different search phases while maintaining high performance and extensibility.