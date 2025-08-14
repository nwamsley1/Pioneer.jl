# Plan: Simplify Iteration Strategy to Fixed 3-Phase Alternating

## Status: ✅ IMPLEMENTED (2025-01-14)

## Overview

Simplify the iteration strategy to always use exactly 3 phases with a fixed alternating pattern:
1. **Phase 1**: Centered at zero bias
2. **Phase 2**: Shifted right by +max_tolerance_ppm
3. **Phase 3**: Shifted left by -max_tolerance_ppm

This removes unnecessary complexity and configuration options while maintaining robust parameter space exploration.

## Current Implementation (To Be Removed)

### Configurable Parameters to Remove
- `max_phases` - Currently allows 1-N phases
- `bias_shift_strategy` - Currently allows alternating/positive_first/negative_first
- `bias_shift_magnitude` - Currently allows custom values or max_tolerance

### Current IterationSettings Structure
```julia
struct IterationSettings
    mass_tolerance_scale_factor::Float32
    iterations_per_phase::Int64
    max_phases::Int64                         # TO BE REMOVED
    bias_shift_strategy::Symbol               # TO BE REMOVED
    bias_shift_magnitude::Union{Float32, Symbol}  # TO BE REMOVED
end
```

## Proposed Simplified Implementation

### New IterationSettings Structure
```julia
struct IterationSettings
    mass_tolerance_scale_factor::Float32     # KEEP: User-configurable scaling
    iterations_per_phase::Int64              # KEEP: User-configurable iterations
    
    function IterationSettings(
        scale_factor::Float32 = 2.0f0,
        iterations_per_phase::Int64 = 3
    )
        @assert scale_factor > 1.0f0 "Scale factor must be greater than 1"
        @assert iterations_per_phase > 0 "Iterations per phase must be positive"
        new(scale_factor, iterations_per_phase)
    end
end
```

### Simplified JSON Configuration

#### Before (Complex)
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,
    "iterations_per_phase": 3,
    "max_phases": 3,                    // REMOVE
    "bias_shift_strategy": "alternating", // REMOVE
    "bias_shift_magnitude": "max_tolerance" // REMOVE
}
```

#### After (Simple)
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,
    "iterations_per_phase": 3
}
```

## Implementation Changes

### 1. Update types.jl

#### Remove fields from IterationSettings
```julia
# OLD
struct IterationSettings
    mass_tolerance_scale_factor::Float32
    iterations_per_phase::Int64
    max_phases::Int64                         # REMOVE
    bias_shift_strategy::Symbol               # REMOVE
    bias_shift_magnitude::Union{Float32, Symbol}  # REMOVE
    ...
end

# NEW
struct IterationSettings
    mass_tolerance_scale_factor::Float32
    iterations_per_phase::Int64
    
    function IterationSettings(
        scale_factor::Float32 = 2.0f0,
        iterations_per_phase::Int64 = 3
    )
        @assert scale_factor > 1.0f0 "Scale factor must be greater than 1"
        @assert iterations_per_phase > 0 "Iterations per phase must be positive"
        new(scale_factor, iterations_per_phase)
    end
end
```

#### Update parameter extraction
```julia
# In ParameterTuningSearchParameters constructor
iteration_settings = if hasproperty(search_params, :iteration_settings)
    iter = search_params.iteration_settings
    IterationSettings(
        hasproperty(iter, :mass_tolerance_scale_factor) ? 
            Float32(iter.mass_tolerance_scale_factor) : 2.0f0,
        hasproperty(iter, :iterations_per_phase) ? 
            Int64(iter.iterations_per_phase) : 3
    )
else
    IterationSettings()  # Use defaults
end
```

### 2. Update ParameterTuningSearch.jl

#### Simplify next_iteration!
```julia
function next_iteration!(state::IterationState, settings::IterationSettings)
    state.total_iterations += 1
    state.current_iteration_in_phase += 1
    
    # Check if we need to transition to a new phase
    if state.current_iteration_in_phase > settings.iterations_per_phase
        state.current_phase += 1
        state.current_iteration_in_phase = 1
    end
    
    # Always exactly 3 phases
    return state.current_phase <= 3
end
```

#### Simplify calculate_phase_bias_shift
```julia
function calculate_phase_bias_shift(phase::Int64, params)::Float32
    max_tol = getMaxTolerancePpm(params)
    
    if phase == 1
        return 0.0f0  # No shift in first phase
    elseif phase == 2
        return max_tol  # Positive shift in second phase
    elseif phase == 3
        return -max_tol  # Negative shift in third phase
    else
        error("Invalid phase: $phase. Only phases 1-3 are supported.")
    end
end
```

#### Update reset_for_new_phase!
```julia
function reset_for_new_phase!(search_context, ms_file_idx, params, phase::Int64, iteration_state::IterationState)
    initial_tolerance = getFragTolPpm(params)
    
    # Calculate bias shift for this phase (simplified)
    bias_shift = calculate_phase_bias_shift(phase, params)
    push!(iteration_state.phase_bias_shifts, bias_shift)
    
    # Reset to initial tolerance with new bias
    new_model = create_capped_mass_model(
        bias_shift,
        initial_tolerance,
        initial_tolerance,
        getMaxTolerancePpm(params)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
    
    @info "Phase $phase: Reset to initial tolerance (±$(round(initial_tolerance, digits=1)) ppm) " *
          "with bias shift $(round(bias_shift, digits=1)) ppm"
end
```

#### Update process_file! logging
```julia
# Simplified logging in process_file!
@info "Processing file: $parsed_fname (index: $ms_file_idx)"
@info "Iteration settings: scale_factor=$(settings.mass_tolerance_scale_factor), " *
      "iterations_per_phase=$(settings.iterations_per_phase)"
# Remove mention of max_phases since it's always 3

# Later in the function
@info "Completed $(iteration_state.total_iterations) total iterations across " *
      "$(iteration_state.current_phase) phases (max 3). Converged: $converged"
```

### 3. Clean Up Documentation

#### Update GENERALIZED_ITERATION_PLAN.md
- Remove references to configurable strategies
- Document fixed 3-phase approach
- Update examples to remove deprecated parameters

#### Update PARAMETER_TUNING_IMPROVEMENTS_SUMMARY.md
- Remove documentation for bias_shift_strategy
- Remove documentation for bias_shift_magnitude
- Remove documentation for max_phases
- Update all examples

#### Remove BIAS_SHIFT_STRATEGIES_EXPLAINED.md
- This file becomes obsolete with the simplified approach
- Could be replaced with a simpler explanation of the fixed 3-phase system

## Benefits of Simplification

### 1. **Reduced Complexity**
- Fewer parameters to configure
- Less chance for user error
- Clearer behavior

### 2. **Predictable Behavior**
- Always exactly 3 phases
- Always same exploration pattern
- Easier to debug and support

### 3. **Maintained Robustness**
- Still explores full parameter space
- Covers zero, positive, and negative bias regions
- Proven effective pattern

### 4. **Simpler Code**
- Less branching logic
- Fewer edge cases
- Easier to maintain

## Migration Path

### For Users with Existing Configurations

#### Old Configuration
```json
{
  "iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,
    "iterations_per_phase": 3,
    "max_phases": 4,
    "bias_shift_strategy": "positive_first",
    "bias_shift_magnitude": 30.0
  }
}
```

#### New Configuration (Automatic)
```json
{
  "iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,  // Preserved
    "iterations_per_phase": 3            // Preserved
    // Other parameters ignored/removed
  }
}
```

### Backward Compatibility
- Code will ignore deprecated parameters if present
- Won't break existing configurations
- Will log info about using fixed 3-phase approach

## Testing Requirements

### Unit Tests to Update
1. Test that exactly 3 phases are always used
2. Test bias shifts are always [0, +max_tol, -max_tol]
3. Remove tests for different strategies
4. Remove tests for custom bias magnitudes

### Integration Tests
1. Verify convergence with simplified system
2. Test with various scale_factor and iterations_per_phase values
3. Ensure backward compatibility with old configs

## Implementation Steps

1. **Update types.jl**
   - Remove fields from IterationSettings
   - Update constructor
   - Update parameter extraction

2. **Update ParameterTuningSearch.jl**
   - Simplify next_iteration!
   - Simplify calculate_phase_bias_shift
   - Update reset_for_new_phase!
   - Clean up logging

3. **Update Documentation**
   - Update GENERALIZED_ITERATION_PLAN.md
   - Update PARAMETER_TUNING_IMPROVEMENTS_SUMMARY.md
   - Remove or replace BIAS_SHIFT_STRATEGIES_EXPLAINED.md

4. **Test Changes**
   - Run existing tests
   - Add new tests for fixed behavior
   - Test backward compatibility

5. **Commit Changes**
   - Clear commit message about simplification
   - Note that functionality is preserved

## Summary

This simplification removes unnecessary configuration complexity while maintaining the robust parameter exploration that makes the system effective. The fixed 3-phase alternating pattern (zero → positive → negative) provides comprehensive coverage of the parameter space without requiring users to understand or configure bias shift strategies.