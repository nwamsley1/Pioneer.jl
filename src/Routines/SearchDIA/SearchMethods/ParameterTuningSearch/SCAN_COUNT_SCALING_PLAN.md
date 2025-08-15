# Plan: Simplified Scan Count Scaling Strategy

## Status: PENDING APPROVAL

## Current Behavior

Currently:
1. Scan count starts at `initial_scan_count` (e.g., 500)
2. Doubles each iteration within a phase (500 → 1000 → 2000 → 4000)
3. Capped at `max_parameter_tuning_scans` (e.g., 8000)
4. Resets to initial when entering a new phase

## Proposed Behavior

New approach:
1. Run all 3 phases with `initial_scan_count`
2. If no convergence, multiply scan count by a scaling factor (e.g., 2.0)
3. Run all 3 phases again with the new scan count
4. Repeat until scan count reaches or exceeds `max_parameter_tuning_scans`
5. Run one final attempt with the maximum scan count
6. Give up if still no convergence

## Example Flow

With `initial_scan_count=500`, `max_parameter_tuning_scans=8000`, `scale_factor=2.0`:

```
Attempt 1: Run 3 phases with 500 scans
  Phase 1: 500 scans (3 iterations)
  Phase 2: 500 scans (3 iterations)  
  Phase 3: 500 scans (3 iterations)
  → No convergence

Attempt 2: Run 3 phases with 1000 scans (500 * 2)
  Phase 1: 1000 scans (3 iterations)
  Phase 2: 1000 scans (3 iterations)
  Phase 3: 1000 scans (3 iterations)
  → No convergence

Attempt 3: Run 3 phases with 2000 scans (1000 * 2)
  Phase 1: 2000 scans (3 iterations)
  Phase 2: 2000 scans (3 iterations)
  Phase 3: 2000 scans (3 iterations)
  → No convergence

Attempt 4: Run 3 phases with 4000 scans (2000 * 2)
  Phase 1: 4000 scans (3 iterations)
  Phase 2: 4000 scans (3 iterations)
  Phase 3: 4000 scans (3 iterations)
  → No convergence

Attempt 5: Run 3 phases with 8000 scans (capped at max)
  Phase 1: 8000 scans (3 iterations)
  Phase 2: 8000 scans (3 iterations)
  Phase 3: 8000 scans (3 iterations)
  → No convergence → GIVE UP
```

## Implementation Changes

### 1. Modify IterationState

Add fields to track scan scaling:
```julia
mutable struct IterationState
    # Existing fields...
    current_phase::Int64
    current_iteration_in_phase::Int64
    total_iterations::Int64
    phase_bias_shifts::Vector{Float32}
    converged::Bool
    collection_tolerance::Float32
    
    # New fields for scan scaling
    current_scan_count::Int64        # Current scan count for this attempt
    scan_attempt::Int64              # Which scaling attempt we're on (1, 2, 3...)
end
```

### 2. New Function: `run_phases_with_scan_count`

```julia
function run_phases_with_scan_count(
    scan_count::Int64,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # Run all 3 phases with fixed scan count
    for phase in 1:3
        # Run phase with fixed scan count
        converged = run_single_phase(phase, scan_count, ...)
        if converged
            return true
        end
    end
    return false
end
```

### 3. Simplified process_file!

```julia
function process_file!(...)
    scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    scale_factor = 2.0  # Hardcoded or configurable
    
    while scan_count <= max_scans
        @info "Attempting parameter tuning with $scan_count scans"
        
        # Run all 3 phases with current scan count
        converged = run_phases_with_scan_count(scan_count, ...)
        
        if converged
            return true
        end
        
        # Scale up for next attempt
        scan_count = min(scan_count * scale_factor, max_scans)
        
        # If we've reached max, do one final attempt
        if scan_count >= max_scans && scan_count != max_scans
            scan_count = max_scans
        elseif scan_count == max_scans
            # Already tried with max, give up
            break
        end
    end
    
    # Failed to converge
    apply_fallback_strategy!(...)
end
```

### 4. Remove Per-Iteration Scan Scaling

Remove the current logic that doubles scans each iteration:
- Remove `add_scans_for_iteration!` 
- Replace with fixed scan count per phase attempt

## Benefits

1. **Simpler Logic**: Each phase runs with consistent scan count
2. **More Predictable**: Clear progression of scan counts
3. **Better Resource Usage**: Don't waste time with many low-scan iterations
4. **Clearer Failure Mode**: Know exactly when we've exhausted options

## Configuration

Keep existing parameters:
- `initial_scan_count`: Starting point
- `max_parameter_tuning_scans`: Maximum allowed

Add (optional, could be hardcoded):
- `scan_scale_factor`: How much to multiply scan count between attempts (default: 2.0)

## Questions for Approval

1. Should the scale factor be configurable or hardcoded to 2.0?
2. Should we run all 3 phases before scaling, or allow early scaling if a phase shows no promise?
3. Should we keep the per-iteration mass tolerance scaling within each phase?
4. Any preference on the exact scaling sequence?