# ParameterTuningSearch Improvements Summary

## Overview

This document summarizes the comprehensive improvements made to the ParameterTuningSearch module in Pioneer.jl. These changes significantly enhance the robustness and flexibility of parameter tuning, particularly for challenging datasets with extreme mass biases or difficult-to-converge parameters.

## Major Changes Implemented

### 1. Configurable Parameters via JSON

All previously hardcoded values are now configurable through the JSON parameters file, providing users with full control over the parameter tuning process.

### 2. Enhanced Iteration Strategy

The iteration strategy has been completely redesigned with a phase-based approach that systematically explores different parameter spaces.

### 3. Improved Error Handling

Graceful fallback mechanisms ensure the pipeline continues even when parameter tuning fails for specific files.

## JSON Parameter Structure

### Complete Parameter Hierarchy

```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "intensity_filter_quantile": 0.25,  // NEW: Quantile for filtering fragments by intensity
      "max_tol_ppm": 50.0,                // NEW: Maximum allowed mass tolerance (moved from search_settings)
      // ... other fragment settings ...
    },
    "search_settings": {
      // Basic settings
      "initial_scan_count": 500,           // NEW: Starting number of scans
      "max_parameter_tuning_scans": 8000,  // RENAMED from expanded_scan_count
      "topn_peaks": 200,                   // NEW: Number of top peaks to consider
      "max_q_value": 0.01,                 // NEW: Q-value threshold for parameter tuning
      "max_frags_for_mass_err_estimation": 5, // RENAMED from max_best_rank
      
      // NEW: Iteration configuration (SIMPLIFIED 2025-01-14)
      "iteration_settings": {
        "mass_tolerance_scale_factor": 1.5,  // Still configurable
        "iterations_per_phase": 3            // Still configurable
        // Fixed 3-phase system with alternating ±max_tolerance shifts
      }
    }
  }
}
```

## Detailed Parameter Descriptions

### Fragment Settings

#### `intensity_filter_quantile` (NEW)
- **Type**: Float (0.0 - 1.0)
- **Default**: 0.25
- **Purpose**: Filters out low-intensity fragments during mass error model fitting
- **Impact**: Lower values (e.g., 0.1) keep only high-intensity fragments, potentially more accurate but fewer data points
- **Example**: 0.25 means only fragments above the 25th percentile intensity are used

#### `max_tol_ppm` (NEW - Moved from search_settings)
- **Type**: Float
- **Default**: 50.0
- **Purpose**: Maximum allowed mass tolerance regardless of scaling
- **Impact**: Prevents tolerance from expanding beyond reasonable limits
- **Safety**: Acts as a hard cap on tolerance expansion
- **Phase Shifts**: Used as the bias shift magnitude in Phase 2 (+max_tol_ppm) and Phase 3 (-max_tol_ppm)

### Search Settings

#### `initial_scan_count` (NEW)
- **Type**: Integer
- **Default**: 500
- **Purpose**: Number of scans to start with in the first iteration
- **Impact**: Lower values speed up initial iterations but may not find enough PSMs
- **Growth**: Doubles each iteration until `max_parameter_tuning_scans` is reached

#### `max_parameter_tuning_scans` (RENAMED)
- **Type**: Integer
- **Default**: 8000
- **Previous Name**: `expanded_scan_count`
- **Purpose**: Maximum number of scans to use in parameter tuning
- **Impact**: Caps the scan count to prevent excessive computation time
- **Behavior**: Scan count progression: 500 → 1000 → 2000 → 4000 → 8000 (capped)

#### `topn_peaks` (NEW)
- **Type**: Integer or null
- **Default**: null (no filtering)
- **Purpose**: Limits the number of peaks considered per spectrum
- **Impact**: Reduces computation time and noise in spectra with many peaks
- **Example**: 200 means only the 200 most intense peaks are considered

#### `max_q_value` (NEW)
- **Type**: Float
- **Default**: 0.01
- **Purpose**: Q-value threshold specifically for parameter tuning
- **Impact**: More stringent values (e.g., 0.001) require higher confidence PSMs
- **Independence**: Separate from the global q-value threshold for final results

#### `max_frags_for_mass_err_estimation` (RENAMED)
- **Type**: Integer
- **Default**: 5
- **Previous Name**: `max_best_rank`
- **Purpose**: Number of top fragments per precursor used for mass error estimation
- **Impact**: More fragments provide better statistics but may include noise
- **Clarity**: Renamed for better understanding of its purpose

### Iteration Settings (SIMPLIFIED 2025-01-14)

#### `mass_tolerance_scale_factor`
- **Type**: Float (must be > 1.0)
- **Default**: 2.0
- **Purpose**: Factor by which mass tolerance is scaled each iteration
- **Impact**: 
  - Conservative (1.2-1.5): Slower convergence, more precise
  - Moderate (1.5-2.0): Balanced speed and precision
  - Aggressive (2.0-3.0): Faster convergence, may overshoot
- **Example Progression** (factor=1.5): 20 ppm → 30 ppm → 45 ppm

#### `iterations_per_phase`
- **Type**: Integer
- **Default**: 3
- **Purpose**: Number of iterations before resetting to a new phase
- **Impact**: More iterations explore current bias region more thoroughly
- **Recommendation**: 2-4 iterations per phase is typically optimal

#### Fixed 3-Phase System (Not Configurable)
The system always uses exactly 3 phases with fixed bias shifts:
- **Phase 1**: Zero bias (explores around 0 ppm)
- **Phase 2**: Positive shift (+max_tol_ppm from fragment_settings)
- **Phase 3**: Negative shift (-max_tol_ppm from fragment_settings)

This fixed pattern ensures comprehensive parameter space coverage while maintaining simplicity.

## Configuration Examples

### Conservative Configuration (High-Accuracy Instruments)
```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "max_tol_ppm": 30.0
    },
    "search_settings": {
      "initial_scan_count": 1000,
      "max_parameter_tuning_scans": 5000,
      "iteration_settings": {
        "mass_tolerance_scale_factor": 1.2,
        "iterations_per_phase": 4
      }
    }
  }
}
```
**Use Case**: Modern Orbitrap instruments with stable calibration
**Note**: Will explore ±30 ppm bias regions with slow scaling

### Balanced Configuration (Default)
```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "max_tol_ppm": 50.0
    },
    "search_settings": {
      "initial_scan_count": 500,
      "max_parameter_tuning_scans": 8000,
      "iteration_settings": {
        "mass_tolerance_scale_factor": 2.0,
        "iterations_per_phase": 3
      }
    }
  }
}
```
**Use Case**: General purpose, works well for most datasets
**Note**: Will explore ±50 ppm bias regions with standard scaling

### Aggressive Configuration (Challenging Datasets)
```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "max_tol_ppm": 75.0
    },
    "search_settings": {
      "initial_scan_count": 250,
      "max_parameter_tuning_scans": 10000,
      "iteration_settings": {
        "mass_tolerance_scale_factor": 2.5,
        "iterations_per_phase": 2
      }
    }
  }
}
```
**Use Case**: Older instruments, poorly calibrated data, or unknown sample characteristics
**Note**: Will explore ±75 ppm bias regions with fast scaling

### Fast Configuration (Quick Processing)
```json
{
  "parameter_tuning": {
    "search_settings": {
      "initial_scan_count": 1000,
      "max_parameter_tuning_scans": 4000,
      "topn_peaks": 150,
      "iteration_settings": {
        "mass_tolerance_scale_factor": 3.0,
        "iterations_per_phase": 2
      }
    }
  }
}
```
**Use Case**: Large-scale processing where speed is prioritized
**Note**: Fewer iterations per phase for faster convergence

## Backward Compatibility

### Default Behavior
If no `iteration_settings` section is provided, the system uses defaults that match the previous behavior:
- Scale factor of 2.0 (doubling)
- 3 iterations per phase
- 3 phases maximum
- Alternating bias strategy

### Legacy Parameter Support
- `expanded_scan_count` is still supported and maps to `max_parameter_tuning_scans`
- `max_best_rank` is still supported and maps to `max_frags_for_mass_err_estimation`
- Missing parameters use sensible defaults

## How the Phase System Works

### Phase 1: Initial Exploration
```
Iteration 1: tolerance = 20 ppm, bias = 0 ppm
Iteration 2: tolerance = 30 ppm, bias = adjusted
Iteration 3: tolerance = 45 ppm, bias = adjusted
```

### Phase 2: Positive Bias Exploration
```
Reset: tolerance = 20 ppm, bias = +50 ppm
Iteration 4: tolerance = 20 ppm, bias = +50 ppm
Iteration 5: tolerance = 30 ppm, bias = adjusted
Iteration 6: tolerance = 45 ppm, bias = adjusted
```

### Phase 3: Negative Bias Exploration
```
Reset: tolerance = 20 ppm, bias = -50 ppm
Iteration 7: tolerance = 20 ppm, bias = -50 ppm
Iteration 8: tolerance = 30 ppm, bias = adjusted
Iteration 9: tolerance = 45 ppm, bias = adjusted
```

## Benefits of the New System

### 1. Robustness
- Handles extreme mass biases that were previously undetectable
- Systematic exploration ensures all reasonable parameter spaces are tested
- Graceful fallback mechanisms prevent pipeline failures

### 2. Flexibility
- Full control over convergence behavior through JSON configuration
- Different strategies for different instrument types and data quality
- Can be tuned for speed vs. accuracy trade-offs

### 3. Transparency
- Clear logging shows exactly what parameters are being tested
- Phase and iteration tracking provides insight into convergence process
- Diagnostic outputs help troubleshoot difficult files

### 4. Efficiency
- Early convergence stops unnecessary iterations
- Scan count management prevents excessive computation
- TopN peak filtering reduces noise in complex spectra

## Migration Guide

### From Old to New Configuration

#### Old Configuration (Minimal)
```json
{
  "parameter_tuning": {
    "search_settings": {
      "expanded_scan_count": 10000
    }
  }
}
```

#### New Configuration (Equivalent)
```json
{
  "parameter_tuning": {
    "search_settings": {
      "max_parameter_tuning_scans": 10000,
      "iteration_settings": {
        "mass_tolerance_scale_factor": 2.0,
        "iterations_per_phase": 3,
        "max_phases": 3
      }
    }
  }
}
```

### Recommended Migration Steps

1. **Start with defaults**: Run with no `iteration_settings` to use defaults
2. **Monitor convergence**: Check logs to see how many phases/iterations are needed
3. **Adjust if needed**: 
   - If converging too slowly: Increase `mass_tolerance_scale_factor`
   - If not converging: Increase `max_phases` or adjust `bias_shift_magnitude`
   - If taking too long: Reduce `iterations_per_phase` or `max_parameter_tuning_scans`

## Troubleshooting Guide

### Problem: Files not converging even after all phases

**Solutions**:
1. Increase `max_tol_ppm` in fragment_settings to allow wider search windows and larger bias shifts
2. Reduce `mass_tolerance_scale_factor` for finer exploration within each phase
3. Increase `iterations_per_phase` for more thorough exploration
4. Check data quality - may be inherently problematic

### Problem: Parameter tuning taking too long

**Solutions**:
1. Reduce `max_parameter_tuning_scans` 
2. Increase `initial_scan_count` to start with more data
3. Increase `mass_tolerance_scale_factor` for faster expansion
4. Enable `topn_peaks` filtering (e.g., 150-200)

### Problem: Poor mass accuracy after convergence

**Solutions**:
1. Decrease `mass_tolerance_scale_factor` for finer steps
2. Increase `iterations_per_phase` for more thorough exploration
3. Decrease `intensity_filter_quantile` to use only high-quality fragments
4. Increase `max_frags_for_mass_err_estimation` for better statistics

## Summary

These improvements transform ParameterTuningSearch from a rigid, sometimes fragile process into a robust, configurable system that can handle a wide variety of challenging datasets. The phase-based approach with configurable scaling provides a systematic way to explore parameter space, while maintaining backward compatibility and providing clear paths for optimization based on specific needs.