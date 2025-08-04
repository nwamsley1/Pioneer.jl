# Parameter Tuning Robustness Implementation Plan

## Executive Summary

This plan outlines a comprehensive approach to make ParameterTuningSearch more robust against edge cases and pathological scenarios. The implementation focuses on adaptive parameter adjustment, diagnostic tracking, and graceful degradation using sensible defaults when optimal parameters cannot be determined.

## Implementation Status (As of 2025-08-02)

### Current Progress
- ✅ **Phase 1: Graceful Fallback System** - COMPLETED
  - Implemented fallback to conservative defaults (±50 ppm, identity RT)
  - Added warning system with file-specific messages
  - Basic diagnostic structure in place
  - **Issue Found**: Type initialization prevents fallback from working (needs fix)
  
- ✅ **Phase 2: Mass Bias Detection Strategy** - COMPLETED
  - Created `bias_detection.jl` with BiasSearchStrategy
  - Implemented intelligent bias search up to ±50 ppm
  - Added adaptive window expansion for bias
  - **Issue Found**: Yields too few PSMs in simple samples (41 vs 100 needed)
  
- ✅ **Phase 3: Diagnostic Infrastructure** - COMPLETED
  - Created `diagnostics.jl` with comprehensive tracking
  - Added ParameterTuningDiagnostics and ParameterTuningStatus types
  - Integrated diagnostic reporting into main flow
  
- ✅ **Boundary Sampling Validation** - COMPLETED
  - Created `boundary_sampling.jl` with validation logic
  - Implemented 5% threshold checking at tolerance boundaries
  - Added automatic tolerance expansion when needed
  
- ✅ **Cross-Run Parameter Learning** - COMPLETED
  - Created `cross_run_learning.jl` with parameter history tracking
  - Implemented intelligent initial parameter selection
  - Added adaptive strategy selection based on variability

- ✅ **Module Loading Issues** - RESOLVED
  - Created `types.jl` to break circular dependencies
  - Fixed loading order in `importScripts.jl`
  - All types now load correctly

### Critical Issues Discovered (2025-08-02)

1. **Type System Bug**: 
   - `rt_to_irt_model` initialized as `Ref{SplineRtConversionModel}()` instead of `Ref{RtConversionModel}()`
   - Prevents fallback to `IdentityModel`, causing crashes
   - **Location**: Line 107 in `ParameterTuningSearch.jl`
   - **Fix Required**: Change initialization to use abstract type

2. **PSM Threshold Too High**:
   - Current minimum of 1000 PSMs unrealistic for simple samples
   - E. coli test only yields 281-350 PSMs per iteration
   - Bias detection only finds 41 PSMs (needs 100)
   - **Fix Required**: Implement adaptive thresholds based on sample complexity

### Immediate Fixes Required

#### Fix 1: Type Initialization Bug (CRITICAL)
```julia
# In init_search_results(), line 107:
# WRONG:
Ref{SplineRtConversionModel}(),

# CORRECT:
Ref{RtConversionModel}(),
```

This allows the field to hold either `SplineRtConversionModel` or `IdentityModel`.

#### Fix 2: Update types.jl Structure Definition
```julia
# In types.jl, line 113:
# WRONG:
rt_to_irt_model::Base.Ref{<:RtConversionModel}

# CORRECT:
rt_to_irt_model::Base.Ref{RtConversionModel}
```

Remove the `<:` to match the initialization.

### Next Steps
1. ✅ Fix type initialization bug in `init_search_results` 
2. ✅ Update struct definition in `types.jl`
3. Implement adaptive PSM thresholds (new Phase 1.3)
4. Enhance low-PSM handling strategies
5. Add sample complexity estimation
6. Create fallback parameter profiles
7. Continue with Phases 4-5 after fixes

## Problem Analysis

### Current Limitations

1. **Hard Failure Mode**: Files that don't converge are marked as failed with no recovery
2. **Fixed Sample Rate**: The 2% sampling rate may be insufficient for sparse data or overly conservative for dense data
3. **Simple Convergence Logic**: Only checks mass offset and total tolerance, missing other failure modes
4. **No TopN Filtering**: When tolerance is large or isolation windows are wide, the search becomes prohibitively slow
5. **Limited Diagnostics**: No tracking of convergence metrics or failure patterns

### Key Concepts: Mass Bias vs Mass Tolerance

**Mass Bias (Accuracy)**: Systematic offset in measured masses
- Example: All masses are 18 ppm lower than expected
- Corrected by `mass_offset` in MassErrorModel
- Applied as: `corrected_mz = observed_mz - (mass_offset * observed_mz / 1e6)`

**Mass Tolerance (Precision)**: Random variation around the bias-corrected value
- Example: ±20 ppm window after bias correction
- Defined by `(left_tolerance, right_tolerance)` in MassErrorModel
- Creates search window: `[theoretical - left_tol, theoretical + right_tol]`

### Critical Edge Cases

#### Case 1: Large Mass Bias
- **Scenario**: True tolerance is ±20 ppm but bias is -18 ppm
- **Problem**: Initial search window (e.g., ±20 ppm around 0 bias) misses most peaks
- **Current Behavior**: Might converge slowly or fail entirely
- **Needed**: Intelligent bias detection and search window adjustment

#### Case 2: Ultra-Wide Isolation Windows
- **Scenario**: DIA with 25+ m/z isolation windows
- **Problem**: Too many precursor candidates per spectrum
- **Current Behavior**: Extremely slow due to exhaustive searching
- **Needed**: TopN peak selection and/or stricter scoring thresholds

#### Case 3: Low Complexity Samples
- **Scenario**: Few peptides present (e.g., purified proteins, QC samples)
- **Problem**: Cannot collect enough PSMs for reliable fitting
- **Current Behavior**: Fails after max iterations
- **Needed**: Graceful fallback to conservative defaults

#### Case 4: Inadequate Tolerance Boundary Sampling
- **Scenario**: Initial tolerance guess too close to true tolerance (e.g., 20 ppm guess, 22 ppm true)
- **Problem**: Insufficient sampling at distribution tails for accurate tolerance estimation
- **Current Behavior**: May underestimate true tolerance, leading to missed matches
- **Needed**: Boundary density check - ensure <5% of matches are within 10% of tolerance boundaries. If violated, expand tolerance and re-search.

#### Case 5: Excessive Initial Tolerance
- **Scenario**: Initial tolerance much larger than true tolerance (e.g., 50 ppm guess, 10 ppm true)
- **Problem**: Too many false positive matches in initial search, poor discrimination
- **Current Behavior**: Slow search, potential for incorrect parameter estimation
- **Needed**: Iterative tolerance reduction or TopN filtering from the start

#### Case 6: Variable Mass Error Across Runs
- **Scenario**: Different MS files have different mass calibrations
- **Problem**: Each file needs independent parameter fitting
- **Current Behavior**: Handles this correctly (per-file fitting)
- **Enhancement**: Better cross-file diagnostic reporting

## Implementation Plan

### Phase 1: Graceful Fallback System (Priority 1)

#### 1.1 Remove Hard Failure Mode ✅ COMPLETED
```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    try
        # Existing parameter tuning logic
        # ...
    catch e
        @warn "Parameter tuning failed for file $ms_file_idx. Using conservative defaults." exception=e
        
        # Set conservative defaults
        results.mass_err_model[] = MassErrorModel(
            0.0f0,  # No bias assumption
            (50.0f0, 50.0f0)  # ±50 ppm tolerance
        )
        
        # Simple RT model or identity mapping
        results.rt_to_irt_model[] = IdentityRTModel()
        
        # Mark in diagnostics but don't fail
        setWarningFlag!(getMSData(search_context), ms_file_idx, 
                       "PARAM_TUNING_FALLBACK")
        
        # Store diagnostic info for user
        storeDiagnostics!(search_context, ms_file_idx, 
                         "Used conservative defaults due to tuning failure")
    end
    
    return results
end
```

#### 1.2 Warning System ✅ COMPLETED
```julia
# Add to SearchContext or MSData
struct ParameterTuningWarnings
    file_warnings::Dict{Int, Vector{String}}
    global_warnings::Vector{String}
end

# Easy detection for users
function hasParameterTuningWarnings(search_context::SearchContext)
    return !isempty(search_context.param_warnings.file_warnings)
end

function getParameterTuningReport(search_context::SearchContext)
    # Generate markdown report of all warnings and fallbacks
end
```

#### 1.3 Adaptive PSM Thresholds (NEW - Priority 1)

Current issue: Fixed 1000 PSM minimum is too high for simple samples like E. coli test data.

```julia
struct AdaptivePSMThresholds
    # Minimum thresholds by sample complexity
    high_complexity_min::Int      # Default: 1000 (human plasma, tissue)
    medium_complexity_min::Int    # Default: 500 (cell lysates)
    low_complexity_min::Int       # Default: 200 (purified proteins, QC)
    
    # Absolute minimums for fitting
    mass_model_min::Int          # Default: 100 (minimum for reliable mass fitting)
    rt_model_min::Int            # Default: 50 (minimum for RT spline)
    
    # Fallback thresholds
    bias_detection_min::Int      # Default: 50 (minimum for bias search)
end

function estimate_sample_complexity(
    initial_psm_count::Int,
    n_unique_peptides::Int,
    n_proteins::Int,
    isolation_window_width::Float32
) :: Symbol
    # Heuristic based on initial results
    if n_proteins > 1000 || n_unique_peptides > 5000
        return :high_complexity
    elseif n_proteins > 100 || n_unique_peptides > 500
        return :medium_complexity
    else
        return :low_complexity
    end
end

function get_adaptive_psm_threshold(
    complexity::Symbol,
    thresholds::AdaptivePSMThresholds
) :: Int
    if complexity == :high_complexity
        return thresholds.high_complexity_min
    elseif complexity == :medium_complexity
        return thresholds.medium_complexity_min
    else
        return thresholds.low_complexity_min
    end
end
```

**Fallback Profiles for Low-PSM Scenarios:**
```julia
struct FallbackProfile
    name::String
    min_psms::Int
    mass_tolerance::Tuple{Float32, Float32}
    rt_model_type::Symbol  # :spline, :linear, :identity
    confidence_level::Symbol  # :high, :medium, :low
end

const FALLBACK_PROFILES = [
    FallbackProfile("Conservative", 0, (50.0f0, 50.0f0), :identity, :low),
    FallbackProfile("Low Sample", 50, (30.0f0, 30.0f0), :linear, :medium),
    FallbackProfile("Medium Sample", 200, (25.0f0, 25.0f0), :spline, :medium),
    FallbackProfile("Normal", 500, (20.0f0, 20.0f0), :spline, :high)
]

function select_fallback_profile(psm_count::Int) :: FallbackProfile
    for profile in reverse(FALLBACK_PROFILES)
        if psm_count >= profile.min_psms
            return profile
        end
    end
    return FALLBACK_PROFILES[1]  # Most conservative
end
```

### Phase 2: Cross-Run Parameter Learning (Priority 1)

#### 2.1 Parameter History Tracking
```julia
# Add to SearchContext
struct ParameterHistory
    file_parameters::Dict{Int, TuningResults}  # Per-file results
    global_stats::GlobalParameterStats
end

struct TuningResults
    mass_offset::Float32
    mass_tolerance::Tuple{Float32, Float32}
    converged::Bool
    psm_count::Int
    warnings::Vector{String}
end

struct GlobalParameterStats
    median_mass_offset::Float32
    mass_offset_mad::Float32
    median_tolerance::Float32
    tolerance_mad::Float32
    n_successful_files::Int
end
```

#### 2.2 Intelligent Initial Parameter Selection
```julia
function get_initial_parameters(search_context::SearchContext, ms_file_idx::Int)
    history = search_context.parameter_history
    
    # First file: use config defaults
    if isempty(history.file_parameters)
        return get_default_parameters()
    end
    
    # Subsequent files: use learned parameters
    stats = history.global_stats
    
    # Start with median values from previous successful runs
    initial_bias = stats.median_mass_offset
    
    # Add safety margin based on observed variability
    bias_uncertainty = 2.0 * stats.mass_offset_mad
    tolerance_buffer = 1.5 * stats.median_tolerance + 2.0 * stats.tolerance_mad
    
    # For bias search, expand range based on observed variability
    bias_search_range = max(50.0, 3.0 * stats.mass_offset_mad)
    
    return InitialParameters(
        bias_estimate = initial_bias,
        bias_search_range = bias_search_range,
        initial_tolerance = tolerance_buffer,
        informed_by_history = true
    )
end
```

#### 2.3 Adaptive Strategy Selection
```julia
function select_search_strategy(history::ParameterHistory)
    # If all previous files had similar parameters, use focused search
    if history.global_stats.mass_offset_mad < 5.0 &&
       history.global_stats.tolerance_mad < 3.0
        return FocusedSearchStrategy()  # Narrow search around expected values
    else
        return BroadSearchStrategy()    # Wide search for variable data
    end
end
```

### Phase 3: Mass Bias Detection Strategy (Priority 1)

#### 2.1 Intelligent Bias Search
```julia
struct BiasSearchStrategy
    initial_tolerance::Float32  # From config (e.g., 20 ppm)
    bias_search_range::Float32  # How far to search for bias (e.g., ±50 ppm)
    bias_search_steps::Vector{Float32}  # Steps to try: [0, -10, 10, -20, 20, -30, 30, ...]
end

function detect_mass_bias(spectra, search_context, params, ms_file_idx)
    strategy = BiasSearchStrategy(
        getFragTolPpm(params),
        50.0f0,  # Search up to ±50 ppm bias
        [0.0, -10.0, 10.0, -20.0, 20.0, -30.0, 30.0, -40.0, 40.0, -50.0, 50.0]
    )
    
    best_psm_count = 0
    best_bias = 0.0f0
    
    for bias_guess in strategy.bias_search_steps
        # Quick search with bias offset
        temp_mass_model = MassErrorModel(
            bias_guess,
            (strategy.initial_tolerance, strategy.initial_tolerance)
        )
        
        # Do limited search (fewer spectra) to estimate PSM yield
        psm_count = estimate_psm_yield(spectra, search_context, temp_mass_model, params)
        
        if psm_count > best_psm_count
            best_psm_count = psm_count
            best_bias = bias_guess
        end
        
        # Early termination if found good yield
        if psm_count > getMinPsms(params) * 0.5
            break
        end
    end
    
    return best_bias, best_psm_count
end
```

#### 2.2 Adaptive Window Expansion for Bias
When large bias is detected, the search window needs special handling:

```julia
function adjust_search_window_for_bias(initial_tolerance, detected_bias)
    # If bias is large relative to tolerance, expand asymmetrically
    if abs(detected_bias) > initial_tolerance * 0.5
        if detected_bias > 0
            # Positive bias: expand right side more
            left_tol = initial_tolerance
            right_tol = initial_tolerance + abs(detected_bias) * 0.5
        else
            # Negative bias: expand left side more
            left_tol = initial_tolerance + abs(detected_bias) * 0.5
            right_tol = initial_tolerance
        end
    else
        # Small bias: symmetric expansion
        left_tol = right_tol = initial_tolerance
    end
    
    return (left_tol, right_tol)
end
```

#### 2.3 Boundary Sampling Validation
```julia
struct BoundarySamplingCheck
    boundary_fraction_threshold::Float32  # e.g., 0.05 (5%)
    boundary_zone_fraction::Float32       # e.g., 0.1 (10% of tolerance)
end

function check_boundary_sampling(ppm_errors::Vector{Float32}, 
                               mass_model::MassErrorModel,
                               checker::BoundarySamplingCheck)
    # Remove bias to analyze tolerance boundaries
    centered_errors = ppm_errors .- getMassOffset(mass_model)
    left_tol, right_tol = getLeftTol(mass_model), getRightTol(mass_model)
    
    # Define boundary zones (outer 10% of tolerance range)
    left_boundary = -left_tol * (1 - checker.boundary_zone_fraction)
    right_boundary = right_tol * (1 - checker.boundary_zone_fraction)
    
    # Count matches in boundary zones
    n_left_boundary = count(err -> err < left_boundary, centered_errors)
    n_right_boundary = count(err -> err > right_boundary, centered_errors)
    n_total = length(centered_errors)
    
    # Calculate fractions
    left_fraction = n_left_boundary / n_total
    right_fraction = n_right_boundary / n_total
    
    # Check if we have adequate sampling
    adequate_left = left_fraction < checker.boundary_fraction_threshold
    adequate_right = right_fraction < checker.boundary_fraction_threshold
    
    return BoundarySamplingResult(
        adequate_sampling = adequate_left && adequate_right,
        left_boundary_fraction = left_fraction,
        right_boundary_fraction = right_fraction,
        expansion_needed = !adequate_left || !adequate_right,
        suggested_expansion_factor = max(
            adequate_left ? 1.0 : 1.5,
            adequate_right ? 1.0 : 1.5
        )
    )
end
```

### Phase 3: Diagnostic Infrastructure (Priority 2)

#### 3.1 Comprehensive Diagnostics
```julia
struct ParameterTuningDiagnostics
    # Bias detection phase
    bias_search_attempted::Bool
    bias_search_results::Vector{Tuple{Float32, Int}}  # (bias, psm_count) pairs
    initial_bias_estimate::Float32
    
    # Convergence tracking
    iteration_count::Int
    iteration_metrics::Vector{IterationMetrics}
    converged::Bool
    convergence_reason::String
    
    # Final parameters
    final_mass_offset::Float32
    final_mass_tolerance::Tuple{Float32, Float32}
    final_psm_count::Int
    used_fallback::Bool
    fallback_reason::String
    
    # Warnings
    warnings::Vector{String}
end

struct IterationMetrics
    psm_count::Int
    mass_offset::Float32
    mass_tolerance::Tuple{Float32, Float32}
    sample_rate::Float32
    search_time::Float64
end
```

#### 3.2 Diagnostic Reporting
Generate clear report for each file:
```
## Parameter Tuning Report - File: sample_001.raw

### Status: ⚠️ Converged with warnings

### Mass Calibration
- Initial bias estimate: -18.5 ppm (detected via bias search)
- Final mass offset: -18.2 ppm
- Final tolerance: ±22.3 ppm
- Search window adjustments: 2

### Convergence
- Iterations: 4
- Final PSM count: 2,847
- Convergence reason: Stable parameters for 2 iterations

### Warnings
- Large mass bias detected (-18.5 ppm)
- Initial search window adjusted for bias compensation

### Recommendations
- Consider mass recalibration for this instrument
- Current parameters should work but efficiency is reduced
```

### Phase 4: Adaptive Sampling Strategy (Priority 2)

#### 4.1 Dynamic Sample Rate
```julia
struct AdaptiveSamplingState
    base_sample_rate::Float32
    current_sample_rate::Float32
    min_psms_per_iteration::Int
    
    function adjust_sample_rate!(state, current_psm_count, target_psm_count)
        if current_psm_count < state.min_psms_per_iteration
            # Increase sampling
            state.current_sample_rate = min(
                state.current_sample_rate * 2.0,
                1.0  # Cap at 100%
            )
            return "Increased sample rate to $(state.current_sample_rate)"
        elseif current_psm_count > target_psm_count * 2
            # Decrease sampling if we have plenty
            state.current_sample_rate = max(
                state.current_sample_rate * 0.7,
                state.base_sample_rate
            )
            return "Decreased sample rate to $(state.current_sample_rate)"
        end
        return nothing
    end
end
```

#### 4.2 TopN Peak Filtering
```julia
struct TopNPeakFilter
    enabled::Bool
    n_peaks::Int
    activation_conditions::TopNActivation
end

struct TopNActivation
    when_sample_rate_above::Float32  # e.g., 0.1
    when_search_time_above::Float64   # e.g., 5.0 seconds
    when_isolation_width_above::Float32  # e.g., 20 m/z
end

function should_use_topn(filter::TopNPeakFilter, state::SearchState)
    return filter.enabled && (
        state.current_sample_rate > filter.activation_conditions.when_sample_rate_above ||
        state.last_search_time > filter.activation_conditions.when_search_time_above ||
        state.avg_isolation_width > filter.activation_conditions.when_isolation_width_above
    )
end
```

### Phase 5: Enhanced Convergence Logic (Priority 3)

#### 5.1 Multi-Criteria Convergence
```julia
struct ConvergenceChecker
    # Stability thresholds
    mass_offset_stability_ppm::Float32  # e.g., 2.0 ppm
    tolerance_stability_factor::Float32  # e.g., 0.1 (10% change)
    min_stable_iterations::Int  # e.g., 2
    
    # Quality thresholds
    min_psm_count::Int
    max_iterations::Int
    
    function check_convergence(checker, history::Vector{IterationMetrics})
        if length(history) < checker.min_stable_iterations
            return false, "Insufficient iterations"
        end
        
        # Check mass offset stability
        recent = history[end-checker.min_stable_iterations+1:end]
        offset_changes = [abs(recent[i].mass_offset - recent[i-1].mass_offset) 
                         for i in 2:length(recent)]
        offset_stable = all(change < checker.mass_offset_stability_ppm 
                           for change in offset_changes)
        
        # Check tolerance stability
        tol_changes = [max(
            abs(recent[i].mass_tolerance[1] - recent[i-1].mass_tolerance[1]),
            abs(recent[i].mass_tolerance[2] - recent[i-1].mass_tolerance[2])
        ) / mean(recent[i].mass_tolerance) for i in 2:length(recent)]
        
        tolerance_stable = all(change < checker.tolerance_stability_factor 
                              for change in tol_changes)
        
        # Check PSM count
        sufficient_psms = recent[end].psm_count >= checker.min_psm_count
        
        if offset_stable && tolerance_stable && sufficient_psms
            return true, "Parameters stable"
        elseif length(history) >= checker.max_iterations
            return true, "Maximum iterations reached"
        else
            reasons = String[]
            !offset_stable && push!(reasons, "Mass offset unstable")
            !tolerance_stable && push!(reasons, "Tolerance unstable")
            !sufficient_psms && push!(reasons, "Insufficient PSMs")
            return false, join(reasons, ", ")
        end
    end
end
```

### Phase 6: Configuration Updates (Priority 3)

```json
{
    "parameter_tuning": {
        "fallback_strategy": {
            "enabled": true,
            "default_mass_tolerance_ppm": 50.0,
            "default_rt_model": "identity",
            "warn_on_fallback": true
        },
        "bias_detection": {
            "enabled": true,
            "search_range_ppm": 50.0,
            "quick_search_sample_rate": 0.005,
            "min_psms_for_detection": 100
        },
        "adaptive_psm_thresholds": {
            "enabled": true,
            "high_complexity_min": 1000,
            "medium_complexity_min": 500,
            "low_complexity_min": 200,
            "mass_model_min": 100,
            "rt_model_min": 50,
            "bias_detection_min": 50
        },
        "adaptive_sampling": {
            "enabled": true,
            "min_sample_rate": 0.02,
            "max_sample_rate": 1.0,
            "rate_multiplier": 2.0,
            "min_psms_per_iteration": 500
        },
        "topn_filtering": {
            "enabled": true,
            "n_peaks": 20,
            "activate_above_sample_rate": 0.1,
            "activate_above_search_time_seconds": 5.0,
            "activate_above_isolation_width_mz": 20.0
        },
        "convergence": {
            "mass_offset_stability_ppm": 2.0,
            "tolerance_stability_factor": 0.1,
            "min_stable_iterations": 2,
            "max_iterations": 10
        },
        "boundary_sampling": {
            "enabled": true,
            "boundary_fraction_threshold": 0.05,
            "boundary_zone_fraction": 0.1,
            "auto_expand_on_violation": true
        },
        "cross_run_learning": {
            "enabled": true,
            "min_files_for_learning": 2,
            "parameter_variability_threshold": 5.0
        }
    }
}
```

## Updated Convergence Loop Integration

The enhanced `process_file!` function would incorporate all these improvements:

```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    diagnostics = ParameterTuningDiagnostics()
    
    try
        # 1. Get initial parameters from cross-run learning
        initial_params = get_initial_parameters(search_context, ms_file_idx)
        
        # 2. Attempt bias detection if needed
        if should_attempt_bias_detection(initial_params)
            best_bias, psm_count = detect_mass_bias(spectra, search_context, params, ms_file_idx)
            initial_params.bias_estimate = best_bias
        end
        
        # 3. Main convergence loop with enhancements
        n_attempts = 0
        sampling_state = AdaptiveSamplingState(params)
        
        while n_attempts < params.max_iterations
            # Set current models
            setMassErrorModel!(search_context, ms_file_idx, current_mass_model)
            
            # Collect PSMs with adaptive sampling
            psms = collect_psms_adaptive(spectra, search_context, params, 
                                        ms_file_idx, sampling_state)
            
            # Fit models
            rt_model = fit_irt_model(params, psms)
            fragments = get_matched_fragments(spectra, psms, results, search_context, 
                                            params, ms_file_idx)
            mass_err_model, ppm_errs = fit_mass_err_model(params, fragments)
            
            # Check boundary sampling
            boundary_result = check_boundary_sampling(ppm_errs, mass_err_model, params)
            if !boundary_result.adequate_sampling
                # Expand tolerance and retry
                expand_tolerance!(current_mass_model, boundary_result.suggested_expansion_factor)
                diagnostics.warnings.push!("Inadequate boundary sampling - expanding tolerance")
                continue
            end
            
            # Check convergence
            converged, reason = check_convergence(convergence_checker, iteration_history)
            if converged
                break
            end
            
            # Adjust sampling if needed
            adjust_sample_rate!(sampling_state, length(psms), params.min_psms)
            
            n_attempts += 1
        end
        
        # 4. Store results for cross-run learning
        store_tuning_results!(search_context, ms_file_idx, final_params)
        update_global_statistics!(search_context.parameter_history)
        
    catch e
        # 5. Graceful fallback
        @warn "Parameter tuning failed for file $ms_file_idx. Using conservative defaults." exception=e
        apply_fallback_parameters!(results, params)
        diagnostics.used_fallback = true
        diagnostics.fallback_reason = string(e)
    end
    
    # 6. Generate diagnostics
    generate_diagnostic_report(diagnostics, ms_file_idx)
    
    return results
end
```

## Implementation Priority

### Phase 1 (Week 1) - Critical ✅ COMPLETED
1. ✅ Implement graceful fallback (remove hard failures)
2. ✅ Add warning system with clear user visibility
3. ✅ Basic diagnostic structure
4. 🔧 **IN PROGRESS**: Adaptive PSM thresholds (Phase 1.3)
   - Need to fix type initialization bug first
   - Then implement adaptive thresholds

### Phase 2 (Week 1-2) - High Priority ✅ COMPLETED
1. ✅ Mass bias detection strategy
2. ✅ Adaptive search window for large bias
3. ✅ Integration with existing convergence loop
4. ✅ Boundary sampling validation
5. ✅ Cross-run parameter learning

### Phase 3 (Week 2) - Medium Priority ✅ COMPLETED
1. ✅ Comprehensive diagnostics
2. ✅ Per-file diagnostic reports
3. ✅ Cross-file summary reports

### Critical Fixes Required (Immediate)
1. 🔴 **Type initialization bug** - Prevents fallback from working
2. 🔴 **Adaptive PSM thresholds** - Current thresholds too high for simple samples

### Phase 4 (Week 3) - Medium Priority ⏳ PENDING
1. Adaptive sampling rate
2. TopN peak filtering
3. Performance optimizations

### Phase 5 (Week 3-4) - Lower Priority ⏳ PENDING
1. Enhanced convergence criteria
2. Configuration system
3. Full testing suite

## Success Metrics

1. **Zero Hard Failures**: All files produce results (either fitted or fallback)
2. **Clear Warnings**: Users always informed when fallbacks used
3. **Bias Handling**: Successfully handle bias up to ±50 ppm
4. **Performance**: <2x slowdown for typical cases
5. **Diagnostic Quality**: Clear actionable reports for problem files

## Testing Strategy

### Unit Tests
- Test bias detection with synthetic data
- Verify fallback mechanisms
- Test convergence detection
- Validate adaptive sampling

### Integration Tests  
- Known difficult datasets:
  - Large mass bias samples
  - Wide isolation windows
  - Low complexity samples
- Verify downstream methods handle all parameter combinations

### Edge Case Tests
- Zero PSMs found
- Extreme mass bias (>50 ppm)
- Very noisy data
- Single peptide samples

## Questions for Review

1. Is 50 ppm a reasonable default fallback tolerance?
2. Should bias detection be enabled by default?
3. What's the preferred warning mechanism (log, file, both)?
4. Should we store raw diagnostic data or just summaries?
5. Any preference on diagnostic report format?
6. Should TopN filtering be conservative (off by default)?