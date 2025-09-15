# Two-Tier FTR Control Strategies Implementation Plan

## Overview

This document outlines the implementation plan for flexible False Transfer Rate (FTR) control strategies for match-between-runs analysis. The system provides three different approaches: single-tier experiment-wide control, or two different two-tier hierarchical strategies.

## FTR Control Strategies

### Strategy A: "experiment_wide" (Single-Tier, Original)
```
Single Stage: Experiment-wide FTR control across all transfer candidates
```

### Strategy B: "experiment_local" (Two-Tier, Current Implementation)
```
Stage 1: Experiment-wide FTR control → Stage 2: Local per-file refinement
```

### Strategy C: "global_experiment" (Two-Tier, New Implementation)
```  
Stage 1: Global precursor-level FTR control → Stage 2: Experiment-wide FTR control
```

## Technical Design

### Strategy A: "experiment_wide" (Single-Tier)
**Algorithm**:
1. Apply existing `get_ftr_threshold()` across all transfer candidates
2. Single threshold used experiment-wide
3. **Effect**: Simple, fast FTR control - original behavior

### Strategy B: "experiment_local" (Current Two-Tier)
**Stage 1 - Experiment-wide**: 
1. Apply `get_ftr_threshold()` across all transfer candidates globally
2. Establish baseline threshold

**Stage 2 - Local per-file**:
1. Group candidates by `ms_file_idx` 
2. Apply per-file refinement using `get_hierarchical_ftr_thresholds()`
3. Never less stringent than Stage 1 threshold
4. **Effect**: File-specific adaptation with global floor

### Strategy C: "global_experiment" (New Two-Tier)
**Stage 1 - Global precursor-level**:
1. Group all MBR transfer candidates by `(precursor_idx, isotopes_captured)`
2. For each precursor group, select the **maximum probability** representative among transfer candidates
3. Apply FDR calculation using these precursor representatives:
   - Targets: precursor representatives that are targets  
   - Transfer decoys: precursor representatives where best cross-run match is a decoy
4. Filter: retain precursor groups where representative passes q-value threshold
5. **Effect**: Eliminates entire precursor groups with poor cross-run evidence

**Stage 2 - Experiment-wide**:
1. Take all transfer candidates from precursor groups that passed Stage 1
2. Apply existing `get_ftr_threshold()` across remaining candidates
3. **Effect**: Standard experiment-wide FTR control on the pre-filtered candidate set

## Configuration Options

### "experiment_wide" 
- Single-tier experiment-wide FTR control
- Original behavior, fastest
- Good baseline for most datasets

### "experiment_local" 
- Two-tier: Experiment-wide → Local per-file  
- Current hierarchical implementation
- Better per-file adaptation

### "global_experiment"
- Two-tier: Global precursor-level → Experiment-wide
- New implementation
- Better handling of poor precursor groups

## Implementation Details

### Configuration Parameter
```json
"global": {
    "scoring": {
        "q_value_threshold": 0.01,
        "mbr_ftr_control": "experiment_local",  // "experiment_wide", "experiment_local", "global_experiment"
        "hierarchical_ftr_min_candidates": 50
    }
}
```

### Core Data Structures
```julia
# Precursor representative for global FTR
struct PrecursorRepresentative
    precursor_idx::UInt32
    isotopes_captured::Tuple{Int8, Int8}
    max_prob::Float32              # Best probability among transfer candidates
    is_target::Bool                # Whether representative is target
    is_transfer_decoy::Bool        # Whether best cross-run match is decoy
    file_idx::Int                  # File containing the representative
    original_indices::Vector{Int}  # Indices of all candidates in this group
end

# FTR control configuration
@enum FTRControlMode begin
    FTR_EXPERIMENT_WIDE = 1
    FTR_EXPERIMENT_LOCAL = 2  
    FTR_GLOBAL_EXPERIMENT = 3
end
```

### New Functions

#### 1. Precursor Representative Selection
```julia
function create_precursor_representatives(
    merged_df::DataFrame,
    candidate_mask::Vector{Bool}
) -> Vector{PrecursorRepresentative}
    # Group candidates by (precursor_idx, isotopes_captured)
    # Select max probability representative per group
    # Create PrecursorRepresentative objects
end
```

#### 2. Global FTR Control  
```julia
function apply_global_ftr_control!(
    representatives::Vector{PrecursorRepresentative},
    α_global::Float64
) -> Set{Tuple{UInt32, Tuple{Int8, Int8}}}
    # Calculate q-values for representatives
    # Return set of passing precursor groups
end
```

#### 3. Global-Experiment FTR Pipeline
```julia
function apply_global_experiment_mbr_filter!(
    merged_df::DataFrame,
    params::ScoringSearchParameters,
    fdr_scale_factor::Float32,
    file_names::Vector{String}
) -> Symbol
    # Stage 1: Global precursor-level control
    # Stage 2: Experiment-wide control on filtered set
    # Return filtered probability column name
end
```

#### 4. Unified FTR Dispatcher
```julia
function apply_mbr_filter_by_mode!(
    merged_df::DataFrame,
    params::ScoringSearchParameters,
    fdr_scale_factor::Float32,
    file_names::Vector{String}
) -> Symbol
    mode = parse_ftr_control_mode(params.mbr_ftr_control)
    
    return if mode == FTR_EXPERIMENT_WIDE
        apply_mbr_filter!(merged_df, params, fdr_scale_factor)  # Original single-tier
    elseif mode == FTR_EXPERIMENT_LOCAL
        apply_hierarchical_mbr_filter!(merged_df, params, fdr_scale_factor, file_names)  # Current two-tier
    elseif mode == FTR_GLOBAL_EXPERIMENT
        apply_global_experiment_mbr_filter!(merged_df, params, fdr_scale_factor, file_names)  # New two-tier
    end
end
```

## Diagnostic Logging Enhancement

### Strategy A: "experiment_wide" Logging
```
=== EXPERIMENT-WIDE FTR CONTROL ===
Total candidates: 144,348
Experiment-wide threshold: 0.6789
Candidates retained: 98,234 (68.0%)
```

### Strategy B: "experiment_local" Logging (Current)
```
=== HIERARCHICAL MBR FILTERING ===
Processing 2,602,286 PSMs across 40 files
MBR transfer candidates: 144,348 / 2,602,286 (5.5%)

STAGE 1: GLOBAL FTR ANALYSIS
  Experiment-wide threshold: 0.6789
  Candidates passing: 120,457 (83.4%)

STAGE 2: PER-FILE FTR REFINEMENT  
  Files with local refinement: 35
  Additional candidates filtered: 3,245
  Final candidates retained: 117,212 (81.2% of original)
```

### Strategy C: "global_experiment" Logging (New)
```
=== GLOBAL-EXPERIMENT FTR CONTROL ===
Processing 144,348 MBR transfer candidates

STAGE 1: GLOBAL PRECURSOR CONTROL
  Precursor groups: 45,123
  Representatives analyzed: 45,123  
  Global q-value threshold: 0.01
  Precursor groups passing: 38,456 (85.2%)
  Candidates from passing groups: 120,457 (83.4%)

STAGE 2: EXPERIMENT-WIDE CONTROL  
  Remaining candidates: 120,457
  Experiment-wide threshold: 0.6789
  Final candidates retained: 98,234 (68.0% of original)
```

## Performance Considerations

### Computational Complexity
- **Tier 1**: O(n_candidates) for grouping + O(n_precursors) for FTR calculation
- **Tier 2**: O(n_remaining * log(n_remaining)) for sorting + O(n_remaining) for FTR
- **Tier 3**: O(n_files * n_candidates_per_file) 

### Memory Usage
- **Additional**: `Vector{PrecursorRepresentative}` (~1KB per precursor group)
- **Temporary**: Precursor grouping dictionaries
- **Overall**: Minimal increase (<5% typical)

### Expected Performance Impact
- **"experiment_wide"**: Baseline (fastest) - original single-tier
- **"experiment_local"**: +10-20% compute time - current two-tier implementation
- **"global_experiment"**: +15-25% compute time - new precursor-level filtering

## Implementation Phase Plan

### Phase 1: Core Infrastructure (Week 1)
1. **Data Structures**: Define `PrecursorRepresentative` and `FTRControlMode`
2. **Configuration**: Add `mbr_ftr_control` parameter parsing
3. **Dispatcher**: Implement `apply_mbr_filter_by_mode!`
4. **Testing**: Unit tests for basic parameter parsing

### Phase 2: Global FTR Control (Week 2)
1. **Representative Selection**: `create_precursor_representatives()`
2. **Global Filtering**: `apply_global_ftr_control!()` 
3. **Global-Only Mode**: Complete `apply_global_only_mbr_filter!()`
4. **Testing**: Validate precursor representative logic

### Phase 3: Three-Tier Integration (Week 3)
1. **Three-Tier Pipeline**: `apply_three_tier_mbr_filter!()`
2. **Integration**: Connect with existing two-tier system
3. **Diagnostics**: Enhanced logging for all three tiers
4. **Testing**: End-to-end validation with real datasets

### Phase 4: Performance & Documentation (Week 4)
1. **Optimization**: Performance tuning and memory optimization
2. **Documentation**: Update user guides and technical documentation
3. **Validation**: Compare FTR control effectiveness across modes
4. **Release**: Integration testing and deployment

## Integration Points

### ScoringSearch.jl Changes
```julia
# Replace current MBR filter call
if params.match_between_runs
    prob_column = apply_mbr_filter_by_mode!(
        merged_df, 
        scoring_params, 
        fdr_scale_factor, 
        file_names
    )
else
    prob_column = :prob
end
```

### Parameter Migration
- **Backward Compatibility**: `use_hierarchical_ftr: true` → `mbr_ftr_control: "experiment_local"`
- **Default**: `mbr_ftr_control: "experiment_local"` (preserve current behavior)
- **Validation**: Error on invalid mode strings

### Default Configuration Update
```json
"global": {
    "scoring": {
        "q_value_threshold": 0.01,
        "mbr_ftr_control": "experiment_local",
        "hierarchical_ftr_min_candidates": 50
    }
}
```

## Expected Benefits

### Scientific Benefits
1. **Improved Specificity**: Global tier eliminates poor precursor groups early
2. **Better Sensitivity**: Hierarchical refinement preserves good candidates
3. **Interpretability**: Clear three-level FTR reporting
4. **Flexibility**: Users can choose appropriate stringency level

### Technical Benefits  
1. **Modularity**: Clean separation of FTR control levels
2. **Performance Options**: Trade-off between speed and stringency
3. **Extensibility**: Framework supports additional FTR strategies
4. **Diagnostics**: Comprehensive logging for method development

## Risk Mitigation

### Potential Issues
1. **Overfiltration**: Three-tier might be too stringent
   - **Mitigation**: Provide "global_only" and "two_tier" fallbacks
2. **Performance Impact**: Additional computation overhead
   - **Mitigation**: Benchmark and optimize, provide "none" option
3. **Complexity**: More configuration options to understand
   - **Mitigation**: Clear documentation, sensible defaults

### Testing Strategy
1. **Unit Tests**: Each tier independently
2. **Integration Tests**: Full three-tier pipeline  
3. **Performance Tests**: Benchmark against current system
4. **Scientific Validation**: Compare FTR effectiveness on real datasets

## Success Metrics

### Technical Metrics
- **Performance**: <40% increase in FTR control time
- **Memory**: <10% increase in peak memory usage
- **Reliability**: Zero regressions in existing two-tier mode

### Scientific Metrics  
- **Specificity**: Improved FTR control effectiveness
- **Sensitivity**: Maintained or improved identification rates
- **Consistency**: Reduced variance in FTR across files/experiments

This plan provides a comprehensive framework for implementing sophisticated, tiered FTR control while maintaining backward compatibility and providing clear performance trade-offs for users.