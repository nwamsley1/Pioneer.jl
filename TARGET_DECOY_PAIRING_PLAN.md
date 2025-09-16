# Target-Decoy Pairing Implementation Plan

## Overview

This document outlines a plan to implement 1:1 precursor-level target-decoy pairing in `sort_of_percolator_in_memory` by adding a `pair_idx` column that ensures balanced target-decoy relationships through strategic decoy cloning.

## Problem Statement

### Current Issues
1. **Imbalanced PSM Representation**: Target precursors may have more PSM instances across runs/isotopes than their paired decoy precursors
2. **Many-to-One Relationships**: Multiple target precursors may be associated with the same decoy precursor
3. **Cross-Run Imbalances**: Target precursors appear in more (run, isotopes_captured) combinations than their paired decoys
4. **Pairing Ambiguity**: No clear 1:1 mapping between target and decoy precursors at the PSM level

### Goals
1. **1:1 Precursor Pairing**: Each target precursor paired with exactly one decoy precursor
2. **Shared pair_idx**: All PSM instances of a target precursor share the same pair_idx with all instances of its paired decoy precursor
3. **Balanced PSM Counts**: Equal number of PSM instances between paired target and decoy precursors
4. **Strategic Decoy Cloning**: Create decoy clones to match target precursor's (run, isotopes_captured) coverage
5. **Clean Removal**: Clone decoys removed after processing, preserving original data integrity

## Technical Analysis

### Data Structure Understanding
```julia
# PSM DataFrame structure (relevant columns):
# - precursor_idx::UInt32        # Unique precursor identifier
# - target::Bool                 # true=target, false=decoy
# - isotopes_captured::Tuple     # Isotope information
# - ms_file_idx::Int            # Run/file identifier
# - [other PSM features...]
```

### Key Relationships
1. **Precursor Level**: Each target precursor should map to exactly one decoy precursor via shared pair_idx
2. **PSM Balancing**: Each target precursor should have equal PSM instance counts as its paired decoy precursor
3. **Library Level**: Target-decoy relationships defined in spectral library via pair information

### Pairing Example
```
Target A precursor_idx=100, appears in:
  - (run=1, isotopes=(0,1)) → pair_idx = 1
  - (run=2, isotopes=(0,1)) → pair_idx = 1  
  - (run=3, isotopes=(1,2)) → pair_idx = 1

Decoy C precursor_idx=200, naturally appears in:
  - (run=1, isotopes=(0,1)) → pair_idx = 1

Target B precursor_idx=101, appears in:
  - (run=1, isotopes=(0,1)) → pair_idx = 2
  - (run=2, isotopes=(1,2)) → pair_idx = 2

Original Decoy C would be paired with Target A, but Target B needs a decoy.
Clone Decoy C (precursor_idx=200_clone) created for:
  - (run=1, isotopes=(0,1)) → pair_idx = 2
  - (run=2, isotopes=(1,2)) → pair_idx = 2

Result:
- pair_idx=1: Target A (3 instances) + Decoy C (1 instance) + 2 cloned instances
- pair_idx=2: Target B (2 instances) + Clone Decoy C (2 instances)
```

## Implementation Strategy

### Phase 1: Precursor-Level Analysis
1. **Build Target-Decoy Mapping**: Extract precursor pairing from spectral library
2. **Count PSM Instances per Precursor**: Group by `precursor_idx` and count (run, isotopes_captured) combinations
3. **Identify Imbalances**: Find target precursors with more instances than their paired decoy precursors
4. **Determine Cloning Requirements**: Calculate missing decoy instances needed for balance

### Phase 2: Strategic Decoy Cloning
```julia
function balance_target_decoy_pairs!(psms::DataFrame, spec_lib)
    # 1. Extract precursor-level pairing information from spectral library
    target_to_decoy_map = build_precursor_pairing_map(spec_lib)
    
    # 2. Count PSM instances per precursor
    precursor_instances = count_precursor_instances(psms)
    
    # 3. Identify cloning requirements for each target-decoy pair
    cloning_plan = plan_decoy_cloning(precursor_instances, target_to_decoy_map)
    
    # 4. Create decoy clones to balance PSM counts
    cloned_decoys = create_balanced_decoy_clones(psms, cloning_plan)
    
    # 5. Add cloned decoys to PSM dataframe
    append!(psms, DataFrame(cloned_decoys))
    
    # 6. Assign pair_idx to all target-decoy precursor pairs
    assign_precursor_pair_indices!(psms, target_to_decoy_map)
    
    return psms
end
```

### Phase 3: Precursor-Level Pair Index Assignment
```julia
function assign_precursor_pair_indices!(psms::DataFrame, target_to_decoy_map::Dict)
    psms[!, :pair_idx] = zeros(UInt32, nrow(psms))
    psms[!, :is_clone] = falses(nrow(psms))  # Track cloned rows
    
    current_pair_idx = UInt32(1)
    
    # Assign same pair_idx to all instances of paired target-decoy precursors
    for (target_precursor, decoy_precursor) in target_to_decoy_map
        # All instances of target precursor get same pair_idx
        target_mask = (psms.precursor_idx .== target_precursor) .& psms.target
        psms.pair_idx[target_mask] .= current_pair_idx
        
        # All instances of paired decoy precursor get same pair_idx
        decoy_mask = (psms.precursor_idx .== decoy_precursor) .& (.!psms.target)
        psms.pair_idx[decoy_mask] .= current_pair_idx
        
        current_pair_idx += 1
    end
end
```

### Phase 4: Clone Removal
```julia
function remove_cloned_decoys!(psms::DataFrame)
    # Remove rows marked as clones
    filter!(:is_clone => (x -> !x), psms)
    # Remove temporary tracking columns
    select!(psms, Not([:is_clone]))
end
```

## Detailed Implementation Plan

### Step 1: Precursor Instance Analysis
```julia
function count_precursor_instances(psms::DataFrame)
    # Count (run, isotopes_captured) combinations per precursor
    instance_counts = Dict{UInt32, Vector{Tuple{Int, Tuple}}}()
    
    for row in eachrow(psms)
        precursor = row.precursor_idx
        instance_key = (row.ms_file_idx, row.isotopes_captured)
        
        if precursor ∉ keys(instance_counts)
            instance_counts[precursor] = []
        end
        
        if instance_key ∉ instance_counts[precursor]
            push!(instance_counts[precursor], instance_key)
        end
    end
    
    return instance_counts
end
```

### Step 2: Cloning Plan Development
```julia
function plan_decoy_cloning(instance_counts::Dict, target_to_decoy_map::Dict)
    cloning_plan = Dict{UInt32, Vector{Tuple{Int, Tuple}}}()  # decoy_precursor => missing_instances
    
    for (target_precursor, decoy_precursor) in target_to_decoy_map
        target_instances = get(instance_counts, target_precursor, [])
        decoy_instances = get(instance_counts, decoy_precursor, [])
        
        # Find instances where target exists but decoy doesn't
        missing_decoy_instances = setdiff(target_instances, decoy_instances)
        
        if length(missing_decoy_instances) > 0
            cloning_plan[decoy_precursor] = missing_decoy_instances
        end
    end
    
    return cloning_plan
end
```

### Step 3: Targeted Decoy Cloning
```julia
function create_balanced_decoy_clones(psms::DataFrame, cloning_plan::Dict)
    cloned_rows = DataFrame()
    
    for (decoy_precursor, missing_instances) in cloning_plan
        # Find best template decoy row for cloning
        template_decoy = find_best_decoy_template(psms, decoy_precursor)
        
        # Create clones for each missing instance
        for (run_idx, isotopes) in missing_instances
            cloned_row = clone_decoy_for_instance(template_decoy, run_idx, isotopes)
            push!(cloned_rows, cloned_row)
        end
    end
    
    # Mark all clones for later removal
    cloned_rows[!, :is_clone] = trues(nrow(cloned_rows))
    
    return cloned_rows
end

function clone_decoy_for_instance(template_row, target_run_idx::Int, target_isotopes::Tuple)
    cloned_row = copy(template_row)
    
    # Update instance-specific fields to match target
    cloned_row.ms_file_idx = target_run_idx
    cloned_row.isotopes_captured = target_isotopes
    
    # Mark as clone
    cloned_row.is_clone = true
    
    # Keep other fields (precursor_idx, target=false, library features, etc.)
    return cloned_row
end
```

## Integration Points

### Function Placement in sort_of_percolator_in_memory
```julia
function sort_of_percolator_in_memory!(...)
    # ... existing preprocessing ...
    
    # NEW: Add precursor-level target-decoy pairing at the beginning
    balance_target_decoy_pairs!(psms, spec_lib)
    
    # ... existing XGBoost training and scoring ...
    # All ML training and scoring operates on balanced target-decoy pairs
    # Each target precursor has equal representation with its paired decoy precursor
    
    # NEW: Remove cloned decoys at the end (before final output)
    remove_cloned_decoys!(psms)
    
    # ... existing postprocessing and output ...
end
```

## Clarifying Questions & Considerations

### Technical Questions

1. **Spectral Library Access**: How do we access the target-decoy pairing information from `spec_lib`? Is there a function like `getPairId(spec_lib, precursor_idx)`?

2. **Cloning Strategy**: When cloning decoys, which features should be identical vs. modified?
   - **Identical**: `precursor_idx`, `target` (false), library-derived features
   - **Modified**: `ms_file_idx`, `isotopes_captured` (to match target)
   - **Uncertain**: Scores, RT values, intensity features

3. **Cross-Run Behavior**: If target A appears in runs 1, 2, 3 but its paired decoy B only appears in run 1, do we:
   - Clone decoy B into runs 2 and 3?
   - Or only pair within runs where both naturally occur?

4. **Memory Impact**: Cloning could significantly increase DataFrame size. Should we implement memory monitoring?

5. **Feature Consistency**: How do we ensure cloned decoys have realistic feature values that don't break downstream ML training?

### Algorithmic Questions

6. **Pairing Priority**: If multiple decoys could pair with a target, what's the selection criteria?
   - Closest in m/z?
   - Best scoring decoy?
   - Library-defined pairs only?

7. **Isotope Handling**: How do `isotopes_captured` values affect pairing?
   - Should clones inherit target's isotopes_captured?
   - Or maintain decoy's original isotopes_captured?

8. **Run Assignment**: When cloning decoys for cross-run pairing:
   - Use target's `ms_file_idx`?
   - Use original decoy's `ms_file_idx`?
   - Create hybrid assignments?

### Validation Questions

9. **Quality Control**: How do we validate that:
   - All targets have exactly one paired decoy?
   - No orphaned decoys exist after pairing?
   - Cloning doesn't bias ML training?

10. **Performance Impact**: What's the expected computational overhead of this approach?

## Potential Risks & Mitigation

### Risk 1: Memory Explosion
**Problem**: Extensive cloning could double or triple DataFrame size
**Mitigation**: 
- Implement memory usage monitoring
- Use streaming processing for very large datasets
- Consider lazy cloning (generate clones on-demand)

### Risk 2: Feature Authenticity
**Problem**: Cloned decoys may have unrealistic feature combinations
**Mitigation**:
- Clone only essential columns, compute others fresh
- Add feature validation for cloned rows
- Consider feature interpolation for realistic values

### Risk 3: ML Bias Introduction
**Problem**: Artificial clones could bias XGBoost training
**Mitigation**:
- Exclude clones from training set
- Weight cloned samples differently
- Validate model performance with/without cloning

### Risk 4: Cross-Run Inconsistencies
**Problem**: Cloning decoys into runs where they weren't observed
**Mitigation**:
- Only clone within runs where original decoy has some evidence
- Use RT prediction models to ensure realistic RT values for clones
- Implement cross-run feature consistency checks

## Success Metrics

1. **Pairing Completeness**: 100% of targets have corresponding decoys
2. **Balance Achievement**: Equal target and decoy counts post-cloning
3. **Memory Efficiency**: <50% memory increase from cloning
4. **Performance Preservation**: <10% increase in processing time
5. **Quality Maintenance**: No degradation in PSM identification rates

## Implementation Timeline

### Week 1: Infrastructure
- Implement spectral library pairing extraction
- Build target-decoy counting functions
- Create basic cloning framework

### Week 2: Core Algorithm
- Implement balanced cloning algorithm
- Add pair_idx assignment logic
- Build clone removal functions

### Week 3: Integration & Testing
- Integrate into sort_of_percolator_in_memory
- Unit tests for each component
- End-to-end testing with real data

### Week 4: Validation & Optimization
- Performance benchmarking
- Memory usage optimization
- ML bias analysis and mitigation

This plan provides a comprehensive approach to implementing 1:1 target-decoy pairing while addressing key technical challenges and potential risks. The modular design allows for iterative implementation and testing of each component.