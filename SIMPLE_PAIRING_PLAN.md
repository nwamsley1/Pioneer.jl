# Simple Target-Decoy Pairing Implementation Plan

## Overview

This document outlines a simplified approach to implement minority-to-majority precursor pairing in `sort_of_percolator_in_memory` through random assignment without any row cloning.

## Problem Statement

### Current Issues
1. **No Pairing Structure**: No systematic 1:1 mapping between target and decoy precursors
2. **Missing pair_idx**: No way to identify which targets and decoys should be compared together  
3. **Unbalanced Comparisons**: ML training may compare targets and decoys that aren't meaningfully paired
4. **Inefficient Use of Data**: No strategy to maximize pairing when target/decoy counts are imbalanced

### Goals
1. **Minority-to-Majority Pairing**: Each precursor from the smaller set (minority) randomly paired with one precursor from the larger set (majority) within iRT bins
   - If more decoys than targets: each target gets paired with one random decoy from same iRT bin
   - If more targets than decoys: each decoy gets paired with one random target from same iRT bin
2. **iRT Stratified Assignment**: Pair within iRT bins containing ~1000 targets each, with overflow to nearby bins
3. **Shared pair_idx**: All PSM instances belonging to paired target-decoy precursors share the same pair_idx
4. **Reproducible Pairing**: Use fixed random seed for consistent results across runs
5. **Type Safety**: pair_idx as `Union{Missing, UInt32}` with proper missing value handling
6. **No Cloning**: No artificial rows created - only use existing PSM data

## Technical Analysis

### Data Structure
```julia
# PSM DataFrame structure (relevant columns):
# - precursor_idx::UInt32        # Unique precursor identifier  
# - target::Bool                 # true=target, false=decoy
# - [other PSM features...]
# 
# NEW COLUMN TO ADD:
# - pair_idx::Union{Missing, UInt32}  # Shared identifier for paired target-decoy precursors (missing for unpaired)
```

### Pairing Example
```
### Example: iRT-Stratified Pairing
Available precursors in PSMs with iRT values:
- Target precursors: [(100, 10.5), (101, 15.2), (102, 22.1), (103, 28.7), (104, 35.3)]  # 5 targets
- Decoy precursors: [(200, 12.1), (201, 16.8), (202, 25.4)]                             # 3 decoys

iRT Bins Created (target: ~1000 per bin, but only 5 total so use 2 bins):
- Bin 1: iRT 10.0 - 20.0  
- Bin 2: iRT 20.0 - 40.0

Bin Assignments:
- Bin 1 targets: [100, 101] (iRT: 10.5, 15.2)
- Bin 1 decoys: [200, 201] (iRT: 12.1, 16.8)  
- Bin 2 targets: [102, 103, 104] (iRT: 22.1, 28.7, 35.3)
- Bin 2 decoys: [202] (iRT: 25.4)

Strategy: More targets than decoys globally, so each decoy gets paired within bin
- Bin 1: Decoy 200 ↔ Target 101, Decoy 201 ↔ Target 100 (pair_idx = 1, 2)
- Bin 2: Decoy 202 ↔ Target 103 (pair_idx = 3)
- Overflow: Target 102 and 104 remain unpaired (pair_idx = missing)

Final Result:
- pair_idx=1: Target 101 ↔ Decoy 200  
- pair_idx=2: Target 100 ↔ Decoy 201
- pair_idx=3: Target 103 ↔ Decoy 202
- Unpaired: Target 102, Target 104 (pair_idx = missing)
```

## Implementation Strategy

### Core Algorithm
```julia
const PAIRING_RANDOM_SEED = 1844  # Fixed seed for reproducible pairing
const TARGET_PRECURSORS_PER_BIN = 1000  # Target number of precursors per iRT bin

function assign_random_target_decoy_pairs!(psms::DataFrame)
    @info "Starting iRT-stratified target-decoy pairing..."
    
    # Set random seed for reproducibility
    Random.seed!(PAIRING_RANDOM_SEED)
    
    # 1. Create iRT bins
    irt_bins = create_irt_bins(psms)
    @info "Created $(length(irt_bins)) iRT bins"
    
    # 2. Extract unique precursors with their iRT values
    target_precursors = get_unique_precursors_with_irt(psms, true)
    decoy_precursors = get_unique_precursors_with_irt(psms, false)
    
    @info "Found $(length(target_precursors)) unique targets, $(length(decoy_precursors)) unique decoys"
    
    # 3. Create stratified random pairings within iRT bins
    pairings = create_stratified_pairings(target_precursors, decoy_precursors, irt_bins)
    
    @info "Created $(length(pairings)) target-decoy pairs across iRT bins"
    
    # 4. Add pair_idx column and assign values
    assign_pair_indices!(psms, pairings)
    
    # 5. Report pairing statistics
    report_pairing_statistics(psms)
    
    return psms
end
```

### Step 1: Extract Unique Precursors
```julia
function get_unique_precursors_with_irt(psms::DataFrame, is_target::Bool)
    mask = psms.target .== is_target
    precursor_irt_pairs = unique([(row.precursor_idx, row.irt) for row in eachrow(psms[mask, [:precursor_idx, :irt]])])
    return precursor_irt_pairs
end

function create_irt_bins(psms::DataFrame)
    # Get all target precursors with their iRT values
    target_precursors = get_unique_precursors_with_irt(psms, true)
    n_targets = length(target_precursors)
    
    if n_targets < TARGET_PRECURSORS_PER_BIN
        @info "Only $n_targets targets available, using single iRT bin"
        return [(-Inf, Inf)]  # Single bin containing all
    end
    
    # Calculate number of bins needed
    n_bins = max(1, div(n_targets, TARGET_PRECURSORS_PER_BIN))
    
    # Sort targets by iRT and determine bin boundaries
    sorted_irts = sort([irt for (_, irt) in target_precursors])
    
    # Create bin boundaries at quantiles
    bin_edges = Vector{Float64}(undef, n_bins + 1)
    bin_edges[1] = -Inf
    bin_edges[end] = Inf
    
    for i in 2:n_bins
        quantile_pos = (i - 1) / n_bins
        idx = max(1, min(length(sorted_irts), round(Int, quantile_pos * length(sorted_irts))))
        bin_edges[i] = sorted_irts[idx]
    end
    
    # Convert to (min, max) tuples
    bins = [(bin_edges[i], bin_edges[i+1]) for i in 1:n_bins]
    
    @info "Created $n_bins iRT bins with ~$(div(n_targets, n_bins)) targets per bin"
    return bins
end
```

### Step 2: Create Random Pairings
```julia
function create_stratified_pairings(target_precursors::Vector{Tuple{UInt32, Float64}}, 
                                   decoy_precursors::Vector{Tuple{UInt32, Float64}}, 
                                   irt_bins::Vector{Tuple{Float64, Float64}})
    n_targets = length(target_precursors)
    n_decoys = length(decoy_precursors)
    
    if n_targets == 0 || n_decoys == 0
        @warn "No target-decoy pairs possible (one set is empty)"
        return Tuple{UInt32, UInt32}[]
    end
    
    # Pre-allocate result array
    all_pairings = Vector{Tuple{UInt32, UInt32}}()
    sizehint!(all_pairings, min(n_targets, n_decoys))
    
    # Assign precursors to bins
    target_bins = assign_precursors_to_bins(target_precursors, irt_bins)
    decoy_bins = assign_precursors_to_bins(decoy_precursors, irt_bins)
    
    # Determine global pairing strategy
    minority_is_target = n_targets <= n_decoys
    @info "Pairing strategy: Each $(minority_is_target ? "target" : "decoy") paired with random $(minority_is_target ? "decoy" : "target") within iRT bins"
    
    # Track unpaired precursors per bin for smart overflow handling
    unpaired_targets_by_bin = Vector{Vector{UInt32}}(undef, length(irt_bins))
    unpaired_decoys_by_bin = Vector{Vector{UInt32}}(undef, length(irt_bins))
    
    # Diagnostic counters
    bins_with_excess_targets = 0
    bins_with_excess_decoys = 0
    
    # Pair within each bin
    for (bin_idx, (min_irt, max_irt)) in enumerate(irt_bins)
        bin_targets = get(target_bins, bin_idx, UInt32[])
        bin_decoys = get(decoy_bins, bin_idx, UInt32[])
        
        @info "Bin $bin_idx analysis: $(length(bin_targets)) targets, $(length(bin_decoys)) decoys (iRT: $(round(min_irt, digits=2)) - $(round(max_irt, digits=2)))"
        
        # Diagnose bin imbalance
        if length(bin_targets) > length(bin_decoys)
            bins_with_excess_targets += 1
            @info "  → Excess targets in bin $bin_idx: $(length(bin_targets) - length(bin_decoys)) targets will remain unpaired"
        elseif length(bin_decoys) > length(bin_targets)
            bins_with_excess_decoys += 1
            @info "  → Excess decoys in bin $bin_idx: $(length(bin_decoys) - length(bin_targets)) decoys need pairing with targets from other bins"
        end
        
        bin_pairings, unpaired_t, unpaired_d = pair_within_bin(bin_targets, bin_decoys)
        append!(all_pairings, bin_pairings)
        
        unpaired_targets_by_bin[bin_idx] = unpaired_t
        unpaired_decoys_by_bin[bin_idx] = unpaired_d
        
        @info "  → Bin $bin_idx result: $(length(bin_pairings)) pairs, $(length(unpaired_t)) unpaired targets, $(length(unpaired_d)) unpaired decoys"
    end
    
    @info "Bin summary: $bins_with_excess_targets bins with excess targets, $bins_with_excess_decoys bins with excess decoys"
    
    # Handle overflow pairing: prioritize pairing excess decoys with nearby unpaired targets
    overflow_pairings = handle_cross_bin_overflow(unpaired_targets_by_bin, unpaired_decoys_by_bin, irt_bins)
    append!(all_pairings, overflow_pairings)
    
    if length(overflow_pairings) > 0
        @info "Overflow pairing across bins: $(length(overflow_pairings)) additional pairs"
    end
    
    total_unpaired = (minority_is_target ? length(unpaired_decoys) : length(unpaired_targets)) - length(overflow_pairings)
    @info "Final result: $(length(all_pairings)) total pairs, $total_unpaired unpaired precursors"
    
    return all_pairings
end

function assign_precursors_to_bins(precursors::Vector{Tuple{UInt32, Float64}}, 
                                 irt_bins::Vector{Tuple{Float64, Float64}})
    bin_assignments = Dict{Int, Vector{UInt32}}()
    
    # Pre-allocate bin vectors
    for i in 1:length(irt_bins)
        bin_assignments[i] = Vector{UInt32}()
    end
    
    for (precursor_idx, irt_val) in precursors
        # Find which bin this precursor belongs to
        bin_idx = findfirst(((min_irt, max_irt),) -> min_irt <= irt_val < max_irt, irt_bins)
        
        if bin_idx === nothing
            # Handle edge case - assign to last bin if >= max value
            bin_idx = length(irt_bins)
        end
        
        push!(bin_assignments[bin_idx], precursor_idx)
    end
    
    return bin_assignments
end

function pair_within_bin(targets::Vector{UInt32}, decoys::Vector{UInt32})
    n_targets = length(targets)
    n_decoys = length(decoys)
    
    if n_targets == 0 || n_decoys == 0
        return Tuple{UInt32, UInt32}[], targets, decoys
    end
    
    # Always try to pair as many as possible within bin
    n_pairs = min(n_targets, n_decoys)
    
    # Use randperm for efficient random selection without copying full arrays
    target_indices = randperm(n_targets)
    decoy_indices = randperm(n_decoys)
    
    # Create pairings
    pairings = [(targets[target_indices[i]], decoys[decoy_indices[i]]) for i in 1:n_pairs]
    
    # Determine unpaired precursors
    unpaired_targets = targets[target_indices[(n_pairs+1):end]]
    unpaired_decoys = decoys[decoy_indices[(n_pairs+1):end]]
    
    return pairings, unpaired_targets, unpaired_decoys
end

function handle_cross_bin_overflow(unpaired_targets_by_bin::Vector{Vector{UInt32}}, 
                                  unpaired_decoys_by_bin::Vector{Vector{UInt32}},
                                  irt_bins::Vector{Tuple{Float64, Float64}})
    cross_bin_pairings = Vector{Tuple{UInt32, UInt32}}()
    n_bins = length(irt_bins)
    
    # Priority: pair excess decoys with unpaired targets in nearby bins
    for bin_idx in 1:n_bins
        excess_decoys = unpaired_decoys_by_bin[bin_idx]
        
        if isempty(excess_decoys)
            continue
        end
        
        @info "Processing $(length(excess_decoys)) excess decoys from bin $bin_idx"
        
        # Find unpaired targets in nearby bins (search outward from current bin)
        for distance in 0:(n_bins-1)
            if isempty(excess_decoys)
                break
            end
            
            # Check bins at this distance
            nearby_bins = Int[]
            if distance == 0
                push!(nearby_bins, bin_idx)
            else
                # Add bins at +/- distance
                if bin_idx - distance >= 1
                    push!(nearby_bins, bin_idx - distance)
                end
                if bin_idx + distance <= n_bins
                    push!(nearby_bins, bin_idx + distance)
                end
            end
            
            for target_bin_idx in nearby_bins
                if isempty(excess_decoys)
                    break
                end
                
                available_targets = unpaired_targets_by_bin[target_bin_idx]
                if isempty(available_targets)
                    continue
                end
                
                # Pair as many as possible
                n_pairs = min(length(excess_decoys), length(available_targets))
                
                if n_pairs > 0
                    @info "  → Cross-bin pairing: $(n_pairs) decoys from bin $bin_idx with targets from bin $target_bin_idx (distance: $distance)"
                    
                    # Use randperm for random selection
                    decoy_indices = randperm(length(excess_decoys))[1:n_pairs]
                    target_indices = randperm(length(available_targets))[1:n_pairs]
                    
                    for i in 1:n_pairs
                        push!(cross_bin_pairings, (available_targets[target_indices[i]], excess_decoys[decoy_indices[i]]))
                    end
                    
                    # Remove paired precursors
                    excess_decoys = excess_decoys[setdiff(1:length(excess_decoys), decoy_indices)]
                    unpaired_targets_by_bin[target_bin_idx] = available_targets[setdiff(1:length(available_targets), target_indices)]
                end
            end
        end
        
        if !isempty(excess_decoys)
            @info "  → Warning: $(length(excess_decoys)) decoys from bin $bin_idx could not be paired (no available targets)"
        end
    end
    
    return cross_bin_pairings
end
```

### Step 3: Assign Pair Indices
```julia
function assign_pair_indices!(psms::DataFrame, pairings::Vector{Tuple{UInt32, UInt32}})
    # Initialize pair_idx column with missing values (proper type handling)
    psms[!, :pair_idx] = Vector{Union{Missing, UInt32}}(undef, nrow(psms))
    fill!(psms.pair_idx, missing)
    
    # Assign pair_idx for each target-decoy pair
    for (pair_id, (target_precursor, decoy_precursor)) in enumerate(pairings)
        # All instances of this target precursor get this pair_idx
        target_mask = (psms.precursor_idx .== target_precursor) .& (psms.target .== true)
        psms.pair_idx[target_mask] .= UInt32(pair_id)
        
        # All instances of this decoy precursor get this pair_idx  
        decoy_mask = (psms.precursor_idx .== decoy_precursor) .& (psms.target .== false)
        psms.pair_idx[decoy_mask] .= UInt32(pair_id)
    end
end
```

### Step 4: Validation and Reporting
```julia
function report_pairing_statistics(psms::DataFrame)
    # Count paired vs unpaired rows
    paired_mask = .!ismissing.(psms.pair_idx)
    n_paired = sum(paired_mask)
    n_unpaired = sum(.!paired_mask)
    
    @info "Pairing complete: $n_paired PSMs paired, $n_unpaired PSMs unpaired"
    
    # Count by target/decoy
    paired_targets = sum(paired_mask .& psms.target)
    paired_decoys = sum(paired_mask .& .!psms.target)
    unpaired_targets = sum(.!paired_mask .& psms.target)
    unpaired_decoys = sum(.!paired_mask .& .!psms.target)
    
    @info "Targets: $paired_targets paired, $unpaired_targets unpaired"
    @info "Decoys: $paired_decoys paired, $unpaired_decoys unpaired"
    
    # Validate pairing balance
    validate_pairing_balance(psms)
end

function validate_pairing_balance(psms::DataFrame)
    paired_psms = psms[.!ismissing.(psms.pair_idx), :]
    
    if nrow(paired_psms) == 0
        @warn "No paired PSMs found"
        return
    end
    
    # Count targets and decoys for each pair_idx
    pair_counts = combine(groupby(paired_psms, :pair_idx)) do group
        (
            n_targets = sum(group.target),
            n_decoys = sum(.!group.target),
            target_precursors = length(unique(group.precursor_idx[group.target])),
            decoy_precursors = length(unique(group.precursor_idx[.!group.target]))
        )
    end
    
    @info "Validation: $(nrow(pair_counts)) pairs created"
    
    # Check for 1:1 precursor pairing
    invalid_pairs = pair_counts[(pair_counts.target_precursors .!= 1) .| (pair_counts.decoy_precursors .!= 1), :]
    if nrow(invalid_pairs) > 0
        @error "Invalid pairing detected: some pairs don't have exactly 1 target and 1 decoy precursor"
        println(invalid_pairs)
    else
        @info "✓ All pairs have exactly 1 target and 1 decoy precursor"
    end
end
```

## Integration Point

```julia
function sort_of_percolator_in_memory!(...)
    # ... existing preprocessing ...
    
    # NEW: Add random target-decoy pairing
    assign_random_target_decoy_pairs!(psms)
    
    # ... existing XGBoost training and scoring ...
    # NOTE: Training/scoring may need to handle missing pair_idx values
    
    # ... existing postprocessing and output ...
    # NOTE: Consider filtering out unpaired PSMs before final output
end
```

## Diagnostic Output & Verification

The implementation includes comprehensive diagnostic prints to verify correct operation:

### Bin Creation Diagnostics
```
@info "Created 3 iRT bins with ~1200 targets per bin"
```

### Per-Bin Analysis
```
@info "Bin 1 analysis: 1180 targets, 890 decoys (iRT: 10.50 - 25.30)"
@info "  → Excess targets in bin 1: 290 targets will remain unpaired"
@info "  → Bin 1 result: 890 pairs, 290 unpaired targets, 0 unpaired decoys"
```

### Cross-Bin Overflow Handling
```
@info "Bin summary: 2 bins with excess targets, 1 bins with excess decoys"
@info "Processing 45 excess decoys from bin 3"
@info "  → Cross-bin pairing: 30 decoys from bin 3 with targets from bin 2 (distance: 1)"
@info "  → Cross-bin pairing: 15 decoys from bin 3 with targets from bin 1 (distance: 2)"
```

### Final Summary
```
@info "Final result: 2340 total pairs, 125 unpaired precursors"
@info "Pairing complete: 4680 PSMs paired, 892 PSMs unpaired"
@info "Targets: 2340 paired, 245 unpaired"
@info "Decoys: 2340 paired, 0 unpaired"
```

### Validation Diagnostics
```
@info "Validation: 2340 pairs created"
@info "✓ All pairs have exactly 1 target and 1 decoy precursor"
```

These diagnostics allow verification that:
1. iRT bins are created with appropriate target counts
2. Bin imbalances are identified and handled correctly
3. Cross-bin overflow pairing works for excess decoys
4. Final pair counts match expectations
5. All paired PSMs have valid target-decoy relationships

## Implementation Questions & Considerations

### Technical Questions

1. **iRT Bin Sizing**: With ~1000 targets per bin, what happens with very small datasets?
   - If <1000 total targets, use single bin containing all precursors
   - Minimum bin size to ensure meaningful pairing within bins

2. **Overflow Handling**: When targets in a bin exceed available decoys (or vice versa):
   - Pair remaining precursors across bins randomly
   - Maintain iRT proximity where possible for overflow pairing

3. **Type Safety**: Using `Union{Missing, UInt32}` for pair_idx column:
   - Proper initialization with `Vector{Union{Missing, UInt32}}(undef, nrow)`
   - Fill with `missing` values initially
   - Downstream functions must handle missing values correctly

4. **Memory Efficiency**: Pre-allocation strategies for O(n) performance:
   - Pre-allocate bin assignment dictionaries
   - Use `sizehint!` for result vectors
   - Avoid quadratic operations in bin assignment

5. **Reproducibility**: Fixed random seed (1844) ensures:
   - Consistent pairing across multiple runs
   - Deterministic results for testing and validation
   - Same target-decoy pairs for identical input data

### Performance Considerations

6. **Random Selection Performance**: Using `randperm()` for efficient random selection
   - `randperm(n)[1:k]` more efficient than `shuffle(array)[1:k]` for partial selection
   - Avoids copying full arrays when only selecting subset

7. **Masking Operations**: Multiple boolean masking operations could be optimized
   - Could pre-compute precursor → row indices mapping

8. **Memory Allocation**: Vector creation and assignments
   - Pre-allocate pair_idx column
   - Use views where possible

### Validation Questions

9. **Quality Control**: How to validate correct pairing?
   - Check each pair_idx has exactly 1 target + 1 decoy precursor
   - Verify no duplicate pairings
   - Confirm unpaired precursors have missing pair_idx

10. **Downstream Impact**: How does this affect existing analysis?
    - Will filtering by pair_idx break existing workflows?
    - Need to update any functions that assume all PSMs are included?

## Potential Issues & Mitigation

### Issue 1: Unpaired Precursors
**Problem**: Some precursors won't be paired and may be excluded from analysis
**Mitigation**: 
- Clear documentation of which precursors are unpaired
- Option to include unpaired PSMs in separate analysis
- Consider strategies to maximize pairing (e.g., abundance-based)

### Issue 2: Random Variation
**Problem**: Different runs produce different pairings, affecting reproducibility
**Mitigation**:
- Optional random seed parameter
- Save pairing results for record-keeping
- Consider deterministic pairing strategies

### Issue 3: Imbalanced Datasets
**Problem**: Highly imbalanced target/decoy ratios leave many precursors unpaired
**Mitigation**:
- Report statistics on pairing efficiency  
- Current strategy maximizes usage of minority set (all get paired)
- Consider abundance-weighted selection from majority set

### Issue 4: Downstream Compatibility
**Problem**: Existing code may not handle missing pair_idx values
**Mitigation**:
- Thorough testing of downstream functions
- Clear documentation of pair_idx usage
- Optional filtering to remove unpaired PSMs

## Success Metrics

1. **Pairing Completeness**: 100% of minority set precursors get paired
2. **Performance**: <10% increase in processing time
3. **Memory Efficiency**: Minimal memory overhead (just one new column)
4. **Validation**: 100% of pairs have exactly 1 target + 1 decoy precursor
5. **Reproducibility**: Consistent results with same random seed

## Implementation Timeline

### Phase 1: Core Implementation (1-2 days)
- Implement precursor extraction functions
- Create random pairing algorithm
- Add pair_idx assignment logic

### Phase 2: Integration & Validation (1 day)
- Integrate into `sort_of_percolator_in_memory`
- Add comprehensive validation and reporting
- Unit tests for each component

### Phase 3: Testing & Optimization (1-2 days)
- Test with real datasets of various sizes
- Performance benchmarking and optimization
- Handle edge cases (empty target/decoy sets, etc.)

This simplified approach eliminates the complexity of row cloning while achieving the core goal of establishing 1:1 target-decoy relationships for downstream analysis.