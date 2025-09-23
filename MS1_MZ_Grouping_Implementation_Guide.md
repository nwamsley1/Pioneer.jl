# MS1 M/Z Grouping Implementation Guide

## Overview

This document provides detailed implementation instructions for combining precursors with identical m/z values into shared design matrix columns during MS1 fitting in Pioneer.jl. This approach solves multicollinearity issues caused by target/decoy pairs and other precursors sharing identical masses.

## Problem Statement

**Current Issue**: In MS1 deconvolution, target and decoy precursors with identical masses create separate columns in the design matrix, leading to multicollinearity that can cause numerical instability and poor coefficient estimation.

**Solution**: Group precursors by m/z value into shared design matrix columns, fit coefficients for the groups, then distribute the fitted coefficients back to individual precursors.

## Implementation Locations

**Primary Files to Modify**:
1. `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl` - New MS1-specific design matrix construction
2. `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl` - Integration in MS1 chromatogram building
3. `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl` - Integration in MS1 PSM identification

**Configuration**: No new parameters required initially - use hardcoded precision.

## Step-by-Step Implementation

### Step 1: Create M/Z Grouping Data Structure

**File**: `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl`

Add the following data structure **at the top of the file**:

```julia
"""
    MzGroupingMap

Data structure for managing m/z-based grouping of precursors in MS1 design matrix construction.
Groups precursors with identical (rounded) m/z values into shared design matrix columns.
"""
mutable struct MzGroupingMap
    # Maps rounded m/z → design matrix column
    mz_to_col::Dict{UInt32, UInt16}

    # Maps precursor ID → rounded m/z (for coefficient retrieval)
    precid_to_mz_group::Dict{UInt32, UInt32}

    # Maps rounded m/z → list of precursor IDs in that group
    mz_group_to_precids::Dict{UInt32, Vector{UInt32}}

    # Precision for m/z rounding (default: 100000 for 5 decimal places)
    mz_precision::UInt32

    # Current column counter
    current_col::UInt16
end

"""
    MzGroupingMap(mz_precision::UInt32 = 100000)

Constructor for MzGroupingMap with configurable m/z precision.
Default precision of 100000 groups masses within ~0.01 ppm at 1000 Da.
"""
function MzGroupingMap(mz_precision::UInt32 = UInt32(100000))
    return MzGroupingMap(
        Dict{UInt32, UInt16}(),
        Dict{UInt32, UInt32}(),
        Dict{UInt32, Vector{UInt32}}(),
        mz_precision,
        UInt16(0)
    )
end

"""
    round_mz_for_grouping(mz::Float32, precision::UInt32) -> UInt32

Round m/z to specified precision for grouping identical masses.
Example: round_mz_for_grouping(1000.12345f0, 100000) = 100012345
"""
function round_mz_for_grouping(mz::Float32, precision::UInt32)::UInt32
    return round(UInt32, mz * precision)
end

"""
    reset!(grouping::MzGroupingMap)

Reset the grouping map for reuse across scans.
"""
function reset!(grouping::MzGroupingMap)
    empty!(grouping.mz_to_col)
    empty!(grouping.precid_to_mz_group)
    empty!(grouping.mz_group_to_precids)
    grouping.current_col = UInt16(0)
    return nothing
end
```

### Step 2: Create MS1-Specific Design Matrix Construction

**File**: `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl`

Add the new function **after the existing `buildDesignMatrix!` function**:

```julia
"""
    buildDesignMatrixMS1!(H::SparseArray{UInt32,Float32},
                          matches::Vector{m}, misses::Vector{m},
                          nmatches::Int64, nmisses::Int64,
                          mz_grouping::MzGroupingMap,
                          precursors::LibraryPrecursors;
                          block_size = 10000) where {m<:MatchIon{Float32}}

MS1-specific design matrix construction that groups precursors by m/z values.
Identical to buildDesignMatrix! except precursors with identical m/z are
assigned to the same design matrix column to avoid multicollinearity.
"""
function buildDesignMatrixMS1!(
    H::SparseArray{UInt32,Float32},
    matches::Vector{m},
    misses::Vector{m},
    nmatches::Int64,
    nmisses::Int64,
    mz_grouping::MzGroupingMap,
    precursors::LibraryPrecursors;
    block_size = 10000
) where {m<:MatchIon{Float32}}

    T = Float32

    # Get precursor m/z values for grouping
    precursor_mzs = getMz(precursors)

    # Calculate number of rows (same logic as original buildDesignMatrix!)
    M = 1
    for i in range(2, nmatches)
        if getPeakInd(matches[i]) != getPeakInd(matches[i - 1])
            M += 1
        end
    end
    M += nmisses

    # Handle array resizing (same logic as original)
    if (nmatches + nmisses) >= length(H.colval) - 1
        block_size = max(block_size, nmatches + nmisses - length(H.colval))
        append!(H.colval, zeros(eltype(H.colval), block_size))
        append!(H.rowval, zeros(eltype(H.rowval), block_size))
        append!(H.nzval, zeros(eltype(H.nzval), block_size))
        append!(H.x, zeros(eltype(H.x), block_size))
        append!(H.matched, zeros(eltype(H.matched), block_size))
        append!(H.isotope, zeros(eltype(H.isotope), block_size))
        append!(H.colptr, Vector{UInt32}(undef, block_size))
    end

    # Initialize tracking variables
    row = 0
    last_peak_ind = -1
    last_row, last_col = -1, -1
    j = 0
    H.n_vals = 0

    # Process matches with m/z grouping
    for i in range(1, nmatches)
        match = matches[i]
        prec_id = getPrecID(match)
        prec_mz = precursor_mzs[prec_id]

        # Get or create m/z group column
        mz_group = round_mz_for_grouping(prec_mz, mz_grouping.mz_precision)

        if !haskey(mz_grouping.mz_to_col, mz_group)
            # Create new column for this m/z group
            mz_grouping.current_col += 1
            mz_grouping.mz_to_col[mz_group] = mz_grouping.current_col
            mz_grouping.mz_group_to_precids[mz_group] = UInt32[]
        end

        # Track precursor → m/z group mapping
        mz_grouping.precid_to_mz_group[prec_id] = mz_group

        # Add precursor to group if not already present
        if prec_id ∉ mz_grouping.mz_group_to_precids[mz_group]
            push!(mz_grouping.mz_group_to_precids[mz_group], prec_id)
        end

        # Use m/z group column instead of individual precursor column
        col = mz_grouping.mz_to_col[mz_group]

        # Handle peak indexing (same as original)
        if getPeakInd(match) != last_peak_ind
            row += 1
            last_peak_ind = getPeakInd(match)
        end

        # Add entry to design matrix (same as original)
        if (col != last_col) | (row != last_row)
            j += 1
            H.n_vals += 1
        end
        H.colval[j] = col
        H.rowval[j] = row
        H.nzval[j] += getPredictedIntensity(match)
        H.x[j] = getIntensity(match)
        H.matched[j] = true
        H.isotope[j] = getIsoIdx(match)

        last_col = col
        last_row = row
    end

    # Process misses with m/z grouping
    i = j + 1
    for j_miss in range(1, nmisses)
        miss = misses[j_miss]
        prec_id = getPrecID(miss)

        # Only add miss if precursor was seen in matches (same logic as original)
        if haskey(mz_grouping.precid_to_mz_group, prec_id)
            mz_group = mz_grouping.precid_to_mz_group[prec_id]
            col = mz_grouping.mz_to_col[mz_group]

            row += 1
            H.colval[i] = col
            H.rowval[i] = row
            H.nzval[i] += getPredictedIntensity(miss)
            H.x[i] = zero(Float32)
            H.matched[i] = false
            H.isotope[i] = getIsoIdx(miss)
            i += 1
            H.n_vals += 1
        end
    end

    # Sort sparse matrix (same as original)
    sortSparse!(H)
end
```

### Step 3: Create Coefficient Distribution Function

**File**: `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl`

Add the coefficient distribution function **after the MS1 design matrix function**:

```julia
"""
    distribute_ms1_coefficients!(individual_weights::Vector{Float32},
                                fitted_weights::Vector{Float32},
                                mz_grouping::MzGroupingMap,
                                precursor_ids::Vector{UInt32})

Distribute fitted coefficients from m/z groups back to individual precursors.
After solving the grouped design matrix, this assigns the fitted coefficient
for each m/z group to all precursors in that group.
"""
function distribute_ms1_coefficients!(
    individual_weights::Vector{Float32},
    fitted_weights::Vector{Float32},
    mz_grouping::MzGroupingMap,
    precursor_ids::Vector{UInt32}
)
    for prec_id in precursor_ids
        if haskey(mz_grouping.precid_to_mz_group, prec_id)
            mz_group = mz_grouping.precid_to_mz_group[prec_id]
            if haskey(mz_grouping.mz_to_col, mz_group)
                group_col = mz_grouping.mz_to_col[mz_group]
                # Assign fitted coefficient to this precursor
                if group_col <= length(fitted_weights)
                    individual_weights[prec_id] = fitted_weights[group_col]
                else
                    individual_weights[prec_id] = 0.0f0
                end
            else
                individual_weights[prec_id] = 0.0f0
            end
        else
            # Precursor not in any group (shouldn't happen for precursors in current scan)
            individual_weights[prec_id] = 0.0f0
        end
    end
end
```

### Step 4: Integrate into IntegrateChromatogramsSearch

**File**: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`

**Location**: In the `build_chromatograms` function with `::MS1CHROM` parameter (around line 425-435)

**Changes needed**:

1. **Add m/z grouping initialization** after line 451 (after `isotopes_dict = getIsotopes(...)`):

```julia
    isotopes_dict = getIsotopes(seqs, pmz, pids, pcharge, QRoots(5), 5)

    # NEW: Create m/z grouping map for MS1
    mz_grouping = MzGroupingMap(UInt32(100000))  # 5 decimal place precision
```

2. **Replace the buildDesignMatrix! call** around line 99-106:

**Find this block**:
```julia
            buildDesignMatrix!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                getIdToCol(search_data)
            )
```

**Replace with**:
```julia
            # Reset grouping for this scan
            reset!(mz_grouping)

            # Use MS1-specific design matrix construction with m/z grouping
            buildDesignMatrixMS1!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                mz_grouping,
                precursors
            )
```

3. **Add coefficient distribution** after the `solveHuber!` call (around line 137):

**Find this block**:
```julia
            solveHuber!(
                Hs,
                residuals,
                weights,
                params.ms1_huber_delta,
                params.ms1_lambda,
                params.max_iter_newton,
                params.max_iter_bisection,
                params.max_iter_outer,
                search_context.deconvolution_stop_tolerance[],
                search_context.deconvolution_stop_tolerance[],
                params.max_diff,
                params.ms1_reg_type
            )
```

**Add immediately after**:
```julia
            # NEW: Distribute grouped coefficients back to individual precursors
            distribute_ms1_coefficients!(
                weights,  # Individual precursor weights (modified in-place)
                weights,  # Fitted group weights (same array, different indexing)
                mz_grouping,
                precs_temp[1:prec_temp_size]
            )
```

### Step 5: Integrate into SecondPassSearch

**File**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

**Location**: In the MS1 deconvolution section (around line 476-491)

**Changes needed**:

1. **Add m/z grouping initialization** after the existing variable declarations in the MS1 section:

```julia
    # NEW: Add after existing variable declarations in MS1 section (around line 445)
    mz_grouping = MzGroupingMap(UInt32(100000))  # 5 decimal place precision
    precursors = getPrecursors(getSpecLib(search_context))
```

2. **Replace the buildDesignMatrix! call** around line 99 in the MS1 section:

**Find this block** (in the MS1 chromatogram function):
```julia
            buildDesignMatrix!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                getIdToCol(search_data)
            )
```

**Replace with**:
```julia
            # Reset grouping for this scan
            reset!(mz_grouping)

            # Use MS1-specific design matrix construction with m/z grouping
            buildDesignMatrixMS1!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                mz_grouping,
                precursors
            )
```

3. **Add coefficient distribution** after the MS1 `solveHuber!` call (around line 491):

**Find this block**:
```julia
            solveHuber!(
                Hs,
                residuals,
                weights,
                params.ms1_huber_delta,
                params.ms1_lambda,
                params.max_iter_newton,
                params.max_iter_bisection,
                params.max_iter_outer,
                search_context.deconvolution_stop_tolerance[],
                search_context.deconvolution_stop_tolerance[],
                params.max_diff,
                params.ms1_reg_type
            )
```

**Add immediately after**:
```julia
            # NEW: Distribute grouped coefficients back to individual precursors
            # Note: Need to get current precursor list for this scan
            current_precursors = UInt32[]
            for (prec_id, col_idx) in pairs(getIdToCol(search_data))
                if col_idx > 0
                    push!(current_precursors, prec_id)
                end
            end

            distribute_ms1_coefficients!(
                weights,  # Individual precursor weights (modified in-place)
                weights,  # Fitted group weights (same array, different indexing)
                mz_grouping,
                current_precursors
            )
```

## Expected Behavior Changes

### Before Implementation
- Target/decoy pairs with identical masses create separate design matrix columns
- Multicollinearity can cause numerical instability in MS1 fitting
- Each precursor gets its own coefficient regardless of mass similarity

### After Implementation
- Precursors with identical m/z (rounded to 5 decimal places) share design matrix columns
- Fitted coefficients for the group are distributed to all precursors in the group
- Reduced multicollinearity improves numerical stability
- Individual precursor tracking is preserved for downstream analysis

## Testing Strategy

### Unit Tests
Create tests in `test/UnitTests/` directory:

```julia
# Test m/z grouping logic
@testset "MZ Grouping" begin
    grouping = MzGroupingMap(100000)

    # Test identical masses get grouped
    mz1 = round_mz_for_grouping(1000.12345f0, 100000)
    mz2 = round_mz_for_grouping(1000.12345f0, 100000)
    @test mz1 == mz2

    # Test different masses don't get grouped
    mz3 = round_mz_for_grouping(1000.12346f0, 100000)
    @test mz1 != mz3
end

# Test coefficient distribution
@testset "Coefficient Distribution" begin
    # Create mock data and test distribution logic
    # Verify that grouped coefficients are correctly assigned
end
```

### Integration Tests
- Run small dataset through both old and new MS1 fitting
- Compare quantification results for precursors with/without identical masses
- Verify that multicollinearity warnings are reduced

## Performance Considerations

### Memory Usage
- `MzGroupingMap` adds ~3 Dict structures per scan
- Overhead should be minimal compared to existing design matrix construction
- Memory scales with number of unique m/z groups (typically << number of precursors)

### Computational Overhead
- Additional m/z rounding and dictionary lookups per precursor
- Coefficient distribution is O(n) where n = precursors per scan
- Overall overhead should be <5% of total deconvolution time

### Scalability
- Grouping efficiency depends on precision choice
- With 100,000 precision, groups ~0.01 ppm at 1000 Da
- Can be made configurable if needed

## Configuration Options (Future)

Consider adding to `PioneerParameters` structure:

```julia
# In global_settings or optimization.deconvolution.ms1
"mz_grouping": {
    "enabled": true,
    "precision": 100000,  # 5 decimal places
    "min_group_size": 2   # Only group if ≥2 precursors
}
```

## Validation

### Success Criteria
1. **Functional**: MS1 quantification produces reasonable coefficients
2. **Numerical**: Reduced condition number in design matrices with identical masses
3. **Accuracy**: Quantification accuracy maintained or improved vs. original method
4. **Performance**: <10% overhead in MS1 deconvolution time

### Potential Issues and Solutions

**Issue**: M/Z precision too coarse, groups distinct precursors
**Solution**: Increase precision or add validation to check actual m/z differences

**Issue**: Memory overhead in large experiments
**Solution**: Clear grouping maps more frequently or use more memory-efficient structures

**Issue**: Coefficient distribution logic errors
**Solution**: Add comprehensive bounds checking and validation

**Issue**: Integration complexity
**Solution**: Start with IntegrateChromatogramsSearch only, add SecondPassSearch after validation

## Implementation Priority

1. **Phase 1**: Implement core data structures and MS1 design matrix function
2. **Phase 2**: Integrate into IntegrateChromatogramsSearch only
3. **Phase 3**: Add comprehensive testing and validation
4. **Phase 4**: Integrate into SecondPassSearch
5. **Phase 5**: Add configuration options and optimization

This approach allows for incremental implementation and testing at each phase.