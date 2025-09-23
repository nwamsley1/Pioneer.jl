# MS1 M/Z Grouping Implementation Plan

## Overview

This document outlines two approaches for combining precursors with identical m/z values into the same column of the design matrix during MS1 fitting. This addresses multicollinearity issues caused by target/decoy pairs and other precursors sharing identical masses while allowing individual coefficient assignment back to each precursor.

## Current State Analysis

### Current MS1 Design Matrix Construction

**Location**: `IntegrateChromatogramsSearch/utils.jl:538-545` and `SecondPassSearch/utils.jl:472-490`

**Current Flow**:
```julia
# 1. Build design matrix with individual precursor columns
buildDesignMatrix!(Hs, ion_matches, ion_misses, nmatches, nmisses, getIdToCol(search_data))

# 2. Each precursor gets its own column based on precursor ID
if iszero(precID_to_col[getPrecID(match)])
    prec_col += one(UInt32)  # New column for each unique precursor ID
    update!(precID_to_col, getPrecID(match), UInt16(prec_col))
end
```

**Problem**: Target/decoy pairs with identical masses create separate but potentially identical columns, causing multicollinearity.

### Data Structures Involved

**Key Components**:
- `precID_to_col::ArrayDict{UInt32, UInt16}` - Maps precursor ID → design matrix column
- `getPrecID(match)` - Returns unique precursor identifier
- `getMz(getPrecursors(spec_lib))[precursor_id]` - Returns m/z for precursor
- `weights[column_idx]` - Fitted coefficients per design matrix column

## Proposed Solutions

### Approach 1: M/Z-Based Column Mapping (Recommended)

**Core Concept**: Group precursors by rounded m/z values into shared design matrix columns, then distribute fitted coefficients back to individual precursors.

#### Implementation Details

**1. Add M/Z Grouping Data Structure**

Create new data structure to manage m/z-based grouping:

```julia
struct MzGroupingMap
    # Maps rounded m/z → design matrix column
    mz_to_col::Dict{UInt32, UInt16}

    # Maps precursor ID → m/z group (for coefficient retrieval)
    precid_to_mz_group::Dict{UInt32, UInt32}

    # Maps m/z group → list of precursor IDs in that group
    mz_group_to_precids::Dict{UInt32, Vector{UInt32}}

    # Precision for m/z rounding (e.g., 100000 for 5 decimal places)
    mz_precision::UInt32

    # Current column counter
    current_col::UInt16
end

function MzGroupingMap(mz_precision::UInt32 = 100000)
    return MzGroupingMap(
        Dict{UInt32, UInt16}(),
        Dict{UInt32, UInt32}(),
        Dict{UInt32, Vector{UInt32}}(),
        mz_precision,
        zero(UInt16)
    )
end

# Round m/z to specified precision for grouping
function round_mz_for_grouping(mz::Float32, precision::UInt32)::UInt32
    return round(UInt32, mz * precision)
end
```

**2. Modified buildDesignMatrix! for MS1**

Create MS1-specific version that uses m/z grouping:

```julia
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

    # Get precursor m/z values for grouping
    precursor_mzs = getMz(precursors)

    # Same row counting logic as original
    M = 1
    for i in range(2, nmatches)
        if getPeakInd(matches[i]) != getPeakInd(matches[i - 1])
            M += 1
        end
    end
    M += nmisses

    # Array size management (same as original)
    if (nmatches + nmisses) >= length(H.colval) - 1
        # ... resize logic ...
    end

    # Modified column assignment using m/z grouping
    row = 0
    last_peak_ind = -1
    last_row, last_col = -1, -1
    j = 0
    H.n_vals = 0

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

        # Peak indexing logic (same as original)
        if getPeakInd(match) != last_peak_ind
            row += 1
            last_peak_ind = getPeakInd(match)
        end

        # Add to design matrix (same as original)
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

    # Miss handling (similar modifications)
    # ... miss processing logic with m/z grouping ...

    sortSparse!(H)
end
```

**3. Coefficient Distribution**

After fitting, distribute fitted coefficients back to individual precursors:

```julia
function distribute_ms1_coefficients!(
    individual_weights::Vector{Float32},
    fitted_weights::Vector{Float32},
    mz_grouping::MzGroupingMap,
    precursor_ids::Vector{UInt32}
)
    for prec_id in precursor_ids
        if haskey(mz_grouping.precid_to_mz_group, prec_id)
            mz_group = mz_grouping.precid_to_mz_group[prec_id]
            group_col = mz_grouping.mz_to_col[mz_group]

            # Assign fitted coefficient to this precursor
            individual_weights[prec_id] = fitted_weights[group_col]
        else
            # Precursor not in any group (shouldn't happen)
            individual_weights[prec_id] = 0.0f0
        end
    end
end
```

**4. Integration Points**

**IntegrateChromatogramsSearch Modification**:
```julia
# In build_chromatograms() for MS1CHROM
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS1CHROM
)
    # ... existing setup code ...

    # NEW: Create m/z grouping map for MS1
    mz_grouping = MzGroupingMap(100000)  # 5 decimal place precision
    precursors = getPrecursors(getSpecLib(search_context))

    for scan_idx in scan_range
        # ... existing scan processing ...

        if nmatches > 2
            # Use MS1-specific design matrix construction
            buildDesignMatrixMS1!(
                Hs,
                ion_matches,
                ion_misses,
                nmatches,
                nmisses,
                mz_grouping,
                precursors
            )

            # ... weight initialization ...

            # Solve with grouped design matrix
            solveHuber!(Hs, residuals, grouped_weights, ...)

            # Distribute coefficients back to individual precursors
            distribute_ms1_coefficients!(
                weights,  # Individual precursor weights
                grouped_weights,  # Fitted group weights
                mz_grouping,
                precs_temp[1:prec_temp_size]
            )

            # ... existing chromatogram recording ...
        end
    end
end
```

#### Advantages
- **Eliminates multicollinearity** by combining identical masses
- **Preserves individual precursor tracking** for downstream analysis
- **Minimal performance impact** - only adds grouping overhead
- **Configurable precision** for m/z rounding tolerance

#### Potential Issues
- **M/Z precision choice**: Too coarse groups distinct precursors, too fine misses identical masses
- **Memory overhead**: Additional data structures for grouping
- **Complexity**: More complex coefficient distribution logic

---

### Approach 2: Post-Fitting Regularization Adjustment

**Core Concept**: Keep individual columns but apply post-fitting regularization specifically to groups of precursors with identical m/z values.

#### Implementation Details

**1. Identify M/Z Groups After Design Matrix Construction**

```julia
function identify_mz_groups(
    precID_to_col::ArrayDict{UInt32, UInt16},
    precursors::LibraryPrecursors,
    mz_precision::Float32 = 1e-5
)::Dict{UInt32, Vector{UInt32}}

    mz_groups = Dict{UInt32, Vector{UInt32}}()
    precursor_mzs = getMz(precursors)

    for (prec_id, col) in pairs(precID_to_col)
        rounded_mz = round(UInt32, precursor_mzs[prec_id] / mz_precision)

        if !haskey(mz_groups, rounded_mz)
            mz_groups[rounded_mz] = UInt32[]
        end
        push!(mz_groups[rounded_mz], prec_id)
    end

    # Return only groups with multiple precursors
    return filter(p -> length(p.second) > 1, mz_groups)
end
```

**2. Post-Fitting Coefficient Averaging**

```julia
function regularize_identical_mz_coefficients!(
    weights::Vector{Float32},
    precID_to_col::ArrayDict{UInt32, UInt16},
    mz_groups::Dict{UInt32, Vector{UInt32}}
)
    for (rounded_mz, prec_ids) in mz_groups
        if length(prec_ids) > 1
            # Calculate average coefficient for this m/z group
            group_columns = [precID_to_col[pid] for pid in prec_ids]
            group_weights = [weights[col] for col in group_columns]
            avg_weight = mean(group_weights)

            # Assign average to all precursors in group
            for col in group_columns
                weights[col] = avg_weight
            end
        end
    end
end
```

**3. Integration**

```julia
# After solveHuber! in MS1 fitting
solveHuber!(Hs, residuals, weights, ...)

# NEW: Apply m/z-based regularization
mz_groups = identify_mz_groups(getIdToCol(search_data), precursors)
regularize_identical_mz_coefficients!(weights, getIdToCol(search_data), mz_groups)

# Continue with existing logic...
```

#### Advantages
- **Minimal code changes** - post-processing approach
- **Preserves existing infrastructure** - no changes to design matrix construction
- **Simple implementation** - straightforward averaging logic

#### Disadvantages
- **Less principled** - post-hoc averaging vs. proper grouped fitting
- **Doesn't address multicollinearity** during fitting - may still have numerical issues
- **Potential information loss** - ignores correlation structure during optimization

## Questions for Clarification

1. **M/Z Precision**: What precision should be used for grouping m/z values? Options:
   - 5 decimal places (10^-5 Da): Groups masses within 0.01 ppm at 1000 Da
   - 4 decimal places (10^-4 Da): Groups masses within 0.1 ppm at 1000 Da
   - PPM-based grouping: Group within X ppm tolerance

2. **Scope of Application**:
   - Apply only to target/decoy pairs with identical masses?
   - Apply to all precursors with identical masses regardless of origin?
   - How to handle cases where >2 precursors share the same mass?

3. **Performance Requirements**:
   - Is the additional computational overhead of Approach 1 acceptable?
   - Should this be optional/configurable via parameters?

4. **Validation Strategy**:
   - How to test that grouping correctly identifies identical masses?
   - How to validate that coefficient distribution preserves quantitative accuracy?

5. **Edge Cases**:
   - What happens if a precursor group has all zero coefficients?
   - Should grouping respect charge states (same mass but different charge)?
   - How to handle isotope patterns within grouped precursors?

## Recommended Implementation Path

**Phase 1**: Implement Approach 1 (M/Z-Based Column Mapping) for MS1 fitting only
- Start with 5 decimal place precision (10^-5 Da)
- Add comprehensive unit tests for grouping logic
- Apply only in `IntegrateChromatogramsSearch` and `SecondPassSearch` MS1 paths

**Phase 2**: Add configuration parameters
- Make m/z precision configurable via `PioneerParameters`
- Add flag to enable/disable m/z grouping
- Add diagnostic output showing grouping statistics

**Phase 3**: Validation and optimization
- Compare quantification accuracy with/without grouping
- Optimize data structures for large-scale experiments
- Consider extending to MS2 if beneficial

This approach provides a principled solution to the multicollinearity problem while maintaining the flexibility to track individual precursors for downstream analysis.