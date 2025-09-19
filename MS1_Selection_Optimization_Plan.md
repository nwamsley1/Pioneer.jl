# MS1 PSM Selection Optimization Plan

## Problem Statement

The current MS1 feature selection logic has both performance and biological accuracy issues:

1. **Quadratic Performance**: Current RT proximity search is O(n²) - for each MS1 PSM group, it searches through all MS2 PSMs
2. **Selection Strategy Conflict**: Two competing approaches exist:
   - `parseMs1Psms`: Selects max intensity MS1 scan per precursor
   - Current RT proximity: Selects MS1 scan closest to MS2 apex RT

## Proposed Hybrid Solution

**Goal**: Combine the benefits of both approaches while achieving O(n log n) performance.

### Strategy Overview
1. **Track max intensity MS1 RT** for comparison with MS2 apex RT
2. **Use RT-proximate MS1 scan** for all other features
3. **Efficient alignment** using pre-built lookup structures

## Detailed Implementation Plan

### Step 1: Pre-process MS2 PSM RT Lookup

**Create efficient MS2 RT lookup structure:**
```julia
# Build precursor_idx → MS2_RT mapping (O(n))
ms2_rt_lookup = Dict{UInt32, Float32}()
for row in eachrow(psms)
    ms2_rt_lookup[row.precursor_idx] = row.rt
end
```

**Benefits:**
- O(1) RT lookup instead of O(n) search
- Eliminates quadratic behavior
- Pre-computed once, used many times

### Step 2: Enhanced MS1 PSM Processing

**Replace `parseMs1Psms` with hybrid selection:**
```julia
function parseMs1PsmsHybrid(
    psms::DataFrame,
    spectra::MassSpecData,
    ms2_rt_lookup::Dict{UInt32, Float32}
)
    if !hasproperty(psms, :precursor_idx) || (size(psms, 1) == 0)
        return DataFrame()
    end

    # Add RT column
    rts = [getRetentionTime(spectra, scan_idx) for scan_idx in psms.scan_idx]
    psms[!, :rt] = rts

    # Group by precursor and apply hybrid selection
    return combine(groupby(psms, :precursor_idx)) do group
        # Find max intensity scan
        max_intensity_idx = argmax(group[!, :weight])
        max_intensity_rt = group[max_intensity_idx, :rt]

        # Find RT-closest scan to MS2 apex
        precursor_idx = group[1, :precursor_idx]
        if haskey(ms2_rt_lookup, precursor_idx)
            ms2_rt = ms2_rt_lookup[precursor_idx]
            rt_diffs = abs.(group[!, :rt] .- ms2_rt)
            closest_rt_idx = argmin(rt_diffs)

            # Start with RT-closest scan features
            result_row = group[closest_rt_idx, :]

            # Add max intensity RT as additional feature
            result_row[!, :rt_max_intensity] = max_intensity_rt
            result_row[!, :rt_diff_max_intensity] = abs(max_intensity_rt - ms2_rt)

            return result_row
        else
            # Fallback to max intensity if no MS2 match
            return group[max_intensity_idx, :]
        end
    end
end
```

### Step 3: Feature Engineering Benefits

**New features available for ML models:**
```julia
# Existing features from RT-closest scan:
:rt_ms1                    # RT of scan closest to MS2 apex
:weight_ms1               # Intensity features from RT-closest scan
:gof_ms1                  # Spectral quality from RT-closest scan
:error_ms1                # Mass accuracy from RT-closest scan

# New features from max intensity scan:
:rt_max_intensity_ms1     # RT of strongest MS1 signal
:rt_diff_max_intensity_ms1 # |RT_max_intensity - RT_ms2_apex|

# Derived features:
:rt_diff_ms1 = abs(:rt_ms1 - :rt_ms2)  # RT alignment quality
:intensity_vs_apex_ms1 = :weight_ms1 / :weight_max_intensity_ms1  # Relative intensity
```

### Step 4: Performance Optimization

**Algorithm Complexity:**
- **Current**: O(n²) - n MS1 groups × n MS2 PSM searches
- **Proposed**: O(n log n) - Dict lookup O(1) + groupby operations O(n log n)

**Memory Efficiency:**
- Pre-compute MS2 RT lookup once: ~4 bytes per MS2 PSM
- Eliminate repeated DataFrame filtering operations
- Single-pass MS1 processing

### Step 5: Biological Interpretation Benefits

**RT-Closest Features (Primary)**:
- Better temporal alignment with MS2 identification
- More accurate mass error estimation
- Correct spectral quality metrics
- Proper intensity quantification

**Max Intensity Features (Secondary)**:
- Peak detection confirmation
- Chromatographic consistency validation
- Signal-to-noise assessment
- RT drift detection

## Implementation Steps

### Phase 1: Core Algorithm
1. **Replace `parseMs1Psms`** with `parseMs1PsmsHybrid`
2. **Build MS2 RT lookup** before MS1 processing
3. **Update feature column names** to reflect dual approach

### Phase 2: Feature Integration
4. **Add new features** to ML model feature lists
5. **Update feature importance** analysis
6. **Modify QC plots** to show both RT alignments

### Phase 3: Performance Validation
7. **Benchmark performance** on large datasets
8. **Profile memory usage** improvements
9. **Validate biological accuracy** with known standards

## Code Location Changes

**Files to Modify:**
1. `SecondPassSearch/utils.jl` - Replace `parseMs1Psms` function
2. `SecondPassSearch/SecondPassSearch.jl` - Update MS1 processing call
3. `ScoringSearch/` - Add new features to ML model
4. `WriteOutputs/qcPlots.jl` - Update QC visualizations

**Key Changes:**
```julia
# In SecondPassSearch.jl (lines 340-344):
# OLD:
ms1_psms = parseMs1Psms(ms1_psms, spectra)

# NEW:
ms2_rt_lookup = Dict{UInt32, Float32}(
    row.precursor_idx => row.rt for row in eachrow(psms)
)
ms1_psms = parseMs1PsmsHybrid(ms1_psms, spectra, ms2_rt_lookup)
```

## Expected Outcomes

### Performance Improvements
- **10-100x faster** MS1 processing for large datasets
- **Reduced memory** allocation from DataFrame operations
- **Better scalability** for high-resolution MS data

### Biological Accuracy
- **Improved feature quality** from proper RT alignment
- **Enhanced peak detection** from max intensity tracking
- **Better ML model performance** from richer feature set

### Analytical Benefits
- **RT consistency validation** across MS1/MS2
- **Chromatographic quality assessment**
- **Peak shape characterization**
- **Cross-run alignment validation**

## Risk Mitigation

**Potential Issues:**
1. **Memory overhead** from RT lookup Dict → Minimal (4 bytes/PSM)
2. **Feature naming conflicts** → Careful column naming convention
3. **Backward compatibility** → Maintain same output schema initially

**Testing Strategy:**
1. **Unit tests** for `parseMs1PsmsHybrid` function
2. **Integration tests** with existing SecondPassSearch pipeline
3. **Performance benchmarks** on ecoli_test dataset
4. **Biological validation** with known retention time standards

## Migration Path

### Phase 1: Core Implementation (Week 1)
- Implement `parseMs1PsmsHybrid` function
- Add RT lookup construction
- Update SecondPassSearch integration

### Phase 2: Feature Enhancement (Week 2)
- Add new features to ML pipeline
- Update feature importance analysis
- Enhance QC plotting

### Phase 3: Optimization & Validation (Week 3)
- Performance profiling and optimization
- Biological validation studies
- Documentation and testing

This approach provides the best of both worlds: **biological accuracy** from RT proximity and **signal quality** from max intensity tracking, while achieving significant **performance improvements** through algorithmic optimization.