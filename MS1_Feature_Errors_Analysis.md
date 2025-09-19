# Critical MS1 Feature Errors Analysis

## Overview

This document identifies significant errors and missed opportunities in Pioneer.jl's MS1 features calculation pipeline, with particular focus on MS1-MBR integration and feature joining logic.

## Critical MS1 Feature Errors Found 🚨

### **1. Major MS1-MBR Integration Gap** ⚠️

**Current MBR features are MS2-only**:
```julia
# Current MBR features (percolatorSortOf.jl:717-744)
MBR_num_runs, MBR_best_irt_diff, MBR_rv_coefficient,
MBR_log2_weight_ratio, MBR_max_pair_prob, MBR_is_missing
```

**Missing MS1 MBR features**:
- **MS1 peak presence across runs** - Strong indicator of true positive
- **MS1 isotope pattern consistency** - Quality metric across runs
- **MS1 RT alignment quality** - Independent RT validation
- **MS1 mass accuracy consistency** - Systematic error detection

### **2. Flawed MS1 Feature Join Logic** 🚨

**Problem** (`SecondPassSearch.jl:395-426`):
```julia
# MS1 features joined ONLY by precursor_idx - no RT consideration!
psms = leftjoin(psms, ms1_psms, on = :precursor_idx)
miss_mask = ismissing.(psms[!, ms1_cols[1]])
```

**What Actually Happens Upstream**:

1. **MS1 PSM Collection** (`SecondPassSearch/utils.jl:309-531`):
   ```julia
   # Each MS1 scan produces individual Ms1ScoredPSM records
   struct Ms1ScoredPSM {
       precursor_idx::UInt32
       scan_idx::UInt32      # Different scans for same precursor!
       weight::Float32       # Individual scan intensities
       spectral_contrast::Float16
       gof::Float16
       rt::Float32           # Scan retention time
       # ... other features
   }
   ```

2. **Multiple Records Per Precursor**: MS1 search creates **one PSM per scan per precursor**
   - Same `precursor_idx` appears in multiple scans across chromatographic peak
   - Each scan has different `scan_idx`, `rt`, `weight`, and spectral scores
   - **No aggregation or apex selection** for MS1 PSMs

3. **Problematic Join** (`SecondPassSearch.jl:397-403`):
   ```julia
   # MS2 PSMs (1 per precursor apex) joined to MS1 PSMs (many per precursor)
   psms = leftjoin(psms, ms1_psms, on = :precursor_idx, makeunique = true)
   ```

**Critical Issues Identified**:

- **Arbitrary MS1 Scan Selection**: `leftjoin` picks the **first** MS1 PSM with matching precursor_idx
- **No Apex Selection**: Unlike MS2 which keeps only `best_scan` apex PSMs (line 393), MS1 has no apex filtering
- **RT Mismatch**: MS1 features may come from scan far from MS2 PSM retention time
- **Peak Quality Ignored**: May select weak MS1 signal instead of chromatographic apex
- **Data Loss**: Discards potentially better MS1 measurements from other scans

**Concrete Example**:
```
MS2 PSM: precursor_456, scan_1050, rt=25.4min (apex)
MS1 PSMs:
  - precursor_456, scan_1045, rt=25.1min, weight=1000  # weak
  - precursor_456, scan_1048, rt=25.3min, weight=5000  # strong
  - precursor_456, scan_1052, rt=25.5min, weight=3000  # decent

Result: MS2 PSM gets MS1 features from scan_1045 (first match, weak signal)
Should get: MS1 features from scan_1048 (strongest signal, close RT)
```

### **3. Inadequate RT Alignment** ⚠️

**Evidence** (`SecondPassSearch.jl:259-263`):
```julia
# Commented out MS1 iRT calculations!
#ms1_psms[!,:irt] = zeros(Float32, size(ms1_psms, 1))
#    ms1_psms[i,:irt] = rt_irt_model(getRetentionTime(spectra, ms1_psms[i,:scan_idx]))
```

**Problems**:
- **No MS1-specific RT calibration** - Uses ad-hoc RT differences
- **Missing iRT validation** - MS1 RT consistency not checked
- **Poor temporal matching** - MS1/MS2 scan timing not validated

### **4. MS1 Missing Data Handling** ❌

**Issue** (`SecondPassSearch.jl:412-415`):
```julia
if size(ms1_psms, 1) == 0
    miss_mask = trues(size(psms, 1))  # ALL PSMs marked as missing MS1!
    psms[!, col] = -1*ones(Float32, size(psms, 1))  # Default -1 values
end
```

**Problems**:
- **Cascading failures** - No MS1 PSMs → all features marked missing
- **Poor fallback values** - Default -1 may confuse ML models
- **No partial recovery** - Can't use even basic MS1 scan information

### **5. Missing MS1 MBR Opportunities** 💡

**Proposed MS1 MBR features**:

```julia
# Proposed MS1 MBR features:
MBR_ms1_peak_consistency    // MS1 peak found across runs at expected RT
MBR_ms1_isotope_similarity  // Isotope pattern correlation across runs
MBR_ms1_mass_accuracy      // MS1 mass error consistency
MBR_ms1_rt_precision       // MS1 RT alignment quality
MBR_ms1_signal_strength    // MS1 peak intensity reliability
```

## **Recommended Fixes**

### **1. Add MS1 MBR Features**
```julia
# In MBR calculation (percolatorSortOf.jl)
MBR_ms1_peak_found = check_ms1_peak_at_rt(precursor, expected_rt, tolerance)
MBR_ms1_isotope_quality = calculate_ms1_pattern_consistency(isotope_matches)
MBR_ms1_rt_accuracy = abs(observed_ms1_rt - predicted_rt)
```

### **2. Improve MS1 Feature Joining**
```julia
# Add RT proximity constraint
psms = leftjoin(psms, ms1_psms,
    on = [:precursor_idx, :rt => :rt_ms1],  # RT proximity join
    matchmissing = :notequal,
    rt_tolerance = 0.5  # minutes
)
```

### **3. Implement MS1 RT Calibration**
```julia
# Enable proper MS1 iRT calculation
for i in 1:size(ms1_psms, 1)
    ms1_psms[i,:irt] = rt_irt_model(getRetentionTime(spectra, ms1_psms[i,:scan_idx]))
    ms1_psms[i,:irt_error] = abs(ms1_psms[i,:irt] - expected_irt)
end
```

### **4. Enhance Missing Data Strategy**
```julia
# Graduated missing data handling
ms1_quality_score = calculate_ms1_completeness(ms1_matches, rt_proximity, mass_accuracy)
ms1_confidence_level = categorize_ms1_quality(ms1_quality_score)  # high/medium/low/missing
```

## **Impact Assessment**

These errors likely cause:
- **Reduced sensitivity** from poor MS1 feature matching
- **Lower specificity** from incorrect MS1-MS2 associations
- **Missed identifications** from inadequate MBR with MS1 validation
- **Suboptimal scoring** from poor MS1 feature quality

The MS1 MBR integration could significantly improve identification confidence, especially for challenging cases where MS2 evidence is ambiguous but MS1 provides strong supporting evidence.

## Code References

**MS1 Implementation**:
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl:395-426` - MS1 feature joining logic
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl:259-263` - Commented RT calibration
- `src/utils/ML/percolatorSortOf.jl:717-744` - Current MBR feature calculation
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl` - Feature configuration

**MS1 Feature Pipeline**:
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:372-387` - MS1 error collection
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl:544-604` - MS1 tolerance estimation
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl:327-886` - MS1 feature integration