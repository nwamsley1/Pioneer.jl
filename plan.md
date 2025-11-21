# Plan: Fix rt_to_irt vs rt_to_refined_irt Usage Throughout Pipeline

## Problem Statement

After implementing iRT refinement, we have **TWO sets of RT↔iRT conversion models**, but many parts of the code are using the **wrong ones**:

1. **Library iRT models** (`rt_irt_map`, `irt_rt_map`): RT ↔ library_iRT (from spectral library)
2. **Refined iRT models** (`rt_to_refined_irt_map`, `refined_irt_to_rt_map`): RT ↔ refined_iRT (library adjusted by refinement)

### Critical Understanding

**Refined iRT is FILE-SPECIFIC and must be calculated on-the-fly:**
- **Observed refined iRT** = `rt_to_refined_irt_map[file_idx](scan_rt)` - converts scan RT to refined iRT space
- **Predicted refined iRT** = `refinement_model[file_idx](sequence, library_irt)` - applies sequence-specific correction to library iRT

**Cannot store a single "refined iRT" value per precursor globally** because:
- Each MS file has its own refinement model
- Same precursor has different refined iRT in different files
- Must be computed per-file using file-specific models

### Current Issues Identified

1. **getRtIrtModel() returns LIBRARY iRT models** everywhere:
   - SecondPassSearch (6+ locations)
   - HuberTuningSearch (1 location)
   - IntegrateChromatogramsSearch (2 locations)
   - **Should use refined iRT models** after FirstPassSearch

2. **No getter for refined models by file index**
   - Have: `getRtToRefinedIrtMap(s)` → returns whole dict
   - Need: `getRtToRefinedIrtModel(s, index)` → returns model for specific file

3. **getPredIrt/setPredIrt are effectively obsolete**
   - Currently store library iRT values per precursor
   - Only used in 2 places (SecondPassSearch/utils.jl lines 905, 909)
   - Should be replaced with direct library iRT access: `getIrt(precursors)[prec_idx]`
   - Or kept as-is since they return library iRT (which is correct)

4. **FirstPassSearch.jl stores library iRT in irt_obs**
   - Lines 577, 587, 589, 595 set library iRT values
   - This is actually **CORRECT** - irt_obs should store library iRT
   - getPredIrt returns library iRT (which is what SecondPassSearch needs for the refinement calculation)

## Solution Approach

**Add `getRtToRefinedIrtModel()` with fallback** + Update all post-FirstPassSearch code to use refined RT→iRT models

**DO NOT change irt_obs or getPredIrt/setPredIrt** - they correctly store/return library iRT

**Files to modify:** 6 files (not 7), ~15 locations (not 20)

---

## Implementation Steps

### Part 1: Add getRtToRefinedIrtModel() Getters

**File:** SearchTypes.jl
**Location:** After line 508

**Add these two functions:**

```julia
"""
    getRtToRefinedIrtModel(s::SearchContext, index::Integer)

Get RT → refined_iRT model for file index.
Falls back to library iRT if refined unavailable.
Returns identity model if neither exists.

# Usage
observed_refined_irt = getRtToRefinedIrtModel(context, file_idx)(scan_rt)
"""
function getRtToRefinedIrtModel(s::SearchContext, index::I) where {I<:Integer}
    if haskey(s.rt_to_refined_irt_map, index)
        return s.rt_to_refined_irt_map[index]
    elseif haskey(s.rt_irt_map, index)
        @debug "Refined iRT model not found for file $index, falling back to library iRT model"
        return s.rt_irt_map[index]
    else
        return IdentityModel()
    end
end

"""
    getRefinedIrtToRtModel(s::SearchContext, index::Integer)

Get refined_iRT → RT model for file index.
Falls back to library iRT if refined unavailable.
Returns identity model if neither exists.

# Usage
predicted_rt = getRefinedIrtToRtModel(context, file_idx)(refined_irt)
"""
function getRefinedIrtToRtModel(s::SearchContext, index::I) where {I<:Integer}
    if haskey(s.refined_irt_to_rt_map, index)
        return s.refined_irt_to_rt_map[index]
    elseif haskey(s.irt_rt_map, index)
        @debug "Refined iRT model not found for file $index, falling back to library iRT model"
        return s.irt_rt_map[index]
    else
        return IdentityModel()
    end
end
```

---

### Part 2: Update SecondPassSearch

**File:** SecondPassSearch/SecondPassSearch.jl

**Line 495:**
```julia
# FROM: getRtIrtModel(search_context, ms_file_idx)
# TO:   getRtToRefinedIrtModel(search_context, ms_file_idx)
```

**Lines 550-553:**
```julia
# FROM:
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                      Float32(-1),
                      abs.(rt_to_irt_model.(psms[!,:rt]) .- rt_to_irt_model.(psms[!,:rt_ms1]))))

# TO:
rt_to_refined_irt_model = getRtToRefinedIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                      Float32(-1),
                      abs.(rt_to_refined_irt_model.(psms[!,:rt]) .- rt_to_refined_irt_model.(psms[!,:rt_ms1]))))
```

**Line 570:**
```julia
# FROM: getRtIrtModel(search_context, ms_file_idx),
# TO:   getRtToRefinedIrtModel(search_context, ms_file_idx),
```

**File:** SecondPassSearch/utils.jl

**Lines 181-183:**
```julia
# FROM:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# TO:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Also update line 190:**
```julia
# FROM: if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) ||
# TO:   if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) ||
# (no change needed - uses irt_start/irt_stop which are different variables)
```

**Lines 386-388:**
```julia
# FROM:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# TO:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Line 834** (parameter name for clarity):
```julia
# FROM: rt_to_irt_interp::RtConversionModel,
# TO:   rt_to_refined_irt_interp::RtConversionModel,
```

**Lines 904, 909** (update variable name):
```julia
# FROM: rt_to_irt_interp
# TO:   rt_to_refined_irt_interp
```

Note: Detailed changes for computing refined iRT predictions are shown in the "Detailed SecondPassSearch/utils.jl Changes" section below.

**Line 1020** (parameter name):
```julia
# FROM: rt_to_irt_interp::RtConversionModel
# TO:   rt_to_refined_irt_interp::RtConversionModel
```

**Line 1072:**
```julia
# FROM: irts = rt_to_irt_interp.(psms.rt)
# TO:   irts = rt_to_refined_irt_interp.(psms.rt)
```

---

### Part 3: Update HuberTuningSearch

**File:** HuberTuningSearch/utils.jl

**Lines 225-227:**
```julia
# FROM:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# TO:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Also update line 234:**
```julia
# FROM: if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
# TO:   if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
# (no change needed - uses irt_start/irt_stop which are different variables)
```

---

### Part 4: Update IntegrateChromatogramsSearch

**File:** IntegrateChromatogramsSearch/utils.jl

**Lines 254-256:**
```julia
# FROM:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# TO:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Also update line 263:**
```julia
# FROM: if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
# TO:   if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
# (no change needed - uses irt_start/irt_stop which are different variables)
```

**Lines 491-493:**
```julia
# FROM:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# TO:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

---

### Part 5: Keep FirstPassSearch as-is

**File:** FirstPassSearch/FirstPassSearch.jl
**Lines 575-597** - NO CHANGES

**Why:** This code stores library iRT values in irt_obs, which is correct. getPredIrt/setPredIrt deal with library iRT, which is the base value needed for refinement.

**File:** FirstPassSearch/utils.jl
**Line 98:** `irt[i] = rt_irt(rt[i])`

**Add clarifying comment:**
```julia
# Line 98:
# Note: Using library iRT here since refinement happens after FirstPassSearch completes
irt[i] = rt_irt(rt[i])
```

---

## Detailed SecondPassSearch/utils.jl Changes

**Location:** add_features! function

**Lines 853, 855** (variable allocation):
```julia
# FROM:
irt_obs = zeros(Float32, N)
irt_pred = zeros(Float32, N)

# TO:
refined_irt_obs = zeros(Float32, N)
refined_irt_pred = zeros(Float32, N)
```

**Before line 895** (after line 894, before tasks_per_thread):
```julia
# Get refinement model for this file (will be used in parallel loop)
refinement_model = getIrtRefinementModel(search_context, ms_file_idx)
```

**Lines 904-913** (inside the parallel loop):

**CURRENT CODE:**
```julia
irt_obs[i] = rt_to_irt_interp(rt[i])
irt_pred[i] = getPredIrt(search_context, prec_idx)
#irt_diff[i] = abs(irt_obs[i] - first(prec_id_to_irt[prec_idx]))
irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_refined_irt)
if !ms1_missing[i]
    ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx))
else
    ms1_irt_diff[i] = 0f0
end
irt_error[i] = abs(irt_obs[i] - irt_pred[i])
```

**CHANGED CODE:**
```julia
# Calculate observed refined iRT from scan RT
refined_irt_obs[i] = rt_to_refined_irt_interp(rt[i])

# Calculate predicted refined iRT using refinement model + library iRT
library_irt = getPredIrt(search_context, prec_idx)
refined_irt_pred[i] = if !isnothing(refinement_model) && refinement_model.use_refinement
    refinement_model(precursor_sequence[prec_idx], library_irt)
else
    library_irt
end

# Difference between observed and best refined iRT from other runs
irt_diff[i] = abs(refined_irt_obs[i] - prec_id_to_irt[prec_idx].best_refined_irt)

# MS1-level iRT difference
if !ms1_missing[i]
    ms1_refined_irt_obs = rt_to_refined_irt_interp(ms1_rt[i])
    ms1_irt_diff[i] = abs(ms1_refined_irt_obs - refined_irt_pred[i])
else
    ms1_irt_diff[i] = 0f0
end

# Error between observed and predicted refined iRT
irt_error[i] = abs(refined_irt_obs[i] - refined_irt_pred[i])
```

**Lines 932-933** (DataFrame column assignment):
```julia
# FROM:
psms[!,:irt_obs] = irt_obs
psms[!,:irt_pred] = irt_pred

# TO:
psms[!,:refined_irt_obs] = refined_irt_obs
psms[!,:refined_irt_pred] = refined_irt_pred
```

---

## Summary Table

| File | Lines | Change | Count |
|------|-------|--------|-------|
| SearchTypes.jl | 508+ | Add getRtToRefinedIrtModel/getRefinedIrtToRtModel | +2 functions |
| FirstPassSearch/utils.jl | 98 | Add comment only | 0 code changes |
| SecondPassSearch.jl | 495, 550-553, 570 | getRtIrtModel → getRtToRefinedIrtModel + rt_to_irt_model → rt_to_refined_irt_model | 3 calls + 1 var |
| SecondPassSearch/utils.jl | 181-183, 386-388 | getRtIrtModel → getRtToRefinedIrtModel + irt → refined_irt | 2 calls + 6 var uses |
| SecondPassSearch/utils.jl | 834, 1020 | Rename rt_to_irt_interp → rt_to_refined_irt_interp (params) | 2 params |
| SecondPassSearch/utils.jl | 894+ | Add refinement_model lookup before loop | +1 line |
| SecondPassSearch/utils.jl | 904-913 | Compute refined iRT predictions in loop | ~15 lines |
| SecondPassSearch/utils.jl | 1072 | Rename rt_to_irt_interp → rt_to_refined_irt_interp (usage) | 1 line |
| HuberTuningSearch/utils.jl | 225-227 | getRtIrtModel → getRtToRefinedIrtModel + irt → refined_irt | 1 call + 3 var uses |
| IntegrateChromatogramsSearch/utils.jl | 254-256, 491-493 | getRtIrtModel → getRtToRefinedIrtModel + irt → refined_irt | 2 calls + 6 var uses |

**Total:** 6 files, 2 new functions, 8 function call updates, ~22 variable renames, 1 major loop rewrite

---

## Expected Impact

1. **More accurate iRT-based features** in SecondPassSearch with clear naming:
   - `refined_irt_obs` = observed refined iRT (using rt_to_refined_irt model)
   - `refined_irt_pred` = predicted refined iRT (using refinement model + library iRT)
   - `irt_error` = |refined_irt_obs - refined_irt_pred| in refined iRT space
   - `irt_diff` = |refined_irt_obs - best_refined_irt_from_other_runs|
   - `ms1_irt_diff` = MS1-level refined iRT error

2. **Improved downstream searches** use refined iRT:
   - HuberTuningSearch
   - IntegrateChromatogramsSearch

3. **Clearer code** with explicit getRtToRefinedIrtModel() calls

4. **Potential ID increase** from more accurate iRT features

5. **getPredIrt/setPredIrt remain unchanged** - correctly return library iRT

6. **irt_obs remains unchanged** - correctly stores library iRT

7. **Fallback safety** if refined model unavailable (uses library iRT)
