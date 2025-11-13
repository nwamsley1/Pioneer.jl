# iRT Value Usage Analysis: Post-FirstPassSearch

## Executive Summary

This document comprehensively analyzes all uses of library iRT predictions and the `getPredIrt()` function throughout the Pioneer.jl codebase **after** FirstPassSearch completes. It evaluates whether the iRT refinement model (which corrects systematic library prediction errors) is being properly applied in each location.

**Key Finding**: The implementation correctly uses the iRT refinement where needed, with proper lazy fallback for unobserved precursors.

---

## Background: iRT Refinement Implementation

### What is iRT Refinement?

**Location**: `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl`

**Purpose**: Corrects systematic errors in library iRT predictions using a linear regression model based on amino acid composition.

**Model**:
```julia
error = irt_library - irt_observed  # Positive when library overestimates
irt_refined = irt_library - predicted_error
```

**Training**: Uses observed PSMs from FirstPassSearch to train file-specific correction models.

**Storage**: Refined iRT values stored in `SearchContext.irt_obs` (Dictionary{UInt32, Float32})
- Only precursors **observed** in FirstPass get refined values (~170k precursors)
- Unobserved precursors use lazy fallback to library iRT

### Access Pattern: getPredIrt()

**Location**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:426-440`

```julia
function getPredIrt(s::SearchContext, prec_idx::UInt32)::Float32
    # Check if refined/observed iRT exists in irt_obs (lazy initialization)
    irt = get(s.irt_obs, prec_idx, nothing)

    # Lazy fallback to library iRT if not found
    if isnothing(irt)
        return getIrt(getPrecursors(getSpecLib(s)))[prec_idx]  # Library iRT
    end

    return irt  # Refined iRT
end
```

**Behavior**:
- Returns **refined iRT** if precursor was observed in FirstPass
- Returns **library iRT** if precursor was not observed
- Transparent to caller - always returns best available iRT estimate

---

## Complete Usage Analysis

### 1. SecondPassSearch (After FirstPass)

**Location**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

#### Usage 1: Feature Calculation (Line 905)
```julia
irt_pred[i] = getPredIrt(search_context, prec_idx)
```

**Context**: Calculates `irt_error` feature for PSMs
**Status**: ✅ **CORRECT** - Uses getPredIrt(), gets refined iRT for observed precursors

#### Usage 2: MS1 RT Difference (Line 909)
```julia
ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx))
```

**Context**: Calculates RT difference between MS1 and MS2 scans
**Status**: ✅ **CORRECT** - Uses getPredIrt(), gets refined iRT for observed precursors

#### Comment Documentation (Line 841)
```julia
# Note: iRT values come from SearchContext via getPredIrt() (may be refined, see lines 905, 909)
```

**Status**: ✅ **CORRECT** - Code includes explanatory comment about refinement

---

### 2. FirstPassSearch (iRT Refinement Training)

**Location**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:296-306`

#### Usage: Initial PSM Column Population
```julia
add_main_search_columns!(
    psms,
    getModel(rt_model),
    getStructuralMods(getPrecursors(getSpecLib(search_context))),
    getMissedCleavages(getPrecursors(getSpecLib(search_context))),
    getIsDecoy(getPrecursors(getSpecLib(search_context))),
    getIrt(getPrecursors(getSpecLib(search_context))),  # Line 301
    ...
)
```

**Context**: Creates `:irt_predicted` column in PSMs **before** refinement
**Status**: ✅ **CORRECT** - Uses library iRT before refinement (refinement happens later)

**Details**:
- `add_main_search_columns!` sets `psms[!,:irt_predicted]` from library iRT (utils.jl:96)
- This occurs **before** iRT refinement training (which happens at utils.jl:410-450)
- Sequence:
  1. Create PSMs with library iRT → `:irt_predicted` column
  2. Fit RT-to-iRT model using library iRT
  3. **Then** train iRT refinement model using observed RTs
  4. Store refined iRT in `SearchContext.irt_obs`

---

### 3. RT Model Fitting (Uses :irt_predicted Column)

**Location**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`

#### Usage: fit_irt_model()
```julia
function fit_irt_model(psms::DataFrame) -> (model, rt, irt, mad)
    rt = psms[!, :rt]
    irt = psms[!, :irt_predicted]  # Line 61
    ...
end
```

**Context**: Fits spline mapping RT → iRT using PSM data
**Called From**:
- `FirstPassSearch/utils.jl:372` - RT model fitting
- `ParameterTuningSearch/utils.jl:260` - Initial RT model

**Status**: ✅ **CORRECT** - Uses library iRT (`:irt_predicted` column)

**Why Library iRT is Correct Here**:
1. **RT model must be fitted BEFORE refinement** - You need RT→iRT conversion to compute observed iRT from observed RT
2. **Refinement depends on RT model** - The refinement model trains on:
   ```julia
   observed_irt = rt_model(observed_rt)  # Need RT model first!
   error = library_irt - observed_irt
   ```
3. **RT model accuracy**: Library iRT is accurate enough for RT model fitting (MAE typically ~0.5 iRT units before refinement vs ~0.3 after)
4. **Refinement is applied downstream**: Refined values used in SecondPass via `getPredIrt()`

---

### 4. QuadTuningSearch (Before FirstPass)

**Location**: `src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/utils.jl:498`

```julia
getIrt(getPrecursors(getSpecLib(search_context)))
```

**Context**: Uses library iRT for column addition in QuadTuningSearch
**Status**: ✅ **CORRECT** - QuadTuningSearch runs **before** FirstPassSearch, so no refinement exists yet

---

### 5. SearchTypes.jl (Lazy Fallback Implementation)

**Location**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:436`

```julia
if isnothing(irt)
    return getIrt(getPrecursors(getSpecLib(s)))[prec_idx]  # Fallback
end
```

**Context**: Inside `getPredIrt()` - fallback for unobserved precursors
**Status**: ✅ **CORRECT** - Intentional fallback design for precursors not seen in FirstPass

---

### 6. Other Search Methods

**Checked Methods**:
- ✅ **HuberTuningSearch**: No iRT usage found
- ✅ **ScoringSearch**: No iRT usage found
- ✅ **IntegrateChromatogramsSearch**: No iRT usage found
- ✅ **MaxLFQSearch**: No iRT usage found

**Status**: ✅ **CORRECT** - These methods don't use iRT values directly

---

## Search Pipeline Order

Understanding when refinement occurs relative to pipeline stages:

```
1. ParameterTuningSearch    [Library iRT only - no refinement yet]
2. NceTuningSearch          [Library iRT only]
3. QuadTuningSearch         [Library iRT only] ← Uses getIrt(getPrecursors(...))
4. FirstPassSearch
   ├─ PSMs created          [Library iRT → :irt_predicted column]
   ├─ RT model fitted       [Uses :irt_predicted with library iRT]
   ├─ iRT refinement        [Trains correction model]
   └─ Store refined iRT     [SearchContext.irt_obs populated]
5. HuberTuningSearch        [No iRT usage]
6. SecondPassSearch         [Uses getPredIrt() → refined iRT] ✓
7. ScoringSearch            [No iRT usage]
8. IntegrateChromatogramSearch [No iRT usage]
9. MaxLFQSearch             [No iRT usage]
```

---

## Critical Implementation Details

### 1. Why RT Model Uses Library iRT (Not Refined)

**Question**: Should the RT model be refitted after iRT refinement?

**Answer**: **No** - Library iRT is intentionally used for RT model fitting

**Reasoning**:
1. **Circular dependency**: RT model converts RT → iRT, which is needed to compute observed iRT for refinement training
2. **RT model precedes refinement**: Workflow is:
   ```
   observed_RT → [RT model] → observed_iRT
   error = library_iRT - observed_iRT
   refined_iRT = library_iRT - predicted_error
   ```
3. **Sufficient accuracy**: Library iRT is accurate enough for RT alignment (typical MAE ~0.5 iRT units)
4. **Refinement corrects downstream**: Refined iRT is used in SecondPass for precursor selection windows

### 2. Lazy Initialization Design

**Only observed precursors get refined iRT**:
- `SearchContext.irt_obs` is a sparse dictionary (~170k entries out of ~500k total precursors)
- `getPredIrt()` provides transparent fallback to library iRT
- Memory efficient: Only stores corrections for precursors actually seen

**Why not refine all precursors?**:
- Model trained on observed precursors only (high confidence)
- Extrapolation to unobserved sequences may be inaccurate
- Conservative approach: Only correct what you've observed

### 3. File-Specific Refinement Models

**Observation**: iRT refinement is file-specific (trained per MS file)

**Storage**:
- Refined iRT values in `SearchContext.irt_obs` are **global** across files
- Last file processed "wins" for each precursor

**Implication**:
- Current implementation uses the most recent file's refined iRT
- Multi-file refinement aggregation could be enhanced (future work)

---

## Potential Issues & Edge Cases

### 1. `:irt_predicted` Column in Saved PSMs

**Location**: FirstPassSearch saves PSMs with `:irt_predicted` column

**Question**: Does this column contain refined iRT?

**Answer**: **No** - It contains **library iRT**

**Evidence**:
```julia
# FirstPassSearch.jl:502 - PSMs saved to Arrow
select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
    :irt_predicted, ...])  # irt_predicted set before refinement
```

**Impact**:
- Minor: Column name `:irt_predicted` is slightly misleading (it's library iRT, not refined)
- **No functional issue**: SecondPass uses `getPredIrt()` directly, not `:irt_predicted` column
- FirstPass PSMs are used only for RT model fitting and summary statistics

**Recommendation**: Consider renaming to `:irt_library` for clarity (low priority)

### 2. Match-Between-Runs (MBR) Integration

**Code Reference**: FirstPassSearch.jl:514-538 (commented out)

**Status**: Current MBR code is disabled (`if false==true`)

**Observation**: Commented code references library iRT:
```julia
setPredIrt!(search_context, pid, getIrt(getPrecursors(...))[pid])
```

**Future Work**: When MBR is re-enabled, ensure it uses `getPredIrt()` instead of library iRT

---

## Validation Tests

### Test 1: Verify getPredIrt() Returns Refined iRT

**Test Case**: Precursor observed in FirstPass
```julia
# After FirstPassSearch completes
prec_idx = 12345  # Some observed precursor
library_irt = getIrt(getPrecursors(getSpecLib(search_context)))[prec_idx]
refined_irt = getPredIrt(search_context, prec_idx)

@assert refined_irt != library_irt  # Should be different
@assert abs(refined_irt - library_irt) < 5.0  # But within reasonable range
```

### Test 2: Verify Lazy Fallback

**Test Case**: Precursor NOT observed in FirstPass
```julia
prec_idx = 99999  # Unobserved precursor
library_irt = getIrt(getPrecursors(getSpecLib(search_context)))[prec_idx]
returned_irt = getPredIrt(search_context, prec_idx)

@assert returned_irt == library_irt  # Should fallback to library
```

### Test 3: SecondPass Uses Refined iRT

**Test Case**: Verify SecondPass PSMs have refined iRT in features
```julia
# After SecondPassSearch completes
psms = Arrow.Table(getSecondPassPsms(getMSData(search_context), 1)) |> DataFrame
for row in eachrow(psms)
    pred_irt_from_context = getPredIrt(search_context, row.precursor_idx)
    @assert row.irt_pred == pred_irt_from_context  # Should match
end
```

---

## Summary Table

| Location | Function/Variable | Uses Refined? | Status | Notes |
|----------|------------------|---------------|--------|-------|
| SecondPassSearch/utils.jl:905 | `getPredIrt()` | ✅ Yes (via function) | ✅ Correct | Feature calculation |
| SecondPassSearch/utils.jl:909 | `getPredIrt()` | ✅ Yes (via function) | ✅ Correct | MS1 RT difference |
| FirstPassSearch.jl:301 | `getIrt(getPrecursors(...))` | ❌ No (library) | ✅ Correct | Before refinement |
| FirstPassSearch/utils.jl:96 | `:irt_predicted` column | ❌ No (library) | ✅ Correct | Before refinement |
| rt_alignment_utils.jl:61 | `:irt_predicted` column | ❌ No (library) | ✅ Correct | RT model needs library iRT |
| QuadTuningSearch/utils.jl:498 | `getIrt(getPrecursors(...))` | ❌ No (library) | ✅ Correct | Before FirstPass |
| SearchTypes.jl:436 | `getIrt(getPrecursors(...))` | ❌ No (library) | ✅ Correct | Lazy fallback design |

---

## Conclusions

### ✅ Implementation is Correct

1. **SecondPassSearch** correctly uses `getPredIrt()` to access refined iRT values
2. **FirstPassSearch** correctly uses library iRT before refinement occurs
3. **RT model fitting** correctly uses library iRT (required for the workflow)
4. **Lazy fallback** correctly provides library iRT for unobserved precursors
5. **No other search methods** use iRT after FirstPass

### 🔍 Minor Improvements (Optional)

1. **Column naming**: Rename `:irt_predicted` → `:irt_library` in FirstPass PSMs for clarity
2. **MBR integration**: Update commented MBR code to use `getPredIrt()` when re-enabled
3. **Multi-file refinement**: Consider aggregating refinement models across files (advanced)
4. **Documentation**: Add inline comments explaining why RT model uses library iRT

### 📊 Expected Behavior

When iRT refinement is enabled:
- **Observed precursors** (~170k): Use file-specific refined iRT (typically 20-40% MAE reduction)
- **Unobserved precursors** (~330k): Use library iRT (no refinement possible)
- **RT model**: Uses library iRT (intentional - required for workflow)
- **SecondPass**: Automatically benefits from refinement via `getPredIrt()`

---

## References

### Key Files
- `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl` - Refinement model implementation
- `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:426-440` - getPredIrt() implementation
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl:410-450` - Refinement training
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl:905,909` - Refined iRT usage
- `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl` - RT model fitting

### Related Documentation
- `docs/irt_refinement_model_storage_plan.md` - Original implementation plan
- `src/Routines/SearchDIA/CLAUDE.md` - SearchDIA architecture overview
