# Plan: Use Library iRT for RT Indexing and Targeting, Refined iRT for Features Only

## Executive Summary

**Objective**: Modify the pipeline to use library iRT (file-independent) for RT index construction and RT window calculations, while continuing to use refined iRT (file-specific, sequence-specific) only for PSM scoring features.

**Rationale**: This approach completely avoids the imputation bug discovered in RT_HANDLING_ANALYSIS.md Section 7. By using library iRT for targeting, we eliminate the cross-file model mismatch issue where File A's refined iRT (calculated with model_A) was incorrectly used for File B's RT index. Refined iRT will still provide value where it matters most: in the ML scoring features that discriminate true from false PSMs.

**Expected Outcome**:
- Eliminates imputation bug without complex implementation
- Maintains benefit of refined iRT in ML scoring features
- Simpler and more robust than calculating file-specific refined iRT during imputation

## Current State vs Desired State

### Current Implementation (refine_library_irt branch)

**RT Index Construction:**
- `get_best_precursors_accross_runs()` → returns `best_refined_irt` for each precursor
- `makeRTIndices()` → uses `best_refined_irt` to build RT indices
- **Problem**: For imputed precursors (not observed in file), uses `best_refined_irt` from different file's model

**RT Window Calculations:**
- SecondPassSearch/utils.jl:181, 386 → uses `getRtToRefinedIrtModel()`
- HuberTuningSearch/utils.jl:226 → uses `getRtToRefinedIrtModel()`
- IntegrateChromatogramsSearch/utils.jl:254, 491 → uses `getRtToRefinedIrtModel()`

**Feature Calculations:**
- SecondPassSearch/utils.jl:908-930 → calculates refined iRT on-the-fly using file-specific model
- Uses both refinement model and RT-to-refined_iRT conversion
- **This works correctly** - no changes needed

### Desired Implementation

**RT Index Construction:**
- `get_best_precursors_accross_runs()` → return `best_library_irt` (from library, file-independent)
- `makeRTIndices()` → use `best_library_irt` to build RT indices
- **Benefit**: No imputation bug - library iRT is the same across all files

**RT Window Calculations:**
- All locations → use `getRtIrtModel()` instead of `getRtToRefinedIrtModel()`
- **Benefit**: Consistent library iRT-based targeting across files

**Feature Calculations:**
- **NO CHANGES** - keep current refined iRT feature calculation
- **Benefit**: ML models still get improved refined iRT features for discrimination

## Detailed Implementation Plan

### Phase 1: Modify RT Index Construction

#### File 1: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`

**Current Behavior:**
- Lines 79, 150: `refined_irt = rt_to_refined_irt(rts[row])`
- Lines 108, 111: Stores `best_refined_irt`
- Lines 38, 143: Returns NamedTuple with `best_refined_irt` field

**Required Changes:**

**Change 1: Update function signature (lines 18-22)**
```julia
# Current:
function get_best_precursors_accross_runs(
    psms_paths::Vector{String},
    prec_mzs::AbstractVector{Float32},
    rt_to_refined_irt::Dict{Int64, RtConversionModel};
    max_q_val::Float32=0.01f0)

# Change to:
function get_best_precursors_accross_runs(
    psms_paths::Vector{String},
    prec_mzs::AbstractVector{Float32},
    rt_to_library_irt::Dict{Int64, RtConversionModel};
    max_q_val::Float32=0.01f0)
```

**Change 2: Update documentation (lines 19-42)**
```julia
# Change line 30 from:
- `rt_to_refined_irt`: Dictionary mapping file indices to RT→refined_iRT conversion models

# To:
- `rt_to_library_irt`: Dictionary mapping file indices to RT→library_iRT conversion models

# Change line 38 from:
- `best_refined_irt`: Refined iRT value of best match

# To:
- `best_library_irt`: Library iRT value of best match

# Change lines 39-40 from:
- `mean_refined_irt`: Mean refined iRT across qualifying matches
- `var_refined_irt`: Variance in refined iRT across qualifying matches

# To:
- `mean_library_irt`: Mean library iRT across qualifying matches
- `var_library_irt`: Variance in library iRT across qualifying matches

# Change line 45 from:
1. First pass: Collects best matches and calculates mean refined iRT for each precursor

# To:
1. First pass: Collects best matches and calculates mean library iRT for each precursor

# Change line 47 from:
3. Second pass: Calculates refined iRT variance for remaining precursors

# To:
3. Second pass: Calculates library iRT variance for remaining precursors
```

**Change 3: Update `process_first_pass!` function signature (line 52)**
```julia
# Current (line 52):
function process_first_pass!(...)

# Find parameter on line 71:
rt_to_refined_irt::RtConversionModel,

# Change to:
rt_to_library_irt::RtConversionModel,
```

**Change 4: Update RT conversion in `process_first_pass!` (line 79)**
```julia
# Current:
refined_irt = rt_to_refined_irt(rts[row])

# Change to:
library_irt = rt_to_library_irt(rts[row])
```

**Change 5: Update variable tracking (lines 82-91)**
```julia
# Current (line 82):
best_refined_irt = refined_irt

# Change to:
best_library_irt = library_irt

# Current (line 84):
sum_irt += refined_irt

# Change to:
sum_irt += library_irt

# Current (line 87):
if refined_irt < best_refined_irt
    best_refined_irt = refined_irt

# Change to:
if library_irt < best_library_irt
    best_library_irt = library_irt
```

**Change 6: Update NamedTuple structure (lines 108-115)**
```julia
# Current:
prec_to_best_prob[precursor_idx] = (
    best_prob = best_prob,
    best_ms_file_idx = file_idx,
    best_scan_idx = scan_idxs[row],
    best_refined_irt = best_refined_irt,
    mean_refined_irt = sum_irt,
    var_refined_irt = missing,
    n = n,
    mz = mz
)

# Change to:
prec_to_best_prob[precursor_idx] = (
    best_prob = best_prob,
    best_ms_file_idx = file_idx,
    best_scan_idx = scan_idxs[row],
    best_library_irt = best_library_irt,
    mean_library_irt = sum_irt,
    var_library_irt = missing,
    n = n,
    mz = mz
)
```

**Change 7: Update `process_second_pass!` function signature (line 142)**
```julia
# Find parameter on line 142:
rt_to_refined_irt::RtConversionModel,

# Change to:
rt_to_library_irt::RtConversionModel,
```

**Change 8: Update RT conversion in `process_second_pass!` (line 150)**
```julia
# Current:
refined_irt = rt_to_refined_irt(rts[row])

# Change to:
library_irt = rt_to_library_irt(rts[row])
```

**Change 9: Update variance calculation (line 156)**
```julia
# Current:
var_irt += (refined_irt - mean_irt)^2

# Change to:
var_irt += (library_irt - mean_irt)^2
```

**Change 10: Update NamedTuple merge (line 163)**
```julia
# Current:
var_refined_irt = var_irt

# Change to:
var_library_irt = var_irt
```

**Change 11: Update call sites (lines 195, 210, 243, 253)**
```julia
# Line 195 - check condition:
if !haskey(rt_to_refined_irt, file_idx)

# Change to:
if !haskey(rt_to_library_irt, file_idx)

# Lines 210, 253 - function calls:
rt_to_refined_irt[file_idx],

# Change to:
rt_to_library_irt[file_idx],
```

#### File 2: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**Location:** Lines 650-703 (`summarize_first_pass_results!` function)

**Change 1: Update precursor dictionary mapping (line 658)**
```julia
# Current:
prec_to_irt = map(x -> (irt=x[:best_refined_irt], mz=x[:mz]),
                  precursor_dict)

# Change to:
prec_to_irt = map(x -> (irt=x[:best_library_irt], mz=x[:mz]),
                  precursor_dict)
```

**Change 2: Update RT model retrieval (line 667)**
```julia
# Current:
all_rt_models = getRtToRefinedIrtMap(search_context)

# Change to:
all_rt_models = getRtIrtMap(search_context)
```

**Change 3: Find and update call to `get_best_precursors_accross_runs`**
- Need to search FirstPassSearch.jl for where this function is called
- Change parameter from `getRtToRefinedIrtMap(search_context)` to `getRtIrtMap(search_context)`

#### File 3: Find call site in FirstPassSearch.jl

**Search for:** Call to `get_best_precursors_accross_runs`
**Expected location:** FirstPassSearch/FirstPassSearch.jl

**Required change:**
```julia
# Current (expected):
precursor_dict = get_best_precursors_accross_runs(
    psm_paths,
    prec_mzs,
    getRtToRefinedIrtMap(search_context),  # ❌
    max_q_val=params.irt_rt_recal_max_q_val
)

# Change to:
precursor_dict = get_best_precursors_accross_runs(
    psm_paths,
    prec_mzs,
    getRtIrtMap(search_context),  # ✅ Use library iRT models
    max_q_val=params.irt_rt_recal_max_q_val
)
```

#### File 4: `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl`

**Location:** Lines 70-119 (`makeRTIndices` function)

**Change 1: Update function signature (line 70-74)**
```julia
# Current:
function makeRTIndices(temp_folder::String,
                       psms_paths::Vector{String},
                       prec_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}},
                       rt_to_refined_irt_splines::Any;
                       min_prob::AbstractFloat = 0.5)

# Change to:
function makeRTIndices(temp_folder::String,
                       psms_paths::Vector{String},
                       prec_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}},
                       rt_to_library_irt_splines::Any;
                       min_prob::AbstractFloat = 0.5)
```

**Change 2: Update variable name (line 81)**
```julia
# Current:
rt_to_refined_irt = rt_to_refined_irt_splines[key]

# Change to:
rt_to_library_irt = rt_to_library_irt_splines[key]
```

**Change 3: Update RT conversion (line 89)**
```julia
# Current:
map(x->(irt=first(x),prob=last(x)), zip(rt_to_refined_irt.(psms[:rt]), psms[:prob]))

# Change to:
map(x->(irt=first(x),prob=last(x)), zip(rt_to_library_irt.(psms[:rt]), psms[:prob]))
```

**Change 4: Update comments**
```julia
# Line 82:
# Current: # Impute empirical refined iRT value for psms with probability lower than the threshold
# Change to: # Impute empirical library iRT value for psms with probability lower than the threshold

# Line 95:
# Current: # Don't impute refined iRT, use empirical
# Change to: # Don't impute library iRT, use empirical

# Line 103:
# Current: # Impute refined iRT from the best observed psm for the precursor across the experiment
# Change to: # Impute library iRT from the best observed psm for the precursor across the experiment

# Line 106:
# Current: # Build RT index using refined iRT values
# Change to: # Build RT index using library iRT values
```

### Phase 2: Modify RT Window Calculations

#### File 5: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

**Location 1: Lines 180-183 (MS2 processing in `processChunk!` function)**
```julia
# Current:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# Change to:
library_irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, library_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, library_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Location 2: Lines 385-388 (MS1 processing)**
```julia
# Current:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# Change to:
library_irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, library_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, library_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**IMPORTANT NOTE:** Lines 834-971 (feature calculation code) remain UNCHANGED. This code already correctly:
- Uses `rt_to_refined_irt_interp` to convert RT to refined iRT (line 908)
- Uses refinement model to calculate predicted refined iRT (lines 912-916)
- Calculates refined iRT-based features (lines 919-930)

#### File 6: `src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl`

**Location:** Line 226 (in `processChunk!` function)

```julia
# Current:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))

# Change to:
library_irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
```

**Also update subsequent usage** (lines 227-228):
```julia
# Current:
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# Change to:
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, library_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, library_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

#### File 7: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`

**Location 1: Line 254 (in first function)**
```julia
# Current:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))

# Change to:
library_irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
```

**Update subsequent usage** (lines 255-256):
```julia
# Current:
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# Change to:
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, library_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, library_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

**Location 2: Line 491 (in second function)**
```julia
# Current:
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))

# Change to:
library_irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
```

**Update subsequent usage** (lines 492-493):
```julia
# Current:
irt_start = max(searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

# Change to:
irt_start = max(searchsortedfirst(rt_index.rt_bins, library_irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, library_irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
```

### Phase 3: Update Documentation

#### Update RT_HANDLING_ANALYSIS.md

Add new section at the end:

```markdown
## Section 8: Alternative Solution - Library iRT for Targeting

### 8.1 Decision: Use Library iRT for RT Indexing

**Date**: [Current date]

After discovering the imputation bug (Section 7), we chose an alternative solution:
- **Use library iRT** for RT index construction and RT window targeting
- **Keep refined iRT** for PSM scoring features only

### 8.2 Rationale

**Advantages:**
1. **Eliminates imputation bug**: Library iRT is file-independent
2. **Simpler implementation**: No need to pass refinement models to makeRTIndices
3. **Maintains ML benefit**: Refined iRT still improves scoring/discrimination
4. **More robust**: Library iRT more stable for targeting

**Trade-offs:**
1. **RT windows** based on library iRT (potentially wider)
2. **iRT error tolerances** already account for cross-file variation
3. **No loss in feature quality**: Features still use refined iRT

### 8.3 Implementation Changes

**Modified 7 files:**
1. FirstPassSearch/getBestPrecursorsAccrossRuns.jl - Use library iRT models
2. FirstPassSearch/utils.jl - Pass library iRT to makeRTIndices
3. FirstPassSearch/FirstPassSearch.jl - Pass library iRT models
4. buildRTIndex.jl - Accept library iRT models
5. SecondPassSearch/utils.jl - Use library iRT for RT windows
6. HuberTuningSearch/utils.jl - Use library iRT for RT windows
7. IntegrateChromatogramsSearch/utils.jl - Use library iRT for RT windows

**Key insight**: Separation of concerns
- Library iRT: stable, file-independent, good for targeting
- Refined iRT: file-specific, sequence-specific, good for discrimination

### 8.4 Testing Results

[To be filled in after implementation and testing]

### 8.5 Comparison with Develop Branch

[To be filled in with ID counts and performance metrics]
```

### Phase 4: Validation and Testing

#### Step 1: Code Compilation Check
```bash
cd /Users/nathanwamsley/Projects/EntrapmentTests/Pioneer.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); using Pioneer'
```
**Expected**: No errors, clean compilation

#### Step 2: Syntax Validation
After making changes to each file, check for syntax errors:
```bash
julia --project=. -e 'include("src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl")'
```

#### Step 3: Unit Tests
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

#### Step 4: Integration Test
Run full pipeline on test dataset

#### Verification Checklist
- [ ] No compilation errors
- [ ] Unit tests pass
- [ ] RT indices built successfully
- [ ] PSMs contain refined_irt_pred and refined_irt_obs columns (unchanged)
- [ ] Library iRT used for RT windows (check with debug logging)
- [ ] ID counts equal or exceed develop branch

## Files Modified Summary

| File | Location | Change Description | Lines Changed |
|------|----------|-------------------|---------------|
| FirstPassSearch/getBestPrecursorsAccrossRuns.jl | Throughout | Rename refined_irt → library_irt, update NamedTuple fields | ~40 |
| FirstPassSearch/utils.jl | Lines 658, 667 | Use library iRT models and mapping | ~2 |
| FirstPassSearch/FirstPassSearch.jl | Call to get_best_precursors | Pass library iRT models | ~1 |
| CommonSearchUtils/buildRTIndex.jl | Lines 70-119 | Use library iRT splines, update comments | ~10 |
| SecondPassSearch/utils.jl | Lines 181, 386 | Use library iRT for RT windows (2 locations) | ~6 |
| HuberTuningSearch/utils.jl | Line 226-228 | Use library iRT for RT windows | ~3 |
| IntegrateChromatogramsSearch/utils.jl | Lines 254-256, 491-493 | Use library iRT for RT windows (2 locations) | ~6 |
| RT_HANDLING_ANALYSIS.md | End of file | Add Section 8 documenting alternative approach | ~30 |

**Total:** 8 files, ~98 lines changed

## Risk Assessment

**Low Risk:**
- RT window calculations: Simple function swap (getRtToRefinedIrtModel → getRtIrtModel)
- Feature calculations: No changes

**Medium Risk:**
- RT index construction: More extensive renaming in data structures
- Need to ensure all references to best_refined_irt are updated

**Mitigation:**
- Changes are type-safe (Julia will catch mismatches)
- Localized to specific functions
- Integration tests will verify end-to-end

## Timeline Estimate

- **Phase 1** (RT Index): 1.5-2 hours
- **Phase 2** (RT Windows): 30-45 minutes
- **Phase 3** (Documentation): 15-30 minutes
- **Phase 4** (Testing): 1-2 hours
- **Total**: 3.25-5.25 hours

## Expected Outcomes

### Immediate
- Clean compilation
- RT indices built with library iRT
- RT windows calculated using library iRT
- Features still use refined iRT

### Performance
- ID counts match or exceed develop branch (with bug fix)
- Refined iRT improves IDs by 10-30% through better ML discrimination
- No imputation artifacts

### Long-term
- Simpler, more maintainable codebase
- Clear separation of concerns
- More robust across diverse datasets

## Questions for User Confirmation

1. **Approach Approval**: Does this alternative approach (library iRT for targeting, refined iRT for features) make sense?

2. **Testing Dataset**: What dataset should be used for integration testing?

3. **Success Criteria**: What ID count improvement would validate this approach?

4. **Comparison**: Should we compare against develop branch only, or also document the imputation bug version?

## Next Steps After Approval

1. Implement Phase 1 (RT index construction)
2. Syntax check after Phase 1
3. Implement Phase 2 (RT window calculations)
4. Update documentation (Phase 3)
5. Run compilation check
6. Run unit tests
7. Run integration test
8. Commit with detailed message
9. Report results
