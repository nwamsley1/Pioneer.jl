# iRT Refinement Implementation Plan

## Overview
Implement file-specific iRT refinement using amino acid composition to improve retention time prediction accuracy. After FirstPassSearch, train a linear model to predict the error between library iRT and observed iRT based on amino acid counts and library iRT. Use this model to refine iRT values for all subsequent searches.

---

## Key Concepts

### Terminology
- **library_irt** (`:irt_predicted`): Original iRT from spectral library
- **observed_irt**: RT converted to iRT via `rt_to_irt` spline
- **refined_irt** (`:refined_irt`): Library iRT corrected using refinement model (`IrtRefinementModel(sequence, library_irt)`)
- **observed_refined_irt**: RT converted via `rt_to_refined_irt` spline

### Data Flow
```
FirstPassSearch Workflow:
  1. Fit rt_to_irt spline (RT → library_irt)
  2. Calculate observed_irt = rt_to_irt(rt) for top PSMs
  3. Train IrtRefinementModel: irt_error = f(library_irt, AA_counts)
  4. Apply model to get refined_irt for each PSM
  5. Fit rt_to_refined_irt spline (RT → refined_irt)
  6. Overwrite rt_irt_map with rt_to_refined_irt spline
  7. Add :refined_irt column to PSM files
  8. Build RT indices using refined_irt

SecondPassSearch+ (all subsequent searches):
  - Use rt_to_refined_irt (stored in rt_irt_map) to convert scan RT → observed_refined_irt
  - Apply IrtRefinementModel to library_irt → refined_irt for search windows
  - Calculate irt_error = |observed_refined_irt - refined_irt|
```

---

## Implementation Steps

### 1. Create IrtRefinementModel Struct

**File:** `src/structs/RetentionTimeConversionModel.jl`

**Location:** After line 58 (after `RtConversionModel() = IdentityModel()`)

**ADD THIS CODE:**

```julia
"""
    IrtRefinementModel

File-specific model to refine library iRT predictions using amino acid composition.

# Fields
- `use_refinement::Bool`: Whether refinement improves validation MAE
- `aa_coefficients::Dict{Char, Float32}`: Per-AA weights (20 standard AAs)
- `intercept::Float32`: Model intercept
- `irt_coefficient::Float32`: Weight for library_irt feature
- `mae_original::Float32`: Validation MAE without refinement
- `mae_refined::Float32`: Validation MAE with refinement
- `r2_train::Float32`: Training R²
- `r2_val::Float32`: Validation R²

# Callable Interface
Model is callable: `refined_irt = model(sequence::String, library_irt::Float32)`

# Algorithm
Predicts error = library_irt - observed_irt, then:
refined_irt = library_irt - predicted_error

# Example
```julia
model = IrtRefinementModel(true, aa_weights, 0.5f0, 0.1f0, ...)
refined = model("PEPTIDE", 50.0f0)  # Returns refined iRT
```
"""
struct IrtRefinementModel
    use_refinement::Bool
    aa_coefficients::Dict{Char, Float32}
    intercept::Float32
    irt_coefficient::Float32
    mae_original::Float32
    mae_refined::Float32
    r2_train::Float32
    r2_val::Float32
end

"""
    (model::IrtRefinementModel)(sequence::String, library_irt::Float32) -> Float32

Apply iRT refinement to a sequence. Zero-allocation via Dict lookup.
"""
function (model::IrtRefinementModel)(sequence::String, library_irt::Float32)::Float32
    if !model.use_refinement
        return library_irt
    end

    # Calculate predicted error
    error_pred = model.intercept + model.irt_coefficient * library_irt

    # Add AA contributions
    for aa in sequence
        if haskey(model.aa_coefficients, aa)
            error_pred += model.aa_coefficients[aa]
        end
    end

    # Return refined iRT (subtract predicted error)
    return library_irt - error_pred
end

# 20 standard amino acids for model training
const STANDARD_AAS = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
```

---

### 2. Update SearchContext

**File:** `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

#### Change 2a: Add field to struct (line 223)

**CURRENT CODE:**
```julia
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    precursor_dict::Base.Ref{Dictionary}
```

**CHANGED CODE:**
```julia
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}  # NEW FIELD
    precursor_dict::Base.Ref{Dictionary}
```

#### Change 2b: Initialize in constructor (line 267)

**CURRENT CODE:**
```julia
            Dict{Int64, RtConversionModel}(),  # line 267: irt_rt_map
            Dict{Int64, RtConversionModel}(),  # line 268: rt_irt_map
            Ref{Dictionary}(),                 # line 269: precursor_dict
```

**CHANGED CODE:**
```julia
            Dict{Int64, RtConversionModel}(),                          # line 267: irt_rt_map
            Dict{Int64, RtConversionModel}(),                          # line 268: rt_irt_map
            Dict{Int64, Union{IrtRefinementModel, Nothing}}(),         # NEW LINE: irt_refinement_models
            Ref{Dictionary}(),                                         # line 269: precursor_dict (now line 270)
```

#### Change 2c: Add getters/setters (after line 524)

**CURRENT CODE:**
```julia
function setRtIrtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.rt_irt_map[index] = rcm
end
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
```

**ADD AFTER setRtIrtMap!:**
```julia
"""
    getIrtRefinementModel(s::SearchContext, index::Integer)

Get iRT refinement model for MS file index. Returns nothing if not found.
"""
function getIrtRefinementModel(s::SearchContext, index::I) where {I<:Integer}
    return get(s.irt_refinement_models, index, nothing)
end

"""
    setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::Integer)

Store iRT refinement model for MS file index.
"""
function setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::I) where {I<:Integer}
    s.irt_refinement_models[index] = model
end
```

#### Change 2d: Update setPrecursorDict! signature (line 525)

**CURRENT CODE:**
```julia
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
```

**CHANGED CODE:**
```julia
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_refined_irt::Float32, mean_refined_irt::Union{Missing, Float32}, var_refined_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
```

**Note:** Field names changed from `:best_irt`, `:mean_irt`, `:var_irt` to `:best_refined_irt`, `:mean_refined_irt`, `:var_refined_irt`

---

### 3. Create Helper Functions for iRT Refinement

**File:** `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl` (NEW FILE)

**CREATE THIS FILE with full contents from plan (lines 158-434 in original plan)**

[Contents remain the same as in original plan - the helper functions prepare_features_dataframe, fit_irt_refinement_model, and add_refined_irt_column!]

---

### 4. Modify map_retention_times!() in FirstPassSearch

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**REPLACE entire function (lines 345-473) with updated version that includes 5-step workflow**

[Contents remain the same as in original plan - lines 445-640]

---

### 5. Update Precursor Dictionary Field Names

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**CURRENT CODE (lines 530-536):**
```julia
PrecToIrtType = Dictionary{UInt32,
    NamedTuple{
        (:best_prob, :best_ms_file_idx, :best_scan_idx, :best_irt, :mean_irt, :var_irt, :n, :mz),
        Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}
    }
}
```

**CHANGED CODE:**
```julia
PrecToIrtType = Dictionary{UInt32,
    NamedTuple{
        (:best_prob, :best_ms_file_idx, :best_scan_idx, :best_refined_irt, :mean_refined_irt, :var_refined_irt, :n, :mz),
        Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}
    }
}
```

**Changes:**
- `:best_irt` → `:best_refined_irt`
- `:mean_irt` → `:mean_refined_irt`
- `:var_irt` → `:var_refined_irt`

---

### 6. Update get_best_precursors_accross_runs

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`

#### Change 6a: Update function documentation (lines 34-42)

**CURRENT CODE:**
```julia
- `best_prob`: Highest probability score
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_irt`: iRT value of best match
- `mean_irt`: Mean iRT across qualifying matches
- `var_irt`: Variance in iRT across qualifying matches
- `n`: Number of qualifying matches
- `mz`: Precursor m/z value
```

**CHANGED CODE:**
```julia
- `best_prob`: Highest probability score
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_refined_irt`: refined iRT value of best match
- `mean_refined_irt`: Mean refined iRT across qualifying matches
- `var_refined_irt`: Variance in refined iRT across qualifying matches
- `n`: Number of qualifying matches
- `mz`: Precursor m/z value
```

#### Change 6b: Update readPSMs! function signature (lines 57-64)

**CURRENT CODE:**
```julia
    function readPSMs!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_irt::Float32,
                                                    mean_irt::Union{Missing, Float32},
                                                    var_irt::Union{Missing, Float32},
                                                    n::Union{Missing, UInt16},
                                                    mz::Float32}},
```

**CHANGED CODE:**
```julia
    function readPSMs!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_refined_irt::Float32,
                                                    mean_refined_irt::Union{Missing, Float32},
                                                    var_refined_irt::Union{Missing, Float32},
                                                    n::Union{Missing, UInt16},
                                                    mz::Float32}},
```

#### Change 6c: Update variable names throughout readPSMs! (lines 74-120+)

**CURRENT CODE (line 79):**
```julia
            irt =  rt_irt(rts[row])
```

**CHANGED CODE:**
```julia
            refined_irt =  rt_irt(rts[row])  # rt_irt now contains rt_to_refined_irt spline
```

**CURRENT CODE (lines 84-86):**
```julia
            n = passed_q_val ? one(UInt16) : zero(UInt16)
            mean_irt = passed_q_val ? irt : zero(Float32)
            var_irt = zero(Float32)
```

**CHANGED CODE:**
```julia
            n = passed_q_val ? one(UInt16) : zero(UInt16)
            mean_refined_irt = passed_q_val ? refined_irt : zero(Float32)
            var_refined_irt = zero(Float32)
```

**CURRENT CODE (line 94):**
```julia
                best_prob, best_ms_file_idx, best_scan_idx, best_irt, old_mean_irt, var_irt, old_n, mz = prec_to_best_prob[precursor_idx]
```

**CHANGED CODE:**
```julia
                best_prob, best_ms_file_idx, best_scan_idx, best_refined_irt, old_mean_refined_irt, var_refined_irt, old_n, mz = prec_to_best_prob[precursor_idx]
```

**Continue replacing all instances of `irt`, `mean_irt`, `var_irt`, `best_irt` with `refined_irt`, `mean_refined_irt`, `var_refined_irt`, `best_refined_irt` throughout the function**

---

### 7. Modify create_rt_indices!() to Use Refined iRT

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

#### Change 7a: Update mapping (lines 569-571)

**CURRENT CODE:**
```julia
    # Create precursor to iRT mapping
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]),
                      precursor_dict)
```

**CHANGED CODE:**
```julia
    # Create precursor to refined iRT mapping
    prec_to_refined_irt = map(x -> (irt=x[:best_refined_irt], mz=x[:mz]),
                               precursor_dict)
```

#### Change 7b: Update makeRTIndices call (lines 601-607)

**CURRENT CODE:**
```julia
    # Make RT indices only for valid files
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        valid_psm_paths,
        prec_to_irt,
        valid_rt_models,
        min_prob=params.max_prob_to_impute
    )
```

**CHANGED CODE:**
```julia
    # Make RT indices using refined iRT
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        valid_psm_paths,
        prec_to_refined_irt,  # Now uses refined iRT
        valid_rt_models,
        min_prob=params.max_prob_to_impute
    )
```

---

### 8. Update buildRTIndex.jl to Use :refined_irt

**File:** `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl`

**Location:** Lines 80-109 inside makeRTIndices function

**CURRENT CODE (lines 87-90):**
```julia
        prec_set = Dict(zip(
            psms[:precursor_idx],
            map(x->(irt=first(x),prob=last(x)), zip(rt_to_irt.(psms[:rt]), psms[:prob]))
        ))
```

**CHANGED CODE:**
```julia
        # Map observed precursors to refined iRT (from rt_to_refined_irt spline) and probability
        prec_set = Dict(zip(
            psms[:precursor_idx],
            map(x->(irt=first(x),prob=last(x)), zip(rt_to_irt.(psms[:rt]), psms[:prob]))
        ))
```

**Note:** No actual code change needed here - the comment clarifies that `rt_to_irt` now contains the `rt_to_refined_irt` spline after our changes in map_retention_times!(). The variable naming convention in makeRTIndices parameter `prec_to_irt` refers to the refined iRT values passed from create_rt_indices!().

**CURRENT CODE (lines 92-105):**
```julia
        Threads.@threads for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
            prec_ids[i] = prec_id
            irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}
            #Don't impute irt, use empirical
            if haskey(prec_set, prec_id)
                _irt_, prob = prec_set[prec_id]
                if (prob >= min_prob)
                    irts[i], mzs[i]  = _irt_, mz
                    continue
                end
            end
            #Impute irt from the best observed psm for the precursor accross the experiment
            irts[i], mzs[i] = irt,mz
        end
```

**CHANGED CODE (update comments for clarity):**
```julia
        Threads.@threads for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
            prec_ids[i] = prec_id
            irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}
            # Use empirical refined iRT if available with high confidence
            if haskey(prec_set, prec_id)
                _irt_, prob = prec_set[prec_id]
                if (prob >= min_prob)
                    irts[i], mzs[i]  = _irt_, mz
                    continue
                end
            end
            # Otherwise use library refined iRT (from prec_to_irt mapping which contains refined values)
            irts[i], mzs[i] = irt,mz
        end
```

---

### 9. Update SecondPassSearch to Use Refined iRT

**File:** `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

**NO CODE CHANGES NEEDED**

SecondPassSearch already uses `rt_irt_map` from SearchContext, which we overwrite with `rt_to_refined_irt` splines in map_retention_times!() (Section 4, line 552). All subsequent searches automatically use refined iRT.

---

### 10. Update add_main_search_columns!() RT Column Names

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**Function:** `add_main_search_columns!` (starts around line 50)

#### Change 10a: Rename local variable (line 69)

**CURRENT CODE:**
```julia
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    irt = zeros(Float32, N);  # <-- LINE 69
    TIC = zeros(Float16, N);
```

**CHANGED CODE:**
```julia
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    observed_refined_irt = zeros(Float32, N);  # <-- RENAMED from irt
    TIC = zeros(Float16, N);
```

#### Change 10b: Update variable assignment (line 98)

**CURRENT CODE:**
```julia
                irt_pred[i] = Float32(prec_irt[prec_idx]);
                rt[i] = Float32(scan_retention_time[scan_idx[i]]);
                irt[i] = rt_irt(rt[i])  # <-- LINE 98
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
```

**CHANGED CODE:**
```julia
                irt_pred[i] = Float32(prec_irt[prec_idx]);
                rt[i] = Float32(scan_retention_time[scan_idx[i]]);
                observed_refined_irt[i] = rt_irt(rt[i])  # rt_irt now contains rt_to_refined_irt spline
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
```

**Note:** After our changes in map_retention_times!(), `rt_irt` parameter contains the `rt_to_refined_irt` spline (we overwrote rt_irt_map at line 552 in Section 4). So `rt_irt(rt[i])` now computes `observed_refined_irt`.

#### Change 10c: Update DataFrame column assignment (line 109)

**CURRENT CODE:**
```julia
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:irt] = irt  # <-- LINE 109
    psms[!,:TIC] = TIC
```

**CHANGED CODE:**
```julia
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:observed_refined_irt] = observed_refined_irt  # <-- RENAMED column
    psms[!,:TIC] = TIC
```

---

### 11. Update Imports and Exports

#### Change 11a: Add include statement

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**Location:** At top of file with other includes

**ADD:**
```julia
include("../../CommonSearchUtils/irt_refinement_utils.jl")
```

#### Change 11b: Add exports

**File:** `src/structs/RetentionTimeConversionModel.jl`

**Location:** After existing exports (if any) or at end of file

**ADD:**
```julia
export IrtRefinementModel, STANDARD_AAS
```

---

### 12. Add Debugging Print Statements

**Already included in code above:**

1. ✅ `fit_irt_refinement_model()` - Lines 328-338 in irt_refinement_utils.jl
2. ✅ `add_refined_irt_column!()` - Lines 410, 423, 430 in irt_refinement_utils.jl
3. ✅ `map_retention_times!()` - Lines 467, 506, 528, 570, 574 in utils.jl
4. **TODO**: Add @user_info in `create_rt_indices!()` to confirm using refined iRT

---

## Summary of All File Changes

| File | Lines | Changes | Before → After |
|------|-------|---------|----------------|
| `RetentionTimeConversionModel.jl` | 58+ | Add IrtRefinementModel struct | N/A → New code |
| `SearchTypes.jl` | 223 | Add field | `rt_irt_map` → `rt_irt_map, irt_refinement_models` |
| `SearchTypes.jl` | 267-269 | Initialize field | 2 Dicts → 3 Dicts |
| `SearchTypes.jl` | 524+ | Add getters/setters | N/A → New functions |
| `SearchTypes.jl` | 525 | Update signature | `:best_irt` → `:best_refined_irt` |
| `irt_refinement_utils.jl` | NEW | Create file | N/A → New file |
| `FirstPassSearch/utils.jl` | 345-473 | Rewrite function | Simple RT mapping → 5-step workflow |
| `FirstPassSearch/utils.jl` | 530-536 | Update type | `:best_irt` → `:best_refined_irt` |
| `FirstPassSearch/utils.jl` | 569-571 | Update mapping | `prec_to_irt` → `prec_to_refined_irt` |
| `FirstPassSearch/utils.jl` | 604 | Update call | `prec_to_irt` → `prec_to_refined_irt` |
| `FirstPassSearch/utils.jl` | 69 | Rename variable | `irt` → `observed_refined_irt` |
| `FirstPassSearch/utils.jl` | 98 | Update assignment | `irt[i]` → `observed_refined_irt[i]` |
| `FirstPassSearch/utils.jl` | 109 | Update column | `:irt` → `:observed_refined_irt` |
| `getBestPrecursorsAccrossRuns.jl` | Multiple | Update all refs | `irt` → `refined_irt` throughout |
| `buildRTIndex.jl` | 87-105 | Update comments | Clarify refined iRT usage |

---

## Testing Strategy

### Unit Tests
1. Test `IrtRefinementModel` callable interface with synthetic sequences
2. Test `prepare_features_dataframe()` creates 21 features correctly
3. Test `fit_irt_refinement_model()` with known error patterns
4. Test `add_refined_irt_column!()` adds correct column

### Integration Tests
1. Run FirstPassSearch on small dataset
2. Verify :refined_irt column exists in PSM files
3. Verify observed_refined_irt column exists
4. Verify refinement model stored in SearchContext
5. Verify RT indices use refined iRT
6. Run SecondPassSearch and verify it uses refined iRT from rt_irt_map

### Validation Checks
1. Compare `irt_error` distributions before/after refinement
2. Plot observed_refined_irt vs refined_irt scatter plots per file
3. Verify refinement only enabled when MAE improves
4. Check files with insufficient PSMs (<20) fallback to library iRT
5. Verify consistent terminology in all output columns

---

## Expected Outcomes

1. **Improved RT Accuracy**: Refined iRT reduces MAE for files where model validates successfully
2. **Fallback Safety**: Files with insufficient PSMs or no improvement use library iRT (refinement_model.use_refinement = false)
3. **Consistent Naming**:
   - `:irt_predicted` = library iRT (unchanged)
   - `:refined_irt` = library iRT with refinement applied
   - `:observed_refined_irt` = RT converted via rt_to_refined_irt spline
4. **Debugging Support**: Print statements confirm each step executes correctly
5. **Backward Compatibility**: Searches work with or without refinement enabled

---

## Next Steps After Approval

1. Create IrtRefinementModel struct in RetentionTimeConversionModel.jl
2. Update SearchContext fields and accessors
3. Create irt_refinement_utils.jl with helper functions
4. Modify map_retention_times!() with 5-step workflow
5. Update all field names from `irt` to `refined_irt`
6. Update get_best_precursors_accross_runs variable names
7. Update create_rt_indices!() to use refined_irt
8. Update add_main_search_columns!() variable and column names
9. Add imports and exports
10. Run integration tests
11. Generate validation plots
12. Document results
